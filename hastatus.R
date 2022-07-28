#-------------------------------------------------#
#---Ecological Niche Modeling of Rumex hastatus---#
#---MSc thesis / Martinee et al. 2022 Nord J Bot---#
#-------------------------------------------------#

####---START---####

rm(list = ls())

# packages
if(!require("sdm")){
  install.packages("sdm")
}
installAll() # only once (done)

library(randomForest)
library(sdm)
library(dismo)
library(rgbif)
library(usdm)
library(ecospat)
library(pROC)

library(raster)
library(countrycode)
library(rnaturalearthdata)
library(rnaturalearth)
library(ggmap)
library(OpenStreetMap)
library(sf)
library(sp)
library(terra)
library(CoordinateCleaner)
library(mapdata)
library(maps)
library(viridis)
library(ggspatial)
library(mapview)
library(rgdal)

library(ggplot2)
library(ggfortify)
library(ggnewscale)
library(patchwork)
library(cluster)
library(lfda)
library(matrixTests)
library(caTools)
library(reshape)
library(cluster)
library(MASS)
library(vegan)
library(mvnormtest)
library(rrcov)
library(precrec)
library(RColorBrewer)
library(cowplot)
library(scales)
library(dplyr)
library(tidyr)

library(plyr)

install.packages("")

citation("rgbif")
citation("raster")
citation("mapview")
citation("scales")
citation("ggnewscale")
help("")
#-------------------#

#------------------------#
#---1 Occurrence data-----

#-----------------------------#
#---1.1 GBIF and FIELD data----

# set R environment to access GBIF database
usethis::edit_r_environ() # GBIF_USER/PW/EMAIL

# get taxonKey for Rumex hastatus
hastatusKey <- name_backbone(name = "Rumex hastatus")
hastatusKey$speciesKey # 7739918
occ_count(taxonKey = hastatusKey$speciesKey) # records = 444

# retrieve occurrences (georeferrenced only)
hastatusGBIF <- occ_download(pred("taxonKey", 7739918), 
                             pred("hasCoordinate", TRUE))

occ_download_meta(hastatusGBIF) # Download key: 0023419-210914110416597

hastatus_occ <- occ_download_get("0023419-210914110416597", 
                                 path = "data/occ/", overwrite = T) %>% 
  occ_download_import()

# cite GBIF download
occ_download_meta(hastatusGBIF) %>% gbif_citation()

# R Open Sci pipeline:
# select columns of interest 
hastatus_occ <- hastatus_occ %>%
  dplyr::select(species, decimalLongitude, decimalLatitude, countryCode, individualCount,
                gbifID, family, taxonRank, coordinateUncertaintyInMeters, year,
                basisOfRecord, institutionCode, datasetName)

# remove records without coordinates
hastatus_occ <- hastatus_occ %>%
  filter(!is.na(decimalLongitude))%>%
  filter(!is.na(decimalLatitude))

# plot data to get an overview
wm <- borders("world", colour = "gray50", fill = "gray50")
ggplot() + coord_fixed() + wm +
  geom_point(data = hastatus_occ, aes(x = decimalLongitude, y = decimalLatitude),
             colour = "darkred", size = 0.75)+
  theme_bw()

# convert country code from ISO2c to ISO3c
hastatus_occ$countryCode <- countrycode(hastatus_occ$countryCode, 
                                        origin = 'iso2c', 
                                        destination = 'iso3c')

# flag problems
hastatus_occ <- data.frame(hastatus_occ)
?clean_coordinates
flags <- clean_coordinates(x = hastatus_occ,
                           lon = "decimalLongitude",
                           lat = "decimalLatitude",
                           countries = "countryCode",
                           species = "species",
                           tests = c("capitals", "centroids", "equal", "gbif", 
                                     "institutions", "outliers", "seas", "zeros"))

summary(flags)
plot(flags, lon = "decimalLongitude", lat = "decimalLatitude")

# flagged records ( 2 outliers)
hastatus_occ_flagged <- hastatus_occ[!flags$.summary,]

# remove flagged records
hastatus_occ_clean <- hastatus_occ[flags$.summary,]

# check for duplicates
names(hastatus_occ_clean)[2:3] <- c("decimallongitude", "decimallatitude") # rename

hastatus_occ_clean <- hastatus_occ_clean %>% # directly removes duplicates
  cc_dupl() # Removed 137 records.

names(hastatus_occ_clean)[2:3] <- c("decimalLongitude", "decimalLatitude") # rerename

# subset df
hastatus_occ_clean <- hastatus_occ_clean[, c(1, 2, 3)]

# write csv file
write.csv2(hastatus_occ_clean, "data/occ/hastatusGBIF_clean.csv", row.names = F)

# occurrences from field work (Georg Miehe & Joachim Schmidt)
hastatusFIELD <- read.csv2("data/occ/hastatusFIELD.csv")

# check for duplicates
names(hastatusFIELD)[2:3] <- c("decimallongitude", "decimallatitude") # rename

hastatusFIELD_clean <- hastatusFIELD %>% # directly removes duplicates
  cc_dupl()

names(hastatusFIELD_clean)[2:3] <- c("decimalLongitude", "decimalLatitude") # rerename

# plot data to get an overview
ggplot() + coord_fixed() + wm +
  geom_point(data = hastatusFIELD_clean, aes(x = decimalLongitude, y = decimalLatitude),
             colour = "darkred", size = 0.75) +
  theme_bw()

# write csv file
write.csv2(hastatusFIELD_clean, "data/occ/hastatusFIELD_clean.csv", row.names = F)

# merge cleaned occurrences from FIELD and GBIF
hastatus_all <- rbind(hastatus_occ_clean, hastatusFIELD_clean)

# write csv file
write.csv2(hastatus_all, "data/occ/hastatus_all.csv", row.names = F)
#---------------------------------------------------#

#------------------------#
#---1.2 Elevation data----

# let's get spatial
sp <- read.csv2("data/occ/hastatus_all.csv")
class(sp)
dim(sp)
sum(is.na(sp))

# manipulate df
sp$species <- 1 # 1 for presence (presence-only data)

# convert df to a spatial points df
coordinates(sp) <- c("decimalLongitude", "decimalLatitude")
#--------------------------------------------------------#

# srtm
srtm1 <- raster::getData("SRTM", lon = 65, lat = 25, path = "data/srtm/")
srtm2 <- raster::getData("SRTM", lon = 70, lat = 25, path = "data/srtm/")
srtm3 <- raster::getData("SRTM", lon = 75, lat = 25, path = "data/srtm/")
srtm4 <- raster::getData("SRTM", lon = 80, lat = 25, path = "data/srtm/")
srtm5 <- raster::getData("SRTM", lon = 85, lat = 25, path = "data/srtm/")
srtm6 <- raster::getData("SRTM", lon = 90, lat = 25, path = "data/srtm/")
srtm7 <- raster::getData("SRTM", lon = 95, lat = 25, path = "data/srtm/")
srtm8 <- raster::getData("SRTM", lon = 100, lat = 25, path = "data/srtm/")
srtm9 <- raster::getData("SRTM", lon = 105, lat = 25, path = "data/srtm/")
srtm10 <- raster::getData("SRTM", lon = 65, lat = 30, path = "data/srtm/")
srtm11 <- raster::getData("SRTM", lon = 70, lat = 30, path = "data/srtm/")
srtm12 <- raster::getData("SRTM", lon = 75, lat = 30, path = "data/srtm/")
srtm13 <- raster::getData("SRTM", lon = 80, lat = 30, path = "data/srtm/")
srtm14 <- raster::getData("SRTM", lon = 85, lat = 30, path = "data/srtm/")
srtm15 <- raster::getData("SRTM", lon = 90, lat = 30, path = "data/srtm/")
srtm16 <- raster::getData("SRTM", lon = 95, lat = 30, path = "data/srtm/")
srtm17 <- raster::getData("SRTM", lon = 100, lat = 30, path = "data/srtm/")
srtm18 <- raster::getData("SRTM", lon = 105, lat = 30, path = "data/srtm/")
srtm19 <- raster::getData("SRTM", lon = 65, lat = 35, path = "data/srtm/")
srtm20 <- raster::getData("SRTM", lon = 70, lat = 35, path = "data/srtm/")
srtm21 <- raster::getData("SRTM", lon = 75, lat = 35, path = "data/srtm/")
srtm22 <- raster::getData("SRTM", lon = 80, lat = 35, path = "data/srtm/")
srtm23 <- raster::getData("SRTM", lon = 85, lat = 35, path = "data/srtm/")
srtm24 <- raster::getData("SRTM", lon = 90, lat = 35, path = "data/srtm/")
srtm25 <- raster::getData("SRTM", lon = 95, lat = 35, path = "data/srtm/")
srtm26 <- raster::getData("SRTM", lon = 100, lat = 35, path = "data/srtm/")
srtm27 <- raster::getData("SRTM", lon = 105, lat = 35, path = "data/srtm/")
srtm28 <- raster::getData("SRTM", lon = 65, lat = 40, path = "data/srtm/")
srtm29 <- raster::getData("SRTM", lon = 70, lat = 40, path = "data/srtm/")
srtm30 <- raster::getData("SRTM", lon = 75, lat = 40, path = "data/srtm/")
srtm31 <- raster::getData("SRTM", lon = 80, lat = 40, path = "data/srtm/")
srtm32 <- raster::getData("SRTM", lon = 85, lat = 40, path = "data/srtm/")
srtm33 <- raster::getData("SRTM", lon = 90, lat = 40, path = "data/srtm/")
srtm34 <- raster::getData("SRTM", lon = 95, lat = 40, path = "data/srtm/")
srtm35 <- raster::getData("SRTM", lon = 100, lat = 40, path = "data/srtm/")
srtm36 <- raster::getData("SRTM", lon = 105, lat = 40, path = "data/srtm/")

srtm_merged <- raster::mosaic(srtm1, srtm2, srtm3, srtm4, srtm5, srtm6,
                              srtm7, srtm8, srtm9, srtm10, srtm11, srtm12,
                              srtm13, srtm14, srtm15, srtm16,
                              srtm17, srtm18, srtm19, srtm20, srtm21,
                              srtm22, srtm23, srtm24, srtm25, srtm26,
                              srtm27, srtm28, srtm29, srtm30, srtm31, srtm32,
                              srtm33, srtm34, srtm35, srtm36, fun = mean)

plot(srtm_merged)
plot(sp, add = T, col = "darkred")

raster::writeRaster(srtm_merged, "data/srtm/srtm_merged.tif", format = "GTiff",
                    overwrite = T) # write raster

# extract elevations from occ points using srtm data
elev <- raster::extract(srtm_merged, sp) # extract elev of pts
elev <- as.data.frame(elev) # make df

# merge with hastatus_all
hastatus_all_elev <- cbind(hastatus_all, elev)

# check
head(hastatus_all_elev)
range(hastatus_all_elev$elev)
mean(hastatus_all_elev$elev)
sd(hastatus_all_elev$elev)

# write csv file
write.csv2(hastatus_all_elev, "data/occ/hastatus_all_elev.csv", row.names = F)
#----------------------------------------------------------#


#-------------------#
#---1.3 Mapping-----

#------------------------#
#---1.3.1 Study area-----

# set theme
theme_set(theme_bw())

# occ, gap & elev data
sp <- read.csv2("data/occ/hastatus_all.csv")
gap.pts <- read.csv2("data/pca/d_gap_elev.csv", header = T)
gap.pts <- gap.pts[, c(1,2,5)]
head(gap.pts)

# river data 
rivers <- readOGR("data/hydro/mrb_shp_zip/mrb_named_rivers.shp")
head(rivers)
rivers@data$RIVER

my_rivers <- c("Salween", "Mekong", "Lancang Jiang", "Mekong (also Lancang Jiang)", "Yangtze", "Brahmaputra")
my_rivers2 <- c("Yarkant He", "Irrawaddy", "Min Jiang", "Yalong Jiang", "Yuan Jiang", "Red")

rivers_sa <- rivers[rivers@data$RIVER %in% my_rivers, ]
rivers.df <- fortify(rivers_sa)

rivers_sa2 <- rivers[rivers@data$RIVER %in% my_rivers2, ]
rivers.df2 <- fortify(rivers_sa2)

# crop elev data to sa extent
srtm_merged <- raster::raster("data/srtm/srtm_merged.tif") # load
extent_sa <- as(extent(66, 106.5, 20.8, 39), 'SpatialPolygons')
projection(extent_sa) <- projection(srtm_merged)
srtm_sa <- raster::crop(srtm_merged, extent_sa)

# extract points to df
srtm_sa_low <- aggregate(srtm_sa, fact = 10) # lower res for faster plotting
raster::writeRaster(srtm_sa_low, "data/srtm/srtm_sa_low.tif", format = "GTiff",
                    overwrite = T) # write raster

# make df
srtm_sa_low <- raster::raster("data/srtm/srtm_sa_low.tif") # load

srtm_pts <- rasterToPoints(srtm_sa_low)
srtm_df <- data.frame(srtm_pts)
colnames(srtm_df) <- c("lon", "lat", "elev")
head(srtm_df)

# world data
world <- ne_countries(scale = "medium", returnclass = "sf")
sf::sf_use_s2(FALSE) # turn off s2 processing

world_pts <- cbind(world, st_coordinates(st_centroid(world$geometry))) # country centroids
world_pts <- world_pts[c(-112), ] # remove Siachen Glacier (Kashmir) label

# labels & rectangle
him <- read.csv2("data/spatial/him.csv")
hdm <- read.csv2("data/spatial/hdm.csv")
tibet <- read.csv2("data/spatial/tibet.csv")
rect_sa <- read.csv2("data/spatial/rect_sa.csv")
rect_gap <- read.csv2("data/spatial/rect_gap.csv")
riv_lab <- read.csv2("data/hydro/rivers.csv")
riv_lab2 <- read.csv2("data/hydro/rivers2.csv")
places <- read.csv2("data/spatial/places.csv")
cities <- read.csv2("data/spatial/cities.csv")
cities_lab <- read.csv2("data/spatial/cities-lab.csv")

# plot sa (Figure 1)
sa <- ggplot() +
      geom_raster(data = srtm_df, aes(lon, lat, fill = elev)) +
      scale_fill_gradient(low = "grey90", high = "grey10") +
      geom_sf(data = world, col = "black", fill = NA, size = 0.1) +
      geom_path(data = rivers.df, aes(x = long, y = lat, group = group), 
                color = "dodgerblue", size = 0.4) +
      geom_text(data= world_pts, aes(x = X, y = Y, label = name),
                color = "black", fontface = "bold", check_overlap = FALSE, 
                size = 2) +
      geom_text(data = riv_lab, aes(lon, lat, label = river),
                color = "white", size = 2, fontface = "bold", angle = -45) +
      geom_text(data = tibet, aes(lon, lat, label = mtn),
                color = "black", size = 2.6, fontface = "bold") +
      geom_point(data = sp, aes(x = decimalLongitude, 
                                y = decimalLatitude),
                                size = 1.25, shape = 21,
                                color = "white", fill = "black") +
      annotation_scale(pad_x = unit(8.3, "cm"), 
                   pad_y = unit(5.9, "cm"),
                   height = unit(0.1, "cm"),
                   text_cex = 0.55) +
      annotation_north_arrow(which_north = "true", 
                         height = unit(.6, "cm"),
                         width = unit(.5, "cm"),
                         pad_x = unit(0.1, "cm"), 
                         pad_y = unit(0.1, "cm"),
                         style = north_arrow_fancy_orienteering) +
      coord_sf(xlim = c(66, 106.5), 
               ylim = c(20.8, 39), 
               expand = FALSE) +
      theme(legend.key.size = unit(0.4, 'cm'), 
            legend.key.height = unit(0.4, 'cm'), 
            legend.key.width = unit(0.4, 'cm'), 
            legend.title = element_text(size = 7.5, face = "bold"), 
            legend.text = element_text(size = 7.25),
            panel.background =  element_rect(fill = "aliceblue"),
            panel.grid.minor = element_blank(),
            panel.grid.major = element_blank(),
            axis.text.x = element_text(size = 8),
            axis.text.y = element_text(size = 8),
            axis.title = element_text(size = 8)) +
      labs(x = "Longitude", y = "Latitude", fill = "Elevation (m)")
sa

ggsave("fig1.pdf", plot = sa,      # save as pdf
       path = "plots/figs", device = "pdf",
       width = 16, height = 8 , units = c("cm"), dpi = 300)

ggsave("fig1.png", plot = sa,      # save as png
       path = "plots/figs/", device = "png",
       width = 16, height = 8 , units = c("cm"), dpi = 330)

# plot sa and gap pts (Figure S1)
sa2 <- ggplot() +
  geom_raster(data = srtm_df, aes(lon, lat, fill = elev)) +
  scale_fill_gradient(low = "grey90", high = "grey10") +
  geom_path(data = rect_sa, aes(x = lon, y = lat), 
            color = "black", size = .4, linetype = "dashed") +
  geom_sf(data = world, col = "black", fill = NA, size = 0.1) +
  geom_path(data = rivers.df, aes(x = long, y = lat, group = group), 
            color = "dodgerblue", size = 0.4) +
  geom_text(data= world_pts, aes(x = X, y = Y, label = name),
            color = "black", fontface = "bold", check_overlap = FALSE, 
            size = 2) +
  geom_text(data = riv_lab, aes(lon, lat, label = river),
            color = "white", size = 2, fontface = "bold", angle = -45) +
  geom_text(data = him, aes(lon, lat, label = mtn),
            color = "black", size = 3.6, fontface = "bold", angle = -45) +
  geom_text(data = hdm, aes(lon, lat, label = mtn),
            color = "black", size = 3.2, fontface = "bold", angle = 30) +
  geom_text(data = tibet, aes(lon, lat, label = mtn),
            color = "black", size = 3, fontface = "bold") +
  geom_path(data = rect_gap, aes(x = lon, y = lat), 
            color = "black", size = .4) +
  geom_point(data = gap.pts, aes(x = lon, 
                            y = lat),
             size = 1.25, shape = 21,
             color = "white", fill = "#e41a1c") +
  annotation_scale(pad_x = unit(8.35, "cm"), 
                   pad_y = unit(5.95, "cm"),
                   height = unit(0.1, "cm"),
                   text_cex = 0.55) +
  annotation_north_arrow(which_north = "true", 
                         height = unit(.6, "cm"),
                         width = unit(.5, "cm"),
                         pad_x = unit(0.1, "cm"), 
                         pad_y = unit(0.1, "cm"),
                         style = north_arrow_fancy_orienteering) +
  coord_sf(xlim = c(66, 106.5), 
           ylim = c(20.8, 39), 
           expand = FALSE) +
  theme(legend.key.size = unit(0.4, 'cm'), 
        legend.key.height = unit(0.4, 'cm'), 
        legend.key.width = unit(0.4, 'cm'), 
        legend.title = element_text(size = 7.5, face = "bold"), 
        legend.text = element_text(size = 7.25),
        panel.background =  element_rect(fill = "aliceblue"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 8)) +
  labs(x = "Longitude", y = "Latitude", fill = "Elevation (m)")
sa2

ggsave("fig_s1.png", plot = sa2,       # save as png
       path = "plots/suppl", device = "png",
       width = 16, height = 8 , units = c("cm"), dpi = 330)

# plot all places mentioned in manuscript (Figure S10)
sa10 <- ggplot() +
  geom_raster(data = srtm_df, aes(lon, lat, fill = elev)) +
  scale_fill_gradient(low = "grey90", high = "grey10") +
  geom_sf(data = world, col = "black", fill = NA, size = 0.1) +
  geom_path(data = rivers.df, aes(x = long, y = lat, group = group), 
            color = "dodgerblue", size = 0.4) +
  geom_path(data = rivers.df2, aes(x = long, y = lat, group = group), 
            color = "dodgerblue", size = 0.4) +
  geom_text(data= world_pts, aes(x = X, y = Y, label = name),
            color = "black", fontface = "bold", check_overlap = FALSE, 
            size = 2) +
  geom_text(data = riv_lab, aes(lon, lat, label = river),
            color = "white", size = 2, fontface = "bold", angle = -45) +
  geom_text(data = riv_lab2, aes(lon, lat, label = river),
            color = "white", size = 2, fontface = "bold", angle = -45) +
  geom_point(data = cities, aes(lon, lat),
            color = "black", fill = "white",  size = 1, shape = 21) +
  geom_text(data = cities_lab, aes(lon, lat, label = city_lab),
            color = "white", size = 1.6, fontface = "bold", angle = -30) +
  geom_text(data = places, aes(lon, lat, label = place),
            color = "black", size = 1.8, angle = -30) +
  annotation_scale(pad_x = unit(8.35, "cm"), 
                   pad_y = unit(5.95, "cm"),
                   height = unit(0.1, "cm"),
                   text_cex = 0.55) +
  annotation_north_arrow(which_north = "true", 
                         height = unit(.6, "cm"),
                         width = unit(.5, "cm"),
                         pad_x = unit(0.1, "cm"), 
                         pad_y = unit(0.1, "cm"),
                         style = north_arrow_fancy_orienteering) +
  coord_sf(xlim = c(66, 106.5), 
           ylim = c(20.8, 39), 
           expand = FALSE) +
  theme(legend.key.size = unit(0.4, 'cm'), 
        legend.key.height = unit(0.4, 'cm'), 
        legend.key.width = unit(0.4, 'cm'), 
        legend.title = element_text(size = 7.5, face = "bold"), 
        legend.text = element_text(size = 7.25),
        panel.background =  element_rect(fill = "aliceblue"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 8)) +
  labs(x = "Longitude", y = "Latitude", fill = "Elevation (m)")
sa10

ggsave("fig_s10.png", plot = sa10,      # save as png
       path = "plots/suppl", device = "png",
       width = 16, height = 8 , units = c("cm"), dpi = 330)
#-----------------------------------------------------------#


#---------------------------------#
#---2 Ecological Niche Modeling----

#-------------------------------------#
#---2.1 Rumex hastatus occurrences----

# load occurrences
sp <- read.csv2("data/occ/hastatus_all.csv")
dim(sp)

# manipulate df
sp$species <- 1 # 1 for presence (presence-only data)

# convert df to a spatial points df
coordinates(sp) <- c("decimalLongitude", "decimalLatitude")

# check coordinate reference system
proj4string(sp) # NA

#-----------------------------------------------------------#
#---2.2 ENM under current climatic conditions----

#---2.2.1 CHELSA climate data preparation----

## BIOCLIM 1979-2013 (CHELSA v1.2, res 30 arc sec)

# load CHELSA bioclim data on my machine
bio1 <- raster("D:/CHELSA/V1.2/present/CHELSA_bio10_01.tif")
bio2 <- raster("D:/CHELSA/V1.2/present/CHELSA_bio10_02.tif")
bio3 <- raster("D:/CHELSA/V1.2/present/CHELSA_bio10_03.tif")
bio4 <- raster("D:/CHELSA/V1.2/present/CHELSA_bio10_04.tif")
bio5 <- raster("D:/CHELSA/V1.2/present/CHELSA_bio10_05.tif")
bio6 <- raster("D:/CHELSA/V1.2/present/CHELSA_bio10_06.tif")
bio7 <- raster("D:/CHELSA/V1.2/present/CHELSA_bio10_07.tif")
bio8 <- raster("D:/CHELSA/V1.2/present/CHELSA_bio10_08.tif")
bio9 <- raster("D:/CHELSA/V1.2/present/CHELSA_bio10_09.tif")
bio10 <- raster("D:/CHELSA/V1.2/present/CHELSA_bio10_10.tif")
bio11 <- raster("D:/CHELSA/V1.2/present/CHELSA_bio10_11.tif")
bio12 <- raster("D:/CHELSA/V1.2/present/CHELSA_bio10_12.tif")
bio13 <- raster("D:/CHELSA/V1.2/present/CHELSA_bio10_13.tif")
bio14 <- raster("D:/CHELSA/V1.2/present/CHELSA_bio10_14.tif")
bio15 <- raster("D:/CHELSA/V1.2/present/CHELSA_bio10_15.tif")
bio16 <- raster("D:/CHELSA/V1.2/present/CHELSA_bio10_16.tif")
bio17 <- raster("D:/CHELSA/V1.2/present/CHELSA_bio10_17.tif")
bio18 <- raster("D:/CHELSA/V1.2/present/CHELSA_bio10_18.tif")
bio19 <- raster("D:/CHELSA/V1.2/present/CHELSA_bio10_19.tif")

# stack layers
biochelsa <- stack(bio1, bio2, bio3, bio4, bio5, bio6,
                   bio7, bio8, bio9, bio10, bio11, bio12,
                   bio13, bio14, bio15, bio16, bio17, bio18, bio19)

# WORLDCLIM data (version 2.1, res 2.5 min) FOR LAYER NAMES
bio <- raster::getData("worldclim", var = "bio", res = 2.5, 
                       path = "data/bioclim/")
bio # NOT FOR MODELLING BUT...

names(bio) # ...to use layer names to rename CHELSA layers

# rename layers
names(biochelsa)
names(biochelsa) <- names(bio)
names(biochelsa)
#--------------------------------------------------------#

# create spatial polygon to crop study area
x_coords <- c(66, 66, 75, 85,   106.5, 106.5, 92.6, 92.6, 75, 75, 66)
y_coords <- c(31, 39, 39, 32.6, 32.6,  20.8,  20.8, 26,   26, 31, 31)

polyg <- sp::Polygon(cbind(x_coords,y_coords)) # make polygon from matrix

polygo <- sp::Polygons(list(polyg), ID = "A") # polygon class
str(polygo, 1)
polygo

polygon <- sp::SpatialPolygons(list(polygo)) # spatial polygon
polygon

plot(polygon)

shapefile(x = polygon, 
          file = "data/spatial/polygon.shp",
          overwrite = T) # save as shapefile

# plot study area
plot(biochelsa[[1]], 
     xlim = c(65, 110), 
     ylim = c(20, 40))

plot(sp, add = T, col = "red") # add spp
plot(polygon, add = T) # add polygon

# save polygon as png
png("plots/map/polygon.png", width = 14.97, height = 8.5, 
    units = "cm", res = 600)
plot(polygon)
dev.off()

#---------------------------------------#

# define extent of study area
studyarea <- shapefile("data/spatial/polygon.shp")
studyarea

# set projection
projection(studyarea) 
projection(studyarea) <- projection(bioc)

# crop study area
e <- extent(studyarea)
bioc <- raster::crop(biochelsa, e) # crop raster object
bioc <- raster::mask(bioc, studyarea) # mask to remove NA values
bioc

plot(bioc[[1]])
plot(bioc[[12]])
#-------------------------------------------------------#

# harmonize layers (see CHELSA v1.2 tech specification)
bioc$bio1<- bioc$bio1/10
bioc$bio2<- bioc$bio2/10
bioc$bio3<- bioc$bio3/10
bioc$bio4<- bioc$bio4/10
bioc$bio5<- bioc$bio5/10
bioc$bio6<- bioc$bio6/10
bioc$bio7<- bioc$bio7/10
bioc$bio8<- bioc$bio8/10
bioc$bio9<- bioc$bio9/10
bioc$bio10<- bioc$bio10/10
bioc$bio11<- bioc$bio11/10

bioc

plot(bioc[[1]])
plot(bioc[[12]])

# write raster
raster::writeRaster(bioc, "data/bioclim/chelsa_bioc_studyarea.tif", format = "GTiff",
                    overwrite)
#-----------------------------------------------------------------------------------#

# read study area with chelsa bioclim variables
bioc <- raster::brick("data/bioclim/chelsa_bioc_studyarea.tif")

# WORLDCLIM data (version 2.1, res 2.5 min) FOR LAYER NAMES
bio <- raster::getData("worldclim", var = "bio", res = 2.5, 
                       path = "data/bioclim/")
bio # NOT FOR MODELLING BUT...
names(bio) # ...to use layer names to rename CHELSA layers

names(bioc) # check layer names
names(bioc) <- names(bio) # rename
names(bioc) 

# extract bioclim variables from occurence points
ex <- raster::extract(bioc, sp)
head(ex)
sum(is.na(ex)) # check for NAs (< 0 means NAs in df; not the case)

# check for multicolinearity
vc <- vifcor(ex)
vc # 12 variables 

vs <- vifstep(ex)
vs # 9 variables

# exclude problematic bioclim var (based on vs) from raster object
bioc <- exclude(bioc, vs)
bioc

# normalize data
bioc.norm <- raster::scale(bioc)
bioc.norm

raster::writeRaster(bioc.norm, "data/bioclim/bioc_normalized.tif", format = "GTiff",
                    overwrite = T) # write raster
#------------------------------------------------------------------#


#---2.2.2 Model fitting----

## Down-sampling RF method by Valavi & al. 2021

# load normalized data
bioc.norm <- raster::brick("data/bioclim/bioc_normalized.tif")

names(bioc.norm) # check layer names
names(bioc.norm) <- names(bioc) # rename
names(bioc.norm)
bioc.norm

# specify predictor variables and create sdmData object + 10000 background points
#?sdmData
set.seed(555)
d <- sdmData(species ~., train = sp, predictors = bioc.norm,
             bg = list(method = "gRandom", n = 10000)) # ~. selects all predictor var in bioc
d

# make sdm data object a df; species 1 for presence, 0 for absence (bg pts)
ddata <- as.data.frame(d)
head(ddata)
ddata <- ddata[, -1] # drop 1st column

# convert response to factor for classification
ddata$species <- as.factor(ddata$species)

# split data in training and test dataset / repeated split sample
percent <- floor((nrow(ddata)/10)*7) # 70 % for training / 30 % for testing

# define two lists containing 25 dfs
ddata.train <- list()
ddata.test <- list()

# for loop to subsample 25 different df for training and testing
set.seed(555) # make reproducible
for(k in 1:25){
  ddata <- ddata[sample(nrow(ddata)), ] # sample rows
  ddata.train[[k]] <- ddata[1:percent, ]
  ddata.test[[k]] <- ddata[(percent+1):nrow(ddata), ]
}

# check dfs -> differ in order of samples (muy bien)
head(ddata.train[[1]])
head(ddata.train[[25]])
head(ddata.test[[1]])
head(ddata.test[[25]])

# calculate sample sizes (always sample same amount of presences and pseudo-absences)
spsize <- list() # define list
for(k in 1:25){
  prNums <- as.numeric(table(ddata.train[[k]]$species)["1"]) # number of presence records
  spsize[[k]] <- c("0" = prNums, "1" = prNums) 
}

spsize # check list

# for loop to fit 25 rf down-sampling models and put them in a list
rf_down <- list()
set.seed(555)
for(k in 1:25){
rf_down[[k]] <- randomForest(species ~ .,
                              data = ddata.train[[k]],
                              ntree = 1000,
                              sampsize = spsize[[k]],
                              importance = TRUE,
                              replace = TRUE) # reduce overfitting
}

rf_down # check

plot(rf_down[[1]]) # plot error against no of trees
plot(rf_down[[25]])
#------------------------------------------------------------------------#

#---2.2.3 Prediction and averaging----

# load normalized data
bioc.norm <- raster::brick("data/bioclim/bioc_normalized.tif")

names(bioc.norm) # check layer names
names(bioc.norm) <- names(bioc) # rename
names(bioc.norm)

# generate predictions using all 25 fitted models (pp = prediction present)
pp <- list()
for(k in 1:25){
  pp[[k]] <- predict(bioc.norm, rf_down[[k]], type = "prob", index = 2)
}
pp

# make raster brick
pp.brick <- raster::brick(pp)

raster::writeRaster(pp.brick, "output/predict/pp_rf-down.tif", format = "GTiff",
                    overwrite = T) # write raster

# calculate mean over all raster layers
pp.avg <- calc(pp.brick, fun = mean)
pp.avg

raster::writeRaster(pp.avg, "output/ensemble/pp_avg_rf-down.tif", format = "GTiff",
                    overwrite = T) # write raster

#---2.2.3 Plotting----

# load raster data
pp.avg <- raster("output/ensemble/pp_avg_rf-down.tif")

# make df
pp.pts <- rasterToPoints(pp.avg)
pp.df <- data.frame(pp.pts)
colnames(pp.df) <- c("lon", "lat", "niche")
head(pp.df)

# enm on sa map
f1 <- ggplot() +
  geom_raster(data = srtm_df, aes(lon, lat, fill = elev)) +
  scale_fill_gradient(low = "grey90", high = "grey10") +
  geom_sf(data = world, col = "black", fill = NA, size = 0.1) +
  guides(fill = "none") +
  new_scale_fill() +
  geom_raster(data = pp.df, aes(lon, lat, fill = niche)) +
  scale_fill_viridis_c(option = "mako") +
  geom_path(data = rivers.df, aes(x = long, y = lat, group = group), 
            color = "dodgerblue", size = 0.3) +
  geom_text(data = riv_lab, aes(lon, lat, label = river),
            color = "white", size = 2, fontface = "bold", angle = -45) +
  annotation_scale(pad_x = unit(8, "cm"), 
                   pad_y = unit(5.7, "cm"),
                   height = unit(0.1, "cm"),
                   text_cex = 0.55) +
  annotation_north_arrow(which_north = "true", 
                         height = unit(.6, "cm"),
                         width = unit(.5, "cm"),
                         pad_x = unit(0.1, "cm"), 
                         pad_y = unit(0.1, "cm"),
                         style = north_arrow_fancy_orienteering) +
  coord_sf(xlim = c(66, 106.5), 
           ylim = c(20.8, 39), 
           expand = FALSE) +
  theme(legend.key.size = unit(0.4, 'cm'), 
        legend.key.height = unit(0.4, 'cm'), 
        legend.key.width = unit(0.4, 'cm'), 
        legend.title = element_text(size = 7.5, face = "bold"), 
        legend.text = element_text(size = 7.25),
        panel.background =  element_rect(fill = "aliceblue"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 8)) +
  labs(x = "Longitude", y = "Latitude", fill = "Habitat\nsuitability")
f1

ggsave("pp_map.pdf", plot = f1,      # save as pdf
       path = "plots/sdm/publ", device = "pdf",
       width = 16, height = 12 , units = c("cm"), dpi = 300)
#----------------------------------------------------------------------------------#


#---------------------------------------------#
#---2.3 ENM under past climatic conditions----

#---2.3.1 CHELSA climate data preparation----

## LGM ~21,000 ya / PMIP3 / CCSM4 (CHELSA v1.2, res 30 arc-sec)
# load CHELSA bioclim LGM data
lgm1 <- raster("D:/CHELSA/V1.2/lgm/CCSM4/CHELSA_PMIP_CCSM4_BIO_01.tif")
lgm2 <- raster("D:/CHELSA/V1.2/lgm/CCSM4/CHELSA_PMIP_CCSM4_BIO_02.tif")
lgm3 <- raster("D:/CHELSA/V1.2/lgm/CCSM4/CHELSA_PMIP_CCSM4_BIO_03.tif")
lgm4 <- raster("D:/CHELSA/V1.2/lgm/CCSM4/CHELSA_PMIP_CCSM4_BIO_04.tif")
lgm5 <- raster("D:/CHELSA/V1.2/lgm/CCSM4/CHELSA_PMIP_CCSM4_BIO_05.tif")
lgm6 <- raster("D:/CHELSA/V1.2/lgm/CCSM4/CHELSA_PMIP_CCSM4_BIO_06.tif")
lgm7 <- raster("D:/CHELSA/V1.2/lgm/CCSM4/CHELSA_PMIP_CCSM4_BIO_07.tif")
lgm8 <- raster("D:/CHELSA/V1.2/lgm/CCSM4/CHELSA_PMIP_CCSM4_BIO_08.tif")
lgm9 <- raster("D:/CHELSA/V1.2/lgm/CCSM4/CHELSA_PMIP_CCSM4_BIO_09.tif")
lgm10 <- raster("D:/CHELSA/V1.2/lgm/CCSM4/CHELSA_PMIP_CCSM4_BIO_10.tif")
lgm11 <- raster("D:/CHELSA/V1.2/lgm/CCSM4/CHELSA_PMIP_CCSM4_BIO_11.tif")
lgm12 <- raster("D:/CHELSA/V1.2/lgm/CCSM4/CHELSA_PMIP_CCSM4_BIO_12.tif")
lgm13 <- raster("D:/CHELSA/V1.2/lgm/CCSM4/CHELSA_PMIP_CCSM4_BIO_13.tif")
lgm14 <- raster("D:/CHELSA/V1.2/lgm/CCSM4/CHELSA_PMIP_CCSM4_BIO_14.tif")
lgm15 <- raster("D:/CHELSA/V1.2/lgm/CCSM4/CHELSA_PMIP_CCSM4_BIO_15.tif")
lgm16 <- raster("D:/CHELSA/V1.2/lgm/CCSM4/CHELSA_PMIP_CCSM4_BIO_16.tif")
lgm17 <- raster("D:/CHELSA/V1.2/lgm/CCSM4/CHELSA_PMIP_CCSM4_BIO_17.tif")
lgm18 <- raster("D:/CHELSA/V1.2/lgm/CCSM4/CHELSA_PMIP_CCSM4_BIO_18.tif")
lgm19 <- raster("D:/CHELSA/V1.2/lgm/CCSM4/CHELSA_PMIP_CCSM4_BIO_19.tif")

# stack layers
lgmchelsa <- stack(lgm1, lgm2, lgm3, lgm4, lgm5, lgm6, lgm7,
                   lgm8, lgm9, lgm10, lgm11, lgm12, lgm13, lgm14,
                   lgm15, lgm16, lgm17, lgm18, lgm19) 

# WORLDCLIM data (version 2.1, res 2.5 min) FOR LAYER NAMES
bio <- raster::getData("worldclim", var = "bio", res = 2.5, 
                       path = "data/bioclim/")
bio # NOT FOR MODELLING BUT...
names(bio) # ...to use layer names to rename CHELSA layers

# change names of biof so the model can understand
names(lgmchelsa)
names(lgmchelsa) <- names(bio)
names(lgmchelsa)
#-------------------------------------------------#

# define extent of study area
studyarea <- shapefile("data/spatial/polygon.shp")
studyarea

# set projection
projection(lgmchelsa) 
projection(studyarea) <- projection(lgmchelsa)

# crop study area
e <- extent(studyarea)
lgmc <- raster::crop(lgmchelsa, e) # crop raster object
lgmc <- raster::mask(lgmc, studyarea) # mask to remove NA values
lgmc

plot(lgmc[[1]])
plot(lgmc[[12]])
#-------------------------------------------------------#

# harmonize layer units
lgmc$bio1 <- (lgmc$bio1/10)-273.15
lgmc$bio2 <- (lgmc$bio2/10)#-273.15
lgmc$bio3 <- (lgmc$bio3/10)#-273.15
lgmc$bio4 <- (lgmc$bio4/10)#-273.15
lgmc$bio5 <- (lgmc$bio5/10)-273.15
lgmc$bio6 <- (lgmc$bio6/10)-273.15
lgmc$bio7 <- (lgmc$bio7/10)#-273.15
lgmc$bio8 <- (lgmc$bio8/10)-273.15
lgmc$bio9 <- (lgmc$bio9/10)-273.15
lgmc$bio10 <- (lgmc$bio10/10)-273.15
lgmc$bio11 <- (lgmc$bio11/10)-273.15

lgmc$bio12 <- lgmc$bio12/10
lgmc$bio13 <- lgmc$bio13/10
lgmc$bio14 <- lgmc$bio14/10
lgmc$bio15 <- lgmc$bio15/10
lgmc$bio16 <- lgmc$bio16/10
lgmc$bio17 <- lgmc$bio17/10
lgmc$bio18 <- lgmc$bio18/10
lgmc$bio19 <- lgmc$bio19/10

lgmc

plot(lgmc[[1]])
plot(lgmc[[12]])

# write raster
raster::writeRaster(lgmc, "data/bioclim/chelsa_lgmc_studyarea.tif", format = "GTiff",
                    overwrite = T)
#------------------------------------------------------------------------------------#

# read study area with chelsa bioclim variables
lgmc <- raster::brick("data/bioclim/chelsa_lgmc_studyarea.tif")
lgmc

names(lgmc) # check layer names
names(lgmc) <- names(bio) # rename
names(lgmc)

# exclude same multicolinear var as in bioc
lgmc <- exclude(lgmc, vs)
lgmc

# normalize data
lgmc.norm <- raster::scale(lgmc)
lgmc.norm

raster::writeRaster(lgmc.norm, "data/bioclim/lgmc_normalized.tif", format = "GTiff",
                    overwrite = T) # write raster
#-----------------------------------------------------------------------------------#

#---2.3.2 Prediction and averaging----

# load normalized data
lgmc.norm <- raster::brick("data/bioclim/lgmc_normalized.tif")

names(lgmc.norm) # check layer names
names(lgmc.norm) <- names(bioc) # rename
names(lgmc.norm)
lgmc.norm

# generate predictions using all 25 fitted models (pp = prediction lgm)
plgm <- list()
for(k in 1:25){
  plgm[[k]] <- predict(lgmc.norm, rf_down[[k]], type = "prob", index = 2)
}
plgm

# make raster brick
plgm.brick <- raster::brick(plgm)

raster::writeRaster(plgm.brick, "output/predict/plgm_rf-down.tif", format = "GTiff",
                    overwrite = T) # write raster

# calculate mean over all raster layers
plgm.avg <- calc(plgm.brick, fun = mean)
plgm.avg

raster::writeRaster(plgm.avg, "output/ensemble/plgm_avg_rf-down.tif", format = "GTiff",
                    overwrite = T) # write raster


#---2.3.3 Plotting----

# load raster data
plgm.avg <- raster("output/ensemble/plgm_avg_rf-down.tif")

# make df
plgm.pts <- rasterToPoints(plgm.avg)
plgm.df <- data.frame(plgm.pts)
colnames(plgm.df) <- c("lon", "lat", "niche")
head(plgm.df)

# enm on sa map
label_lgm <- read.csv2("data/spatial/label_lgm.csv")

f2 <- ggplot() +
  geom_raster(data = srtm_df, aes(lon, lat, fill = elev)) +
  scale_fill_gradient(low = "grey90", high = "grey10") +
  geom_sf(data = world, col = "black", fill = NA, size = 0.1) +
  guides(fill = "none") +
  new_scale_fill() +
  geom_raster(data = plgm.df, aes(lon, lat, fill = niche)) +
  scale_fill_viridis_c(option = "mako") +
  guides(fill = "none") +
  geom_path(data = rivers.df, aes(x = long, y = lat, group = group), 
            color = "dodgerblue", size = 0.3) +
  geom_text(data = label_lgm, aes(lon, lat, label = period),
            color = "white", size = 2.75, fontface = "bold") +
  coord_sf(xlim = c(66, 106.5), 
           ylim = c(20.8, 39), 
           expand = FALSE) +
  theme(legend.key.size = unit(0.3, 'cm'), 
        legend.key.height = unit(0.3, 'cm'), 
        legend.key.width = unit(0.3, 'cm'), 
        legend.title = element_text(size = 6, face = "bold"), 
        legend.text = element_text(size = 6),
        panel.background =  element_rect(fill = "aliceblue"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title = element_blank(),
        axis.ticks = element_line(size = .3)) +
  labs(x = "", y = "", fill = "Habitat\nsuitability") 
f2

ggsave("plgm_map.pdf", plot = f2,      # save as pdf
       path = "plots/sdm/publ", device = "pdf",
       width = 16, height = 12 , units = c("cm"), dpi = 300)
#---------------------------------------------#

#---2.4 ENM under future climatic conditions----

#---2.4.1 RCP 26----

## 2061-2080 / CMIP5 / CCSM4 (CHELSA v1.2, res 30 arc sec)

# load CHELSA bioclim data 
bio261 <- raster("D:/CHELSA/V1.2/future/CCSM4/CHELSA_bio_mon_CCSM4_rcp26_r1i1p1_g025.nc_1_2061-2080_V1.2.tif")
bio262 <- raster("D:/CHELSA/V1.2/future/CCSM4/CHELSA_bio_mon_CCSM4_rcp26_r1i1p1_g025.nc_2_2061-2080_V1.2.tif")
bio263 <- raster("D:/CHELSA/V1.2/future/CCSM4/CHELSA_bio_mon_CCSM4_rcp26_r1i1p1_g025.nc_3_2061-2080_V1.2.tif")
bio264 <- raster("D:/CHELSA/V1.2/future/CCSM4/CHELSA_bio_mon_CCSM4_rcp26_r1i1p1_g025.nc_4_2061-2080_V1.2.tif")
bio265 <- raster("D:/CHELSA/V1.2/future/CCSM4/CHELSA_bio_mon_CCSM4_rcp26_r1i1p1_g025.nc_5_2061-2080_V1.2.tif")
bio266 <- raster("D:/CHELSA/V1.2/future/CCSM4/CHELSA_bio_mon_CCSM4_rcp26_r1i1p1_g025.nc_6_2061-2080_V1.2.tif")
bio267 <- raster("D:/CHELSA/V1.2/future/CCSM4/CHELSA_bio_mon_CCSM4_rcp26_r1i1p1_g025.nc_7_2061-2080_V1.2.tif")
bio268 <- raster("D:/CHELSA/V1.2/future/CCSM4/CHELSA_bio_mon_CCSM4_rcp26_r1i1p1_g025.nc_8_2061-2080_V1.2.tif")
bio269 <- raster("D:/CHELSA/V1.2/future/CCSM4/CHELSA_bio_mon_CCSM4_rcp26_r1i1p1_g025.nc_9_2061-2080_V1.2.tif")
bio2610 <- raster("D:/CHELSA/V1.2/future/CCSM4/CHELSA_bio_mon_CCSM4_rcp26_r1i1p1_g025.nc_10_2061-2080_V1.2.tif")
bio2611 <- raster("D:/CHELSA/V1.2/future/CCSM4/CHELSA_bio_mon_CCSM4_rcp26_r1i1p1_g025.nc_11_2061-2080_V1.2.tif")
bio2612 <- raster("D:/CHELSA/V1.2/future/CCSM4/CHELSA_bio_mon_CCSM4_rcp26_r1i1p1_g025.nc_12_2061-2080_V1.2.tif")
bio2613 <- raster("D:/CHELSA/V1.2/future/CCSM4/CHELSA_bio_mon_CCSM4_rcp26_r1i1p1_g025.nc_13_2061-2080_V1.2.tif")
bio2614 <- raster("D:/CHELSA/V1.2/future/CCSM4/CHELSA_bio_mon_CCSM4_rcp26_r1i1p1_g025.nc_14_2061-2080_V1.2.tif")
bio2615 <- raster("D:/CHELSA/V1.2/future/CCSM4/CHELSA_bio_mon_CCSM4_rcp26_r1i1p1_g025.nc_15_2061-2080_V1.2.tif")
bio2616 <- raster("D:/CHELSA/V1.2/future/CCSM4/CHELSA_bio_mon_CCSM4_rcp26_r1i1p1_g025.nc_16_2061-2080_V1.2.tif")
bio2617 <- raster("D:/CHELSA/V1.2/future/CCSM4/CHELSA_bio_mon_CCSM4_rcp26_r1i1p1_g025.nc_17_2061-2080_V1.2.tif")
bio2618 <- raster("D:/CHELSA/V1.2/future/CCSM4/CHELSA_bio_mon_CCSM4_rcp26_r1i1p1_g025.nc_18_2061-2080_V1.2.tif")
bio2619 <- raster("D:/CHELSA/V1.2/future/CCSM4/CHELSA_bio_mon_CCSM4_rcp26_r1i1p1_g025.nc_19_2061-2080_V1.2.tif")

# stack layers
bio26 <- stack(bio261, bio262, bio263, bio264, bio265, bio266,
               bio267, bio268, bio269, bio2610, bio2611, bio2612,
               bio2613, bio2614, bio2615, bio2616, bio2617, bio2618, bio2619)

# change names of biof so the model can understand
names(bio26)
names(bio26) <- names(bio)
names(bio26)
#------------------------------------------------#

# define extent of study area
studyarea <- shapefile("data/spatial/polygon.shp")
studyarea

# set projection
projection(bio26) 
projection(bio26) <- projection(bioc)

# crop study area
e <- extent(studyarea)
bioc26 <- raster::crop(bio26, e) # crop raster object
bioc26 <- raster::mask(bioc26, studyarea) # mask to remove NA values
bioc26

plot(bioc26[[1]])
plot(bioc26[[12]])
#-----------------------------------------------------------#

# harmonize layers (see CHELSA v1.2 tech specification)
bioc26$bio1<- bioc26$bio1/10
bioc26$bio2<- bioc26$bio2/10
bioc26$bio3<- bioc26$bio3/10
bioc26$bio4<- bioc26$bio4/10
bioc26$bio5<- bioc26$bio5/10
bioc26$bio6<- bioc26$bio6/10
bioc26$bio7<- bioc26$bio7/10
bioc26$bio8<- bioc26$bio8/10
bioc26$bio9<- bioc26$bio9/10
bioc26$bio10<- bioc26$bio10/10
bioc26$bio11<- bioc26$bio11/10

bioc26

plot(bioc26[[1]])
plot(bioc26[[12]])

# write raster
raster::writeRaster(bioc26, "data/bioclim/chelsa_bioc26_studyarea.tif", format = "GTiff",
                    overwrite = T)
#------------------------------------------------------------------------------------#

# read study area with chelsa bioclim variables
bioc26 <- raster::brick("data/bioclim/chelsa_bioc26_studyarea.tif")
bioc26

names(bioc26) # check layer names
names(bioc26) <- names(bio) # rename
names(bioc26)

# exclude same multicolinear var as in bioc
bioc26 <- exclude(bioc26, vs)
bioc26

# normalize data
bioc26.norm <- raster::scale(bioc26)
bioc26.norm

raster::writeRaster(bioc26.norm, "data/bioclim/bioc26_normalized.tif", format = "GTiff",
                    overwrite = T) # write raster
#----------------------------------------------------------------------------------------#

#---2.4.1.1 Prediction and averaging----

# load normalized data
bioc26.norm <- raster::brick("data/bioclim/bioc26_normalized.tif")

names(bioc26.norm) # check layer names
names(bioc26.norm) <- names(bioc) # rename
names(bioc26.norm)
bioc26.norm

# generate predictions using all 25 fitted models (p26 = prediction rcp 26)
p26 <- list()
for(k in 1:25){
  p26[[k]] <- predict(bioc26.norm, rf_down[[k]], type = "prob", index = 2)
}
p26

# make raster brick
p26.brick <- raster::brick(p26)

raster::writeRaster(p26.brick, "output/predict/p26_rf-down.tif", format = "GTiff",
                    overwrite = T)

# calculate mean over all raster layers
p26.avg <- calc(p26.brick, fun = mean)
p26.avg

raster::writeRaster(p26.avg, "output/ensemble/p26_avg_rf-down.tif", format = "GTiff",
                    overwrite = T) # write raster


#---2.4.1.2 Plotting----

# load raster data
p26.avg <- raster("output/ensemble/p26_avg_rf-down.tif")

# make df
p26.pts <- rasterToPoints(p26.avg)
p26.df <- data.frame(p26.pts)
colnames(p26.df) <- c("lon", "lat", "niche")
head(p26.df)

# enm on sa map
label_26 <- read.csv2("data/spatial/label_26.csv")

f3 <- ggplot() +
  geom_raster(data = srtm_df, aes(lon, lat, fill = elev)) +
  scale_fill_gradient(low = "grey90", high = "grey10") +
  geom_sf(data = world, col = "black", fill = NA, size = 0.1) +
  guides(fill = "none") +
  new_scale_fill() +
  geom_raster(data = p26.df, aes(lon, lat, fill = niche)) +
  scale_fill_viridis_c(option = "mako") +
  geom_path(data = rivers.df, aes(x = long, y = lat, group = group), 
            color = "dodgerblue", size = 0.3) +
  geom_text(data = riv_lab, aes(lon, lat, label = river),
            color = "white", size = 2, fontface = "bold", angle = -45) +
  geom_text(data = label_26, aes(lon, lat, label = period),
            color = "black", size = 5, fontface = "bold") +
  annotation_scale(pad_x = unit(8.7, "cm"), 
                   pad_y = unit(6.2, "cm"),
                   height = unit(0.1, "cm"),
                   text_cex = 0.55) +
  annotation_north_arrow(which_north = "true", 
                         height = unit(.6, "cm"),
                         width = unit(.5, "cm"),
                         pad_x = unit(0.1, "cm"), 
                         pad_y = unit(0.1, "cm"),
                         style = north_arrow_fancy_orienteering) +
  coord_sf(xlim = c(66, 106.5), 
           ylim = c(20.8, 39), 
           expand = FALSE) +
  theme(legend.key.size = unit(0.4, 'cm'), 
        legend.key.height = unit(0.4, 'cm'), 
        legend.key.width = unit(0.4, 'cm'), 
        legend.title = element_text(size = 7.5, face = "bold"), 
        legend.text = element_text(size = 7.25),
        panel.background =  element_rect(fill = "aliceblue"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 8)) +
  labs(x = "Longitude", y = "Latitude", fill = "Habitat\nsuitability") 
f3

ggsave("fig_s4.png", plot = f3,      # save as png
       path = "plots/suppl", device = "png",
       width = 16, height = 8 , units = c("cm"), dpi = 330)
#----------------------------------------------#

#---2.4.2 RCP 45----

## 2061-2080 / CMIP5 / CCSM4 (CHELSA v1.2, res 30 arc sec)

# load CHELSA bioclim data 
bio451 <- raster("D:/CHELSA/V1.2/future/CCSM4/CHELSA_bio_mon_CCSM4_rcp45_r1i1p1_g025.nc_1_2061-2080_V1.2.tif")
bio452 <- raster("D:/CHELSA/V1.2/future/CCSM4/CHELSA_bio_mon_CCSM4_rcp45_r1i1p1_g025.nc_2_2061-2080_V1.2.tif")
bio453 <- raster("D:/CHELSA/V1.2/future/CCSM4/CHELSA_bio_mon_CCSM4_rcp45_r1i1p1_g025.nc_3_2061-2080_V1.2.tif")
bio454 <- raster("D:/CHELSA/V1.2/future/CCSM4/CHELSA_bio_mon_CCSM4_rcp45_r1i1p1_g025.nc_4_2061-2080_V1.2.tif")
bio455 <- raster("D:/CHELSA/V1.2/future/CCSM4/CHELSA_bio_mon_CCSM4_rcp45_r1i1p1_g025.nc_5_2061-2080_V1.2.tif")
bio456 <- raster("D:/CHELSA/V1.2/future/CCSM4/CHELSA_bio_mon_CCSM4_rcp45_r1i1p1_g025.nc_6_2061-2080_V1.2.tif")
bio457 <- raster("D:/CHELSA/V1.2/future/CCSM4/CHELSA_bio_mon_CCSM4_rcp45_r1i1p1_g025.nc_7_2061-2080_V1.2.tif")
bio458 <- raster("D:/CHELSA/V1.2/future/CCSM4/CHELSA_bio_mon_CCSM4_rcp45_r1i1p1_g025.nc_8_2061-2080_V1.2.tif")
bio459 <- raster("D:/CHELSA/V1.2/future/CCSM4/CHELSA_bio_mon_CCSM4_rcp45_r1i1p1_g025.nc_9_2061-2080_V1.2.tif")
bio4510 <- raster("D:/CHELSA/V1.2/future/CCSM4/CHELSA_bio_mon_CCSM4_rcp45_r1i1p1_g025.nc_10_2061-2080_V1.2.tif")
bio4511 <- raster("D:/CHELSA/V1.2/future/CCSM4/CHELSA_bio_mon_CCSM4_rcp45_r1i1p1_g025.nc_11_2061-2080_V1.2.tif")
bio4512 <- raster("D:/CHELSA/V1.2/future/CCSM4/CHELSA_bio_mon_CCSM4_rcp45_r1i1p1_g025.nc_12_2061-2080_V1.2.tif")
bio4513 <- raster("D:/CHELSA/V1.2/future/CCSM4/CHELSA_bio_mon_CCSM4_rcp45_r1i1p1_g025.nc_13_2061-2080_V1.2.tif")
bio4514 <- raster("D:/CHELSA/V1.2/future/CCSM4/CHELSA_bio_mon_CCSM4_rcp45_r1i1p1_g025.nc_14_2061-2080_V1.2.tif")
bio4515 <- raster("D:/CHELSA/V1.2/future/CCSM4/CHELSA_bio_mon_CCSM4_rcp45_r1i1p1_g025.nc_15_2061-2080_V1.2.tif")
bio4516 <- raster("D:/CHELSA/V1.2/future/CCSM4/CHELSA_bio_mon_CCSM4_rcp45_r1i1p1_g025.nc_16_2061-2080_V1.2.tif")
bio4517 <- raster("D:/CHELSA/V1.2/future/CCSM4/CHELSA_bio_mon_CCSM4_rcp45_r1i1p1_g025.nc_17_2061-2080_V1.2.tif")
bio4518 <- raster("D:/CHELSA/V1.2/future/CCSM4/CHELSA_bio_mon_CCSM4_rcp45_r1i1p1_g025.nc_18_2061-2080_V1.2.tif")
bio4519 <- raster("D:/CHELSA/V1.2/future/CCSM4/CHELSA_bio_mon_CCSM4_rcp45_r1i1p1_g025.nc_19_2061-2080_V1.2.tif")

# stack layers
bio45 <- stack(bio451, bio452, bio453, bio454, bio455, bio456,
               bio457, bio458, bio459, bio4510, bio4511, bio4512,
               bio4513, bio4514, bio4515, bio4516, bio4517, bio4518, bio4519)

# change names of biof so the model can understand
names(bio45)
names(bio45) <- names(bio)
names(bio45)
#------------------------------------------------#

# define extent of study area
studyarea <- shapefile("data/spatial/polygon.shp")
studyarea

# set projection
projection(bio45) 
projection(bio45) <- projection(bioc)

# crop study area
e <- extent(studyarea)
bioc45 <- raster::crop(bio45, e) # crop raster object
bioc45 <- raster::mask(bioc45, studyarea) # mask to remove NA values
bioc45

plot(bioc45[[1]])
plot(bioc45[[12]])
#-------------------------------------------------------#

# harmonize layers (see CHELSA v1.2 tech specification)
bioc45$bio1<- bioc45$bio1/10
bioc45$bio2<- bioc45$bio2/10
bioc45$bio3<- bioc45$bio3/10
bioc45$bio4<- bioc45$bio4/10
bioc45$bio5<- bioc45$bio5/10
bioc45$bio6<- bioc45$bio6/10
bioc45$bio7<- bioc45$bio7/10
bioc45$bio8<- bioc45$bio8/10
bioc45$bio9<- bioc45$bio9/10
bioc45$bio10<- bioc45$bio10/10
bioc45$bio11<- bioc45$bio11/10

bioc45

plot(bioc45[[1]])
plot(bioc45[[12]])

# write raster
raster::writeRaster(bioc45, "data/bioclim/chelsa_bioc45_studyarea.tif", format = "GTiff",
                    overwrite = T)
#------------------------------------------------------------------------------------#

# read study area with chelsa bioclim variables
bioc45 <- raster::brick("data/bioclim/chelsa_bioc45_studyarea.tif")
bioc45

names(bioc45) # check layer names
names(bioc45) <- names(bio) # rename
names(bioc45) 

# exclude same multicolinear var as in bioc
bioc45 <- exclude(bioc45, vs)
bioc45

# normalize data
bioc45.norm <- raster::scale(bioc45)
bioc45.norm

raster::writeRaster(bioc45.norm, "data/bioclim/bioc45_normalized.tif", format = "GTiff",
                    overwrite = T) # write raster
#----------------------------------------------------------------------------------------#

#---2.4.2.1 Prediction and averaging----

# load normalized data
bioc45.norm <- raster::brick("data/bioclim/bioc45_normalized.tif")

names(bioc45.norm) # check layer names
names(bioc45.norm) <- names(bioc) # rename
names(bioc45.norm)
bioc45.norm

# generate predictions using all 25 fitted models (p45 = prediction rcp 45)
p45 <- list()
for(k in 1:25){
  p45[[k]] <- predict(bioc45.norm, rf_down[[k]], type = "prob", index = 2)
}
p45

# make raster brick
p45.brick <- raster::brick(p45)

raster::writeRaster(p45.brick, "output/predict/p45_rf-down.tif", format = "GTiff",
                    overwrite = T)

# calculate mean over all raster layers
p45.avg <- calc(p45.brick, fun = mean)
p45.avg

raster::writeRaster(p45.avg, "output/ensemble/p45_avg_rf-down.tif", format = "GTiff",
                    overwrite = T) # write raster

#---2.4.2.2 Plotting----

# load raster data
p45.avg <- raster("output/ensemble/p45_avg_rf-down.tif")

# make df
p45.pts <- rasterToPoints(p45.avg)
p45.df <- data.frame(p45.pts)
colnames(p45.df) <- c("lon", "lat", "niche")
head(p45.df)

# enm on sa map
label_45 <- read.csv2("data/spatial/label_45.csv")

f4 <- ggplot() +
  geom_raster(data = srtm_df, aes(lon, lat, fill = elev)) +
  scale_fill_gradient(low = "grey90", high = "grey10") +
  geom_sf(data = world, col = "black", fill = NA, size = 0.1) +
  guides(fill = "none") +
  new_scale_fill() +
  geom_raster(data = p45.df, aes(lon, lat, fill = niche)) +
  scale_fill_viridis_c(option = "mako") +
  geom_path(data = rivers.df, aes(x = long, y = lat, group = group), 
            color = "dodgerblue", size = 0.3) +
  geom_text(data = label_45, aes(lon, lat, label = period),
            color = "white", size = 2.75, fontface = "bold") +
coord_sf(xlim = c(66, 106.5), 
         ylim = c(20.8, 39), 
         expand = FALSE) +
  theme(legend.key.size = unit(0.3, 'cm'), 
        legend.key.height = unit(0.3, 'cm'), 
        legend.key.width = unit(0.3, 'cm'), 
        legend.title = element_text(size = 6, face = "bold"), 
        legend.text = element_text(size = 6),
        panel.background =  element_rect(fill = "aliceblue"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title = element_blank(),
        axis.ticks = element_line(size = .3)) +
  labs(x = "", y = "", fill = "Habitat\nsuitability")
f4

ggsave("p45_map.pdf", plot = f4,      # save as pdf
       path = "plots/sdm/publ", device = "pdf",
       width = 16, height = 12 , units = c("cm"), dpi = 300)
#------------------------------------------#

#---2.4.3 RCP 60----

## 2061-2080 / CMIP5 / CCSM4 (CHELSA v1.2, res 30 arc sec)

# load CHELSA bioclim data 
bio601 <- raster("D:/CHELSA/V1.2/future/CCSM4/CHELSA_bio_mon_CCSM4_rcp60_r1i1p1_g025.nc_1_2061-2080_V1.2.tif")
bio602 <- raster("D:/CHELSA/V1.2/future/CCSM4/CHELSA_bio_mon_CCSM4_rcp60_r1i1p1_g025.nc_2_2061-2080_V1.2.tif")
bio603 <- raster("D:/CHELSA/V1.2/future/CCSM4/CHELSA_bio_mon_CCSM4_rcp60_r1i1p1_g025.nc_3_2061-2080_V1.2.tif")
bio604 <- raster("D:/CHELSA/V1.2/future/CCSM4/CHELSA_bio_mon_CCSM4_rcp60_r1i1p1_g025.nc_4_2061-2080_V1.2.tif")
bio605 <- raster("D:/CHELSA/V1.2/future/CCSM4/CHELSA_bio_mon_CCSM4_rcp60_r1i1p1_g025.nc_5_2061-2080_V1.2.tif")
bio606 <- raster("D:/CHELSA/V1.2/future/CCSM4/CHELSA_bio_mon_CCSM4_rcp60_r1i1p1_g025.nc_6_2061-2080_V1.2.tif")
bio607 <- raster("D:/CHELSA/V1.2/future/CCSM4/CHELSA_bio_mon_CCSM4_rcp60_r1i1p1_g025.nc_7_2061-2080_V1.2.tif")
bio608 <- raster("D:/CHELSA/V1.2/future/CCSM4/CHELSA_bio_mon_CCSM4_rcp60_r1i1p1_g025.nc_8_2061-2080_V1.2.tif")
bio609 <- raster("D:/CHELSA/V1.2/future/CCSM4/CHELSA_bio_mon_CCSM4_rcp60_r1i1p1_g025.nc_9_2061-2080_V1.2.tif")
bio6010 <- raster("D:/CHELSA/V1.2/future/CCSM4/CHELSA_bio_mon_CCSM4_rcp60_r1i1p1_g025.nc_10_2061-2080_V1.2.tif")
bio6011 <- raster("D:/CHELSA/V1.2/future/CCSM4/CHELSA_bio_mon_CCSM4_rcp60_r1i1p1_g025.nc_11_2061-2080_V1.2.tif")
bio6012 <- raster("D:/CHELSA/V1.2/future/CCSM4/CHELSA_bio_mon_CCSM4_rcp60_r1i1p1_g025.nc_12_2061-2080_V1.2.tif")
bio6013 <- raster("D:/CHELSA/V1.2/future/CCSM4/CHELSA_bio_mon_CCSM4_rcp60_r1i1p1_g025.nc_13_2061-2080_V1.2.tif")
bio6014 <- raster("D:/CHELSA/V1.2/future/CCSM4/CHELSA_bio_mon_CCSM4_rcp60_r1i1p1_g025.nc_14_2061-2080_V1.2.tif")
bio6015 <- raster("D:/CHELSA/V1.2/future/CCSM4/CHELSA_bio_mon_CCSM4_rcp60_r1i1p1_g025.nc_15_2061-2080_V1.2.tif")
bio6016 <- raster("D:/CHELSA/V1.2/future/CCSM4/CHELSA_bio_mon_CCSM4_rcp60_r1i1p1_g025.nc_16_2061-2080_V1.2.tif")
bio6017 <- raster("D:/CHELSA/V1.2/future/CCSM4/CHELSA_bio_mon_CCSM4_rcp60_r1i1p1_g025.nc_17_2061-2080_V1.2.tif")
bio6018 <- raster("D:/CHELSA/V1.2/future/CCSM4/CHELSA_bio_mon_CCSM4_rcp60_r1i1p1_g025.nc_18_2061-2080_V1.2.tif")
bio6019 <- raster("D:/CHELSA/V1.2/future/CCSM4/CHELSA_bio_mon_CCSM4_rcp60_r1i1p1_g025.nc_19_2061-2080_V1.2.tif")

# stack layers
bio60 <- stack(bio601, bio602, bio603, bio604, bio605, bio606,
               bio607, bio608, bio609, bio6010, bio6011, bio6012,
               bio6013, bio6014, bio6015, bio6016, bio6017, bio6018, bio6019)

# change names of biof so the model can understand
names(bio60)
names(bio60) <- names(bio)
names(bio60)
#------------------------------------------------#

# define extent of study area
studyarea <- shapefile("data/spatial/polygon.shp")
studyarea

# set projection
projection(bio60) 
projection(bio60) <- projection(bioc)

# crop study area
e <- extent(studyarea)
bioc60 <- raster::crop(bio60, e) # crop raster object
bioc60 <- raster::mask(bioc60, studyarea) # mask to remove NA values
bioc60

plot(bioc60[[1]])
plot(bioc60[[12]])
#-------------------------------------------------------#

# harmonize layers (see CHELSA v1.2 tech specification)
bioc60$bio1<- bioc60$bio1/10
bioc60$bio2<- bioc60$bio2/10
bioc60$bio3<- bioc60$bio3/10
bioc60$bio4<- bioc60$bio4/10
bioc60$bio5<- bioc60$bio5/10
bioc60$bio6<- bioc60$bio6/10
bioc60$bio7<- bioc60$bio7/10
bioc60$bio8<- bioc60$bio8/10
bioc60$bio9<- bioc60$bio9/10
bioc60$bio10<- bioc60$bio10/10
bioc60$bio11<- bioc60$bio11/10

bioc60

plot(bioc60[[1]])
plot(bioc60[[12]])

# write raster
raster::writeRaster(bioc60, "data/bioclim/chelsa_bioc60_studyarea.tif", format = "GTiff",
                    overwrite = T)
#------------------------------------------------------------------------------------#

# read study area with chelsa bioclim variables
bioc60 <- raster::brick("data/bioclim/chelsa_bioc60_studyarea.tif")
bioc60

names(bioc60) # check layer names
names(bioc60) <- names(bio) # rename
names(bioc60)

# exclude same multicolinear var as in bioc
bioc60 <- exclude(bioc60, vs)
bioc60

# normalize data
bioc60.norm <- raster::scale(bioc60)
bioc60.norm

raster::writeRaster(bioc60.norm, "data/bioclim/bioc60_normalized.tif", format = "GTiff",
                    overwrite = T) # write raster
#------------------------------------------#

#---2.4.3.1 Prediction and averaging----

# load normalized data
bioc60.norm <- raster::brick("data/bioclim/bioc60_normalized.tif")

names(bioc60.norm) # check layer names
names(bioc60.norm) <- names(bioc) # rename
names(bioc60.norm)
bioc60.norm

# generate predictions using all 25 fitted models (p26 = prediction rcp 26)
p60 <- list()
for(k in 1:25){
  p60[[k]] <- predict(bioc60.norm, rf_down[[k]], type = "prob", index = 2)
}
p60

# make raster brick
p60.brick <- raster::brick(p60)

raster::writeRaster(p60.brick, "output/predict/p60_rf-down.tif", format = "GTiff",
                    overwrite = T)

# calculate mean over all raster layers
p60.avg <- calc(p60.brick, fun = mean)
p60.avg

raster::writeRaster(p60.avg, "output/ensemble/p60_avg_rf-down.tif", format = "GTiff",
                    overwrite = T) # write raster

#---2.4.3.2 Plotting----

# load raster data
p60.avg <- raster("output/ensemble/p60_avg_rf-down.tif")

# make df
p60.pts <- rasterToPoints(p60.avg)
p60.df <- data.frame(p60.pts)
colnames(p60.df) <- c("lon", "lat", "niche")
head(p60.df)

# enm on sa map
label_60 <- read.csv2("data/spatial/label_60.csv")

f5 <- ggplot() +
  geom_raster(data = srtm_df, aes(lon, lat, fill = elev)) +
  scale_fill_gradient(low = "grey90", high = "grey10") +
  geom_sf(data = world, col = "black", fill = NA, size = 0.1) +
  guides(fill = "none") +
  new_scale_fill() +
  geom_raster(data = p60.df, aes(lon, lat, fill = niche)) +
  scale_fill_viridis_c(option = "mako") +
  geom_path(data = rivers.df, aes(x = long, y = lat, group = group), 
            color = "dodgerblue", size = 0.3) +
  geom_text(data = riv_lab, aes(lon, lat, label = river),
            color = "white", size = 2, fontface = "bold", angle = -45) +
  geom_text(data = label_60, aes(lon, lat, label = period),
            color = "black", size = 5, fontface = "bold") +
  annotation_scale(pad_x = unit(8.7, "cm"), 
                   pad_y = unit(6.2, "cm"),
                   height = unit(0.1, "cm"),
                   text_cex = 0.55) +
  annotation_north_arrow(which_north = "true", 
                         height = unit(.6, "cm"),
                         width = unit(.5, "cm"),
                         pad_x = unit(0.1, "cm"), 
                         pad_y = unit(0.1, "cm"),
                         style = north_arrow_fancy_orienteering) +
  coord_sf(xlim = c(66, 106.5), 
           ylim = c(20.8, 39), 
           expand = FALSE) +
  theme(legend.key.size = unit(0.4, 'cm'), 
        legend.key.height = unit(0.4, 'cm'), 
        legend.key.width = unit(0.4, 'cm'), 
        legend.title = element_text(size = 7.5, face = "bold"), 
        legend.text = element_text(size = 7.25),
        panel.background =  element_rect(fill = "aliceblue"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 8)) +
  labs(x = "Longitude", y = "Latitude", fill = "Habitat\nsuitability") 
f5

ggsave("fig_s5.png", plot = f5,      # save as png
       path = "plots/suppl", device = "png",
       width = 16, height = 8 , units = c("cm"), dpi = 330)
#---------------------------------------------------------#

#---2.4.4 RCP 85----

## 2061-2080 / CMIP5 / CCSM4 (CHELSA v1.2, res 30 arc sec)

# load CHELSA bioclim data 
bio851 <- raster("D:/CHELSA/V1.2/future/CCSM4/CHELSA_bio_mon_CCSM4_rcp85_r1i1p1_g025.nc_1_2061-2080_V1.2.tif")
bio852 <- raster("D:/CHELSA/V1.2/future/CCSM4/CHELSA_bio_mon_CCSM4_rcp85_r1i1p1_g025.nc_2_2061-2080_V1.2.tif")
bio853 <- raster("D:/CHELSA/V1.2/future/CCSM4/CHELSA_bio_mon_CCSM4_rcp85_r1i1p1_g025.nc_3_2061-2080_V1.2.tif")
bio854 <- raster("D:/CHELSA/V1.2/future/CCSM4/CHELSA_bio_mon_CCSM4_rcp85_r1i1p1_g025.nc_4_2061-2080_V1.2.tif")
bio855 <- raster("D:/CHELSA/V1.2/future/CCSM4/CHELSA_bio_mon_CCSM4_rcp85_r1i1p1_g025.nc_5_2061-2080_V1.2.tif")
bio856 <- raster("D:/CHELSA/V1.2/future/CCSM4/CHELSA_bio_mon_CCSM4_rcp85_r1i1p1_g025.nc_6_2061-2080_V1.2.tif")
bio857 <- raster("D:/CHELSA/V1.2/future/CCSM4/CHELSA_bio_mon_CCSM4_rcp85_r1i1p1_g025.nc_7_2061-2080_V1.2.tif")
bio858 <- raster("D:/CHELSA/V1.2/future/CCSM4/CHELSA_bio_mon_CCSM4_rcp85_r1i1p1_g025.nc_8_2061-2080_V1.2.tif")
bio859 <- raster("D:/CHELSA/V1.2/future/CCSM4/CHELSA_bio_mon_CCSM4_rcp85_r1i1p1_g025.nc_9_2061-2080_V1.2.tif")
bio8510 <- raster("D:/CHELSA/V1.2/future/CCSM4/CHELSA_bio_mon_CCSM4_rcp85_r1i1p1_g025.nc_10_2061-2080_V1.2.tif")
bio8511 <- raster("D:/CHELSA/V1.2/future/CCSM4/CHELSA_bio_mon_CCSM4_rcp85_r1i1p1_g025.nc_11_2061-2080_V1.2.tif")
bio8512 <- raster("D:/CHELSA/V1.2/future/CCSM4/CHELSA_bio_mon_CCSM4_rcp85_r1i1p1_g025.nc_12_2061-2080_V1.2.tif")
bio8513 <- raster("D:/CHELSA/V1.2/future/CCSM4/CHELSA_bio_mon_CCSM4_rcp85_r1i1p1_g025.nc_13_2061-2080_V1.2.tif")
bio8514 <- raster("D:/CHELSA/V1.2/future/CCSM4/CHELSA_bio_mon_CCSM4_rcp85_r1i1p1_g025.nc_14_2061-2080_V1.2.tif")
bio8515 <- raster("D:/CHELSA/V1.2/future/CCSM4/CHELSA_bio_mon_CCSM4_rcp85_r1i1p1_g025.nc_15_2061-2080_V1.2.tif")
bio8516 <- raster("D:/CHELSA/V1.2/future/CCSM4/CHELSA_bio_mon_CCSM4_rcp85_r1i1p1_g025.nc_16_2061-2080_V1.2.tif")
bio8517 <- raster("D:/CHELSA/V1.2/future/CCSM4/CHELSA_bio_mon_CCSM4_rcp85_r1i1p1_g025.nc_17_2061-2080_V1.2.tif")
bio8518 <- raster("D:/CHELSA/V1.2/future/CCSM4/CHELSA_bio_mon_CCSM4_rcp85_r1i1p1_g025.nc_18_2061-2080_V1.2.tif")
bio8519 <- raster("D:/CHELSA/V1.2/future/CCSM4/CHELSA_bio_mon_CCSM4_rcp85_r1i1p1_g025.nc_19_2061-2080_V1.2.tif")

# stack layers
bio85 <- stack(bio851, bio852, bio853, bio854, bio855, bio856,
               bio857, bio858, bio859, bio8510, bio8511, bio8512,
               bio8513, bio8514, bio8515, bio8516, bio8517, bio8518, bio8519)

# change names of biof so the model can understand
names(bio85)
names(bio85) <- names(bio)
names(bio85)
#------------------------------------------------#

# define extent of study area
studyarea <- shapefile("data/spatial/polygon.shp")
studyarea

# set projection
projection(bio85) 
projection(bio85) <- projection(bioc)

# crop study area
e <- extent(studyarea)
bioc85 <- raster::crop(bio85, e) # crop raster object
bioc85 <- raster::mask(bioc85, studyarea) # mask to remove NA values
bioc85

plot(bioc85[[1]])
plot(bioc85[[12]])
#-------------------------------------------------------#

# harmonize layers (see CHELSA v1.2 tech specification)
bioc85$bio1<- bioc85$bio1/10
bioc85$bio2<- bioc85$bio2/10
bioc85$bio3<- bioc85$bio3/10
bioc85$bio4<- bioc85$bio4/10
bioc85$bio5<- bioc85$bio5/10
bioc85$bio6<- bioc85$bio6/10
bioc85$bio7<- bioc85$bio7/10
bioc85$bio8<- bioc85$bio8/10
bioc85$bio9<- bioc85$bio9/10
bioc85$bio10<- bioc85$bio10/10
bioc85$bio11<- bioc85$bio11/10

bioc85

plot(bioc85[[1]])
plot(bioc85[[12]])

# write raster
raster::writeRaster(bioc85, "data/bioclim/chelsa_bioc85_studyarea.tif", format = "GTiff",
                    overwrite = T)
#------------------------------------------------------------------------------------#

# read study area with chelsa bioclim variables
bioc85 <- raster::brick("data/bioclim/chelsa_bioc85_studyarea.tif")
bioc85

names(bioc85) # check layer names
names(bioc85) <- names(bio) # rename
names(bioc85)

# exclude same multicolinear var as in bioc
bioc85 <- exclude(bioc85, vs)
bioc85

# normalize data
bioc85.norm <- raster::scale(bioc85)
bioc85.norm

raster::writeRaster(bioc85.norm, "data/bioclim/bioc85_normalized.tif", format = "GTiff",
                    overwrite = T) # write raster
#------------------------------------------#

#---2.4.4.1 Prediction and averaging----

# load normalized data
bioc85.norm <- raster::brick("data/bioclim/bioc85_normalized.tif")

names(bioc85.norm) # check layer names
names(bioc85.norm) <- names(bioc) # rename
names(bioc85.norm)
bioc85.norm

# generate predictions using all 25 fitted models (p26 = prediction rcp 26)
p85 <- list()
for(k in 1:25){
  p85[[k]] <- predict(bioc85.norm, rf_down[[k]], type = "prob", index = 2)
}
p85

# make raster brick
p85.brick <- raster::brick(p85)

raster::writeRaster(p85.brick, "output/predict/p85_rf-down.tif", format = "GTiff",
                    overwrite = T)

# calculate mean over all raster layers
p85.avg <- calc(p85.brick, fun = mean)
p85.avg

raster::writeRaster(p85.avg, "output/ensemble/p85_avg_rf-down.tif", format = "GTiff",
                    overwrite = T) # write raster

#---2.4.4.2 Plotting----

# load raster data
p85.avg <- raster("output/ensemble/p85_avg_rf-down.tif")

# make df
p85.pts <- rasterToPoints(p85.avg)
p85.df <- data.frame(p85.pts)
colnames(p85.df) <- c("lon", "lat", "niche")
head(p85.df)

# enm on sa map
label_85 <- read.csv2("data/spatial/label_85.csv")

f6 <- ggplot() +
  geom_raster(data = srtm_df, aes(lon, lat, fill = elev)) +
  scale_fill_gradient(low = "grey90", high = "grey10") +
  geom_sf(data = world, col = "black", fill = NA, size = 0.1) +
  guides(fill = "none") +
  new_scale_fill() +
  geom_raster(data = p85.df, aes(lon, lat, fill = niche)) +
  scale_fill_viridis_c(option = "mako") +
  geom_path(data = rivers.df, aes(x = long, y = lat, group = group), 
            color = "dodgerblue", size = 0.3) +
  geom_text(data = riv_lab, aes(lon, lat, label = river),
            color = "white", size = 2, fontface = "bold", angle = -45) +
  geom_text(data = label_85, aes(lon, lat, label = period),
            color = "black", size = 5, fontface = "bold") +
  annotation_scale(pad_x = unit(8.7, "cm"), 
                   pad_y = unit(6.2, "cm"),
                   height = unit(0.1, "cm"),
                   text_cex = 0.55) +
  annotation_north_arrow(which_north = "true", 
                         height = unit(.6, "cm"),
                         width = unit(.5, "cm"),
                         pad_x = unit(0.1, "cm"), 
                         pad_y = unit(0.1, "cm"),
                         style = north_arrow_fancy_orienteering) +
  coord_sf(xlim = c(66, 106.5), 
           ylim = c(20.8, 39), 
           expand = FALSE) +
  theme(legend.key.size = unit(0.4, 'cm'), 
        legend.key.height = unit(0.4, 'cm'), 
        legend.key.width = unit(0.4, 'cm'), 
        legend.title = element_text(size = 7.5, face = "bold"), 
        legend.text = element_text(size = 7.25),
        panel.background =  element_rect(fill = "aliceblue"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 8)) +
  labs(x = "Longitude", y = "Latitude", fill = "Habitat\nsuitability") 
f6

ggsave("fig_s6.png", plot = f6,      # save as png
       path = "plots/suppl", device = "png",
       width = 16, height = 8 , units = c("cm"), dpi = 330)
#---------------------------------------------------------#


#---------------------------#
#---2.5 Model evaluation----

#---2.5.1 ROC AUC----

## ROC AUC / test data
# relative likelihood predicted on test data
pred_down.test <- list()
for(k in 1:25){
  pred_down.test[[k]] <- predict(rf_down[[k]], ddata.test[[k]], type = "prob")[, "1"]
}

head(pred_down.test[[1]])
head(pred_down.test[[25]])

# get AUC
auc.test <- vector()
for(k in 1:25){
  auc.test[[k]] <- as.vector(auc(roc(ddata.test[[k]]$species, pred_down.test[[k]])))
}

auc.test
mean(auc.test)
sd(auc.test)

# plot ROC curves
roc.test <- list()
for(k in 1:25){
  roc.test[[k]] <- roc(ddata.test[[k]]$species, pred_down.test[[k]])
}

roc.test.plot <- ggroc(roc.test, legacy.axes = TRUE) +
  theme_classic() +
  geom_abline(size = .75, col = "darkgrey") +
  geom_line(size = .5) +
  xlab("1-Specificity (false positives)") +
  ylab("Sensitivity (true positives)") +
  theme(legend.position = "none",
        axis.title = element_text(size = 9))
roc.test.plot

## ROC AUC / training data
# relative likelihood predicted on train data
pred_down.train <- list()
for(k in 1:25){
  pred_down.train[[k]] <- predict(rf_down[[k]], ddata.train[[k]], type = "prob")[, "1"]
}

head(pred_down.train[[1]])
head(pred_down.train[[25]])

# get AUC for training data
auc.train <- vector()
for(k in 1:25){
  auc.train[[k]] <- as.vector(auc(roc(ddata.train[[k]]$species, pred_down.train[[k]])))
}

auc.train
mean(auc.train)
sd(auc.train)

# plot ROC curves
roc.train <- list()
for(k in 1:25){
  roc.train[[k]] <- roc(ddata.train[[k]]$species, pred_down.train[[k]])
}

roc.train.plot <- ggroc(roc.train, legacy.axes = TRUE) + 
  theme_classic() +
  geom_abline(size = .75, col = "darkgrey") +
  geom_line(size = .5) +
  xlab("1-Specificity (false positives)") +
  ylab("Sensitivity (true positives)") +
  theme(legend.position = "none",
        axis.title = element_text(size = 9))
roc.train.plot

## graphical display of both roc and auc of test and train data
# plot AUC boxplot of train and test data
auc.df <- data.frame(auc.test, auc.train) # make df
auc.df$test <- c(rep("test", 25)) # add column with character
auc.df$train <- c(rep("train", 25)) # add column with character

auc.test.plot <- ggplot(data = auc.df, mapping = aes(x = test, y = auc.test)) +
  geom_jitter(alpha = .9, width = .15, size = 1.5, color = "darkgrey") +
  geom_boxplot(alpha = 0.1) +
  theme_minimal() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title = element_text(size = 9)) +
  xlab("Test data (30%)") +
  ylab("AUC")
auc.test.plot # test data

auc.train.plot <- ggplot(data = auc.df, mapping = aes(x = train, y = auc.train)) +
  geom_jitter(alpha = .9, width = .15, size = 1.5, color = "darkgrey") +
  geom_boxplot(alpha = 0.1) +
  theme_minimal() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title = element_text(size = 9)) +
  xlab("Training data (70%)") +
  ylab("AUC")
auc.train.plot # train data

# make panel plot combining roc curves and auc boxplots
roc.auc.plot <- plot_grid(roc.train.plot, auc.train.plot,
                          roc.test.plot, auc.test.plot,
                          ncol = 2,
                          labels = c("a", "", "b", ""),
                          align = "v")
roc.auc.plot

ggsave("roc.auc.rf-down.pdf", plot = roc.auc.plot,      # save as pdf
       path = "plots/evaluation/", device = "pdf",
       width = 16, height = 12 , units = c("cm"), dpi = 600)
#--------------------------------------------------------#

#---2.5.2 TSS----

## TSS / test data
max.tss.test <- list()
for(k in 1:25){
  max.tss.test[[k]] <- ecospat.max.tss(pred_down.test[[k]],ddata.test[[k]]$species)
}

max.tss.test

## TSS / training data
max.tss.train <- list()
for(k in 1:25){
  max.tss.train[[k]] <- ecospat.max.tss(pred_down.train[[k]],ddata.train[[k]]$species)
}

max.tss.train

## make overview table
# vectors with max tss values and their respective thresholds for train and test data
tss.test <- vector()
th.test <- vector()
tss.train <- vector()
th.train <- vector()
for(k in 1:25){
  tss.test[[k]] <- as.vector(max.tss.test[[k]]$max.TSS)
  th.test[[k]] <- as.vector(max.tss.test[[k]]$max.threshold)
  tss.train[[k]] <- as.vector(max.tss.train[[k]]$max.TSS)
  th.train[[k]] <- as.vector(max.tss.train[[k]]$max.threshold)
}

tss.th.df <- data.frame(tss.train, th.train, tss.test, th.test) # make df
View(tss.th.df)

mean(tss.test) # descriptive tss
sd(tss.test)
mean(tss.train)
sd(tss.train)

mean(th.test) # descriptive th
sd(th.test)
mean(th.train)
sd(th.train)
#-------------------------------------------------------------------------#

#---2.5.3 Variable importance----

# calculate mean decrease of accuracy for all variables of fitted models
important <- list()
for(k in 1:25){
  important[[k]] <- as.data.frame(importance(rf_down[[k]], type = 1))
}
important
varImpPlot(rf_down[[20]], type = 1)

# merge mean decrease accuracy of vars of 25 models in one df
important.var <- rep(c("bio2", "bio3", "bio5", "bio8", "bio9", "bio12", "bio14", "bio15", "bio18"), 25) # colnames

important.mda <- vector()
important.mda <- as.vector(sapply(important, "[[", "MeanDecreaseAccuracy")) # get mdas to write as columns in new df

important.df <- data.frame(important.var, important.mda) # bind

var.imp.box <- ggplot(data = important.df, 
                      aes(x = reorder(important.var, important.mda), # plot
                          y = important.mda)) +
  geom_jitter(alpha = .9, width = .15, size = .75, color = "darkgrey") +
  geom_boxplot(alpha = 0.1, outlier.colour = "black") +
  theme_minimal() +
  theme(axis.title = element_text(size = 9)) +
  coord_flip() +
  xlab("Bioclimatic variable") +
  ylab("Mean decrease in accuracy")

var.imp.box # boxplot

ggsave("varImp.pdf", plot = var.imp.box,      # save as pdf
       path = "plots/evaluation/", device = "pdf",
       width = 16, height = 6, units = c("cm"), dpi = 600)
#--------------------------------------------------------#


#----------------------------------------#
#---2.6 Range shift in space and time----

# load raster layers with predictions
pp.avg   <- raster("output/ensemble/pp_avg_rf-down.tif")
plgm.avg <- raster("output/ensemble/plgm_avg_rf-down.tif")
p26.avg  <- raster("output/ensemble/p26_avg_rf-down.tif")
p45.avg  <- raster("output/ensemble/p45_avg_rf-down.tif")
p60.avg  <- raster("output/ensemble/p60_avg_rf-down.tif")
p85.avg  <- raster("output/ensemble/p85_avg_rf-down.tif")

#---2.6.1 Binarized range shift (threshold set to maximizing tss)----

# set threshold for max TSS (for train data)
th1 <- 0.63

# convert past, current and future prob of occurrences to TSS-thresholded presences and absences
pa1 <- raster(pp.avg) # create empty raster with same extent and resolution than en1
pa1[] <- ifelse(pp.avg[] >= th1, 1, 0) # which pixels within pp.avg >= than th, if yes 1, if no 0
plot(pa1) 

pa2 <- raster(p26.avg) # create empty raster with same extent and resolution than en2
pa2[] <- ifelse(p26.avg[] >= th1, 1, 0) # which pixels within p26.avg >= than th, if yes 1, if no 0
plot(pa2) 

pa3 <- raster(p45.avg) # create empty raster with same extent and resolution than en2
pa3[] <- ifelse(p45.avg[] >= th1, 1, 0) # which pixels within p45.avg >= than th, if yes 1, if no 0
plot(pa3)

pa4 <- raster(p60.avg) # create empty raster with same extent and resolution than en2
pa4[] <- ifelse(p60.avg[] >= th1, 1, 0) # which pixels within p60.avg >= than th, if yes 1, if no 0
plot(pa4)

pa5 <- raster(p85.avg) # create empty raster with same extent and resolution than en2
pa5[] <- ifelse(p85.avg[] >= th1, 1, 0) # which pixels within p85.avg >= than th, if yes 1, if no 0
plot(pa5)

pa6 <- raster(plgm.avg) # create empty raster with same extent and resolution than en1
pa6[] <- ifelse(plgm.avg[] >= th1, 1, 0) # which pixels within plgm.avg >= than th, if yes 1, if no 0
plot(pa6) 

#---2.6.1.1 LGM - RCP 4.5----

# stabilization from lgm to future 45
paletti <- colorRampPalette(c("chocolate", "yellow", "lightgreen"))
stabil <- pa6+pa1+pa3 # cells with value = 3 have always been stable (suitable)
plot(stabil, col = paletti(3))

# change from lgm to future
changelgm45b <- pa3 - pa6
plot(changelgm45b, col = paletti(3))

# gain and loss 
changex.pts <- rasterToPoints(changelgm45b)
changex.df <- data.frame(changex.pts)
colnames(changex.df) <- c("lon", "lat", "niche")
head(changex.df)

# extract gains
gain <- "1"
change.df <- filter(changex.df, niche %in% gain)

change.df <- mutate(change.df, 
                    niche = factor(case_when(niche %in% gain ~ "Gain")))
head(change.df)

# extract losses
loss <- "-1"
change1.df <- filter(changex.df, niche %in% loss)

change1.df <- mutate(change1.df, 
                    niche = factor(case_when(niche %in% loss ~ "Loss")))
head(change1.df)

# make df
stabil.pts <- rasterToPoints(stabil)
stabil.df <- data.frame(stabil.pts)
colnames(stabil.df) <- c("lon", "lat", "niche")
head(stabil.df)

# extract stable grid cells (3)
stab <- "3"
stabil.df <- filter(stabil.df, niche %in% stab)

stabil.df <- mutate(stabil.df, 
                    niche = factor(case_when(niche %in% stab ~ "Stabilization")))
head(stabil.df)

# combine dfs
changes.df <- rbind(change.df, stabil.df, change1.df)

# plot
change.f <- ggplot() +
  geom_raster(data = srtm_df, aes(lon, lat, fill = elev)) +
  scale_fill_gradient(low = "grey90", high = "grey10") +
  geom_sf(data = world, col = "black", fill = NA, size = 0.1) +
  guides(fill = "none") +
  new_scale_fill() +
  geom_raster(data = changes.df, aes(lon, lat, fill = niche)) +
  scale_fill_manual(values =  c("#ef6548", "#80cdc1", "#d73027")) +
  geom_path(data = rivers.df, aes(x = long, y = lat, group = group), 
            color = "dodgerblue", size = 0.3) +
  geom_text(data = riv_lab, aes(lon, lat, label = river),
            color = "white", size = 2, fontface = "bold", angle = -45) +
  annotation_scale(pad_x = unit(7.7, "cm"), 
                   pad_y = unit(5.5, "cm"),
                   height = unit(0.1, "cm"),
                   text_cex = 0.55) +
  annotation_north_arrow(which_north = "true", 
                         height = unit(.6, "cm"),
                         width = unit(.5, "cm"),
                         pad_x = unit(0.1, "cm"), 
                         pad_y = unit(0.1, "cm"),
                         style = north_arrow_fancy_orienteering) +
  coord_sf(xlim = c(66, 106.5), 
           ylim = c(20.8, 39), 
           expand = FALSE) +
  theme(legend.key.size = unit(0.4, 'cm'), 
        legend.key.height = unit(0.4, 'cm'), 
        legend.key.width = unit(0.4, 'cm'), 
        legend.title = element_text(size = 7.5, face = "bold"), 
        legend.text = element_text(size = 7.25),
        panel.background =  element_rect(fill = "aliceblue"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 8)) +
  labs(x = "Longitude", y = "Latitude", fill = "Habitat") 
change.f

ggsave("changes_bin.pdf", plot = change.f,      # save as pdf
       path = "plots/sdm/publ", device = "pdf",
       width = 16, height = 12 , units = c("cm"), dpi = 300)

#---2.6.1.2 LGM - PP----

# change from lgm to present (Fig 3c)
changelgmppb <- pa1 - pa6

# gain and loss 
changelgmpp.pts <- rasterToPoints(changelgmppb)
changelgmpp.df <- data.frame(changelgmpp.pts)
colnames(changelgmpp.df) <- c("lon", "lat", "niche")
head(changelgmpp.df)

# extract gains
gain <- "1"
lgmppG.df <- filter(changelgmpp.df, niche %in% gain)

lgmppG.df <- mutate(lgmppG.df, 
                    niche = factor(case_when(niche %in% gain ~ "Gain")))
head(lgmppG.df)

# extract losses
loss <- "-1"
lgmppL.df <- filter(changelgmpp.df, niche %in% loss)

lgmppL.df <- mutate(lgmppL.df, 
                     niche = factor(case_when(niche %in% loss ~ "Loss")))
head(lgmppL.df)

# combine dfs
lgmppGL.df <- rbind(lgmppG.df, lgmppL.df)

# plot
lgmpp.f <- ggplot() +
  geom_raster(data = srtm_df, aes(lon, lat, fill = elev)) +
  scale_fill_gradient(low = "grey90", high = "grey10") +
  geom_sf(data = world, col = "black", fill = NA, size = 0.1) +
  guides(fill = "none") +
  new_scale_fill() +
  geom_raster(data = lgmppGL.df, aes(lon, lat, fill = niche)) +
  scale_fill_manual(values =  c("#80cdc1", "#d73027")) +
  geom_path(data = rivers.df, aes(x = long, y = lat, group = group), 
            color = "dodgerblue", size = 0.3) +
  geom_text(data = riv_lab, aes(lon, lat, label = river),
            color = "white", size = 2, fontface = "bold", angle = -45) +
  annotation_scale(pad_x = unit(7.7, "cm"), 
                   pad_y = unit(5.5, "cm"),
                   height = unit(0.1, "cm"),
                   text_cex = 0.55) +
  annotation_north_arrow(which_north = "true", 
                         height = unit(.6, "cm"),
                         width = unit(.5, "cm"),
                         pad_x = unit(0.1, "cm"), 
                         pad_y = unit(0.1, "cm"),
                         style = north_arrow_fancy_orienteering) +
  coord_sf(xlim = c(66, 106.5), 
           ylim = c(20.8, 39), 
           expand = FALSE) +
  theme(legend.key.size = unit(0.4, 'cm'), 
        legend.key.height = unit(0.4, 'cm'), 
        legend.key.width = unit(0.4, 'cm'), 
        legend.title = element_text(size = 7.5, face = "bold"), 
        legend.text = element_text(size = 7.25),
        panel.background =  element_rect(fill = "aliceblue"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 8)) +
  labs(x = "Longitude", y = "Latitude", fill = "Habitat") 
lgmpp.f


#---2.6.1.3 PP - RCP 4.5----

# change from present to future 45 (Fig 3d)
changepp45b <- pa3 - pa1

# gain and loss 
changepp45.pts <- rasterToPoints(changepp45b)
changepp45.df <- data.frame(changepp45.pts)
colnames(changepp45.df) <- c("lon", "lat", "niche")
head(changepp45.df)

# extract gains
gain <- "1"
pp45G.df <- filter(changepp45.df, niche %in% gain)

pp45G.df <- mutate(pp45G.df, 
                    niche = factor(case_when(niche %in% gain ~ "Gain")))
head(pp45G.df)

# extract losses
loss <- "-1"
pp45L.df <- filter(changepp45.df, niche %in% loss)

pp45L.df <- mutate(pp45L.df, 
                    niche = factor(case_when(niche %in% loss ~ "Loss")))
head(pp45L.df)

# combine dfs
pp45GL.df <- rbind(pp45G.df, pp45L.df)

# plot
pp45.f <- ggplot() +
  geom_raster(data = srtm_df, aes(lon, lat, fill = elev)) +
  scale_fill_gradient(low = "grey90", high = "grey10") +
  geom_sf(data = world, col = "black", fill = NA, size = 0.1) +
  guides(fill = "none") +
  new_scale_fill() +
  geom_raster(data = pp45GL.df, aes(lon, lat, fill = niche)) +
  scale_fill_manual(values =  c("#80cdc1", "#d73027")) +
  geom_path(data = rivers.df, aes(x = long, y = lat, group = group), 
            color = "dodgerblue", size = 0.3) +
  geom_text(data = riv_lab, aes(lon, lat, label = river),
            color = "white", size = 2, fontface = "bold", angle = -45) +
  annotation_scale(pad_x = unit(7.7, "cm"), 
                   pad_y = unit(5.5, "cm"),
                   height = unit(0.1, "cm"),
                   text_cex = 0.55) +
  annotation_north_arrow(which_north = "true", 
                         height = unit(.6, "cm"),
                         width = unit(.5, "cm"),
                         pad_x = unit(0.1, "cm"), 
                         pad_y = unit(0.1, "cm"),
                         style = north_arrow_fancy_orienteering) +
  coord_sf(xlim = c(66, 106.5), 
           ylim = c(20.8, 39), 
           expand = FALSE) +
  theme(legend.key.size = unit(0.4, 'cm'), 
        legend.key.height = unit(0.4, 'cm'), 
        legend.key.width = unit(0.4, 'cm'), 
        legend.title = element_text(size = 7.5, face = "bold"), 
        legend.text = element_text(size = 7.25),
        panel.background =  element_rect(fill = "aliceblue"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 8)) +
  labs(x = "Longitude", y = "Latitude", fill = "Habitat") 
pp45.f

#---2.6.1.4 PP - RCP 2.6----

# change from present to future 26 (Fig S5b)
changepp26b <- pa2 - pa1

# gain and loss 
changepp26.pts <- rasterToPoints(changepp26b)
changepp26.df <- data.frame(changepp26.pts)
colnames(changepp26.df) <- c("lon", "lat", "niche")
head(changepp26.df)

# extract gains
gain <- "1"
pp26G.df <- filter(changepp26.df, niche %in% gain)

pp26G.df <- mutate(pp26G.df, 
                   niche = factor(case_when(niche %in% gain ~ "Gain")))
head(pp26G.df)

# extract losses
loss <- "-1"
pp26L.df <- filter(changepp26.df, niche %in% loss)

pp26L.df <- mutate(pp26L.df, 
                   niche = factor(case_when(niche %in% loss ~ "Loss")))
head(pp26L.df)

# combine dfs
pp26GL.df <- rbind(pp26G.df, pp26L.df)

# plot
pp26.f <- ggplot() +
  geom_raster(data = srtm_df, aes(lon, lat, fill = elev)) +
  scale_fill_gradient(low = "grey90", high = "grey10") +
  geom_sf(data = world, col = "black", fill = NA, size = 0.1) +
  guides(fill = "none") +
  new_scale_fill() +
  geom_raster(data = pp26GL.df, aes(lon, lat, fill = niche)) +
  scale_fill_manual(values =  c("#80cdc1", "#d73027")) +
  geom_path(data = rivers.df, aes(x = long, y = lat, group = group), 
            color = "dodgerblue", size = 0.3) +
  geom_text(data = riv_lab, aes(lon, lat, label = river),
            color = "white", size = 2, fontface = "bold", angle = -45) +
  annotation_scale(pad_x = unit(7.7, "cm"), 
                   pad_y = unit(5.5, "cm"),
                   height = unit(0.1, "cm"),
                   text_cex = 0.55) +
  annotation_north_arrow(which_north = "true", 
                         height = unit(.6, "cm"),
                         width = unit(.5, "cm"),
                         pad_x = unit(0.1, "cm"), 
                         pad_y = unit(0.1, "cm"),
                         style = north_arrow_fancy_orienteering) +
  coord_sf(xlim = c(66, 106.5), 
           ylim = c(20.8, 39), 
           expand = FALSE) +
  theme(legend.key.size = unit(0.4, 'cm'), 
        legend.key.height = unit(0.4, 'cm'), 
        legend.key.width = unit(0.4, 'cm'), 
        legend.title = element_text(size = 7.5, face = "bold"), 
        legend.text = element_text(size = 7.25),
        panel.background =  element_rect(fill = "aliceblue"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 8)) +
  labs(x = "Longitude", y = "Latitude", fill = "Habitat") 
pp26.f

#---2.6.1.5 PP - RCP 6.0----

# change from present to future 60 (Fig Fig S6b)
changepp60b <- pa4 - pa1

# gain and loss 
changepp60.pts <- rasterToPoints(changepp60b)
changepp60.df <- data.frame(changepp60.pts)
colnames(changepp60.df) <- c("lon", "lat", "niche")
head(changepp60.df)

# extract gains
gain <- "1"
pp60G.df <- filter(changepp60.df, niche %in% gain)

pp60G.df <- mutate(pp60G.df, 
                   niche = factor(case_when(niche %in% gain ~ "Gain")))
head(pp60G.df)

# extract losses
loss <- "-1"
pp60L.df <- filter(changepp60.df, niche %in% loss)

pp60L.df <- mutate(pp60L.df, 
                   niche = factor(case_when(niche %in% loss ~ "Loss")))
head(pp60L.df)

# combine dfs
pp60GL.df <- rbind(pp60G.df, pp60L.df)

# plot
pp60.f <- ggplot() +
  geom_raster(data = srtm_df, aes(lon, lat, fill = elev)) +
  scale_fill_gradient(low = "grey90", high = "grey10") +
  geom_sf(data = world, col = "black", fill = NA, size = 0.1) +
  guides(fill = "none") +
  new_scale_fill() +
  geom_raster(data = pp60GL.df, aes(lon, lat, fill = niche)) +
  scale_fill_manual(values =  c("#80cdc1", "#d73027")) +
  geom_path(data = rivers.df, aes(x = long, y = lat, group = group), 
            color = "dodgerblue", size = 0.3) +
  geom_text(data = riv_lab, aes(lon, lat, label = river),
            color = "white", size = 2, fontface = "bold", angle = -45) +
  annotation_scale(pad_x = unit(7.7, "cm"), 
                   pad_y = unit(5.5, "cm"),
                   height = unit(0.1, "cm"),
                   text_cex = 0.55) +
  annotation_north_arrow(which_north = "true", 
                         height = unit(.6, "cm"),
                         width = unit(.5, "cm"),
                         pad_x = unit(0.1, "cm"), 
                         pad_y = unit(0.1, "cm"),
                         style = north_arrow_fancy_orienteering) +
  coord_sf(xlim = c(66, 106.5), 
           ylim = c(20.8, 39), 
           expand = FALSE) +
  theme(legend.key.size = unit(0.4, 'cm'), 
        legend.key.height = unit(0.4, 'cm'), 
        legend.key.width = unit(0.4, 'cm'), 
        legend.title = element_text(size = 7.5, face = "bold"), 
        legend.text = element_text(size = 7.25),
        panel.background =  element_rect(fill = "aliceblue"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 8)) +
  labs(x = "Longitude", y = "Latitude", fill = "Habitat") 
pp60.f

#---2.6.1.6 PP - RCP 8.5----

# change from present to future 85 (Fig Fig S7b)
changepp85b <- pa5 - pa1

# gain and loss 
changepp85.pts <- rasterToPoints(changepp85b)
changepp85.df <- data.frame(changepp85.pts)
colnames(changepp85.df) <- c("lon", "lat", "niche")
head(changepp85.df)

# extract gains
gain <- "1"
pp85G.df <- filter(changepp85.df, niche %in% gain)

pp85G.df <- mutate(pp85G.df, 
                   niche = factor(case_when(niche %in% gain ~ "Gain")))
head(pp85G.df)

# extract losses
loss <- "-1"
pp85L.df <- filter(changepp85.df, niche %in% loss)

pp85L.df <- mutate(pp85L.df, 
                   niche = factor(case_when(niche %in% loss ~ "Loss")))
head(pp85L.df)

# combine dfs
pp85GL.df <- rbind(pp85G.df, pp85L.df)

# plot
pp85.f <- ggplot() +
  geom_raster(data = srtm_df, aes(lon, lat, fill = elev)) +
  scale_fill_gradient(low = "grey90", high = "grey10") +
  geom_sf(data = world, col = "black", fill = NA, size = 0.1) +
  guides(fill = "none") +
  new_scale_fill() +
  geom_raster(data = pp85GL.df, aes(lon, lat, fill = niche)) +
  scale_fill_manual(values =  c("#80cdc1", "#d73027")) +
  geom_path(data = rivers.df, aes(x = long, y = lat, group = group), 
            color = "dodgerblue", size = 0.3) +
  geom_text(data = riv_lab, aes(lon, lat, label = river),
            color = "white", size = 2, fontface = "bold", angle = -45) +
  annotation_scale(pad_x = unit(7.7, "cm"), 
                   pad_y = unit(5.5, "cm"),
                   height = unit(0.1, "cm"),
                   text_cex = 0.55) +
  annotation_north_arrow(which_north = "true", 
                         height = unit(.6, "cm"),
                         width = unit(.5, "cm"),
                         pad_x = unit(0.1, "cm"), 
                         pad_y = unit(0.1, "cm"),
                         style = north_arrow_fancy_orienteering) +
  coord_sf(xlim = c(66, 106.5), 
           ylim = c(20.8, 39), 
           expand = FALSE) +
  theme(legend.key.size = unit(0.4, 'cm'), 
        legend.key.height = unit(0.4, 'cm'), 
        legend.key.width = unit(0.4, 'cm'), 
        legend.title = element_text(size = 7.5, face = "bold"), 
        legend.text = element_text(size = 7.25),
        panel.background =  element_rect(fill = "aliceblue"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 8)) +
  labs(x = "Longitude", y = "Latitude", fill = "Habitat") 
pp85.f
#------------------------------------------------------------------------------------------#

#------------------------------------#
#---2.7 Visualize bio12 over time----

## make bio12 dfs
# lgm
lgmc.pts <- rasterToPoints(lgmc$bio12)
lgmc.df <- data.frame(lgmc.pts)
colnames(lgmc.df) <- c("lon", "lat", "bio12")
head(lgmc.df)

# present
bioc.pts <- rasterToPoints(bioc$bio12)
bioc.df <- data.frame(bioc.pts)
colnames(bioc.df) <- c("lon", "lat", "bio12")
head(bioc.df)

# future 26
bioc26.pts <- rasterToPoints(bioc26$bio12)
bioc26.df <- data.frame(bioc26.pts)
colnames(bioc26.df) <- c("lon", "lat", "bio12")
head(bioc26.df)

# future 45
bioc45.pts <- rasterToPoints(bioc45$bio12)
bioc45.df <- data.frame(bioc45.pts)
colnames(bioc45.df) <- c("lon", "lat", "bio12")
head(bioc45.df)

# future 60
bioc60.pts <- rasterToPoints(bioc60$bio12)
bioc60.df <- data.frame(bioc60.pts)
colnames(bioc60.df) <- c("lon", "lat", "bio12")
head(bioc60.df)

# future 85
bioc85.pts <- rasterToPoints(bioc85$bio12)
bioc85.df <- data.frame(bioc85.pts)
colnames(bioc85.df) <- c("lon", "lat", "bio12")
head(bioc85.df)
#--------------------------------------------------#

## preparation
# color palette
rain <-  colorRampPalette(c("#e0f3f8", "#2c7bb6", "#fee090", "#fdae61", "#f46d43", "#d73027", "#b2182b", "#a50026", "#67001f"))

# define polygon to crop gap area
x_coords_zoom <- c(76.8, 76.8, 83.3, 83.3, 76.8)
y_coords_zoom <- c(28.5, 31.3, 31.3, 28.5, 28.5)

polyZ <- sp::Polygon(cbind(x_coords_zoom,y_coords_zoom)) # make polygon from matrix

polyZOO <- sp::Polygons(list(polyZ), ID = "A") # polygon class
str(polyZOO, 1)

polyZOOM <- sp::SpatialPolygons(list(polyZOO)) # spatial polygon
polyZOOM

shapefile(x = polyZOOM, 
          file = "data/spatial/polyZOOM.shp",
          overwrite = T) # save as shapefile

# define extent of zoom
zoom <- shapefile("data/spatial/polyZOOM.shp")
zoom

# set projection
projection(zoom) 
projection(zoom) <- projection(bioc)

# crop study area
e_zoom <- extent(zoom)
zoom_sa <- raster::crop(bioc$bio12, e_zoom) # crop raster object
zoom_sa <- raster::mask(zoom_sa, zoom) # mask to remove NA values
zoom_sa

# present zoom
zoom.pts <- rasterToPoints(zoom_sa)
zoom.df <- data.frame(zoom.pts)
colnames(zoom.df) <- c("lon", "lat", "bio12")
head(zoom.df)

# crop sp spatial dataframe
sp_zoom <- raster::crop(sp, e_zoom)
sp_zoom <- fortify(sp_zoom)
#-------------------------------------------------------------------#

## plots
# lgm
label_lgm <- read.csv2("data/spatial/label_lgm.csv")

lgmc_12 <- ggplot() +
  geom_raster(data = srtm_df, aes(lon, lat, fill = elev)) +
  scale_fill_gradient(low = "grey90", high = "grey10") +
  geom_sf(data = world, col = "black", fill = NA, size = 0.1) +
  guides(fill = "none") +
  new_scale_fill() +
  geom_raster(data = lgmc.df, aes(lon, lat, fill = bio12)) +
  scale_fill_gradientn(colors = rain(100)) +
  guides(fill = "none") +
  geom_path(data = rivers.df, aes(x = long, y = lat, group = group), 
            color = "dodgerblue4", size = 0.3) +
  geom_text(data = label_lgm, aes(lon, lat, label = period),
            color = "white", size = 2.75, fontface = "bold") +
  coord_sf(xlim = c(66, 106.5), 
           ylim = c(20.8, 39), 
           expand = FALSE) +
  theme(legend.key.size = unit(0.4, 'cm'), 
        legend.key.height = unit(0.4, 'cm'), 
        legend.key.width = unit(0.4, 'cm'), 
        legend.title = element_text(size = 7.5, face = "bold"), 
        legend.text = element_text(size = 7.25),
        panel.background =  element_rect(fill = "aliceblue"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title = element_blank(),
        axis.ticks = element_line(size = .3)) +
  labs(x = "", y = "", fill = "bio12") +
  expand_limits(fill = c(NA, 9999))
lgmc_12

ggsave("lgmc_12_map.pdf", plot = lgmc_12,      # save as pdf
       path = "plots/bio12/publ/", device = "pdf",
       width = 16, height = 12 , units = c("cm"), dpi = 300)

# present
rect <- read.csv2("data/spatial/rect_zoom.csv") # for zoom rectangle 
sp <- read.csv2("data/occ/hastatus_all.csv")
arrow <- data.frame(lon1 = 76.8, lat1 = 28.5, lon2 = 76.8, lat2 = 21)

bioc_12 <- ggplot() +
  geom_raster(data = srtm_df, aes(lon, lat, fill = elev)) +
  scale_fill_gradient(low = "grey90", high = "grey10") +
  geom_sf(data = world, col = "black", fill = NA, size = 0.1) +
  guides(fill = "none") +
  new_scale_fill() +
  geom_raster(data = bioc.df, aes(lon, lat, fill = bio12)) +
  scale_fill_gradientn(colors = rain(100)) +
  geom_path(data = rivers.df, aes(x = long, y = lat, group = group), 
            color = "dodgerblue4", size = 0.3) +
  geom_text(data = riv_lab, aes(lon, lat, label = river),
            color = "black", size = 2, fontface = "bold", angle = -45) +
  geom_path(data = rect, aes(x = lon, y = lat), 
           color = "black", size = .4) +
  geom_segment(aes(x = lon1, y = lat1, xend = lon2, yend = lat2), data = arrow, 
               arrow = arrow(length = unit(0.25, "cm")), 
               size = .4, color = "black", lineend = "round") +
  geom_point(data = sp, aes(x = decimalLongitude, 
                          y = decimalLatitude),
         size = 1, shape = 21,
         color = "white", fill = "black") +
  annotation_scale(pad_x = unit(7.6, "cm"), 
                   pad_y = unit(5.4, "cm"),
                   height = unit(0.1, "cm"),
                   text_cex = 0.55) +
  annotation_north_arrow(which_north = "true", 
                         height = unit(.6, "cm"),
                         width = unit(.5, "cm"),
                         pad_x = unit(0.1, "cm"), 
                         pad_y = unit(0.1, "cm"),
                         style = north_arrow_fancy_orienteering) +
  coord_sf(xlim = c(66, 106.5), 
           ylim = c(20.8, 39), 
           expand = FALSE) +
  theme(legend.key.size = unit(0.4, 'cm'), 
        legend.key.height = unit(0.4, 'cm'), 
        legend.key.width = unit(0.4, 'cm'), 
        legend.title = element_text(size = 7.5, face = "bold"), 
        legend.text = element_text(size = 7.25),
        panel.background =  element_rect(fill = "aliceblue"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 8)) +
  labs(x = "Longitude", y = "Latitude", fill = "bio12") +
  expand_limits(fill = c(NA, 9999))
bioc_12

ggsave("bioc_12_map.pdf", plot = bioc_12,      # save as pdf
       path = "plots/bio12/publ/", device = "pdf",
       width = 16, height = 12 , units = c("cm"), dpi = 300)

# present zoom 
zoom_12 <- ggplot() +
  geom_raster(data = srtm_df, aes(lon, lat, fill = elev)) +
  scale_fill_gradient(low = "grey90", high = "grey10") +
  geom_sf(data = world, col = "black", fill = NA, size = 0.1) +
  guides(fill = "none") +
  new_scale_fill() +
  geom_raster(data = zoom.df, aes(lon, lat, fill = bio12)) +
  scale_fill_gradientn(colors = rain(100)) +
  guides(fill = "none") +
  geom_point(data = sp_zoom, aes(x = long, 
                            y = lat),
             size = 1, shape = 21,
             color = "white", fill = "black") +
  annotation_scale(pad_x = unit(0.1, "cm"), 
                   pad_y = unit(0.1, "cm"),
                   height = unit(0.1, "cm"),
                   text_cex = 0.45) +
  coord_sf(xlim = c(76.8, 83.3), 
           ylim = c(28.5, 31.3), 
           expand = FALSE) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title = element_blank(),
        axis.ticks = element_line(size = .3)) +
  labs(x = "", y = "", fill = "bio12") +
  expand_limits(fill = c(NA, 9999))
zoom_12

ggsave("zoom_12_map.pdf", plot = zoom_12,      # save as pdf
       path = "plots/bio12/publ/", device = "pdf",
       width = 8, height = 6 , units = c("cm"), dpi = 300)

# future 26
label_26 <- read.csv2("data/spatial/label_26.csv")

bioc26_12 <- ggplot() +
  geom_raster(data = srtm_df, aes(lon, lat, fill = elev)) +
  scale_fill_gradient(low = "grey90", high = "grey10") +
  geom_sf(data = world, col = "black", fill = NA, size = 0.1) +
  guides(fill = "none") +
  new_scale_fill() +
  geom_raster(data = bioc26.df, aes(lon, lat, fill = bio12)) +
  scale_fill_gradientn(colors = rain(100)) +
  geom_path(data = rivers.df, aes(x = long, y = lat, group = group), 
            color = "dodgerblue4", size = 0.4) +
  geom_text(data = riv_lab, aes(lon, lat, label = river),
            color = "black", size = 2, fontface = "bold", angle = -45) +
  geom_text(data = label_26, aes(lon, lat, label = period),
            color = "black", size = 5, fontface = "bold") +
  annotation_scale(pad_x = unit(8.75, "cm"), 
                   pad_y = unit(6.25, "cm"),
                   height = unit(0.1, "cm"),
                   text_cex = 0.55) +
  annotation_north_arrow(which_north = "true", 
                         height = unit(.6, "cm"),
                         width = unit(.5, "cm"),
                         pad_x = unit(0.1, "cm"), 
                         pad_y = unit(0.1, "cm"),
                         style = north_arrow_fancy_orienteering) +
  coord_sf(xlim = c(66, 106.5), 
           ylim = c(20.8, 39), 
           expand = FALSE) +
  theme(legend.key.size = unit(0.4, 'cm'), 
        legend.key.height = unit(0.4, 'cm'), 
        legend.key.width = unit(0.4, 'cm'), 
        legend.title = element_text(size = 7.5, face = "bold"), 
        legend.text = element_text(size = 7.25),
        panel.background =  element_rect(fill = "aliceblue"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 8)) +
  labs(x = "Longitude", y = "Latitude", fill = "bio12") +
  expand_limits(fill = c(NA, 9999))
bioc26_12

ggsave("fig_s7.png", plot = bioc26_12,      # save as pdf
       path = "plots/suppl/", device = "png",
       width = 16, height = 8 , units = c("cm"), dpi = 330)

# future 45
label_45 <- read.csv2("data/spatial/label_45.csv")

bioc45_12 <- ggplot() +
  geom_raster(data = srtm_df, aes(lon, lat, fill = elev)) +
  scale_fill_gradient(low = "grey90", high = "grey10") +
  geom_sf(data = world, col = "black", fill = NA, size = 0.1) +
  guides(fill = "none") +
  new_scale_fill() +
  geom_raster(data = bioc45.df, aes(lon, lat, fill = bio12)) +
  scale_fill_gradientn(colors = rain(100)) +
  guides(fill = "none") +
  geom_path(data = rivers.df, aes(x = long, y = lat, group = group), 
            color = "dodgerblue4", size = 0.3) +
  geom_text(data = label_45, aes(lon, lat, label = period),
            color = "white", size = 2.75, fontface = "bold") +
  coord_sf(xlim = c(66, 106.5), 
           ylim = c(20.8, 39), 
           expand = FALSE) +
  theme(legend.key.size = unit(0.4, 'cm'), 
        legend.key.height = unit(0.4, 'cm'), 
        legend.key.width = unit(0.4, 'cm'), 
        legend.title = element_text(size = 7.5, face = "bold"), 
        legend.text = element_text(size = 7.25),
        panel.background =  element_rect(fill = "aliceblue"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        axis.title = element_blank(),
        axis.ticks = element_line(size = .3)) +
  labs(x = "", y = "", fill = "bio12") +
  expand_limits(fill = c(NA, 9999))
bioc45_12

ggsave("bioc45_12_map.pdf", plot = bioc45_12,      # save as pdf
       path = "plots/bio12/publ/", device = "pdf",
       width = 16, height = 12 , units = c("cm"), dpi = 300)

# future 60
label_60 <- read.csv2("data/spatial/label_60.csv")

bioc60_12 <- ggplot() +
  geom_raster(data = srtm_df, aes(lon, lat, fill = elev)) +
  scale_fill_gradient(low = "grey90", high = "grey10") +
  geom_sf(data = world, col = "black", fill = NA, size = 0.1) +
  guides(fill = "none") +
  new_scale_fill() +
  geom_raster(data = bioc60.df, aes(lon, lat, fill = bio12)) +
  scale_fill_gradientn(colors = rain(100)) +
  geom_path(data = rivers.df, aes(x = long, y = lat, group = group), 
            color = "dodgerblue4", size = 0.4) +
  geom_text(data = riv_lab, aes(lon, lat, label = river),
            color = "black", size = 2, fontface = "bold", angle = -45) +
  geom_text(data = label_60, aes(lon, lat, label = period),
            color = "black", size = 5, fontface = "bold") +
  annotation_scale(pad_x = unit(8.75, "cm"), 
                   pad_y = unit(6.25, "cm"),
                   height = unit(0.1, "cm"),
                   text_cex = 0.55) +
  annotation_north_arrow(which_north = "true", 
                         height = unit(.6, "cm"),
                         width = unit(.5, "cm"),
                         pad_x = unit(0.1, "cm"), 
                         pad_y = unit(0.1, "cm"),
                         style = north_arrow_fancy_orienteering) +
  coord_sf(xlim = c(66, 106.5), 
           ylim = c(20.8, 39), 
           expand = FALSE) +
  theme(legend.key.size = unit(0.4, 'cm'), 
        legend.key.height = unit(0.4, 'cm'), 
        legend.key.width = unit(0.4, 'cm'), 
        legend.title = element_text(size = 7.5, face = "bold"), 
        legend.text = element_text(size = 7.25),
        panel.background =  element_rect(fill = "aliceblue"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 8)) +
  labs(x = "Longitude", y = "Latitude", fill = "bio12") +
  expand_limits(fill = c(NA, 9999))
bioc60_12

ggsave("fig_s8.png", plot = bioc60_12,      # save as pdf
       path = "plots/suppl/", device = "png",
       width = 16, height = 8 , units = c("cm"), dpi = 330)

# future 85
label_85 <- read.csv2("data/spatial/label_85.csv")

bioc85_12 <- ggplot() +
  geom_raster(data = srtm_df, aes(lon, lat, fill = elev)) +
  scale_fill_gradient(low = "grey90", high = "grey10") +
  geom_sf(data = world, col = "black", fill = NA, size = 0.1) +
  guides(fill = "none") +
  new_scale_fill() +
  geom_raster(data = bioc85.df, aes(lon, lat, fill = bio12)) +
  scale_fill_gradientn(colors = rain(100)) +
  geom_path(data = rivers.df, aes(x = long, y = lat, group = group), 
            color = "dodgerblue4", size = 0.4) +
  geom_text(data = riv_lab, aes(lon, lat, label = river),
            color = "black", size = 2, fontface = "bold", angle = -45) +
  geom_text(data = label_85, aes(lon, lat, label = period),
            color = "black", size = 5, fontface = "bold") +
  annotation_scale(pad_x = unit(8.75, "cm"), 
                   pad_y = unit(6.25, "cm"),
                   height = unit(0.1, "cm"),
                   text_cex = 0.55) +
  annotation_north_arrow(which_north = "true", 
                         height = unit(.6, "cm"),
                         width = unit(.5, "cm"),
                         pad_x = unit(0.1, "cm"), 
                         pad_y = unit(0.1, "cm"),
                         style = north_arrow_fancy_orienteering) +
  coord_sf(xlim = c(66, 106.5), 
           ylim = c(20.8, 39), 
           expand = FALSE) +
  theme(legend.key.size = unit(0.4, 'cm'), 
        legend.key.height = unit(0.4, 'cm'), 
        legend.key.width = unit(0.4, 'cm'), 
        legend.title = element_text(size = 7.5, face = "bold"), 
        legend.text = element_text(size = 7.25),
        panel.background =  element_rect(fill = "aliceblue"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 8)) +
  labs(x = "Longitude", y = "Latitude", fill = "bio12") +
  expand_limits(fill = c(NA, 9999))
bioc85_12

ggsave("fig_s9.png", plot = bioc85_12,      # save as pdf
       path = "plots/suppl/", device = "png",
       width = 16, height = 8 , units = c("cm"), dpi = 330)
#--------------------------------------------------------------------------#


#-------------------------------------------#
#---3 Statistical analysis-------------------

#---3.1 PCA ------------------------

# prepare (normalized) bioclim var for actual occurrences
occ_pca <- sdmData(species ~ ., train = sp, predictors = bioc.norm) # bioc var of duplicate-cleaned sp occ
occ_pca

occ_pca <- as.data.frame(occ_pca) # make it a df
head(occ_pca)
occ_pca <- occ_pca[, 3:11] # drop first two columns
head(occ_pca)

occ_pca$Group <- c("occ") # add column with ID "occ" for occurrence
head(occ_pca)

write.csv2(occ_pca, "data/pca/occ_pca_bioc.norm.csv", row.names = F)
#------------------------------------------#

# prepare bioclim var for pseudo-occurrences within the gap

# elev data from srtm
srtm_gap1 <- raster::raster("data/srtm/srtm_54_07.tif") 
srtm_gap2 <- raster::raster("data/srtm/srtm_55_07.tif")
srtm_gap3 <- raster::raster("data/srtm/srtm_56_07.tif")

srtm_gap <- raster::mosaic(srtm_gap1, srtm_gap2, srtm_gap3,  fun = mean) # merge srtm layers

# define polygon to crop gap area
x_coords_gap <- c(85, 85, 98, 98, 85)
y_coords_gap <- c(26, 30, 30, 26, 26)

polyG <- sp::Polygon(cbind(x_coords_gap,y_coords_gap)) # make polygon from matrix

polyGA <- sp::Polygons(list(polyG), ID = "A") # polygon class
str(polyGA, 1)

polyGAP <- sp::SpatialPolygons(list(polyGA)) # spatial polygon
polyGAP

plot(srtm_gap)
plot(sp, add = T)
plot(polyGAP, add = T)

shapefile(x = polyGAP, 
          file = "data/spatial/polyGAP.shp",
          overwrite = T) # save as shapefile
#-------------------------------------------#

# define extent of study area
gaparea <- shapefile("data/spatial/polyGAP.shp")
gaparea

# set projection
projection(gaparea) 
projection(gaparea) <- projection(bioc)

# crop study area
e_gap <- extent(gaparea)
gap <- raster::crop(srtm_gap, e_gap) # crop raster object
gap <- raster::mask(gap, gaparea) # mask to remove NA values
gap
plot(gap)
#--------------------------------------------------------------#

# select random points in gap area and extract elevation
set.seed(555)
pts_gap <- randomPoints(gap, 500)
points(pts_gap,col = "blue")

gap_df <- as.data.frame(pts_gap) # make df with coordinates of pts

pts_elev <- raster::extract(gap, pts_gap) # extract elev of pts
pts_elev <- as.data.frame(pts_elev) # make df

d_gap <- cbind(gap_df, pts_elev) # merge dfs with coord and respective elev of pts 
head(d_gap)

# prepare df as input for sdmData object
names(d_gap)[1] <- "lon"
names(d_gap)[2] <- "lat"
names(d_gap)[3] <- "elev"

d_gap$species <- 1 # 1 for presence (presence-only data)

d_gap$Group <- c("gap") # add column with ID "occ" for occurrence

d_gap <- d_gap %>%
  filter(elev >= 600 ) %>%
  filter(elev <= 3200) # select only points with elev range of rumex hastatus
head(d_gap)

range(d_gap$elev) # 603 - 3157
mean(d_gap$elev) # 1743.75
sd(d_gap$elev) # 753.63

write.csv2(d_gap, "data/pca/d_gap_elev.csv", row.names = F) # write csv file

# plot selected points within gap region
coordinates(d_gap) <- c("lon", "lat") # convert df to a spatial points df

plot(gap)
plot(d_gap, add = T, col = "blue")
#-----------------------------------------------------------------------#

# read pseudo-occurrences within gap and make spatial df
d_gap <- read.csv2("data/pca/d_gap_elev.csv", header = T)
coordinates(d_gap) <- c("lon", "lat")

# create sdmData object by extracting bioc variables from pseudo-occ
gap_pca <- sdmData(species ~ ., train = d_gap, predictors = bioc.norm)
gap_pca

gap_pca <- as.data.frame(gap_pca) # make it a df
head(gap_pca)
gap_pca <- gap_pca[, 4:13] # drop first three columns
head(gap_pca)

write.csv2(gap_pca, "data/pca/gap_pca_bioc.norm.csv", row.names = F) # write csv file
#------------------------------------------------------------------------------------------#

# read dfs with bioc.norm data
occ_pca_norm <- read.csv2("data/pca/occ_pca_bioc.norm.csv", header = T)
gap_pca_norm <- read.csv2("data/pca/gap_pca_bioc.norm.csv", header = T)

# merge dfs vertically
pca_og_norm <- rbind(occ_pca_norm, gap_pca_norm)
colnames(pca_og_norm) <- c("2", "3", " ", "8", "  ", "12", "14", "15", "18", "Group")
head(pca_og_norm)

# write csv file
#write.csv2(pca_og_norm, "data/pca/occ+gap_pca_bioc.csv", row.names = F) # write csv file

# perfom pca and plot biplot using ggfortify
pca_data_norm <- pca_og_norm[, -10]
pca <- prcomp(pca_data_norm, scale. = F, center = T) # data already scaled
summary(pca)
pca$rotation

# biplot
pca.plot <- autoplot(pca, size = 1, alpha = .8,
            loadings = TRUE, loadings.label = T,
            data = pca_og_norm, colour = "Group",
            loadings.colour = "black",
            loadings.label.size = 2.75,
            loadings.label.font = "bold",
            loadings.label.colour = "#e41a1c",
            loadings.label.vjust = 0, 
            loadings.label.hjust = -0.1,
            scale = F) +
            theme_minimal() +
            theme(legend.key.size = unit(.5, "cm"),
                  legend.title = element_blank(),
            text = element_text(size = 8.5),
            axis.ticks.x = element_blank(),
            axis.title = element_text(size = 7.5)) +
            scale_color_manual(values = c("black", "#80cdc1"))
pca.plot

# save biplot
ggsave("pca_plot.pdf", plot = pca.plot,      # save as pdf
       path = "plots/pca_lda/publ", device = "pdf",
       width = 8, height = 6, units = c("cm"), dpi = 600)
#--------------------------------------------#

#-------------------------------#
#---3.2 MANOVA & LDA------------

## CODE from Numerical Ecology with R (Borcard & al. 2018) + custom function (David Zélény)

## Unscaled data
# read dfs with unscaled bioc data
occ_pca_unscaled <- read.csv2("data/pca/occ_pca_bioc.csv", header = T)
gap_pca_unscaled <- read.csv2("data/pca/gap_pca_bioc.csv", header = T)

# merge dfs vertically
pca_og_unscaled <- rbind(occ_pca_unscaled, gap_pca_unscaled)

# select groups
gr <- pca_og_unscaled[, 10]
gr

# select environmental variables
envi_unscaled <- as.matrix(pca_og_unscaled[1:258, 1:9])

# calculate euclidean distance matrix
envi.d_unscaled <- dist(envi_unscaled)

#  analysis of multivariate homogeneity of group dispersions (variances)
envidat_unscaled <- betadisper(envi.d_unscaled, gr)
envidat_unscaled

# permutation-based test of multivariate homogeneity of group dispersions (variances)
permutest(envidat_unscaled) # significant, within-group covariance matrices are not homogeneous

## Normalized data
# read dfs with bioc.norm data
occ_pca_norm <- read.csv2("data/pca/occ_pca_bioc.norm.csv", header = T)
gap_pca_norm <- read.csv2("data/pca/gap_pca_bioc.norm.csv", header = T)

# merge dfs vertically
pca_og_norm <- rbind(occ_pca_norm, gap_pca_norm)

# select environmental variables
envi_norm <- as.matrix(pca_og_norm[1:258, 1:9])

# calculate euclidean distance matrix
envi.d_norm <- dist(envi_norm)

#  analysis of multivariate homogeneity of group dispersions (variances)
envidat_norm <- betadisper(envi.d_norm, gr)
envidat_norm

# permutation-based test of multivariate homogeneity of group dispersions (variances)
permutest(envidat_norm) # significant, within-group covariance matrices are not homogeneous
#------------------------------------------------------------------------------------------#

# perform manova (wilks test)
Wilks.test(envi_unscaled, gr) # with no assumed equal variances (like in Li et al. 2019, J. Biogeogr.)
Wilks.test(envi_norm, gr) 

# make envia df
envi_unscaled.df <- data.frame(envi_unscaled)
envi_norm.df <- data.frame(envi_norm) # make df

# compute lda
lda1 <- lda(gr ~ ., data = envi_unscaled.df)
lda1
lda2 <- lda(gr ~ ., data = envi_norm.df)
lda2

# predict and a posteriori groupings and post. probs.
lda1.values <- predict(lda1, envi_unscaled.df)
lda1.values$posterior
lda1.values$class
lda1.values$x

plot(lda1)

# make df for ggplot
occ_lda <- data.frame(lda1.values$x[1:130])
gap_lda <- data.frame(lda1.values$x[131:258])

colnames(occ_lda) <- "LD1"
colnames(gap_lda) <- "LD1"

occ_lda$Group <- "occ"
gap_lda$Group <- "gap"

lda.df <- rbind(occ_lda, gap_lda)

# plot lda histogram
dashline <- ddply(lda.df, "Group", summarise, grp.mean=mean(LD1))

ld1 <- ggplot(lda.df, aes(LD1 , fill = Group)) + 
       geom_density(alpha = .8) +
       geom_vline(data = dashline, 
              aes(xintercept = grp.mean),
                linetype = "dashed", size = .4) +
       xlab("LD1") +
       ylab("Density") +
       theme_minimal() +
       theme(legend.key.size = unit(.5, "cm"),
             legend.title = element_blank(),
       text = element_text(size = 8.5),
       axis.ticks.x=element_blank(),
       axis.title = element_text(size = 7.5)) +
       scale_fill_manual(values = c("black", "#80cdc1"))
ld1
       
ggsave("lda.pdf", plot = ld1,      # save as pdf
       path = "plots/pca_lda/publ", device = "pdf",
       width = 8, height = 6, units = c("cm"), dpi = 600)

# make panel plot combining roc curves and auc boxplots
pca.lda.plot <-  plot_grid(pca.plot, ld1,
                          ncol = 2,
                          labels = c("a", "b"),
                          align = "h")
pca.lda.plot

ggsave("pca.lda.panel.pdf", plot = pca.lda.plot,      # save as pdf
       path = "plots/pca_lda/publ", device = "pdf",
       width = 16, height = 6 , units = c("cm"), dpi = 600)
#-------------------------------------------------------------------#

#------------------------------------#
#---3.3 bio12 whisker plot------------
box12 <- data.frame(pca_og_unscaled[, c(6,10)])
head(box12)

#plot
plot12 <- ggplot(data = box12, aes(x = Group, y = bio12, color = Group)) +
  geom_jitter(alpha = .9, width = .25, size = 1) +
  geom_boxplot(alpha = .1, color = "black", size = .4) +
  scale_color_manual(values = c("#d73027", "#2c7bb6")) +
  guides(color = "none") +
  theme_minimal() +
  theme(axis.ticks.x = element_blank(),
        legend.key.size = unit(.5, "cm"),
        legend.title = element_blank(),
        axis.text.y = element_text(size = 6),
        axis.text.x = element_text(size = 7),                         
        axis.title = element_text(size = 7)) +
  labs(x = "", y = "bio12", color = "")
plot12 # test data
#-----------------------------------------------------#


#---------------------------#
#---3.4 Wilcoxon------------

# nonparametrical test for differing tendencies of bioclim variables
wilcox.test(envi_unscaled.df$bio2 ~ gr)
wilcox.test(envi_unscaled.df$bio3 ~ gr)
wilcox.test(envi_unscaled.df$bio5 ~ gr)
wilcox.test(envi_unscaled.df$bio8 ~ gr)
wilcox.test(envi_unscaled.df$bio9 ~ gr)
wilcox.test(envi_unscaled.df$bio12 ~ gr)
wilcox.test(envi_unscaled.df$bio14 ~ gr)
wilcox.test(envi_unscaled.df$bio15 ~ gr)
wilcox.test(envi_unscaled.df$bio18 ~ gr)
#--------------------------------------#


#--------------------------------#
#---4 Figure assembly-------------

# Figure 1 --> study_area_map

# Figure 2
fig2 <- f1 / (pca.plot | ld1) +
  plot_layout(heights = unit(c(6, 8), c("cm", "null"))) +
  plot_annotation(tag_levels = "a") & 
  theme(plot.tag = element_text(size = 12, face = "bold"))

ggsave("fig2.pdf", plot = fig2,      # save as pdf
       path = "plots/figs", device = "pdf",
       width = 16, height = 14, units = c("cm"), dpi = 300)

ggsave("fig2.png", plot = fig2,      # save as png
       path = "plots/figs", device = "png",
       width = 16, height = 14 , units = c("cm"), dpi = 330)

# Figure 3
layout3 <- "
AABB
CCCC
"

fig3 <- f2 + f4 + change.f +
  plot_layout(design = layout3, heights = c(1, 2)) +
  plot_annotation(tag_levels = "a") & 
  theme(plot.tag = element_text(size = 12, face = "bold"))

ggsave("fig3.pdf", plot = fig3,      # save as pdf
       path = "plots/figs", device = "pdf",
       width = 16, height = 12, units = c("cm"), dpi = 300)

ggsave("fig3.png", plot = fig3,      # save as png
       path = "plots/figs", device = "png",
       width = 16, height = 12 , units = c("cm"), dpi = 330)

# Figure 4
layout4 <- "
AABB
CCCC
DDEE
"

fig4 <- lgmc_12 + bioc45_12 + bioc_12 + zoom_12 + plot12 +
  plot_layout(design = layout4, heights = c(1, 2, 1)) +
  plot_annotation(tag_levels = "a") & 
  theme(plot.tag = element_text(size = 12, face = "bold"))

ggsave("fig4.pdf", plot = fig4,      # save as pdf
       path = "plots/figs", device = "pdf",
       width = 16, height = 16, units = c("cm"), dpi = 300)

ggsave("fig4.png", plot = fig4,      # save as png
       path = "plots/figs", device = "png",
       width = 16, height = 16 , units = c("cm"), dpi = 330)

# Figure 3 new
layout3 <- "
AABB
CCCC
DDDD
"

fig3new <- f2 + f4 + lgmpp.f + pp45.f +
  plot_layout(design = layout3, heights = c(1, 2, 2)) +
  plot_annotation(tag_levels = "a") & 
  theme(plot.tag = element_text(size = 12, face = "bold"))

ggsave("fig3new.pdf", plot = fig3new,      # save as pdf
       path = "plots/figs", device = "pdf",
       width = 16, height = 20, units = c("cm"), dpi = 300)

ggsave("fig3new.png", plot = fig3new,      # save as png
       path = "plots/figs", device = "png",
       width = 16, height = 20 , units = c("cm"), dpi = 330)
#--------------------------------------------------------------------------#
####---STOP---####