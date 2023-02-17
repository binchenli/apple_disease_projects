library(sdm)
#installAll()
library(dismo)
library(dplyr)
library(tidyr)
library(mapview)
library(raster)
library(parallel)
library(usdm)
# library(ggplot2)
library(maptools)
library(shiny)
library(kernlab)
library(sf)
library(rgdal)
library(rgeos)
library(sp)
library(corrplot)
library(ncdf4)
library(mapview)
library(lattice)

getwd()
setwd('../apple_disease_projects/datas/environment/')

#----------read data---------
#species

diseases_pest <- read.csv("diseases_8pests.csv",header = TRUE, check.names = F)

#spatial points data frame

coordinates(diseases_pest) ~ lon + lat
# proj4string(diseases_pest$lon, diseases_pest$lat) <- projection(raster())
mapview(bioc[[1]])
mapview(diseases_pest)

#shp
country <- readOGR("D:/datas/border/chian-country-provinces/2.shp")
crs(country) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

plot(country)
# country <- spTransform(country, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
# plot(country)

getwd()
#human footprint
# rastlist <- list.files(path = "D:/datas/hfp2000_2018/", pattern='.tif', all.files=TRUE, full.names=T)
humanfootprint<- stack("D:/datas/hfp2000_2018/humanfootprint_2000_2017_mean.tif")

# newproj<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
# humanfootprint <- projectRaster(humanfootprint,crs = newproj)

# humanfootprint_2000_2017_mean<- stackApply(humanfootprint , indices = 1, fun = mean)
# 
# rf <- writeRaster(humanfootprint_2000_2017_mean, 
#                   filename="D:/datas/hfp2000_2018/humanfootprint_2000_2017_mean.tif",overwrite=TRUE)
plot(mask(humanfootprint,country))
plot(humanfootprint )

# humanfootprint_2000_2017_mean <- mask(humanfootprint_2000_2017_mean,country)

humanfootprint <- resample(humanfootprint,clim,method= 'bilinear')


# bioc <- dropLayer(bioc,5)

#environment
rastlist <- list.files(path = "D:/datas/2020_2017_mean/", pattern='.nc', all.files=TRUE, full.names=T)
bioc <- stack(rastlist)
bioc <- resample(bioc, clim, method='bilinear')
# 
# #elev
# elev <- raster::getData('alt', country='CHN', mask=FALSE)
# #重采样
# elev <- resample(elev, bioc)
# palette<-colorRampPalette(c("blue","white","red"))
# plot(mask(elev,country),col=palette(10))
# # str(elev)


#-------------------------------------

bioc_hfp <- stack(bioc,mask(humanfootprint,country))
# bioc <- stack(bioc,mask(elev,country))
names(bioc_hfp) <- c('pre','tmn','tmp','tmx','hfp')
bioc_hfp
# rf <- writeRaster(bioc_hfp, filename="D:/datas/bioc_hfp.tif", format="GTiff", overwrite=TRUE)


#----------------soil and topographic-------------
fileslist <- list.files(path = "./soils/", pattern='.tif', all.files=TRUE, full.names=T)
soil <- stack(fileslist)

#重采样
soil <- resample(soil, bioc_hfp, method = 'bilinear')

topographic_files_list <- list.files('./Topographics/', pattern='.tif', all.files=TRUE, full.names=T)
topographic_files <- stack(topographic_files_list)

#重采样
topographic_files <- resample(topographic_files, bioc_hfp, method = 'bilinear')


envs <- stack(bioc_hfp, soil, topographic_files)
names(envs) <- c('pre','tmn','tmp','tmx','hfp','soc','prop','sand','silt','aspect','slope','curvature','dem')
envs
rf <- writeRaster(envs, filename="./envs.tif", format="GTiff", overwrite=TRUE)



