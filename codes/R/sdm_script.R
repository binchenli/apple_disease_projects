#---------- Script for ENM 2020 online course:

# In this demonstration, we are going to demonstrate the sdm (and usdm) packagge 
# to fit species distribution models for "Acinonyx jubatus" 
# and predict/project the potential distribution in the current and future times
# and measure the changes in the potential distribution (range shift) due to climate change


# install.packages('sdm')
# or
# devtools::install_github('babaknaimi/sdm')

# setwd

library(sdm)
# installAll() # only firest time after installing the sdm package


library(dismo)
library(dplyr)
library(tidyr)
library(mapview)

# Acinonyx jubatus
gbif("Acinonyx","jubatus",download= F)

sp <- gbif("Acinonyx","jubatus",download= T)

class(sp)

dim(sp)
table(sp$basisOfRecord)

sp <- sp %>% 
  filter(basisOfRecord %in% c("HUMAN_OBSERVATION","OBSERVATION","PRESERVED_SPECIMEN"))
nrow(sp)

spg <- sp %>% select(lon,lat)
head(spg)
spg$species <- 1
spg <- spg %>% drop_na()
nrow(spg)
#------------
class(spg)

coordinates(spg) <- c('lon','lat')
##################

# download the bioclim data:

bio <- raster::getData('worldclim', var='bio', res=10)
bio
names(bio)

plot(bio[[1]])  
points(spg)

e <- drawExtent()

spg <- crop(spg, e)

points(spg,col='red')
bioc <- crop(bio, e)
plot(bioc[[1]])

#------------
library(usdm)

vif(bioc)
ex <- raster::extract(bioc,spg)
head(ex)

v <- vifstep(ex)

#vifcor
v
bioc <- exclude(bioc, v)
#--------------------
library(sdm)


d <- sdmData(species~., spg, predictors= bioc, bg = list(method='gRandom',n=1000))
d

getmethodNames()

m <- sdm(species~., d, methods=c('glm','brt','rf','fda'), replication=c('sub','boot'),
         test.p=30,n=3, parallelSetting=list(ncore=4,method='parallel'))

m
#m@models$species$rf$`13`@object

gui(m)

p1 <- predict(m, bio,filename='pr.img')
p1
names(p1)

plot(p1[[c(1,7,13,23)]])

#en1 <- ensemble(m, bio, filename='en.img',setting=list(method='weighted',stat='tss',opt=2))
en1 <- ensemble(m, p1, filename='en.img',setting=list(method='weighted',stat='tss',opt=2))

plot(en1)
##################
biof <- raster::getData('CMIP5', var='bio', res=10, rcp=85, model='CN', year=70)

biof
plot(biof[[1]])
names(biof) <- names(bio)

en2 <- ensemble(m, biof, filename='enf.img',setting=list(method='weighted',stat='tss',opt=2))
#--------------
plot(en2)
cl <- colorRampPalette(c('#3E49BB','#3498DB','yellow','orange','red','darkred'))
#------
plot(en1, col=cl(200))
plot(en2, col=cl(200))

proj4string(spg) <- projection(en1)

library(mapview)
mapview(en1,col.regions=cl(200)) + spg

#-----------
ch <- en2 - en1
cl2 <- colorRampPalette(c('red','orange','yellow','gray','green','blue'))
plot(ch,col=cl2(200))
#----
df <- as.data.frame(d)
df <- data.frame(species=df$species,coordinates(d))
xy <- as.matrix(df[,c('lon','lat')])
head(xy)
p <- raster::extract(en1,xy)
head(p)
nrow(df)
length(p)
ev <- evaluates(df$species,p)
ev@statistics

th <- ev@threshold_based$threshold[2]

pa1 <- raster(en1)

pa1[] <- ifelse(en1[] >= th, 1, 0)
plot(pa1)

pa2 <- raster(en1)

pa2[] <- ifelse(en2[] >= th, 1, 0)
plot(pa2)

chp <- pa2 - pa1
plot(chp,col=c('red','gray','blue'))

#---------------

rcurve(m,id=7:12)

plot(getVarImp(m,method='rf'))

niche(bio,d,n=c('bio4','bio18'),col=cl(200))


getwd()





