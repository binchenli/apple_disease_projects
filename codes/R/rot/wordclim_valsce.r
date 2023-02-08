library(sdm)
library(dismo)
library(dplyr)
library(tidyr)
library(mapview)
library(raster)
library(parallel)
library(usdm)
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
#----------------------------
worldclim_bio <- raster::getData('worldclim', var='bio', res=10)
names(worldclim_bio)
plot(worldclim_bio[[1]]) 


country <- readOGR("../../../datas/environment/Province/省级行政区.shp")
country <- spTransform(country, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

#crop background to reserch area

worldclim_bioc <- crop(worldclim_bio, country)
plot(worldclim_bioc[[1]])

#check cor

vif_values <- vif(worldclim_bioc)
v <- vifstep(worldclim_bioc)
v
#Collinear variables were excluded

worldclim_bioc <- exclude(worldclim_bioc, v)
plot(worldclim_bioc)


#----------read data---------
disease_point<- read.csv("valsce.csv",header = TRUE, check.names = F)

#spatial points data frame

coordinates(disease_point) <- ~ lon + lat


#-----------------fit modes-----------------

worldclim_d <- sdmData(infestation_percentage
             ~., disease_point,
             predictors = worldclim_bioc,bg = list(method='gRandom',n=1000))
worldclim_d
worldclim_df <- as.data.frame(worldclim_d) # train datas
write.csv(worldclim_df,"./worldclim_data_valsce.csv", row.names = FALSE)

#getmethodNames()
worldclim_m <- sdm(infestation_percentage ~., 
         data=worldclim_d, methods=c('glmnet','maxlike','mlp','rbf',
                           'gam','gbm','rpart','svm','rf')
         ,replications='boot',test.p=30,n=10,
         parallelSetting=list(ncore=12,method='parallel')) #replications : sub , cv, boot
worldclim_m




#m@models$species$rf$`13`@object
# roc(m)
gui(worldclim_m)

predictions <- predict(worldclim_m, worldclim_bioc,
                       mean = F, filename ='all.img',overwrite=TRUE)

# predictions <- predict(m,bioc_hfp,species = 'rot',
#                        mean = T,filename ='rot.gri',overwrite=TRUE)


cl <- colorRampPalette(c('#3E49BB','#3498DB','yellow','orange','red','darkred'))

plot(predictions,col=cl(200))                   
# plot(predictions, breaks = c(seq(0,1,0.5)))    
# if this preds is a dataframe ,than output can be a csv
plot(predictions[[1:7]],)





#---------------ensemble methord----------

worldclim_en1 <- ensemble(m,worldclim_bioc,filename='worldclim_en1.img',
                setting=list(method='weighted',stat='tss',opt=2),overwrite = T)


cl <- colorRampPalette(c('#3E49BB','#3498DB','yellow','orange','red','darkred'))

worldclim_en_auc <-ensemble(m,bioc,filename = 'en_auc.img',
#                   setting = list(stat = 'auc', opt = 2),overwrite = T)


# mapview(en1)
plot(worldclim_en1) 
points(rot,cex=1,pch='*',col='black')



#--------------------forecast---------

biof <- raster::getData('CMIP5',var='bio',res= 10, rcp =85, year= 70, model = 'AC')
names(biof) <- names(bio)
biof <- crop(biof, country)
predictions_frost <- predict(m,biof,species = 'rot',
                             mean = T,filename ='rot_f.gri',overwrite=TRUE)
plot(predictions_frost)
enf <- calc(predictions_frost[[c(1,3,7)]],mean)
mapview(enf)
plot(stack(en,enf))
mapview(stack(en,enf))

#------------change-------

ch <- enf - en
cl <- colorRampPalette(c('red','orange','yellow','gray','green','blue'))
plot(ch, col=cl(100))

#-----------get threshould---------
th <- getEvaluation(m,stat='threshold')
mean(th[1:3,2])

df <- data.frame(as.data.frame(d),coordinates(d))
head(df)

pr <- raster::extract(en, df[,c('lon','lat')])
head(pr)
ev <- evaluates(df$rot, pr)
ev@threshold_based
th <- 0.6370913 

pa <- en

pa[] <- ifelse(pa[] >= th, 1, 0)

paf <- enf

paf[] <- ifelse(paf[] >= th, 1, 0)

plot(pa)
plot(paf)


pa.ch <- paf - pa
plot(pa.ch)


cl <- colorRampPalette(c('red','gray80','darkgreen'))
plot(pa.ch,col=cl(3))


#---------------
rcurve(m)

plot(getVarImp(m))

niche(bioc_hfp,en1,n=c('hfp','tmx'),col=cl(200))

