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

getwd()
#----------read data---------
disease_point<- read.csv("../../../datas/apple_diseases/apple_diseases_sdm/fu_lan_bing.csv",
                         header = TRUE, check.names = F)
nrow(disease_point)
presence_data <- filter(disease_point, infestation_percentage == 1)
write.csv(presence_data,"./presence_data.csv", row.names = FALSE)
absence_data <- filter(disease_point, infestation_percentage == 0)
write.csv(absence_data,"./absence_data.csv", row.names = FALSE)
#spatial points data frame

coordinates(disease_point) <- ~ lon + lat

#env
bioc_hfp <-stack('../../../datas/environment/2020_2017_mean/bioc_hfp.tif')
names(bioc_hfp) <- c('pre','tmn','tmp','tmx','hfp')
plot(bioc_hfp[[2]])
points(presence_data,cex=1,pch='*',col='red')
points(absence_data,cex=1,pch='*',col='black')


#-----------------fit modes-----------------

#only presence
nrow(presence_data)
d <- sdmData(infestation_percentage
             ~pre + tmn + tmx + tmp + hfp, disease_point,
             predictors = bioc_hfp,bg = list(method='gRandom',n=1000))
d
#presence-absence

# d <- sdmData(infestation_percentage
#              ~pre + tmn + tmx + tmp + hfp, disease_point,
#              predictors = bioc_hfp)
# d

# saving datas

df <- as.data.frame(d) # train datas
write.csv(df,"./op_data_fu_lan_bing.csv", row.names = FALSE)

# getmethodNames()
m <- sdm(infestation_percentage ~pre + tmn + tmx + tmp + hfp, 
         data=d, methods=c('glmnet','maxlike',
                           'gam','gbm','rpart','svm','rf')
         ,replications='boot',test.p=30,n=10,
         parallelSetting=list(ncore=12,method='parallel')) #replications : sub , cv, boot
m




#m@models$species$rf$`13`@object
# roc(m)
gui(m)

predictions <- predict(m,bioc_hfp,
                       mean = T,filename ='predictions_only_preesence.img',overwrite=TRUE)

# predictions <- predict(m,bioc_hfp,species = 'rot',
#                        mean = T,filename ='rot.gri',overwrite=TRUE)


cl <- colorRampPalette(c('#3E49BB','#3498DB','yellow','orange','red','darkred'))

plot(predictions,col=cl(200))                   
# plot(predictions, breaks = c(seq(0,1,0.5)))    
# if this preds is a dataframe ,than output can be a csv
plot(predictions[[1:7]],)





#---------------ensemble methord----------

ensemble_only_presence <- ensemble(m,bioc_hfp,filename='ensemble_only_presence.img',
                setting=list(method='unweighted',stat='tss',opt=2),overwrite = T)


cl <- colorRampPalette(c('#3E49BB','#3498DB','yellow','orange','red','darkred'))

# en_auc <-ensemble(m,bioc,filename = 'en_auc.img',
#                   setting = list(stat = 'auc', opt = 2),overwrite = T)


# mapview(en1)
plot(ensemble_only_presence,col=cl(200))  
points(presence_data,cex=1,pch='*',col='black')



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

