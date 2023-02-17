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
library(car)
library(corrplot)



getwd()
#----------read data---------
disease_point<- read.csv("../../../datas/apple_diseases/datas/fu_lan_bing.csv",header = TRUE, check.names = F)

nrow(disease_point)
presence_data <- filter(disease_point, infestation_percentage == 1)
coordinates(presence_data ) <- ~ lon + lat
# write.csv(presence_data,"./presence_data.csv", row.names = FALSE)
absence_data <- filter(disease_point, infestation_percentage == 0)
# write.csv(absence_data,"./absence_data.csv", row.names = FALSE)

#spatial points data frame

coordinates(disease_point) <- ~ lon + lat

#env
envs <-stack('../../../datas/environment/envs.tif')
names(envs) <- c('pre','tmn','tmp','tmx','hfp','soc','prop','sand','silt','aspect','slope','curvature','dem')
plot(envs[[8]])


#-----------------fit modes-----------------

d <- sdmData(infestation_percentage
             ~., disease_point,
             predictors = envs,bg = list(method='gRandom',n=1000))
d

df <- as.data.frame(d) # train datas
write.csv(df,"./data_valsce.csv", row.names = FALSE)

#getmethodNames()
m <- sdm(infestation_percentage ~., 
         data=d, methods=c('glmnet','maxlike','mlp','rbf',
                           'gam','gbm','rpart','svm','rf')
         ,replications=c('sub','boot'),test.p=30,n=10,
         parallelSetting=list(ncore=12,method='parallel')) #replications : sub , cv, boot
m




#m@models$species$rf$`13`@object
# roc(m)
gui(m)

predictions <- predict(m,envs,
                       mean = T,filename ='predictions.img',overwrite=TRUE)
#predictions <- predict(bioc, m ,type = "revalonse",
#                     method = c('glm','gam','gbm','glmnet','maxlike','svm','rf'))

# predictions <- predict(m,envs,species = 'rot',
#                        mean = T,filename ='rot.gri',overwrite=TRUE)

#p1 <- predict(m, bio,filename='pr.img')

cl <- colorRampPalette(c('#3E49BB','#3498DB','yellow','orange','red','darkred'))

plot(predictions,col=cl(200))                   
# plot(predictions, breaks = c(seq(0,1,0.5)))    
# if this preds is a dataframe ,than output can be a csv
plot(predictions[[1:7]],)





#---------------ensemble methord----------

en1 <- ensemble(m,envs,filename='fu_lan_bing.img',
                setting=list(method='weighted',stat='tss',opt=2),overwrite = T)


# cl <- colorRampPalette(c('#3E49BB','#3498DB','yellow','orange','red','darkred'))
cl <- colorRampPalette(c('gray','blue','yellow','red'))
# en_auc <-ensemble(m,bioc,filename = 'en_auc.img',
#                   setting = list(stat = 'auc', opt = 2),overwrite = T)

delete_en1 <- en1
delete_en1[delete_en1 < 0] <- NA
# mapview(en1)
plot(delete_en1,col=cl2(200)) 
points(presence_data,cex=1,pch='*',col='black')

rf <- writeRaster(delete_en1,
                  filename="../../../results/fu_lan_bing_ensemble.tif",overwrite=TRUE)

#------climate change
biof <- raster::getData('CMIP5', var='bio', res=10, rcp=85, model='CN', year=70)

biof
plot(biof[[1]])
names(biof) <- names(bio)

en2 <- ensemble(m, biof, filename='enf.img',setting=list(method='weighted',stat='tss',opt=2))
#--------------
plot(en2)
cl2 <- colorRampPalette(c('#3498DB','yellow','orange','red','darkred'))
#------
plot(en1, col=cl(200))
plot(en2, col=cl(200))

#spg species distribution
proj4string(spg) <- projection(en1)

library(mapview)
mapview(en1,col.regions=cl(200)) + spg

#-----------
change <- en2 - en1
cl3 <- colorRampPalette(c('red','orange','yellow','gray','green','blue'))
plot(change,col=cl2(200))

#----
# presence_data
df
df <- data.frame(species=df$infestation_percentage,coordinates(d))

xy <- as.matrix(df[,c('lon','lat')])
head(xy)
nrow(xy)
p <- raster::extract(en1,xy)
head(p)
nrow(df)
length(p)

#病害发生点和不存在点的模型预测值

evaluates_point_data <- cbind(df,p)
write.csv(evaluates_point_data,"./evaluates_point_data.csv", row.names = FALSE)

ev <- evaluates(df$species,p)
ev@statistics

th <- ev@threshold_based$threshold[2]


#-----------------
plot(en1 > th, col=cl2(500)) 
proj4string(presence_data) <- projection(en1)
points(presence_data,cex=1,pch='*',col='black')

#---------------


pa1 <- raster(en1)

pa1[] <- ifelse(en1[] >= th, 1, 0)
plot(pa1)

pa2 <- raster(en1)

pa2[] <- ifelse(en2[] >= th, 1, 0)
plot(pa2)

chp <- pa2 - pa1
plot(chp,col=c('red','gray','blue'))

#---------------other analyse--------------

rcurve(m,id=1:4)

plot(getVarImp(m,method=c('glmnet','maxlike','mlp','rbf',
                          'gam','gbm','rpart','svm','rf')),
     main='relative variable importance(averaged)')

# Mapping Ecological Niche using selected two variables

niche(envs,h=en1,n=c('hfp','tmn'),col=cl(10),rnd=2)


getwd()



envs


