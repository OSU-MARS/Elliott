# Code for Species identification from Planet imagery
# authors: Liviu and subsemantu'
# May 5, 2022
rm(list=ls())

library(sp)
library(sf)
library(terra)
library(raster)
library(glcm)
library(ranger) # used for random forest
library(stars) #compatibility between sf and raster
library(caret)


#---FUNCTIONS---------------------------------------------------------------------#

multiFocal <- function(x, w=matrix(1, nr=3, nc=3), ...) {
  if(is.character(x)) {
    x <- brick(x)
  }
  # The function to be applied to each individual layer
  fun <- function(ind, x, w, ...){
    focal(x[[ind]], w=w, na.rm=TRUE,...)
  }
  n    <- seq(nlayers(x))
  list <- lapply(X=n, FUN=fun, x=x, w=w, ...)
  out  <- stack(list)
  names(out) <- names(x)
  return(out)
}

setwd("k:\\People\\MARS\\Strimbu\\TreeSeg2\\")
getwd()

#### ------ INPUT DATA ----#

#chm=rast("ULC_CHM_2fts.tif") #import CHM to be classified
#chm=as(chm, "Raster")
#chm=chm[[1]]
#crs(chm)=crs(rast("ULC_CHM_2fts.tif"))

loff =rast("k:\\RawData\\Imagery\\Planet\\composite_2022_02_09_rect.tif")
lon =rast("k:\\RawData\\Imagery\\Planet\\composite_2022_08_16_rect.tif")
names(loff)[1:4]=paste0("loff", substr(names(loff)[1:4], nchar(names(loff)[1:4])-1, nchar(names(loff)[1:4])))
names(lon)[1:4]=paste0("lon", substr(names(lon)[1:4], nchar(names(lon)[1:4])-1, nchar(names(lon)[1:4])))
chm=rast("k:\\People\\MARS\\Strimbu\\TreeSeg2\\ESSRF_CHM_60cm.tif")

planet=c(loff,lon)

#names(planet)[1:4]=paste0(names(planet)[1:4],"_","lon")
#names(planet)[-c(1:4)]=paste0(names(planet)[-c(1:4)],"_","loff")
planet$loff_5=(planet$loff_4-planet$loff_3) / (planet$loff_4 + planet$loff_3)
planet$lon_5=(planet$lon_4-planet$lon_3) / (planet$lon_4 + planet$lon_3)
planet = planet[[c(1:4,9,5:8,10)]] # order the layers
planet
planet=project(planet,"epsg:32610") # re-project the image coordinate system using terra package

planet1=resample(planet, rast(ext=ext(planet), resolution=1, crs=crs(planet)), method="near")

#train.shp=read_sf("d:\\Work\\ULC\\GIS\\ULC_GroundTreesTrain_DBHCrown.shp")
train.shp=st_read("./GIS/Species_aoi.shp")
plot(st_geometry(train.shp), col="red")
aggregate(train.shp$Area_m2 ~ train.shp$Species, FUN = sum)

#### - RANDOM FOREST TRAINING ---- #

#-----Create Training data-------------------------------------------------------------------#

xycrw.df = vector('list',dim(train.shp)[1])

for(i in 1:dim(train.shp)[1]) {
  #i=11
  print(dim(train.shp)[1]-i)
  polyg.tmp = train.shp[i,]
  #st_crs(polyg.tmp)
  #plot(st_geometry(polyg.tmp))
 
  rsp.tmp = extract(planet1, polyg.tmp, fun=NULL, method="simple",  exact=TRUE) #select only pixels that have the centroid intersected by the polygon
  
  #create a new dataset with two predictors mean and Std.dev.
  
  rsc.tmp.mu=as.data.frame(t(apply(rsp.tmp[,1:10],2,mean, na.rm=TRUE)))
  rsc.tmp.mad=as.data.frame(t(apply(rsp.tmp[,1:10],2,mad, na.rm=TRUE)))
  names(rsc.tmp.mu)=paste0(names(rsc.tmp.mu),"_mu")
  names(rsc.tmp.mad)=paste0(names(rsc.tmp.mad),"_mad")
  rsc.tmp=cbind(rsc.tmp.mu,rsc.tmp.mad)
  
  rsc.tmp$ID = polyg.tmp$Id
  rsp.tmp$ID = polyg.tmp$Id
  rsc.tmp$Species = polyg.tmp$Species
  rsp.tmp$Species = polyg.tmp$Species
  xycrw.df[[i]] = rsc.tmp 
 
  
}
xycrw.df.rf = do.call('rbind',xycrw.df) #creates the training data set to be used in RF for crown
xycrw.df.rf$Species = as.factor(xycrw.df.rf$Species)
xycrw.df.rf2= xycrw.df.rf
xycrw.df.rf2$Species=as.character(xycrw.df.rf2$Species)
#xycrw.df.rf2$Species[xycrw.df.rf2$Species=="WC"] = "DF" 
xycrw.df.rf2$Species = as.factor(xycrw.df.rf2$Species)
#head(xycrw.df.rf2)
#xycrw.df.rf2 = xycrw.df.rf2[,-which(names(xycrw.df.rf2)=="ID")]
#summary(xycrw.df.rf)

fit.crw.rf = ranger(Species~.,data=xycrw.df.rf2)

xycrw.df.rf.bin = xycrw.df.rf2
xycrw.df.rf.bin$Species=ifelse(xycrw.df.rf.bin$Species=="DF",1,0)

planet.dif=xycrw.df.rf.bin[,c(1:5)]-xycrw.df.rf.bin[,c(6:10)]
planet.dif[,c("Species")]=xycrw.df.rf.bin[,c("Species")]

#Logistic regression on the difference between seasons-alternative to Random Forest
fit.glm.crw = glm(Species~.,data=planet.dif, family = binomial)
summary(fit.glm.crw)
yhat = predict(fit.glm.crw, type="response", newdata=planet.dif)
summary(yhat)
hist(yhat)

yhat=ifelse(yhat>0.5,1,0)
table(yhat,planet.dif$Species)

## -- Cross validation
yhat=rep(0,dim(xycrw.df.rf2)[1])

for (i in 1:dim(xycrw.df.rf2)[1]) {
  #i=1
  test.data = xycrw.df.rf2[i,]
  fit.crw.rf.cv = ranger(Species~.,data=xycrw.df.rf2[-i,])
  yhat[i]=as.character(predict(fit.crw.rf.cv, data=test.data)$predictions)
  
}

cf.mat=confusionMatrix(as.factor(yhat),as.factor(xycrw.df.rf2$Species))
cf.mat$byClass


fit.crw.rf = ranger(Species~.,data=xycrw.df.rf2)
saveRDS(fit.crw.rf,"./Results/RandForest.RDS")


RF.model = readRDS("./Results/RandForest.RDS")
rm()

#--CLASSIFY Species Crown Level----------------------------------------------------------------------#


#shp.list =list.files("k:\\People\\MARS\\Strimbu\\TreeSeg2\\Results",pattern=".shp") #read the stands: trees and crowns

#shp.crown = shp.list[which(substr(shp.list,1,3) =="ITC")] #select only the crowns

# stands = st_read(".\\GIS\\ESRF_Stands062022Fixed.shp")
# plot(st_geometry(stands))
# cat(crs(stands))
# #head(stands)
# stands=st_transform(stands, 32610)

#crown.shp = st_read("d:\\Work\\ULC\\Results\\CrownsInvent.shp")
#cat(crs(crown.shp))

trees1 = st_read("./Results/Compiled/TreesInventCompl500lt290.shp")
trees2 = st_read("./Results/Compiled/TreesInventCompl1000lt290.shp")
trees3 = st_read("./Results/Compiled/TreesInventCompl1500lt290.shp")
trees4 = st_read("./Results/Compiled/TreesInventCompl1970lt290.shp")
trees.nn =rbind(trees1, trees2, trees3, trees4)
trees.nn=st_transform(trees.nn,32610)
cat(crs(trees.nn))
stand.list=vector("list", length(unique(trees.nn$OSU_SID)))

osu_sid = unique(trees.nn$OSU_SID)

 for (sid in 1:100){ #length(shp.crown)) {
   #sid=2
  #c=1
  #crown.shp = st_read(paste0("k:/People/MARS/Strimbu/TreeSeg2/Results/",shp.crown[sid]))
  #crown.shp=st_transform(crown.shp,32610)
  #crs(crown.shp)
  #crown.shp$OSU_SID = as.integer(substr(shp.crown[sid],17,nchar(shp.crown[sid])-4))
  crown.shp =  st_buffer(trees.nn[which(trees.nn$OSU_SID == osu_sid[sid]),],3) #buffer of 3 m
  
  #xycrw.df = vector('list',dim(crown.shp)[1])
  xycrw.df = vector('list', length(osu_sid))
   
#create the data to be classified with Random Forests
#for(i in 1:dim(crown.shp)[1]) {
for(i in 1:dim(crown.shp)[1]) {
  #i=2
  #print(paste0("Crowns left to be executed for OSU_SID ",substr(shp.crown[sid],17,nchar(shp.crown[sid])-4) ," (stand#",sid,") : ",dim(crown.shp)[1]-i))
  print(paste0("Trees left from stand ", osu_sid[sid], " (#", sid ,")"," for which species is computed # ", dim(crown.shp)[1]-i))
  #polyg.tmp = crown.shp[i,] #for using crown
  polyg.tmp = crown.shp[i,]
  #plot(st_geometry(polyg.tmp))
  rsp.tmp = extract(planet1, polyg.tmp, fun=NULL, method="simple",  exact=TRUE) #select only pixels that have the centroid intersected by the polygon
  
  #create a new dataset with two predictors mean and mad
  
  rsc.tmp.mu=as.data.frame(t(apply(rsp.tmp[,1:10],2,mean, na.rm=TRUE)))
  rsc.tmp.mad=as.data.frame(t(apply(rsp.tmp[,1:10],2,mad, na.rm=TRUE)))
  names(rsc.tmp.mu)=paste0(names(rsc.tmp.mu),"_mu")
  names(rsc.tmp.mad)=paste0(names(rsc.tmp.mad),"_mad")
  rsc.tmp=cbind(rsc.tmp.mu,rsc.tmp.mad)
  
  rsc.tmp$ID = polyg.tmp$treeID
  
  xycrw.df[[i]] = rsc.tmp 
 }
xycrw.df.rf = do.call('rbind',xycrw.df) #creates the  data set to be classified with RF for crown
#names(xycrw.df.rf)

#table(is.na(xycrw.df.rf[1]))

yhat=as.character(predict(RF.model, data=xycrw.df.rf[,1:20])$predictions)
#yhat = ifelse(yhat == 0,"HW","DF")

crown.shp$Species = 0
crown.shp$Species=yhat
table(crown.shp$Species)


#trees.tmp = intersect(vect(trees.nn), vect(crown.shp[,c("treeID","Species")])) #use terra package for intersection

trees.tmp = merge(st_drop_geometry(crown.shp[,c("OSU_SID","treeID","Species")]),st_drop_geometry(trees.nn), by=c("OSU_SID", "treeID"))
length(trees.tmp$OSU_SID) # check validity of the intersection
trees.tmp$Ht_ft = trees.tmp$zmax
trees.tmp$Area_ft2 = trees.tmp$area
trees.tmp=trees.tmp[,-c(4,25)]
trees.tmp = trees.tmp[,c(1:3,27, 4:26)]

#trees.tmp$Ht2BC_ft = trees.tmp$zq65
#trees.tmp$Species = "DF"


trees.tmp$DBH_in=0
trees.tmp$DBH_in[which(trees.tmp$Species=="DF")] = -47.2424 + 46.2518 * exp(0.00336 * trees.tmp$Ht_ft[which(trees.tmp$Species=="DF" )])
trees.tmp$DBH_in[which(trees.tmp$Species=="DF" & trees.tmp$Ht_ft >= 280)] = 90
#trees.tmp$DBH_in[which(trees.tmp$Species=="WH")] = -105.1 +104.3*exp(0.00186* trees.tmp$Ht_ft[which(trees.tmp$Species=="WH" )])
#trees.tmp$DBH_in[which(trees.tmp$Species=="WH" & trees.tmp$Ht_ft >= 260)] = 70
trees.tmp$DBH_in[which(trees.tmp$Species=="HW")] = -24.2932+ 28.4565*exp(0.00457* trees.tmp$Ht_ft[which(trees.tmp$Species=="HW" )])
trees.tmp$DBH_in[which(trees.tmp$Species=="HW" & trees.tmp$Ht_ft >= 230)] = 60

#summary(trees.tmp$Ht_ft)
#summary(trees.tmp$DBH_in)
#table(trees.tmp$Species)

stand.list[[sid]]<-trees.tmp
}


tree.final=do.call("rbind", stand.list)
table(tree.final$Species)
tree.final=st_as_sf(tree.final, coords = c("X","Y"),crs=6557) #the coordinates of of the original trees were in 6557

st_write(tree.final,"./Results/Compiled/TreesFinal100.shp", append=TRUE)




