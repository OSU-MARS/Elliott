library(dplyr)
library(terra)

# reference data frame
#library(readr)
#xycrw.df.rf = tibble(read_rds(file.path(getwd(), "GIS/Trees/segmentation/SpeciesTrainingData.rds")))

trainingPolygons = vect(file.path(getwd(), "GIS/Trees/segmentation/SpeciesITC.shp"))

dataSourcePath = "D:/Elliott/GIS/DOGAMI/2021 OLC Coos County"
chm = rast(file.path(dataSourcePath, "CHM", "CHM.vrt"))
orthoimagery = rast(file.path(dataSourcePath, "orthoimage", "orthoimage.vrt"))
orthoimageCellSize = max(res(orthoimagery))

for (trainingPolygonIndex in 1:nrow(trainingPolygons))
{
  # extract training data for this polygon
  # For now, assume orthoimagery and CHM are in the same CRS and that orthoimagery is higher resolution.
  trainingPolygon = trainingPolygons[trainingPolygonIndex]
  polygonOrthoimage = crop(orthoimagery, trainingPolygon, touches = FALSE) # take only cells with centroids within the training polygon
  polygonChm = resample(crop(chm, buffer(trainingPolygon, orthoimageCellSize)), polygonOrthoimage)
  
  c(polygonChm, polygonOrthoimage) # stack CHM Z and RGB+NIR raster bands
}

for (i in 1:length(orthoimageFileNamesbypol)) {
  polygonChm = disagg(polygonChm, fact=2)
  bands.cz = c(ortho, polygonChm)
  #plot(ortho$NIR)
  #plot(z.var, add=TRUE)
  #rgb2g = (0.3 * ortho$R  + 0.59 *ortho$G + 0.11 * ortho$B)/ 65535 #RGB to gray
  #plot(rgb2g)
  
  x.var = terra::extract(bands.cz, trainingPolygons[i,], fun = NULL, method="simple",  exact=TRUE) #select only pixels that have the centroid intersected by the polygon
  
  #plot(z.var)
  
  idx.z = which(!is.na(x.var$Z))
  if (length(idx.z) > 0) {
    z.tr = myOtsu(x.var$Z[idx.z]) 
  } else {
    z.tr = myOtsu(x.var$Z) 
  }
  
  x.var.tr = x.var[x.var$Z>z.tr,]
  #class(x.var)
  #rgb2g = (0.3 * x.var$R  + 0.59 *x.var$G + 0.11 * x.var$B)/ 65535 #RGB to gray
  #intens.tr = myOtsu(rgb2g)
  #x.var.tr = x.var[rgb2g>intens.tr,]
  #sum(rgb2g>intens.tr)/length(rgb2g)
  #hist(rgb2g,30)
  
  
  
  tmp = matrix(apply(x.var.tr,2,function(x) quantile(x, seq(0,1,by = 0.10)))[,c("R","G","B","NIR","Z")], nrow=1)
  colnames(tmp) = c(paste0("R",seq(0,100, by = 10)), 
                    paste0("G",seq(0,100, by = 10)),
                    paste0("B",seq(0,100, by = 10)), 
                    paste0("NIR",seq(0,100, by = 10)),
                    paste0("Z",seq(0,100, by = 10)))
  tmp = as.data.frame(tmp)
  tmp[,c("Rm","Gm","Bm","NIRm","Zm")] = as.numeric(apply(x.var.tr[,c("R","G","B","NIR","Z")], 2, mean))
  tmp[,c("Rsd","Gsd","Bsd","NIRsd", "Zsd")] = as.numeric(apply(x.var.tr[,c("R","G","B","NIR","Z")], 2, sd))
  #indices computed using the formulas from "A visible band index for remote sensing leaf chlorophyll content at the canopy scale" by Hunt et al 2013
  tmp$NDVIm = (tmp$NIRm-tmp$Rm) / (tmp$NIRm+tmp$Rm) #NDVI
  tmp$NGRDIm = (tmp$Gm-tmp$Rm) / (tmp$Gm+tmp$Rm) #Normalized Green-Red Difference Index: Tucker 1979
  tmp$CIGm = tmp$NIRm/tmp$Gm - 1 #Gitelson et al 2003
  tmp$CVIm = tmp$NIRm * tmp$Rm / tmp$Gm^2 #Vinciniet al 2008
  tmp2 = matrix(apply(tmp, 2, function(x) quantile(x, seq(0,1,by = 0.10)))[,c("NDVIm","NGRDIm","CIGm","CVIm")], nrow=1)
  colnames(tmp2) = c(paste0("NDVI",seq(0,100, by = 10)), 
                     paste0("NGRDI",seq(0,100, by = 10)),
                     paste0("CIG",seq(0,100, by = 10)), 
                     paste0("CVI",seq(0,100, by = 10)))
  tmp2 = as.data.frame(tmp2)
  tmp=cbind(tmp,tmp2)
  tmp$Species = trainingPolygons[i,]$Species
  
  xycrw.df[[i]] = tmp
  
}  

print(paste0("Time to create Input data: ", difftime(Sys.time(), datain.t, units = "mins"), " min"))

xycrw.df.rf = do.call('rbind',xycrw.df) #creates the training data set to be used in RF for crown
xycrw.df.rf$Species = as.factor(xycrw.df.rf$Species)
saveRDS(xycrw.df.rf, "k:/People/MARS/Strimbu/TreeSeg2/TrainingData.RDS")

names(xycrw.df.rf)
xycrw.df.rf2= xycrw.df.rf #[,-c(45:55,60,65)] # no heights
#xycrw.df.rf2= xycrw.df.rf[,-c(1:44, 56:64)] #eliminate the absolute values for the reflectance


#xycrw.df.rf2$Species=as.factor(xycrw.df.rf2$Species)

#Random Forest using ranger
# set.seed(1969)
# fit.crw.rf = ranger(Species~., data=xycrw.df.rf2, importance="permutation")
# fit.crw.rf$confusion.matrix
# barplot(importance(fit.crw.rf))
# 
# ctrl <- trainControl(method = 'cv', 
#                     number = 10,
#                     classProbs = TRUE,
#                     savePredictions = TRUE,
#                     verboseIter = TRUE)
# 
# rfFit <- train(Species ~ ., 
#                data = xycrw.df.rf2, 
#                method = "ranger",
#                importance = "permutation", #***
#                trControl = ctrl,
#                verbose = T)
# 
# x.importance = varImp(rfFit)$importance
# x.importance$Var.names = row.names(x.importance)
# 
# x.importance = x.importance[order(x.importance$Overall, decreasing = TRUE),]
# barplot(x.importance$Overall[1:12], names=x.importance$Var.names[1:12])
# 
# saveRDS(fit.crw.rf,"k:/People/MARS/Strimbu/TreeSeg2/TrainingFIT.RDS")

## -- Cross validation
yhat.rf=rep(0,dim(xycrw.df.rf2)[1])
yhat.svm.r=rep(0,dim(xycrw.df.rf2)[1])
yhat.svm.l=rep(0,dim(xycrw.df.rf2)[1])
yhat.nn=rep(0,dim(xycrw.df.rf2)[1])

for (i in 1:dim(xycrw.df.rf2)[1]) {
  #i=1
  test.data = xycrw.df.rf2[i,]
  train_control <- trainControl(method="repeatedcv", number=10, repeats=3)
  ### Random Forest
  #fit.rf = train(Species ~., data = xycrw.df.rf2[-i,], method = "ranger", 
  #               trControl = train_control, tuneGrid = expand.grid(mtry=floor(sqrt(dim(xycrw.df.rf)[2]-1)),
  #                                                                 splitrule = 'gini',
  #                                                                 min.node.size = c(1:7)))
  #yhat.rf[i]=as.character(predict(fit.rf, newdata=test.data))
  
  ###### Neural Networks
  
  #fit.nnet <- train(Species ~., data = xycrw.df.rf2[-i,], method = "nnet", 
  #                  trControl = train_control, preProcess = c("center","scale"), verbose = FALSE,
  #                  tuneGrid = expand.grid(size=c(2, 3), decay = seq(0.05, 1, length.out=10)))
  #yhat.nn[i]=as.character(predict(fit.s, newdata = test.data))
  
  ##### Support Vector Machine
  
  fit.svm.r <- train(Species ~., data = xycrw.df.rf2[-i,], method = "svmRadial", #radial Basis Function Kernel
                     trControl = train_control, preProcess = c("center","scale"), tuneLength = 10)
  
  fit.svm.l <- train(Species ~., data = xycrw.df.rf2[-i,], method = "svmLinear", #linear kernel
                     trControl = train_control, preProcess = c("center","scale"), 
                     tuneGrid =expand.grid(C = c(0.01, 0.05,0.1, 1)))
  
  yhat.svm.r[i]=as.character(predict(fit.svm.r, newdata = test.data))
  yhat.svm.l[i]=as.character(predict(fit.svm.l, newdata = test.data))
}

cf.mat.nn=confusionMatrix(as.factor(yhat.nn),as.factor(xycrw.df.rf2$Species))
cf.mat.nn$byClass

cf.mat.rf=confusionMatrix(as.factor(yhat.rf),as.factor(xycrw.df.rf2$Species))
cf.mat.rf$byClass

cf.mat.svm.r=confusionMatrix(as.factor(yhat.svm.r),as.factor(xycrw.df.rf2$Species))
cf.mat.svm.r$byClass

cf.mat.svm.l=confusionMatrix(as.factor(yhat.svm.l),as.factor(xycrw.df.rf2$Species))
cf.mat.svm.l$byClass


#--CLASSIFY Species Crown Level----------------------------------------------------------------------#
RF.model = readRDS("k:/People/MARS/Strimbu/TreeSeg2/TrainingFIT.RDS")
rm()

#RF
#fit.rf <- train(Species ~., data = xycrw.df.rf2, method = "ranger", 
#                trControl = train_control)
#saveRDS(fit.rf,"k:/People/MARS/Strimbu/TreeSeg2/TrainingRF.RDS")

#+++++++++SVM
# --- Radial Kernell
# fit.svm.r <- train(Species ~., data = xycrw.df.rf2, method = "svmRadial", 
#                  trControl = train_control, preProcess = c("center","scale"), tuneLength = 10)
# whichTwoPct <- tolerance(fit.svm.r$results, metric = "Kappa", 
#                          tol = 1, maximize = TRUE) 
# 
# fit.svm.rf = fit.svm.r$results[whichTwoPct,] #identify best SVM parametrization
# 
# fit.svm.r <- train(Species ~., data = xycrw.df.rf2, method = "svmRadial", 
#                  trControl = train_control, preProcess = c("center","scale"), 
#                  tuneGrid = expand.grid(C = fit.svm.rf$C, sigma=fit.svm.rf$sigma) )
#saveRDS(fit.svm.r,"k:/People/MARS/Strimbu/TreeSeg2/TrainingSVMradial.RDS")

# ---- Linear kernel
fit.svm.l <- train(Species ~., data = xycrw.df.rf2, method = "svmLinear", 
                   trControl = train_control, preProcess = c("center","scale"), 
                   tuneGrid =expand.grid(C = c(0.01, 0.05,0.1, 1)))
whichTwoPct <- tolerance(fit.svm.l$results, metric = "Kappa", 
                         tol = 1, maximize = TRUE) 

fit.svm.lf = fit.svm.l$results[whichTwoPct,] #identify best SVM parametrization

fit.svm.l <- train(Species ~., data = xycrw.df.rf2, method = "svmLinear", 
                   trControl = train_control, preProcess = c("center","scale"), 
                   tuneGrid = expand.grid(C = fit.svm.lf$C) )

saveRDS(fit.svm.l,"k:/People/MARS/Strimbu/TreeSeg2/TrainingSVMlinear.RDS")

#shp.list =list.files("k:\\People\\MARS\\Strimbu\\TreeSeg2\\Results",pattern=".shp") #read the stands: trees and crowns

#shp.crown = shp.list[which(substr(shp.list,1,3) =="ITC")] #select only the crowns

# stands = st_read(".\\GIS\\ESRF_Stands062022Fixed.shp")
# plot(st_geometry(stands))
# cat(crs(stands))
# #head(stands)
# stands=st_transform(stands, 32610)

#crown.shp = st_read("d:\\Work\\ULC\\Results\\CrownsInvent.shp")
#cat(crs(crown.shp))

#===============================================================================
#Connect stands with orthophotos
stand.list = as.data.frame(list.files("k:\\People\\MARS\\Strimbu\\TreeSeg2\\Results\\", pattern = "ITCDalp_OSU_SID"))
names(stand.list)=c("stand")
stand.list$rec.id = seq.int(nrow(stand.list))
stand.list = stand.list[which(substr(stand.list$stand, nchar(stand.list$stand)-2,nchar(stand.list$stand))=="shp"),]

#trees = vect("k:\\ProcessedData\\ForestInvent\\Invent062023\\TreesFinal_062023.shp")

stand.list$osu_sid = substr(stand.list$stand, 17, nchar(stand.list$stand)-4)

ij.stand = matrix(0, ncol=length(orthoimageFileNames), nrow = dim(stand.list)[1])

time2ortho=Sys.time()
for (t in  1:nrow(stand.list)) {
  #t=2
  stand.tmp = stand.b[stand.b$OSU_SID == stand.list$osu_sid[t],]
  #plot(stand.tmp)
  for (tt in 1:length(orthoimageFileNames)) {
    ortho =rast(paste0("k:\\RawData\\Imagery\\NV5_4bands_2021_RayTrace\\MaxZSmooth\\",orthoimageFileNames[[tt]])) 
    ext.ortho = ext(ortho)
    tmp=try(terra::intersect(ext.ortho, stand.tmp), silent=TRUE)
    #tmp=try(crop(ortho,stand.tmp),silent=TRUE)
    if(class(tmp)=="try-error") {
      next
    } else {
      print(c(t,tt))
      ij.stand[t,tt]=1
    }  
    
  }
}
print(paste0("Time to id the orthophotos needed for each polygon: ", 
             difftime(Sys.time(),time2ortho, units = "mins"), " min"))
saveRDS(ij.stand, "k:\\People\\MARS\\Strimbu\\TreeSeg2\\Results\\ij_stand.RDS")

#feature extraction and species prediction 
#================================================
#FOREST LEVEL CLASSIFICATION
#run thru stands
#ij.stand =readRDS("k:\\People\\MARS\\Strimbu\\TreeSeg2\\Results\\ij_stand.RDS")

for (sid in 97:200) {#nrow(stand.list)) {
  #sid=2
  print(paste0("Execute stand #", sid, " ->OSU_SID = ", stand.list$stand[sid]))
  trainingPolygons = vect(paste0("k:\\People\\MARS\\Strimbu\\TreeSeg2\\Results\\", stand.list$stand[sid])) #select all crowns from stand osu_SID
  trainingPolygons = project(trainingPolygons, "epsg:6557")
  #cat(crs(trainingPolygons))
  #plot(trainingPolygons, add=TRUE)
  
  #================================================================================================
  id.ortho = which(ij.stand[sid,]==1)  
  rlist = list()
  for (ii in 1:length(id.ortho)) {
    rlist[[ii]] = rast(paste0("k:\\RawData\\Imagery\\NV5_4bands_2021_RayTrace\\MaxZSmooth\\",
                              orthoimageFileNames[id.ortho[ii]]))
  }
  rsrc <- sprc(rlist)
  ortho <- mosaic(rsrc)
  #plot(ortho$NIR)
  
  #================================================================================================
  #Execute classification INSIDE the stands
  xycrw.df = vector('list',dim(trainingPolygons)[1])
  #print(paste0("Number of crowns in the stand: ", dim(trainingPolygons)[1]))
  for (i in 1:nrow(trainingPolygons)) {
    #i=50
    z.var = terra::crop(chm,ext(trainingPolygons[i,]))
    z.var = disagg(z.var, fact=2)
    if (ext(ortho) != ext(z.var)) {
      ortho.c=terra::crop(ortho, ext(z.var))
    } 
    #plot(ortho.c$NIR)
    #plot(z.var, add=TRUE)
    
    bands.cz = c(ortho.c, z.var)
    
    x.var = terra::extract(bands.cz, trainingPolygons[i,], fun = NULL, method="simple",  exact=TRUE) #select only pixels that have the centroid intersected by the polygon
    
    if(sum(is.na(x.var))>0) {
      tmp=matrix(0,1,55)
      colnames(tmp) = c(paste0("R",seq(0,100, by = 10)), 
                        paste0("G",seq(0,100, by = 10)),
                        paste0("B",seq(0,100, by = 10)), 
                        paste0("NIR",seq(0,100, by = 10)),
                        paste0("Z",seq(0,100, by = 10)))
      tmp = as.data.frame(tmp)
      tmp[,c("Rm","Gm","Bm","NIRm","Zm")] = 0
      tmp[,c("Rsd","Gsd","Bsd","NIRsd", "Zsd")] = 0
      tmp$NDVIm = 0
      tmp$NGRDIm = 0
      tmp$CIGm = 0
      tmp$CVIm = 0
      tmp2 = matrix(apply(tmp, 2, function(x) quantile(x, seq(0,1,by = 0.10), na.rm=TRUE))[,c("NDVIm","NGRDIm","CIGm","CVIm")], nrow=1)
      colnames(tmp2) = c(paste0("NDVI",seq(0,100, by = 10)), 
                         paste0("NGRDI",seq(0,100, by = 10)),
                         paste0("CIG",seq(0,100, by = 10)), 
                         paste0("CVI",seq(0,100, by = 10)))
      tmp2 = as.data.frame(tmp2)
      tmp=cbind(tmp,tmp2)
      
      tmp$ID = trainingPolygons[i,]$treeID
      
      tmp$SpeciesHat.svm = "UN"
      xycrw.df[[i]] = tmp
      next
    } 
    
    idx.z = which(!is.na(x.var$Z))
    if (length(idx.z) > 0) {
      z.tr = myOtsu(x.var$Z[idx.z]) 
    } else {
      z.tr = myOtsu(x.var$Z)
    }
    
    x.var.tr = x.var[x.var$Z>z.tr & x.var$fraction>0.4,]
    
    tmp = matrix(apply(x.var.tr,2,function(x) quantile(x, seq(0,1,by = 0.10),na.rm=TRUE))[,c("R","G","B","NIR","Z")], nrow=1)
    colnames(tmp) = c(paste0("R",seq(0,100, by = 10)), 
                      paste0("G",seq(0,100, by = 10)),
                      paste0("B",seq(0,100, by = 10)), 
                      paste0("NIR",seq(0,100, by = 10)),
                      paste0("Z",seq(0,100, by = 10)))
    tmp = as.data.frame(tmp)
    tmp[,c("Rm","Gm","Bm","NIRm","Zm")] = as.numeric(apply(x.var.tr[,c("R","G","B","NIR","Z")], 2, mean))
    tmp[,c("Rsd","Gsd","Bsd","NIRsd", "Zsd")] = as.numeric(apply(x.var.tr[,c("R","G","B","NIR","Z")], 2, sd))
    #indices computed using the formulas from "A visible band index for remote sensing leaf chlorophyll content at the canopy scale" by Hunt et al 2013
    tmp$NDVIm = (tmp$NIRm-tmp$Rm) / (tmp$NIRm+tmp$Rm) #NDVI
    tmp$NGRDIm = (tmp$Gm-tmp$Rm) / (tmp$Gm+tmp$Rm) #Normalized Green-Red Difference Index: Tucker 1979
    tmp$CIGm = tmp$NIRm/tmp$Gm - 1 #Gitelson et al 2003
    tmp$CVIm = tmp$NIRm * tmp$Rm / tmp$Gm^2 #Vinciniet al 2008
    tmp2 = matrix(apply(tmp, 2, function(x) quantile(x, seq(0,1,by = 0.10), na.rm=TRUE))[,c("NDVIm","NGRDIm","CIGm","CVIm")], nrow=1)
    colnames(tmp2) = c(paste0("NDVI",seq(0,100, by = 10)), 
                       paste0("NGRDI",seq(0,100, by = 10)),
                       paste0("CIG",seq(0,100, by = 10)), 
                       paste0("CVI",seq(0,100, by = 10)))
    tmp2 = as.data.frame(tmp2)
    tmp=cbind(tmp,tmp2)
    
    tmp$ID = trainingPolygons[i,]$treeID
    #Random Forest
    #yhat.svm=as.character(predict(fit.rf, newdata = tmp))
    #tmp$SpeciesHat.rf = yhat.rf
    
    #Support Vector Machine
    if(sum(is.na(tmp))>0) {
      tmp$SpeciesHat.svm = "UN"
      xycrw.df[[i]] = tmp
      next
    } else {
      yhat.svm=predict(fit.svm.l, newdata = tmp)
      tmp$SpeciesHat.svm = yhat.svm
      xycrw.df[[i]] = tmp
    }
    
  }  
  xycrw.df.all =  do.call('rbind',xycrw.df) #creates the training data set to be used in SVM for crown
  xycrw.df.all$OSU_SID = stand.list$osu_sid[sid]
  #=====================================================================================
  
  trainingPolygons$SpeciesHat.svm = xycrw.df.all$SpeciesHat.svm
  writeVector(trainingPolygons, paste0("k:/People/MARS/Strimbu/TreeSeg2/Species/PredSVM/",
                                str_remove(stand.list$stand[sid],".shp"),"SVM.shp"), overwrite = TRUE)
  rm(ortho, rlist, tmp, tmp2);gc()
  
  saveRDS(xycrw.df.all,paste0("k:/People/MARS/Strimbu/TreeSeg2/Species/RDS_SVM/",
                              str_remove(stand.list$stand[sid],".shp"),".RDS"))
}


## one time setup
# Since gdalbuildvrt doesn't flow metadata (https://github.com/OSGeo/gdal/issues/3627) VRTs need manual editing
# after creation to add <VRTRasterBand/Description> elements naming raster bands. This most likely isn't important 
# to the CHM (or DSM or DTM) as terra defaults the band name to the file name but is needed for orthoimagery
# to set the band names to R, G, B, and NIR.
# create .vrt for CHM after chmJob.R has processed all tiles
chmSourcePath = "D:/Elliott/GIS/DOGAMI/2021 OLC Coos County/CHM"
chmFilePaths = file.path(chmSourcePath, list.files(orthoimageSourcePath, "\\.tif"))
vrt(chmFilePaths, file.path(chmSourcePath, "CHM.vrt"), overwrite = TRUE)

# create .vrt for orthoimages after orthoImageJob.R has processed all tiles
orthoimageSourcePath = "D:/Elliott/GIS/DOGAMI/2021 OLC Coos County/orthoimage"
orthoimageFilePaths = file.path(orthoimageSourcePath, list.files(orthoimageSourcePath, "\\.tif"))
vrt(orthoimageFilePaths, file.path(orthoimageSourcePath, "orthoimage.vrt"), overwrite = TRUE)
