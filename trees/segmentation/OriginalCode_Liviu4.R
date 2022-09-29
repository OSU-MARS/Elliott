# R code for landscape level individual tree crown segmentation 
# Original code with all section explicitly written with the help from Dr. Liviu Ene, Skogforsk SE.
library(stringr)
library(sp)
library(sf)
library(terra)
library(raster)
library(stars) #compatibility between sf and raster
library(lidR)
library(rgeos)
library(FNN)


setwd("k:\\People\\MARS\\strimbu\\Inventory\\Liviu") 
getwd()
esrf.inv15 = st_read("StandsJun29_Fixed.shp") #boundary stands ESRF
#plot(st_geometry(esrf.inv15))
#head(esrf.inv15)
#esrf.inv15$ac = esrf.inv15$Shape_Area / 43560 # check the surface of the polygons
esrf.inv15$T4I = esrf.inv15$TPA_Total * esrf.inv15$Acres
esrf.inv15 = st_transform(esrf.inv15, 6557)
#st_crs(esrf.inv15)

esrf.inv15.id = esrf.inv15[,-(3:44)]
esrf.inv15.id = esrf.inv15.id[,-(1)]
#head(esrf.inv15.id)
#st_crs(esrf.inv15.id)

tree.files=dir(path="k:\\People\\MARS\\Strimbu\\Inventory\\Liviu\\Results", pattern=".shp") #read the las files from the working directory into a list
tree.files.id = str_remove_all(tree.files,".shp")
tree.files.id = str_sub(tree.files.id,start=18,end=-1)#nchar(tree.files.id))

tmp2 = vector("list",length(tree.files.id))

for (i in 1:length(tree.files.id)) {
  #i=875
  tmp=tree.files.id[i]
  if (str_detect(tmp,"_")) {
    tmp2[[i]] = str_split_fixed(tmp,"_",2)[1]
  } else {
    tmp2[[i]]=tmp
  }
}

tmp2 = unlist(tmp2)
tmp2 = sort(unique(tmp2))
tmp2.df=as.data.frame(tmp2)
not.exe=esrf.inv15$StandID[which(!esrf.inv15$StandID %in% tmp2)]


#-------execute inventory at plot level ----#
assess.dif=data.frame("StandID"=numeric(),"T4Lidar"=numeric(),  "T4Inv"=numeric(), 
                      "Diff.Lidar.Inv.abs"=numeric(), "Diff.Lidar.Inv.rel"=numeric())
for (i in not.exe[1:150]) {
  #i=1661
  print(paste0("Execute Stand: ", esrf.inv15.id$StandID[which(esrf.inv15$StandID==i,)]))
  poly.inv0=esrf.inv15.id[which(esrf.inv15.id$StandID==i),]
  
  #plot(st_geometry(poly.inv))
  time.norm.s=Sys.time()
  if (dim(poly.inv0)[1]>1) {
    print(paste0("Number of polygons with the same StandID: ",dim(poly.inv0)[1]))
    standid.t=poly.inv0$StandID
    standid.t=paste0(standid.t,"_",1:length(standid.t))
    poly.inv0$StandID=standid.t
    las = clip_roi(las.files, poly.inv0)
    for (k in 1:length(standid.t)) {
      # k=1
      poly.inv=st_buffer(poly.inv0[k,], 30) #creates a buffer to  eliminate boundary effects
      las.k=las_update(las[[k]]) #update the header
      dtm=crop(dtm.nv5,ext(las.k)) #crop the large DTM to the las extent
      las.n=normalize_height(las.k, dtm)
      las.n = filter_poi(las.n, Classification==1L, Z>= 0, Z<= 400)
      #rm(las);gc()
      chm <- grid_canopy(las.n, res=2, algorithm=p2r(subcircle=0)) #compute the canopy height model using a pixel of 0.5 units, in this case meters
      #plot(chm)
      time.norm.e = Sys.time()
      print(paste0("Normalization time: ", difftime(time.norm.e, time.norm.s, units="secs")))
      
      
      #smooth the CHM using Gaussian filtering
      gauss.param=focalMat(chm, c(1,5),"Gauss")
      #chm.smooth=terra::focal(chm[[1]],w=gauss.param, na.rm=TRUE)
      chm.smooth=chm
      #plot(chm.smooth)
      poly.inv.h=quantile(chm, probs = c(0.1, 0.85, 0.95), type=7,names = FALSE)
      #poly.inv.h[[2]]
      
      lmf.f <- function(x) { x * log((poly.inv.h[[3]]-poly.inv.h[[2]]))/25 + 
          (poly.inv.h[[3]]-poly.inv.h[[2]])/10 }
      #lmf.f <- function(x) { x * 0.01 + 6 }
      local.max=find_trees(chm.smooth,lmf(lmf.f, hmin=poly.inv.h[[1]]), uniqueness = "bitmerge")  #identifies the tree tops from the CHM using the local maximum filter with a window of 6 feets and height >8 ft
      #plot(local.max, pch=20, add=TRUE)
      
      time.seg.s=Sys.time()
      ht.thres=poly.inv.h[[1]]
      cr.max=poly.inv.h[[3]] / 8
      algo <- dalponte2016(chm.smooth, local.max, th_tree=ht.thres, max_cr=cr.max) #create an R object that stores the algorithm to be used for tree segmentation
      #algo=watershed(chm, th_tree=15)
      #algo=li2012(dt1=1, dt2=2, Zu=30, hmin=10, speed_up = 8)
      #algo=silva2016(chm.smooth,local.max, max_cr_factor = 0.5, exclusion = 0.3, ID="treeID")
      las.t=segment_trees(las.n,algorithm=algo, attribute="treeID") # segment point cloud: identifies individual tree crowns
      
      time.seg.e = Sys.time()
      print(paste0("Segmentation time: ", difftime(time.seg.e, time.seg.s, units="secs")))
      
      
      #Remove  las.n to free memory
      #rm(las.n);gc()
      
      itc <- delineate_crowns(las.t, type = c("bbox"))  #create a polygon for each crown (convex hull)
      itc@data$area=as.numeric(st_area(st_as_sf(itc)))
      #hist(itc@data$area,100)
      
      tree.est=  tree_metrics(las.t, func=~list(X=mean(X), Y=mean(Y), Z=max(Z), # compute the location and height of each tree
                                                z60=quantile(Z,probs=0.60), z70=quantile(Z,probs=0.70), #compute percentiles
                                                z80=quantile(Z,probs=0.80), z90=quantile(Z,probs=0.90), 
                                                z95=quantile(Z,probs=0.95), z99=quantile(Z,probs=0.99)), 
                              attribute="treeID")
      #join the lidar metrics with crown info
      tree.est=merge(tree.est,itc,by=c("treeID"))
      
      #remove las.t file to free memory
      rm(las.t);gc()
      
      #--------Neighboring elimination process ----#
      #print("Start the neigboring elimination process ")
      time.nn.s=Sys.time()
      min.ht = poly.inv.h[[1]]
      min.cr.area = 15
      tree.est.t=subset(tree.est, tree.est@data$area > min.cr.area,
                        tree.est@data$Z99 > min.ht) #eliminate trees with crown radius smaller than 10 feet or 50 sq.ft
      
      ######eliminate trees that have the top less than a threshold distance apart (like 7 ft)
      #compute nearest neighbor distance, up to the farthest 3rd tree
      #nnwhich identifies the index of the tree in a list not the tree number. The matching of index with treeID should be done explicitly
      
      trees.df=st_drop_geometry(st_as_sf(tree.est))
      
      trees.nn3=get.knn(trees.df[,c("X","Y")], k=3, algorithm = "brute")
      trees.nn3=as.data.frame(cbind(trees.nn3$nn.index,trees.nn3$nn.dist))
      names(trees.nn3)=c("NN1", "NN2", "NN3", "NN1Dist","NN2Dist", "NN3Dist")
      
      
      trees.nnz=matrix(0, nrow(trees.df),3)
      trees.nnx=matrix(0, nrow(trees.df),3)
      trees.nny=matrix(0, nrow(trees.df),3)
      trees.nnindex=matrix(0, nrow(trees.df),3)
      trees.nndist=matrix(0, nrow(trees.df),3)
      
      #for each tree compute the xyz, the index and the distance to nearest 3 neighbors
      for (j in 1:nrow(trees.df)){
        trees.nnz[j,]=trees.df$Z[as.numeric(trees.nn3[j,c("NN1","NN2","NN3")])]
        trees.nnx[j,]=trees.df$X[as.numeric(trees.nn3[j,c("NN1","NN2","NN3")])]
        trees.nny[j,]=trees.df$Y[as.numeric(trees.nn3[j,c("NN1","NN2","NN3")])]
        trees.nnindex[j,]=as.numeric(trees.nn3[j,c("NN1","NN2","NN3")])
        trees.nndist[j,]=as.numeric(trees.nn3[j,c("NN1Dist","NN2Dist","NN3Dist")])
      }
      
      
      #rename the columns within the matrix
      colnames(trees.nnz)=c("Z1","Z2","Z3")
      colnames(trees.nnx)=c("X1","X2","X3")
      colnames(trees.nny)=c("Y1","Y2","Y3")
      colnames(trees.nnindex)=c("NN1","NN2","NN3")
      colnames(trees.nndist)=c("NN1Dist","NN2Dist","NN3Dist")
      
      trees.nn=cbind(trees.nnindex, trees.nndist, trees.nnz, trees.nnx, trees.nny, trees.df)
      
      
      min.d2nt=6.5
      trees.nn$FP=0
      for (j in 1:nrow(trees.nn)){
        #j=2
        tmp=trees.nn[j,c("NN1Dist", "NN2Dist","NN3Dist", "Z1", "Z2","Z3")]
        z.tmp=tmp[4:6][which(tmp[1:3]<=min.d2nt)]
        if (length(z.tmp)==0) {next}
        if(trees.nn$Z[j] <= max(z.tmp)) {
          trees.nn$FP[j] = 1
        }
      }
      
      trees.final = trees.nn[which(trees.nn$FP==0 & trees.nn$area>20),]
      trees.final$StandID=i
      
      time.nn.e = Sys.time()
      print(paste0("Neighbors Elimination time: ", difftime(time.nn.e, time.nn.s, units="secs")))
      
      #sum(trees.nn$NN1Dist<min.d2nt)
      #sum(trees.nn$FP)
      
      #create the final tree measurements file
      trees.shp=trees.nn
      trees.shp=subset(trees.shp,trees.shp$FP==0)
      coordinates(trees.shp)=~X+Y
      trees.shp=st_as_sf(trees.shp)
      st_crs(trees.shp)=st_crs(itc)
      trees.shp=st_transform(trees.shp, 6557)
      #crs(trees.shp)
      poly.inv.f=st_intersection(trees.shp, 
                                 st_geometry(poly.inv))
      poly.inv.f$StandID=poly.inv$StandID
      #st_crs(poly.inv.f)
      #plot(st_geometry(trees.shp), add=TRUE)
      #plot(st_geometry(poly.inv.f)), add=TRUE)
      
      st_write(poly.inv.f, paste0(getwd(),"/Results/TSegDalp_StandID_", poly.inv$StandID ,".shp"), append=FALSE)#write the segmented trees as a shapefile  
      
      diff.Lidar.Inv =as.data.frame(cbind(poly.inv$StandID, 
                                          nrow(poly.inv.f), poly.inv$T4I,
                                          floor(nrow(poly.inv.f)- poly.inv$T4I),
                                          -1+nrow(poly.inv.f) / poly.inv$T4I))
      names(diff.Lidar.Inv)=c("StandID","T4Lidar", "T4Inv", "Diff.Lidar.Inv.abs", "Diff.Lidar.Inv.rel")                         
      assess.dif=rbind(assess.dif, diff.Lidar.Inv)
      
      
    }
  } else {
    
    las = clip_roi(las.files, poly.inv0)
    plot(st_geometry(poly.inv0))
    las=las_update(las) #update the header
    #dtm=rasterize_terrain(las,res=2,algorithm=tin())
    dtm=crop(dtm.nv5,ext(las)) #crop the large DTM to the las extent
    las.n=normalize_height(las, dtm)
    las.n = filter_poi(las.n, Classification==1L, Z>= 0, Z<= 400)
    #rm(las);gc()
    chm <- grid_canopy(las.n, res=2, algorithm=p2r(subcircle=0)) #compute the canopy height model using a pixel of 0.5 units, in this case meters
    #plot(chm)
    time.norm.e = Sys.time()
    print(paste0("Normalization time: ", difftime(time.norm.e, time.norm.s, units="secs")))
    
    
    #smooth the CHM using Gaussian filtering
    gauss.param=focalMat(chm, c(1,5),"Gauss")
    #chm.smooth=terra::focal(chm[[1]],w=gauss.param, na.rm=TRUE)
    chm.smooth=chm
    #plot(chm.smooth)
    poly.inv.h=quantile(chm, probs = c(0.1, 0.85, 0.95), type=7,names = FALSE)
    #poly.inv.h[[2]]
    
    lmf.f <- function(x) { x * log((poly.inv.h[[3]]-poly.inv.h[[2]]))/25 + 
        (poly.inv.h[[3]]-poly.inv.h[[2]])/10 }
    #lmf.f <- function(x) { x * 0.01 + 6 }
    local.max=find_trees(chm.smooth,lmf(lmf.f, hmin=poly.inv.h[[1]]), uniqueness = "bitmerge")  #identifies the tree tops from the CHM using the local maximum filter with a window of 6 feets and height >8 ft
    #plot(local.max, pch=20, add=TRUE)
    
    time.seg.s=Sys.time()
    ht.thres=poly.inv.h[[1]]
    cr.max=poly.inv.h[[3]] / 8
    algo <- dalponte2016(chm.smooth, local.max, th_tree=ht.thres, max_cr=cr.max) #create an R object that stores the algorithm to be used for tree segmentation
    #algo=watershed(chm, th_tree=15)
    #algo=li2012(dt1=1, dt2=2, Zu=30, hmin=10, speed_up = 8)
    #algo=silva2016(chm.smooth,local.max, max_cr_factor = 0.5, exclusion = 0.3, ID="treeID")
    las.t=segment_trees(las.n,algorithm=algo, attribute="treeID") # segment point cloud: identifies individual tree crowns
    
    time.seg.e = Sys.time()
    print(paste0("Segmentation time: ", difftime(time.seg.e, time.seg.s, units="secs")))
    
    
    #Remove  las.n to free memory
    #rm(las.n);gc()
    
    itc <- delineate_crowns(las.t, type = c("bbox"))  #create a polygon for each crown (convex hull)
    itc@data$area=as.numeric(st_area(st_as_sf(itc)))
    #hist(itc@data$area,100)
    
    tree.est=  tree_metrics(las.t, func=~list(X=mean(X), Y=mean(Y), Z=max(Z), # compute the location and height of each tree
                                              z60=quantile(Z,probs=0.60), z70=quantile(Z,probs=0.70), #compute percentiles
                                              z80=quantile(Z,probs=0.80), z90=quantile(Z,probs=0.90), 
                                              z95=quantile(Z,probs=0.95), z99=quantile(Z,probs=0.99)), 
                            attribute="treeID")
    #join the lidar metrics with crown info
    tree.est=merge(tree.est,itc,by=c("treeID"))
    
    #remove las.t file to free memory
    rm(las.t);gc()
    
    #--------Neighboring elimination process ----#
    #print("Start the neigboring elimination process ")
    time.nn.s=Sys.time()
    min.ht = poly.inv.h[[1]]
    min.cr.area = 15
    tree.est.t=subset(tree.est, tree.est@data$area > min.cr.area,
                      tree.est@data$Z99 > min.ht) #eliminate trees with crown radius smaller than 10 feet or 50 sq.ft
    
    ######eliminate trees that have the top less than a threshold distance apart (like 7 ft)
    #compute nearest neighbor distance, up to the farthest 3rd tree
    #nnwhich identifies the index of the tree in a list not the tree number. The matching of index with treeID should be done explicitly
    
    trees.df=st_drop_geometry(st_as_sf(tree.est))
    
    trees.nn3=get.knn(trees.df[,c("X","Y")], k=3, algorithm = "brute")
    trees.nn3=as.data.frame(cbind(trees.nn3$nn.index,trees.nn3$nn.dist))
    names(trees.nn3)=c("NN1", "NN2", "NN3", "NN1Dist","NN2Dist", "NN3Dist")
    
    
    trees.nnz=matrix(0, nrow(trees.df),3)
    trees.nnx=matrix(0, nrow(trees.df),3)
    trees.nny=matrix(0, nrow(trees.df),3)
    trees.nnindex=matrix(0, nrow(trees.df),3)
    trees.nndist=matrix(0, nrow(trees.df),3)
    
    #for each tree compute the xyz, the index and the distance to nearest 3 neighbors
    for (j in 1:nrow(trees.df)){
      trees.nnz[j,]=trees.df$Z[as.numeric(trees.nn3[j,c("NN1","NN2","NN3")])]
      trees.nnx[j,]=trees.df$X[as.numeric(trees.nn3[j,c("NN1","NN2","NN3")])]
      trees.nny[j,]=trees.df$Y[as.numeric(trees.nn3[j,c("NN1","NN2","NN3")])]
      trees.nnindex[j,]=as.numeric(trees.nn3[j,c("NN1","NN2","NN3")])
      trees.nndist[j,]=as.numeric(trees.nn3[j,c("NN1Dist","NN2Dist","NN3Dist")])
    }
    
    
    #rename the columns within the matrix
    colnames(trees.nnz)=c("Z1","Z2","Z3")
    colnames(trees.nnx)=c("X1","X2","X3")
    colnames(trees.nny)=c("Y1","Y2","Y3")
    colnames(trees.nnindex)=c("NN1","NN2","NN3")
    colnames(trees.nndist)=c("NN1Dist","NN2Dist","NN3Dist")
    
    trees.nn=cbind(trees.nnindex, trees.nndist, trees.nnz, trees.nnx, trees.nny, trees.df)
    
    
    min.d2nt=6.5
    trees.nn$FP=0
    for (j in 1:nrow(trees.nn)){
      #j=2
      tmp=trees.nn[j,c("NN1Dist", "NN2Dist","NN3Dist", "Z1", "Z2","Z3")]
      z.tmp=tmp[4:6][which(tmp[1:3]<=min.d2nt)]
      if (length(z.tmp)==0) {next}
      if(trees.nn$Z[j] <= max(z.tmp)) {
        trees.nn$FP[j] = 1
      }
    }
    
    trees.final = trees.nn[which(trees.nn$FP==0 & trees.nn$area>20),]
    trees.final$StandID=i
    
    time.nn.e = Sys.time()
    print(paste0("Neighbors Elimination time: ", difftime(time.nn.e, time.nn.s, units="secs")))
    
    #sum(trees.nn$NN1Dist<min.d2nt)
    #sum(trees.nn$FP)
    
    #create the final tree measurements file
    trees.shp=trees.nn
    trees.shp=subset(trees.shp,trees.shp$FP==0)
    coordinates(trees.shp)=~X+Y
    trees.shp=st_as_sf(trees.shp)
    st_crs(trees.shp)=st_crs(itc)
    trees.shp=st_transform(trees.shp, 6557)
    #crs(trees.shp)
    poly.inv.f=st_intersection(trees.shp, 
                               st_geometry(poly.inv0))
    poly.inv.f$StandID=poly.inv0$StandID
    #plot(st_geometry(trees.shp), add=TRUE)
    #plot(st_geometry(plot.inv.f), add=TRUE)
    
    st_write(poly.inv.f, paste0(getwd(),"/Results/TSegDalp_StandID_", poly.inv0$StandID ,".shp"), append=FALSE)#write the segmented trees as a shapefile  
    
    diff.Lidar.Inv =as.data.frame(cbind(poly.inv0$StandID, 
                                        nrow(poly.inv.f), poly.inv0$T4I,
                                        floor(nrow(poly.inv.f)- poly.inv0$T4I),
                                        -1+nrow(poly.inv.f) / poly.inv0$T4I))
    names(diff.Lidar.Inv)=c("StandID","T4Lidar", "T4Inv", "Diff.Lidar.Inv.abs", "Diff.Lidar.Inv.rel")                         
    assess.dif=rbind(assess.dif, diff.Lidar.Inv)
    
  }
  
}





