#################################################################################
######################### CREATE UK_NE SURFACE IN HPC ###########################
#################################################################################
library(readxl)
library(mgcv)
library(raster)
library(rgdal)
library(ncdf4)
library(lubridate)

############################ 1. FIT NO2 GAM MODEL #############################
load("/rds/general/user/ww3414/ephemeral/RWorkspace/GAMmodel.RData")

############################ 2. EXTRACT POINT #############################
#### 2.0 IMPORT RASTER #####
roadlenmaj_1000 <- raster("2.raster/roadlenmaj_1000.tif")
trafload_50 <- raster("2.raster/trafload_50.tif")
ldurban_5000 <- raster("2.raster/ldurban_5000.tif")
grespace_500 <- raster("2.raster/grespace_500.tif")
natural_5000 <- raster("2.raster/natural_5000.tif")

##### 2.1 create newsites dataset #####
# London
pts <- read.csv("3.fishnet25_LONDON.csv");pts[1] <- list(NULL)

pts$ID <- 1:nrow(pts)
# Add lat-long to pts
pts.xy <- cbind(Easting = as.numeric(as.character(pts$POINT_X)),Northing = as.numeric(as.character(pts$POINT_Y))) # Create coordinates variable
colnames(pts.xy) <- c("X","Y")
sp.xy <- SpatialPointsDataFrame(pts.xy, data=data.frame(pts$ID,pts$ID), proj4string = CRS("+init=epsg:27700")) # Create the SpatialPointsDataFrame
sp.latlon <- spTransform(sp.xy,crs("+init=epsg:4326"))
colnames(sp.latlon@coords)[colnames(sp.latlon@coords) == "Easting"] <- "lon"
colnames(sp.latlon@coords)[colnames(sp.latlon@coords) == "Northing"] <- "lat"
pts.latlon <- data.frame(sp.latlon@coords)
colnames(pts.latlon) <- c("lon","lat")
pts.xy <- as.data.frame(pts.xy)

roadlenmaj_1000_pt <- extract(roadlenmaj_1000,pts.xy)
trafload_50_pt <- extract(trafload_50,pts.xy)
ldurban_5000_pt <- extract(ldurban_5000,pts.xy)
grespace_500_pt <- extract(grespace_500,pts.xy)
natural_5000_pt <- extract(natural_5000,pts.xy)

siteID <- c(1:nrow(pts.xy))
newsites <- cbind.data.frame(siteID,
                             pts.xy$X,
                             pts.xy$Y,
                             roadlenmaj_1000_pt,
                             trafload_50_pt,
                             ldurban_5000_pt,
                             grespace_500_pt,
                             natural_5000_pt)
names(newsites) <- c("siteID","X","Y","roadlenmaj_1000","trafload_50","ldurban_5000","grespace_500","natural_5000")

rm(grespace_500,grespace_500_pt,ldurban_5000,ldurban_5000_pt,natural_5000,natural_5000_pt,roadlenmaj_1000,
   roadlenmaj_1000_pt,trafload_50,trafload_50_pt,siteID)
gc()
##### 2.2 Create newdaily dataset #####
# Extract daily MACC values
# Command brick reads all layers (time slices) in the file
no22015 <- brick("4.ENSa.2015.NO2.yearlyrea.daymean.nc", varname = "NO2")
no22015 <- crop(no22015, extent(-9, 2.1, 49.5, 61)) # UK extent

# Get date index from the file
idxs <- as.integer(getZ(no22015))

# Subset the data - one month (try one month first)
no22015lon1 <- no22015[[1:31]]
no22015lon2 <- no22015[[32:59]]
no22015lon3 <- no22015[[60:90]]
no22015lon4 <- no22015[[91:120]]
no22015lon5 <- no22015[[121:151]]
no22015lon6 <- no22015[[152:181]]
no22015lon7 <- no22015[[182:212]]
no22015lon8 <- no22015[[213:243]]
no22015lon9 <- no22015[[244:273]]
no22015lon10 <- no22015[[274:304]]
no22015lon11 <- no22015[[305:334]]
no22015lon12 <- no22015[[335:365]]

# Put coordinates and extract values
a <- pts.latlon$lat
b <- pts.latlon$lon
c <- c(b,a); rm(a,b)

location <- matrix(c, ncol = 2)

vars1 <- extract(no22015lon1, location)
vars2 <- extract(no22015lon2, location)
vars3 <- extract(no22015lon3, location)
vars4 <- extract(no22015lon4, location)
vars5 <- extract(no22015lon5, location)
vars6 <- extract(no22015lon6, location)
vars7 <- extract(no22015lon7, location)
vars8 <- extract(no22015lon8, location)
vars9 <- extract(no22015lon9, location)
vars10 <- extract(no22015lon10, location)
vars11 <- extract(no22015lon11, location)
vars12 <- extract(no22015lon12, location)

vars <- cbind(vars1,vars2,vars3,vars4,vars5,vars6,vars7,vars8,vars9,vars10,vars11,vars12)
rm(no22015,no22015lon1,no22015lon2,no22015lon3,no22015lon4,no22015lon5,no22015lon6,no22015lon7,
   no22015lon8,no22015lon9,no22015lon10,no22015lon11,no22015lon12,
   vars1,vars2,vars3,vars4,vars5,vars6,vars7,vars8,vars9,vars10,vars11,vars12); gc()

# Merge dates and values and fix data frame names
siteID <- c(1:nrow(pts.xy))
df.no22015 <- cbind(siteID,vars)
rm(vars);gc()
df.no22015 <- data.frame(df.no22015)
names(df.no22015)  <- c("siteID", idxs)

ndays <- 365 # user choose: depending on how many months are included!!!
idxs <- as.Date(as.character(idxs), "%Y%m%d")
c1 <- rep(siteID, each = ndays) #siteID
c2 <- rep(pts.xy$X, each = ndays) #X
c3 <- rep(pts.xy$Y, each = ndays)#Y
#c4 <- rep(pts.latlon$lat, each = ndays) #lat
#c5 <- rep(pts.latlon$lon, each = ndays) #lon
c6 <- rep(idxs, nrow(pts.xy)) #date
c7 <- yday(c6) #dayID
c8 <- rep("2015", ndays*nrow(pts.xy)) #year
c9 <- c(t(df.no22015[2:ncol(df.no22015)])) #MACC

newdaily <- data.frame(c1,c2,c3,c6,c7,c8,c9)

names(newdaily) <- c("siteID","X","Y","date","dayID","year","MACC")
rm(c1,c2,c3,c6,c7,c8,c9)
#gc()

save.image("/rds/general/user/ww3414/ephemeral/RWorkspace/SurfaceDev_HPC_London_PredReady.RData")

############################ 3. MODEL PREDICTION #############################
# set up new sites 
valid <- newdaily
validSites <- newsites

# free memory space
rm(newdaily, newsites,d,daily,df.no22015,location,sites,
   subset,trainSites,trainSites2,trainSitesTmp,c,nTrainSites,
   no22015,uPredBySite,uSEBySite,uSites)
gc()
###

valid$sValidPred<-rep(0,nrow(valid)) # set up storage for validation predictions and prediction standard errors
valid$sValidSE<-rep(0,nrow(valid))

for (day in c(1:365)){ #use 31 for now, should be 1:365
  tmp<-predict(sMod[[day]],newdata=valid[valid$dayID==day,],se.fit=TRUE)
  valid$sValidSE[valid$dayID==day]<-tmp$se.fit  # prediction standard errors for spatial smooth component, validation
} 

# prediction and standard errors for second stage
tmp<-predict(uMod,newdata=validSites,se.fit=TRUE)
uValidPredBySite<-tmp$fit
uValidSEBySite<-tmp$se.fit
uValidTmp<-as.data.frame(as.factor(validSites$siteID))
uValidTmp<-cbind(uValidTmp,uValidPredBySite,uValidSEBySite)
names(uValidTmp)<-c("siteID","uValidPred","uValidSE")
valid<-merge(valid,uValidTmp,by.x="siteID",by.y="siteID",sort=FALSE)
#uvalidtmp.tmp <- as.data.frame(cbind(rep(uValidTmp$siteID,each=ndays),rep(uValidTmp$uValidPred,each=59),rep(uValidTmp$uValidSE,each=59))) 
#names(uvalidtmp.tmp) <- c("siteID","uValidPred","uValidSE")
#valid <- as.data.frame(cbind(valid, uvalidtmp.tmp))
#valid[,59] <- NULL #the four steps above achieve the same result as the merge step, but quicker
valid$validSE<-NA

# prediction and standard errors for first stage covariate smooths
valid$siteFactor<-train$siteFactor[1] # just need a dummy value
tmp<-predict(gMod,newdata=valid,type="terms",se.fit=TRUE)
if(ncol(tmp$fit)>2){
  valid$gValidPred<-apply(tmp$fit[,(2:ncol(tmp$fit))],1,sum)
  valid$gValidVar<-apply(tmp$se.fit[,(2:ncol(tmp$se.fit))]^2,1,sum)
} else{
  valid$gValidPred<-tmp$fit[,2]
  valid$gValidVar<-tmp$se.fit[,2]^2
}

# total prediction, summed over components
valid$validPred<-valid$gValidPred+valid$sValidPred+valid$uValidPred

# total standard error from summing variance terms over components
for(day in days){   # day-specific variance for first stage spatial smooths
  valid$validSE[valid$dayID==day]<-valid$sValidSE[valid$dayID==day]^2+sModSig2[day]  # actually on var scale for the moment
}
valid$validSE<-sqrt(valid$validSE+uMod$sig2+valid$uValidSE^2+valid$gValidVar)

# standard error components for each stage: second stage components are common to all months for a given site
valid$validSEstage1<-NA
for(day in days){
  valid$validSEstage1[valid$dayID==day]<-valid$sValidSE[valid$dayID==day]^2+sModSig2[day]  # actually on var scale for the moment
}
valid$validSEstage1<-sqrt(valid$validSEstage1+valid$gValidVar)
valid$validSEstage2<-sqrt(uMod$sig2+valid$uValidSE^2)

############################ 4. MODEL OUTPUT #############################
write.csv(valid,"fullValid_london.csv")
gamPred_LON <- cbind.data.frame(valid$X, valid$Y, valid$dayID, valid$validPred)

save.image("/rds/general/user/ww3414/ephemeral/RWorkspace/SurfaceDev_HPC_London_Result.RData")

# Subset valid by day
# valid1 <- valid[which(valid$dayID==1),]
# valid2 <- valid[which(valid$dayID==2),]

# Create netCDF from valid
