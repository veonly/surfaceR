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
load("/rds/general/user/ww3414/ephemeral/RWorkspace/SurfaceDev_HPC_London_PredReady.RData")

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
