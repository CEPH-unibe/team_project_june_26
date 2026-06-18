# Title: 2-stage analyses exposure-response function for daily temperature-related mortality 2014-2014 -------
# Group 1: Coral Salvador, Apolline Saucy, Meret Haldemann, Jacopo Vanoli
# Date created: 16.06.2026
# last update: 18.06.2026

# Prepare environment -----------------------------------------------------------
rm(list=ls(all=TRUE))
library(data.table) # HANDLE LARGE DATASETS
library(dlnm) ; library(gnm) ; library(splines) # MODELLING TOOLS
library(sf) ; library(terra) # HANDLE SPATIAL DATA
library(exactextractr) # FAST EXTRACTION OF AREA-WEIGHTED RASTER CELLS
library(dplyr) ; library(tidyr) # DATA MANAGEMENT TOOLS
library(ggplot2) ; library(patchwork) # PLOTTING TOOLS;
library(mixmeta)

## Load data ---------------------------------------------------------------------

### Load lookup table for communities and districts
lookup <- fread("/Volumes/FS/_ISPM/CCH/AnnualTeamProject2026/Boundaries_and_shapefiles/Gemeindestand_lookup_districts.csv")
listdistrictname <- unique(lookup$`Bezirks-nummer`) ## there are 144 districts
head(lookup)

### Load mortality data and limit to 2014 onward
mort <- readRDS("/Volumes/FS/_ISPM/CCH/AnnualTeamProject2026/death data/death6924.RDS")
mort.14 <- mort%>%
  filter(year(date)>=2014)
summary(mort.14)
length(unique(mort.14$muncode)) # 2246 unique communities

### Load temperature data and limit to 2014 onward
temp <- fread( "/Volumes/FS/_ISPM/CCH/AnnualTeamProject2026/Historical_temp_data/historical_temp_popw_muni_2000_2024.csv")
temp14 <- temp%>%filter(year(time)>=2014)
temp14 <- temp14%>%
  rename("date"=time,
         "muncode" = GDENR,
         "tmean" = mean_value)

### Load min mortality function
source('01_exposure_response_functions/findmin.R')

## Prepare dataset for analyses ---------------------------------------------------------
mort.14.m <- merge(mort.14, temp14, by=c("muncode", "date"), all.x=T)

summary(mort.14.m[is.na(tmean)])
table(mort.14.m[is.na(tmean), .(muncode)])
nrow(unique(mort.14.m[is.na(tmean), .(muncode)]))  ## 141 Municipalities are dropped from the analyses
                                                   ## (for different time extents)
                                                   ## because the codes do not match the most recent definition (2026 shapefile)
### Visualize
mort.14[muncode%in%c(unique(mort.14.m[is.na(tmean)]$muncode))]%>%
  group_by(muncode)%>%
  summarize(mindate = min(date),
            meandate = mean(date),
            mediandate= median(date),
            maxdate= max(date))%>%
  ggplot(aes(y=muncode))+
  geom_segment(aes(x = mindate, xend=maxdate))  ## they are missing for the full extent,
                                                ##but their data could be available for future years if covered by a new community code that did stay until 2026

mort.14 <- mort.14.m ## For now, we ignore this problem and drop these observations

mort.14 <- mort.14%>%arrange(muncode, date) ## Order dataset
mort.14$dow <- wday(mort.14$date) ## add day of the week


### Remove useless data
rm(temp, mort, mort.14.m)

# Analyses ----------------------------------------------------------------------------

## Model specifications --------------------------------------

### Specification of the exposure response function
varfun <- "ns"
vardegree <- NULL
varper <- c(50,90)

### Specification of the lag function
lag <- 3
lagnk <- 2

ncoef <- length(varper) + ifelse(varfun=="bs",vardegree,1)  # for 2 knots, these are 3 coefficients
# coeff <- array(NA,dim=list(ncoef,length(listdistrictname)),
#                dimnames=list(paste0("coeff", 1:ncoef)
#                              ,listdistrictname))
# coeff <- t(coeff)
# vcovv <- array(NA,dim=list(ncoef,ncoef,length(listdistrictname)),
#                dimnames=list(paste0("coeffh", 1:ncoef),paste0("coeffv", 1:ncoef)
#                              ,listdistrictname))
# head(vcovv)
#avertmean_district <- matrix(NA, nrow=length(listdistrictname), ncol=1)
#rangetmean_district <- matrix(NA, nrow=length(listdistrictname), ncol=1)

# #These are the percentile distribution for which we want to predict in the second stage analyses
# predper <- c(seq(0,1,0.1),2:98,seq(99,100,0.1))
# average_dist_district <- matrix(NA, nrow=length(predper), ncol=length(listdistrictname))

### List municipalities by district element
dlist<- split(lookup$`BFS Gde-nummer`, lookup$`Bezirks-nummer`)

### Empty vectors and matrices to store results for district-specific models
firststage<-list()
cp_list <- list()
coefall <- matrix(NA, nrow=length(dlist), ncol=ncoef)
vcovall <- list()


## Run first stage (small area analysis at the community level, aggregate by district) ---------------------------

### Loop through districts
for (i in seq_along(dlist)) {
  dist <- dlist[[i]]
  datafull<-mort.14[muncode%in%c(dlist[[i]]),]
  if(nrow(datafull>=1)){
  datafull$year <- year(datafull$date)
  datafull$month <- month(datafull$date)
  datafull <- subset(datafull, month%in%6:9)
  datafull$doy <- yday(datafull$date)


  # DEFINE SPLINES OF DAY OF THE YEAR
spldoy <- onebasis(datafull$doy, "ns", df=3)

argvar <- list(fun="ns", knots=quantile(datafull$tmean, varper/100, na.rm=T),
               Boundary.knots=range(datafull$tmean)) #Ana suggested inckuding two knots
# arglag <- list(fun="ns", knots=lagnk) # to discuss this but I think that it is reasonable (or maybe we can use the strata fuction)
arglag <- list(fun="integer") # to discuss this but I think that it is reasonable (or maybe we can use the strata fuction)

datafull$group <- factor(paste(datafull$muncode, datafull$year, sep="-"))
group <- factor(paste(datafull$muncode, datafull$year, sep="-"))
group <- with(datafull, factor(paste(muncode, year, sep="-")))
cbtmean <- crossbasis(datafull$tmean, lag=lag, argvar=argvar, arglag=arglag, #Check how 7 lags look, and if the results are too noisy/inaccurate, reduce them (e.g., 4 lags)
                      group=group)

# DEFINE THE STRATA
datafull[, stratum:=factor(paste(muncode, year, month, dow, sep=":"))]

# RUN THE MODEL
# NB: EXCLUDE EMPTY STRATA, OTHERWISE BIAS IN gnm WITH quasipoisson # we can control the seasonal patterns using the approach from Gasparrini et al--> https://github.com/gasparrini/CTS-smallarea
datafull[,  keep:=sum(dcount)>0, by=stratum]
modfull <- gnm(dcount ~ cbtmean ,
               eliminate=stratum, data=datafull, family=quasipoisson, subset=keep)

mmti<-findmin(cbtmean, model=modfull, from=quantile(datafull$tmean, 0.25), to=quantile(datafull$tmean, 0.90)) #check the function

cp_list[[i]] <- crossreduce(cbtmean, modfull, cen=mmti)

coefall[i,] <-  cp_list[[i]]$coefficients
vcovall[[i]] <-  cp_list[[i]]$vcov
  }
  else{
    cp_list[[i]] <- NA

    coefall[i,] <-   NA
    vcovall[[i]] <-  NA
  }
}

names(cp_list) <- names(dlist)
rownames(coefall) <- names(dlist)
names(vcovall) <- names(dlist)

### Store coefficients and variance-covariance matrices

firststage <-list(coefall=coefall, vcovall=vcovall)

### Plot results first stage

pdf(paste0("01_exposure_response_functions/firststageplots_bydistricts_", lag, "_days.pdf"))
for(x in 1:(length(cp_list)-1)){
  plot(cp_list[[x]],
       main=paste("First stage",  names(cp_list)[x])
  )
}
dev.off()

### Save results first stage

saveRDS(cp_list, "01_exposure_response_functions/crosspreds_stage1.rds")
saveRDS(firststage, "01_exposure_response_functions/coeffs_vcov_stage1.rds")

## Run second stage analyses -----------------------------------------------------

### add district and canton information to temperature data
temp14.m <- merge(temp14, lookup%>%
                    select(Kanton, `BFS Gde-nummer`, `Bezirks-nummer`), by.x="muncode", by.y="BFS Gde-nummer")
summary(temp14.m)
temp14.m <- temp14.m[lubridate::month(date)%in%6:9,] ## restrict to summer months

# View(temp14.m[is.na(`Bezirks-nummer`)])

### Define fixed effects predictors
avgtmean <- temp14.m[, mean(tmean, na.rm=T), by=c("Bezirks-nummer", "Kanton")]
rangetmean <- temp14.m[, c("min","max"):= as.list(range(tmean, na.rm=T)), by=c("Bezirks-nummer", "Kanton")]
rangetmean <- unique(rangetmean%>%select(`Bezirks-nummer`, min, max))
rangetmean[, rangetmean:=(max-min)]

metavarALL <- cbind(rangetmean, avgtmean[,-1])
metavarALL <-
  metavarALL%>%
  rename("avgtmean" = V1)

# avgtmean   <- sapply(temp14.m,function(x) mean(x$tmean,na.rm=TRUE)) #average of mean temperature (ºC)
# rangetmean <- sapply(dlist,function(x) diff(range(x$tmean,na.rm=TRUE))) #range of mean temperature (ºC)
#
# metavarALL<-data.frame(avgtmean=avgtmean, rangetmean=rangetmean, district=district, district=district)
#
coefmeta <- coefall
vcovmeta <- vcovall

#### Run meta-analysis ----------------------------------------------------

mvall <- mixmeta(coefmeta~rangetmean+avgtmean,vcovmeta, metavarALL,
                 control=list(showiter=T), random=~1|`Kanton`/`Bezirks-nummer`, method="reml")
# mvall <- mixmeta::mixmeta(coefmeta~rangetmean+avgtmean,vcovmeta, metavarALL,
#                  control=list(showiter=T), random=~1|Kanton/`Bezirks-nummer`, method="reml")

# ## BLUPS AT district LEVEL FROM TWO-LEVEL MODEL
#
# districtblup <- exp(blup(mvall))
# rownames(districtblup) <- names(cp_list)[1:143]
# saveRDS(blup, "01_exposure_response_functions/secondstage.rds")

### Pooled estimates ----------------------------------------------------

temp14.m <- temp14.m[lubridate::month(date)%in%6:9,]
avgtmean <- temp14.m[, mean(tmean, na.rm=T),
                     by=c("Bezirks-nummer", "Kanton")]
avgrange_district <- temp14.m[, .(avgtmean = mean(tmean, na.rm = TRUE),
                                  rangetmean  = max(tmean, na.rm = TRUE)-min(tmean, na.rm = TRUE)),
                              by = .(`Bezirks-nummer`, Kanton)]
predper <- c(seq(0,1,0.1),2:98,seq(99,100,0.1))
qt_district <- temp14.m[, .(percentile = predper,
                            quantile   = quantile(tmean, probs = predper / 100)),
                        by = .(`Bezirks-nummer`)]
qt_district_wide <- dcast(data=qt_district, percentile ~ `Bezirks-nummer`,
                          value.var = "quantile")


meanaverage_dist_district <- rowMeans(qt_district_wide[,-1])
knots <- meanaverage_dist_district [predper %in% c(50,90)]
bvar <- onebasis(meanaverage_dist_district , fun="ns", knots=knots)


newdata <- data.frame(avgtmean=mean(avgrange_district$avgtmean),
                      rangetmean=mean(avgrange_district$rangetmean))

predfit <- predict(mvall, newdata=newdata, vcov=T)
predpool <- crosspred(bvar, coef=predfit$fit, vcov=predfit$vcov,
                      model.link="log", at=meanaverage_dist_district, cen=16) ## I changed to 16 (visual minimum)

#### Plot

pdf("01_exposure_response_functions/secondstage_pooled.pdf",width=9,height=7)
plot(predpool)
dev.off()

### BLUPs -----------------------------------------------

blup <- mixmeta::blup(mvall, vcov=T)
qt_district_wide <- qt_district_wide[,-1]

pdf("01_exposure_response_functions/blups_districts.pdf",width=9,height=13)
layout(matrix(seq(5*3),nrow=5,byrow=T))
par(mar=c(4,3.8,3,2.4),mgp=c(2.5,1,0),las=1)


mintempdistrict <- rep(NA,143)
for(i in 1:143) {

  predvar <- qt_district_wide[[i]]

  # REDEFINE THE FUNCTION USING ALL THE ARGUMENTS (BOUNDARY KNOTS INCLUDED)
  argvar <- list(x=predvar,fun=varfun,
                 knots=qt_district_wide[predper %in% varper,i, with = FALSE][[1]],
                 Bound=qt_district_wide[c(1,length(predper)),i, with = FALSE][[1]])
  if(!is.null(vardegree)) argvar$degree <- vardegree
  bvar <- do.call(dlnm::onebasis,argvar)

  ### ERC PLOTS
  minperdistrict <- (50:90)[which.min((bvar%*%blup[[i]]$blup)[50:90])]
  mintempdistrict[i] <- qt_district_wide[predper==minperdistrict,i, with = FALSE]
  predblup <- crosspred(bvar,coef=blup[[i]]$blup,vcov=blup[[i]]$vcov,
                                 model.link="log",by=0.1,cen=mintempdistrict[i])

  plot(predblup,ylim=c(0,2.5),lab=c(6,5,7),xlab="Mean Temperature",
       ylab="Relative Risk", main=names(cp_list)[i], col="red",
       ci.arg=list(col=alpha("red",0.1)),lwd=1.5)

}
dev.off()


