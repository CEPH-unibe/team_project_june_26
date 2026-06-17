
################################################################################
# MAIN MODELS
################################################################################
rm(list=ls(all=TRUE))
library(data.table) # HANDLE LARGE DATASETS
library(dlnm) ; library(gnm) ; library(splines) # MODELLING TOOLS
library(sf) ; library(terra) # HANDLE SPATIAL DATA
library(exactextractr) # FAST EXTRACTION OF AREA-WEIGHTED RASTER CELLS
library(dplyr) ; library(tidyr) # DATA MANAGEMENT TOOLS
library(ggplot2) ; library(patchwork) # PLOTTING TOOLS
################################################################################
### save environment with everything
# SET DIRECTORY OF DATA INPUT

dir <- "/Volumes/FS/_ISPM/CCH/AnnualTeamProject2026/death data"
dirout <- "/Volumes/FS/_ISPM/CCH/AnnualTeamProject2026/01_exposure_response_functions/"

lookup <- fread("/Volumes/FS/_ISPM/CCH/AnnualTeamProject2026/Boundaries_and_shapefiles/Gemeindestand_lookup_districts.csv")

listdistrictname <- unique(lookup$Bezirksname)
head(lookup)

# SPECIFICATION OF THE EXPOSURE FUNCTION
varfun <- "ns"
vardegree <- NULL
varper <- c(50,90)

# SPECIFICATION OF THE LAG FUNCTION
lag <- 5
lagnk <- 2

ncoef <- length(varper) + ifelse(varfun=="bs",vardegree,1)
coeff <- array(NA,dim=list(ncoef,length(listdistrictname)),
               dimnames=list(paste0("coeff", 1:3)
                             ,listdistrictname))
coeff <- t(coeff)
vcovv <- array(NA,dim=list(ncoef,ncoef,length(listdistrictname)),
               dimnames=list(paste0("coeffh", 1:3),paste0("coeffv", 1:3)
                             ,listdistrictname))
head(vcovv)

#avertmean_district <- matrix(NA, nrow=length(listdistrictname), ncol=1)
#rangetmean_district <- matrix(NA, nrow=length(listdistrictname), ncol=1)

predper <- c(seq(0,1,0.1),2:98,seq(99,100,0.1))
average_dist_district <- matrix(NA, nrow=length(predper), ncol=length(listdistrictname))

# Load min mortality function
source("01_exposure_response_functions/findmin.R")

#each element include municipalities within each district
dlist<- split(lookup$`BFS Gde-nummer`, lookup$`Bezirks-nummer`)


firststage<-list()
cp_list <- list()
coefall <- matrix(NA, nrow=length(dlist), ncol=9)
vcovall <- list()



mort <- readRDS("/Volumes/FS/_ISPM/CCH/AnnualTeamProject2026/death data/death6924.RDS")
mort.14 <- mort%>%
  filter(year(date)>=2014)
summary(mort.14)

temp <- fread( "/Volumes/FS/_ISPM/CCH/AnnualTeamProject2026/Historical_temp_data/historical_temp_popw_muni_2000_2024.csv")
temp14 <- temp%>%filter(year(time)>=2014)
temp14 <- temp14%>%
  rename("date"=time,
         "muncode" = GDENR,
         "tmean" = mean_value)
# mort.14$tmean <- rnorm(nrow(mort.14))

mort.14.m <- merge(mort.14, temp14, by=c("muncode", "date"), all.x=T)

summary(mort.14.m[is.na(tmean)])
View(mort.14.m[is.na(tmean)])
table(mort.14.m[is.na(tmean), .(muncode)])
nrow(unique(mort.14.m[is.na(tmean), .(muncode)]))

mort.14 <- mort.14.m

setkey(mort.14, muncode, date)
mort.14 <- mort.14%>%arrange(muncode, date)
mort.14$dow <- wday(mort.14$date)

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

argvar <- list(fun="ns", knots=quantile(datafull$tmean, c(50,90)/100, na.rm=T),
               Boundary.knots=range(datafull$tmean)) #Ana suggested inckuding two knots
arglag <- list(fun="ns", knots=2) # to discuss this but I think that it is reasonable (or maybe we can use the strata fuction)

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

cp_list[[i]] <- crosspred(cbtmean, modfull, cen=mmti)

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

#Store coefficients and variance-covariance matrices
firststage <-list(coefall=coefall, vcovall=vcovall)

pdf(paste0("01_exposure_response_functions/firststageplots_bydistricts.pdf"))
for(x in 1:length(cp_list)){
  plot(cp_list[[x]], "overall")
}
dev.off()

saveRDS(cp_list, "01_exposure_response_functions/crosspreds_stage1.rds")
saveRDS(firststage, "01_exposure_response_functions/coeffs_vcov_stage1.rds")

# second stage analyses

# add district numbers to temperature data
temp14.m <- merge(temp14, lookup%>%
                    select(Kanton, `BFS Gde-nummer`, `Bezirks-nummer`), by.x="muncode", by.y="BFS Gde-nummer")
# summary(temp14.m)
# View(temp14.m[is.na(`Bezirks-nummer`)])
#
# #predictors
temp14.m <- temp14.m[month(date)%in%6:9,]
avgtmean <- temp14.m[, mean(tmean, na.rm=T), by=c("Bezirks-nummer", "Kanton")]
rangetmean <- temp14.m[, c("min","max"):= as.list(range(tmean, na.rm=T)), by=c("Bezirks-nummer", "Kanton")]


rangetmean <- unique(rangetmean%>%select(`Bezirks-nummer`, min, max))
rangetmean[, rangetmean:=(max-min)]

metavarALL <- cbind(rangetmean, avgtmean[,-1])
metavarALL$avgtmean <- metavarALL$V1
# avgtmean   <- sapply(temp14.m,function(x) mean(x$tmean,na.rm=TRUE)) #average of mean temperature (ºC)
# rangetmean <- sapply(dlist,function(x) diff(range(x$tmean,na.rm=TRUE))) #range of mean temperature (ºC)
#
# metavarALL<-data.frame(avgtmean=avgtmean, rangetmean=rangetmean, district=district, district=district)
#
coefmeta <- coefall
vcovmeta <- vcovall

# mvall <- mixmeta(coefmeta~rangetmean+avgtmean,vcovmeta, metavarALL,
#                  control=list(showiter=T), random=~1|`Kanton`/`Bezirks-nummer`, method="reml")
mvall <- mixmeta(coefmeta~rangetmean+avgtmean,vcovmeta, metavarALL,
                 control=list(showiter=T), random=~1|`Bezirks-nummer`, method="reml")

# blup(mvall, se=T, pi=T)
#
# # BLUPS AT district LEVEL FROM TWO-LEVEL MODEL
# districtblup <- exp(blup(mvall))

# rownames(districtblup) <- names(cp_list)[1:143]




saveRDS(mvall, "01_exposure_response_functions/secondstage.rds")



#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
