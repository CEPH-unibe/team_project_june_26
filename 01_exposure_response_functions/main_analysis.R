
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
dirout <- "/Volumes/FS/_ISPM/CCH/AnnualTeamProject2026/exp_response/"

lookup <- fread("/Volumes/FS/_ISPM/CCH/AnnualTeamProject2026/Boundaries_and_shapefiles/Gemeindestand_lookup_districts.csv")

listdistrictname <- lookup$Bezirksname
head(lookup)

# SPECIFICATION OF THE EXPOSURE FUNCTION
varfun <- "ns"
vardegree <- NULL
varper <- c(50,90)

# SPECIFICATION OF THE LAG FUNCTION
lag <- 7
lagnk <- 2

ncoef <- length(varper) + ifelse(varfun=="bs",vardegree,1)
coeff <- array(NA,dim=list(ncoef,length(listdistrictname)),
               dimnames=list(paste0("coeff", 1:3)
                             ,listdistrictname))
vcovv <- array(NA,dim=list(ncoef,ncoef,length(listdistrictname)),
               dimnames=list(paste0("coeffh", 1:3),paste0("coeffv", 1:3)
                             ,listdistrictname))

#avertmean_district <- matrix(NA, nrow=length(listdistrictname), ncol=1)
#rangetmean_district <- matrix(NA, nrow=length(listdistrictname), ncol=1)

predper <- c(seq(0,1,0.1),2:98,seq(99,100,0.1))
average_dist_district <- matrix(NA, nrow=length(predper), ncol=length(listdistrictname))


# CREATE LIST WITH ONE DATASET PER DISTRICT
dlist_all <- list()

#each element include municipalities within each district
dlist<- split(lookup$muncode, lookup$BEZNR)

for (i in seq_along(dlist)) {
  datafull<-dlist[[i]]
  datafull$year <- year(datafull$date)
  datafull$month <- month(datafull$date)
  datafull <- subset(datafull, month%in%6:9)
  datafull$doy <- yday(datafull$date)
  
  
  # DEFINE SPLINES OF DAY OF THE YEAR
spldoy <- onebasis(datafull$doy, "ns", df=3)

argvar <- list(fun="ns", knots=quantile(datafull$tmean, c(50,90)/100, na.rm=T)) #Ana suggested inckuding two knots
arglag <- list(fun="ns", knots=2) # to discuss this but I think that it is reasonable (or maybe we can use the strata fuction)
group <- factor(paste(datafull$muncode, datafull$year, sep="-"))
cbtmean <- crossbasis(datafull$tmean, lag=7, argvar=argvar, arglag=arglag, #Check how 7 lags look, and if the results are too noisy/inaccurate, reduce them (e.g., 4 lags)
                      group=group)

# DEFINE THE STRATA 
datafull[, stratum:=factor(paste(BEZNR, year, month, dow, sep=":"))]

# RUN THE MODEL
# NB: EXCLUDE EMPTY STRATA, OTHERWISE BIAS IN gnm WITH quasipoisson # we can control the seasonal patterns using the approach from Gasparrini et al--> https://github.com/gasparrini/CTS-smallarea
datafull[,  keep:=sum(dtot)>0, by=stratum]
modfull <- gnm(dtot ~ cbtmean , 
               eliminate=stratum, data=datafull, family=quasipoisson, subset=keep)

mmti<-findmin(cbtmean, model=modfull, from=quantile(datafull$tmean, 0.25), to=quantile(datafull$tmean, 0.90)) #check the function

cp_list[[i]] <- crosspred(cbtmean, modfull, cen=mmti)

coefall[i,] <-  cp_list[[i]]$coefficients
vcovall[[i]] <-  cp_list[[i]]$vcov
}
  
#Store coefficients and variance-covariance matrices 
firststage <-list(coefall=coefall, vcovall=vcovall)

#predictors
avgtmean   <- sapply(dlist,function(x) mean(x$tmean_CRU,na.rm=TRUE)) #average of mean temperature (ºC)
rangetmean <- sapply(dlist,function(x) diff(range(x$tmean_CRU,na.rm=TRUE))) #range of mean temperature (ºC)

metavarALL<-data.frame(avgtmean=avgtmean, rangetmean=rangetmean, district=district, district=district)

coefmeta <- coefall
vcovmeta <- vcovall


















