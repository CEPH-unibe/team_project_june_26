
###################################################################################
#Preprocess mortality data into clean time series
###################################################################################


#Load libraries and functions
source("~/Library/CloudStorage/OneDrive-UniversitaetBern/my projects/team_project/team_project_june_26/01_exposure_response_functions/mortality/00.pkg.R")

#Define directories
deathdir <- "/Volumes/FS/_ISPM/CCH/01Data/Mortality_CH/"
savedir <- "~/Library/CloudStorage/OneDrive-UniversitaetBern/my projects/team_project/team_project_june_26/01_exposure_response_functions/mortality/"

#Load mortality data
death1 <- fread(paste0(deathdir,"/mortality_CH_1969-2018/mort_6918.csv"))
death2 <- fread(paste0(deathdir,"/mortality_CH_2017_2024/Datenlieferung_UNIBE_Vicedo_TU ab 2017_20251118.csv"))
death3 <- fread(paste0(deathdir,"/mortality_CH_2017_2024/Datenlieferung_UNIBE_Vicedo_TU 2024_20260105.csv"))

###PROCESS DEATHS IN THE PERIOD 1969-2018: DEATH1
#Subset relevant columns
death1 <- death1[,c("comm_resi","dod")]

#Rename columns
setnames(death1, new = c("muncode","date"), old = c("comm_resi","dod"))

#Reformat dates
death1[, date := as.IDate(date, format = "%d%b%Y")][, date := as.Date(date, format = "%Y-%m-%d")]

#Aggregate counts
death1 <- death1[, .(dcount = .N), by = .(muncode,  date)]

#Reorder by municipality code and date
setorder(death1, muncode, date)

################################################
###PROCESS DEATHS IN THE PERIOD 2012-2023, SUBSETTING TO 2019-2023 ONLY: DEATH2

death2[, date := as.Date(paste(EREIGNIS_JJJJ_GES_N, EREIGNIS_MM_GES_N, EREIGNIS_TT_GES_N, sep = "-"))]
death2[, muncode := WORT_AKT_GEM_N]

#Aggregate counts
death2 <- death2[, .(dcount = .N), by = .(muncode,  date)]

#Subset data from 2019 to 2023 
death2 <- death2[year(date)>=2019 & year(date)<=2023,]

#Reorder by municipality code and date
setorder(death2, muncode, date)


################################################
###PROCESS DEATHS IN THE PERIOD 2020-2024, SUBSETTING TO 2024 ONLY: DEATH3

death3[, date := as.Date(paste(EREIGNIS_JJJJ_GES_N, EREIGNIS_MM_GES_N, EREIGNIS_TT_GES_N, sep = "-"))]
death3[, muncode := WORT_AKT_GEM_N]

#Aggregate counts
death3 <- death3[, .(dcount = .N), by = .(muncode,  date)]

#Subset data from 2019 to 2023 
death3 <- death3[year(date)==2024,]

#Reorder by municipality code and date
setorder(death3, muncode, date)

################################################
###APPEND ALL PERIODS TO SINGLE PERIOD: 1969-2024
death <- bind_rows(death1,death2,death3)

################################################
###SAVE MORTALITY DATA
saveRDS(death, paste0(savedir,"death6924.RDS"))

