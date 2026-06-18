@ -1,47 +0,0 @@
  # add temperature data to mortality

  # First, prepare lookup table municipality to districts and verify shapefiles

  # libraries
  library(data.table)
library(dplyr)
library(sf)

# load data

# lookup table
lookup <- fread("/Volumes/FS/_ISPM/CCH/AnnualTeamProject2026/Boundaries_and_shapefiles/Gemeindestand_lookup_districts.csv")
summary(lookup)
sum(duplicated(lookup$Gemeindename))
sum(duplicated(lookup$`BFS Gde-nummer`))


#shapefiles
communes.shp <- st_read("/Volumes/FS/_ISPM/CCH/AnnualTeamProject2026/Boundaries_G1_Commune_20260101/Boundaries_G1_Commune_20260101.shp")
districts.shp <- st_read("/Volumes/FS/_ISPM/CCH/AnnualTeamProject2026/Boundaries_G1_District_20260101/Boundaries_G1_District_20260101.shp")

plot(districts.shp)

# mortality

mort <- readRDS("/Volumes/FS/_ISPM/CCH/AnnualTeamProject2026/death data/death6924.RDS")
mort.14 <- mort%>%
  filter(year(date)>=2014)
summary(mort.14)
# check overlap between commune shapefile and lookup table

match <- merge(communes.shp, lookup, by.x="GDENR", by.y="BFS Gde-nummer", all=T)
summary(match)

# check overlap lookup table and mortality data

match2 <- merge(mort.14, lookup, by.x="muncode", by.y="BFS Gde-nummer", all=T)
summary(match2)
View(match2[is.na(`Bezirks-nummer`)])

temp <- fread("/Volumes/FS/_ISPM/CCH/AnnualTeamProject2026/historical_temp_popw_2000_2024.csv")

# mortality data is already pop averaged at the district level

head(mort)
head(temp)
