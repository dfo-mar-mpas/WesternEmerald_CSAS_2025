library(Mar.datawrangling)

get_data('rv', data.dir = "C:/Users/HarbinJ/Documents/data/rv") # Note you'll have to change your data.dir


GSDET$latitude <- 0
GSDET$longitude <- 0
GSDET$year <- 0
missions <- unique(GSDET$MISSION)

GSINF <-GSINF[-which(is.na(GSINF$SDATE)),]
for (i in seq_along(missions)) {
  GSDET$latitude[which(GSDET$MISSION == missions[i])] <- GSINF$LATITUDE[which(GSINF$MISSION == missions[i])][1]
  GSDET$longitude[which(GSDET$MISSION == missions[i])]  <- GSINF$LONGITUDE[which(GSINF$MISSION == missions[i])][1]
  GSDET$year[which(GSDET$MISSION == missions[i])]  <- unique(as.numeric(substr(GSINF$SDATE[which(GSINF$MISSION == missions[i])],1,4)))
}


gsdet <- GSDET
