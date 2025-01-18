#Code for formatting the RV data for comparison to the eDNA

#load libraries ----------
library(sf)
library(tidyverse)
library(rnaturalearth)
library(rnaturalearthhires)
library(lubridate)
library(taxize)
library(worrms)
library(MarConsNetData)

#load functions -----------
source("code/ClassifyFunction.R")

#Load projections
latlong <- "+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"
utm <- "+proj=utm +zone=19N +datum=WGS84 +units=m +no_defs"
utmkm <- "+proj=utm +zone=20 +datum=NAD83 +units=km +no_defs +ellps=GRS80 +towgs84=0,0,0"

# #load RV data
# get_data('rv',data.dir="R:/Science/CESD/HES_MPAGroup/Data/RVdata/RVdata")
# rvdat <- summarize_catches()

#Basemaps
bioregion <- data_bioregion()%>%
  st_transform(latlong)

#bounding box for the bioregion
bioregion_box <- bioregion%>%
  st_transform(utmkm)%>%
  st_buffer(50)%>% #50km buffer
  st_transform(latlong)%>%
  st_bbox()%>%st_as_sfc()
#plotlimits
plotlims <- bioregion%>%
  st_transform(utmkm)%>%
  st_buffer(20)%>% #50km buffer
  st_transform(latlong)%>%
  st_bbox()

#maritimes network -- note that this cannot be shared or made available online to the public
maritimes_network <- data_draft_areas()%>%
  st_transform(latlong)%>%
  st_make_valid()%>%
  dplyr::select(Classification_E,SiteName_E)%>%
  rename(status=Classification_E,name=SiteName_E)

#load the webca buffer
webca_buffer <- read_sf("data/Shapefiles/WEBCA_10k_85k.shp")%>%
                st_transform(latlong)%>%
                mutate(status="buffer",name="webca_buffer")%>%
                rename(geoms=geometry)%>%
                st_difference(maritimes_network%>%filter(name=="Western/Emerald Banks Marine Refuge"))%>%
                dplyr::select(status,name)


#bioclassification polygons
bioclass <- read_sf("data/Shapefiles/bioclassification_clusters.shp")%>%
  st_transform(latlong)

#Create basemap intersected with the bounding box. 
basemap_atlantic <- rbind(ne_states(country = "Canada",returnclass = "sf")%>%
                            dplyr::select(name_en,geometry)%>%
                            st_as_sf()%>%
                            st_union()%>%
                            st_transform(latlong)%>%
                            st_as_sf()%>%
                            mutate(country="Canada"),
                          ne_states(country = "United States of America",returnclass = "sf")%>%
                            dplyr::select(name_en,geometry)%>%
                            st_as_sf()%>%
                            st_union()%>%
                            st_transform(latlong)%>%
                            st_as_sf()%>%
                            mutate(country="USA"))%>%
  st_intersection(.,bioregion_box)

#Get RV data
# rvopen <- data_MarRV_survey_open()
# 
# #assign a common naming system (this will help to get around year conventions)
# gs_names <- c("GSCAT","GSINF","GSMISSIONS","GSSPECIES")
# 
# for(i in gs_names){assign(i,rvopen[[which(grepl(i,names(rvopen)))]])}
# 
# #Now identify a coordinate for each station. For some stations there is only the 'start' others have the 'end'
# #for those with a start and an end we can calculate the mid-point
# 
# GSINF <- GSINF%>%
#   rowwise()%>%
#   mutate(
#     # Check if both start and end points are available
#     MLAT = if (!is.na(SLAT) & !is.na(SLONG) & !is.na(ELAT) & !is.na(ELONG)) {
#       geosphere::midPoint(c(SLONG, SLAT), c(ELONG, ELAT))[2] # Extract latitude
#     } else if (!is.na(SLAT) & !is.na(SLONG)) {
#       SLAT # Only start points available
#     } else if (!is.na(ELAT) & !is.na(ELONG)) {
#       ELAT # Only end points available
#     } else {
#       NA_real_ # Neither points are available
#     },
#     MLONG = if (!is.na(SLAT) & !is.na(SLONG) & !is.na(ELAT) & !is.na(ELONG)) {
#       geosphere::midPoint(c(SLONG, SLAT), c(ELONG, ELAT))[1] # Extract longitude
#     } else if (!is.na(SLAT) & !is.na(SLONG)) {
#       SLONG # Only start points available
#     } else if (!is.na(ELAT) & !is.na(ELONG)) {
#       ELONG # Only end points available
#     } else {
#       NA_real_ # Neither points are available
#     },
#     ID = case_when(
#       !is.na(SLAT) & !is.na(ELAT) ~ "Middle", # Both start and end points available
#       !is.na(SLAT) ~ "Start", # Only start points are available
#       !is.na(ELAT) ~ "End", # Only end points are available
#       TRUE ~ NA_character_ # If neither are available
#     )
#   ) %>%
#   ungroup()
# 
# rvdata <- GSCAT%>% #catch data
#           left_join(.,GSMISSIONS)%>% #general info on the cruises
#           left_join(.,GSSPECIES%>%rename(SPEC = CODE))%>% #species IDs
#           left_join(.,GSINF)%>% #set specific data
#           mutate(std_count = TOTNO*1.75/DIST, #standardized to a set tow distance of 1.75 nm
#                  std_wgt = TOTWGT*1.75/DIST,
#                  id=paste(MISSION,SETNO,sep="-"))%>%
#           st_as_sf(coords=c("MLONG","MLAT"),crs=latlong)  #convert to sf object
            
#save(rvdata,file="data/rvdata.RData")\
#save(GSCAT,GSINF,GSMISSIONS,GSSPECIES,file="data/RVSurveyPull.RData")

load("data/rvdata.RData") #only need to run the above commmented code once. 
load("data/RVSurveyPull.RData")

rv_stations <- GSINF %>%
              mutate(date = as.POSIXct(SDATE), year = year(date), decade = paste0(floor(year / 10) * 10, "'s")) %>%
              st_as_sf(coords = c("MLONG", "MLAT"), crs = latlong) %>%
              st_join(., maritimes_network%>%rbind(webca_buffer), join = st_intersects)%>%
              filter(name %in% c("webca_buffer","Western/Emerald Banks Marine Refuge"))%>%
              mutate(id=paste(MISSION,SETNO,sep="-"))

rv_df <- rvdata%>%
         filter(id %in% rv_stations$id)%>%
         left_join(.,rv_stations%>%data.frame()%>%dplyr::select(id,name))

#id the species captured in the RV survey within the WEBCA + buffer

webca_species_rv <- rv_df%>%
                 filter(name == "Western/Emerald Banks Marine Refuge")%>%
                 data.frame()%>%
                 dplyr::select(SCI_NAME)%>%
                 distinct(SCI_NAME)%>%
                 mutate(type="webca")%>%
                 rename(name=1)

webca_species_buffer <- rv_df%>%
                        filter(name == "webca_buffer")%>%
                        data.frame()%>%
                        pull(SCI_NAME)%>%
                        unique()%>%
                        setdiff(.,webca_species)%>% #the additional species
                        data.frame()%>%
                        rename(name=1)%>%
                        mutate(type="buffer")

webca_species <- rbind(webca_species_rv,webca_species_buffer)

#Classify the species now

id_names <-webca_species%>%
              mutate(latin = name, #clean up the names for entry into various online taxonomic databases via APIs and the taxise package.#basic data cleaning
                     latin = gsub("SPP.","G.",latin,fixed=T),#this specifies that these are to the genus level
                     latin = gsub("SP.","G.",latin,fixed=T),
                     latin = gsub(" SP","G.",latin,fixed=T),
                     latin = gsub("S.P.","G.",latin,fixed=T),
                     latin = gsub("S.C.","",latin,fixed=T), #Not sure what is meant here (superclass?) but we can ID based on the db's to class.

                     latin = gsub("COLUS G.","BUCCINIDAE F.",latin), #COLUS SP are a genus of Buccinidae but this doesn't return anything in itis. This (like Hayus) will have to modified after to add the family name.
                     latin = gsub("HYAS G.","OREGONIIDAE F.",latin),#these are specified as toad crabs (genus hyas)# this will fix the gsub error caused by the first hyas sub but for some reason the online db's return plants. Have to add 'hyas' to the tax id manually

                     latin = gsub("BRYOZOANS BRACHIOPODA P.","BRACHIOPODA P.",latin,fixed=T), #for some reason COMM 'Lampshells' are listed as under a grouped phylum id. In this case we expect them to be part of the brachiopoda and not bryozoa
                     latin = gsub("OREGONIIDAE F. ARANEUS","HYAS ARANEUS",latin),# this will fix the gsub error caused by the first hyas sub
                     latin = gsub("OREGONIIDAE F. COARCTATUS","HYAS COARCTATUS",latin),# this will fix the gsub error caused by the first hyas sub
                     latin = gsub("PENNATULACEA","PENNATULACEA O.",latin), #this will specify that this is to the order level.
                     latin = gsub("GORGONOCEPHALIDAE,ASTERONYCHIDAE F.","PHRYNOPHIURIDA O.",latin), #there are two families grouped for basket stars, so we will go up one taxonomic level to order Phrynophiurida
                     latin = ifelse(grepl("GEPHYREA",latin),"SIPUNCULA P.",latin), #Gephyrea is a former taxon group which now is split into three phyla including sipuncula. This is a worm
                     latin = gsub("\\s*\\([^()]+\\)", "", latin), #get rid of the (names) for some species 
                     
                     latin = ifelse(name =="ICELUS SPATULA","icelus spatula",latin),
                     latin = ifelse(name == "SPIRONTOCARIS SPINUS","spirontocaris spinus",latin),
                     latin = ifelse(name == "EUMICROTREMUS SPINOSUS","eumicrotremus spinosus",latin),
                     latin = gsub("OPHIUROIDEA","ophiuroidea",latin), #this will specify that this is to the order level.
                     latin = ifelse(name == "PORANIA (PORANIA) PULVILLUS","Porania (Porania) pulvillus",latin), #this one needs the bracket to be keyed out in worrms
                     
                     latin = trimws(latin))%>% #get rid of trailing whitespace

              filter(!grepl("eggs",tolower(latin)),# we have skate, whelk and snail/slug eggs - these are removed.
                      latin != "ORGANIC DEBRIS", #remove
                      latin != "UNIDENTIFIED", #remove
                      latin != "ASCOPHYLLUM NODOSUM", #remove - seaweed would just be captured on the upcast since the depths are much deeper than any marine algae/plants would be expected from
                      latin !="STONES AND ROCKS", #remove
                      latin !="PURSE LITTLE SKATE")%>%#this is an egg case (purse)
                mutate(latin = tolower(latin))

#taxonomy levels for the final table
PhyloNames <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")

trad_species <- id_names$latin

tax_list_trad <- data.frame()

for(i in 1:length(trad_species)){
  
  sp <- trad_species[i]
  
  message(paste0("Working on ",sp," ",i," of ",length(trad_species)))
  
  temp <- classification(sp,db="worms") #run classification
  
  temp2 <- temp[[1]]%>% #unpack classification
    data.frame()
  
  if(nrow(temp2[1])>1){
    
    temp2 <- temp2%>%
      select(rank,name)%>%
      spread(rank,name)%>%
      mutate(aphiaID = temp[[1]]%>%data.frame()%>%slice(n())%>%pull(id))
    
    temp3 <- temp2%>% #trim classification
      select(all_of(c(names(temp2)[names(temp2)%in%PhyloNames],"aphiaID")))
    
    #if data is missing or a certain taxonomic level isn't identified.
    missing_cols <- setdiff(c(PhyloNames,"aphiaID"),names(temp3))
    
    if(length(missing_cols)>0){
      
      temp3[,missing_cols] <- NA
      
      temp3 <- temp3[,c(PhyloNames,"aphiaID")]
      
    }
  }else{temp3 <- data.frame(matrix(NA, nrow = 1, ncol = length(PhyloNames)))
  names(temp3) = PhyloNames
  temp3$aphiaID = NA}
  
  temp3$species_filter <- sp # this is for linking back in the original code
  
  tax_list_trad <- rbind(tax_list_trad,temp3)
  
}                     

trad_formatted <- tax_list_trad%>%
  dplyr::select(all_of(c(PhyloNames,"aphiaID","species_filter")))%>%
  rename(latin = species_filter)%>%
  left_join(.,id_names)%>%
  rename(SCI_NAME = name)

rv_formatted <- rv_df%>%
                mutate(id = 1:n())%>%
                left_join(.,trad_formatted)%>%
                filter(!is.na(aphiaID))

save(rv_formatted,file="data/rv_taxonomy.RData")

## eDNA Data ----

