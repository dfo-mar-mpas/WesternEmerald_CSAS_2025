#Script to build a table of all eDNA taxonomy from the 2020 survey

#load libraries ---------
library(tidyverse)
library(rnaturalearth)
library(MarConsNetData)
library(ROracle)
library(Mar.datawrangling)
library(sf)

#Load projections --------
latlong <- "+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"
utm <- "+proj=utm +zone=19N +datum=WGS84 +units=m +no_defs"
utmkm <- "+proj=utm +zone=20 +datum=NAD83 +units=km +no_defs +ellps=GRS80 +towgs84=0,0,0"
CanProj <- "+proj=lcc +lat_1=49 +lat_2=77 +lat_0=63.390675 +lon_0=-91.86666666666666 +x_0=6200000 +y_0=3000000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"

#Load the cleaned datasets -----------
load("data/rv_taxonomy.RData")
load("data/edna_taxonomy.RData")
load("data/edna_taxonomy_other.RData") #those not in the buffered area around WEBMR

PhyloNames <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")


#Load base mapping polygons --------------
bioregion <- data_planning_areas()%>% 
              st_transform(CanProj)

#plotlimits
plotlims <- bioregion%>%
            st_transform(utmkm)%>%
            st_buffer(20)%>% #50km buffer
            st_transform(CanProj)%>%
            st_bbox()

#maritimes network -- note that this cannot be shared or made available online to the public
maritimes_network <- data_draft_areas()%>%
                    st_transform(CanProj)%>%
                    st_make_valid()%>%
                    dplyr::select(Classification_E,SiteName_E)%>%
                    rename(status=Classification_E,name=SiteName_E)


#load the webca buffer
webca_buffer <- read_sf("data/Shapefiles/WEBCA_10k_85k.shp")%>%
                st_transform(CanProj)%>%
                mutate(status="buffer",name="webca_buffer")%>%
                rename(geoms=geometry)%>%
                st_difference(maritimes_network%>%filter(name=="Western/Emerald Banks Marine Refuge"))%>%
                dplyr::select(status,name)

edna_buffer <- maritimes_network%>% #this will grab the eDNA stations within 50 km of webca
              filter(name=="Western/Emerald Banks Marine Refuge")%>%
              st_transform(utmkm)%>%
              st_buffer(50)%>%
              st_transform(CanProj)%>%
              mutate(name='buffer_50km')

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
                            st_transform(CanProj)

#Groundfish diet data --- 
diet_df <- read.csv("data/Diet_Data_Allregions.csv")

haddock_diet <- diet_df%>%filter(grepl("aeglefinus",predator_species))


## create a list of species detected by eDNA across all surveys




edna_webmr <- edna_data%>%
              distinct(latin,.keep_all=TRUE)%>%
              dplyr::select(all_of(c(PhyloNames,"latin")))

edna_outside_webmr <- edna_data_other%>%
                      distinct(latin,.keep_all=TRUE)%>%
                      dplyr::select(all_of(c(PhyloNames,"latin")))%>%
                      filter(!latin %in% edna_webmr$latin)


#load diet databsae

focal_sp <- c("HALIBUT(ATLANTIC)","AMERICAN PLAICE","COD(ATLANTIC)",
             "WINTER SKATE","HADDOCK","WHITE HAKE") #based on the Murillo et al. Res Doc indicator table of focal groundfish species

# un = "" #set this with the username
# pw = "" #set this with the password
# 
# cxn <- ROracle::dbConnect(DBI::dbDriver("Oracle"), un, pw, "PTRAN")
# 
# get_data(db="rv",cxn=cxn)
# 
# 
# diet_db <- STOMACH_DATA_VW%>%
#            rename(CODE = SPEC)%>%
#            left_join(GSSPECIES%>%select(SPEC,COMM,CODE))%>%
#            rename(preycode = PREYSPECCD)%>%
#            left_join(GSSPECIES%>%
#                        rename(preycode = CODE,prey_comm=COMM,prey=SPEC)%>%
#                        select(prey,prey_comm,preycode))%>%
#            mutate(date=as.POSIXct(SDATE),
#                   year=year(date),
#                   month=month(date))%>%
#            filter(!is.na(SLONGDD),!is.na(SLATDD))%>%
#            st_as_sf(coords = c("SLONGDD","SLATDD"),crs=latlong)%>%
#            st_transform(CanProj)%>%
#            st_join(., webca_buffer, join = st_intersects)

#if you want more details about the size of fish you can connect to GSDET 

#from Mike McMahon
#"I think you can hope to join to GSDET using FSHNO (in addition to MISSION, SETNO and SPEC), but we found in the last week or so that FSHNO is not a great/reliable linkage.  Manon just discovered it, and is not happy at all, but if you don't need to join prey to a spcific predator, you should be able to generate patterns"
                     
#save(diet_db,file="data/diet_db.RData") 

load("data/diet_db.RData")

