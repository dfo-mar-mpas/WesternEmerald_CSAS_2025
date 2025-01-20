#make a map of eDNA stations
library(sf)
library(tidyverse)
library(rnaturalearth)
library(rnaturalearthhires)
library(MarConsNetData)
library(ggspatial)

#Load projections
latlong <- "+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"
utm <- "+proj=utm +zone=19N +datum=WGS84 +units=m +no_defs"
utmkm <- "+proj=utm +zone=20 +datum=NAD83 +units=km +no_defs +ellps=GRS80 +towgs84=0,0,0"
CanProj <- "+proj=lcc +lat_1=49 +lat_2=77 +lat_0=63.390675 +lon_0=-91.86666666666666 +x_0=6200000 +y_0=3000000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"

palette_cust <- c("#154360", "#FF5733", "#1ABC9C")

#Basemaps
bioregion <- read_sf("r:/Science/CESD/HES_MPAGroup/Data/Shapefiles/MaritimesPlanningArea.shp")%>%
             st_transform(CanProj)

# bioregion <- data_bioregion()%>% #isn't quite the planning region
#   st_transform(latlong)

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

maritimes_network_trimmed <- maritimes_network%>%
                              st_make_valid()%>%
                              st_intersection(bioregion)

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

#Load the cleaned data
load("data/RVSurveyPull.RData")
load("data/edna_taxonomy.RData") #'edna_data'

edna_comp_stations <- edna_data%>% #this is the data with the stations close to or inside WEBCA
                      mutate(stations=paste(gsub("-","",cruise),setno,sep="-"))%>%
                      distinct(stations)%>%
                      pull(stations)

edna_paired_staions <- read.csv("data/merged_edna_samples.csv")%>% #this is the eDNA sample metadata from the 2020 survey where it was collected 
                        st_as_sf(coords=c("longitude","latitude"),crs=latlong)%>%
                        mutate(MISSION = gsub("-","",cruise),
                               stations = paste(MISSION,setno,sep="-"))%>%
                        distinct(stations)%>%
                        pull(stations)%>%
                        setdiff(edna_comp_stations)

edna_rv_df  <- GSINF %>%
               filter(MISSION == "NED2020025")%>%
               mutate(stations = paste(MISSION,SETNO,sep="-"))%>%
               distinct(stations,.keep_all = TRUE)%>%
               mutate(type = case_when(stations %in% edna_paired_staions ~ "Paired samples",
                                       stations %in% edna_comp_stations ~ "WEBMR comparison",
                                       TRUE ~ "RV sets"),
                      type = factor(type,levels=c("RV sets","Paired samples","WEBMR comparison")))%>%
               st_as_sf(coords = c("MLONG", "MLAT"), crs = latlong)%>%
               st_transform(CanProj)
     
#make the map
map_plot <- ggplot()+
            geom_sf(data=bioregion,fill=NA)+
            geom_sf(data=basemap_atlantic)+
            geom_sf(data=maritimes_network_trimmed,fill=NA)+
            geom_sf(data=maritimes_network%>%filter(name == "Western/Emerald Banks Marine Refuge"),fill="coral")+
            geom_sf(data=webca_buffer,lty=2,col="grey75",lwd=0.5,fill=NA)+
            geom_sf(data=edna_rv_df%>%filter(type=="RV sets"),aes(fill=type),pch=21,size=1.2)+
            geom_sf(data=edna_rv_df%>%filter(type!="RV sets"),aes(fill=type),size=3,pch=21)+
            theme_bw()+
            coord_sf(expand = 0,xlim=plotlims[c(1,3)],ylim=plotlims[c(2,4)])+
            theme(legend.position="inside",
                  legend.position.inside = c(0.88,0.08),
                  legend.title = element_blank(),
                  legend.background = element_blank(),
                  legend.key = element_blank())+
            annotation_scale(location="tr")+
            scale_fill_manual(values=palette_cust)

ggsave("output/eDNA_sample_map.png",map_plot,height=8,width=8,units="in",dpi=300)
