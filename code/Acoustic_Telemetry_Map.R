### Map of acoustic reciever locations

#load libraries ----
  library(tidyverse)
  library(sf)
  library(rnaturalearth)
  library(ggspatial)
  library(arcpullr)
  library(MarConsNetData)
  library(marmap)

  s2_as_sf = FALSE

#projections ------
  latlong <- "+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"
  utm <- "+proj=utm +zone=20 +datum=NAD83 +units=km +no_defs +ellps=GRS80 +towgs84=0,0,0"
  CanProj <- "+proj=lcc +lat_1=49 +lat_2=77 +lat_0=63.390675 +lon_0=-91.86666666666666 +x_0=6200000 +y_0=3000000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"

#load network shapefiles

  #Load the Scotian Shelf-Bay of Fundy Planning Region
  bioregion <- data_planning_areas()%>%
    st_transform(CanProj)
  
  #load the network polygons 
  network <- data_draft_areas()%>%
    st_transform(CanProj)
  
#set up plotting limits
  plot_lims <- bioregion%>%
    st_bbox()%>% #get the bounding box
    st_as_sfc()%>%
    st_transform(utm)%>% #convert to a planar (km) projection
    st_buffer(25)%>% #create a buffer on that bounding box of 25 km - this is faster than doing a buffer on the polygon
    st_transform(CanProj)%>%
    st_bbox()

#get a basemap for the land
basemap <- rbind(ne_states(country = "Canada",returnclass = "sf")%>%
                   dplyr::select(name_en,geometry)%>%
                   st_as_sf()%>%
                   st_union()%>%
                   st_transform(CanProj)%>%
                   st_as_sf()%>%
                   mutate(country="Canada"),
                 ne_states(country = "United States of America",returnclass = "sf")%>%
                   dplyr::select(name_en,geometry)%>%
                   st_as_sf()%>%
                   st_union()%>%
                   st_transform(CanProj)%>%
                   st_as_sf()%>%
                   mutate(country="USA"),
                 ne_states(country = "Greenland",returnclass = "sf")%>%
                   dplyr::select(name_en,geometry)%>%
                   st_as_sf()%>%
                   st_union()%>%
                   st_transform(CanProj)%>%
                   st_as_sf()%>%
                   mutate(country="Denmark"))


#buffered polygon considered in the analysis
buffer_poly <- read_sf("data/WEBCA_10k_85k.shp")%>%st_transform(CanProj)

#Query OTN for reciever locations
proj_long_upp <- -40.00 #encompassing the regional network
proj_long_low <- -70.00
proj_lat_upp <- 60.00
proj_lat_low <- 40.00

geoserver_receivers <- readr::read_csv('https://members.oceantrack.org/geoserver/otn/ows?service=WFS&version=1.0.0&request=GetFeature&typeName=otn:stations_receivers&outputFormat=csv', guess_max = 13579)

otn_stations <- geoserver_receivers %>%
  filter(!is.na(deploy_date)) %>% # Remove rows with missing deploy_date values
  filter(stn_lat >= proj_lat_low & stn_lat <= proj_lat_upp &# Filter stations based on latitude and longitude bounds
           stn_long >= proj_long_low & stn_long <= proj_long_upp)%>%
  st_as_sf(coords=c("stn_long","stn_lat"),crs=latlong)%>%
  st_transform(CanProj)%>%
  mutate(marine = is.na(as.numeric(st_intersects(., basemap))),
         halifax_line = ifelse(grepl("HFX",station_name),"Halifax line","OTN"))%>%
  filter(marine)

#make the map

p1 <- ggplot()+
  geom_sf(data=bioregion,fill=NA)+
  geom_sf(data=basemap)+
  geom_sf(data=basemap%>%filter(country=="Canada"),fill="grey60")+
  geom_sf(data=network)+
  #geom_sf(data=buffer_poly,fill=NA,lty=2)+
  geom_sf(data=network%>%filter(SiteName_E == "Western/Emerald Banks Marine Refuge"),fill="coral2")+
  geom_sf(data=network%>%filter(SiteName_E %in% c("Musquash Estuary Marine Protected Area",
                                                  "St. Anns Bank Marine Protected Area ",
                                                  "Gully Marine Protected Area",
                                                  "Eastern Shore Islands")),fill="aquamarine")+
  geom_sf(data=otn_stations%>%filter(halifax_line == "OTN"),size=0.15,shape=20)+
  geom_sf(data=otn_stations%>%filter(halifax_line == "Halifax line"),shape=21,fill="deepskyblue2",size=0.75,col="black",lwd=0.25)+
  coord_sf(expand=0,xlim=plot_lims[c(1,3)],ylim=plot_lims[c(2,4)])+
  theme_bw()+
  annotation_scale(location="br")

ggsave("output/acoustic_map.png",p1,height=6,width=6,units="in",dpi=300)


