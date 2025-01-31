## plot up the Scotian Shelf Ichthyoplankton Project (SSPI) data

#load libraries
library(tidyverse)
library(sf)
library(rnaturalearth)
library(MarConsNetData)

#Load projections
latlong <- "+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"
utmkm <- "+proj=utm +zone=20 +datum=NAD83 +units=km +no_defs +ellps=GRS80 +towgs84=0,0,0"
CanProj <- "+proj=lcc +lat_1=49 +lat_2=77 +lat_0=63.390675 +lon_0=-91.86666666666666 +x_0=6200000 +y_0=3000000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"

#Basemaps
bioregion <- data_planning_areas()%>%
  st_transform(CanProj)%>%
  st_make_valid()

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
  rename(status=Classification_E,name=SiteName_E)%>%
  st_make_valid()%>%
  st_intersection(bioregion) #

#load the SSPI shapefiles
larval_poly <- read_sf("data/SSIP/SSIPdiversity_poly.shp")%>%
               st_transform(CanProj)%>%
               st_intersection(bioregion)

larval_cons <- larval_poly%>%st_intersection(maritimes_network)

#load land basemap
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

ggplot()+
  geom_sf(data=bioregion,fill=NA)+
  geom_sf(data=basemap_atlantic)+
  geom_sf(data=larval_cons,fill="coral2")+
  geom_sf(data=maritimes_network,fill=NA,linewidth=0.5)+
  geom_sf(data=maritimes_network%>%filter(name=="Western/Emerald Banks Marine Refuge"),linewidth=1.2,fill=NA,col="black")+
  theme_bw()+
  coord_sf(expand=0,xlim=plotlims[c(1,3)],ylim=plotlims[c(2,4)])
