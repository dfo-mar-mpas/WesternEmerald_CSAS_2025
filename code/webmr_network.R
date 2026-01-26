##MCT project call map

## select stations for the 2024 summer RV survey

#load libraries ----
library(sf)
library(tidyverse)
library(rnaturalearth)
library(viridis)
library(ggnewscale)
library(patchwork)
library(ggspatial)
library(marmap)
library(MarConsNetData)

curdir <- getwd()

#projections ----
latlong <- "+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"
utm <- "+proj=utm +zone=20 +datum=NAD83 +units=km +no_defs +ellps=GRS80 +towgs84=0,0,0"
CanProj <- "+proj=lcc +lat_1=49 +lat_2=77 +lat_0=63.390675 +lon_0=-91.86666666666666 +x_0=6200000 +y_0=3000000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"

#load the network shapefile ----

#load bioregion
bioregion <- data_planning_areas()%>%
  st_transform(CanProj)%>%
  st_make_valid()

network <- data_draft_areas()%>%
            st_transform(CanProj)%>%
            st_make_valid()%>%
            mutate(class = ifelse(grepl("Tier",Classification_E),"Draft","Existing"))%>%
            st_intersection(.,bioregion)%>%
            mutate(type=case_when(grepl("Marine Protected",SiteName_E) ~ "MPA",
                                  grepl("Marine Refuge",SiteName_E) ~ "MR",
                                  Classification_E == "Areas of Interest (AOI)" ~ "AOI",
                                  LeadAgency_E %in% c("Canadian Wildlife Service","Parks Canada") ~ "OECM",
                                  TRUE ~ "Draft"),
                   type_fr=case_when(type=="MPA" ~ "ZPM", #french translations 
                                     type=="MR" ~ "RM",
                                     type=="AOI" ~ "SI",
                                     type=="OECM" ~ "AMCE",
                                     type=="Draft" ~ "Provisoire",
                                     TRUE ~ "NA"),
                   type=factor(type,levels=c("MPA","MR","AOI","OECM","Draft")),
                   type_fr=factor(type_fr,levels=c("ZPM","RM","SI","AMCE","Provisoire")))

##Basemap ----
plot_bounds <- bioregion%>%
  st_transform(utm)%>%
  st_buffer(25)%>%
  st_transform(CanProj)%>%
  st_bbox()


basemap <- ne_states(country = "Canada",returnclass = "sf")%>%
  dplyr::select(name_en,geometry)%>%
  mutate(country = "Canada")%>%
  st_union()%>%
  st_as_sf()%>%
  rbind(.,
        ne_states(country = "United States of America",returnclass = "sf")%>%
          dplyr::select(name_en,geometry)%>%
          mutate(country = "United States of America")%>%
          st_union()%>%
          st_as_sf())%>%
  st_transform(CanProj)%>%
  st_as_sf()

#download the bathymetric data

bathy <- read_sf("data/bathymetry/contour_250.shp")%>%
         st_transform(CanProj)%>%
         st_intersection(bioregion)


p1 <- ggplot()+
  geom_sf(data=bathy,fill=NA,lwd=0.25,col="grey80")+
  geom_sf(data=bioregion,fill=NA,col="grey30")+
  geom_sf(data=basemap,fill="white",col="black")+
  geom_sf(data=network,aes(fill=type))+
  geom_sf(data=network%>%filter(SiteName_E ==  "Western/Emerald Banks Marine Refuge"),,col="black",linewidth=1.3,fill=NA)+
  scale_fill_viridis(discrete=TRUE,option="D")+
  labs(fill="")+
  theme_bw()+
  theme(legend.position="inside",
        legend.position.inside = c(0.9,0.1),
        legend.background = element_blank(),
        legend.title = element_blank())+
  coord_sf(expand=0,xlim=plot_bounds[c(1,3)],ylim=plot_bounds[c(2,4)])+
  annotation_scale(location="bl")

#french translantion 

p1_fr <- ggplot()+
  geom_sf(data=bathy,fill=NA,lwd=0.25,col="grey80")+
  geom_sf(data=bioregion,fill=NA,col="grey30")+
  geom_sf(data=basemap,fill="white",col="black")+
  geom_sf(data=network,aes(fill=type_fr))+
  geom_sf(data=network%>%filter(SiteName_E ==  "Western/Emerald Banks Marine Refuge"),,col="black",linewidth=1.3,fill=NA)+
  scale_fill_viridis(discrete=TRUE,option="D")+
  scale_x_continuous(
    labels = function(x) paste0(abs(x), "° O")
  )+
  scale_y_continuous(
    labels = function(y) paste0(abs(y), "° N")
  )+
  labs(fill="")+
  theme_bw()+
  theme(legend.position="inside",
        legend.position.inside = c(0.9,0.1),
        legend.background = element_blank(),
        legend.title = element_blank())+
  coord_sf(expand=0,xlim=plot_bounds[c(1,3)],ylim=plot_bounds[c(2,4)])+
  annotation_scale(location="bl")

ggsave("output/webmr_network.png",p1 ,width=8,height=7,units="in",dpi=300)
ggsave("output/webmr_network_fr.png",p1_fr ,width=8,height=7,units="in",dpi=300)
