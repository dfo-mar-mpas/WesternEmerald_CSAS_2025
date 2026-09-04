##O'brien classifiation plot

#load libraries ----
library(sf)
library(tidyverse)
library(rnaturalearth)
library(rnaturalearthhires)
library(MarConsNetData)
library(ggspatial)
library(terra)
library(tidyterra)
library(rasterVis)
library(viridis)
library(scales)
library(patchwork)

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
  st_intersection(bioregion) #clean up the edges

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

#load bioclassification 

lab_df <- data.frame(cl=1:6,
                     name=c("Slope","Laurentian Channel/Shelf Break","ESS",
                            "ESS:Banks","WSS/Outer BoF","WSS:Banks/Inner BoF"),
                     name_fr=c(
                       "Pente",
                       "Chenal Laurentien/Rupture de plateau",
                       "PNEO",
                       "PNEO : Bancs",
                       "PNEOcc/Exterieur BdF",
                       "PNEOcc : Bancs/Interieur BdF"
                     )
                    )

bioclass <- read_sf("data/Shapefiles/bioclassification_clusters.shp")%>%
            filter(region=="Maritimes")%>%
            st_transform(CanProj)%>%
            left_join(lab_df)%>%
            mutate(name=factor(name,levels=name),
                   name_fr=factor(name_fr,levels=name_fr))

webca_box <- maritimes_network%>%
             filter(name=="Western/Emerald Banks Marine Refuge")%>%
            st_transform(utmkm)%>%
            st_buffer(10)%>%
            st_transform(CanProj)%>%
            st_bbox()

webca_box_poly <- webca_box%>%st_as_sfc()

target_order <- lab_df$name[c(2, 3, 4, 1, 5, 6)] #for the legend
target_order_fr <- lab_df$name_fr[c(2, 3, 4, 1, 5, 6)] #for the legend

p1 <- ggplot()+
      geom_sf(data=bioregion,fill=NA)+
      geom_sf(data=bioclass,aes(fill=name))+
      geom_sf(data=basemap_atlantic)+
      geom_sf(data=maritimes_network,fill=NA)+
      geom_sf(data=maritimes_network%>%filter(name=="Western/Emerald Banks Marine Refuge"),linewidth=1.2,fill=NA,col="black")+
      geom_sf(data=webca_box_poly,lty=2,lwd=0.5,fill=NA)+    
      theme_bw()+
      theme(axis.title = element_blank(),
            legend.title = element_blank())+
  scale_fill_discrete(breaks = target_order)+
      coord_sf(xlim=plotlims[c(1,3)],ylim=plotlims[c(2,4)])+
      annotation_scale()

p1_fr <- ggplot()+
  geom_sf(data=bioregion,fill=NA)+
  geom_sf(data=bioclass,aes(fill=name_fr))+
  geom_sf(data=basemap_atlantic)+
  geom_sf(data=maritimes_network,fill=NA)+
  geom_sf(data=maritimes_network%>%filter(name=="Western/Emerald Banks Marine Refuge"),linewidth=1.2,fill=NA,col="black")+
  geom_sf(data=webca_box_poly,lty=2,lwd=0.5,fill=NA)+    
  theme_bw()+
  theme(axis.title = element_blank(),
        legend.title = element_blank())+
  scale_fill_discrete(breaks = target_order_fr)+
  coord_sf(xlim=plotlims[c(1,3)],ylim=plotlims[c(2,4)])+
  annotation_scale()

p2 <- ggplot()+
      geom_sf(data=bioclass,aes(fill=factor(cl)))+
      geom_sf(data=maritimes_network,fill=NA)+
      geom_sf(data=maritimes_network%>%filter(name=="Western/Emerald Banks Marine Refuge"),linewidth=1.2,fill=NA,col="black")+
      theme_bw()+
      theme(axis.title = element_blank(),
            legend.position = "none",
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            plot.margin = unit(c(0, 0, 0, 0), "lines"))+
      coord_sf(xlim=webca_box[c(1,3)],ylim=webca_box[c(2,4)])

#combination plots
combo_plot <- p1 + inset_element(p2, left = 0.55, bottom = 0.05, right = 1, top = 0.5)
combo_plot_fr <- p1_fr + inset_element(p2, left = 0.55, bottom = 0.05, right = 1, top = 0.5)

#save plots
ggsave("output/Fig15b.png",combo_plot,height=5,width=7,units="in",dpi=600)
ggsave("output/Fig15b_fr.png",combo_plot_fr,height=5,width=7.5,units="in",dpi=600)
