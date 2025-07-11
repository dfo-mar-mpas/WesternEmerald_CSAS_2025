## 2024 Ecosystem Survey map

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
library(basemaps)
library(viridis)
library(scales)
library(patchwork)
library(magick)

#Load projections
latlong <- "+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"
utmkm <- "+proj=utm +zone=20 +datum=NAD83 +units=km +no_defs +ellps=GRS80 +towgs84=0,0,0"
CanProj <- "+proj=lcc +lat_1=49 +lat_2=77 +lat_0=63.390675 +lon_0=-91.86666666666666 +x_0=6200000 +y_0=3000000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"

#Basemaps
bioregion <- data_planning_areas()%>%
  st_transform(CanProj)%>%
  st_make_valid()

#load in the 2024 survey
rv_df <- read.csv("c:/Users/stanleyr/Documents/Github/mar_network/data/eDNA sampling/2024_Summer_Stations.csv")%>%
         st_as_sf(coords=c("LON_DD","LAT_DD"),crs=latlong)%>%
         st_transform(CanProj)%>%
         mutate(set_type = ifelse(TYPE=="PRIMARY","Primary set","Secondary set"))

rv_strata <- read_sf("c:/Users/stanleyr/Documents/Github/mar_network/data/shapefiles/MaritimesRegionEcosystemAssessmentStrata(2014-).shp")%>%
             st_transform(CanProj)%>%
            st_intersection(bioregion)



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

webca_box <- maritimes_network%>%
  filter(name=="Western/Emerald Banks Marine Refuge")%>%
  st_transform(utmkm)%>%
  st_buffer(10)%>%
  st_transform(CanProj)%>%
  st_bbox()

webca_box_poly <- webca_box%>%st_as_sfc()


p1 <- ggplot()+
  geom_sf(data=rv_strata,fill=NA,col="grey80",lwd=0.3)+
  geom_sf(data=bioregion,fill=NA)+
  geom_sf(data=basemap_atlantic)+
  geom_sf(data=maritimes_network%>%filter(name!="Western/Emerald Banks Marine Refuge"),fill="grey90",alpha=0.5)+
  geom_sf(data=maritimes_network%>%filter(name=="Western/Emerald Banks Marine Refuge"),linewidth=1.2,fill="grey90",alpha=0.5,col="black")+
  geom_sf(data=rv_df%>%filter(TYPE=="PRIMARY"),aes(shape=set_type),fill="cornflowerblue",col="black",size=1)+
  geom_sf(data=rv_df%>%filter(TYPE!="PRIMARY"),aes(shape=set_type),col="black",size=0.6)+
  scale_shape_manual(values=c(21,3))+
  geom_sf(data=webca_box_poly,lty=2,lwd=0.5,fill=NA)+    
  theme_bw()+
  theme(axis.title = element_blank(),
        legend.position="inside",
        legend.position.inside = c(0.84,0.92),
        legend.title = element_blank(),
        legend.background = element_blank())+
  coord_sf(xlim=plotlims[c(1,3)],ylim=plotlims[c(2,4)])+
  annotation_scale()

p2 <- ggplot()+
  geom_sf(data=rv_strata,fill=NA,col="grey80",lwd=0.3)+
  geom_sf(data=maritimes_network%>%filter(name!="Western/Emerald Banks Marine Refuge"),fill="grey90",alpha=0.5)+
  geom_sf(data=maritimes_network%>%filter(name=="Western/Emerald Banks Marine Refuge"),fill="grey90",alpha=0.5,linewidth=1.2,col="black")+
  geom_sf(data=rv_df%>%filter(TYPE=="PRIMARY"),aes(shape=set_type),fill="cornflowerblue",col="black",size=2.5)+
  geom_sf(data=rv_df%>%filter(TYPE!="PRIMARY"),aes(shape=set_type),col="black",size=1)+
  scale_shape_manual(values=c(21,3))+
  theme_bw()+
  theme(axis.title = element_blank(),
        legend.position = "none",
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "lines"))+
  coord_sf(xlim=webca_box[c(1,3)],ylim=webca_box[c(2,4)])

combo_plot <- p1 + inset_element(p2, left = 0.55, bottom = 0.05, right = 1, top = 0.5)

ggsave("output/webca_rv_combo.png",combo_plot,height=5,width=5,units="in",dpi=600)

##make a plot with the RV strata with full extent

webmr <- maritimes_network%>%filter(name=="Western/Emerald Banks Marine Refuge")

webmr_strata <- rv_strata%>%
                st_intersection(webmr)%>%
                mutate(area=round(st_area(.)/1000/1000,2))

rv_area_df <- rv_strata%>%
              filter(StrataID %in% unique(webmr_strata$StrataID))%>%
              mutate(area=round(st_area(.)/1000/1000,2))

webmr_strata <- webmr_strata%>%
                mutate(prop_area = round(webmr_strata$area/rv_area_df$area,3)*100,
                       lab = paste0(StrataID," (",prop_area,"%)"))%>%
                arrange(-prop_area)%>%
                mutate(lab = factor(lab,levels=lab))%>%
                filter(as.numeric(prop_area)>5) #cut out the small areas

webmr_strata_full <- rv_strata%>%
                     filter(StrataID %in% unique(webmr_strata$StrataID))%>%
                     left_join(.,webmr_strata%>%data.frame()%>%dplyr::select(StrataID,lab))

plot_lims2 <- rv_strata%>%
              filter(StrataID %in% unique(webmr_strata$StrataID))%>%
              st_transform(utmkm)%>%
              st_buffer(5)%>%
              st_transform(CanProj)%>%
              st_bbox()
        


p3 <- ggplot()+
  geom_sf(data=basemap_atlantic)+
  geom_sf(data=maritimes_network%>%filter(name!="Western/Emerald Banks Marine Refuge"),fill="grey90",alpha=0.5)+
  geom_sf(data=rv_strata,fill=NA,col="grey80",lwd=0.3)+
  geom_sf(data=webmr_strata_full,aes(fill=lab),alpha=0.7)+
  geom_sf(data=webmr_strata,aes(fill=lab))+
  geom_sf(data=webmr,fill=NA,lwd=1.2)+
  geom_sf(data=rv_df%>%filter(TYPE=="PRIMARY"),aes(shape=set_type),fill="cornflowerblue",col="black",size=2.5,show.legend = FALSE)+
  geom_sf(data=rv_df%>%filter(TYPE!="PRIMARY"),aes(shape=set_type),col="black",size=1,show.legend = FALSE)+
    scale_shape_manual(values=c(21,3))+
    theme_bw()+
    theme(axis.title = element_blank(),
          #legend.position = "none",
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          plot.margin = unit(c(0, 0, 0, 0), "lines"))+
  coord_sf(expand=0,xlim=plot_lims2[c(1,3)],ylim=plot_lims2[c(2,4)])+
  labs(fill="RV Strata")+
  annotation_scale(location = "br")


p4 <- ggplot()+
  geom_sf(data=rv_strata,fill=NA,col="grey80",lwd=0.3)+
  geom_sf(data=bioregion,fill=NA)+
  geom_sf(data=basemap_atlantic)+
  geom_sf(data=maritimes_network%>%filter(name!="Western/Emerald Banks Marine Refuge"),fill="grey90",alpha=0.5)+
  geom_sf(data=maritimes_network%>%filter(name=="Western/Emerald Banks Marine Refuge"),linewidth=1.2,fill="grey90",alpha=0.5,col="black")+
  geom_sf(data=rv_df%>%filter(TYPE=="PRIMARY"),aes(shape=set_type),fill="cornflowerblue",col="black",size=0.85)+
  geom_sf(data=rv_df%>%filter(TYPE!="PRIMARY"),aes(shape=set_type),col="black",size=0.5)+
  scale_shape_manual(values=c(21,3))+
  geom_sf(data=plot_lims2%>%st_as_sfc(),lty=2,lwd=0.5,fill=NA)+    
  theme_bw()+
  theme(axis.title = element_blank(),
        legend.position="inside",
        legend.position.inside = c(0.80,0.92),
        legend.title = element_blank(),
        legend.background = element_blank())+
  coord_sf(xlim=plotlims[c(1,3)],ylim=plotlims[c(2,4)])+
  annotation_scale()+
  guides(shape = guide_legend(override.aes = list(size = 2.2)))
  
#assemble and save
combo_plot2 <- p4 + p3 + plot_layout(ncol=2)


combo_path <- "output/rv_strata_webmr_plot.png"

ggsave(combo_path,combo_plot2,height=5,width=9,units="in",dpi=300)

image <- image_read(combo_path)

image_write(image_trim(image), path = combo_path)
