
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
library(RColorBrewer)

#Load projections
latlong <- "+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"
utmkm <- "+proj=utm +zone=20 +datum=NAD83 +units=km +no_defs +ellps=GRS80 +towgs84=0,0,0"
CanProj <- "+proj=lcc +lat_1=49 +lat_2=77 +lat_0=63.390675 +lon_0=-91.86666666666666 +x_0=6200000 +y_0=3000000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"

#Basemaps
bioregion <- data_planning_areas()%>%
  st_transform(CanProj)%>%
  st_make_valid()

bioregions <- read_sf("r:/Science/CESD/HES_MPAGroup/Data/Shapefiles/DFO_Marine_Bioregions_Clipped_1M_CAEAC_2012_05_31.shp")%>%
              st_transform(CanProj)%>%
              filter(!Legend %in% c("Saint-Pierre et Miquelon / Saint-Pierre et Miquelon","13. Great Lakes / Grands Lacs"))%>%
              mutate(name = gsub("^[0-9]+\\.\\s*|\\s*/.*", "", Legend))

eez <- read_sf("r:/Science/CESD/HES_MPAGroup/Data/Shapefiles/Canada_EEZ.shp")%>%st_transform(CanProj)

#plotting limits-----
bioregion_lims <- bioregions%>%
                  st_bbox()%>%
                  st_as_sfc()%>%
                  st_buffer(10)%>%
                  st_bbox()


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
                            mutate(country="USA"),
                         ne_states(country = "Greenland",returnclass = "sf")%>%
                            dplyr::select(name_en,geometry)%>%
                            st_as_sf()%>%
                            st_union()%>%
                            st_transform(latlong)%>%
                            st_as_sf()%>%
                            mutate(country="Greenland"),
                          ne_states(country = "France",returnclass = "sf")%>%
                            dplyr::select(name_en,geometry)%>%
                            st_as_sf()%>%
                            st_union()%>%
                            st_transform(latlong)%>%
                            st_as_sf()%>%
                            mutate(country="France"),
                         ne_states(country = "Iceland",returnclass = "sf")%>%
                           dplyr::select(name_en,geometry)%>%
                           st_as_sf()%>%
                           st_union()%>%
                           st_transform(latlong)%>%
                           st_as_sf()%>%
                           mutate(country="Iceland"))%>%
                    st_transform(CanProj)

#bathymetric layer -----
mar_bathy <- rast("data/bathymetry/gebco_cropped.tif")


contour_300m <- as.contour(mar_bathy, levels = -300)%>%
                st_as_sf()%>%
                st_transform(CanProj)%>%
                st_intersection(bioregion)

#select sites for the EBM slide ----
target_sites <- c("Gully Marine Protected Area","Eastern Shore Islands","Musquash Estuary Marine Protected Area",
                  "Fundian Channel-Browns Bank","St. Anns Bank Marine Protected Area ","Eastern Canyons Marine Refuge",
                  "Western/Emerald Banks Marine Refuge")

#set colours for the target sitse
site_colors <- setNames(brewer.pal(length(target_sites), "Set2"), target_sites)


##\MPAs and OECMs in Canada -----
mpas <- read_sf("c:/Users/stanleyr/Documents/Github/CanMonitoringFramework/data/shapefiles/MPAs.shp")%>%st_transform(CanProj)
tht <- read_sf("c:/Users/stanleyr/Documents/Github/CanMonitoringFramework/data/shapefiles/ThT_MPA.shp")%>%st_transform(CanProj)

mpas[grepl("offshore",tolower(mpas$NAME_E)),"geometry"] <- tht%>%dplyr::select(geometry)


CCAs <- mpas%>%
  mutate(name=ifelse(grepl("Eastport",NAME_E),"Eastport Marine Protected Area",NAME_E))%>%
  filter(OWNER_T == 1)%>%
  st_transform(CanProj)%>%
  st_intersection(bioregions%>%dplyr::select(Legend))#government of conservation areas


#Plot 1 (Regional) _ Scotian SHelf-Bay of Fundy Bioregion with target sites fill coloured
p1 <- ggplot()+
  geom_sf(data=bioregion,fill=NA)+
  geom_sf(data=basemap_atlantic%>%filter(country=="Canada"),fill="grey55")+
  geom_sf(data=contour_300m,lwd=0.6,col="grey")+
  geom_sf(data=maritimes_network,fill=NA)+
  geom_sf(data=maritimes_network%>%filter(name %in% target_sites),aes(fill=name))+
  geom_sf(data=maritimes_network%>%filter(name=="Western/Emerald Banks Marine Refuge"),linewidth=1.2,alpha=0.5,col="black")+ 
  theme_bw()+
  scale_fill_manual(values = site_colors)+
  theme(axis.title = element_blank(),
        legend.position="none")+
  coord_sf(xlim=plotlims[c(1,3)],ylim=plotlims[c(2,4)])+
  annotation_scale()

ggsave("output/network_polys/network_polys.png",p1,height=6,width=6,units="in",dpi=600)


#loop to make images for each site with corresponding colours to P1 for the figures. 
    # for (site in target_sites) {
    #   
    #   site_data <- maritimes_network %>% filter(name == site)
    #   
    #   p <- ggplot() +
    #     geom_sf(data = site_data, fill = site_colors[site], color = "black",lwd=1.1) +  # Fill with assigned color
    #     theme_void() 
    #   
    #   if(site == "Western/Emerald Banks Marine Refuge"){site = gsub("/"," ",site)}
    #   
    #   # Save the plot
    #   ggsave(filename = paste0("output/network_polys/", site, ".png"), 
    #          plot = p, 
    #          width = 4, height = 4, dpi = 300)
    # }

#Canadian bioregions map with focal extent of p1 --------
p2 <- ggplot()+
  geom_sf(data=bioregions,aes(fill=Legend))+
  geom_sf(data=basemap_atlantic)+
  geom_sf(data=basemap_atlantic%>%filter(country=="Canada"),fill="grey55")+
  geom_sf(data=plotlims%>%st_as_sfc(),fill=NA,linetype=2)+
  theme_bw()+
  theme(legend.position="none")+
  scale_fill_viridis(discrete=TRUE)+
  coord_sf(expand=0,xlim=bioregion_lims[c(1,3)],ylim=bioregion_lims[c(2,4)])

ggsave("output/network_polys/Canada_bioregions.png",p2,height=6,width=6,units="in",dpi=600)

#Canadian Network map with sites coloured to match bioregions (p2)
p3 <- ggplot()+
  geom_sf(data=eez,fill=NA)+
  geom_sf(data=basemap_atlantic)+
  geom_sf(data=basemap_atlantic%>%filter(country=="Canada"),fill="grey55")+
  geom_sf(data=CCAs,aes(fill=Legend))+
  geom_sf(data=plotlims%>%st_as_sfc(),fill=NA,linetype=2)+
  theme_bw()+
  theme(legend.position="none")+
  scale_fill_viridis(discrete=TRUE)+
  coord_sf(expand=0,xlim=bioregion_lims[c(1,3)],ylim=bioregion_lims[c(2,4)])
  
ggsave("output/network_polys/Canada_ccas.png",p3,height=6,width=6,units="in",dpi=600)
