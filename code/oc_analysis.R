## Blue Carbon View

#make a map of eDNA stations
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

#Load projections
latlong <- "+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"
utmkm <- "+proj=utm +zone=20 +datum=NAD83 +units=km +no_defs +ellps=GRS80 +towgs84=0,0,0"
CanProj <- "+proj=lcc +lat_1=49 +lat_2=77 +lat_0=63.390675 +lon_0=-91.86666666666666 +x_0=6200000 +y_0=3000000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"

source("code/random_raster_gen_functions.R")

palette_cust <- c("#154360", "#FF5733", "#1ABC9C")

#Basemaps
bioregion <- read_sf("r:/Science/CESD/HES_MPAGroup/Data/Shapefiles/MaritimesPlanningArea.shp")%>%
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
  rename(status=Classification_E,name=SiteName_E)%>%
  st_make_valid()%>%
  st_intersection(bioregion) #clean up the edges

#load the priority areas from Epstien et al. 2025
bc_areas <- read_sf("data/Shapefiles/PA_Shapefile.shp")%>%
            st_transform(CanProj)

oc_rast <- rast("data/SSOCDen_Rock.tif")[[1]]%>% #the first is the 'mean' OC for the Scotian Shelf
           project(CanProj)

if (!st_crs(maritimes_network) == crs(oc_rast)) {
  maritimes_network <- st_transform(maritimes_network, crs(oc_rast))
}


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

#make map

p1 <- ggplot()+
  geom_sf(data=bioregion,fill=NA)+
  geom_sf(data=basemap_atlantic)+
  geom_spatraster(data=oc_rast)+
  geom_sf(data=maritimes_network,fill=NA,linewidth=0.8)+
  geom_sf(data=maritimes_network%>%filter(name=="Western/Emerald Banks Marine Refuge"),linewidth=1.2,fill=NA)+
  #geom_sf(data=bc_areas,fill="red")+
  theme_bw()+
  coord_sf(expand=0,xlim=plotlims[c(1,3)],ylim=plotlims[c(2,4)])+
  scale_fill_viridis_b(na.value = "transparent")+
  labs(fill="OC Density")

ggsave("output/WEBMR_OC.png",p1,width=6,height=6,units="in",dpi=300)

#coverage analysis
masked_oc <- mask(oc_rast,vect(maritimes_network))

total_oc_conservation <- sum(values(masked_oc),na.rm=T)
total_oc_region <- sum(values(oc_rast),na.rm=T)

total_oc_conservation / total_oc_region


maritimes_oc <- NULL

for (i in 1:nrow(maritimes_network)) {
  
  message("Working on ",maritimes_network[i,]%>%pull(name)," site ",i," of ",nrow(maritimes_network))
  
  # Extract individual site as a subset
  site <- maritimes_network[i, ]
  
  # Mask raster with the individual site
  site_masked_raster <- mask(oc_rast, vect(site))
  
  # Calculate total OC density for the site
  total_oc_site <- sum(values(site_masked_raster),na.rm = TRUE)
  
  # Calculate proportion of OC density for this site
  proportion_site <- total_oc_site / total_oc_region
  
  # Store results
  
  out <- data.frame(site = maritimes_network[i,]%>%pull(name),
                    total_oc = total_oc_site,
                    prop_oc = proportion_site)
  
  maritimes_oc <- rbind(maritimes_oc,out)
  
  
}

plot_oc <- maritimes_oc%>%
           arrange(prop_oc)%>%
           mutate(name_oc = factor(site,levels=site))

p2 <- ggplot(plot_oc%>%filter(prop_oc>0.005),aes(x=prop_oc,y=name_oc,fill=prop_oc))+
    geom_bar(stat="identity",col="black")+
    theme_bw()+
    scale_x_continuous(labels = percent, expand = c(0, 0.001))+ 
    scale_fill_viridis_b(labels = percent)+
    labs(y="",x="Proportion coverage of total regional OC",fill="")+
    theme(legend.position = "inside",
          legend.position.inside = c(0.9,0.15),
          legend.background = element_blank())

ggsave("output/OC_analysis.png",p2,width=6,height=6,units="in",dpi=300)

  
