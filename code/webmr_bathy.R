#make a map of eDNA stations
library(sf)
library(tidyverse)
library(rnaturalearth)
library(rnaturalearthhires)
library(MarConsNetData)
library(ggspatial)
library(httr)
library(terra)
library(viridis)
library(tidyterra)

#Load projections
latlong <- "+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"
utm <- "+proj=utm +zone=19N +datum=WGS84 +units=m +no_defs"
utmkm <- "+proj=utm +zone=20 +datum=NAD83 +units=km +no_defs +ellps=GRS80 +towgs84=0,0,0"
CanProj <- "+proj=lcc +lat_1=49 +lat_2=77 +lat_0=63.390675 +lon_0=-91.86666666666666 +x_0=6200000 +y_0=3000000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"

palette_cust <- c("#154360", "#FF5733", "#1ABC9C")

#Basemaps
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
                      rename(status=Classification_E,name=SiteName_E)%>%
                      st_make_valid()%>%
                      st_intersection(bioregion)

plotlims_small <- maritimes_network%>%
                  filter(name=="Western/Emerald Banks Marine Refuge")%>%
                  st_transform(utmkm)%>%
                  st_buffer(10)%>%
                  st_transform(CanProj)%>%
                  st_bbox()

#load bathymetry layer
mar_bathy <- rast("data/bathymetry/gebco_cropped.tif")

contour_200m <- as.contour(mar_bathy, levels = -200)%>%
                st_as_sf()%>%
                st_transform(CanProj)

contour_300m <- as.contour(mar_bathy, levels = -300)%>%
  st_as_sf()%>%
  st_transform(CanProj)

webmr_bathy <- mar_bathy%>%
               crop(.,maritimes_network%>%
                      filter(name=="Western/Emerald Banks Marine Refuge")%>%
                      st_transform(4326)%>%
                      st_bbox()%>%
                      st_as_sfc()%>%
                      st_as_sf()%>%
                      ext())%>%
                mask(.,maritimes_network%>%
                       filter(name=="Western/Emerald Banks Marine Refuge")%>%
                       st_transform(4326))%>%
                project(CanProj)

p1 <- ggplot()+
  geom_sf(data=contour_300m,lwd=0.7,col="grey")+
  geom_sf(data=maritimes_network,fill=NA)+
  geom_spatraster(data=webmr_bathy)+
  geom_sf(data=maritimes_network%>%filter(name=="Western/Emerald Banks Marine Refuge"),col="black",lwd=1.2,fill=NA)+
  theme_bw()+
  scale_fill_viridis(na.value = "transparent")+
  coord_sf(expand=0,xlim=plotlims_small[c(1,3)],ylim=plotlims_small[c(2,4)])+
  theme(legend.position="inside",
        legend.position.inside = c(0.1,0.85),
        legend.background = element_blank())+
  labs(fill="Depth (m)")

ggsave("output/webmr_bathy.png",p1,height=6,width=6,units="in",dpi=600)





# #download bathymetry data
# gebco_extent <- bioregion%>%
#                 st_transform(4326)%>%
#                 st_bbox()%>%
#                 st_as_sfc()%>%
#                 st_as_sf()%>%
#                 ext()
# 
# gebco_file <- "data/bathymetry/gebco_2023_n66.9608_s23.607_w-83.4739_e-46.0687.nc" #this is on the R Drive data/bathymetry/gebco
# gebco_raster <- rast(gebco_file)
# 
# gebco_raster <- gebco_raster[[1]] 
# 
# gebco_cropped <- crop(gebco_raster, gebco_extent)
# 
# writeRaster(gebco_cropped, "data/bathymetry/gebco_cropped.tif", overwrite = TRUE)
