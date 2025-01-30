## make a raster that is more managable in size based on the national one for the bioregion 

#load libraries
library(tidyverse)
library(raster)

#full raster -> https://borealisdata.ca/dataverse/bluecarboncanada

oc_rast <- raster("data/OCDen.tif") #ntoe this is part of the .gitignore becausae it is so large

rast_proj <- crs(oc_rast) #reaster projection 

#crop region
bioregion <- data_planning_areas()%>%
             st_transform(utmkm)%>%
             st_buffer(20)%>% #50km buffer
             st_transform(rast_proj)%>%
             st_bbox()%>%
             st_as_sfc()%>%
             as_Spatial()

oc_rast_crop <- crop(oc_rast,bioregion)

writeRaster(oc_rast_crop,"data/OCDEN_Maritimes.tif") #for processing later
