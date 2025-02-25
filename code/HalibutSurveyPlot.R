#Load libraries ---
library(tidyverse)
library(MarConsNetData)
library(Mar.datawrangling)
library(ROracle)
library(sf)
library(patchwork)
library(viridis)

#load halibut survey pull
load("data/halibut_survey_isdb_pull.RData") #see code at end for generation.

#Load projections
latlong <- "+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"
utmkm <- "+proj=utm +zone=20 +datum=NAD83 +units=km +no_defs +ellps=GRS80 +towgs84=0,0,0"
CanProj <- "+proj=lcc +lat_1=49 +lat_2=77 +lat_0=63.390675 +lon_0=-91.86666666666666 +x_0=6200000 +y_0=3000000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"

#Basemaps
bioregion <- data_planning_areas()%>%
  st_transform(CanProj)%>%
  st_make_valid()

#bathymetry
bathy <- read_sf("data/Shapefiles/contour_250.shp")%>%st_transform(CanProj)

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

plot_df <- bycat_dat%>%
           mutate(type=ifelse(SETCD_ID == 4, "Fixed","Random"),
                  id=paste0(id=paste(STATION,YEAR,sep="-")))%>%
           distinct(id,.keep_all=T)%>% # just he stations. 
           st_as_sf(coords=c("LONG1","LAT1"),crs=latlong)%>%
           st_transform(CanProj)%>%
           st_intersection(bioregion)

p1 <- ggplot()+
  geom_sf(data=bioregion,fill=NA)+
  geom_sf(data=basemap_atlantic)+
  geom_sf(data=bathy,lwd=0.25,col="grey80",lty=2)+
  geom_sf(data=maritimes_network%>%filter(name!="Western/Emerald Banks Marine Refuge"),fill="grey90",alpha=0.5)+
  geom_sf(data=maritimes_network%>%filter(name=="Western/Emerald Banks Marine Refuge"),linewidth=1.2,fill="grey90",alpha=0.5,col="black")+
  #geom_sf(data=plot_df%>%filter(YEAR %in% c(2016,2024)),aes(shape=type,fill=type),col="black",size=1)+ #2016 last year full fixed
  geom_sf(data=plot_df%>%filter(YEAR==2023),aes(shape=factor(YEAR)),size=0.5)+
  geom_sf(data=plot_df%>%filter(YEAR==2024),aes(shape=factor(YEAR)),fill="cornflowerblue",size=1.5)+
  scale_shape_manual(values=c(3,21))+
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
  geom_sf(data=maritimes_network%>%filter(name!="Western/Emerald Banks Marine Refuge"),fill="grey90",alpha=0.5)+
  geom_sf(data=maritimes_network%>%filter(name=="Western/Emerald Banks Marine Refuge"),fill="grey90",alpha=0.5,linewidth=1.2,col="black")+
  geom_sf(data=bathy,lwd=0.25,col="grey80",lty=2)+
  geom_sf(data=plot_df%>%filter(YEAR==2023),aes(shape=factor(YEAR)),size=0.75)+
  geom_sf(data=plot_df%>%filter(YEAR==2024),aes(shape=factor(YEAR)),fill="cornflowerblue",size=2)+
  scale_shape_manual(values=c(3,21))+
  theme_bw()+
  theme(axis.title = element_blank(),
        legend.position = "none",
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "lines"))+
  coord_sf(xlim=webca_box[c(1,3)],ylim=webca_box[c(2,4)])

combo_plot <- p1 + inset_element(p2, left = 0.55, bottom = 0.05, right = 1, top = 0.5)

ggsave("output/halibut_isdb_combo.png",combo_plot,height=5,width=5,units="in",dpi=600)



#code to do the isdb pull ----
# use.pkg = 'roracle'
# 
# get_data('isdb', 
#          data.dir ="data/isdb",
#          fn.oracle.username = "stanleyr",
#          fn.oracle.password = "homerun7",
#          fn.oracle.dsn = "PTRAN") 

# get_data('isdb', data.dir = "data/isdb")
# data_tweaks2('isdb', data.dir = "data/isdb")
# 
# ISGEARS = ISGEARS[ISGEARS$GEARCD_ID %in% c(50, 51), ] # longline and longlines set near bottom
# ISTRIPTYPECODES <- ISTRIPTYPECODES[ISTRIPTYPECODES$TRIPCD_ID %in% c(7057),]
# ISSETTYPECODES <- ISSETTYPECODES[ISSETTYPECODES$SETCD_ID %in% c(4,5),] # 5 = stratified Random, 4 fixed (2022-present)
# self_filter()
# 
# 
# ISSETPROFILE_WIDE2 <- ISSETPROFILE_WIDE %>% dplyr::select(FISHSET_ID, SET_NO, LAT1, LONG1, YEAR, DEP1, DATE_TIME1)
# ISFISHSETS2 <- ISFISHSETS %>% dplyr::select(FISHSET_ID, TRIP_ID, HAULCCD_ID, SETCD_ID, GEAR_ID, STATION, STRATUM_ID, NAFAREA_ID, LEN_LONGLINE, NUM_HOOK_HAUL)
# ISCATCHES2 <- ISCATCHES %>% dplyr::select(FISHSET_ID, CATCH_ID, SPECCD_ID, EST_NUM_CAUGHT, EST_COMBINED_WT, EST_KEPT_WT, EST_DISCARD_WT, EST_REDUCTION_WT, SPECSCD_ID)
# ISTRIPS2 <- ISTRIPS %>% dplyr::select(TRIP_ID, TRIPCD_ID, OBSCD_ID, COMMENTS)
# 
# temp_merge <- merge(ISSETPROFILE_WIDE2, ISFISHSETS2, by = "FISHSET_ID")
# temp_merge2 <- merge(x = temp_merge, y = ISCATCHES2, by = "FISHSET_ID", all = F)
# temp_merge2 <- merge(temp_merge2, ISGEARS[, c("GEAR_ID", "GEARCD_ID", "HOOKCD_ID", "HOOKSIZE")], by = "GEAR_ID")
# temp_merge2 <- merge(temp_merge2, ISTRIPS[, c("TRIP_ID", "TRIP", "TRIPCD_ID", "OBSCD_ID")])
# temp_merge2 <- merge(temp_merge2, ISSPECIESCODES[, c("SPECCD_ID", "COMMON", "SCIENTIFIC")], by = "SPECCD_ID", all.x = T)
# temp_merge2$QUARTER_YEAR <- quarter(x = temp_merge2$DATE_TIME1)
# 
# bycat_dat <- temp_merge2 %>% dplyr::select(FISHSET_ID, YEAR, DATE_TIME1, QUARTER_YEAR, TRIP_ID, TRIP, 
#                                            OBSCD_ID, TRIPCD_ID, SET_NO, LAT1, LONG1, DEP1, STATION, STRATUM_ID, 
#                                            NAFAREA_ID, SETCD_ID, HAULCCD_ID, GEARCD_ID, LEN_LONGLINE, NUM_HOOK_HAUL, 
#                                            HOOKSIZE, HOOKCD_ID, CATCH_ID, SPECCD_ID, COMMON, EST_COMBINED_WT, EST_NUM_CAUGHT, 
#                                            EST_KEPT_WT, EST_DISCARD_WT, EST_REDUCTION_WT, SCIENTIFIC)
# 
# save(bycat_dat,file="data/halibut_survey_isdb_pull.RData")
