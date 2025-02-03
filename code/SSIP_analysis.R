## plot up the Scotian Shelf Ichthyoplankton Project (SSPI) data

#load libraries
library(tidyverse)
library(sf)
library(rnaturalearth)
library(MarConsNetData)
library(purrr)
library(scales)
library(viridis)
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
  st_intersection(bioregion) #

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

#load the SSPI shapefiles
larval_poly <- read_sf("data/SSIP/SSIPdiversity_poly.shp")%>%
               st_transform(CanProj)%>%
               st_intersection(bioregion)

larval_cons <- larval_poly%>%st_intersection(maritimes_network)

ssip_data <- data.frame(species=c("Haddock",
                                  "Mackerel",
                                  "American plaice",
                                  "Pollock",
                                  "Redfish",
                                  "Silver hake",
                                  "Witch flounder",
                                  "Yellowtail flounder",
                                  "Fish diversity"),
                        shp = c("data/SSIP/SSIP_Haddock_80p.shp",
                                "data/SSIP/SSIP_Mackerel_80p.shp",
                                "data/SSIP/SSIP_Plaice_80p.shp",
                                "data/SSIP/SSIP_Pollock_80p.shp",
                                "data/SSIP/SSIP_Redfish_80p.shp",
                                "data/SSIP/SSIP_SilverHake_80p.shp",
                                "data/SSIP/SSIP_WitchFlounder_80p.shp",
                                "data/SSIP/SSIP_YellowtailFlounder_80p.shp",
                                "data/SSIP/SSIPdiversity_poly.shp"))

ssip_sf <- map2_dfr(ssip_data$shp, ssip_data$species, function(shp_path, sp) {
  st_read(shp_path, quiet = TRUE) %>%
    mutate(species = sp)%>%
    st_transform(CanProj)%>%
    st_make_valid()%>%
    st_intersection(bioregion)%>%
    group_by(species)%>%
    summarise(geometry = st_union(geometry),
              area = as.numeric(st_area(geometry))/1000000 )
})


ssip_webca_sf <- ssip_sf%>%
  st_intersection(maritimes_network%>%filter(name=="Western/Emerald Banks Marine Refuge"))%>%
  dplyr::select(species,geometry)%>%
  group_by(species)%>%
  summarise(
    geometry = st_union(geometry),
    area = as.numeric(st_area(geometry))/1000000  # Area in km²
  )

diveristy_target <- ssip_sf%>%filter(species=="Fish diversity")%>%pull(area)*0.2

ssip_network_sf <- ssip_sf%>%
                   st_intersection(maritimes_network%>%st_union())%>%
                   dplyr::select(species,geometry)%>%
                   group_by(species)%>%
                    summarise(
                      geometry = st_union(geometry),
                      area = as.numeric(st_area(geometry))/1000000  # Area in km²
                    )%>%
                    left_join(ssip_webca_sf%>%
                              data.frame()%>%
                              dplyr::select(species,area)%>%
                              rename(webca_area=area))%>%
                    left_join(ssip_sf%>%
                             data.frame()%>%
                             dplyr::select(species,area)%>%
                             rename(bioregion_area=area))%>%
                    mutate(webca_prop = round((webca_area/bioregion_area*100),2))

ordered_species <- ssip_network_sf %>%
  arrange(desc(webca_prop)) %>%  # Sort by webca_prop
  pull(species)  # Extract species names
                  
ordered_species <- c(
  "Fish diversity",
  ordered_species[ordered_species != "Fish diversity"]
)

ssip_network_sf <- ssip_network_sf %>%
  mutate(species = factor(species, levels = ordered_species))%>%
  data.frame()

ssip_sf <- ssip_sf %>%
  mutate(species = factor(species, levels = ordered_species))

ssip_sf_other <- ssip_sf%>%
                 st_difference(maritimes_network)%>% #what's not in the network
                  mutate(species = factor(species, levels = ordered_species))


p1 <- ggplot()+
  geom_sf(data=bioregion,fill=NA)+
  geom_sf(data=basemap_atlantic)+
  geom_sf(data=ssip_sf,fill="coral2",alpha=0.3,col=NA)+
  geom_sf(data=ssip_network_sf,fill="coral2")+
  geom_sf(data=maritimes_network,fill=NA,linewidth=0.5)+
  geom_sf(data=maritimes_network%>%filter(name=="Western/Emerald Banks Marine Refuge"),linewidth=1.2,fill=NA,col="black")+
  geom_text(data = ssip_network_sf,
            aes(x = Inf, y = -Inf, label = paste0("WEBMR % coverage: ",webca_prop)),
            hjust = 1, vjust = -0.5,  # Position in bottom right
            size = 3) +
  theme_bw()+
  coord_sf(expand=0,xlim=plotlims[c(1,3)],ylim=plotlims[c(2,4)])+
  facet_wrap(~species,nrow=3)+
  theme(strip.background = element_rect(fill="white"),
        axis.title = element_blank())

ggsave("output/ssip_all.png",p1,height=12,width=12,units="in",dpi=600)


#just the haddock 

haddock_area <- ssip_sf%>%
                filter(species=="Haddock")%>%
                mutate(area = as.numeric(st_area(geometry))/1000000)%>%
                pull(area)

haddock_area_network <- ssip_sf%>%
                        filter(species=="Haddock")%>%
                        st_intersection(maritimes_network%>%st_union())%>%
                        mutate(area = as.numeric(st_area(geometry))/1000000)%>%
                        pull(area)

ssip_network_haddock <- ssip_sf%>%
  filter(species=="Haddock")%>%
  st_intersection(maritimes_network)%>%
  dplyr::select(name,geometry)%>%
  group_by(name)%>%
  summarise(
    geometry = st_union(geometry),
    area = as.numeric(st_area(geometry))/1000000# Area in km²
  )%>%
  ungroup()%>%
  data.frame()%>%
  mutate(prop_area_total = round(area/haddock_area,3),
         prop_area_network = round(area/haddock_area_network,3))%>%
  dplyr::select(-geometry)%>%
  arrange(prop_area_total)%>%
  mutate(name = factor(name, levels=name),
         abbrev = case_when(name=="Western/Emerald Banks Marine Refuge" ~ "WEBMR",
                            name=="Fundian Channel-Browns Bank" ~ "Fundian",
                            TRUE ~ name),
         abbrev = factor(abbrev,levels = abbrev))




p2 <- ggplot()+
  geom_sf(data=bioregion,fill=NA)+
  geom_sf(data=basemap_atlantic)+
  geom_sf(data=ssip_sf%>%filter(species == "Haddock"),fill="coral2",alpha=0.3,col=NA)+
  geom_sf(data=ssip_network_sf%>%filter(species == "Haddock"),fill="coral2")+
  geom_sf(data=maritimes_network,fill=NA,linewidth=0.5)+
  geom_sf(data=maritimes_network%>%filter(name=="Western/Emerald Banks Marine Refuge"),linewidth=1.2,fill=NA,col="black")+
  theme_bw()+
  coord_sf(expand=0,xlim=plotlims[c(1,3)],ylim=plotlims[c(2,4)])+
  theme(strip.background = element_rect(fill="white"),
        axis.title = element_blank());p2

p2_bar <- ggplot(ssip_network_haddock,aes(y=abbrev))+
          geom_bar(stat="identity",aes(x=prop_area_total,fill=prop_area_total),col="black")+
          theme_bw()+
          scale_x_continuous(labels = percent, expand = c(0, 0.005))+ 
          scale_fill_viridis_b(labels = percent,n.breaks=5,option="D")+
          #geom_text(
          #   aes(x=prop_area_total, label = scales::percent(prop_area_total, accuracy = 0.1)), 
          #   hjust = 0,  # Align to the left side of the text
          #   vjust = -0.5,  # Center vertically
          #   angle = 270,  # Rotate 90 degrees
          #   nudge_x = 0.005,  # Slight nudge to prevent overlap
          #   size = 2.5
          # ) +
          theme(axis.title = element_blank(),
                legend.position = "none",
                axis.text.y = element_text(size = rel(0.9)));p2_bar
          
combo_haddock <- p2 + p2_bar + plot_layout(ncol=2, widths = c(0.7, 0.3),guides = 'collect') & theme(plot.margin = unit(c(0.05, 0.05, 0.05, 0.05), "lines"))
                  # inset_element(p2_bar, 
                  #               left = 0.65,   # Adjust these values to position the inset
                  #               bottom = 0.05, 
                  #               right = 0.95,  # Adjust these values to control size
                  #               top = 0.3)
ggsave("output/haddock_inset.png",combo_haddock,height=6,width=12,units="in",dpi=600)
