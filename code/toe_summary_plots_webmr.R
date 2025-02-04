#load libraries
library(dplyr)
library(tidyr)
library(sf)
library(ggplot2)
library(rnaturalearth)
library(scales)
library(patchwork)
library(ggspatial)
library(viridis)
library(ggimage)
library(rphylopic)

sf_use_s2(FALSE)

#Projections ------------
latlong <- "+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"

#focal species

focal_sp <- c("Melanogrammus aeglefinus","Gadus morhua","Hippoglossoides platessoides","Leucoraja ocellata")


focal_sp_df <- NULL
for(i in focal_sp){
  
  focal_sp_df  <- rbind(focal_sp_df,data.frame(species=i,uuid = get_uuid(i))) #playing with phylopic
  
  
}

#load the time of emergence summaries
load("c:/Users/stanleyr/Documents/Github/MAR_thermal_emerg/output/toe_summaries/all_toe_summaries.RData")

toe_dat <- toe_summaries%>%
                 filter(mod != "GFDL",
                        species %in% focal_sp,
                        NAME == "Western Emerald Bank Conservation Area")%>%
                 mutate(climate_proj = gsub("\\.","-",climate_proj))# the . was left in from the ensemble calculation


#calculate the total area for each species in the network
    agg_network_area <- toe_dat%>%
                      group_by(climate_proj,mod,species)%>% #slightly different for each mod and projection 
                      summarise(total_area=sum(cell_area))%>%
                      ungroup()%>%
                      data.frame()%>%
                      arrange(species)
    
    #so now we will add together grid cells that emerged in the same year
    agg_annual_toe <- toe_dat%>%
                      group_by(climate_proj,mod,species,ToE)%>%
                      summarise(area_lost=sum(cell_area))%>%
                      ungroup()%>%
                      data.frame()
                        

#extract the summaries of habitat loss (e.g., when a cell becomes too warm) ----------------------------

species <- unique(agg_annual_toe$species)
years <- 2015:2100
mods <- unique(agg_annual_toe$mod)
projs <- unique(agg_annual_toe$climate_proj)

    #Big 'for loop' == this could probably be done using 'do' in dplyr but this doesn't take long and is easy to follow. 
    habitat_loss <- NULL #will grow each loop
    for(i in species){
      message(paste0("Working on ",i))
      for(p in projs){
        for(m in mods){
         
            temp <- agg_annual_toe%>%
                        filter(species==i,climate_proj == p,mod==m,!is.na(ToE))%>%
                        rbind(.,data.frame(climate_proj=p,mod=m,species=i,ToE=setdiff(years,.$ToE),area_lost=0))%>% #add in the years that are missing as 0 loss years
                        arrange(ToE)%>% #sort them
                        mutate(cum_lost = cumsum(area_lost))%>%
                        left_join(.,agg_network_area)%>%# add in the total area for a given mod, projection and species
                        mutate(prop_lost = cum_lost/total_area)%>%
                        data.frame()%>%suppressMessages()
            
            habitat_loss <- rbind(habitat_loss,temp)
    
        } #end of 'm' mods loop
      } #end of 'p' climate_proj loop
    } #end of 'i' species loop

## Species by site 'emergence' based on a threshold of habitat loss summary data -----------------------

    #total area in each site that is occupied by each species
      agg_site_area <- toe_dat%>%
                          group_by(climate_proj,mod,species,NAME)%>% #slightly different for each mod and projection 
                          summarise(total_area=sum(cell_area))%>%
                          ungroup()%>%
                          data.frame()
      
      habitat_loss_site <- toe_dat%>%
                           mutate(ToE = ifelse(is.na(ToE),2500,ToE))%>% #2500 is a placeholder for 'NA' or 'not emerged'
                           group_by(climate_proj,mod,species,NAME,ToE)%>%
                           summarise(area_lost=sum(cell_area))%>%
                           ungroup()%>%
                           left_join(agg_site_area)%>% # add in the total area within each site. 
                           mutate(prop_area=area_lost/total_area)%>%
                           arrange(climate_proj,mod,species,NAME,ToE)%>% #make sure everything is ordered so that ToE's are sequential
                           group_by(climate_proj,mod,species,NAME)%>%
                           mutate(cum_sum=cumsum(prop_area))%>%
                           ungroup()%>%
                           data.frame()%>%
                           mutate(climate_proj_fact=factor(ifelse(climate_proj == "2-6","RCP 2.6","RCP 8.5")),
                                  mod = factor(mod, levels=c("AWI","IPSL","HAD","Ensemble")))%>%
                           left_join(.,focal_sp_df)


#construct the plot
webca_loss <- ggplot(habitat_loss_site, aes(x = ToE, y = cum_sum, color = mod,linetype = mod)) +
              geom_line(linewidth = 1,show.legend = FALSE) +
              geom_point(shape=21,aes(fill=mod),col="black",size=2.1, alpha = 0.6) +
              facet_grid(species~climate_proj_fact) +
              scale_linetype_manual(values = c("AWI" = "dashed", "IPSL" = "dashed", "HAD" = "dashed", "Ensemble" = "solid")) +
              scale_color_manual(values = c("AWI" = "coral2", "IPSL" = "aquamarine4", "HAD" = "cornflowerblue", "Ensemble" = "black")) +
              scale_fill_manual(values = c("AWI" = "coral2", "IPSL" = "aquamarine4", "HAD" = "cornflowerblue", "Ensemble" = "black")) +
              theme_bw() +
              labs(
                x = "Time of Emergence (ToE)",
                y = "% Habitat extent lost",
                fill = "Climate model"
              ) +
              theme(
                legend.position = "bottom",
                panel.grid.minor.y = element_blank(),
                strip.background = element_rect(fill="white"),
                strip.text.y = element_text(size = rel(0.7)),
                panel.grid.minor.x = element_line(color = "lightgrey", linetype = "dashed")
              ) +
              scale_size_continuous(range = c(1, 5)) +
              scale_x_continuous(
                limits = c(2015, 2100),
                breaks = seq(2020, 2100, by = 10)
              ) +
              scale_y_continuous(labels=percent)+
              # Add horizontal line at 1 to show complete habitat loss
              geom_hline(yintercept = 0.5, linetype = "dashed", color = "red", alpha = 0.5)+
              guides(
                color = guide_legend(override.aes = list(size = 5, shape = 21, stroke = 1)),
                linetype = "none",  # Hide linetype legend
                fill = guide_legend(override.aes = list(size = 5, shape = 21, stroke = 1))
              )

ggsave("output/webca_focal_habitat_loss.png",webca_loss,height=7.1,width=8,units="in",dpi=600)