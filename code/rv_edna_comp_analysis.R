## Comparision of rv data and eDNA data for webca

#load libraries ----------
library(sf)
library(tidyverse)
library(rnaturalearth)
library(rnaturalearthhires)
library(lubridate)
library(MarConsNetData)
library(vegan)
library(ggnewscale)

source("code/venn_bar.R")

#for plotting and data collation 
PhyloNames <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")

#Load the cleaned data
load("data/rv_taxonomy.RData")
load("data/edna_taxonomy.RData")

#Bar plots - the prettier venn diagrams

#trim the RV data to just those species detected in the paired samples
rv_comp_direct <- rv_formatted%>%
                  data.frame()%>%
                  filter(MISSION == gsub("-","",unique(edna_data$cruise)),
                         SETNO %in% unique(edna_data$setno))%>%
                  distinct(aphiaID,.keep_all=TRUE)%>%
                  rename(species_filter = latin)%>%
                  dplyr::select(all_of(c(PhyloNames,"aphiaID","species_filter")))%>%
                  mutate(method='rv')%>%
                  rbind(.,edna_data%>%
                          data.frame()%>%
                          distinct(aphiaID,.keep_all=TRUE)%>%
                          rename(species_filter = latin)%>%
                          dplyr::select(all_of(c(PhyloNames,"aphiaID","species_filter")))%>%
                          mutate(method='edna'))

bar_direct <- rv_comp_direct%>%
              group_by(aphiaID)%>%
              mutate(Status = case_when(
                sum(method == "rv") == 1 & sum(method == "edna") == 0 ~ "Traditional Only",
                sum(method == "rv") == 0 & sum(method == "edna") == 1 ~ "eDNA Only",
                TRUE ~ "Shared"
              ))%>%
              ungroup()%>%
              distinct(aphiaID,.keep_all=TRUE)%>%
              group_by(Phylum,Status)%>%
              summarise(Count=n())%>%
              ungroup()%>%
              mutate(Status = factor(Status,levels=c("Shared","Traditional Only","eDNA Only")))%>%
              data.frame()

bar_plot_direct <- venn_bar(bar_direct)

ggsave("output/bar_plot_direct.png",bar_plot_direct,height=5,width=7,units="in",dpi=300)
              

#trim the RV data to just the species detected near WEBMR in the cruise where eDNA was collected 
rv_comp_survey <- rv_formatted%>%
                  data.frame()%>%
                  filter(MISSION == gsub("-","",unique(edna_data$cruise)))%>%
                  rename(species_filter = latin)%>%
                  dplyr::select(all_of(c(PhyloNames,"aphiaID","species_filter")))%>%
                  mutate(method='rv')%>%
                  rbind(.,edna_data%>%
                          data.frame()%>%
                          rename(species_filter = latin)%>%
                          dplyr::select(all_of(c(PhyloNames,"aphiaID","species_filter")))%>%
                          mutate(method='edna'))

#comparison of rv data to all the surveys near WEBMR since 1970
rv_comp_all <- rv_formatted%>%
               data.frame()%>%
                distinct(aphiaID,.keep_all=TRUE)%>%
                rename(species_filter = latin)%>%
                dplyr::select(all_of(c(PhyloNames,"aphiaID","species_filter")))%>%
                mutate(method='rv')%>%
                rbind(.,edna_data%>%
                        data.frame()%>%
                        distinct(aphiaID,.keep_all=TRUE)%>%
                        rename(species_filter = latin)%>%
                        dplyr::select(all_of(c(PhyloNames,"aphiaID","species_filter")))%>%
                        mutate(method='edna'))

bar_all <- rv_comp_all%>%
            group_by(aphiaID)%>%
            mutate(Status = case_when(
              sum(method == "rv") == 1 & sum(method == "edna") == 0 ~ "Traditional Only",
              sum(method == "rv") == 0 & sum(method == "edna") == 1 ~ "eDNA Only",
              TRUE ~ "Shared"
            ))%>%
            ungroup()%>%
            distinct(aphiaID,.keep_all=TRUE)%>%
            group_by(Phylum,Status)%>%
            summarise(Count=n())%>%
            ungroup()%>%
            mutate(Status = factor(Status,levels=c("Shared","Traditional Only","eDNA Only")))%>%
            data.frame()

bar_plot_all <- venn_bar(bar_all)

ggsave("output/bar_plot_all.png",bar_plot_all,height=5,width=7,units="in",dpi=300)

#Species accumulation curve comparisons 

edna_species_df <- read.csv("data/2020eDNA_COI12S16S_FinalTaxaTable.csv")%>%
                  dplyr::select(c("OTU_ID",edna_data$stations))%>%
                  rowwise() %>%
                  mutate(otu_count = sum(c_across(2:8)))%>%
                  data.frame()%>%
                  filter(otu_count>0)%>%
                  dplyr::select(-otu_count)%>%
                  pivot_longer(
                    cols = starts_with("eDNA20"),  # Select site columns
                    names_to = "Site",                # Name for the new 'site' column
                    values_to = "Count"               # Name for the new values column
                  ) %>%
                  pivot_wider(
                    names_from = OTU_ID,  # Use OTU_ID as the new column names
                    values_from = Count   # Fill values from the Count column
                  )%>%
                  column_to_rownames(var = "Site")

rv_species_df <- rv_formatted%>%
                 data.frame()%>%
                 filter(MISSION == gsub("-","",unique(edna_data$cruise)),
                       SETNO %in% unique(edna_data$setno))%>%
                 mutate(stations = paste(MISSION,SETNO,sep="-"))%>%
                 dplyr::select(stations,latin,std_wgt)%>%
                  pivot_wider(
                    names_from = latin,           # Species names become columns
                    values_from = std_wgt,        # Biomass fills the values
                    values_fill = 0               # Fill missing values with 0
                  ) %>%
                  column_to_rownames(var = "stations")

#get the data for plotting 
edna_specpool <- species_pool_estimate(edna_species_df)
edna_plot_data <- (edna_specpool$plot_data)%>%mutate(method="eDNA")
edna_accum_data <- (edna_specpool$accum_data)%>%mutate(method="eDNA")

rv_specpool <- species_pool_estimate(rv_species_df)
rv_plot_data <- (rv_specpool$plot_data)%>%mutate(method="RV survey")
rv_accum_data <- (rv_specpool$accum_data)%>%mutate(method="RV survey")

plot_data <- rbind(edna_plot_data,rv_plot_data)
accum_data <- rbind(edna_accum_data,rv_accum_data)

y_intercepts <- plot_data%>%
                group_by(method)%>%
                summarise(y_intercept = unique(asymptote))%>%
                ungroup()


p1 <- ggplot(plot_data, aes(x = Sites)) +
  geom_line(aes(y = Fitted, color = Observed), linetype="dotted",size = 1) +  # Fitted curve
  geom_line(data=plot_data%>%filter(Observed == "Observed"),aes(y = Fitted, color = Observed),size=1)+
  geom_point(data = accum_data, aes(y = Richness), shape=21,fill = "cornflowerblue", size = 2.5) +  # Observed richness
  geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.2, fill = "red") +  # Error bounds
  scale_color_manual(values = c("Observed" = "blue", "Extrapolated" = "red")) +
  labs(
    x = "Number of Sites",
    y = "Species Richness",
    color = ""
  ) +
  theme_bw() +
  facet_wrap(~method,ncol=2)+
  geom_hline(data = y_intercepts, aes(yintercept = y_intercept),colour="black", linetype = "longdash") +
  geom_hline(data = y_intercepts%>%mutate(method=rev(method)), aes(yintercept = y_intercept),colour="grey75", linetype = "longdash")+
  theme(strip.background = element_rect(fill="white"),
        legend.position = "none")

ggsave("output/speccum_plots.png",p1,height=5,width=5,units="in",dpi=300)



