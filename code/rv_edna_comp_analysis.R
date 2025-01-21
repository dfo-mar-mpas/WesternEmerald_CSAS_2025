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
source("code/species_pool_estimate.R")

common_name <- function(x) {
  
  require(worrms)
  
  # Attempt to get common name records
  cn <- tryCatch(
    wm_common_id(x),
    error = function(e) NULL  # Return NULL if there's an error
  )
  
  # Check if there are any records
  if (is.null(cn) || nrow(cn) == 0) {
    return(NA)  # Return NA if no records
  }
  
  # Filter for English language and format the common name
  common <- cn %>%
    filter(language == "English") %>%
    pull(vernacular) %>%
    sub("^(\\w)(\\w*)", "\\U\\1\\L\\2", ., perl = TRUE)
  
  # Return the first common name, or NA if none match the filter
  if (length(common) > 1) {
    return(common[1])
  } else {
    return(common)
  }
}

#for plotting and data collation 
PhyloNames <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")

#Load the cleaned data
load("data/rv_taxonomy.RData")
load("data/edna_taxonomy.RData")

#Bar plots - the prettier venn diagrams

#there is a mismatch for yellowtail flounder with some records in the rv_formatted with the Species name (Limanda ferruginea) instead of the the now accepted myzopsetta ferruginea

edna_data <- data.frame(edna_data)

rows_to_update <- which(edna_data$Species == "Limanda ferruginea")

edna_data[rows_to_update,c("latin",PhyloNames,"aphiaID")] <- rv_formatted%>%
                                                             data.frame()%>%
                                                             filter(latin=="myzopsetta ferruginea")%>%dplyr::select(latin,all_of(PhyloNames),aphiaID)%>%
                                                             slice(1)%>%as.vector()

edna_data[rows_to_update,"OTU_ID"] <- "Myzopsetta ferruginea"


                  
#filter out just the species               
rv_comp_direct <- rv_formatted%>%
                  data.frame()%>%
                  filter(MISSION == gsub("-","",unique(edna_data$cruise)),
                         SETNO %in% unique(edna_data$setno))%>%
                  distinct(aphiaID,.keep_all=TRUE)%>%
                  filter(!is.na(Species))%>%
                  rename(species_filter = latin)%>%
                  dplyr::select(all_of(c(PhyloNames,"aphiaID","species_filter")))%>%
                  mutate(method='rv')%>%
                  rbind(.,edna_data%>%
                          data.frame()%>%
                          distinct(aphiaID,.keep_all=TRUE)%>%
                          filter(!is.na(Species))%>%
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
               filter(!is.na(Species))%>%
               rename(species_filter = latin)%>%
               dplyr::select(all_of(c(PhyloNames,"aphiaID","species_filter")))%>%
               mutate(method='rv')%>%
               rbind(.,edna_data%>%
                       data.frame()%>%
                       distinct(aphiaID,.keep_all=TRUE)%>%
                       filter(!is.na(Species))%>%
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

#doesn't render well (titles can't be aligned vertically with the legend). Will use Illustrator instead
# bar_combo <- (bar_plot_direct + labs(title = "2020 paired RV-eDNA comparison")) + 
#              (bar_plot_all + labs(title = "1970-2024 RV-eDNA comparison") + theme(legend.position="none"))+
#               plot_layout(nrow=2)
# 
# ggsave("output/bar_combo.png",bar_combo,height=8,width=7,units = "in",dpi=300)

#Species level accumulation curve comparisons 

edna_meta <- edna_data%>%
             dplyr::select(stations,cruise,setno)%>%
             mutate(MISSION = gsub("-","",cruise),
                    sample_site = paste(MISSION,setno,sep="-"))%>%
             distinct(sample_site,.keep_all=TRUE)%>%
             dplyr::select(stations,sample_site)%>%
             rename(Site=stations)
             

edna_species_df <- read.csv("data/2020eDNA_COI12S16S_FinalTaxaTable.csv")%>%
                  dplyr::select(c("OTU_ID",edna_data$stations))%>%
                  left_join(.,edna_data%>%
                              dplyr::select(OTU_ID,Species,aphiaID)%>%
                              distinct(aphiaID,.keep_all=TRUE))%>%
                  rowwise() %>%
                  mutate(otu_count = sum(c_across(2:8)))%>%
                  data.frame()%>%
                  filter(otu_count>0,
                         !is.na(Species))%>%
                  dplyr::select(-otu_count,-OTU_ID,-Species)%>%
                  rename(OTU_ID = aphiaID)%>%
                  pivot_longer(
                    cols = starts_with("eDNA20"),  # Select site columns
                    names_to = "Site",                # Name for the new 'site' column
                    values_to = "Count"               # Name for the new values column
                  ) %>%
                  pivot_wider(
                    names_from = OTU_ID,  # Use OTU_ID as the new column names
                    values_from = Count   # Fill values from the Count column
                  )%>%
                  left_join(.,edna_meta)%>%
                  dplyr::select(-Site)%>%
                  column_to_rownames(var = "sample_site") #to match with the RV data


rv_species_df <- rv_formatted%>%
                 data.frame()%>%
                 filter(MISSION == gsub("-","",unique(edna_data$cruise)),
                        SETNO %in% unique(edna_data$setno))%>%
                 mutate(stations = paste(MISSION,SETNO,sep="-"))%>%
                 filter(!is.na(Species))%>%
                 dplyr::select(stations,aphiaID,std_wgt)%>%
                  pivot_wider(
                    names_from = aphiaID,           # Species names become columns
                    values_from = std_wgt,        # Biomass fills the values
                    values_fill = 0               # Fill missing values with 0
                  ) %>%
                  column_to_rownames(var = "stations")

rv_species_df_all <- rv_formatted%>%
                    data.frame()%>%
                    mutate(stations = paste(MISSION,SETNO,sep="-"))%>%
                    filter(!is.na(Species))%>%
                    dplyr::select(stations,latin,std_wgt)%>%
                    pivot_wider(
                      names_from = latin,           # Species names become columns
                      values_from = std_wgt,        # Biomass fills the values
                      values_fill = 0               # Fill missing values with 0
                    ) %>%
                    column_to_rownames(var = "stations")

#combination eDNA + RV 2020
combo_matrix <- cbind(
  
  edna_species_df%>%
    dplyr::select(intersect(colnames(edna_species_df),colnames(rv_species_df))), #shared
  
  edna_species_df%>%
    dplyr::select(setdiff(colnames(edna_species_df),colnames(rv_species_df))), #eDNA only
  
  rv_species_df[match(rownames(edna_species_df),rownames(rv_species_df)),]%>%
    dplyr::select(setdiff(colnames(rv_species_df),colnames(edna_species_df))) #RV only
  
)

# Replace NAs with 0 (indicating absence)
combo_matrix[is.na(combo_matrix)] <- 0

# Convert non-zero values to 1 (binary presence/absence matrix)
combo_matrix[combo_matrix > 0] <- 1

#calculate species accumulation curves

    #eDNA only (edna within WEBMR)
    edna_specpool <- species_pool_estimate(edna_species_df)
    edna_plot_data <- (edna_specpool$plot_data)%>%mutate(method="eDNA")
    edna_accum_data <- (edna_specpool$accum_data)%>%mutate(method="eDNA")
    
    #rv only (2020 paried with eDNA)
    rv_specpool <- species_pool_estimate(rv_species_df)
    rv_plot_data <- (rv_specpool$plot_data)%>%mutate(method="RV survey")
    rv_accum_data <- (rv_specpool$accum_data)%>%mutate(method="RV survey")
    
    #rv all only - 1970-2024
    rv_specpool_all <- species_pool_estimate(rv_species_df_all)
    rv_plot_data_all <- (rv_specpool_all$plot_data)%>%mutate(method="RV survey")
    rv_accum_data_all <- (rv_specpool_all$accum_data)%>%mutate(method="RV survey")
    
    #combined eDNA + RV (2020)
    combo_specpool <- species_pool_estimate(combo_matrix)
    combo_plot_data <- (combo_specpool$plot_data)%>%mutate(method="eDNA + RV")
    combo_accum_data <- (combo_specpool$accum_data)%>%mutate(method="eDNA + RV")
    
    #combination data
    plot_data <- rbind(edna_plot_data,rv_plot_data)
    accum_data <- rbind(edna_accum_data,rv_accum_data)

#plot data assembly
    plot_data_comp <- plot_data%>%
                      mutate(method=gsub("RV survey","RV paired",method))%>%
                      rbind(.,rv_plot_data_all,combo_plot_data)%>%
                      mutate(method=factor(method,levels=c("RV survey","RV paired","eDNA","eDNA + RV")))
    
    accum_data_comp <- accum_data%>%
                       mutate(method=gsub("RV survey","RV paired",method))%>%
                       rbind(.,rv_accum_data_all,combo_accum_data)%>%
                       mutate(method=factor(method,levels=c("RV survey","RV paired","eDNA","eDNA + RV")))
    
    y_intercepts <- plot_data%>%
                    group_by(method)%>%
                    summarise(y_intercept = unique(asymptote))%>%
                    ungroup()
    
    y_intercepts_all <- plot_data_comp%>%
                        group_by(method)%>%
                        summarise(y_intercept = unique(asymptote))%>%
                        ungroup()%>%
                        mutate(method=factor(method,levels=c("RV survey","RV paired","eDNA","eDNA + RV")))

#make plots
    
    #direct comparison
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

    #just RV survey (total)
    p2 <- ggplot(rv_plot_data_all, aes(x = Sites)) +
      geom_line(aes(y = Fitted, color = Observed), linetype="dotted",size = 1) +  # Fitted curve
      geom_line(data=rv_plot_data_all%>%filter(Observed == "Observed"),aes(y = Fitted, color = Observed),size=1)+
      geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.2, fill = "red") +  # Error bounds
      scale_color_manual(values = c("Observed" = "blue", "Extrapolated" = "red")) +
      labs(
        x = "Number of Sites",
        y = "Species Richness",
        color = ""
      ) +
      theme_bw() +
      geom_hline(yintercept = unique(rv_plot_data_all$asymptote),colour="black", linetype = "longdash") +
      theme(strip.background = element_rect(fill="white"),
            legend.position = "none")

ggsave("output/speccum_plots_all.png",p2,height=5,width=5,units="in",dpi=300)

#rv survey (total) vs 2020 eDNA
  p3 <- ggplot(plot_data_comp, aes(x = Sites)) +
    geom_line(aes(y = Fitted, color = Observed), linetype="dotted",linewidth = 1) +  # Fitted curve
    geom_line(data=plot_data_comp%>%filter(Observed == "Observed"),aes(y = Fitted, color = Observed),size=1)+
    geom_point(data = accum_data_comp%>%filter(method !="RV survey"), aes(y = Richness), shape=21,fill = "cornflowerblue", size = 2.5) +  # Observed richness
    geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.2, fill = "red") +  # Error bounds
    scale_color_manual(values = c("Observed" = "blue", "Extrapolated" = "red")) +
    labs(
      x = "Number of Sites",
      y = "Species Richness",
      color = ""
    ) +
    theme_bw() +
    facet_wrap(~method,ncol=4,scales="free_x")+
    geom_hline(data = y_intercepts_all, aes(yintercept = y_intercept),colour="black", linetype = "longdash") +
    #geom_hline(data = y_intercepts_all%>%mutate(method=rev(method)), aes(yintercept = y_intercept),colour="grey75", linetype = "longdash")+
    theme(strip.background = element_rect(fill="white"),
          legend.position = "none")

ggsave("output/speccum_plots_all.png",p3,height=5,width=7,units="in",dpi=300)


##create the summary table with all 90 species detected

PhyloNames <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")

rv_data_table <- rv_formatted%>%
                 data.frame()%>%
                 filter(MISSION == gsub("-","",unique(edna_data$cruise)),
                         SETNO %in% unique(edna_data$setno),
                        !is.na(Species))%>%
                 distinct(latin,.keep_all=TRUE)%>%
                 dplyr::select(all_of(PhyloNames),aphiaID)

edna_data_table <- edna_data%>%
                   distinct(latin,.keep_all = TRUE)%>%
                   dplyr::select(all_of(PhyloNames),aphiaID)

species_table <- rv_data_table%>%
                 filter(aphiaID %in% setdiff(rv_data_table$aphiaID,edna_data_table$aphiaID))%>%
                 mutate(type = "RV")%>%
                 rbind(.,
                       rv_data_table%>%
                         filter(aphiaID %in% intersect(rv_data_table$aphiaID,edna_data_table$aphiaID))%>%
                         mutate(type = "Shared")
                       )%>%
                 rbind(.,
                       edna_data_table%>%
                         filter(aphiaID %in% setdiff(edna_data_table$aphiaID,rv_data_table$aphiaID))%>%
                         mutate(type = "eDNA"))%>%
                 mutate(type=factor(type,levels=c("RV","Shared","eDNA")))%>%
                 arrange(type,Phylum,Class,Order,Family,Genus,Species)

write.csv(species_table,"output/comparison_species_table.csv",row.names=FALSE)

#common WEBCA species

grouped_data <- edna_data%>%
                data.frame()%>%
                filter(!is.na(Species))%>%
                dplyr::select(all_of(PhyloNames),aphiaID,latin)%>%
                rbind(.,
                rv_formatted%>%
                  data.frame()%>%
                  filter(!is.na(Species))%>%
                  dplyr::select(all_of(PhyloNames),aphiaID,latin)
                )%>%
                distinct(latin,.keep_all = TRUE)
                
sets_per_year <- rv_formatted %>%
                 data.frame()%>%
                 filter(!is.na(Species))%>%
                 group_by(YEAR) %>%
                 summarise(total_sets = n_distinct(SETNO))

# Summarize total count and weight by species and year, normalized by the number of sets per year
fish_summary <- rv_formatted %>%
  data.frame()%>%
  filter(!is.na(Species),
         Phylum == "Chordata")%>%
  group_by(YEAR, latin) %>%
  summarise(
    total_count = sum(std_count, na.rm = TRUE),  # Total count for the species in the year
    total_weight = sum(std_wgt, na.rm = TRUE)   # Total weight for the species in the year
  ) %>%
  left_join(sets_per_year) %>%      # Join with the sets_per_year data
  mutate(
    count_per_set = total_count / total_sets,    # Normalize count by number of sets
    weight_per_set = total_weight / total_sets, # Normalize weight by number of sets
    count_rank = rank(-count_per_set),          # Rank species by normalized count
    weight_rank = rank(-weight_per_set)         # Rank species by normalized weight
  ) %>%
  arrange(YEAR, weight_rank) 
                         
          
top_10_fish <- fish_summary %>%
  filter(count_rank <= 10 | weight_rank <= 10) %>%
  mutate(
    in_top_10_count = ifelse(count_rank <= 10, 1, 0),  # 1 if in top 10 by count
    in_top_10_weight = ifelse(weight_rank <= 10, 1, 0) # 1 if in top 10 by weight
  )

# Summarize the number of occurrences in the top 10 for each aphiaID
top_10_summary <- top_10_fish %>%
  group_by(latin) %>%
  summarise(
    top_10_count_occurrences = sum(in_top_10_count),   # Times in top 10 by count
    top_10_weight_occurrences = sum(in_top_10_weight) # Times in top 10 by weight
  ) %>%
  arrange(desc(top_10_count_occurrences), desc(top_10_weight_occurrences))%>%
  left_join(.,grouped_data)%>%# Sort by most frequent
  mutate(rank = 1:n())
# View the summary

#missed species
top_10_summary%>%
  slice(1:10)%>%
  filter(aphiaID %in% setdiff(aphiaID,edna_data%>%pull(aphiaID)%>%unique()))%>%
  dplyr::select(Species,aphiaID,rank)


rv_comp_direct%>%filter(aphiaID %in% setdiff(top_10_summary%>%slice(1:10)%>%pull(aphiaID),edna_data%>%pull(aphiaID)%>%unique()))

top_10_df <- top_10_summary%>%
             slice(1:10)%>%
             rowwise()%>%
             mutate(Common = common_name(aphiaID),
                    type = ifelse(aphiaID %in% setdiff(aphiaID,edna_data%>%pull(aphiaID)%>%unique()),"RV only","Shared"))%>%
             dplyr::select(Common,Species,aphiaID,rank,type)
             

#Manual fixes
top_10_df[top_10_df$Species=="Myzopsetta ferruginea","Common"] <- "Yellowtail flounder"
top_10_df[top_10_df$Species=="Pollachius virens","Common"] <- "Pollock"
top_10_df[top_10_df$Species=="Pseudopleuronectes americanus","Common"] <- "Winter flounder"

write.csv(top_10_df,file="output/top_10_fish.csv",row.names = FALSE)


  tt <- read.csv("data/2020eDNA_COI12S16S_FinalTaxaTable.csv")%>%
        rowwise() %>%
        mutate(otu_count = sum(c_across(starts_with("EDNA"))))%>%
        data.frame()%>%
        #filter(otu_count>0)%>%
        gather(key="station","count",2:8)%>%
        mutate(latin=gsub(" g\\.$","",OTU_ID),
               latin=gsub(" cf.","",latin),
               latin=gsub(" f\\.$","",latin))%>%
        filter(grepl("gadus",tolower(latin))| grepl("pollac",tolower(latin)))%>%
        dplyr::select(-OTU_ID,-station,-count,-otu_count)%>%
          pivot_longer(
            cols = starts_with("eDNA20"),  # Select site columns
            names_to = "Site",                # Name for the new 'site' column
            values_to = "Count"               # Name for the new values column
          ) %>%
        filter(Count>0)%>%
        mutate(id=paste(latin,Site,Count,sep="-"))%>% #this code is crudgy AF but works
        distinct(id,.keep_all=T)
  
  #was gadus observed in webca
  tt%>%filter(Site %in% unique(edna_data$stations))
  

## general summary stats

#total distance
rv_formatted%>%mutate(id=paste(MISSION,SETNO,sep="-"))%>%distinct(id,.keep_all=TRUE)%>%pull(DIST)%>%sum()
