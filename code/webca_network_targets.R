#Design Target plot

#load libraries
library(tidyverse)
library(viridis)
library(scales)

#load the target summary from the app

abbreviation_table <- read.csv("data/abbrev_table_network_features.csv")%>%
                      mutate(feature = trimws(feature))%>%
                      dplyr::select(feature,abbrev)


webca_targets <- read.csv("data/targetstable.csv")%>%
                 filter(Achieved > 0,
                        type != "Commercial fishery landing")%>%
                 rename(feature=Variable)%>%
                 mutate(feature = trimws(feature))%>%
                 left_join(.,abbreviation_table)%>%
                 mutate(type = gsub(" classification","",type),
                        type_f = factor(type,levels=c("Biodiversity hotspot","Depleted species","Biogenic habitat",
                                                    "Seabird functional group","Invertebrate functional group",
                                                    "Fish functional group","Scope for growth","Natural disturbance",
                                                    "Geomorphic unit","Biophysical unit")),
                        abbrev = ifelse(feature=="Medium Benthic Benthivores (East)" & type=="Invertebrate functional group","Inv Benthic Med 4VW",abbrev), #overlap in names
                        abbrev = ifelse(feature=="Small Benthic Benthivores (East)" & type=="Invertebrate functional group","Inv Benthic Sm 4VW",abbrev))%>%
                  arrange(type_f,Percent.to.Goal)     

# Create blank rows for each type_f level
blank_rows <- webca_targets %>%
              group_by(type_f) %>%
              summarise(across(everything(), ~NA), .groups = "drop") 

# Bind and arrange
webca_targets_padded1 <- webca_targets %>%
                        bind_rows(blank_rows) %>%
                        arrange(type_f, desc(Percent.to.Goal))

webca_targets_padded <- webca_targets_padded1 %>%
                        slice(1) %>%
                        mutate(type_f = "Biodiversity hotspot") %>%
                        mutate(across(-feature, ~NA))%>%
                        rbind(.,webca_targets_padded1)%>%
                        mutate(id = n():1,
                               abbrev = ifelse(is.na(abbrev), "", abbrev))
          

# Create named vector for labels
id_labels <- setNames(webca_targets_padded$abbrev, webca_targets_padded$id)

# Plot
webca_plot <- ggplot(webca_targets_padded, aes(y = factor(id), x = Percent.to.Goal/100, fill = type)) +
  geom_vline(xintercept=0.5,lty=2,col="grey70")+
  geom_bar(stat = "identity", col = "black") +
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        legend.position = "none",
        axis.text.y = element_text(size = rel(0.7))) +
  scale_y_discrete(labels = id_labels)+
  scale_x_continuous(expand=(c(0.01,0)),labels=percent)+
  labs(y="",
       x="% of design target")+
  scale_fill_manual(values=c("darkmagenta","darkgoldenrod1","firebrick4",
                             "deepskyblue4","darkslateblue","cadetblue2","darkgreen",
                             "darkorange2","deeppink2","seagreen3"))

ggsave("output/webca_design_targets.png",webca_plot,width=7,height=12,dpi = 1200)
  



## Network analysis with the WEBMR addition 
network_targets <- read.csv("data/targetstable_all.csv")%>%
                    filter(Achieved > 0,
                           type != "Commercial fishery landing")%>%
                    rename(feature=Variable)%>%
                    mutate(feature = trimws(feature))%>%
                    left_join(.,abbreviation_table)%>%
                    mutate(type = gsub(" classification","",type),
                           type_f = factor(type,levels=c("Biodiversity hotspot","Depleted species","Biogenic habitat",
                                                         "Seabird functional group","Invertebrate functional group",
                                                         "Fish functional group","Scope for growth","Natural disturbance",
                                                         "Geomorphic unit","Biophysical unit")),
                           abbrev = ifelse(feature=="Medium Benthic Benthivores (East)" & type=="Invertebrate functional group","Inv Benthic Med 4VW",abbrev), #overlap in names
                           abbrev = ifelse(feature=="Small Benthic Benthivores (East)" & type=="Invertebrate functional group","Inv Benthic Sm 4VW",abbrev))%>%
                    arrange(type_f,Percent.to.Goal)     

network_targets_sum <- network_targets%>%
                       group_by(feature)%>%
                       summarise(Achieved = sum(Achieved,na.rm=T))%>%
                       ungroup()%>%
                       left_join(.,network_targets%>%
                                   distinct(feature,.keep_all=TRUE)%>%
                                   dplyr::select(feature,type,filter,Minimum.Target,abbrev,type_f))%>%
                        mutate( #translate with help from Gemeni
                          abbrev_fr = case_when(
                            # Geomorphology & Habitat Features
                            trimws(abbrev) == "Abyssal Plain"                 ~ "Plaine abyssale",
                            trimws(abbrev) == "Continental Rise"              ~ "Glacis continental",
                            trimws(abbrev) == "Laurentian Slope"              ~ "Pente laurentienne",
                            trimws(abbrev) == "Shelf Bank"                    ~ "Banc du plateau",
                            trimws(abbrev) == "Shelf Basin"                   ~ "Bassin du plateau",
                            trimws(abbrev) == "Shelf Channel"                 ~ "Chenal du plateau",
                            trimws(abbrev) == "Shelf Flat"                    ~ "Plat du plateau",
                            trimws(abbrev) == "Shelp Top"                     ~ "Haut du plateau",
                            trimws(abbrev) == "Shelf Top Bank"                ~ "Banc du haut du plateau",
                            trimws(abbrev) == "Shelf Top Basin"               ~ "Bassin du haut du plateau",
                            trimws(abbrev) == "Slope"                         ~ "Pente",
                            trimws(abbrev) == "Slope Channel"                 ~ "Chenal de pente",
                            trimws(abbrev) == "Slope-Rise-Abyss"              ~ "Pente-Glacis-Abysse",
                            
                            # Regions & Locations
                            trimws(abbrev) == "Baccro Lehave"                 ~ "Baccro LaHave",
                            trimws(abbrev) == "ESS"                           ~ "PNEO",  # Plateau néo-écossais oriental
                            trimws(abbrev) == "GOM"                           ~ "GDM",   # Golfe du Maine
                            trimws(abbrev) == "Lehave Emerald"                ~ "LaHave Émeraude",
                            trimws(abbrev) == "Western Sable"                 ~ "Sable Ouest",
                            
                            # Rank / Target Priority Categories
                            trimws(abbrev) == "V. High"                       ~ "Très élevé",
                            trimws(abbrev) == "High"                          ~ "Élevé",
                            trimws(abbrev) == "Moderate"                      ~ "Modéré",
                            trimws(abbrev) == "Medium"                        ~ "Moyen",
                            trimws(abbrev) == "Low"                           ~ "Faible",
                            trimws(abbrev) == "V. Low"                        ~ "Très faible",
                            
                            # Species & Critical Habitat (CH)
                            trimws(abbrev) == "A Cod 4Vn"                     ~ "Morue de l'Atlantique 4Vn",
                            trimws(abbrev) == "A Cod 4VsW"                    ~ "Morue de l'Atlantique 4VsW",
                            trimws(abbrev) == "A Cod 4X"                      ~ "Morue de l'Atlantique 4X",
                            trimws(abbrev) == "Am Plaice 4VW"                 ~ "Plie canadienne 4VW",
                            trimws(abbrev) == "Am Plaice 4X"                  ~ "Plie canadienne 4X",
                            trimws(abbrev) == "Atl Wolffish"                  ~ "Loup atlantique",
                            trimws(abbrev) == "Cusk"                          ~ "Brosme",
                            trimws(abbrev) == "N Bottlenose CH"               ~ "HC Baleine à bec commune", # Habitat critical
                            trimws(abbrev) == "Ocean Pout"                    ~ "Loup d'Amérique",
                            trimws(abbrev) == "Redfish"                       ~ "Sébaste",
                            trimws(abbrev) == "Right Whale CH"                ~ "HC Baleine noire",
                            trimws(abbrev) == "Roughhead Gren"                ~ "Grenadier berglax",
                            trimws(abbrev) == "Roundnose Gren"                ~ "Grenadier de roche",
                            trimws(abbrev) == "Sand Dollar"                   ~ "Dollar de sable",
                            trimws(abbrev) == "Sm. Skate 4VsW"                ~ "Raie tachetée 4VsW",
                            trimws(abbrev) == "Sm. Skate 4X"                   ~ "Raie tachetée 4X",
                            trimws(abbrev) == "Spiny Dogfish"                 ~ "Aiguillat commun",
                            trimws(abbrev) == "Th. Skate 4VSW"                ~ "Raie épineuse 4VsW",
                            trimws(abbrev) == "Th. Skate 4X"                   ~ "Raie épineuse 4X",
                            trimws(abbrev) == "White Hake 4VW"                ~ "Merluche blanche 4VW",
                            trimws(abbrev) == "White Hake 4X"                 ~ "Merluche blanche 4X",
                            trimws(abbrev) == "Wint. Skate 4VsW"              ~ "Raie tachetée 4VsW", # Winter Skate -> Raie tachetée (ou Raie d'hiver)
                            
                            # Benthic / Invertebrates / Corals / Sponges
                            trimws(abbrev) == "Boltenia"                      ~ "Boltenia",
                            trimws(abbrev) == "Lg Gorgonian KDE"              ~ "Gorgones t.g. CDE", # Colonnes / Estimation de la densité par noyau
                            trimws(abbrev) == "Lg Gorgonian SDM"              ~ "Gorgones t.g. MDH", # Modèle de distribution de l'habitat
                            trimws(abbrev) == "Sea Pens KDE"                  ~ "Platules CDE",
                            trimws(abbrev) == "Sea Pens SDM"                  ~ "Platules MDH",
                            trimws(abbrev) == "Sm Gorgonian SDM"              ~ "Pty. gorgones MDH",
                            trimws(abbrev) == "Soft Coral"                    ~ "Mousse de mer / Corail mou",
                            trimws(abbrev) == "Sponges Gen"                   ~ "Éponges gén.",
                            trimws(abbrev) == "Vazella"                       ~ "Vazella",
                            
                            # Diversity Measures
                            trimws(abbrev) == "Fish Diversity 4VW"            ~ "Diversité des poissons 4VW",
                            trimws(abbrev) == "Fish Diversity 4X"             ~ "Diversité des poissons 4X",
                            trimws(abbrev) == "Invert Diversity 4VW"          ~ "Diversité des invertébrés 4VW",
                            trimws(abbrev) == "Invert Diversity 4X"           ~ "Diversité des invertébrés 4X",
                            trimws(abbrev) == "Larval Diversity"              ~ "Diversité larvaire",
                            trimws(abbrev) == "Small Fish Diversity Fish 4VW" ~ "Div. pty. poissons 4VW",
                            trimws(abbrev) == "Small Fish Diversity Fish 4X"  ~ "Div. pty. poissons 4X",
                            trimws(abbrev) == "Small Invert Diversity 4VW"    ~ "Div. pty. invertébrés 4VW",
                            trimws(abbrev) == "Small Invert Diversity 4X"     ~ "Div. pty. invertébrés 4X",
                            
                            # Guilds & Functional Groups (Size / Habitat / Feeding)
                            trimws(abbrev) == "Benthic Lg 4VW"                ~ "Benth. G 4VW",
                            trimws(abbrev) == "Benthic Lg 4X"                 ~ "Benth. G 4X",
                            trimws(abbrev) == "Benthic Med 4X"                ~ "Benth. M 4X",
                            trimws(abbrev) == "Inv Ben Colon FF 4VW"          ~ "Invert Benth Colon FF 4VW",
                            trimws(abbrev) == "Inv Ben n-Colon FF 4VW"        ~ "Invert Benth non-Colon FF 4VW",
                            trimws(abbrev) == "Inv Ben n-Colon FF 4X"         ~ "Invert Benth non-Colon FF 4X",
                            trimws(abbrev) == "Inv Benthic Med 4VW"           ~ "Invert Benth M 4VW",
                            trimws(abbrev) == "Inv Benthic Sm 4VW"            ~ "Invert Benth P 4VW",
                            trimws(abbrev) == "Inv Detritivore 4VW"           ~ "Invert Détritivore 4VW",
                            trimws(abbrev) == "Inv Detritivore 4X"            ~ "Invert Détritivore 4X",
                            trimws(abbrev) == "Inv Zoopisc SML 4VW"           ~ "Invert Zoopisc PET 4VW",
                            trimws(abbrev) == "Inv Zoopisc SML 4X"            ~ "Invert Zoopisc PET 4X",
                            trimws(abbrev) == "Benthic SM 4X"                 ~ "Benth. P 4X",
                            trimws(abbrev) == "Pisc Benthic Lg 4VW"           ~ "Pisc Benth. G 4VW",
                            trimws(abbrev) == "Pisc Benthic Lg 4X"            ~ "Pisc Benth. G 4X",
                            trimws(abbrev) == "Pisc Benthic SM 4VW"           ~ "Pisc Benth. P 4VW",
                            trimws(abbrev) == "Pisc Benthic SM 4X"            ~ "Pisc Benth. P 4X",
                            trimws(abbrev) == "Pisc Pelagic SML 4VW"          ~ "Pisc Pélag. PET 4VW",
                            trimws(abbrev) == "Plank Pelagic SML 4VW"         ~ "Planc Pélag. PET 4VW",
                            trimws(abbrev) == "Plank Pelagic SML 4X"          ~ "Planc Pélag. PET 4X",
                            trimws(abbrev) == "Zoopisc Benthic SML 4VW"       ~ "Zoopisc Benth. PET 4VW",
                            trimws(abbrev) == "Zoopisc Benthic SML 4X"        ~ "Zoopisc Benth. PET 4X",
                            trimws(abbrev) == "Zoopisc Pelagic SML 4VW"       ~ "Zoopisc Pélag. PET 4VW",
                            trimws(abbrev) == "Zoopisc Pelagic SML 4X"        ~ "Zoopisc Pélag. PET 4X",
                            
                            # Seabird Foraging Guilds
                            trimws(abbrev) == "Plunge Dive Pisc"              ~ "Pisc plongeur",
                            trimws(abbrev) == "Pursuit Dive Pisc"             ~ "Pisc poursuite sous-marine",
                            trimws(abbrev) == "Pursuit Dive Plank"            ~ "Planc poursuite sous-marine",
                            trimws(abbrev) == "Shall Pursuit General"         ~ "Généraliste poursuite peu profonde",
                            trimws(abbrev) == "Ship-Follow General"           ~ "Généraliste suiveur de navire",
                            trimws(abbrev) == "Surface Seize Plank"           ~ "Planc saisie en surface",
                            trimws(abbrev) == "Surface Shall Dive Coast Pisc" ~ "Pisc côtier plongée peu prof./surface",
                            trimws(abbrev) == "Surface Shallo Dive Pisc"      ~ "Pisc plongée peu prof./surface",
                            
                            TRUE ~ as.character(abbrev) # Keeps unmatched strings intact
                          )
                        )

# Create blank rows for each type_f level
blank_rows_network <- network_targets_sum %>%
                      group_by(type_f) %>%
                      summarise(across(everything(), ~NA), .groups = "drop") 

# Bind and arrange
network_targets_padded1 <- network_targets_sum %>%
                          bind_rows(blank_rows_network) %>%
                          arrange(type_f, desc(Achieved))

#add blank row to the start to make this semetrical
network_targets_padded <- network_targets_padded1 %>%
                          slice(1) %>%
                          mutate(type_f = "Biodiversity hotspot") %>%
                          mutate(across(-feature, ~NA))%>%
                          rbind(.,network_targets_padded1)%>%
                          mutate(id = n():1,
                                 abbrev = ifelse(is.na(abbrev), "", abbrev),
                                 abbrev_fr = ifelse(is.na(abbrev_fr), "", abbrev_fr),
                                 target_met = ifelse(Achieved>=Minimum.Target,"yes","no"))


# Create named vector for labels
id_labels_network <- setNames(network_targets_padded$abbrev, network_targets_padded$id)
id_labels_network_fr <- setNames(network_targets_padded$abbrev_fr, network_targets_padded$id)

linedata <- network_targets_padded%>%
            dplyr::select(id,Achieved,type,target_met)%>%
            rbind(.,
                  network_targets_padded%>%
                    filter(target_met == "no")%>%
                    dplyr::select(id,Minimum.Target,type,target_met)%>%
                    rename(Achieved = Minimum.Target))

network_webca <- network_targets_padded%>%
                 dplyr::select(abbrev,abbrev_fr,type,id,target_met)%>%
                 left_join(webca_targets_padded%>%
                             filter(abbrev !="")%>%
                             dplyr::select(abbrev,Achieved))


#which features does WEBMR have more than 50% of in the network
webca_threshold_features <- webca_targets_padded%>%
                            filter(Percent.to.Goal>=50)%>%
                            pull(feature)


network_target_plot <- ggplot() +
                        geom_line(data = linedata,aes(y = factor(id),x=Achieved,group=factor(id)),lwd=0.5,col="grey70")+
                        geom_bar(data=network_targets_padded,aes(y = factor(id),  fill = type,x = Achieved),stat = "identity", col = NA,alpha=0.75) +
                        geom_bar(data=network_webca,aes(y = factor(id),fill = type, x = Achieved),stat="identity",col=NA)+
                        geom_bar(data=network_targets_padded%>%filter(feature %in% webca_threshold_features),aes(y = factor(id),x = Achieved),fill=NA,stat = "identity", col = "black") +
                        geom_point(data=network_targets_padded,aes(y = factor(id),  fill = type,x = Minimum.Target))+
                        theme_bw() +
                        theme(strip.background = element_blank(),
                              strip.text = element_blank(),
                              legend.position = "none",
                              axis.text.y = element_text(size = rel(0.7)),
                              panel.grid = element_blank()) +
                        scale_y_discrete(labels = id_labels_network) +
                        scale_x_continuous(expand = c(0.01, 0)) +
                        scale_fill_manual(values=c("darkmagenta","darkgoldenrod1","firebrick4",
                                                   "deepskyblue4","darkslateblue","cadetblue2","darkgreen",
                                                   "darkorange2","deeppink2","seagreen3"))+
                        labs(y = "", x = "% of design target")


network_target_plot_fr <- ggplot() +
                        geom_line(data = linedata,aes(y = factor(id),x=Achieved,group=factor(id)),lwd=0.5,col="grey70")+
                        geom_bar(data=network_targets_padded,aes(y = factor(id),  fill = type,x = Achieved),stat = "identity", col = NA,alpha=0.75) +
                        geom_bar(data=network_webca,aes(y = factor(id),fill = type, x = Achieved),stat="identity",col=NA)+
                        geom_bar(data=network_targets_padded%>%filter(feature %in% webca_threshold_features),aes(y = factor(id),x = Achieved),fill=NA,stat = "identity", col = "black") +
                        geom_point(data=network_targets_padded,aes(y = factor(id),  fill = type,x = Minimum.Target))+
                        theme_bw() +
                        theme(strip.background = element_blank(),
                              strip.text = element_blank(),
                              legend.position = "none",
                              axis.text.y = element_text(size = rel(0.7)),
                              panel.grid = element_blank()) +
                        scale_y_discrete(labels = id_labels_network_fr) +
                        scale_x_continuous(expand = c(0.01, 0)) +
                        scale_fill_manual(values=c("darkmagenta","darkgoldenrod1","firebrick4",
                                                   "deepskyblue4","darkslateblue","cadetblue2","darkgreen",
                                                   "darkorange2","deeppink2","seagreen3"))+
                        labs(y = "", x = "% de l'objectif de conception")

#save plots
ggsave("output/network_target_plot.png",network_target_plot,width=7,height=12,dpi = 1200)
ggsave("output/network_target_plot_fr.png",network_target_plot_fr,width=7,height=12,dpi = 1200)



#summaries

network_targets%>%
  distinct(feature,.keep_all=T)%>%
  pull(filter)%>%
  table()%>%
  data.frame()

network_targets%>%
  distinct(feature,.keep_all=T)%>%
  filter(feature %in% webca_threshold_features)%>%
  arrange(filter,type)%>%
  dplyr::select(filter,type,feature)
