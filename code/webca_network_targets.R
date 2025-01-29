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
       x="% of design target")

ggsave("output/webca_design_targets.png",webca_plot,height=8,width=4,units="in",dpi=300)
  



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
                                   dplyr::select(feature,type,filter,Minimum.Target,abbrev,type_f))

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
                                 target_met = ifelse(Achieved>=Minimum.Target,"yes","no"))


# Create named vector for labels
id_labels_network <- setNames(network_targets_padded$abbrev, network_targets_padded$id)

linedata <- data.frame(abbrev_lab=rep(plotdata$abbrev_lab,2),plan=c(plotdata$plan,plotdata$target))

linedata <- network_targets_padded%>%
            dplyr::select(id,Achieved,type,target_met)%>%
            rbind(.,
                  network_targets_padded%>%
                    filter(target_met == "no")%>%
                    dplyr::select(id,Minimum.Target,type,target_met)%>%
                    rename(Achieved = Minimum.Target))

network_webca <- network_targets_padded%>%
                 dplyr::select(abbrev,type,id,target_met)%>%
                 left_join(webca_targets_padded%>%
                             filter(abbrev !="")%>%
                             dplyr::select(abbrev,Achieved))
              
            


network_target_plot <- ggplot() +
                        geom_line(data = linedata,aes(y = factor(id),x=Achieved,group=factor(id)),lwd=0.5,col="grey70")+
                        geom_bar(data=network_targets_padded,aes(y = factor(id),  fill = type,x = Achieved,alpha = target_met),stat = "identity", col = NA) +
                        geom_bar(data=network_webca,aes(y = factor(id),fill = type, x = Achieved,alpha = target_met),stat="identity",col="black")+
                        geom_point(data=network_targets_padded,aes(y = factor(id),  fill = type,x = Minimum.Target))+
                        scale_alpha_manual(values = c("yes" = 1, "no" = 0.5)) +  # Assign transparency manually
                        theme_bw() +
                        theme(strip.background = element_blank(),
                              strip.text = element_blank(),
                              legend.position = "none",
                              axis.text.y = element_text(size = rel(0.7)),
                              panel.grid = element_blank()) +
                        scale_y_discrete(labels = id_labels_network) +
                        scale_x_continuous(expand = c(0.01, 0)) +
                        labs(y = "", x = "% of design target")

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

ggsave("output/network_target_plot.png",network_target_plot,width=7,height=12,dpi = 1200)
