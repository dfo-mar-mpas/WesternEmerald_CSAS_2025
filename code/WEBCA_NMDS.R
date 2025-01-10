#WEBCA data

library(vegan)
library(dplyr)

edna.mat <- read.csv("data/2020eDNA_COI12S16S_FinalTaxaTable_WEBCAonly.csv", header = T) %>% glimpse()
edna.mat2 <- t(edna.mat[,-1])
colnames(edna.mat2)<-edna.mat$X.OTU.ID
edna.mat3 <- edna.mat2[rowSums(edna.mat2[])>0,]

web.nmds <- metaMDS(edna.mat3,distance = "jaccard",k=3, trymax = 100)


#extract nmds scores for ggplot
data.scores = as.data.frame(scores(web.nmds)$sites)

data.scores$Sample <- rownames(data.scores)
data.scores$Inside <- c(rep("inside",3),rep("outside",4))


species.scores <- as.data.frame(scores(web.nmds, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
species.scores$species <- rownames(species.scores) 

hull.webca <- data.scores %>%
  as.data.frame() %>%
  group_by(Inside) %>%
  slice(chull(x=NMDS1,y=NMDS2)) #for this dataset, there's no differentiation between bottom and surface

p2 = ggplot() + 
  geom_polygon(data=hull.webca,aes(x=NMDS1,y=NMDS2,fill=Inside),color="black",alpha=0.30) +
  scale_fill_manual(values=c("firebrick","cornflowerblue"))+
  #ggnewscale::new_scale_fill()+
  geom_point(data=data.scores, aes(x = NMDS1, y = NMDS2, fill=Inside),size = 4, shape=21, colour="black")+
  #scale_fill_manual(values=colorScales)+
  #geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species), alpha=0.5)+
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 14), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2)) + 
  labs(x = "NMDS1", y = "NMDS2")  + 
  geom_text(aes(x=Inf, y=Inf, vjust=48,hjust=1.2,label=paste("Stress =",round(web.nmds$stress,3),"k =",web.nmds$ndim)));p2


# Make a list of species present per sample number
edna.list<-edna.mat %>%
  pivot_longer(cols=-X.OTU.ID, names_to="Sample",values_to = "count")

specs.list <- edna.list %>% 
  filter(count>0) %>% 
  rename(Species=X.OTU.ID) %>%
  mutate(Inside = ifelse(Sample %in% c("EDNA20_478631","EDNA20_478635","EDNA20_478639"),"inside","outside")) %>%
  group_by(Sample)

save.image("data/WEBCASpeciesList.RData")
