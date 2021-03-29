## FUNCTIONS ##

#Function for abundance
abundance = function(x){return(sum(x,na.rm=TRUE))}

#Grouping the ranks, conversion of ASV data/table into abundance table for various taxonomic ranks
group_by_rank = function(otu_table,taxonomy,metadata,rank){
  n_sample = nrow(otu_table)
  otu_table = cbind(t(otu_table),(taxonomy))
  if (rank>1){
    otu_table$taxa = apply(otu_table[,(n_sample+1):(n_sample+rank)],1,paste,collapse="")
  }else{
    otu_table$taxa = otu_table[,n_sample+1]
  }
  otu_table = group_by(otu_table,taxa)
  grouped_table = summarise_at(otu_table,1:n_sample,sum)
  taxa_names = grouped_table$taxa
  grouped_table = grouped_table[,-1]
  grouped_table = t(grouped_table)
  colnames(grouped_table) = taxa_names
  taxa_cols = 1:ncol(grouped_table)
  grouped_table = cbind(grouped_table,metadata)
  taxa_rel_abundance = colSums(grouped_table[,taxa_cols]/rowSums(grouped_table[,taxa_cols],na.rm=T),na.rm=T)
  grouped_table = group_by(grouped_table,SampleID)
  grouped_table = pivot_longer(grouped_table,names_to="taxa",values_to="abundance",cols=all_of(taxa_cols))
  grouped_table$taxa = factor(grouped_table$taxa,levels=names(sort(taxa_rel_abundance,decreasing=TRUE)))
  return(grouped_table)}

## WORKSPACE DATA ##
#Installing Packages
install.packages("zoo")
install.packages("tidyverse")
install.packages("vegan")
install.packages("viridis")
install.packages("data.table")
install.packages("dplyr")
install.packages("tidyr")
install.packages("stringr")
install.packages("ggplot2")

#Loading Packages
library(zoo)
library(tidyverse)
library(vegan)
library(viridis)
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2);theme_set(theme_bw()+
                             theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                                   text=element_text(size = 15, face= "bold"),
                                   panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                   legend.text = element_text(size=12))) 


#Set Working Directory
setwd("~/Desktop/BIOL 403/Project Data")

#Loading Data

#ASV data
#reading in ASV data
seagrass_ASV=read.table("seagrass_diatom_removal_asv.txt",header=TRUE,sep="\t")
#filtering the data to exclude unecessary data, we are only focusing on control and GeO2 samples (Sample IDs below are those that are NOT)
#Control or GeO2
seagrass_ASV.filtered<-subset(seagrass_ASV, SampleID=="T21A"|SampleID=="T22A"|SampleID=="T23A"|SampleID=="T27A"
                              |SampleID=="T28A"|SampleID=="T29A"|SampleID=="T21B"|SampleID=="T22B"|SampleID=="T23B"
                              |SampleID=="T24B"|SampleID=="T25B"|SampleID=="T26B"|SampleID=="T27B"|SampleID=="T28B"
                              |SampleID=="T29B"|SampleID=="T210B"|SampleID=="T311A"|SampleID=="T312A"|SampleID=="T313A"
                              |SampleID=="T314A"|SampleID=="T316A"|SampleID=="T318A"|SampleID=="T311B"|SampleID=="T312B"
                              |SampleID=="T313B"|SampleID=="T314B"|SampleID=="T315B"|SampleID=="T316B"|SampleID=="T317B"
                              |SampleID=="T318B"|SampleID=="T319B"|SampleID=="T320B")

#taxonomy data
#reading in taxonomy data, replacing blanks with NA
taxonomy_original = read.table("seagrass_diatom_removal_taxonomy.txt",header=TRUE,sep="\t",na.strings=c("","NA"))
#seaparation code so taxonomy does not appear in one cell (i.e.phylumfamilygenus)
taxonomy = taxonomy_original %>%
  t() %>% 
  na.locf() %>% 
  t() 
write.csv(taxonomy, "taxonomy_propagated.csv")
taxonomy = read.csv("taxonomy_propagated.csv")
#removing ASV and accession, these are not relevant to the taxonomy data frame 
seagrass_taxonomy=taxonomy[3:9]

#metadata
#reading in metadata, replacing blanks with NA
seagrass_metadata=read.table("seagrass_diatom_removal_metadata.txt",header=TRUE,sep="\t",na.strings=c("","NA"))
#changing each sample (i.e. Control 1) to the type of treatment rather than treatment and sample number
#(treatment sample 1 -> treatment)
seagrass_metadata[seagrass_metadata=="Control 1"]<-"Control"
seagrass_metadata[seagrass_metadata=="Control 2"]<-"Control"
seagrass_metadata[seagrass_metadata=="Control 3"]<-"Control"
seagrass_metadata[seagrass_metadata=="Control 4"]<-"Control"
seagrass_metadata[seagrass_metadata=="Control 5"]<-"Control"
seagrass_metadata[seagrass_metadata=="GeO2 1"]<-"GeO2"
seagrass_metadata[seagrass_metadata=="GeO2 2"]<-"GeO2"
seagrass_metadata[seagrass_metadata=="GeO2 3"]<-"GeO2"
seagrass_metadata[seagrass_metadata=="GeO2 4"]<-"GeO2"
seagrass_metadata[seagrass_metadata=="GeO2 5"]<-"GeO2"
seagrass_metadata[seagrass_metadata=="Recol 1"]<-"Recolonization"
seagrass_metadata[seagrass_metadata=="Recol 2"]<-"Recolonization"
seagrass_metadata[seagrass_metadata=="Recol 3"]<-"Recolonization"
seagrass_metadata[seagrass_metadata=="Recol 4"]<-"Recolonization"
seagrass_metadata[seagrass_metadata=="Recol 5"]<-"Recolonization"
seagrass_metadata[seagrass_metadata=="Antibiotics 1"]<-"Antibiotic"
seagrass_metadata[seagrass_metadata=="Antibiotics 2"]<-"Antibiotic"
seagrass_metadata[seagrass_metadata=="Antibiotics 3"]<-"Antibiotic"
seagrass_metadata[seagrass_metadata=="Antibiotics 4"]<-"Antibiotic"
seagrass_metadata[seagrass_metadata=="Antibiotics 5 "]<-"Antibiotic"
seagrass_metadata[seagrass_metadata=="A + GeO2 sterile 1"]<-"A + GeO2 Sterile"
seagrass_metadata[seagrass_metadata=="A + GeO2 sterile 2"]<-"A + GeO2 Sterile"
seagrass_metadata[seagrass_metadata=="A + GeO2 sterile 3"]<-"A + GeO2 Sterile"
seagrass_metadata[seagrass_metadata=="A + GeO2 sterile 4"]<-"A + GeO2 Sterile"
seagrass_metadata[seagrass_metadata=="A + GeO2 sterile 5"]<-"A + GeO2 Sterile"
#filtering data to contain only the treatments we want to look at, Control and GeO2
seagrass_metadata.filtered<-subset(seagrass_metadata, Treatment=="Control"|Treatment=="GeO2")

## DATA PROCESSING ##

#Figure 1#
#Merge ASV and Metadata Filtered (specifically look at Control and GeO2) based on SampleID
metaasv = inner_join(seagrass_metadata.filtered,seagrass_ASV.filtered[,-(2:3)],by="SampleID")
#set seed for reproducibility
set.seed(4)
#NMDS performed on ASV Filtered (removing the Sample ID column, not needed)
seagrass_ASV.filtered<-seagrass_ASV.filtered[,2:2875] 
#NMDS reduced to 2 dimensions, measured with Bray-Curtis
nmds_s4=metaMDS(seagrass_ASV.filtered,k=2,distance="bray") 
nmds_s4
#Initial NMDS plot
plot(nmds_s4)
#The X and Y coordinates/NMDS scores extracted
data.scores<-as.data.frame(scores(nmds_s4))
#Adding original data information to data.scores, changing the RHS to specifically look at the treatment
#and trial
data.scores$Treatment=as.character(metaasv$Treatment)
data.scores$Trial=as.character(metaasv$Trial)
#Confirmation that columns have been added
head(data.scores)

#Figure 2#


#Figure 3#
#Joining metadata and ASV data by SampleID
asvmeta = inner_join(x=seagrass_metadata,y=seagrass_ASV,
                     by="SampleID")

#removing duplicate metadata columns 
asv_cols=c(7:2881)

#calculating the ASV total abundance of each sample/SampleID
asvmeta$total_abundance=apply(asvmeta[,asv_cols],1,abundance)
#Creation of separate row for each ASV per sample/SampleID
asvmeta_long = pivot_longer(asvmeta,names_to="OTU",values_to="abundance",cols=all_of(asv_cols))
#Calculating and creating new column for relative abundance
asvmeta_long$relative_abundance = asvmeta_long$abundance/asvmeta_long$total_abundance
#Calculation for total relative abundance of each ASV
asv_rel_abundance = colSums(asvmeta[,asv_cols]/asvmeta$total_abundance)
#reordering ASVs, most abundant ASVs appear as the top bar
asvmeta_long$OTU = factor(asvmeta_long$OTU,levels=names(sort(asv_rel_abundance,decreasing=TRUE)))

#Preparing for separating the taxonomy original by inserting underscores 
seagrass_taxonomy[]<-lapply(seagrass_taxonomy,function(x)paste("_",x,sep=""))

#adding taxonomy data to metadata, sorts based on various ranks
asvmeta_order = group_by_rank(otu_table = asvmeta[,asv_cols],
                              taxonomy  = seagrass_taxonomy,
                              metadata  = asvmeta[,-asv_cols],
                              rank      = 6) 

#Creating a new column for and calculating the realtive abundances for asvmeta_order data frame
asvmeta_order$relative_abundance = asvmeta_order$abundance/asvmeta_order$total_abundance
#filtering the data to include the top most 15 abundant orders
top_orders = levels(asvmeta_order$taxa)[1:15]
asvmeta_top_orders = filter(asvmeta_order,taxa %in% top_orders)

#displaying all "Other" taxa appropriately and create new column for information 
taxa_levels = c(levels(asvmeta_top_orders$taxa),"Other")
for (s in 1:nrow(asvmeta)){
  sample = asvmeta$SampleID[s]
  sample_only = filter(asvmeta_top_orders,SampleID==sample)
  other_abundance = sample_only$total_abundance[1] - sum(sample_only$abundance)
  other_row = sample_only[1,]
  other_row$taxa = "Other"
  other_row$abundance = other_abundance
  other_row$relative_abundance = other_abundance/other_row$total_abundance
  asvmeta_top_orders = rbind(asvmeta_top_orders,other_row)
}
#separating/sorting the taxonomy to appropriate ranks from asv_meta_top_orders data frame (top most 15 abundant orders)
top_split_taxa = separate(data = asvmeta_top_orders, 
                          col = taxa, 
                          into = c("Rank_blank","Rank_ASV","Rank1","Rank2","Rank3", "Rank4","Rank5","Rank6"), 
                          sep = "_") 

#Creating subsets of top_split_taxa based on treatment type: Control and GeO2 (treatment of interest)
Control_top_taxa = subset(top_split_taxa, top_split_taxa$Treatment=="Control")
GeO2_top_taxa = subset(top_split_taxa, top_split_taxa$Treatment=="GeO2")

## STATISTICS/GRAPH ##

#Figure 1#
#Plotting and formatting NMDS plot
ggplot(data.scores, aes(x=NMDS1, y=NMDS2))+
  geom_point(size=3, aes(shape=Treatment, color=Trial))+
  theme_classic()+
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"),
  axis.text.x = element_text(colour = "black", face = "bold", size = 12),
  legend.text = element_text(size = 12, face ="bold", colour ="black"),
  legend.position = "right", axis.title.y = element_text(face = "bold", size = 14),
  axis.title.x = element_text(face = "bold", size = 14, colour = "black"),
  legend.title = element_text(size = 14, colour = "black", face = "bold"),
  legend.key=element_blank())+
  scale_colour_manual(values=c("#000000","#CCCCCC"))+
  scale_shape_manual(values=c(16,17))+
  ggtitle("NMDS Plot")+theme(plot.title=element_text(face="bold",size=18,hjust=0.5))+
  coord_cartesian(xlim=c(-2,2),ylim=c(-2,2))

#Statistic 2#
#Performing Krusal. test
kruskal.test(richness~treatment, data=treatment_and_richness)
#p-value is greater than 0.05 (0.3173)
#Fail to reject the null hypothesis of equal species richness among the two treatment groups

#Figure 3#
#Colour Palette
my_colors = c("blue","steelblue4","steelblue2","skyblue",
              "turquoise","springgreen","springgreen4","olivedrab",
              "yellowgreen","yellow2","gold3",
              "darkgoldenrod3","darkgoldenrod4","mistyrose4",
              "darkgrey","grey40")
#Plotting the taxa plot for GeO2 treatment (Rank 5)
ggplot(GeO2_top_taxa,aes(x=SampleID,y=relative_abundance, fill=Rank5)) +
  geom_bar(stat="identity",position="stack")+
  scale_fill_manual(values = my_colors, na.value="black")+
  xlab("Sample ID")+
  ylab("Relative Abundance")+
  labs(fill="Genus")+
  facet_grid(.~Treatment,scales = "free_x", space="free")+
  ggtitle("Taxa Plot - GeO2")+theme(plot.title=element_text(face="bold",size=20,hjust=0.5))

#Plotting the taxa plot for Control treatment (Rank 5)
ggplot(Control_top_taxa,aes(x=SampleID,y=relative_abundance, fill=Rank5)) +
  geom_bar(stat="identity",position="stack")+
  scale_fill_manual(values = my_colors, na.value="black")+
  xlab("Sample ID")+
  ylab("Relative Abundance")+
  labs(fill="Genus")+
  facet_grid(.~Treatment,scales = "free_x", space="free")+
  ggtitle("Taxa Plot - Control")+theme(plot.title=element_text(face="bold",size=20,hjust=0.5))
