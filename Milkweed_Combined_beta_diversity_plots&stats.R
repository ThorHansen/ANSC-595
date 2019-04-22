#Set up packages needed

library(reshape2)
library(phyloseq)
library(ggplot2)
library(RColorBrewer)
library(data.table)
library(Rmisc)
library(broom)
library(vegan)
library(gplots)
library(ggpubr)
library(dplyr)
library(plotrix)
library(tidyr)


source("C:/Users/warriorroad/Desktop/Insect_Microbial_Ecology_and_Evolution_Lab/Sequencing_data/code/DenefLab-MicrobeMiseq/R/miseqR.R") 
##From Denef lab microbiome analysis resources. http://deneflab.github.io/MicrobeMiseq/

set.seed(1)



##########################################################
#Normalizing ASVs by proportions and average library size#
##########################################################

phy_scale_mon <- phy.trim.mon %>%
  scale_reads(round = "round") 


phy_scale_rhizo <- phy.filter.rhizo %>%
  scale_reads(round = "round") 


#Save normalized phy objects in case you need them again
saveRDS(phy_scale_mon, "phy_scale_mon.rds")
saveRDS(phy_scale_rhizo, "phy_scale_rhizo.rds")


# If you close and re-open R
phy_scale_mon <- readRDS("phy_scale_mon.rds") #how to read RD files back in
phy_scale_rhizo <- readRDS("phy_scale_rhizo.rds") #how to read RD files back in











####################################
#Calculating beta diversity -insect#
####################################


#Switching colors between plants for plots
sample_data(phy_scale_mon)$PlantSpecies <- factor(sample_data(phy_scale_mon)$PlantSpecies, levels = c("Asyriaca","Acurassavica"))



###################################
#PCOA Beta diversity bray -insect

phy_scale_mon_pcoa_bray <- ordinate(
  physeq = phy_scale_mon, 
  method = "PCoA", 
  distance = "bray"
)

colourCount = 14
getPalette = colorRampPalette(brewer.pal(8, "Set1"))


phy_scale_mon_pcoa_bray_plot <- plot_ordination(
  physeq = phy_scale_mon,
  ordination = phy_scale_mon_pcoa_bray,
  color = "PlantSpecies",
  shape = "PlantSpecies",
  axes = 1:2,
  title = "Bray-Curtis dissimilarity PCoA of Monarch\nbacterial communities by Milkweed species") + 
  scale_color_manual(values = getPalette(colourCount)) +
  geom_point(aes(color = PlantSpecies), alpha = 0.7, size = 4) +
  geom_point(colour = "grey90", size = 1.5)


phy_scale_mon_pcoa_bray_plot




############################################################################
#Permanova and permdisp results of overall Bray-Curtis Dissimilarity -insect

# Calculate bray curtis distance matrix
all.bray.mon <- phyloseq::distance(phy_scale_mon, method = "bray")

# make a data frame from the sample_data
sampled.mon <- data.frame(sample_data(phy_scale_mon))

# Adonis test
adonis(all.bray.mon ~ PlantSpecies, data = sampled.mon)

# Homogeneity of dispersion test
beta.bray.mon <- betadisper(all.bray.mon, sampled.mon$PlantSpecies)
permutest(beta.bray.mon)




####################################
#PCOA Beta diversity jaccard -insect

phy_scale_mon.PA <- (otu_table(phy_scale_mon)>0)*1 #data>0 returns logical array with T or F
phy_scale_mon.PA <- phyloseq(otu_table(phy_scale_mon.PA, taxa_are_rows=TRUE), 
                             sample_data(phy_scale_mon), 
                             tax_table(tax_table(phy_scale_mon)))

phy_scale_mon_pcoa_jaccard.pa <- ordinate(
  physeq = phy_scale_mon.PA, 
  method = "PCoA", 
  trymax=100,
  #k=3,
  distance = "jaccard"
)

colourCount = 14
getPalette = colorRampPalette(brewer.pal(8, "Set1"))

phy_scale_mon_pcoa_jaccard_plot.pa <-plot_ordination(
  physeq = phy_scale_mon.PA,
  ordination = phy_scale_mon_pcoa_jaccard.pa,
  color = "PlantSpecies",
  shape = "PlantSpecies",
  axes = c(1,2),
  title = "Jaccard dissimilarity PCoA of Monarch\nbacterial communities by Milkweed species") + 
  scale_color_manual(values = getPalette(colourCount)) +
  geom_point(aes(color = PlantSpecies), alpha = 0.7, size = 4) +
  geom_point(colour = "grey90", size = 1.5) 

phy_scale_mon_pcoa_jaccard_plot.pa




########################################################################
#Permanova and permdisp results of overall Jaccard Dissimilarity -insect

# Calculate bray curtis distance matrix
all.jaccard.mon <- phyloseq::distance(phy_scale_mon.PA, method = "jaccard")

# make a data frame from the sample_data
sampled.mon <- data.frame(sample_data(phy_scale_mon.PA))

# Adonis test
adonis(all.jaccard.mon ~ PlantSpecies, data = sampled.mon)

# Homogeneity of dispersion test
beta.jaccard.mon <- betadisper(all.jaccard.mon, sampled.mon$PlantSpecies)
permutest(beta.jaccard.mon)




###############################################
#PCOA Beta diversity unweighted unifrac -insect

phy_scale_mon_unifrac <- ordinate(
  physeq = phy_scale_mon, 
  method = "PCoA", 
  distance = "unifrac"
)

colourCount = 14
getPalette = colorRampPalette(brewer.pal(8, "Set1"))

phy_scale_mon_pcoa_unifrac_plot <- plot_ordination(
  physeq = phy_scale_mon,
  ordination = phy_scale_mon_unifrac,
  color = "PlantSpecies",
  shape = "PlantSpecies",
  axes = 1:2,
  title = "Unweighted Unifrac PCoA of Monarch\nbacterial communities by Milkweed species") + 
  scale_color_manual(values = getPalette(colourCount)) +
  geom_point(aes(color = PlantSpecies), alpha = 0.7, size = 4) +
  geom_point(colour = "grey90", size = 1.5) 

phy_scale_mon_pcoa_unifrac_plot




#####################################################################
#Permanova and permdisp results of overall unweighted unifrac -insect

# Calculate bray curtis distance matrix
all.unifrac.mon <- phyloseq::distance(phy_scale_mon, method = "unifrac")

# make a data frame from the sample_data
sampled.mon <- data.frame(sample_data(phy_scale_mon))

# Adonis test
adonis(all.unifrac.mon ~ PlantSpecies, data = sampled.mon)

# Homogeneity of dispersion test
beta.unifrac.mon <- betadisper(all.unifrac.mon, sampled.mon$PlantSpecies)
permutest(beta.unifrac.mon)




#############################################
#PCOA Beta diversity weighted unifrac -insect

phy_mon_scale_wunifrac <- ordinate(
  physeq = phy_scale_mon, 
  method = "PCoA", 
  distance = "wunifrac"
)

colourCount = 14
getPalette = colorRampPalette(brewer.pal(8, "Set1"))

phy_scale_mon_pcoa_wunifrac_plot <-plot_ordination(
  physeq = phy_scale_mon,
  ordination = phy_mon_scale_wunifrac,
  color = "PlantSpecies",
  shape = "PlantSpecies",
  axes = 1:2,
  title = "Weighted Unifrac PCoA of Monarch\nbacterial communities by Milkweed species") + 
  scale_color_manual(values = getPalette(colourCount)) +
  geom_point(aes(color = PlantSpecies), alpha = 0.7, size = 4) +
  geom_point(colour = "grey90", size = 1.5) 

phy_scale_mon_pcoa_wunifrac_plot




###################################################################
#Permanova and permdisp results of overall weighted unifrac -insect

all.wunifrac.mon <- phyloseq::distance(phy_scale_mon, method = "wunifrac")

# make a data frame from the sample_data
sampled.mon <- data.frame(sample_data(phy_scale_mon))

# Adonis test
adonis(all.wunifrac.mon ~ PlantSpecies, data = sampled.mon)

# Homogeneity of dispersion test
beta.wunirac.mon <- betadisper(all.wunifrac.mon, sampled.mon$PlantSpecies)
permutest(beta.wunifrac.mon)










#########################################
#Calculating beta diversity -rhizosphere#
#########################################


#need to switch colors between plants
sample_data(phy_scale_rhizo)$PlantSpecies <- factor(sample_data(phy_scale_rhizo)$PlantSpecies, levels = c("Asyriaca","Acurassavica"))



################################
#PCOA Beta diversity bray -rhizo

phy_scale_rhizo_bray_pcoa <- ordinate(
  physeq = phy_scale_rhizo, 
  method = "PCoA", 
  distance = "bray"
)

colourCount = 14
getPalette = colorRampPalette(brewer.pal(8, "Set1"))

phy_scale_rhizo_pcoa_bray_plot <- plot_ordination(
  physeq = phy_scale_rhizo,
  ordination = phy_scale_rhizo_bray_pcoa,
  color = "PlantSpecies",
  shape = "InsectPresence",
  axes = 1:2,
  title = "Bray-Curtis dissimilarity PCoA of Milkweed\nrhizosphere bacterial communities by\nMilkweed species and Monarch presence") + 
  scale_color_manual(values = getPalette(colourCount)) +
  geom_point(aes(color = PlantSpecies), alpha = 0.7, size = 4) +
  geom_point(colour = "grey90", size = 1.5) 

phy_scale_rhizo_pcoa_bray_plot




###########################################################################
#Permanova and permdisp results of overall Bray-Curtis Dissimilarity -rhizo

# Calculate bray curtis distance matrix
all.bray.rhizo <- phyloseq::distance(phy_scale_rhizo, method = "bray")

# make a data frame from the sample_data
sampled.rhizo <- data.frame(sample_data(phy_scale_rhizo))

# Adonis test
adonis(all.bray.rhizo ~ PlantSpecies*InsectPresence, data = sampled.rhizo)

# Homogeneity of dispersion test
beta.bray.plant.rhizo <- betadisper(all.bray.rhizo, sampled.rhizo$PlantSpecies)
permutest(beta.bray.plant.rhizo)

beta.bray.insect.rhizo <- betadisper(all.bray.rhizo, sampled.rhizo$InsectPresence)
permutest(beta.bray.insect.rhizo)




###################################
#PCOA Beta diversity jaccard -rhizo

phy_scale_rhizo.PA <- (otu_table(phy_scale_rhizo)>0)*1 #data>0 returns logical array with T or F
phy_scale_rhizo.PA <- phyloseq(otu_table(phy_scale_rhizo.PA, taxa_are_rows=TRUE), 
                               sample_data(phy_scale_rhizo), 
                               tax_table(tax_table(phy_scale_rhizo)))

phy_scale_rhizo_pcoa_jaccard.pa <- ordinate(
  physeq = phy_scale_rhizo.PA, 
  method = "PCoA", 
  trymax=100,
  #k=3,
  distance = "jaccard"
)

phy_scale_rhizo_pcoa_jaccard_plot.pa <-plot_ordination(
  physeq = phy_scale_rhizo.PA,
  ordination = phy_scale_rhizo_pcoa_jaccard.pa,
  color = "PlantSpecies",
  shape = "InsectPresence",
  axes = c(1,2),
  title = "Jaccard dissimilarity PCoA of Milkweed\nrhizosphere bacterial communities by\nMilkweed species and Monarch presence") + 
  scale_color_manual(values = getPalette(colourCount)) +
  geom_point(aes(color = PlantSpecies), alpha = 0.7, size = 4) +
  geom_point(colour = "grey90", size = 1.5) 

phy_scale_rhizo_pcoa_jaccard_plot.pa




#######################################################################
#Permanova and permdisp results of overall Jaccard Dissimilarity -rhizo

# Calculate jaccard distance matrix
all.jaccard.rhizo <- phyloseq::distance(phy_scale_rhizo.PA, method = "jaccard")

# make a data frame from the sample_data
sampled.rhizo <- data.frame(sample_data(phy_scale_rhizo.PA))

# Adonis test
adonis(all.jaccard.rhizo ~ PlantSpecies*InsectPresence, data = sampled.rhizo)

# Homogeneity of dispersion test
beta.jaccard.plant.rhizo <- betadisper(all.jaccard.rhizo, sampled.rhizo$PlantSpecies)
permutest(beta.jaccard.plant.rhizo)

beta.jaccard.insect.rhizo <- betadisper(all.jaccard.rhizo, sampled.rhizo$InsectPresence)
permutest(beta.jaccard.insect.rhizo)




##############################################
#PCOA Beta diversity unweighted unifrac -rhizo

phy_scale_rhizo_unifrac <- ordinate(
  physeq = phy_scale_rhizo, 
  method = "PCoA", 
  distance = "unifrac"
)

colourCount = 14
getPalette = colorRampPalette(brewer.pal(8, "Set1"))

phy_scale_rhizo_pcoa_unifrac_plot <- plot_ordination(
  physeq = phy_scale_rhizo,
  ordination = phy_scale_rhizo_unifrac,
  color = "PlantSpecies",
  shape = "InsectPresence",
  axes = 1:2,
  title = "Unweighted Unifrac PCoA of Milkweed\nrhizosphere bacterial communities by\nMilkweed species and Monarch presence") + 
  scale_color_manual(values = getPalette(colourCount)) +
  geom_point(aes(color = PlantSpecies), alpha = 0.7, size = 4) +
  geom_point(colour = "grey90", size = 1.5) 

phy_scale_rhizo_pcoa_unifrac_plot




####################################################################
#Permanova and permdisp results of overall unweighted unifrac -rhizo

# Calculate unweighted unifrac distance matrix
all.unifrac.rhizo <- phyloseq::distance(phy_scale_rhizo, method = "unifrac")

# make a data frame from the sample_data
sampled.rhizo <- data.frame(sample_data(phy_scale_rhizo))

# Adonis test
adonis(all.unifrac.rhizo ~ PlantSpecies*InsectPresence, data = sampled.rhizo)

# Homogeneity of dispersion test
beta.unifrac.plant.rhizo <- betadisper(all.unifrac.rhizo, sampled.rhizo$PlantSpecies)
permutest(beta.unifrac.plant.rhizo)

beta.unifrac.insect.rhizo <- betadisper(all.unifrac.rhizo, sampled.rhizo$InsectPresence)
permutest(beta.unifrac.insect.rhizo)




############################################
#PCOA Beta diversity Weighted unifrac -rhizo

phy_rhizo_scale_wunifrac <- ordinate(
  physeq = phy_scale_rhizo, 
  method = "PCoA", 
  distance = "wunifrac"
)

colourCount = 14
getPalette = colorRampPalette(brewer.pal(8, "Set1"))

phy_scale_rhizo_pcoa_wunifrac_plot <-plot_ordination(
  physeq = phy_scale_rhizo,
  ordination = phy_rhizo_scale_wunifrac,
  color = "PlantSpecies",
  shape = "InsectPresence",
  axes = 1:2,
  title = "Weighted Unifrac PCoA of Milkweed\nrhizosphere bacterial communities by\nMilkweed species and Monarch presence") + 
  scale_color_manual(values = getPalette(colourCount)) +
  geom_point(aes(color = PlantSpecies), alpha = 0.7, size = 4) +
  geom_point(colour = "grey90", size = 1.5) 

phy_scale_rhizo_pcoa_wunifrac_plot




#####################################################
#Permanova results of overall Weighted unifrac -rhizo

# Calculate weighted unifrac distance matrix
all.wunifrac.rhizo <- phyloseq::distance(phy_scale_rhizo, method = "wunifrac")

# make a data frame from the sample_data
sampled.rhizo <- data.frame(sample_data(phy_scale_rhizo))

# Adonis test
adonis(all.wunifrac.rhizo ~ PlantSpecies*InsectPresence, data = sampled.rhizo)

# Homogeneity of dispersion test
beta.wunifrac.plant.rhizo <- betadisper(all.wunifrac.rhizo, sampled.rhizo$PlantSpecies)
permutest(beta.wunifrac.plant.rhizo)

beta.wunifrac.insect.rhizo <- betadisper(all.wunifrac.rhizo, sampled.rhizo$InsectPresence)
permutest(beta.wunifrac.insect.rhizo)


