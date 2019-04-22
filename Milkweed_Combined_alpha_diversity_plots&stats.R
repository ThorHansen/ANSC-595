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











#############################
#Calculating alpha diversity#
#############################


###################################
#Calculate alpha diversity - insect

min_lib_mon <- min(sample_sums(phy_scale_mon)) #filtered data

# Initialize matrices to store richness and evenness estimates
nsamp_mon = nsamples(phy_scale_mon)
trials = 100

richness_mon <- matrix(nrow = nsamp_mon, ncol = trials)
row.names(richness_mon) <- sample_names(phy_scale_mon)

evenness_mon <- matrix(nrow = nsamp_mon, ncol = trials)
row.names(evenness_mon) <- sample_names(phy_scale_mon)

shannon_mon <- matrix(nrow = nsamp_mon, ncol = trials)
row.names(shannon_mon) <- sample_names(phy_scale_mon)

simpson_mon <- matrix(nrow = nsamp_mon, ncol = trials)
row.names(simpson_mon) <- sample_names(phy_scale_mon)





for (i in 1:100) {
  # Subsample
  r_mon <- rarefy_even_depth(phy_scale_mon, sample.size = min_lib_mon, verbose = FALSE, replace = T)
  
  # Calculate richness
  rich_mon <- as.numeric(as.matrix(estimate_richness(r_mon, measures = "Observed")))
  richness_mon[ ,i] <- rich_mon
  
  # Calculate evenness
  even_mon <- as.numeric(as.matrix(estimate_richness(r_mon, measures = "InvSimpson")))
  evenness_mon[ ,i] <- even_mon
  
  # calculate shannon diversity
  shan_mon <- as.numeric(as.matrix(estimate_richness(r_mon, measures = "Shannon")))
  shannon_mon[ ,i] <- shan_mon
  
  # calculate simpson diversity
  simp_mon <- as.numeric(as.matrix(estimate_richness(r_mon, measures = "Simpson")))
  simpson_mon[ ,i] <- simp_mon
  
}

# Create a new dataframe to hold the means and standard deviations of richness estimates
SampleID <- row.names(richness_mon)
mean <- apply(richness_mon, 1, mean)
sd <- apply(richness_mon, 1, sd)
measure <- rep("Richness", nsamp_mon)
rich_stats_mon <- data.frame(SampleID, mean, sd, measure)

# Create a new dataframe to hold the means and standard deviations of evenness estimates
SampleID <- row.names(evenness_mon)
mean <- apply(evenness_mon, 1, mean)
sd <- apply(evenness_mon, 1, sd)
measure <- rep("Inverse\nSimpson", nsamp_mon)
even_stats_mon <- data.frame(SampleID, mean, sd, measure)

# Create a new dataframe to hold the means and standard deviations of shannon diversity estimates
SampleID <- row.names(shannon_mon)
mean <- apply(shannon_mon, 1, mean)
sd <- apply(shannon_mon, 1, sd)
measure <- rep("Shannon\nDiversity", nsamp_mon)
shan_stats_mon <- data.frame(SampleID, mean, sd, measure)

# Create a new dataframe to hold the means and standard deviations of simpson diversity estimates
SampleID <- row.names(simpson_mon)
mean <- apply(simpson_mon, 1, mean)
sd <- apply(simpson_mon, 1, sd)
measure <- rep("Simpson\nDiversity", nsamp_mon)
simp_stats_mon <- data.frame(SampleID, mean, sd, measure)


#combine richness and evenness into one dataframe
alpha_mon <- rbind(rich_stats_mon, even_stats_mon, shan_stats_mon, simp_stats_mon)
s_mon <- data.frame(sample_data(phy_scale_mon))
s_mon$SampleID <- SampleID
alphadivmon <- merge(alpha_mon, s_mon, by = "SampleID")
alphadivmon$PlantSpecies <- factor(alphadivmon$PlantSpecies, levels = c("Asyriaca","Acurassavica"))




###########################################
#Visualize alpha div with barplots - insect

colourCount = 14
getPalette = colorRampPalette(brewer.pal(8, "Set1"))


alphadivmon.plant.plot <- ggplot(alphadivmon, aes(x = PlantSpecies, y = mean, color = PlantSpecies, group = PlantSpecies)) +
  geom_boxplot() +
  facet_grid(measure ~ PlantSpecies, scales = "free", space = "free_x") +
  theme(strip.text.y = element_text(size = 8)) +
  geom_point(aes(fill = PlantSpecies), size = 1, shape = 21) +
  scale_color_manual(values = getPalette(colourCount)) +
  scale_fill_manual(values = getPalette(colourCount)) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_blank())

alphadivmonplot.plant.plot




##########################################
#Alpha div signif by plant species -insect

erich <- estimate_richness(phy_scale_mon, measures = c("Observed", "InvSimpson", "Shannon", "Simpson"))
ttest <- t(sapply(erich, function(x) unlist(t.test(x~sample_data(phy_scale_mon)$PlantSpecies )[c("estimate","p.value","statistic","conf.int")])))
ttest










#########################################
#Calculate alpha diversity - rhizo

min_lib_rhizo <- min(sample_sums(phy_scale_rhizo)) #filtered data

# Initialize matrices to store richness and evenness estimates
nsamp_rhizo = nsamples(phy_scale_rhizo)
trials = 100

richness_rhizo <- matrix(nrow = nsamp_rhizo, ncol = trials)
row.names(richness_rhizo) <- sample_names(phy_scale_rhizo)

evenness_rhizo <- matrix(nrow = nsamp_rhizo, ncol = trials)
row.names(evenness_rhizo) <- sample_names(phy_scale_rhizo)

shannon_rhizo <- matrix(nrow = nsamp_rhizo, ncol = trials)
row.names(shannon_rhizo) <- sample_names(phy_scale_rhizo)

simpson_rhizo <- matrix(nrow = nsamp_rhizo, ncol = trials)
row.names(simpson_rhizo) <- sample_names(phy_scale_rhizo)





for (r in 1:100) {
  # Subsample
  r_rhizo <- rarefy_even_depth(phy_scale_rhizo, sample.size = min_lib_rhizo, verbose = FALSE, replace = T)
  
  # Calculate richness
  rich_rhizo <- as.numeric(as.matrix(estimate_richness(r_rhizo, measures = "Observed")))
  richness_rhizo[ ,r] <- rich_rhizo
  
  # Calculate evenness
  even_rhizo <- as.numeric(as.matrix(estimate_richness(r_rhizo, measures = "InvSimpson")))
  evenness_rhizo[ ,r] <- even_rhizo
  
  # calculate shannon diversity
  shan_rhizo <- as.numeric(as.matrix(estimate_richness(r_rhizo, measures = "Shannon")))
  shannon_rhizo[ ,r] <- shan_rhizo
  
  # calculate simpson diversity
  simp_rhizo <- as.numeric(as.matrix(estimate_richness(r_rhizo, measures = "Simpson")))
  simpson_rhizo[ ,r] <- simp_rhizo
  
}

# Create a new dataframe to hold the means and standard deviations of richness estimates
SampleID <- row.names(richness_rhizo)
mean <- apply(richness_rhizo, 1, mean)
sd <- apply(richness_rhizo, 1, sd)
measure <- rep("Richness", nsamp_rhizo)
rich_stats_rhizo <- data.frame(SampleID, mean, sd, measure)

# Create a new dataframe to hold the means and standard deviations of evenness estimates
SampleID <- row.names(evenness_rhizo)
mean <- apply(evenness_rhizo, 1, mean)
sd <- apply(evenness_rhizo, 1, sd)
measure <- rep("Inverse\nSimpson", nsamp_rhizo)
even_stats_rhizo <- data.frame(SampleID, mean, sd, measure)

# Create a new dataframe to hold the means and standard deviations of shannon diversity estimates
SampleID <- row.names(shannon_rhizo)
mean <- apply(shannon_rhizo, 1, mean)
sd <- apply(shannon_rhizo, 1, sd)
measure <- rep("Shannon\nDiversity", nsamp_rhizo)
shan_stats_rhizo <- data.frame(SampleID, mean, sd, measure)

# Create a new dataframe to hold the means and standard deviations of simpson diversity estimates
SampleID <- row.names(simpson_rhizo)
mean <- apply(simpson_rhizo, 1, mean)
sd <- apply(simpson_rhizo, 1, sd)
measure <- rep("Simpson\nDiversity", nsamp_rhizo)
simp_stats_rhizo <- data.frame(SampleID, mean, sd, measure)


#combine richness and evenness into one dataframe
alpha_rhizo <- rbind(rich_stats_rhizo, even_stats_rhizo, shan_stats_rhizo, simp_stats_rhizo)
s_rhizo <- data.frame(sample_data(phy_scale_rhizo))
s_rhizo$SampleID <- SampleID
alphadivrhizo <- merge(alpha_rhizo, s_rhizo, by = "SampleID")
alphadivrhizo$PlantSpecies <- factor(alphadivrhizo$PlantSpecies, levels = c("Asyriaca","Acurassavica"))




#######################################
#Plot alpha div by plant species -rhizo

colourCount = 14
getPalette = colorRampPalette(brewer.pal(8, "Set1"))

alphadivrhizo.plant.plot <- ggplot(alphadivrhizo, aes(x = PlantSpecies, y = mean, color = PlantSpecies, group = PlantSpecies)) +
  geom_boxplot() +
  facet_grid(measure ~ PlantSpecies, scales = "free", space = "free_x") +
  theme(strip.text.y = element_text(size = 8)) +
  geom_point(aes(fill = PlantSpecies), size = 1, shape = 21) +
  scale_color_manual(values = getPalette(colourCount)) +
  scale_fill_manual(values = getPalette(colourCount)) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_blank())

alphadivrhizo.plant.plot




#########################################
#Plot alpha div by insect presence -rhizo

colourCount = 14
getPalette = colorRampPalette(brewer.pal(8, "Set1"))

alphadivrhizo.insect.plot <- ggplot(alphadivrhizo, aes(x = InsectPresence, y = mean, color = InsectPresence, group = InsectPresence)) +
  geom_boxplot() +
  facet_grid(measure ~ InsectPresence, scales = "free", space = "free_x") +
  theme(strip.text.y = element_text(size = 8)) +
  geom_point(aes(fill = InsectPresence), size = 1, shape = 21) +
  scale_color_manual(values = getPalette(colourCount)) +
  scale_fill_manual(values = getPalette(colourCount)) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_blank())

alphadivrhizo.insect.plot




#############################################################
#Alpha div signif by plant species and insect presence -rhizo

#Adjusting data frame 
alphadivrhizo.split <- split(alphadivrhizo, f = alphadivrhizo$measure)
alphadivrhizo.richness <- alphadivrhizo.split[["Richness"]]
alphadivrhizo.richness <- alphadivrhizo.richness[-c(3,4,7)]
colnames(alphadivrhizo.richness)[colnames(alphadivrhizo.richness) == "mean"] <- "Richness"

alphadivrhizo.inverse.simpson <- alphadivrhizo.split[["Inverse\nSimpson"]]
alphadivrhizo.inverse.simpson <- alphadivrhizo.inverse.simpson[-c(3,4,7)]
colnames(alphadivrhizo.inverse.simpson)[colnames(alphadivrhizo.inverse.simpson) == "mean"] <- "Inverse\nSimpson"

alphadivrhizo.shannon <- alphadivrhizo.split[["Shannon\nDiversity"]]
alphadivrhizo.shannon <- alphadivrhizo.shannon[-c(3,4,7)]
colnames(alphadivrhizo.shannon)[colnames(alphadivrhizo.shannon) == "mean"] <- "Shannon\nDiversity"

alphadivrhizo.simpson <- alphadivrhizo.split[["Simpson\nDiversity"]]
alphadivrhizo.simpson <- alphadivrhizo.simpson[-c(3,4,7)]
colnames(alphadivrhizo.simpson)[colnames(alphadivrhizo.simpson) == "mean"] <- "Simpson\nDiversity"

alphadivrhizo.merge.1 <- merge(alphadivrhizo.richness, alphadivrhizo.inverse.simpson, by.x=c("SampleID", "PlantSpecies", "InsectPresence"), by.y=c("SampleID", "PlantSpecies", "InsectPresence"))

alphadivrhizo.merge.2 <- merge(alphadivrhizo.shannon, alphadivrhizo.simpson, by.x=c("SampleID", "PlantSpecies", "InsectPresence"), by.y=c("SampleID", "PlantSpecies", "InsectPresence"))

alphadivrhizo.whole <- merge(alphadivrhizo.merge.1, alphadivrhizo.merge.2, by.x=c("SampleID", "PlantSpecies", "InsectPresence"), by.y=c("SampleID", "PlantSpecies", "InsectPresence"))







#Assigning unique list for my variables
InsectPresence <- unique(alphadivrhizo$InsectPresence)
PlantSpecies <- unique(alphadivrhizo$PlantSpecies)


ad_metrics.rhizo <- c("Richness", "Inverse\nSimpson", "Shannon\nDiversity", "Simpson\nDiversity")
for(m in ad_metrics.rhizo){
  print(m)
  aov_temp.rhizo.InsectPresenceOnly <- aov(get(m) ~ InsectPresence, data = alphadivrhizo.whole)
  summary(aov_temp.rhizo.InsectPresenceOnly)
  anova_summary.rhizo.InsectPresenceOnly <- as.data.frame(summary(aov_temp.rhizo.InsectPresenceOnly)[[1]])
  write.table(anova_summary.rhizo.InsectPresenceOnly, file = paste0("anova_rhizo_InsectPresenceOnly", ".txt"), sep = "\t", quote = FALSE)
  if (summary(aov_temp.rhizo.InsectPresenceOnly)[[1]][["Pr(>F)"]][[1]]){
    tukey_out.rhizo.InsectPresenceOnly <- TukeyHSD(aov_temp.rhizo.InsectPresenceOnly)
    tukey_out_df.rhizo.InsectPresenceOnly <- as.data.frame(tukey_out.rhizo.InsectPresenceOnly$InsectPresence)
    tukey_out_df.rhizo.InsectPresenceOnly$ad_metric.rhizo <- m
    if (exists("tukey_summary.rhizo.InsectPresenceOnly")) {
      tukey_summary.rhizo.InsectPresenceOnly <- rbind(tukey_summary.rhizo.InsectPresenceOnly, tukey_out_df.rhizo.InsectPresenceOnly)
    } else {
      tukey_summary.rhizo.InsectPresenceOnly <- tukey_out_df.rhizo.InsectPresenceOnly
    }
  }
}
tukey_summary.rhizo.InsectPresenceOnly$q.value <- p.adjust(tukey_summary.rhizo.InsectPresenceOnly$`p adj`, method = "BH")
write.table(tukey_summary.rhizo.InsectPresenceOnly, file = "tukey_summary.rhizo.InsectPresenceOnly.txt", sep = "\t", quote = FALSE)





ad_metrics.rhizo <- c("Richness", "Inverse\nSimpson", "Shannon\nDiversity", "Simpson\nDiversity")
for(m in ad_metrics.rhizo){
  print(m)
  aov_temp.rhizo.PlantSpeciesOnly <- aov(get(m) ~ PlantSpecies, data = alphadivrhizo.whole)
  summary(aov_temp.rhizo.PlantSpeciesOnly)
  anova_summary.rhizo.PlantSpeciesOnly <- as.data.frame(summary(aov_temp.rhizo.PlantSpeciesOnly)[[1]])
  write.table(anova_summary.rhizo.PlantSpeciesOnly, file = paste0("anova_rhizo_PlantSpeciesOnly", ".txt"), sep = "\t", quote = FALSE)
  if (summary(aov_temp.rhizo.PlantSpeciesOnly)[[1]][["Pr(>F)"]][[1]]){
    tukey_out.rhizo.PlantSpeciesOnly <- TukeyHSD(aov_temp.rhizo.PlantSpeciesOnly)
    tukey_out_df.rhizo.PlantSpeciesOnly <- as.data.frame(tukey_out.rhizo.PlantSpeciesOnly$PlantSpecies)
    tukey_out_df.rhizo.PlantSpeciesOnly$ad_metric.rhizo <- m
    if (exists("tukey_summary.rhizo.PlantSpeciesOnly")) {
      tukey_summary.rhizo.PlantSpeciesOnly <- rbind(tukey_summary.rhizo.PlantSpeciesOnly, tukey_out_df.rhizo.PlantSpeciesOnly)
    } else {
      tukey_summary.rhizo.PlantSpeciesOnly <- tukey_out_df.rhizo.PlantSpeciesOnly
    }
  }
}
tukey_summary.rhizo.PlantSpeciesOnly$q.value <- p.adjust(tukey_summary.rhizo.PlantSpeciesOnly$`p adj`, method = "BH")
write.table(tukey_summary.rhizo.PlantSpeciesOnly, file = "tukey_summary.rhizo.PlantSpeciesOnly.txt", sep = "\t", quote = FALSE)





ad_metrics.rhizo <- c("Richness", "Inverse\nSimpson", "Shannon\nDiversity", "Simpson\nDiversity")
for(i in InsectPresence){
  print(i)
  for(m in ad_metrics.rhizo){
    print(m)
    aov_temp.rhizo.InsectPresence <- aov(get(m) ~ PlantSpecies, data = subset(alphadivrhizo.whole, InsectPresence == i))
    summary(aov_temp.rhizo.InsectPresence)
    anova_summary.rhizo.InsectPresence <- as.data.frame(summary(aov_temp.rhizo.InsectPresence)[[1]])
    write.table(anova_summary.rhizo.InsectPresence, file = paste0("anova_rhizo_InsectPresence", ".txt"), sep = "\t", quote = FALSE)
    if (summary(aov_temp.rhizo.InsectPresence)[[1]][["Pr(>F)"]][[1]]){
      tukey_out.rhizo.InsectPresence <- TukeyHSD(aov_temp.rhizo.InsectPresence)
      tukey_out_df.rhizo.InsectPresence <- as.data.frame(tukey_out.rhizo.InsectPresence$PlantSpecies)
      tukey_out_df.rhizo.InsectPresence$InsectPresence <- i
      tukey_out_df.rhizo.InsectPresence$ad_metric.rhizo <- m
      if (exists("tukey_summary.rhizo.InsectPresence")) {
        tukey_summary.rhizo.InsectPresence <- rbind(tukey_summary.rhizo.InsectPresence, tukey_out_df.rhizo.InsectPresence)
      } else {
        tukey_summary.rhizo.InsectPresence <- tukey_out_df.rhizo.InsectPresence
      }
    }
  }
}
tukey_summary.rhizo.InsectPresence$q.value <- p.adjust(tukey_summary.rhizo.InsectPresence$`p adj`, method = "BH")
write.table(tukey_summary.rhizo.InsectPresence, file = "tukey_summary.rhizo.InsectPresence.txt", sep = "\t", quote = FALSE)





ad_metrics.rhizo <- c("Richness", "Inverse\nSimpson", "Shannon\nDiversity", "Simpson\nDiversity")
for(p in PlantSpecies){
  print(p)
  for(m in ad_metrics.rhizo){
    print(m)
    aov_temp.rhizo.PlantSpecies <- aov(get(m) ~ InsectPresence, data = subset(alphadivrhizo.whole, PlantSpecies == p))
    summary(aov_temp.rhizo.PlantSpecies)
    anova_summary.rhizo.PlantSpecies <- as.data.frame(summary(aov_temp.rhizo.PlantSpecies)[[1]])
    write.table(anova_summary.rhizo.PlantSpecies, file = paste0("anova_rhizo_PlantSpecies", ".txt"), sep = "\t", quote = FALSE)
    if (summary(aov_temp.rhizo.PlantSpecies)[[1]][["Pr(>F)"]][[1]]){
      tukey_out.rhizo.PlantSpecies <- TukeyHSD(aov_temp.rhizo.PlantSpecies)
      tukey_out_df.rhizo.PlantSpecies <- as.data.frame(tukey_out.rhizo.PlantSpecies$InsectPresence)
      tukey_out_df.rhizo.PlantSpecies$PlantSpecies <- p
      tukey_out_df.rhizo.PlantSpecies$ad_metric.rhizo <- m
      if (exists("tukey_summary.rhizo.PlantSpecies")) {
        tukey_summary.rhizo.PlantSpecies <- rbind(tukey_summary.rhizo.PlantSpecies, tukey_out_df.rhizo.PlantSpecies)
      } else {
        tukey_summary.rhizo.PlantSpecies <- tukey_out_df.rhizo.PlantSpecies
      }
    }
  }
}
tukey_summary.rhizo.PlantSpecies$q.value <- p.adjust(tukey_summary.rhizo.PlantSpecies$`p adj`, method = "BH")
write.table(tukey_summary.rhizo.PlantSpecies, file = "tukey_summary.rhizo.PlantSpecies.txt", sep = "\t", quote = FALSE)


