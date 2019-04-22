#Load necessary packages and set seed for reproducibility

library(qiime2R)
library(vegan)
library(reshape2)
library(phyloseq)
library(ggplot2); packageVersion("ggplot2")
library(RColorBrewer)
library(data.table)
library(magrittr)
library(readr)

set.seed(1)




#Load phyloseq data set

setwd("C:/Users/warriorroad/Desktop/Insect_Microbial_Ecology_and_Evolution_Lab/Sequencing_data/MilkweedCombined_output/Exported_data")

#Use qiime2R to import artifacts and convert to phyloseqs object
metadata <- read_tsv("Milkweed_prelim_metadata.tsv")
metadata

ASVs<-read_qza("milkweed-combined-completely-filtered-table.qza")
names(ASVs)

taxonomy<-read_qza("milkweed-combined-taxonomy-without-spaces.qza")
taxonomy$uuid

tree<-read_qza("milkweed-combined-rooted-tree.qza")
tree$uuid

phy<-qza_to_phyloseq("milkweed-combined-completely-filtered-table.qza", "milkweed-combined-rooted-tree.qza", "milkweed-combined-taxonomy-without-spaces.qza","Milkweed_prelim_metadata.tsv", tmp="C:/tmp")
phy




#Subset feature table by insect and rhizosphere

otu_table(phy)[1:5,]

otumon<-otu_table(phy)[, 1:13]
otu_table(otumon)[1:5,]

oturhizo<-otu_table(phy)[, 14:41]
otu_table(oturhizo)[1:5,]




#Make monarch and rhizosphere phyloseq objects

physampledata<-sample_data(phy)

phytaxtable<-tax_table(phy)

phytree<-phy_tree(phy)

phymon <- phyloseq(otu_table(otumon, taxa_are_rows=TRUE), 
                   sample_data(physampledata), 
                   tax_table(phytaxtable),
                   phy_tree(phytree))


phyrhizo <- phyloseq(otu_table(oturhizo, taxa_are_rows=TRUE), 
                     sample_data(physampledata), 
                     tax_table(phytaxtable),
                     phy_tree(phytree))




#Examine ASV abundance

taxa_sum_df_mon <- data.frame(sum = taxa_sums(phymon))

head(taxa_sum_df_mon, 15)

smin <- min(taxa_sums(phymon))
smean <- mean(taxa_sums(phymon))
smax <- max(taxa_sums(phymon))

smin
smean
smax


taxa_sum_df_rhizo <- data.frame(sum = taxa_sums(phyrhizo))

head(taxa_sum_df_rhizo, 15)

smin <- min(taxa_sums(phyrhizo))
smean <- mean(taxa_sums(phyrhizo))
smax <- max(taxa_sums(phyrhizo))

smin
smean
smax




#Keep taxa found more than 3 times, for insects, and 10 times, for rhizo, in at least 5% of the samples
phy.filter.mon <- filter_taxa(phymon, function(x) sum(x > 3) > (0.05*length(x)), TRUE)
phy.filter.mon
phymon

phy.filter.rhizo <- filter_taxa(phyrhizo, function(x) sum(x > 10) > (0.05*length(x)), TRUE)
phy.filter.rhizo
phyrhizo




#Next we need to filter out samples with low counts- Histogram of sample sizes
sample_sum_df_mon <- data.frame(sum = sample_sums(phymon))

# Histogram of sample read counts -monarch
sample.read.mon <- ggplot(sample_sum_df_mon, aes(x = sum)) + 
  geom_histogram(color = "black", fill = "indianred", binwidth = 500) +
  ggtitle("Distribution of sample sequencing depth") + 
  xlab("Read counts") +
  theme(axis.title.y = element_blank())

sample.read.mon

smin <- min(sample_sums(phy.filter.mon))
smean <- mean(sample_sums(phy.filter.mon))
smax <- max(sample_sums(phy.filter.mon))

smin
smean
smax



sample_sum_df_rhizo <- data.frame(sum = sample_sums(phyrhizo))

# Histogram of sample read counts -rhizosphere
sample.read.rhizo <- ggplot(sample_sum_df_rhizo, aes(x = sum)) + 
  geom_histogram(color = "black", fill = "indianred", binwidth = 1000) +
  ggtitle("Distribution of sample sequencing depth") + 
  xlab("Read counts") +
  theme(axis.title.y = element_blank())

sample.read.rhizo

smin <- min(sample_sums(phy.filter.rhizo))
smean <- mean(sample_sums(phy.filter.rhizo))
smax <- max(sample_sums(phy.filter.rhizo))

smin
smean
smax




#Getting rid of samples with fewer than 500 reads for insects

phy.trim.mon <- prune_samples(sample_sums(phy.filter.mon)>=500, phy.filter.mon)
sum(otu_table(phy.filter.mon))
sum(otu_table(phy.trim.mon))
saveRDS(phy.trim.mon, "phy.trim.mon.rds")
sample_sum_df2_mon <- data.frame(sum = sample_sums(phy.trim.mon))

# Histogram of sample read counts
trim.read.mon <- ggplot(sample_sum_df2_mon, aes(x = sum)) + 
  geom_histogram(color = "black", fill = "indianred", binwidth = 500) +
  ggtitle("Distribution of sample sequencing depth") + 
  xlab("Read counts") +
  theme(axis.title.y = element_blank())

trim.read.mon




#Sample min/mean/max after removing low count samples for insects
smin <- min(sample_sums(phy.trim.mon))
smean <- mean(sample_sums(phy.trim.mon))
smax <- max(sample_sums(phy.trim.mon))

smin
smean
smax




#Save phyloseq objects and sample files

#Save dada2/phyloseq objects in case you need them again
saveRDS(phy.trim.mon, "phy.trim.mon.rds")
saveRDS(phy.filter.rhizo, "phy.filter.rhizo.rds")


# If you close and re-open R
phy.trim.mon <- readRDS("phy.trim.mon.rds") #how to read RD files back in
phy.filter.rhizo <- readRDS("phy.filter.rhizo.rds") #how to read RD files back in