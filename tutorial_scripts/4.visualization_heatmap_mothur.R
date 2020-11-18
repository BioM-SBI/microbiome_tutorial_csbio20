#!/usr/bin/env Rscript

####################################
# Author: Pantakan Puengrang
# Corresponding Author: Kanthida Kusonmano
# Affiliation: Systems Biology and Bioinformatics Research Group, KMUTT, Thailand
# Used for: Beginner's guide to microbiome analysis tutorial, CSBio2020
# Reference: Uthaipaisanwong, P., et al. (2020). “Beginner’s guide to microbiome analysis: Bioinformatics guidelines and practical concepts for amplicon-based microbiome analysis” In CSBio’20: Proceedings of the Eleventh International Conference on Computational Systems-Biology and Bioinformatics (CSBio2020), November 19– 21, 2020, Bangkok, Thailand. ACM, New York, NY, USA, 3 pages. https: //doi.org/10.1145/3429210.3429211
# Last modified Date: November 17, 2020
####################################

####################################
# Heatmap Visualization 
####################################

# -----------------------------------
# Load library for visualization
# -----------------------------------
library(pheatmap)
library(gplots)
library(tidyr)
library(vegan)
library(RColorBrewer)
library(stringr)

# -----------------------------------
# Set visulaization theme
# -----------------------------------
theme_set(theme_bw())

## Call the functions for microbiome analysis
source("tutorial_scripts/microbiome_scripts.R")

# Set path of input files
shared_file = 'tutorial_inputs/tutorial.0.03.norm.shared'
taxonomy_file = 'tutorial_inputs/tutorial.taxonomy'

# Retrieve abundance table file 
transposed_shared = import_sharedTable(shared_file)

# Retrieve taxonomic data file
taxo = import_taxonomyFile(taxonomy_file)

# Set abundance cut-off value and remove taxa containing abundance in all samples less than 1%
cutoff = 0.01
domin_data = filter_dominance(transposed_shared, cutoff)

# Create dominant taxa annotation 
domin_taxa = subset_taxa(domin_data, taxo)
rownames(domin_data) = heatmap_label(domin_taxa$Genus, domin_data, "_")

# Create list of samples and dominant taxa profile with 'Bray-Curtis disimilarity index'
# Distance metric could be changed as ("bray", "kulczynski", "jaccard", "gower", "altGower", "morisita", "horn", "mountford", "raup", "binomial", "chao", "cao" or "mahalanobis")
dendrogram = anno_distance(transposed_shared, domin_data, 'bray')

# Create color scale for heatmap (yellow to red)
scaleColor <- colorRampPalette(c("lightyellow","red"),space ="rgb")(100)

# Plot Heatmap
png('tutorial.heatmap.png',width=2400, height=2000, res=300)
pheatmap(domin_data,
        cellwidth = 20,
        cellheight = 5,
        fontsize_row = 5,
        fontsize_col = 7,
        angle_col = 45,
        legend = T,
        fontsize = 5,
        color= scaleColor,
        clustering_method = 'ward.D',
        clustering_distance_rows = dendrogram$tax_dist,
        clustering_distance_cols = dendrogram$sam_dist
        )
dev.off()
