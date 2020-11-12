#!/usr/bin/env Rscript

####################################
# Author: Pantakan Puengrang
# Modified Date: 21st October 2020
####################################

####################################
# Heatmap Visualization by pheatmap
####################################

# -----------------------------------
# Import library for visualization
# -----------------------------------
library(pheatmap)
library(gplots)
library(tidyr)
library(vegan)
library(RColorBrewer)
library(stringr)

# -----------------------------------
# Set theme
# -----------------------------------
theme_set(theme_bw())

## Get microbiome transformation and visualization functions
source("tutorial_scripts/microbiome_scripts.R")

# Define input files
shared_file = 'tutorial_input/tutorial.0.03.norm.shared'
taxonomy_file = 'tutorial_input/tutorial.taxonomy'

# Import shared table from Mothur
transposed_shared = import_sharedTable(shared_file)

# Import Taxonomic data
taxo = import_taxonomyFile(taxonomy_file)

# Remove taxa which abundance in all samples less than 1%
cutoff = 0.01
domin_data = filter_dominance(transposed_shared, cutoff)

# Create dominant taxa annotation
domin_taxa = subset_taxa(domin_data, taxo)

# Rename row-names of dominant data frame for heatmap visualization
# Note: Change Taxonomic Rank to (Domain, Phylum, Class, Order, Family, Genus)
#       At domina_taxa$____ >>>>>>>>>>>>>>>>>>>>>>V
rownames(domin_data) = heatmap_label(domin_taxa$Genus, domin_data, "_")

# Create list of samples and dominant taxa profile with 'Bray-Curtis disimilarity index'
# Distance metric could be change as 
# ("bray", "kulczynski", "jaccard", "gower", "altGower", "morisita", "horn", "mountford", "raup", "binomial", "chao", "cao" or "mahalanobis")
dendrogram = anno_distance(transposed_shared, domin_data, 'bray')

# Write grapphically heatmap by shade the color (yellow to red)
scaleColor <- colorRampPalette(c("lightyellow","red"),space ="rgb")(100)

# Plot Heatmap
#png('workshop.heatmap.genusAnno.png',width=2400, height=2000, res=300)
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
#dev.off()

