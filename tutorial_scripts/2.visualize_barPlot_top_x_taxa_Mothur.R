#!/usr/bin/env Rscript

####################################
# Author: Pantakan Puengrang
# Corresponding Author: Kanthida Kusonmano
# Affiliation: Systems Biology and Bioinformatics Research Group, KMUTT, Thailand
# Used for: Beginner's guide to microbiome analysis tutorial, CSBio2020
# Reference: Uthaipaisanwong, P., et al. (2020). “Beginner’s guide to microbiome analysis: Bioinformatics guidelines and practical concepts for amplicon-based microbiome analysis” In CSBio’20: Proceedings of the Eleventh International Conference on Computational Systems-Biology and Bioinformatics (CSBio2020), November 19– 21, 2020, Bangkok, Thailand. ACM, New York, NY, USA, 3 pages. https: //doi.org/10.1145/3429210.3429211
# Last modified Date: November 17, 2020
####################################

########################################################################
# Microbial profiles (barplots) visualization of Top N dominant taxa 
########################################################################

# -----------------------------------
# Load library for visualization
# -----------------------------------
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(RColorBrewer)

# -----------------------------------
# Set visulaization theme
# -----------------------------------
theme_set(theme_bw())

## Call the functions for microbiome analysis
source("tutorial_script/microbiome_scripts.R")

# Set path of input files (from mothur)
shared_file = 'tutorial_inputs/tutorial.0.03.norm.shared'
taxonomy_file = 'tutorial_inputs/tutorial.taxonomy'

# Retrieve abundance table file 
transposed_shared = import_sharedTable(shared_file)

# Retrieve taxonomic data file
taxo = import_taxonomyFile(taxonomy_file)

# Merged taxonomy (at a specific rank) and abundance tables
# Taxonomic rank: Domain, Phylum, Class, Order, Family, Genus
rank = 'Genus'
anno_rank = taxo$Genus

# Create abundance table of the selected rank
group_taxa_data = combine_df(transposed_shared, anno_rank)

# Get top N taxa of each sample to plot as dominant microbial profiles
dominant_data = filter_top_n_taxa(group_taxa_data, N=10)
# Add non-top N taxa as 'others' to the calculated profiles
dominant_data = calculate_other(dominance=dominant_data)

# Prepare data for visualization
tabular_df = transform_df(dominant_data)

# Set color palatte
taxCount = length(rownames(dominant_data))
getPallete = colorRampPalette(brewer.pal(12, "Paired"))
color_set = getPallete(taxCount)

png('tutorial.barPlot_top_10_taxa_Mothur.png', width=3000, height=2800, res=300)
plot_bar(tabular_df, rank, color_set)
dev.off()
