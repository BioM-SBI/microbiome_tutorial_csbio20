#!/usr/bin/env Rscript

####################################
# Author: Pantakan Puengrang
# Modified Date: 21st October 2020
####################################

########################################################################
# Bar plot Visualization with Top N taxa for dominant microbial profile
########################################################################

# -----------------------------------
# Import library for visualization
# -----------------------------------
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(RColorBrewer)

# -----------------------------------
# Set theme
# -----------------------------------
theme_set(theme_bw())

## Get microbiome transformation and visualization functions
source("tutorial_script/microbiome_scripts.R")

# File inputs
shared_file = 'tutorial_inputs/tutorial.0.03.norm.shared'
taxonomy_file = 'tutorial_inputs/tutorial.taxonomy'

# Import shared table from Mothur
transposed_shared = import_sharedTable(shared_file)

# Import Taxonomic data
taxo = import_taxonomyFile(taxonomy_file)

# Merged Taxonomy (at specific rank) and abundance table into taxa table
# Preferred taxa rank
rank = 'Genus'
anno_rank = taxo$Genus

# Create Preferred rank abundance (For example Phylum abundance)
group_taxa_data = combine_df(transposed_shared, anno_rank)

# Find top N taxa of each sample and revealed as dominant microbial profile
# Input is dataframe of specific rank taxa and top (N) taxa were define as input
dominant_data = filter_top_n_taxa(group_taxa_data, N=10)
# Add others at the bottom of the data.frame
dominant_data = calculate_other(dominance=dominant_data)

# Transform data for visualization
tabular_df = transform_df(dominant_data)

# Set color palatte
taxCount = length(rownames(dominant_data))
getPallete = colorRampPalette(brewer.pal(12, "Paired"))
color_set = getPallete(taxCount)

png('tutorial.barPlot_top_10_taxa_Mothur.png', width=3000, height=2800, res=300)
plot_bar(tabular_df, rank, color_set)
dev.off()
####### -------- Finish -------- #######