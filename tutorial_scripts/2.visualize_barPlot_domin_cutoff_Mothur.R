#!/usr/bin/env Rscript

####################################
# Author: Pantakan Puengrang
# Modified Date: 21st October 2020
####################################

####################################
# Bar plot Visualization from Mothur results
####################################

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
source("tutorial_scripts/microbiome_scripts.R")

# File inputs
shared_file = 'tutorial_inputs/tutorial.0.03.norm.shared'
taxonomy_file = 'tutorial_inputs/tutorial.taxonomy'

# Import shared table from Mothur
transposed_shared = import_sharedTable(shared_file)

# Import Taxonomic data
taxo = import_taxonomyFile(taxonomy_file)

# Merged Taxonomy (at specific rank) and abundance table into taxa table
# Preferred taxa rank
rank = 'Phylum'
anno_rank = taxo$Phylum

# Create Preferred rank abundance (For example Phylum abundance)
group_taxa_data = combine_df(transposed_shared, anno_rank)

# Filter to receive dominant taxa
cutoff = 0.01
dominant_data = calculate_other(filter_dominance(group_taxa_data, cutoff))

# Transform data for visualization
tabular_df = transform_df(dominant_data)

# Set color palatte
taxCount = length(rownames(dominant_data))
getPallete = colorRampPalette(brewer.pal(12, "Paired"))
color_set = getPallete(taxCount)

png('tutorial.barPlot_Domain_cutoff_Mothur.png', width=3000, height=2800, res=300)
plot_bar(tabular_df, rank, color_set)
dev.off()
####### -------- Finish -------- #######