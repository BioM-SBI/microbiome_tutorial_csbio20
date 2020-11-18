#!/usr/bin/env Rscript

####################################
# Author: Pantakan Puengrang
# Corresponding Author: Kanthida Kusonmano
# Affiliation: Systems Biology and Bioinformatics Research Group, KMUTT, Thailand
# Used for: Beginner's guide to microbiome analysis tutorial, CSBio2020
# Reference: Uthaipaisanwong, P., et al. (2020). “Beginner’s guide to microbiome analysis: Bioinformatics guidelines and practical concepts for amplicon-based microbiome analysis” In CSBio’20: Proceedings of the Eleventh International Conference on Computational Systems-Biology and Bioinformatics (CSBio2020), November 19– 21, 2020, Bangkok, Thailand. ACM, New York, NY, USA, 3 pages. https: //doi.org/10.1145/3429210.3429211
# Last modified Date: November 17, 2020
####################################

###########################################################################
# Principal Coordinate Analysis Visualization using Mothur Resulting File
###########################################################################

# -----------------------------------
# Load library for visualization
# -----------------------------------
library(ggplot2)
library(dplyr)
library(RColorBrewer)

# -----------------------------------
# Set visualization theme
# -----------------------------------
theme_set(theme_bw())

## Call the functions for microbiome analysis
source("tutorial_scripts/microbiome_scripts.R")

# Set path of input files
coordinate = "tutorial_inputs/tutorial.0.03.norm.braycurtis.0.03.lt.pcoa.axes"
var_file = "tutorial_inputs/tutorial.0.03.norm.braycurtis.0.03.lt.pcoa.loadings"

# Retrieve dissimilarity matrix for plotting PCoA (input from Mothur)
bray_dist = mothur_pcoa(coordinate, var_file)

# Mothur PCoA Plot
## Define a color set for each sample (manually)
colors = c("red", "blue", "green", "orange",  "pink","turquoise")

## Or use color RampPalette to automatically get a color set
sample_count = length(bray_dist$coordinate$group)
getPallete = colorRampPalette(brewer.pal(12, "Paired"))
colors = getPallete(sample_count)

png('tutorial.pcoa.png',units="in", width=6, height=4, res = 300)
plot_pcoa(bray_dist, colors)
dev.off()
