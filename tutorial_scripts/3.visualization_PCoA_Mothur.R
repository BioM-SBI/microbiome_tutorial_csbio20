#!/usr/bin/env Rscript

####################################
# Author: Pantakan Puengrang
# Modified Date: 21st October 2020
####################################

########################################################################
# Principal Coordinate Analysis plot Visualization using Mothur Results
########################################################################

# -----------------------------------
# Download Libraries
# -----------------------------------
install.packages('ggplot2')

# -----------------------------------
# Import library for visualization
# -----------------------------------
library(ggplot2)
library(dplyr)
library(RColorBrewer)

# -----------------------------------
# Set theme
# -----------------------------------
theme_set(theme_bw())

## Get microbiome transformation and visualization functions
source("tutorial_scripts/microbiome_scripts.R")

# Input files
coordinate = "tutorial_input/tutorial.pcoa.axes"
var_file = "tutorial_input/tutorial.pcoa.loadings"

# Prepare data for plotting (input from Mothur)
bray_dist = mothur_pcoa(coordinate, var_file)

# Mothur PCoA Plot
## Define manually the color-set as equal as the number of samples
colors = c("red", "blue", "green", "orange",  "pink","turquoise")

## Or use color RampPalette to automatically get the color-set
sample_count = length(bray_dist$coordinate$group)
getPallete = colorRampPalette(brewer.pal(12, "Paired"))
colors = getPallete(sample_count)

## If user would like to export as PNG de-hash line 46 (png()) and line 48 (dev.off())
# png('workshop_pcoa.png',units="in", width=6, height=4, res = 300)
plot_pcoa(bray_dist, colors)
# dev.off()

# Finish