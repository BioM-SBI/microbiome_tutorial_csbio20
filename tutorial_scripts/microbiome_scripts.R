#!/usr/bin/env Rscript

####################################
# Author: Pantakan Puengrang
# Modified Date: 21st October 2020
####################################

############################################################################################################
# <<Source Code>> functions for Mothur's restuls data transformation and visualization
############################################################################################################

# Function for Transform mothur rarefaction plot
import_mothur_rarefaction = function(file){
  ### ---------------------------------------------------------------------------------- ####
  ### This is a function to import rarefaction table from mothur (rarefaction.single())  ####
  ### ---------------------------------------------------------------------------------- ####
  
  # Import to R dataframe
  data = read.table(file, header = T)
  
  # Pick out the mean column that are not hci or lci
  data = data %>% select('numsampled', matches("^X"))
  
  # Find and replace X0.03. at the beginning of the colnames with nothing
  names = str_replace_all(colnames(data), '^X0.03.','')
  colnames(data) = names
  
  # Melt from matrix format to tabular format with melt() function
  transformed_tbl <- melt(data, id.vars = 'numsampled')
  return(transformed_tbl)
}

# Function for Transform QIIME2 rarefaction plot


# Function for Mothur shared file importing
import_sharedTable <- function(input_shared){
	### ------------------------------------------------------------ ####
	### This is a function to import OTU table (mothur.shared) file #####
	### and to remove column $label and $numOtus ###
	### ------------------------------------------------------------ ####

	# Import to R data.frame 
	shared = read.table(shared_file, header=TRUE)
	# Remove column 1st $label
	shared = shared[-1]
	# Remove column 2st $numOtus
	shared = shared[-2]

	# Transpose dataframe representings features as row and sampleID as column
	transposeTbl = as.data.frame(t(shared[-1]))
	# Rename column name as sample_IDs
	colnames(transposeTbl) = shared$Group
	return(transposeTbl)
}

# Function for Mothur taxonomy file importing
import_taxonomyFile <- function(input_taxonomy){
	### ---------------------------------------------------------------------------- #####
	### This is a function to import taxonomy annotation (mothur.cons.taxonomy) file #####
	### Remove all confident values in parenthesis (xxx) of taxonomic lineage in column 3rd ####
	### which xxx represents 0 -> 100
	### Remove semicolon (;) at the end of string in column 3rd
	### Split taxonomic lineage into individual rank (Domain to Genus)
	### ---------------------------------------------------------------------------- #####

	# Import to R data.frame 
	taxo = read.csv(input_taxonomy, sep='\t')
	
	# Prepare Mothur data before plotting
	taxRank = c('Domain','Phylum','Class','Order','Family','Genus')
	
	# Remove confident number after taxa name, 
	# for example: 
	# Bacteria(100);Chloroflexi(100) to Bacteria;Chloroflexi
	taxo$Taxonomy = str_remove_all(taxo$Taxonomy, pattern = "\\([^()]+\\)")
	
	# Remove the last semicolon sign of str, 
	# for example: 
	# Bacteria;chloroflexi; to Bacteria;chloroflexi
	taxo$Taxonomy = str_remove(taxo$Taxonomy, pattern = ";$")

	# Split Taxonomic lineage into indivisual rank
	taxo = separate(data=taxo, col=Taxonomy, into=taxRank, sep=';')
	return(taxo)
}
# Merge two dataframe 
combine_df <- function(transposeTbl, rank_anno){
	### ----------------------------------------------------------------- ####
	### Combine preferred taxonomy and shared table into single data frame ###
	### ----------------------------------------------------------------- ####

	###  Merged data
	combine_df = as.data.frame(cbind(transposeTbl, taxRank=rank_anno))
	### Group by phylum and summarize values
	taxa_data = as.data.frame(combine_df %>%
		group_by(taxRank) %>%
		summarize_all(funs(sum)))

	# Set row names as Phylum name and drop column $taxRank
	rownames(taxa_data) = taxa_data$taxRank
	taxa_data = taxa_data[-1]

	return(taxa_data)
}
# Filter dominant taxa with User-defined cutoff
filter_dominance <- function(tbl, cutoff){
	### ------------------------------------------------------------------------------------ ####
	### This function calculate relative abundance and filter dominant features (OTU or ASV) ###
	### ------------------------------------------------------------------------------------ ####
	##### Herein, proportion of each feature in each sample were calculate and filter
	##### Dominant features are filtered which their proportion greater than defined cutoff
	### ------------------------------------------------------------------------------------ ####
	# Calculate relative abundance of each 
	relative_data = data.matrix(tbl, rownames.force= NA) %>% prop.table(margin=2)
	# Calculate each maximum abundance (from all samples) of each taxa
	max = apply(relative_data, 1, max)
	# Filter taxa which maximum abundance greater than cutoff
	max_name = names(which(max >= cutoff))
	# Filter dominant taxa to new data frame using name from cutoff
	dominance = relative_data[which(rownames(relative_data) %in% max_name), ]
	return(dominance)
}

filter_top_n_taxa = function(tbl,N=10){
	### ------------------------------------------------------------------------------------ ####
	### This function calculate relative abundance and filter dominant features (OTU or ASV) ####
	### ------------------------------------------------------------------------------------ ####
	##### Herein, proportion of each feature in each sample were calculate and
	##### Dominant features are filtered which observed as top highest N taxa in each sample 
	##### were included
	### ------------------------------------------------------------------------------------ ####
	# Calculate relative abundance of each 
	relative_data = data.matrix(tbl, rownames.force= NA) %>% prop.table(margin=2)
	# Filter taxa (from row names) which observed as top N taxa in each sample
	top_taxa = c()
	# Loop though each sample
	for (sample in 1:ncol(relative_data)) {
	  #print(sample)
	  # use index (variable sample) to select column and order descending (order(-df[,index])) 
	  top_taxa = append(top_taxa, rownames(relative_data[order(-relative_data[,sample]),])[1:N])
	}
	# Find identical taxa name
	top_taxa = unique(top_taxa)
	print(paste("There are",length(unique(top_taxa)),"taxa found as top", N ,"Taxa in this dataset"))
	top_n_domin = relative_data[which(rownames(relative_data) %in% top_taxa), ]
	return(top_n_domin)
}

# Calculate Other and Add to buttom row of domiannt table
calculate_other = function(dominance){
	### ------------------------------------------------------------------------------------ ####
	### Calculate Other 
	### ------------------------------------------------------------------------------------ ####
	##### Also, the remain gap were calculated which minus by 1 (q = 1-p) to show as Other
	##### It represents as all non-dominant taxa 
	# Calculate composition of Other taxa which defined as non-dominant taxa
	Other_taxa = 1 - colSums(dominance)
	dominance = rbind(dominance, Others=Other_taxa)
	# Change data type from matrix to datafrmae
	dominance = as.data.frame(dominance)
	print("Row Others were added at the bottom of the dataframe")
	return(dominance)
}

# Transform data for visualization
transform_df <- function(domin_df){
	### ------------------------------------------------------------------------------------ ####
	### This function transforms data for barplot visualization ###
	### ------------------------------------------------------------------------------------ ####
	
	# Transform table from matrix like format into tarbular format
	tabular_df = domin_df %>% gather('Sample_ID','Relative_abundance')
	# Add Columns 'Phylum' from rownames multiply by number of samples
	tabular_df = cbind(tabular_df, taxRank=rep(rownames(domin_df), ncol(domin_df)))
	# Set factor to fix order of plot
	tabular_df$taxRank = factor(tabular_df$taxRank, levels = rownames(domin_df))
	return(tabular_df)
}

# Plot bar 
plot_bar <- function(tabular_df, rank, color_set){
	###  -----------------------------------------------------------------------  ####
	###  This is a function to visualize Microbial community profile as bar plot  ####
	###  -----------------------------------------------------------------------  ####

	# Get ordered taxa labels
	group_ord = aggregate(tabular_df, list(tabular_df$taxRank), mean)
	group_ord = group_ord[-2]
	group_ord = group_ord[-3]
	group_ord = group_ord[order(-group_ord[2]),]
	ord_tax = array(group_ord[1:nrow(group_ord),]$Group.1)
	tmp = append(ord_tax[!grepl('Others|uncultured',ord_tax)],ord_tax[grepl('uncultured',ord_tax)])
	ord_tax = append(tmp, 'Others')

	# Set factor of order 
	tabular_df$taxRank = factor(tabular_df$taxRank, levels=rev(ord_tax))

	# plot
	plot = ggplot(tabular_df, aes(Sample_ID, Relative_abundance))
	plot + geom_bar(stat="identity", aes(fill=taxRank), color='black') +
	theme(axis.text.x = element_text(angle=45, vjust=1,hjust=1, size=10), 
        legend.position = 'right',
        legend.justification = 'top') +
	xlab('Samples') +
	ylab('Relative Abundance') +
	labs(fill = rank) +
	scale_fill_manual(values = color_set)
}

# Subset taxa annotation from dominant feature ID (OTU / ASV)
subset_taxa = function(domin_df, taxonomy, rowID = F){
        ### ------------------------------------------------------------------------------------ ####
        ### Subset taxa tables with dominant OTU
        ### ------------------------------------------------------------------------------------ ####
        taxo_tmp = taxonomy
        if(rowID == F){
        # if taxonomy table is not add OTU_ID as row name, then we add it
	  	row.names(taxo_tmp) <- taxo_tmp$OTU
	  	}
        Domin_otu = rownames(domin_df)
        dominOtuRepTaxa = subset(taxo_tmp, rownames(taxo_tmp) %in% Domin_otu)
        print(paste('There are',length(Domin_otu),'dominant features'))
        return(dominOtuRepTaxa)
}


# Create heatmap labels in the subsequent figure
heatmap_label = function(taxa_array, domin_df, separator){
        ### ------------------------------------------------------------------------------------ ####
        # Get features representative name with include featureID (OTU/ASV) with Taxa name
        ### ------------------------------------------------------------------------------------ ####
        
        # Get OTU name from data frame
        otuName <- rownames(domin_df)
        # Combined features and taxa name
        df <- data.frame(otuName, taxa_array)
        # Create new label
        # For example; OTU0001_genusABC, OTU0002_genusXYZ
        fullNewName <- paste(df$otuName, df$taxa_array, sep=separator)
        return(fullNewName)
}

# Calculate ecological distance using Bray-Curtis
anno_distance = function(original_shared, dominant_df, vegan_distance){
        ### ------------------------------------------------------------------------------------ ####
        ### Calculate Ecological Distance (dissimilarity) using Vegan package
        ### to create dendro gram in heatmap visualization
        ### ------------------------------------------------------------------------------------ ####
        
        # Calculate distance of taxa using their sample profile
        taxa_dist = vegdist(dominant_df, vegan_distance)
        # Calculate distance od samples using their microbial profile 
        sample_dist = vegdist(t(original_shared), vegan_distance)
        # Export as a list 
        dist_anno = list(tax_dist=taxa_dist, 
                    sam_dist=sample_dist)
        return(dist_anno)
}

# Mothur PCoA Plot
mothur_pcoa = function(file_coor, file_var){
	### -------------------------------------------------------------------- ####
	### This is a function to prepare data mothur.pcoa) and (mothur.loading) ####
	### to visualize as principal coordinate analysis plot 					 ####
	### -------------------------------------------------------------------- ####

	# Import Axes coordinate of PCo1 and PCo2
	bray_pcoa <- read.table(file=file_coor, header=T)
	bray <- bray_pcoa[1:3]

	# Import Variance capture of axis 1 and axis 2
	var_axes = read.table(file = file_var, sep='\t', header=T)
	var_axis1 = round(var_axes[1,2], 2)
	var_axis2 = round(var_axes[2,2], 2)
	mothur_pcoa_mat = list(coordinate=bray, 
						variance1=var_axis1,
						variance2=var_axis2)
	print(paste("The number of samples is",length(rownames(bray)), "samples"))
	return(mothur_pcoa_mat)
}

plot_pcoa = function(pcoa_mat, color_set){
	###  ---------------------------------------------------------------------  ####
	###  This is a function to visualize as principal coordinate analysis plot  ####
	###  From Mothur results												    ####
	###  ---------------------------------------------------------------------  ####
	plot.bray <- ggplot(pcoa_mat$coordinate, aes(axis1,axis2))
	plot.bray + geom_point(size=5, alpha=1, aes(col=group))+
		theme(legend.position = "right") +
		xlab(paste("PCo1 (", pcoa_mat$variance1, "% variation)")) + 
		ylab(paste("PCo2 (", pcoa_mat$variance2,"% variation)")) +
		labs(col = "Sources of sludge") +
		scale_color_manual(values = color_set)
}

####### ----------------------------------------------- ####
#######            QIIME2 data transformation			####
####### ----------------------------------------------- ####

## # Prepare QIIME2 data 
import_qiimeTable = function(qiime2_file, transpose_tbl=F){
  ### ------------------------------------------------------------ ####
  ### This is a function to import QIIME2 collapse file #####
  ### and to remove row 1 representing file title ###
  ### ------------------------------------------------------------ ####
  table = read.csv(qiime2_file, sep ='\t', row.names = 'X.OTU.ID', skip=1)
  # Check Transpose for different purpose
  if(transpose_tbl == F){
  	# Not transpose for Barplot and heatmap visualization
  	return(table)
  	} else {
  	# set transpose_tbl = T; Transpose table for PCoA calculation
  	return(t(table))
  	}
}

# Import metadata design for QIIME2
import_qiimeMeta = function(met, sampleTitle){
  ### ------------------------------------------------------------ ####
  ### This is a function to import metadata for plot annotation    ####
  ### ------------------------------------------------------------ ####
  # Import data
  meta_table = read.table(qiime2_metadata, sep='\t', header = TRUE, row.names = sampleTitle)
  return(meta_table)
}

# Calculate dissimilarity matrix
calculate_bray = function(table, met, method = 'bray'){
  ### ---------------------------------------------------------------------- ####
  ### This is a function to calculate bray-curtis dissimilarity using vegan  ####
  ### ---------------------------------------------------------------------- ####
  # Input is OTU table and preferred method
  dist_mat = vegdist(table, dist = method)
  # Ordination calculation
  betaDist = betadisper(dist_mat, met)
  # Calculate eigenvector
  eigenperc = round(as.vector(eigenvals(betaDist) / sum(eigenvals(betaDist)[eigenvals(betaDist) > 0]))* 100, digits=2)
  veg = list(matrix = dist_mat,
             beta_distance = betaDist,
             variance = eigenperc)
  return(veg)
}

plot_qiime2_pcoa = function(veganInfo, plotTitle = 'Principal Coordinate Analysis Plot',
                            elp = F, cvHull = F, seg = T){
  ### ---------------------------------------------------------------------- ####
  ### This is a function to visualize PCoA  ####
  ### ---------------------------------------------------------------------- ####
  # Parameters
  # elp = ellipse - Radius among group
  # cvHull = hull - Convex hull of the group
  # seg = segment - Linkage between centroid group and individual sample
  # Plot
  plot(veganInfo$beta_distance,
       ellipse = elp, # drop radius among group
       hull = cvHull, # F drop convex hull
       cex = 1,
       main = plotTitle,
       segment = seg, # Show linked of individual to the centroid of the group
       xlab = paste('PCoA1 (',veganInfo$variance[1],'% variation)'),
       ylab = paste('PCoA1 (',veganInfo$variance[2],'% variation)'))
}

# Function for Mothur taxonomy file importing
import_qiimeTaxonomy <- function(input_qiime2_collapse){
  ### ---------------------------------------------------------------------------- #####
  ### This is a function to import taxonomy annotation (collapse data) file        #####
  ### filter and replace
  ### Split taxonomic lineage into individual rank (Domain to Genus)               #####
  ### 
  ### ---------------------------------------------------------------------------- #####
  
  # Import to R data.frame 
  table = read.csv(input_qiime2_collapse, sep ='\t', skip=1)
  # Set Taxonomic lineage as row name
  rownames(table) = table$X.OTU.ID
  
  # Remove taxonomic assigned prefix in the name, 
  # for example: 
  # d__Bacteria;p__Chloroflexi to Bacteria;Chloroflexi
  table$X.OTU.ID = str_remove_all(table$X.OTU.ID, pattern = "[dpcofg]__")
  
  # Find and Replace __ with 'Unassigned, 
  # for example: 
  # Bacteria;__; to Bacteria;Unassigned
  table$X.OTU.ID = str_replace_all(table$X.OTU.ID, pattern = '__', replacement = 'Unassigned')
  
  # Split Taxonomic lineage into indivisual rank
  ## Not remove original
  # Prepare Mothur data before plotting
  taxRank = c('Domain','Phylum','Class','Order','Family','Genus')
  table = separate(data=table, col=X.OTU.ID, into=taxRank, sep=';', remove=FALSE)
  
  ## Extract only taxonomic information (Taxonomic Lineage, Phylum, Class , ... , Genus)
  taxo = table %>% select(X.OTU.ID,taxRank)
  ## Change column name 'X.OTU.ID' to 'OTU'
  colnames(taxo) = c('OTU',taxRank)
  return(taxo)
}


