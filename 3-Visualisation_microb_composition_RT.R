# 14/04/2024
# Katerina Katirtzoglou 
# Rainbow trout (RT) microbiome analysis on 16S rRNA gene amplicon data

#------------------------------------------------------------------------------
# After quality control reads were processed with dada2 and then curated with 
# lulu. Here we use R packages "phyloseq" and "microbiome" to explore fundamental 
# diversity measures and microbiome community composition
#------------------------------------------------------------------------------

# Load packages
library(phyloseq)
library(ggplot2)
library(tidyr)
library(microbiome)
library(VennDiagram)
library(scales)
library(plotly)

# Set working directory
setwd("/Users/aikaterinikatirtzoglou/Desktop/Alaskan_RT/")

## SORT AND PREPARE COMPONENTS FOR PHYLOSEQ

# Import metadata
metadata_trout <- read.csv("metadata.csv", header = TRUE, row.names =  1)

# Modify sample names in metadata so they are the same with the taxonomy table
row.names(metadata_trout) <- gsub("X","RT", rownames(metadata_trout))

# Remove sample column (1st) because the row names will be the sample names in the metadata
metadata_trout <- metadata_trout[,-1]

# Import curated tables of counts and taxonomy
curated_lulu_counts_df <- read.csv("Curated_Table.csv", header = TRUE, row.names =  1)

# Naming samples "RT"+index as in the metadata
colnames(curated_lulu_counts_df) <- gsub("X", "RT", colnames(curated_lulu_counts_df))

# Sort samples in ascending order RT1-RT20, extract numeric part from colnames of the ASV count table
numeric_part_cols <- as.numeric(sub("^RT([0-9]+)$", "\\1", colnames(curated_lulu_counts_df)))
sorted_cols <- colnames(curated_lulu_counts_df)[order(numeric_part_cols)]

# Sort ASV count table colnames based on numeric part
curated_lulu_counts_df <- curated_lulu_counts_df[, order(numeric_part_cols)]

# Reorder the columns of the ASV counts table
curated_lulu_counts_df <- curated_lulu_counts_df[,sorted_cols]

# Import taxonomy table of ASVs
curated_lulu_taxonomy_df <- read.csv("Curated_Tax.csv", header = TRUE, row.names =  1)

# Make the two tables into matrices
curated_lulu_counts_matrix <- as.matrix(curated_lulu_counts_df)
curated_lulu_taxonomy_matrix <- as.matrix(curated_lulu_taxonomy_df)

# PHYLOSEQ

# We are ready to create a phyloseq object
otu = otu_table(curated_lulu_counts_matrix, taxa_are_rows = TRUE)
taxonomy = tax_table(curated_lulu_taxonomy_matrix)
metadata = sample_data(metadata_trout)

RT_physeq = phyloseq(otu, taxonomy, metadata)



#---- Exploring our dataset ----#

# Summarize our dataset
microbiome::summarize_phyloseq(RT_physeq)

# Explore what percent of the total reads accounts for each ASV 
percent((sort(taxa_sums(RT_physeq))/(sum(taxa_sums(RT_physeq)))))


#----- Table 1 -------#

# Extract a table with the alpha diversity measures
alpha_div <- estimate_richness(RT_physeq, measures = c("Shannon", "Simpson", "Observed"))

# Calculate the number of ASV reads per sample
Number_of_reads <- data.frame(colSums(curated_lulu_counts_df))
colnames(Number_of_reads) <- c("Number of reads")

# Merge information on indices and number of reads in one dataframe
alpha_div <- cbind(Number_of_reads,alpha_div)

# Round the decimals
alpha_div$Simpson <- round(alpha_div$Simpson, digits = 3)
alpha_div$Shannon <- round(alpha_div$Shannon, digits = 3)

# Convert alpha diversity to a data frame
alpha_div_df <- as.data.frame(alpha_div)

# Extract table in a .csv file
write.csv(alpha_div_df, file = "alpha_diversity_table.csv")

# Also extracting stats from sequencing results and ASV table
mean(colSums(curated_lulu_counts_df))
sd(colSums(curated_lulu_counts_df))


#------- FIGURE 1 ------- #

# Create a simple rank abundace plot
par(mar = c(10, 4, 4, 2) + 0.1) # make more room on bottom margin
barplot(sort(taxa_sums(RT_physeq), TRUE)[1:10]/nsamples(RT_physeq), las=2)


#------- FIGURE 2 ------- #

# Extract the desired taxa from the tax_table
desired_taxa <- tax_table(RT_physeq)[, "Genus"]

# Sort by the taxa sums
sorted_taxa_sums <- sort(taxa_sums(RT_physeq), decreasing = TRUE)

# Extract the top 10 taxa
top_taxa <- sorted_taxa_sums[1:10]

# Calculate the relative abundance
relative_abundance <- top_taxa / nsamples(RT_physeq)

# Replace ASV names with Genus names
names(relative_abundance) <- desired_taxa[names(relative_abundance)]

# Define colors
custom_colors <- c("blue", "green", "red", "orange", "purple", "yellow", "cyan", "magenta", "brown", "gray")

# Create the bar plot 
top10_genus_barplot <- ggplot(data = data.frame(taxa = names(relative_abundance), abundance = relative_abundance), aes(x = reorder(taxa, -abundance), y = abundance, fill = taxa)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = custom_colors, name = "Genus") +  # Custom colors
  #theme_bw() +  # Change theme
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Angle x-axis labels
  labs(x = "Genera", y = "Abundance")


#------- Figure 3a - compostion plot of ALL samples at Phylum level---------#

# Collapsing ASVs to Phylum level counts table
composition_table_total <- RT_physeq  %>% aggregate_taxa(level = "Phylum")

# Creating a dataframe with the 15 phyla and their count data
composition_for_pie <- as.data.frame(otu_table(composition_table_total))

# Removing non biological samples
composition_for_pie<- composition_for_pie[,-c(21:24)]

# Creating a dataframe adding the counts of all samples per taxon
composition_for_pie <- as.data.frame(rowSums(composition_for_pie))

# Calculate total number of reads
total_reads <- sum(composition_for_pie)

# Calculate relative abundance for each phylum
composition_for_pie$Relative_Abundance <- composition_for_pie[,1] / total_reads

# Identify the five most abundant phyla
top_phyla <- head(composition_for_pie[order(-composition_for_pie$Relative_Abundance),], 5)

# Removing unromalized counts
top_phyla = top_phyla[-1]

# Calculate the abundance of the "Other" category
other_abundance <- sum(composition_for_pie$Relative_Abundance) - sum(top_phyla$Relative_Abundance)

# Create a new dataframe with the top phyla and the "Other" category
pie_data <- rbind(top_phyla, other_abundance)

# Set the row name for the "Other" category
rownames(pie_data)[nrow(pie_data)] <- "Other"


pie_phylum <- plot_ly(data = pie_data, 
                      labels = ~rownames(pie_data), 
                      values = ~Relative_Abundance,
                      textposition = c("inside", "inside", "inside", "outside", "outside", "outside"), # Adjust the position for each label
                      textinfo = "percent",
                      hoverinfo = "label",
                      outsidetextfont = list(color = "black"),
                      marker = list(colors = c("orange", "blue", "#FCE205", "green", "black", "red")),
                      type = "pie") %>%
  layout(legend = list(orientation = "v", x = 1, y = 0.5)) # Align legend to the right of the pie chart


#-------Figure 3b - compostion plot of ALL samples at the Genus level---------#

# Collapsing ASVs to Phylum level
composition_table_total <- RT_physeq  %>% aggregate_taxa(level = "Genus")

# Creating a dataframe with the 15 phyla and their count data
composition_for_pie <- as.data.frame(otu_table(composition_table_total))

# Removing non biological samples
composition_for_pie<- composition_for_pie[,-c(21:24)]

# Creating a dataframe adding the counts of all samples per taxon
composition_for_pie <- as.data.frame(rowSums(composition_for_pie))

# Calculate total number of reads
total_reads <- sum(composition_for_pie)

# Calculate relative abundance for each phylum
composition_for_pie$Relative_Abundance <- composition_for_pie[,1] / total_reads

# Identify the five most abundant phyla
top_phyla <- head(composition_for_pie[order(-composition_for_pie$Relative_Abundance),], 5)

# Removing unromalized counts
top_phyla = top_phyla[-1]

# Calculate the abundance of the "Other" category
other_abundance <- sum(composition_for_pie$Relative_Abundance) - sum(top_phyla$Relative_Abundance)

# Create a new dataframe with the top phyla and the "Other" category
pie_data_genus <- rbind(top_phyla, other_abundance)

# Set the row name for the "Other" category
rownames(pie_data_genus)[nrow(pie_data)] <- "Other"

# make pie chart
pie_genus <- plot_ly(data = pie_data_genus, 
                     labels = ~rownames(pie_data_genus), 
                     values = ~Relative_Abundance,
                     textposition = c("inside", "inside", "inside", "outside", "outside", "outside"), # Adjust the position for each label
                     textinfo = "percent",
                     hoverinfo = "label",
                     outsidetextfont = list(color = "black"),
                     marker = list(colors = c("orange", "blue", "#FCE205","black" , "green", "red")),
                     type = "pie") %>%
  
  layout(legend = list(orientation = "v", x = 1, y = 0.5)) # Align legend to the right of the pie chart



#------- FIGURE 4 ------- #

## Make a venn plot that shows how many ASVs arw shared between biological
##  groups and controls

# Identify column indices for biological and control samples
biological_indices <- 1:20
control_indices <- (ncol(curated_lulu_counts_df) - 3):ncol(curated_lulu_counts_df)

# Extract ASV counts for biological and control samples
biological_counts <- curated_lulu_counts_df[, biological_indices]
control_counts <- curated_lulu_counts_df[, control_indices]

# Make a list of ASVs present in controls (Identify ASVs with zero counts in all control samples)
control_zero_asvs <- rownames(control_counts)[rowSums(control_counts) != 0]

# Make a list of ASVs present in biological samples (Identify ASVs with zero counts in all biological samples)
biological_zero_asvs <- rownames(biological_counts)[rowSums(biological_counts) != 0]


# Define custom colors for the sets
set_colors <- c("#FDBF6F","palegreen4")

# Create Venn diagram
venn.plot <- venn.diagram(
  x = list(Set1 = control_zero_asvs, Set2 = biological_zero_asvs),
  category.names = c("Control", "Biological"),
  fill = set_colors,  # Set custom colors
  title = "",         # Remove title
  cat.pos = c(-20,52),  # Align category names at the same height
  filename = NULL
)


#------- Table S4 ------- #

# Finding which ASVs are only present in controls 
control_specific_ASVs <- setdiff(control_zero_asvs,biological_zero_asvs)

# Finding taxa only present in controls
control_specific_taxa <- curated_lulu_taxonomy_df[rownames(curated_lulu_taxonomy_df) %in% control_specific_ASVs, ]

# Save in a .csv file
write.csv(control_specific_taxa , file = "control_specific_taxa.csv")


#------- Figure 5 ------- #

control_specific_taxa <- control_specific_taxa %>%
  group_by(Phylum, Genus) %>%
  summarise(ASV_Count = n(), .groups = 'drop')

# Summarize ASV_Count by Phylum
phylum_order <- control_specific_taxa %>%
  group_by(Phylum) %>%
  summarise(total_ASV_Count = sum(ASV_Count)) %>%
  arrange(total_ASV_Count) %>%
  pull(Phylum)

# Reorder Phylum in the original dataframe
control_specific_taxa$Phylum <- factor(control_specific_taxa$Phylum, levels = phylum_order)

# Make a palette
fig5_palette <- c("orange", "blue", "#FCE205", "green", "black", "red", "green4", "brown", "blue4", "purple3")

# Plot the data
ggplot(control_specific_taxa, aes(x = ASV_Count, y = reorder(Phylum, ASV_Count), fill = Genus)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  scale_fill_manual(values = fig5_palette) +
  labs(
    x = "",
    y = "") +
  theme(
    axis.text = element_text(size = 10),
    legend.text = element_text(size = 8)
  )


