# Figure 1: Abundance peak
# load libraries
library(ggplot2)
library(tidyverse)

### Metagenome - DNA ###
# load all DNA classifications from the different depths
ten <- read.delim("classifications10.tsv", sep = '\t', header = TRUE) 
onehundred <- read.delim("classifications100.tsv", sep = '\t', header = TRUE)
onetwenty <- read.delim("classifications120.tsv", sep = '\t', header = TRUE)
onethirtyfive <- read.delim("classifications135.tsv", sep = '\t', header = TRUE)
onefifty <- read.delim("classifications150.tsv", sep = '\t', header = TRUE)
onesixtyfive <- read.delim("classifications165.tsv", sep = '\t', header = TRUE)
twohundred <- read.delim("classifications200.tsv", sep = '\t', header = TRUE)

# merge classification.tsv from all depths 
merged_data <- rbind(ten, onehundred, onetwenty, onethirtyfive, onefifty, onesixtyfive, twohundred)

# Plot Sample Column
merged_data <- merged_data %>%
  separate(Sample, into = c("Sample", "Depth", "extra"), sep = "_", remove = FALSE) %>%
  select(-extra) %>%  # Remove the 'extra' column
  mutate(Depth = as.numeric(gsub("m", "", Depth))) %>% # Remove 'm' and convert to numeric
  separate(Taxonomy, into = c("Root", "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";", fill = "right") %>%
  select(-Root)  # Remove the 'Root' column

NorB_data <- merged_data %>%
  filter(Marker == "NorB")

data_taxa <- NorB_data %>%
  group_by(Sample, Depth) %>%
  summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop")

ggplot(data_taxa, aes(x = Depth, y = Abundance)) +
  geom_bar(stat = "identity") + 
  labs(
    title = "NorB DNA Abundance Across Depth",
    x = "Depth (m)",
    y = "Abundance"
  )

ggsave("DNA_NorB_Abundance.jpg")


### Metatranscriptomics - RNA ###

# load all DNA classifications from the different depths
tenRNA <- read.delim("classifications10RNA.tsv", sep = '\t', header = TRUE) 
onehundredRNA <- read.delim("classifications100RNA.tsv", sep = '\t', header = TRUE)
onetwentyRNA <- read.delim("classifications120RNA.tsv", sep = '\t', header = TRUE)
onethirtyfiveRNA <- read.delim("classifications135RNA.tsv", sep = '\t', header = TRUE)
onefiftyRNA <- read.delim("classifications150RNA.tsv", sep = '\t', header = TRUE)
onesixtyfiveRNA <- read.delim("classifications165RNA.tsv", sep = '\t', header = TRUE)
twohundredRNA <- read.delim("classifications200RNA.tsv", sep = '\t', header = TRUE)

# merge classification.tsv from all depths 
merged_dataRNA <- rbind(tenRNA, onehundredRNA, onetwentyRNA, onethirtyfiveRNA, onefiftyRNA, onesixtyfiveRNA, twohundredRNA)

# Split the Sample column
merged_dataRNA <- merged_dataRNA %>%
  separate(Sample, into = c("Sample", "Depth", "extra"), sep = "_", remove = FALSE) %>%
  select(-extra) %>%  # Remove the 'extra' column
  mutate(Depth = as.numeric(gsub("m", "", Depth))) %>% # Remove 'm' and convert to numeric
  separate(Taxonomy, into = c("Root", "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";", fill = "right") %>%
  select(-Root)  # Remove the 'Root' column

NorB_RNAdata <- merged_dataRNA %>%
  filter(Marker == "NorB")

data_taxa_RNA <- NorB_RNAdata %>%
  group_by(Sample, Depth) %>%
  summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop")

ggplot(data_taxa_RNA, aes(x = Depth, y = Abundance)) +
  geom_bar(stat = "identity") + 
  labs(
    title = "NorB RNA Abundance Across Depth",
    x = "Depth (m)",
    y = "Abundance"
  )

ggsave("RNA_NorB_Abundance.jpg")



# Figure 2a + 2b: Alpha diversity
library("tidyverse")
library("vegan")

## Metatranscriptome(RNA) ##
# Load data
alpha_data_RNA <- read.csv("SI_TS_NorB_alpha_diversiy_RNA.csv")

# Split the placerun column
alpha_data_RNA <- alpha_data_RNA %>%
  separate(placerun, into = c("sample", "depth", "gene", "extra"), sep = "_", remove = FALSE) %>%
  select(-extra) %>%  # Remove the 'extra' column
  mutate(depth = as.numeric(gsub("m", "", depth)))  # Remove 'm' and convert to numeric
view(alpha_data_RNA)

# Select relevant columns
alpha_data_RNA_mod <- alpha_data_RNA %>%
  select(depth, phylo_entropy, quadratic, unrooted_pd, rooted_pd, bwpd)

# Pivot data to a longer format for faceting
alpha_data_RNA_longer <- alpha_data_RNA_mod %>%
  pivot_longer(cols = c("phylo_entropy", "quadratic", "unrooted_pd", "rooted_pd", "bwpd"),
               names_to = "metric",
               values_to = "value")

# Define your desired order of metrics
desired_order <- c("unrooted_pd", "rooted_pd", "phylo_entropy", "quadratic", "bwpd")

# Convert 'metric' column to a factor with specified levels
alpha_data_RNA_longer$metric <- factor(alpha_data_RNA_longer$metric, levels = desired_order)

ggplot(alpha_data_RNA_longer, aes(y = depth, x = value)) +
  geom_point(aes(colour = metric, shape = metric), size = 4) +
  scale_y_reverse() +
  labs(title = "Alpha Diversity Metric against Depth (Metatranscriptome)",
       x = "value") +
  theme_minimal() +
  facet_grid(. ~ metric, scales = "free_x")


## Metagenome(DNA) ##
# Load data
alpha_data_DNA <- read.csv("SI_TS_NorB_alpha_diversiy_DNA.csv")

# Split the placerun column
alpha_data_DNA <- alpha_data_DNA %>%
  separate(placerun, into = c("sample", "depth", "gene", "extra"), sep = "_", remove = FALSE) %>%
  select(-extra) %>%  # Remove the 'extra' column
  mutate(depth = as.numeric(gsub("m", "", depth)))  # Remove 'm' and convert to numeric
view(alpha_data_DNA)

# Select relevant columns
alpha_data_DNA_mod <- alpha_data_DNA %>%
  select(depth, phylo_entropy, quadratic, unrooted_pd, rooted_pd, bwpd)

# Pivot data to a longer format for faceting
alpha_data_DNA_longer <- alpha_data_DNA_mod %>%
  pivot_longer(cols = c("phylo_entropy", "quadratic", "unrooted_pd", "rooted_pd", "bwpd"),
               names_to = "metric",
               values_to = "value")

# Define your desired order of metrics
desired_order <- c("unrooted_pd", "rooted_pd", "phylo_entropy", "quadratic", "bwpd")

# Convert 'metric' column to a factor with specified levels
alpha_data_DNA_longer$metric <- factor(alpha_data_DNA_longer$metric, levels = desired_order)

ggplot(alpha_data_DNA_longer, aes(y = depth, x = value)) +
  geom_point(aes(colour = metric, shape = metric), size = 4) +
  scale_y_reverse() +
  labs(title = "Alpha Diversity Metric against Depth (Metagenome)",
       x = "value") +
  theme_minimal() +
  facet_grid(. ~ metric, scales = "free_x")



# Figure 2c + 2d: Beta diversity
# load packages
library(tidyverse)
library(pheatmap)

beta_data_DNA <- read.csv("SI_TS_NorB_beta_diversity_DNA.csv")
beta_data_RNA <- read.csv("SI_TS_NorB_beta_diversity_RNA.csv")

#### Beta Diversity for DNA ###
# Split 'sample_1' and 'sample_2' into their components 
beta_data_DNA <- beta_data_DNA %>%
  mutate(
    sample_1 = gsub("_NorB_complete_profile", "", sample_1),
    sample_2 = gsub("_NorB_complete_profile", "", sample_2)
  ) %>%
  mutate(
    sample_1 = gsub("SI072_", "", sample_1),
    sample_2 = gsub("SI072_", "", sample_2)
  ) %>%
  mutate(
    sample_1 = as.numeric(gsub("m", "", sample_1)),
    sample_2 = as.numeric(gsub("m", "", sample_2))
  )

# Generate a list of all unique samples from both columns
all_samples <- sort(unique(c(beta_data_DNA$sample_1, beta_data_DNA$sample_2)))

# Create an empty matrix filled with NA
beta_matrix_DNA_complete <- matrix(NA,
                                   nrow = length(all_samples),
                                   ncol = length(all_samples),
                                   dimnames = list(all_samples, all_samples))

# Fill in the matrix symmetrically
for (i in seq_len(nrow(beta_data_DNA))) {
  row_name <- as.character(beta_data_DNA$sample_1[i])
  col_name <- as.character(beta_data_DNA$sample_2[i])
  value <- beta_data_DNA$Z_1[i]
  
  # Assign value to both [row, col] and [col, row] to make it symmetrical
  beta_matrix_DNA_complete[row_name, col_name] <- value
  beta_matrix_DNA_complete[col_name, row_name] <- value
}

# Replace NA with 0 if desired (or leave as NA for missing values)
beta_matrix_DNA_complete[is.na(beta_matrix_DNA_complete)] <- 0

# Convert to a data frame for compatibility with pheatmap
beta_matrix_clean <- as.data.frame(beta_matrix_DNA_complete)

# View the completed matrix
print(beta_matrix_clean)

# Sort rows by their names or a specific column
sorted_matrix_DNA <- beta_matrix_clean[order(rownames(beta_matrix_clean)), , drop = FALSE]

# Cluster columns
col_clust <- hclust(dist(t(sorted_matrix_DNA)))  # Create column dendrogram
# Flip dendrogram branches if desired
col_clust$order <- order(rownames(beta_matrix_clean))

# Plot the heatmap
pheatmap(as.matrix(sorted_matrix_DNA),
         cluster_rows = FALSE,
         cluster_cols = col_clust,
         scale = "none",
         color = colorRampPalette(c("black", "gray", "white"))(100),
         border_color = NA,
         main = "KR Distance Heatmap with Dendrogram for DNA",
         fontsize_main = 25,
         filename = "Beta_DNA.jpg",
         width = 12)


#### Beta Diversity for RNA ###

# Split 'sample_1' and 'sample_2' into their components 
beta_data_RNA <- beta_data_RNA %>%
  mutate(
    sample_1 = gsub("_NorB_complete_profile", "", sample_1),
    sample_2 = gsub("_NorB_complete_profile", "", sample_2)
  ) %>%
  mutate(
    sample_1 = gsub("SI072_", "", sample_1),
    sample_2 = gsub("SI072_", "", sample_2)
  ) %>%
  mutate(
    sample_1 = as.numeric(gsub("m", "", sample_1)),
    sample_2 = as.numeric(gsub("m", "", sample_2))
  )

# Generate a list of all unique samples from both columns
all_samples_RNA <- sort(unique(c(beta_data_RNA$sample_1, beta_data_RNA$sample_2)))

# Create an empty matrix filled with NA
beta_matrix_RNA_complete <- matrix(NA,
                                   nrow = length(all_samples_RNA),
                                   ncol = length(all_samples_RNA),
                                   dimnames = list(all_samples_RNA, all_samples_RNA))

# Fill in the matrix symmetrically
for (i in seq_len(nrow(beta_data_RNA))) {
  row_name_RNA <- as.character(beta_data_RNA$sample_1[i])
  col_name_RNA <- as.character(beta_data_RNA$sample_2[i])
  value_RNA <- beta_data_RNA$Z_1[i]
  
  # Assign value to both [row, col] and [col, row] to make it symmetrical
  beta_matrix_RNA_complete[row_name_RNA, col_name_RNA] <- value_RNA
  beta_matrix_RNA_complete[col_name_RNA, row_name_RNA] <- value_RNA
}

# Replace NA with 0 if desired (or leave as NA for missing values)
beta_matrix_RNA_complete[is.na(beta_matrix_RNA_complete)] <- 0

# Convert to a data frame for compatibility with pheatmap
beta_matrix_clean_RNA <- as.data.frame(beta_matrix_RNA_complete)

# View the completed matrix
print(beta_matrix_clean_RNA)

# Sort rows by their names or a specific column
sorted_matrix_RNA <- beta_matrix_clean_RNA[order(rownames(beta_matrix_clean_RNA)), , drop = FALSE]

# Cluster columns
col_clust_RNA <- hclust(dist(t(sorted_matrix_RNA)))  # Create column dendrogram
# Flip dendrogram branches if desired
col_clust_RNA$order <- order(rownames(beta_matrix_clean_RNA))

# Plot the heatmap
pheatmap(as.matrix(sorted_matrix_RNA),
         cluster_rows = FALSE,
         cluster_cols = col_clust_RNA,
         scale = "none",
         color = colorRampPalette(c("black", "gray", "white"))(100),
         border_color = NA,
         main = "KR Distance Heatmap with Dendrogram for RNA", 
         fontsize_main = 25,
         filename = "Beta_RNA.jpg",
         width = 12)



# Figure 3: Bubble plot
library(tidyverse)
library(dplyr)

### Metagenome - DNA ###
# load all DNA classifications from the different depths
ten <- read.delim("classifications10.tsv", sep = '\t', header = TRUE) 
onehundred <- read.delim("classifications100.tsv", sep = '\t', header = TRUE)
onetwenty <- read.delim("classifications120.tsv", sep = '\t', header = TRUE)
onethirtyfive <- read.delim("classifications135.tsv", sep = '\t', header = TRUE)
onefifty <- read.delim("classifications150.tsv", sep = '\t', header = TRUE)
onesixtyfive <- read.delim("classifications165.tsv", sep = '\t', header = TRUE)
twohundred <- read.delim("classifications200.tsv", sep = '\t', header = TRUE)

# merge classification.tsv from all depths 
merged_data <- rbind(ten, onehundred, onetwenty, onethirtyfive, onefifty, onesixtyfive, twohundred)

# Split the Sample column
merged_data <- merged_data %>%
  separate(Sample, into = c("Sample", "Depth", "extra"), sep = "_", remove = FALSE) %>%
  select(-extra) %>%  # Remove the 'extra' column
  mutate(Depth = as.numeric(gsub("m", "", Depth))) %>% # Remove 'm' and convert to numeric
  separate(Taxonomy, into = c("Root", "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";", fill = "right") %>%
  select(-Root)  # Remove the 'Root' column

NorB_data <- merged_data %>%
  filter(Marker == "NorB")

## Order ##
# Sum Abundance by Sample, Order, and Depth
data_taxa_order <- NorB_data %>%
  group_by(Sample, Order, Depth) %>%
  summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
  mutate(Order = gsub("o__", "", Order))

# Make Bubble Plot
ggplot(data_taxa_order, aes(x = Order, y = Depth)) +
  geom_point(aes(size = Abundance, color = Abundance), alpha = 0.7) +  # Map color to Abundance
  scale_y_reverse(limits = c(210, 0)) +  # Shallow at top, deep at bottom
  scale_color_gradient(low = "red", high = "blue") +  # Gradient from red (low) to blue (high)
  scale_size(range = c(1, 10)) +
  theme_light() +
  labs(title = "DNA Abundance by Depth for NorB in Saanich Inlet Resolved to Order",
       x = "Order",
       y = "Depth",
       size = "Abundance") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# save plot
ggsave("NorB_BubblePlot_Order.jpg")


### Metatranscriptomics - RNA ###

# load all DNA classifications from the different depths
tenRNA <- read.delim("classifications10RNA.tsv", sep = '\t', header = TRUE) 
onehundredRNA <- read.delim("classifications100RNA.tsv", sep = '\t', header = TRUE)
onetwentyRNA <- read.delim("classifications120RNA.tsv", sep = '\t', header = TRUE)
onethirtyfiveRNA <- read.delim("classifications135RNA.tsv", sep = '\t', header = TRUE)
onefiftyRNA <- read.delim("classifications150RNA.tsv", sep = '\t', header = TRUE)
onesixtyfiveRNA <- read.delim("classifications165RNA.tsv", sep = '\t', header = TRUE)
twohundredRNA <- read.delim("classifications200RNA.tsv", sep = '\t', header = TRUE)

# merge classification.tsv from all depths 
merged_dataRNA <- rbind(tenRNA, onehundredRNA, onetwentyRNA, onethirtyfiveRNA, onefiftyRNA, onesixtyfiveRNA, twohundredRNA)

# Split the Sample column
merged_dataRNA <- merged_dataRNA %>%
  separate(Sample, into = c("Sample", "Depth", "extra"), sep = "_", remove = FALSE) %>%
  select(-extra) %>%  # Remove the 'extra' column
  mutate(Depth = as.numeric(gsub("m", "", Depth))) %>% # Remove 'm' and convert to numeric
  separate(Taxonomy, into = c("Root", "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";", fill = "right") %>%
  select(-Root)  # Remove the 'Root' column

NorB_RNAdata <- merged_dataRNA %>%
  filter(Marker == "NorB")

## Order ##
# Sum Abundance by Sample, Order, and Depth
RNAdata_taxa_order <- NorB_RNAdata %>%
  group_by(Sample, Order, Depth) %>%
  summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
  mutate(Order = gsub("o__", "", Order))

# Make Bubble Plot
ggplot(RNAdata_taxa_order, aes(x = Order, y = Depth)) +
  geom_point(aes(size = Abundance, color = Abundance), alpha = 0.7) +  # Map color to Abundance
  scale_y_reverse(limits = c(210, 0)) +  # Shallow at top, deep at bottom
  scale_color_gradient(low = "red", high = "blue") +  # Gradient from red (low) to blue (high)
  scale_size(range = c(1, 10)) +
  theme_light() +
  labs(title = "RNA Abundance by Depth for NorB Resolved to Order",
       x = "Order",
       y = "Depth",
       size = "RNA Abundance",
       color = "RNA Abundance") +  # Label for color legend
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("NorB_BubblePlotRNA_Order.jpg")

## GENUS ##
# Sum Abundance by Sample, Phylum, and Depth
RNAdata_taxa_genus <- NorB_RNAdata %>%
  group_by(Sample, Genus, Depth) %>%
  summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
  mutate(Genus = gsub("g__", "", Genus))

# Make Bubble Plot
ggplot(RNAdata_taxa_genus, aes(x = Genus, y = Depth)) +
  geom_point(aes(size = Abundance, color = Abundance), alpha = 0.7) +  # Map color to Abundance
  scale_y_reverse(limits = c(210, 0)) +  # Shallow at top, deep at bottom
  scale_color_gradient(low = "red", high = "blue") +  # Gradient from red (low) to blue (high)
  scale_size(range = c(1, 10)) +
  theme_light() +
  labs(title = "RNA Abundance by Depth for NorB Resolved to Genus",
       x = "Genus",
       y = "Depth",
       size = "RNA Abundance",
       color = "RNA Abundance") +  # Label for color legend
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("NorB_BubblePlotRNA_Genus.jpg")



# Figure 4: iTol plots
# generated from iTol website: https://itol.embl.de/

# Figure 5: Biogeochemistry plot
library(tidyverse)
library(ggplot2)
library(dbplyr)

# Load data
saanich_data <- read.csv("Saanich_Data.csv")

# Convert Depth from kilometers to meters
saanich_data <- saanich_data %>%
  mutate(Depth = Depth * 1000)

# 1. Summarize the data to calculate mean and standard deviation for each depth and chemical
summary_data <- saanich_data %>%
  group_by(Depth) %>%
  summarise(
    mean_WS_H2S = mean(WS_H2S, na.rm = TRUE),
    sd_WS_H2S = sd(WS_H2S, na.rm = TRUE),
    mean_WS_NO3 = mean(WS_NO3, na.rm = TRUE),
    sd_WS_NO3 = sd(WS_NO3, na.rm = TRUE),
    mean_WS_O2 = mean(WS_O2, na.rm = TRUE),
    sd_WS_O2 = sd(WS_O2, na.rm = TRUE)
  )

# 2. Pivot longer to get standard deviations and means
mean_data <- summary_data %>%
  pivot_longer(
    cols = starts_with("mean_WS_"),
    names_to = "Chemical",
    names_prefix = "mean_WS_",
    values_to = "Concentration_uM"
  )

sd_data <- summary_data %>%
  pivot_longer(
    cols = starts_with("sd_WS_"),
    names_to = "Chemical",
    names_prefix = "sd_WS_",
    values_to = "sd_Concentration_uM"
  )

# 3. Combine mean and standard deviation data
combined_data <- left_join(mean_data, sd_data, by = c("Depth", "Chemical"))

# 4. Filter out NaN values
combined_data <- combined_data %>%
  filter(!is.nan(Concentration_uM))

# 5. Order data by Chemical and Depth
combined_data <- combined_data %>%
  arrange(Chemical, Depth)

# 6. Plot all chemicals (O2, NO3, H2S)
ggplot(combined_data, aes(x = Concentration_uM, y = Depth, colour = Chemical)) +
  geom_point() +
  geom_line(aes(group = Chemical), orientation = "y") +
  facet_grid(. ~ Chemical, scales = "free_x") +
  scale_y_reverse() +
  geom_errorbar(aes(xmin = Concentration_uM - sd_Concentration_uM,
                    xmax = Concentration_uM + sd_Concentration_uM)) +
  labs(y = "Depth (m)", x = "Concentration (µM)") +
  theme_bw()

# Plot with GAM smoothing
ggplot(combined_data, aes(x = Concentration_uM, y = Depth, colour = Chemical)) +
  geom_point() +
  geom_line(aes(group = Chemical), orientation = "y") +
  facet_grid(. ~ Chemical, scales = "free_x") +
  scale_y_reverse() +
  geom_errorbar(aes(xmin = Concentration_uM - sd_Concentration_uM,
                    xmax = Concentration_uM + sd_Concentration_uM)) +
  geom_smooth(method = "gam", se = TRUE) +
  labs(y = "Depth (m)", x = "Concentration (µM)") +
  theme_bw()



# Supplementary: Retention plot
library(tidyverse)   # Includes ggplot2, dplyr, readr, etc.
library(pheatmap)
# abundance of Na 

rna_class_data_10m <- read.table(file = 'RNA_classification/s1072_10_rna_classifications.tsv', sep = '\t', header = TRUE)
rna_class_data_100m <- read.table(file = 'RNA_classification/s107_rna_100_classifications.tsv', sep = '\t', header = TRUE)
rna_class_data_120m <- read.table(file = 'RNA_classification/s1072_rna_120_classifications.tsv', sep = '\t', header = TRUE)
rna_class_data_135m <- read.table(file = 'RNA_classification/s1072_135_rna_classifications.tsv', sep = '\t', header = TRUE)
rna_class_data_150m <- read.table(file = 'RNA_classification/s1072_150_rna_classifications.tsv', sep = '\t', header = TRUE)
rna_class_data_165m <- read.table(file = 'RNA_classification/s1072_rna_165_classifications.tsv', sep = '\t', header = TRUE)
rna_class_data_200m <- read.table(file = 'RNA_classification/s1702_rna_200classifications.tsv', sep = '\t', header = TRUE)

all_rna_data <- rbind(rna_class_data_10m, 
                      rna_class_data_100m, 
                      rna_class_data_120m, 
                      rna_class_data_135m, 
                      rna_class_data_150m, 
                      rna_class_data_165m, 
                      rna_class_data_200m)


all_rna_data <- all_rna_data %>%
  separate(Sample, into = c("Sample", "Depth", "extra"), sep = "_", remove = FALSE) %>%
  select(-extra) %>%  # Remove the 'extra' column
  mutate(Depth = as.numeric(gsub("m", "", Depth))) %>% # Remove 'm' and convert to numeric
  separate(Taxonomy, into = c("Root", "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";", fill = "right") %>%
  select(-Root)  # Remove the 'Root' column

rna_norb_data <- all_rna_data %>%
  filter(Marker == "NorB")

taxonomic_levels <- c("Domain","Phylum", "Class", "Order", "Family", "Genus", "Species")

rna_na_percentages <- data.frame(Taxonomic_Level = character(), Percentage = numeric(), stringsAsFactors = FALSE)

# Loop through each taxonomic level
for (level in taxonomic_levels) {
  percentage <- rna_norb_data %>%
    summarise(
      total_abundance = sum(Abundance, na.rm = TRUE),                                  # Total Abundance across all rows
      na_level_abundance = sum(Abundance[is.na(.data[[level]])], na.rm = TRUE),        # Sum of Abundance where current taxonomic level is NA
      percentage = (na_level_abundance / total_abundance) * 100                        # Calculate percentage
    ) %>%
    pull(percentage)
  
  
  rna_na_percentages <- rbind(rna_na_percentages, data.frame(Taxonomic_Level = level, Percentage = percentage))
}

print(na_percentages)

#make rarefaction like curve 

rna_retained <- data.frame ( 
  Taxonomic_level = rna_na_percentages$Taxonomic_Level,
  Percentage = 100 - rna_na_percentages$Percentage)

rna_rarefaction_plot <- ggplot(rna_retained, aes(x = factor(Taxonomic_level, levels = unique(Taxonomic_level)),
                                                 y = Percentage)) +
  geom_line(group = 1, color = "grey") +  
  geom_point(size = 3, color = "lightblue") +  
  labs(title = "RNA Taxonomic Level Retention Curve",
       x = "Taxonomic Level",
       y = "Percentage of RNA Abundance Retained (%)") +
  theme_minimal() 
scale_y_continuous(limits = c(0, 100))

print(rna_rarefaction_plot)


library(tidyverse)   # Includes ggplot2, dplyr, readr, etc.
library(pheatmap)

dna_class_data_10m <- read.table(file = 'S1072_10m_classifications.tsv', sep = '\t', header = TRUE)
dna_class_data_100m <- read.table(file = 'S1072_100m_classifications.tsv', sep = '\t', header = TRUE)
dna_class_data_120m <- read.table(file = 'S1072_120m_classifications.tsv', sep = '\t', header = TRUE)
dna_class_data_135m <- read.table(file = 'S1072_135m_classifications.tsv', sep = '\t', header = TRUE)
dna_class_data_150m <- read.table(file = 'S1072_150m_classifications.tsv', sep = '\t', header = TRUE)
dna_class_data_165m <- read.table(file = 'S1072_165m_classifications.tsv', sep = '\t', header = TRUE)
dna_class_data_200m <- read.table(file = 'S1072_200m_classifications.tsv', sep = '\t', header = TRUE)

dna_all_class_data <- rbind(dna_class_data_10m, 
                            dna_class_data_100m, 
                            dna_class_data_120m, 
                            dna_class_data_135m, 
                            dna_class_data_150m, 
                            dna_class_data_165m, 
                            dna_class_data_200m)

dna_all_class_data <- dna_all_class_data %>%
  separate(Sample, into = c("Sample", "Depth", "extra"), sep = "_", remove = FALSE) %>%
  select(-extra) %>%  # Remove the 'extra' column
  mutate(Depth = as.numeric(gsub("m", "", Depth))) %>% # Remove 'm' and convert to numeric
  separate(Taxonomy, into = c("Root", "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";", fill = "right") %>%
  select(-Root)  # Remove the 'Root' column

dna_norb_data <- all_class_data %>%
  filter(Marker == "NorB")

taxonomic_levels <- c("Domain","Phylum", "Class", "Order", "Family", "Genus", "Species")

dna_na_percentages <- data.frame(Taxonomic_Level = character(), Percentage = numeric(), stringsAsFactors = FALSE)

# Loop through each taxonomic level
for (level in taxonomic_levels) {
  percentage <- dna_norb_data %>%
    summarise(
      total_abundance = sum(Abundance, na.rm = TRUE),                                  # Total Abundance across all rows
      na_level_abundance = sum(Abundance[is.na(.data[[level]])], na.rm = TRUE),        # Sum of Abundance where current taxonomic level is NA
      percentage = (na_level_abundance / total_abundance) * 100                        # Calculate percentage
    ) %>%
    pull(percentage)
  
  dna_na_percentages <- rbind(dna_na_percentages, data.frame(Taxonomic_Level = level, Percentage = percentage))
}

print(dna_na_percentages)

#make rarefaction like curve 

dna_retained <- data.frame ( 
  Taxonomic_level = dna_na_percentages$Taxonomic_Level,
  Percentage = 100 - dna_na_percentages$Percentage)

dna_rarefaction_plot <- ggplot(dna_retained, aes(x = factor(Taxonomic_level, levels = unique(Taxonomic_level)),
                                                 y = Percentage)) +
  geom_line(group = 1, color = "grey") +  
  geom_point(size = 3, color = "lightgreen") +  
  labs(title = "DNA Taxonomic Level Retention Curve",
       x = "Taxonomic Level",
       y = "Percentage of DNA Abundance Retained (%)") +
  theme_minimal() 
scale_y_continuous(limits = c(0, 100))

print(dna_rarefaction_plot)