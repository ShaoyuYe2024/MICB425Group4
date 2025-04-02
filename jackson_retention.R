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
