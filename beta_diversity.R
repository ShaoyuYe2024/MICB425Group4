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





