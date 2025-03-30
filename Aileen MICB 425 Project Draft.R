library(tidyverse)
library(pheatmap)
library(dplyr)
NorB_alpha <- read.csv("SI_TS_NorB_alpha_diversiy.csv")
NorB_beta <- read.csv("SI_TS_NorB_beta_diversity.csv")

# placerun is ID for the file analyzed

# phylo_entropy is phylogenetic entropy
# where high entropy indicates placements are spread widely across the tree
# low entropy means placements are concentrated in specific areas

# quadratic entropy is similar to phylogenetic entropy but calculated
# using a quadratic (squared) measure of differences between placement distributions

# unrooted_pd is equivalent to Faith's PD

# rooted_pd is rooted phylogenetic diversity, similar to unrooted PD but 
# considers higher level branch lengths

# bwpd is generalized PD, incorporates abundance information for placements



# Split the placerun column
NorB_alpha <- NorB_alpha %>%
  separate(placerun, into = c("sample", "depth", "gene", "extra"), sep = "_", remove = FALSE) %>%
  select(-extra) %>%  # Remove the 'extra' column
  mutate(depth = as.numeric(gsub("m", "", depth)))  # Remove 'm' and convert to numeric
# Display the first few rows of the dataset
head(NorB_alpha)

# Select relevant columns
NorB_alpha_wide <- NorB_alpha %>%
  select(depth, phylo_entropy, quadratic, unrooted_pd, rooted_pd, bwpd)

# Pivot data to a longer format for faceting
NorB_alpha_longer <- NorB_alpha_wide %>%
  pivot_longer(cols = c("phylo_entropy", "quadratic", "unrooted_pd", "rooted_pd", "bwpd"),
               names_to = "metric",
               values_to = "value")

# Define your desired order of metrics
desired_order <- c("unrooted_pd", "rooted_pd", "phylo_entropy", "quadratic", "bwpd")

# Convert 'metric' column to a factor with specified levels
NorB_alpha_longer$metric <- factor(NorB_alpha_longer$metric, levels = desired_order)

ggplot(NorB_alpha_longer, aes(y = depth, x = value)) +
  geom_point(aes(colour = metric, shape = metric), size = 4) +
  scale_y_reverse() +
  labs(title = "Alpha Diversity Metric against Depth",
       x = "value") +
  theme_minimal() +
  facet_grid(. ~ metric, scales = "free_x")

# List of classification files
file_list <- c("SI072_100m_classifications.tsv",
               "SI072_120m_classifications.tsv",
               "SI072_150m_classifications.tsv",
               "SI072_200m_classifications.tsv",
               "SI072_10m_classifications.tsv",
               "SI072_135m_classifications.tsv",
               "SI072_165m_classifications.tsv"
)


# Merge all data frames by row names (assuming all files have a common columns)
class_data <- file_list %>%
  lapply(read_tsv, col_types = cols(.default = "c")) %>%
  bind_rows()
# Split the Sample column
class_data <- class_data %>%
  separate(Sample, into = c("Sample", "Depth", "extra"), sep = "_", remove = FALSE) %>%
  select(-extra) %>%  # Remove the 'extra' column
  mutate(Depth = as.numeric(gsub("m", "", Depth))) %>% # Remove 'm' and convert to numeric
  separate(Taxonomy, into = c("Root", "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";", fill = "right") %>%
  select(-Root)  # Remove the 'Root' column
sub_data <- class_data %>%
  filter(Marker == "NorB")
sub_data %>% filter(is.na(as.numeric(Abundance)))

# Sum Abundance by Sample, Phylum, and Depth
data_taxa <- sub_data %>%
  mutate(Abundance = as.numeric(Abundance)) %>%  
  group_by(Sample, Class, Depth) %>%
  summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
  mutate(Class = gsub("c__", "", Class))

# Bubble plot
ggplot(data_taxa, aes(x = Class, y = Depth)) +
  geom_point(aes(size = Abundance), color = "black", alpha = 0.7) +  # Fixed black color
  scale_y_reverse(limits = c(210, 0)) +  # Shallow at top, deep at bottom
  scale_size(range = c(1, 10)) +
  theme_minimal() +
  labs(title = "Bubble Plot of Class Abundance by Depth for NorB",
       x = "Class",
       y = "Depth",
       size = "Abundance") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Bubble plot with bisexual lighting
ggplot(data_taxa, aes(x = Class, y = Depth)) +
  geom_point(aes(size = Abundance, color = Abundance), alpha = 0.7) +  # Map color to Abundance
  scale_y_reverse(limits = c(210, 0)) +  # Shallow at top, deep at bottom
  scale_color_gradient(low = "red", high = "blue") +  # Gradient from red (low) to blue (high)
  scale_size(range = c(1, 10)) +
  theme_light() +
  labs(title = "Bubble Plot of Class Abundance by Depth for NorB",
       x = "Class",
       y = "Depth",
       size = "Abundance",
       color = "Abundance") +  # Label for color legend
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

