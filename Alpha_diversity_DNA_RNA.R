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
