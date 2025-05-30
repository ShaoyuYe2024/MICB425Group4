Analysis_SY
================
Shaoyu Ye
2025-03-28

## R Markdown

``` r
# Load the libraries
library(vegan)
```

    ## 载入需要的程序包：permute

    ## 载入需要的程序包：lattice

``` r
library(tidyverse)   # Includes ggplot2, dplyr, readr, etc.
```

    ## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ## ✔ dplyr     1.1.4     ✔ readr     2.1.5
    ## ✔ forcats   1.0.0     ✔ stringr   1.5.1
    ## ✔ ggplot2   3.5.1     ✔ tibble    3.2.1
    ## ✔ lubridate 1.9.4     ✔ tidyr     1.3.1
    ## ✔ purrr     1.0.4

    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()
    ## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

``` r
library(pheatmap)

# Load the data from the CSV file
alpha_data <- read_csv("../alpha_diversity/SI_TS_NorB_alpha_diversiy.csv")  
```

    ## Rows: 7 Columns: 6
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: ","
    ## chr (1): placerun
    ## dbl (5): phylo_entropy, quadratic, unrooted_pd, rooted_pd, bwpd
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
beta_data <- read.csv("../beta_diversity/SI_TS_NorB_beta_diversity.csv")

# Display the first few rows of the dataset
head(alpha_data)
```

    ## # A tibble: 6 × 6
    ##   placerun                   phylo_entropy quadratic unrooted_pd rooted_pd  bwpd
    ##   <chr>                              <dbl>     <dbl>       <dbl>     <dbl> <dbl>
    ## 1 SI072_100m_NorB_complete_…          3.81      1.84        25.2      26.3  4.56
    ## 2 SI072_10m_NorB_complete_p…          2.66      1.25        16.1      18.4  2.91
    ## 3 SI072_120m_NorB_complete_…          4.28      2.03        31.3      32.4  5.11
    ## 4 SI072_135m_NorB_complete_…          4.54      2.20        29.0      30.1  5.72
    ## 5 SI072_150m_NorB_complete_…          4.75      2.34        31.3      32.4  6.51
    ## 6 SI072_165m_NorB_complete_…          4.78      2.45        26.2      27.3  7.34

``` r
head(beta_data)
```

    ##                           sample_1                         sample_2      Z_1
    ## 1 SI072_100m_NorB_complete_profile  SI072_10m_NorB_complete_profile 2.027050
    ## 2 SI072_100m_NorB_complete_profile SI072_120m_NorB_complete_profile 0.643971
    ## 3 SI072_100m_NorB_complete_profile SI072_135m_NorB_complete_profile 1.038210
    ## 4 SI072_100m_NorB_complete_profile SI072_150m_NorB_complete_profile 1.730160
    ## 5 SI072_100m_NorB_complete_profile SI072_165m_NorB_complete_profile 2.386000
    ## 6 SI072_100m_NorB_complete_profile SI072_200m_NorB_complete_profile 2.074650

``` r
## Saanich TreeSAPP NorB
# Split the placerun column
alpha_data <- alpha_data %>%
  separate(placerun, into = c("sample", "depth", "gene", "extra"), sep = "_", remove = FALSE) %>%
  select(-extra) %>%  # Remove the 'extra' column
  mutate(depth = as.numeric(gsub("m", "", depth)))  # Remove 'm' and convert to numeric
```

    ## Warning: Expected 4 pieces. Additional pieces discarded in 7 rows [1, 2, 3, 4,
    ## 5, 6, 7].

``` r
# Display the first few rows of the dataset
head(alpha_data)
```

    ## # A tibble: 6 × 9
    ##   placerun      sample depth gene  phylo_entropy quadratic unrooted_pd rooted_pd
    ##   <chr>         <chr>  <dbl> <chr>         <dbl>     <dbl>       <dbl>     <dbl>
    ## 1 SI072_100m_N… SI072    100 NorB           3.81      1.84        25.2      26.3
    ## 2 SI072_10m_No… SI072     10 NorB           2.66      1.25        16.1      18.4
    ## 3 SI072_120m_N… SI072    120 NorB           4.28      2.03        31.3      32.4
    ## 4 SI072_135m_N… SI072    135 NorB           4.54      2.20        29.0      30.1
    ## 5 SI072_150m_N… SI072    150 NorB           4.75      2.34        31.3      32.4
    ## 6 SI072_165m_N… SI072    165 NorB           4.78      2.45        26.2      27.3
    ## # ℹ 1 more variable: bwpd <dbl>

# Alpha Diversity

This plot shows all the alpha diversity metrics for DrsAB across the 7
depths of cruise 72. What does this tell us based on the definition of
these metrics? (see Project_TreeSAPP_DsrAB_diversity for definitions)

``` r
# Select relevant columns
alpha_data_wide <- alpha_data %>%
  select(depth, phylo_entropy, quadratic, unrooted_pd, rooted_pd, bwpd)

# Pivot data to a longer format for faceting
alpha_data_longer <- alpha_data_wide %>%
  pivot_longer(cols = c("phylo_entropy", "quadratic", "unrooted_pd", "rooted_pd", "bwpd"),
               names_to = "metric",
               values_to = "value")

# Define your desired order of metrics
desired_order <- c("unrooted_pd", "rooted_pd", "phylo_entropy", "quadratic", "bwpd")

# Convert 'metric' column to a factor with specified levels
alpha_data_longer$metric <- factor(alpha_data_longer$metric, levels = desired_order)

ggplot(alpha_data_longer, aes(y = depth, x = value)) +
  geom_point(aes(colour = metric, shape = metric), size = 4) +
  scale_y_reverse() +
  labs(title = "Alpha Diversity Metric against Depth",
       x = "value") +
  theme_minimal() +
  facet_grid(. ~ metric, scales = "free_x")
```

![](Analysis_Code_SY_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

\#Beta Diversity Now we can look at samples compared to eachother using
KR distance. But first the data needs some cleaning.

``` r
# Split 'sample_1' and 'sample_2' into their components (optional but useful for clarity)
beta_data <- beta_data %>%
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
all_samples <- sort(unique(c(beta_data$sample_1, beta_data$sample_2)))

# Create an empty matrix filled with NA
beta_matrix_complete <- matrix(NA,
                               nrow = length(all_samples),
                               ncol = length(all_samples),
                               dimnames = list(all_samples, all_samples))

# Fill in the matrix symmetrically
for (i in seq_len(nrow(beta_data))) {
  row_name <- as.character(beta_data$sample_1[i])
  col_name <- as.character(beta_data$sample_2[i])
  value <- beta_data$Z_1[i]
  
  # Assign value to both [row, col] and [col, row] to make it symmetrical
  beta_matrix_complete[row_name, col_name] <- value
  beta_matrix_complete[col_name, row_name] <- value
}

# Replace NA with 0 if desired (or leave as NA for missing values)
beta_matrix_complete[is.na(beta_matrix_complete)] <- 0

# Convert to a data frame for compatibility with pheatmap
beta_matrix_clean <- as.data.frame(beta_matrix_complete)

# View the completed matrix
print(beta_matrix_clean)
```

    ##          10      100      120      135      150      165      200
    ## 10  0.00000 2.027050 2.308510 2.649310 3.159560 3.640060 3.415370
    ## 100 2.02705 0.000000 0.643971 1.038210 1.730160 2.386000 2.074650
    ## 120 2.30851 0.643971 0.000000 0.599217 1.278250 1.937730 1.593780
    ## 135 2.64931 1.038210 0.599217 0.000000 0.874400 1.487100 1.171910
    ## 150 3.15956 1.730160 1.278250 0.874400 0.000000 0.971246 0.764753
    ## 165 3.64006 2.386000 1.937730 1.487100 0.971246 0.000000 0.687874
    ## 200 3.41537 2.074650 1.593780 1.171910 0.764753 0.687874 0.000000

Now the plotting:

``` r
# Sort rows by their names or a specific column
sorted_matrix <- beta_matrix_clean[order(rownames(beta_matrix_clean)), , drop = FALSE]

# Cluster columns
col_clust <- hclust(dist(t(sorted_matrix)))  # Create column dendrogram
# Flip dendrogram branches if desired
col_clust$order <- order(rownames(beta_matrix_clean))

# Plot the heatmap
pheatmap(as.matrix(sorted_matrix),
         cluster_rows = FALSE,
         cluster_cols = col_clust,
         scale = "none",
         color = colorRampPalette(c("black", "gray", "white"))(100),
         border_color = NA,
         main = "KR Distance Heatmap with Dendrogram")
```

![](Analysis_Code_SY_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

# Bubble Plots

We have a good picture of the overall diversity for DsrAB across the
samples. Let’s explore the taxonomy of a sample and associate that with
the abundance of DrsAB orthologs.

In order to look at all the samples together we need to contatenate the
TreeSAPP classifications.tsv from each depth. Be sure to carefully
replace the below paths for the files with the paths to the files you’d
like to combine.

Load Data

``` r
# List of classification files
file_list <- c("../classifications/classifications_10m.tsv",
               "../classifications/classifications_100m.tsv",
               "../classifications/classifications_120m.tsv",
               "../classifications/classifications_135m.tsv",
               "../classifications/classifications_150m.tsv",
               "../classifications/classifications_165m.tsv",
               "../classifications/classifications_200m.tsv"
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
```

    ## # A tibble: 0 × 19
    ## # ℹ 19 variables: Sample <chr>, Depth <dbl>, Query <chr>, Marker <chr>,
    ## #   Start_pos <chr>, End_pos <chr>, Domain <chr>, Phylum <chr>, Class <chr>,
    ## #   Order <chr>, Family <chr>, Genus <chr>, Species <chr>, Abundance <chr>,
    ## #   iNode <chr>, E-value <chr>, LWR <chr>, EvoDist <chr>, Distances <chr>

``` r
# Display the first few rows of the dataset
head(sub_data)
```

    ## # A tibble: 6 × 19
    ##   Sample Depth Query   Marker Start_pos End_pos Domain Phylum Class Order Family
    ##   <chr>  <dbl> <chr>   <chr>  <chr>     <chr>   <chr>  <chr>  <chr> <chr> <chr> 
    ## 1 SI072     10 k147_1… NorB   1         211     " d__… " p__… " c_… " o_… " f__…
    ## 2 SI072     10 k147_4… NorB   1         218     " d__… " p__… " c_… " o_… " f__…
    ## 3 SI072     10 k147_3… NorB   1         111     " d__… " p__… " c_… " o_… " f__…
    ## 4 SI072     10 k147_4… NorB   1         222     " d__… " p__… " c_… " o_… " f__…
    ## 5 SI072     10 k147_3… NorB   3         477     " d__… " p__… " c_… " o_… " f__…
    ## 6 SI072     10 k147_1… NorB   1         44      " d__… " p__… " c_… " o_… " f__…
    ## # ℹ 8 more variables: Genus <chr>, Species <chr>, Abundance <chr>, iNode <chr>,
    ## #   `E-value` <chr>, LWR <chr>, EvoDist <chr>, Distances <chr>

# Bubbleplot for NorB

``` r
# Sum Abundance by Sample, Phylum, and Depth
data_taxa <- sub_data %>%
  mutate(Abundance = as.numeric(Abundance)) %>%  
  group_by(Sample, Class, Depth) %>%
  summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
  mutate(Class = gsub("c__", "", Class))

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
```

![](Analysis_Code_SY_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->
