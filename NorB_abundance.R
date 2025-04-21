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
