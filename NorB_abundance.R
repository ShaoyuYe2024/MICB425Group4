


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