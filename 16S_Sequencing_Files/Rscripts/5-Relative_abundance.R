## Load library
library(dplyr)
library(vegan)

## Import the rarefied file with taxonomy info
bac.asv.tax <- read.table("Data/SBAR_16S_Box_Rarefied_09262022.txt", 
                          sep = "\t", header = T, row.names = 1)
dim(bac.asv.tax) #9250    47

# Extract only the sequencing information into asv 
bac.asv <- t(bac.asv.tax[,-c(41:47)])
dim(bac.asv) #40 9250

## Taxonomy 
taxa <- t(bac.asv.tax[,c(41:47)])
dim(taxa) #7 9250

## Import mapping file
metadata <- read.table("Data/SBAR_Incubation_MappingFile_05142022.txt", 
                       sep = "\t", header = T, row.names = 1, check.names = F) %>%
            filter(Sample_type == "Box")
dim(metadata) #40 14


## Match asv with samples with metadata
all(rownames(bac.asv)==rownames(metadata))

#### Calculate relative abundance ####
apply(bac.asv, 1, sum) #all samples should have the same number if rarefied
bac.rel <- decostand(bac.asv, 
                     method = "total")
apply(bac.rel, 1, sum) #1 for all samples

## Transpose the relative abundance dataframe before adding taxonomy info
# bac.rel.t <- as.data.frame(t(bac.rel))
# all(rownames(bac.rel.t) == rownames(taxa)) #should be TRUE

box.rel <- cbind(as.data.frame(t(bac.rel)), t(taxa))
dim(box.rel) #9250 47

#### Add the taxonomy column ####
## sep = "|" is used here
# bac.rel.t$taxonomy <- paste(taxa$Kingdom, 
#                           taxa$Phylum, 
#                           taxa$Class,
#                           taxa$Order, 
#                           taxa$Family, 
#                           taxa$Genus,
#                           taxa$Species, 
#                           sep = "|")

#### Append a column with ASV numbers ####
## First add column 'asv' with sequential numbers
# bac.oak.t <- cbind(bac.oak.t, 
#                    "asv" = 1:nrow(bac.oak.t)) 
# 
# ## Create a column 'ID' with text 'bac' appended with numbers from asv column
# bac.oak.t$ID <- paste('bac', 
#                       bac.oak.t$asv, 
#                       sep = "_") 



write.table(box.rel, 
            file = "Data/SBAR_16S_Abundance_Box_09262022.txt", 
            sep = "\t", 
            quote = F, 
            row.names = T, 
            col.names = NA)