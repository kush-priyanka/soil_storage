## Load library
library(dplyr)
library(vegan)

## Import the rarefied file with taxonomy info
asv.tax <- read.table("Data/SBAR_ITS_Box_Rarefied_09232022.txt", 
                          sep = "\t", header = T, row.names = 1)
dim(asv.tax) #1365 45

# Extract only the sequencing information into asv 
asv <- t(asv.tax[,-c(39:45)])
dim(asv) #38  1365

## Taxonomy 
taxa <- t(asv.tax[,c(39:45)])
dim(taxa) #7 1365

## Import mapping file
metadata <- read.table("Data/SBAR_Incubation_MappingFile_08232022.txt", 
                       sep = "\t", header = T, row.names = 1, check.names = F) %>%
            filter(Sample_type == "Box")
dim(metadata) #40 13
metadata <- metadata[-which(rownames(metadata)=="S_035"),]
md <- metadata[-which(rownames(metadata)=="S_045"),]

## Match asv with samples with metadata
all(rownames(asv)==rownames(md))

asv <- subset(asv, select = colSums(asv)!= 0)
dim(asv) #38 1365

## Match taxa with asv
tax <- taxa[,colnames(asv)]
dim(tax) # 7 1365

#### Calculate relative abundance ####
apply(asv, 1, sum) #all samples should have the same number if rarefied
fun.rel <- decostand(asv, 
                     method = "total")
apply(fun.rel, 1, sum) #1 for all samples

## Transpose the relative abundance dataframe before adding taxonomy info
box.rel <- cbind(as.data.frame(t(fun.rel)), t(tax))
dim(box.rel) #1365 45

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


write.table(box.rel, 
            file = "Data/SBAR_ITS_Abundance_Rar_Box_09232022.txt", 
            sep = "\t", 
            quote = F, 
            row.names = T, 
            col.names = NA)

#### Rel abundance using normalized file###
## Import the rarefied file with taxonomy info
asv.tax.norm <- read.table("Data/SBAR_ITS_Box_Normalized_09232022.txt", 
                      sep = "\t", header = T, row.names = 1)
dim(asv.tax.norm) #1444 46

# Extract only the sequencing information into asv 
nasv <- t(asv.tax.norm[,-c(40:46)])
dim(nasv) #39 1444

## Taxonomy 
ntaxa <- t(asv.tax.norm[,c(40:46)])
dim(ntaxa) #7 1444

## Match asv with samples with metadata
all(rownames(nasv)==rownames(metadata))

#### Calculate relative abundance ####
apply(nasv, 1, sum) #all samples should have the same number if rarefied
nfun.rel <- decostand(nasv, 
                     method = "total")
apply(nfun.rel, 1, sum) #1 for all samples

## Transpose the relative abundance dataframe before adding taxonomy info
nbox.rel <- cbind(as.data.frame(t(nfun.rel)), t(ntaxa))
dim(nbox.rel) #1444 46

write.table(nbox.rel, 
            file = "Data/SBAR_ITS_Abundance_Norm_Box_09232022.txt", 
            sep = "\t", 
            quote = F, 
            row.names = T, 
            col.names = NA)