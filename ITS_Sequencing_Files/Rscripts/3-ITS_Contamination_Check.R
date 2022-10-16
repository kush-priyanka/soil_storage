setwd("C:/Users/priya/Box/AZ_UA_Jan2022/2022_Ongoing_Projects/SBAR2022/Soil_Storage_Sequencing//ITS_Sequencing_Files/")

#Import asv table with sequence and taxonomy information (w blank info)
#File generated from dada2: asv_table_Taxa_ITS_sbar_may2022_0922022.txt 
asv.tax <- read.table("paired_read_post_dada2_files/asv_table_Taxa_ITS_sbar_may2022_09222022.txt",
                    sep = "\t", header = T, row.names = 1)
dim(asv.tax) #2272 82

#Subset asv reads & taxonomy into two variables
asv <- t(asv.tax[,-c(76:82)])
dim(asv) #75 2272
asv <- asv[c(1:74),] #74 6393

tax <- asv.tax[,c(76:82)]
dim(tax) #2272 7

#ensure both variables are data frames
asv <- as.data.frame(asv)
tax <- as.data.frame(tax)

#Read metadata
metadata <- read.table("Data/SBAR_Incubation_MappingFile_08232022.txt", 
                     sep = "\t", header = T, row.names = 1, check.names = F)
dim(metadata) #89 13

#Remove samples not related to the project by matching samples IDs from mapping file
metadata <- metadata[rownames(asv),] 
dim(metadata) #74 13

#Check if all sample names in ASV table match mapping file, should be TRUE for all
all(rownames(metadata) == rownames(asv))

#As you have removed samples,you may have empty ASV columns 
asv <- subset(asv, select = colSums(asv)!= 0)
dim(asv) #74 2236

#Remove any samples from taxonomy table that are not in ASV table 
tax <- as.data.frame(tax[colnames(asv),])
dim(tax) #2236    7
##Check if all sequences in ASV table match taxonomy sequences
all(rownames(tax) == colnames(asv))

#Only keep Fungal sequences
#Remove contaminants:any other Eukaryota
asv <- asv[,grep("k__Fungi", tax$Kingdom)]
tax <- tax[grep("k__Fungi", tax$Kingdom),]


##Merge and taxonomic information including blanks
asv.tax <- cbind(t(asv), tax)
dim(asv.tax) #2236 81

#Save the file
write.table(asv.tax, 
            file = "Data/SBAR_ITS_Taxa_wBlanks_Post-Process_09222022.txt", 
            sep = "\t", quote = F, row.names = T, col.names = NA)

#Check your blank DNA extracts and remove potential lab contamination
#Assign a variable to each of your blanks
ext01.control <- asv[which(rownames(asv) == "S_012"),]
ext03.control <- asv[which(rownames(asv) == "S_031"),]

zymo <- asv[which(rownames(asv) == "S_089"),]


#Check the number of ASVs in each blank
sum(ext01.control)  #173
sum(ext03.control) #1775
sum(zymo)  #24737


#Now check, which ASV is greater than 1% of the total ASVs in that blank
which(ext01.control > 0.000001* sum(ext01.control)) #4
which(ext03.control > 0.000000001* sum(ext03.control)) #5

which(zymo > 0.00000000001* sum(zymo)) #6

#remove asv
tax.wocon <- tax[-c(22,331,333,635,913,1048,184,565,911),]
dim(tax.wocon) #2227 7

asv.wocon <- as.data.frame(asv[,rownames(tax.wocon)])
dim(asv.wocon) #74 2227

#Remove any empty ASVs
asv.wocon <- subset(asv.wocon, select = colSums(asv.wocon)!= 0)
dim(asv.wocon) #74 2227

all(rownames(tax.wocon) == colnames(asv.wocon))

##Merge and taxonomic information
final.tax <- cbind(t(asv.wocon), tax.wocon)
dim(final.tax) #2227 81

#Save the file
write.table(final.tax, 
            file = "Data/SBAR_ITS_ASV_table_Final_woCon_09222022.txt", 
            sep = "\t", quote = F, row.names = T, col.names = NA)

