

#Import asv table with sequence and taxonomy information (w blank info)
#File generated from dada2: asv_table_Taxa_ITS_oak_11.18.2021.txt 
asv.tax <- read.table("Sequencing_Files/Data/asv_table_species_sbar_may2022_05162022.txt",
                    sep = "\t", header = T, row.names = 1)
dim(asv.tax) #18176 97

#Subset asv reads & taxonomy into two variables
asv <- t(asv.tax[,-c(91:97)])
dim(asv) #90 18176

tax <- asv.tax[,c(91:97)]
dim(tax) #18176 7

#ensure both variables are data frames
asv <- as.data.frame(asv)
tax <- as.data.frame(tax)

#Read metadata
metadata <- read.table("Sequencing_Files/Data/SBAR_Incubation_MappingFile_05142022.txt", 
                     sep = "\t", header = T, row.names = 1, check.names = F)
dim(metadata) #89 12

#Remove samples not related to the project by matching samples IDs from mapping file
asv <- asv[rownames(metadata),] 
dim(asv) #89 18176

#Check if all sample names in ASV table match mapping file, should be TRUE for all
all(rownames(metadata) == rownames(asv))

#As you have removed samples,you may have empty ASV columns 
asv <- subset(asv, select = colSums(asv)!= 0)
dim(asv) #89 17900

#Remove any samples from taxonomy table that are not in ASV table 
tax <- as.data.frame(tax[colnames(asv),])
dim(tax) #17900     7
##Check if all sequences in ASV table match taxonomy sequences
all(rownames(tax) == colnames(asv))

#Only keep Bacterial sequences
#Remove contaminants:chloroplasts, mitochondria and Eukaryota
grep("Chloroplast", tax$Order)
grep("Mitochondria", tax$Family)
grep("Eukaryota", tax$Kingdom)

#Remove chloroplasts from ASV and taxonomy table
asv <- asv[,-grep("Chloroplast", tax$Order)]
tax <- tax[-grep("Chloroplast", tax$Order),]
dim(asv) #89 17886
dim(tax) #17886 7

#Remove Mitochondria from ASV and taxonomy table
asv <- asv[,-grep("Mitochondria", tax$Family)]
tax <- tax[-grep("Mitochondria", tax$Family),]
dim(asv) #89 17298
dim(tax) #17298 7

#Remove Eukaryota from ASV and taxonomy table if there is any other kingdom present
#this step was skipped here because no Eukaryota was found
#asv<-asv[,-grep("Eukaryota", tax$Kingdom)]
#tax<-tax[-grep("Eukaryota", tax$Kingdom),]
#dim(asv)
#dim(tax)

##Merge and taxonomic information including blanks
asv.tax <- cbind(t(asv), tax)
dim(asv.tax) #17298    96

#Save the file
write.table(asv.tax, 
            file = "Sequencing_Files/Data/SBAR_16S_Taxa_wBlanks_Post-Process_05182022.txt", 
            sep = "\t", quote = F, row.names = T, col.names = NA)

#Check your blank DNA extracts and remove potential lab contamination
#Assign a variable to each of your blanks
ext01.control <- asv[which(rownames(asv)=="S_012"),]
ext02.control <- asv[which(rownames(asv)=="S_021"),]
ext03.control <- asv[which(rownames(asv)=="S_031"),]
ext04.control <- asv[which(rownames(asv)=="S_039"),]
ext05.control <- asv[which(rownames(asv)=="S_053"),]
ext06.control <- asv[which(rownames(asv)=="S_061"),]
ext09.control <- asv[which(rownames(asv)=="S_082"),]

seq.blank <- asv[which(rownames(asv)=="S_088"),]
zymo <- asv[which(rownames(asv)=="S_089"),]


#Check the number of ASVs in each blank
sum(ext01.control)  #158
sum(ext02.control)  #87
sum(ext03.control) #13
sum(ext04.control) #25
sum(ext05.control)  #14
sum(ext06.control)  #804
sum(ext09.control)  #1090
sum(seq.blank)  #40
sum(zymo)  #87557


#Now check, which ASV is greater than 1% of the total ASVs in that blank
which(ext01.control > 0.01* sum(ext01.control)) #3
which(ext02.control > 0.01* sum(ext02.control)) #3
which(ext03.control > 0.01* sum(ext03.control)) #2
which(ext04.control > 0.01* sum(ext04.control)) #2
which(ext05.control > 0.01* sum(ext05.control)) #2
which(ext06.control > 0.00001* sum(ext06.control)) #27
which(ext09.control > 0.00001* sum(ext09.control)) #17

which(seq.blank > 0.01* sum(seq.blank)) #3
which(zymo > 0.01* sum(zymo)) #9

#remove asv
tax.wocon <- tax[-c(69, 128, 137,  156,  170 , 234,  271, 400,1577,7163,11487, 5735, 16612,6726,12095,156,5253,5735, 6503, 
      7994, 9259,13434, 14285, 15827, 15828, 15829, 15830,16959, 
      16960, 16961, 16962, 16963, 3518, 3560, 4989, 5147,5674, 6128, 
      6726, 7696, 8450, 9092,9778, 12925, 13590, 8303, 13631),]
dim(tax.wocon) #17254 7

asv.wocon <- as.data.frame(asv[,rownames(tax.wocon)])
dim(asv.wocon) #89 17254

# # tax[c(75,7163, 11487),6]


#Remove any empty ASVs
asv.wocon <- subset(asv.wocon, select = colSums(asv.wocon)!= 0)
dim(asv.wocon) #89 17254

all(rownames(tax.wocon) == colnames(asv.wocon))

##Merge and taxonomic information
final.tax <- cbind(t(asv.wocon), tax.wocon)
dim(final.tax) #17254 96

#Save the file
write.table(final.tax, 
            file = "Sequencing_Files/Data/SBAR_16S_ASV_table_Final_05182022.txt", 
            sep = "\t", quote = F, row.names = T, col.names = NA)

