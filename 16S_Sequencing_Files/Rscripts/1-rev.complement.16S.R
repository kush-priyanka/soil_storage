install.packages("seqinr")
library(seqinr)

setwd('D:/Amplicon_Sequencing_Datasets/SBAR_SoilStorage_May2022/Bac_arc__16S')

barcode <- read.table('SBAR_16S_Barcode.txt', sep='\t', header=T)
barcode$Barcode <- as.character(barcode$Barcode)

barcode$Barcode_revcomp <- sapply(barcode$Barcode, 
                                  function(x){toupper(c2s(rev(comp(s2c(x)))))})

barcode_revcomp <- data.frame(barcode$Barcode_revcomp, barcode$SampleID)
colnames(barcode_revcomp) <- c("Barcode", "SampleID")

write.table(barcode_revcomp, "SBAR_16S_barcoderevcomp.txt", 
            sep = '\t', quote = F, row.names = F)