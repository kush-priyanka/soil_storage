install.packages("seqinr")
library(seqinr)

setwd('D:/Amplicon_Sequencing_Datasets/SBAR_SoilStorage_May2022/Fungi_ITS')

barcode <- read.table('SBAR_ITS_Barcode.txt', sep='\t', header=T)
barcode$Barcode <- as.character(barcode$Barcode)

barcode$Barcode_revcomp <- sapply(barcode$Barcode, 
                                  function(x){toupper(c2s(rev(comp(s2c(x)))))})

barcode_revcomp <- data.frame(barcode$Barcode_revcomp, barcode$SampleID)
colnames(barcode_revcomp) <- c("Barcode", "SampleID")

write.table(barcode_revcomp, "SBAR_ITS_barcoderevcomp.txt", 
            sep = '\t', quote = F, row.names = F)