#Install dada2 package
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("dada2")
library(dada2)
packageVersion("dada2") #confirm that DADA2 version is 1.16 or later

#Set working directory to the dataset folder, change to your directory
path <- "/home/u22/pkushwaha/SBAR_May2022/demultiplexed"

#### Define variable path ###
list.files(path)
fnFs <- sort(list.files(path, pattern = "_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_R2_001.fastq", full.names = TRUE))

#### Read in the samples files ####
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 7)
sample.names <- sapply(strsplit(basename(sample.names), "\\.f"), `[`, 1)
sample.names <- paste("S",sample.names,sep = "_")

head(sample.names)

#### Inspect Read Quality Profiles ####
## Forward reads
plotQualityProfile(fnFs[14:24])

## Reverse reads
plotQualityProfile(fnRs[14:24])


#### Filter and trim ####
filt_path <- file.path("/home/u22/pkushwaha/SBAR_May2022", 'filtered_16S')
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, trimLeft = 10, truncLen = c(140,140),
                     maxN = 0, maxEE = c(2,2), truncQ = 2, rm.phix = TRUE,
                     compress = TRUE, multithread = TRUE)

head(out)
View(out)


#### Learn the Error Rates ####
errF <- learnErrors(filtFs[4:10], multithread=TRUE)
errR <- learnErrors(filtRs[4:10], multithread=TRUE)

## Plot error rates
plotErrors(errF, nominalQ = TRUE)
plotErrors(errR, nominalQ = TRUE)


#### DEREPLICATION ####
derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)
#names(derepFs) <- sample.names
#names(derepRs) <- sample.names


#### Sample Inference ####
dadaFs <- dada(derepFs, err = errF, multithread = TRUE)
dadaRs <- dada(derepRs, err = errR, multithread = TRUE)

#describe data from the forward reads of the first sample
dadaFs[[1]]

#### Merge paired reads ####
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose = TRUE)
head(mergers[1])

#### Construct sequence table ####
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))

#### Remove chimeras ####
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", 
                                    multithread = TRUE, verbose = TRUE)
dim(seqtab.nochim)

## percentage of total sequence reads that are NOT chimeras 
sum(seqtab.nochim)/sum(seqtab)

## Export ASV file with reads
write.table(seqtab.nochim, "All_Seqs_16S_reads_sbar_may2022_05162022.txt", 
            sep = "\t", quote = F, row.names = T, col.names = NA)

saveRDS(seqtab.nochim, "All_Seqs_16S_reads_sbar_may2022_05162022.rds")

#### Track reads through the pipeline ####
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
rownames(track) <- sample.names
head(track)
View(track)

## Export a table with reads that have made through the pipeline
write.table(track, "All_Seqs_Filtered_ReadsCounts_16S_sbar_may2022_05162022.txt", 
            sep = "\t", quote = F, row.names = T, col.names = NA)

#### ASSIGN TAXONOMY ####
taxa.genus <- assignTaxonomy(seqtab.nochim, 
                             "/home/u22/pkushwaha/SBAR_May2022/SILVA/silva_nr99_v138.1_train_set.fa.gz", 
                             multithread = TRUE)

## Export table with genus taxonomy only
write.table(taxa.genus, "All_Seqs_silva_16S_genus_sbar_may2022_05162022.txt", 
            sep = "\t", quote = F, row.names = T, col.names = NA)

asv.tax <- cbind(t(seqtab.nochim), taxa.genus)
dim(asv.tax)

## Save the file with asv reads and taxonomy info
write.table(asv.tax, file="asv_table_genus_sbar_may2022_05162022.txt", 
            sep = "\t", quote = F, row.names = T, col.names = NA)


#### SPECIES ASSIGNMENT ####
taxa.species <- addSpecies(taxa.genus, 
                           "/home/u22/pkushwaha/SBAR_May2022/SILVA/silva_species_assignment_v138.1.fa.gz")

## View the assigned taxonomic table
taxa.print <- taxa.species
rownames(taxa.print) <- NULL
head(taxa.print)
View(taxa.print)

## Export table with taxonomy only
write.table(taxa.species, "All_Seqs_silva_16S_species_sbar_may2022_05162022.txt", 
            sep = "\t", quote = F, row.names = T, col.names = NA)

## Merge reads and taxonomic information
bac.asv.tax <- cbind(t(seqtab.nochim), taxa.species)
dim(bac.asv.tax)

## Save the file with asv reads and taxonomy info
write.table(bac.asv.tax, file="asv_table_species_sbar_may2022_05162022.txt", 
            sep = "\t", quote = F, row.names = T, col.names = NA)

