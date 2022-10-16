#Install libraries
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("dada2", version = "3.14")
BiocManager::install("ShortRead")
BiocManager::install("Biostrings")

library(dada2) 
packageVersion("dada2") #1.14.1 #confirm that DADA2 version is 1.8 or later (1.10)
library(ShortRead)
packageVersion("ShortRead") #1.46.0
library(Biostrings)
packageVersion("Biostrings") #2.56.0

#Set working directory to the dataset folder
path <- "~/home/pkushwaha/SBAR_May2022_ITS/demux"

#Set working directory to the unzipped mock dataset folder, here located in downloads
#Here this is done by defining the variable path
list.files(path)
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))

#### Read in the samples files ####
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 7)
sample.names <- sapply(strsplit(basename(sample.names), "\\.f"), `[`, 1)
sample.names <- paste("S",sample.names,sep = "_")

head(sample.names)

#Identify primers
FWD <- "CTTGGTCATTTAGAGGAAGTAA"  #forward primer sequence
REV <- "GCTGCGTTCTTCATCGATGC"  #reverse primer sequence

allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}

FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)

FWD.orients
REV.orients


fnFs.filtN <- file.path("/home/rstudio", "filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
fnRs.filtN <- file.path("/home/rstudio", "filtN", basename(fnRs))

for(i in seq_along(fnFs)) {
  cat("Sample", i, "\n")
  cat("fnF:", fnFs[[i]], "-- fnR:", fnRs[[i]], "\n")
  filterAndTrim(fnFs[[i]], fnFs.filtN[[i]], fnRs[[i]], fnRs.filtN[[i]], 
                maxN = 0, multithread = TRUE, matchIDs = TRUE)
}

##Keep updating as erros appear
for(i in 26:89) {
  cat("Sample", i, "\n")
  cat("fnF:", fnFs[[i]], "-- fnR:", fnRs[[i]], "\n")
  filterAndTrim(fnFs[[i]], fnFs.filtN[[i]], fnRs[[i]], fnRs.filtN[[i]], 
                maxN = 0, multithread = TRUE, matchIDs = TRUE)
}
#filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE)

#Count the number of times the primers appear in the forward and reverse read, 
#while considering all possible primer orientations
primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))

#Remove primers (cutadapt tool)
#CHANGE the cutadapt path on your machine
#Run shell commands from R (version 1.18)
cutadapt <- "/usr/local/bin/cutadapt" 
system2(cutadapt, args = "--version")

path.cut <- file.path("/home/rstudio", "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))


FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC) 
# Run Cutadapt
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i])) # input files
}

#Count the presence of primers in the first cutadapt-ed samples
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))

#Primers should not be detected in the cutadapted reads now.

cutFs <- sort(list.files(path.cut, pattern="_R1_001.fastq", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern="_R2_001.fastq", full.names = TRUE))

#Extract sample names, assuming filenames have same format:
sample.names <- sapply(strsplit(basename(cutFs), "_"), `[`, 6)
sample.names <- sapply(strsplit(basename(sample.names), "\\.f"), `[`, 1)
head(sample.names)

#Inspect read quality profiles
plotQualityProfile(cutFs[1:15])
plotQualityProfile(cutRs[1:20])

filtFs <- file.path(path.cut, "filtered", basename(cutFs))
filtRs <- file.path(path.cut, "filtered", basename(cutRs))

out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, maxN = 0, maxEE = c(2, 2), 
                     truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = TRUE)
head(out)
View(out)

#Learn the error rates
errF <- learnErrors(filtFs[61:74], multithread = TRUE)
errR <- learnErrors(filtRs[61:74], multithread = TRUE)

#Plot error rates
plotErrors(errF, nominalQ = TRUE)
plotErrors(errR, nominalQ = TRUE)


#Dereplicate identical reads
exists <- file.exists(filtFs)
derepFs <- derepFastq(filtFs[exists], verbose = TRUE)
derepRs <- derepFastq(filtRs[exists], verbose = TRUE)

#Name the derep-class objects by the sample names
names(derepFs) <- sample.names[exists]
names(derepRs) <- sample.names[exists]

#Sample inference
dadaFs <- dada(derepFs, err = errF, multithread = TRUE)
dadaRs <- dada(derepRs, err = errR, multithread = TRUE)

#Merge paired reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
head(mergers[1])

#Construct sequence (ASV) table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))

#Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
table(nchar(getSequences(seqtab.nochim)))

#Determine the percentage of total sequence reads that are NOT chimeras 
#Percentage of reads that are chimeras should not be too high
sum(seqtab.nochim)/sum(seqtab)

#Export ASV file with reads
write.table(seqtab.nochim, "All_Seqs_ITS_reads_sbar_may2022_08172022.txt", sep="\t", quote=F, row.names=T, col.names=NA)
saveRDS(seqtab.nochim, "All_Seqs_ITS_reads_sbar_may2022_08172022.rds")

#Track reads through the pipeline
#trimmed, filtered, chimeria removal, etc.
getN <- function(x) sum(getUniques(x))
#track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
track <- cbind(sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
#colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
colnames(track) <- c("denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

#Export a table with reads that have made through the pipeline
write.table(track, "All_Seqs_Denoised_nonchim_ReadsCounts_ITS_sbar_may2022_08172022.txt", sep="\t", quote=F, row.names=T, col.names=NA)
write.table(out, "All_Seqs_filtered_only_ReadsCounts_ITS_sbar_may2022_08172022.txt", sep="\t", quote=F, row.names=T, col.names=NA)

#Assign taxonomy using UNITE ITS database
#IMPORTANT - go to this site https://benjjneb.github.io/dada2/training.html and download the training set
#Place it in the same working directory 'path' before starting the command. 
taxa <- assignTaxonomy(seqtab.nochim, "~/home/pkushwaha/SBAR_May2022_ITS/sh_general_release_10.05.2021_PK", multithread=TRUE, tryRC=T)

#Removing sequence rownames for display only
taxa.print <- taxa  
rownames(taxa.print) <- NULL
head(taxa.print)
View(taxa.print)

#Export ASV and taxonomy tables to path for upload and analysis
write.table(taxa, "All_Seqs_unite_ITS_sbar_may2022_08172022.txt", sep="\t", quote=F, row.names=T, col.names=NA)

##Merge and taxonomic information
fun.asv.tax<-cbind(t(seqtab.nochim), taxa)
dim(fun.asv.tax) #2365 31

#Save the file
write.table(fun.asv.tax, file="asv_table_Taxa_ITS_sbar_may2022_08172022.txt", sep="\t", quote=F, row.names=T, col.names=NA)
