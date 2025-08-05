
#                              ITS_ANALYSE MICROBOME DADA2 & PHYLOSEQ

# First install a recent version of Rstudio: https://www.rstudio.com/products/rstudio/download/#download

# Install last version of R : https://cran.r-project.org/bin/windows/base/
#Select it in Tools Global Options Rversion in Rstudio.

#Installation of packages for Microbiome analyses usig DADA2

#command to install using BioManager
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")


BiocManager::install("phyloseq")
BiocManager::install("DECIPHER")
BiocManager::install("Biostrings")
BiocManager::install("ggplot2")
BiocManager::install("dplyr")
BiocManager::install("dada2")
BiocManager::install("DESeq2")
BiocManager::install("msa")
BiocManager::install("phangorn")
BiocManager::install("devtools")
BiocManager::install("Tax4Fun")
BiocManager::install("Rtools")
install.packages("microbial")
BiocManager::install("edgeR")
install.packages("metagMisc")

#Load libraries

library(phyloseq)
library(DECIPHER)
library(dada2)
library(vegan)
library(Biostrings)
library(ggplot2)
library(dplyr)
library(DESeq2)
library(msa)
library(phangorn)
library(microbial)
#install Tax4Fun manually after downloading from http://tax4fun.gobics.de/
#(need to install first Rtools from https://cran.r-project.org/bin/windows/Rtools/rtools42/rtools.html)
library(Tax4Fun)
library(qiime2R)
library("MicEco")

setwd("your_path/ITS/fastq")

#starting the dada2 pipeline
path <- "your_path/ITS/fastq"  #CHANGE ME to the directory containing the fastq files after unzipping.
list.files(path)

#DETECT AND REMOVE PRIMERS
# Forward and reverse fastq filenames have format: SAMPLENAME-16S_1.fastq and SAMPLENAME-16S_2.fastq
fnFs <- sort(list.files(path, pattern="_1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`,1)
sample.names


#Inspect read quality profiles
plotQualityProfile(fnFs[33:34])

plotQualityProfile(fnRs[1:2])

#Filter and trim: Assign the filenames for the filtered fastq.gz files.

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

#Filter:
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(280,205), maxN=0, maxEE=c(2,2), trimLeft = c(20, 20), truncQ=2, rm.phix=TRUE,compress=TRUE,  multithread=FALSE)
# On Windows set multithread=FALSE; attention ici on trime les primers de 20 pb en amont de chaque sequence for et Rev (trimleft)


head(out) 
write.csv(out, "Stat_filteringITS_SBL.csv")


#Learn the Error Rates

errF <- learnErrors(filtFs, randomize= TRUE,  multithread=TRUE)
errR <- learnErrors(filtRs, randomize= TRUE, multithread=TRUE)

plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

#Sample Inference
dadaFs <- dada(filtFs, err=errF, pool=FALSE, multithread=TRUE)

dadaRs <- dada(filtRs, err=errR, pool=FALSE, multithread=TRUE)

#Inspecting the returned dada-class object:
dadaFs[[1]]

#Merge paired reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

# Inspect the merger data.frame from the first sample

head(mergers[[1]])

#Construct sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)


## [1]  20 293
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
## 
## 251 252 253 254 255 
##   1  88 196   6   2

#Remove non-target-length sequences from your sequence table (eg. seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 250:256]). This is analogous to "cutting a band" in-silico to get amplicons of the targeted length. 
seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 260:433]


#Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

## [1]  20 232
sum(seqtab.nochim)/sum(seqtab)
## [1] 0.964263

#Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

write.csv(track, "Stat_assemblyDADA2ITS.csv")
saveRDS(seqtab.nochim, "seqtabnochim-ITS-Lea.RDS")


seqtab.nochim <- readRDS(file="seqtabnochim-ITS-Lea.RDS")
taxa <- readRDS(file="taxaITS-SBLnewtax.RDS")

#Assign taxonomy (download UNITE database_all_eukaryotes)
taxa <- assignTaxonomy(seqtab.nochim, "yourpath../Taxonomic_databases/sh_general_release_dynamic_s_all_19.02.2025_dev.fasta", multithread=TRUE)


#inspect the taxonomic assignments:
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

View(sample.names)


saveRDS(taxa, "SBLtaxaITS.RDS")

write.csv(taxa, "SBL_ITS.csv")
write.csv(seqtab.nochim, "seqtab.nochim_SBL.csv")



# extract info in a single table (without metadata)
seqtab6 <-as.data.frame(t(seqtab.nochim))
seqtab6$ID <- rownames(seqtab6)
taxa6<-as.data.frame(taxa)
taxa6$ID <- rownames(taxa)
all_data5  <- merge(taxa6,seqtab6,by='ID')
write.csv(all_data5, "alldata_ITS_SBL_UNITE2025.csv")
