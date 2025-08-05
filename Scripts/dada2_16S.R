
#                             16S_ANALYSE MICROBOME DADA2 & PHYLOSEQ

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

n#for Tax4Fun, install manually (download archive from http://tax4fun.gobics.de/)
#Tax4Fun_0.3.1.tar.gz
#install qiime2R
devtools::install_github("jbisanz/qiime2R")

#Barplots & Venn diagram packages
install.packages("microbiome")

install.packages("git2r")
install.packages("remotes")
remotes::install_github("Russel88/MicEco")




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

setwd("your_path/16S")

#starting the dada2 pipeline
path <- "your_path/16S"  #CHANGE ME to the directory containing the fastq files after unzipping.
list.files(path)



# Forward and reverse fastq filenames have format: SAMPLENAME-16S_1.fastq and SAMPLENAME-16S_2.fastq
fnFs <- sort(list.files(path, pattern="_1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_2.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
sample.names
#Inspect read quality profiles
plotQualityProfile(fnFs[1:2])

plotQualityProfile(fnRs[1:2])

#Filter and trim: Assign the filenames for the filtered fastq.gz files.

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

#Filter:  Macrogen 16S primers are 341F :  CCTACGGGNGGCWGCAG (17 nt) and  805R GACTACHVGGGTATCTAATCC (21 nt)

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(280,205), maxN=0, maxEE=c(2,2), trimLeft = c(17, 21), truncQ=2, rm.phix=TRUE,compress=TRUE,  multithread=FALSE)
# On Windows set multithread=FALSE; attention ici on trime les primers de 20 pb en amont de chaque sequence for et Rev (trimleft)


head(out) 
#Learn the Error Rates

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

#Sample Inference
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)

dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

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
seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 380:435]

saveRDS(seqtab, "16SLeaCambodia.RDS")

#Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab2, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
View(seqtab.nochim)


## [1]  20 232
sum(seqtab.nochim)/sum(seqtab2)
## [1] 0.964263

write.csv(seqtab.nochim, "seqtab.nochim16SLeaCambodia.csv")
saveRDS(seqtab.nochim, "seqtab.nochim16SLeaCambodia.RDS")


#Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

write.csv(track, "Stat_assemblyDADA2_16SLeaCambodia.csv")


#Assign taxonomy (download first silva_nr_v132_train_set.fa.gz from silva website)
taxa <- assignTaxonomy(seqtab.nochim, "your_path/Taxonomic_databases/silva_nr_v138_train_set.fa.gz", multithread=TRUE)

#To add species information (download first silva_species_assignment_v132.fa.gz from silva website)
taxa <- addSpecies(taxa, "your_path/Taxonomic_databases/silva_species_assignment_v138.fa.gz")


#Assign taxonomy (different databses options given, most recent and updated is silva138)
#taxa<- assignTaxonomy(seqtab.nochim, "C:/Users/moulinl/Documents/Scripts R/Taxonomic_databases/silva_nr_v132_train_set.fa.gz", multithread=TRUE)
#taxa2 <- assignTaxonomy(seqtab.nochim, "C:/Users/moulinl/Documents/Scripts R/Taxonomic_databases/silva_nr_v138_train_set.fa.gz", multithread=TRUE)
#taxa3 <- assignTaxonomy(seqtab.nochim, "C:/Users/moulinl/Documents/Scripts R/Taxonomic_databases/rdp_train_set_18.fa.gz", multithread=TRUE)

#To add species information (download first silva_species_assignment_v132.fa.gz from silva website)
#taxa1 <- addSpecies(taxa, "C:/Users/moulinl/Documents/Scripts R/Taxonomic_databases/silva_species_assignment_v132.fa.gz")
#taxa2 <- addSpecies(taxa2, "C:/Users/moulinl/Documents/Scripts R/Taxonomic_databases/silva_species_assignment_v138.fa.gz")
#taxa3 <- addSpecies(taxa3, "C:/Users/moulinl/Documents/Scripts R/Taxonomic_databases/rdp_species_assignment_18.fa.gz")

#inspect the taxonomic assignments:
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

#save seqtab and taxa as RDS files (easier to reload later)
saveRDS(taxa, file="taxa16S_LeaCambodge.RDS")
saveRDS(seqtab.nochim, file="16Sseqtab.nochimLea.RDS")

#load metadata file
samdf <- read.csv2("your_path/metadata16S.csv")
samdf

all(rownames(seqtab.nochim) %in% samdf$sample_id) # TRUE pour verifier que les noms concordent

rownames(samdf) <- samdf$sample_id
keep.cols <- c("sample_id", "colname1", "colname2", "colname3", "colname4", "colname5") # to be changed with your colnames
samdf <- samdf[rownames(seqtab.nochim), keep.cols]

samples.out <- rownames(seqtab.nochim)
rownames(samdf) <- samples.out


# check correespondance between samplenames between seqtab and samdf
all(rownames(seqtab.nochim) %in% samdf$sample_id) # TRUE
all(samdf$sample_id %in% rownames(seqtab.nochim)) # TRUE

# check correespondance between taxanames between seqtab and taxa
all(colnames(seqtab.nochim) %in% rownames(taxa)) # TRUE
all(rownames(taxa) %in% colnames(seqtab.nochim)) # TRUE

#save in RDS format the seqtab 

saveRDS(seqtab.nochim, "seqtab.nochim16S.RDS")
saveRDS(taxa, "taxa16S.RDS")
write.csv(taxa, "taxa16S.csv")
write.csv(seqtab.nochim, "seqtab.nochim16S.csv")



# extract info in a single table (without metadata)
seqtab6 <-as.data.frame(t(seqtab.nochim))
seqtab6$ID <- rownames(seqtab6)
taxa6<-as.data.frame(taxa)
taxa6$ID <- rownames(taxa)
all_data5  <- merge(taxa6,seqtab6,by='ID')
write.csv(all_data5, "Alldata_16S_LeaCambodia_silva138.csv")


#phyloseq object without phylogenetic tree

ps_16S <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
                   sample_data(samdf), 
                   tax_table(taxa))

ps_16S
