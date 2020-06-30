# REU DADA2 Workflow
# June 2020
# Author: Kyle Spitler

#clear R brain
rm(list=ls())

# Call libraries
library(dada2); packageVersion("dada2")

# set the working directory to the directory
path <- "~/Desktop/DukeREU/FastQ_Files"
list.files("~/Desktop/DukeREU/FastQ_Files")

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files("~/Desktop/DukeREU/FastQ_Files", pattern = "_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files("~/Desktop/DukeREU/FastQ_Files", pattern = "_R2_001.fastq.gz", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# visualizing the quality profiles of the forward reads
plotQualityProfile(fnFs[1:2])

# visualizing the quality profiles of the reverse reads
plotQualityProfile(fnRs[1:2])

# Place filtered files in filtered/ subdirectory
filtFs <- file.path("~/Desktop/DukeREU/FastQ_Files", "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path("~/Desktop/DukeREU/FastQ_Files", "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(230,160),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
head(out)

# learned error model from the data, by alternating estimation of the error rates and inference of sample composition until they converge on a jointly consistent solution
# forward

errF <- learnErrors(filtFs, multithread=TRUE)

#reverse

errR <- learnErrors(filtRs, multithread=TRUE)

# visualize

plotErrors(errF, nominalQ=TRUE)

# apply core sample interference algorithm

dadaFs <- dada(filtFs, err=errF, multithread=TRUE)

dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

dadaFs[[1]]

# merge the forward and reverse reads together to obtain the full denoised sequences

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

# Inspect the merger data.frame from the first sample

head(mergers[[1]])

# construct Amplicon Sequence Varairnt Table

seqtab <- makeSequenceTable(mergers)
dim(seqtab)



# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

# remove  chimeras

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

sum(seqtab.nochim)/sum(seqtab)

# look at the number of reads that made it through each step in the pipeline

getN <- function(x) sum(getUniques(x))

track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))

# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)

colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")

rownames(track) <- sample.names

head(track)

# assign taxonomy

taxa <- assignTaxonomy(seqtab.nochim, "~/Desktop/DukeREU/Taxonomy/silva_nr_v138_train_set.fa.gz", multithread=TRUE)

taxa <- addSpecies(taxa, "~/Desktop/DukeREU/Taxonomy/silva_species_assignment_v138 (1).fa.gz")

taxa.print <- taxa 

# Removing sequence rownames for display only

rownames(taxa.print) <- NULL

head(taxa.print)

#### Import into phyloseq: ######

library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")

theme_set(theme_bw())

# construct a simple sample data.frame from the information encoded in the filenames

samples.out <- rownames(seqtab.nochim)
subject <- sapply(strsplit(samples.out, "D"), `[`, 1)
gender <- substr(subject,1,1)
subject <- substr(subject,2,999)
day <- as.integer(sapply(strsplit(samples.out, "D"), `[`, 2))
samdf <- data.frame(Subject=subject, Gender=gender, Day=day)
samdf$When <- "Early"
samdf$When[samdf$Day>100] <- "Late"
rownames(samdf) <- samples.out

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))
ps <- prune_samples(sample_names(ps) != "Mock", ps) # Remove mock sample

dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps

saveRDS(ps, file = "ps.rds")

loaded.ps = readRDS("ps.rds")
print(loaded.ps)

rank_names(ps)

table(tax_table(ps)[, "Phylum"], exclude = NULL)

# remove NA reads from taxa table
ps <- subset_taxa(ps, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))

# Compute prevalence of each feature, store as data.frame
prevdf = apply(X = otu_table(ps),
               MARGIN = ifelse(taxa_are_rows(ps), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})

# Add taxonomy and total read counts to this data.frame
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(ps),
                    tax_table(ps))

plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})

# Subset to the remaining phyla
ggplot(prevdf, aes(TotalAbundance, Prevalence / nsamples(ps),color=Phylum)) +
# Include a guess for parameter
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) +  geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")

