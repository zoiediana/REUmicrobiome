# REU DADA2 Workflow
# June 2020
# Author: Kyle Spitler

#clear R brain
rm(list=ls())

knitr::opts_chunk$set(echo = TRUE)

library(readr)
library(fs)
library(R.utils)
library(tidyverse)
library(stringr)
library(phyloseq)
library(ggplot2)
library(dada2)
library(magrittr)

baseFolder <- "~/Desktop/DukeREU/FastQ_Files/Micro_F/High_7wks"

outputFiltered <- "~/Desktop/DukeREU/FastQ_Files/Micro_F/High_7wks/filtered_1"
output.dir <- "~/Desktop/DukeREU/FastQ_Files/Micro_F/High_7wks/output_1"
ps.file <- file.path(output.dir,"ps_MP_high.rds")
track.file <- file.path(output.dir,"track_MP_high.rds")

R1.fileName <- list.files(baseFolder, pattern = "R1_001.fastq.gz")

R2.fileName <- str_replace(R1.fileName,"R1","R2")
sampleInfo <- data.frame(R1.fileName, R2.fileName)
sampleInfo %<>%
  separate(R1.fileName, into = "sampleName", sep = "_", remove=FALSE)

i<-7
plotQualityProfile(file.path(baseFolder, sampleInfo$R1.fileName[i]))
plotQualityProfile(file.path(baseFolder, sampleInfo$R2.fileName[i]))

filtFs <- file.path(outputFiltered , paste0(sampleInfo$sampleName, "_F_filt.fastq.gz"))
filtRs <- file.path(outputFiltered , paste0(sampleInfo$sampleName, "_R_filt.fastq.gz"))
fnFs <- file.path(baseFolder,as.character(sampleInfo$R1.fileName))
fnRs <- file.path(baseFolder,as.character(sampleInfo$R2.fileName))

filt.out <- filterAndTrim(fwd = fnFs, filt = filtFs, rev = fnRs, filt.rev = filtRs, trimLeft=10, truncLen=c(230,230),
                          maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                          compress=TRUE, multithread=FALSE)
filt.out

filtFs <- filtFs[file_exists(filtFs)]
filtRs <- filtRs[file_exists(filtRs)]
sample.names <- filtFs  %>% 
  basename %>%
  str_replace("_F_filt.fastq.gz","") 

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

names(derepFs) <- sampleInfo$sampleName
names(derepRs) <- sampleInfo$sampleName

dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

seqtab <- makeSequenceTable(mergers)

dim(seqtab)

table(nchar(getSequences(seqtab)))

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

getN <- function(x) sum(getUniques(x))
track <- filt.out %>%
  as.data.frame %>%
  rownames_to_column %>%
  mutate(rowname=str_replace(rowname, "_R1.fastq.gz","")) %>%
  rename(sample=rowname, input=reads.in, filtered=reads.out)

dadaFs %>%
  sapply(getN) %>%
  as_tibble(rownames = NA) %>%
  rownames_to_column %>%
  rename(sample=rowname, denoised=value) ->
  denoised

track %<>% full_join(denoised, by=c("sample"))

mergers %>%
  sapply(getN) %>%
  data.frame %>%
  rownames_to_column %>%
  rename(sample=rowname, merged='.') ->
  merged

track %<>% full_join(merged, by=c("sample"))

seqtab %>%
  rowSums %>%
  data.frame %>%
  rownames_to_column() %>%
  rename(sample=rowname, tabled = '.') -> 
  tabled

track %<>% full_join(tabled, by=c("sample"))

seqtab.nochim %>%
  rowSums %>%
  data.frame %>%
  rownames_to_column() %>%
  rename(sample=rowname, nochim='.') -> 
  nochim

track %<>% full_join(nochim, by=c("sample"))

track %>%
  arrange(as.numeric(sample))

write_rds(track, track.file)

taxa <- assignTaxonomy(seqtab.nochim, "~/Desktop/DukeREU/Taxonomy/silva_nr_v138_train_set.fa.gz", multithread=TRUE)

taxa <- addSpecies(taxa, "~/Desktop/DukeREU/Taxonomy/silva_species_assignment_v138 (1).fa.gz")

otus = otu_table(seqtab.nochim, taxa_are_rows=FALSE)

sampleInfo %>%
  select(sampleName,R1.fileName) %>%
  column_to_rownames(var = "sampleName") ->
  sd

ps <- phyloseq(otus, sample_data(sd), tax_table(taxa))

print(ps)

write_rds(ps, ps.file)

print(read_rds(ps.file))

plot_bar(ps)

plot_bar(ps, fill="Kingdom")

plot_bar(ps, fill="Phylum")

p = plot_bar(ps, fill="Phylum")
p + geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack")

