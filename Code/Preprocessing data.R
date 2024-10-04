### IMPORT LIBRARIES
library(biomformat)
library(tidyverse)
library(microeco)
library(magrittr)
library(ape)
library(seqinr)

### IMPORT DATA

setwd("../Data/")

## IMPORT BIOM DATA
biomdata <- biomformat::read_biom("otu_table.blast_NCBI_16S.biom")

## OTU TABLE
otu_table <- as(biom_data(biomdata), "matrix")
otu_table <- as.data.frame(otu_table)

## SAMPLE TABLE
sample <- sample_metadata(biomdata)
sample %<>% mutate(Group = case_when(Group1 == "Case1" ~ "G1", 
                                     Group1 == "Case2" ~ "G2",
                                     Group1 == "Case3" ~ "G3",
                                     TRUE ~ "G4"))

## TAX TABLE
tax_table <- observation_metadata(biomdata)

tax <- data.frame(tax_table[[names(tax_table)[1]]])
tax <- t(tax)
data(taxonomy_table_16S)
colnames(tax) <- colnames(taxonomy_table_16S)
rownames(tax) <- names(tax_table)[1]

for (i in names(tax_table)){
  tab <- t(data.frame(tax_table[[i]]))
  rownames(tab) <- i
  if (ncol(tab) == 7){
    tax <- rbind(tax, tab)
  }
  
}

tax_table <- as.data.frame(tax)
tax_table %<>% tidy_taxonomy

## FASTA DATA
s <- read.fasta("otus_rep.fasta")
names(s) <- rownames(otu_table)

### CREATE MICROTABLE
mt <- microtable$new(otu_table = otu_table,
                     tax_table = tax_table,
                     sample_table = sample,
                     rep_fasta = s)
mt$tidy_dataset()

## SAVE MICROTABLE
saveRDS(mt, file = "metagenome.Rds")
