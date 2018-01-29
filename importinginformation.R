setwd("~/Documents")

##to create biom... 
## python.call(make_otu_table.py -i fileascsv.txt -o file.biom -t rep_set tax_assignments.txt)
## python.call(biom convert -i file.tsv -o file.biom --table-type="OTU table" --to-json) )

library(cooccur)
library(phyloseq)
library(ggplot2)
library("DESeq2")
library(vegan)
library(reshape2)
library(psych)
library(biom)
library(wesanderson)
library(knitr)
library(sp)
library(ade4)
library(igraph)
library(SpiecEasi)

#My OTU table- old method using OTU clustering
pscooccurold <- import_biom("16S_otu_table_newjson.biom",
                         parseFunction=parse_taxonomy_greengenes)

pscooccurold1 <- prune_samples(names(which(sample_sums(pscooccurold) >= 0)),pscooccurold)
pscooccurold1 <- prune_taxa(taxa_sums(pscooccurold1) > 11, pscooccurold1)
sd <- import_qiime_sample_data("corn_rox.tsv")
pscooccurold1 <- merge_phyloseq(pscooccurold1, sd)
pscooccurold1 <- DESeq_varstab(pscooccurold1, ~1)
rename_otus(pscooccurold1)

#save the file as something that can be recalled when doing other fuctions
saveRDS(pscooccurold1, file = "rds/pscooccurold1.rds")
#to restore: readRDS(file = "rds/pscooccurold1.rds")

##########################################################################


# Import sequence table generated from dada2, to be used as otu_table in phyloseq
psnew <- readRDS('rds/16S_seqtab.Rds')
# Rename samples names to include "Sample" in front of each sample number
rownames(psnew) <- paste("Sample",rownames(psnew),sep="")

# Import taxa table and create phyloseq object
taxtabnew <- readRDS('16S_silva_taxtab_spp.Rds')
ps <- phyloseq(otu_table(psnew, taxa_are_rows=FALSE),
               tax_table(taxtabnew)
               )

# Add mapping file, which contains sample data. The 'merge_phyloseq()' function automatically subsets
# the phyloseq object to only include samples in the mapping file
sd <- import_qiime_sample_data("corn_rox.tsv")
pscooccurnew <- merge_phyloseq(ps, sd)

#Accounting for error within the samples
pscooccur1 <- prune_samples(names(which(sample_sums(pscooccurnew) >= 0)),pscooccurnew)
#pscooccur1 <- prune_taxa(taxa_sums(pscooccur1) > 11, pscooccur1)
#sd <- import_qiime_sample_data("corn_rox.tsv")
#pscooccur1 <- merge_phyloseq(pscooccur1, sd)
pscooccur1 <- DESeq_varstab(pscooccur1, ~1)

rename_otus(pscooccur1)
#save the file as something that can be recalled when doing other fuctions
saveRDS(pscooccur1, file = "rds/pscooccur1.rds")
#to restore: readRDS(file = "rds/pscooccur1.rds")


ps.both.raw<- readRDS(file = "rds/combined.both_crn.ps.Rds")
ps.bothnew<- DESeq_varstab(ps.bothnew, ~1)
saveRDS(ps.bothnew, file = "rds/combined.both_crn.ps.deseq.Rds")

ps.ITS.raw <- readRDS(file = "rds/ITS_crn.ps.Rds")
ps.ITSnew <- DESeq_varstab(ps.ITS.raw, ~1)
saveRDS(ps.ITSnew, file = "rds/ITS_crn.ps.deseq.Rds")



