setwd("~/Desktop/Github/seed-microbes")

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


#old methods############################
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

#ITS, 16S, and both - final##########################################
prune.ps <- function(ps, cutoff=0.01,prop_samp=0.1){
  threshold <- prop_samp*nsamples(ps)
  ps <-filter_taxa(ps,function(x)sum(x>0)>=threshold,TRUE)
  ps <- prune_taxa(names(which(taxa_sums(ps)*100/sum(taxa_sums(ps)) > cutoff)),ps)
  return(ps)
}

ps.16S.raw<- readRDS(file = "rds/16S_crn.ps-new.Rds")
ps.16S.new<- DESeq_varstab(ps.16S.raw, ~1)
ps.16S.new10 <- prune.ps(ps.16S.raw)
saveRDS(ps.16S.new, file = "rds/16S_crn.ps.deseq-new.Rds")
saveRDS(ps.16S.new10, file = "rds/16S_crn.ps.deseq-new10.Rds")

ps.both.raw<- readRDS(file = "rds/combined.both_crn.ps-new.Rds")
ps.bothnew<- DESeq_varstab(ps.both.raw, ~1)
ps.both.new10 <- prune.ps(ps.both.raw)
saveRDS(ps.bothnew, file = "rds/combined.both_crn.ps.deseq-new.Rds")
saveRDS(ps.bothnew10, file = "rds/combined.both_crn.ps.deseq-new10.Rds")

ps.ITS.raw <- readRDS(file = "rds/ITS_crn.ps-new.Rds")
ps.ITS.new <- DESeq_varstab(ps.ITS.raw, ~1)
ps.ITS.new10 <- prune.ps(ps.ITS.raw)
saveRDS(ps.ITS.new, file = "rds/ITS_crn.ps.deseq-new.Rds")



