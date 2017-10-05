setwd("~/Documents")

##to create biom... 
## python.call(make_otu_table.py -i fileascsv.txt -o file.biom -t rep_set tax_assignments.txt)
## python.call(biom convert -i file.tsv -o file.biom --table-type="OTU table" --to-json) )

library(cooccur)
library(phyloseq)
library(ggplot2)
library("DESeq2")
library(vegan)
library(psych)
library(biom)
library(wesanderson)
library(knitr)

#My OTU table
pscooccur <- import_biom("16S_otu_table_newjson.biom",
                         parseFunction=parse_taxonomy_greengenes)

#Accounting for error within the samples
pscooccur1 <- prune_samples(names(which(sample_sums(pscooccur) >= 0)),pscooccur)
pscooccur1 <- prune_taxa(taxa_sums(pscooccur1) > 11, pscooccur1)
sd <- import_qiime_sample_data("corn_rox.tsv")
pscooccur1 <- merge_phyloseq(pscooccur1, sd)
pscooccur1 <- DESeq_varstab(pscooccur1, ~1)

#save the file as something that can be recalled when doing other fuctions
saveRDS(pscooccur1, file = "pscooccur1.rds")
#to restore: readRDS(file = "pscooccur1.rds")
