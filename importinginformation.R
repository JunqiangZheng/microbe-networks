##to create biom... create biom -i___ -o____ --table-type="OTU table" --to-json
##this can also be done with pythoninr

setwd("~/Documents")

library(phyloseq)
library(ggplot2)
library("DESeq2")
library(vegan)
library(psych)
library(biom)

###### Test for "sequencing run" as a factor determining sample variation
#
ps_combined = import_biom(16S_otu_table_newjson.biom)
ps_combined_rar = rarefy_even_depth(ps_combined, sample.size = 1000)
ps_combined_DESeq = DESeq_varstab(ps_combined, ~Plate)
ps_combined_DESeq
compare_runs = subset_samples(ps_combined_rar,Compare_runs == 1)
compare_runs = prune_taxa(taxa_sums(compare_runs) > 0, compare_runs)
compare_runs_bray = phyloseq::distance(compare_runs, method = "bray")
compare_runsdf = data.frame(sample_data(compare_runs))
adonis(compare_runs_bray ~ Illumina_run,data=compare_runsdf)


#Prepare phyloseq object
ps_16S = import_biom("16S_otu_table_newjson.biom",
                     #Lucas removed: refseqfilename = "16S_rep_set.fasta",
                     treefilename = "16S_otus.tre",
                     parseFunction=parse_taxonomy_greengenes)