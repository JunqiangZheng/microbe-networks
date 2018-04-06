##16S############################
###importing in the data set
ps.16S <- readRDS(file = "rds/16S_crn.ps.deseq-new10.Rds")

#Running the cooccur function for family:
pscooccur16S_pruned <- prune_samples(names(which(sample_sums(ps.16S) >= 0)), ps.16S)
pscooccur16S_pruned <-subset_samples(pscooccur16S_pruned, type != "DNA_control")
pscooccur16S_pruned_pa <- transform_sample_counts(pscooccur16S_pruned,function(x)1*(x>0))


#Creating the otu table from the pruned phyloseq object for family:
otutab16S <- otu_table(pscooccur16S_pruned_pa) 

#Running the cooccur function for family:
cooccur.species <- cooccur(mat = otutab16S,
                           spp_names = TRUE)
saveRDS(cooccur.species, file = "rds/cooccur.species.10.rds")
           
##ITS##########################
###importing in the data set
ps.ITS <- readRDS(file = "rds/ITS_crn.ps.deseq-new10.Rds")

#Running the cooccur function for family:
pscooccurITS_pruned <- prune_samples(names(which(sample_sums(ps.ITS) >= 0)), ps.ITS)
pscooccurITS_pruned <-subset_samples(pscooccurITS_pruned, type != "DNA_control")
pscooccurITS_pruned_pa <- transform_sample_counts(pscooccurITS_pruned,function(x)1*(x>0))


#Creating the otu table from the pruned phyloseq object for family:
otutabITS <- otu_table(pscooccurITS_pruned_pa) 

#Running the cooccur function for family:
cooccur_ITS <- cooccur(mat = otutabITS,
                       spp_names = TRUE)
saveRDS(cooccur_ITS, file = "rds/cooccur.ITS.10.rds")

##both######################
###importing in the data set
psboth <- readRDS(file = "rds/combined.both_crn.ps.deseq-new10.Rds")

#Running the cooccur function for family:
pscooccurboth_pruned <- prune_samples(names(which(sample_sums(psboth) >= 0)), psboth)
pscooccurboth_pruned <-subset_samples(pscooccurboth_pruned, type != "DNA_control")
pscooccurboth_pruned_pa <- transform_sample_counts(pscooccurboth_pruned,function(x)1*(x>0))


#Creating the otu table from the pruned phyloseq object for family:
otutab.both <- otu_table(pscooccurboth_pruned_pa) 

#Running the cooccur function for family:
cooccur_both <- cooccur(mat = otutab.both,
                        spp_names = TRUE)

saveRDS(cooccur_both, file = "rds/cooccur.both.10.rds")

