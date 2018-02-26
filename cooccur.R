##16S############################
###importing in the data set
pscooccur1 <- readRDS(file = "rds/pscooccur1.rds")

#Running the cooccur function for family:
pscooccur16S_pruned <- prune_samples(names(which(sample_sums(pscooccur1) >= 0)), pscooccur1)
pscooccur16S_pruned <-subset_samples(pscooccur16S_pruned, type != "DNA_control")
pscooccur16S_pruned_pa <- transform_sample_counts(pscooccur16S_pruned,function(x)1*(x>0))


#Creating the otu table from the pruned phyloseq object for family:
otutab16S <- otu_table(pscooccur16S_pruned_pa) 

#Running the cooccur function for family:
cooccur.species <- cooccur(mat = otutab16S,
                      type = "spp_site")
saveRDS(cooccur.species, file = file = "rds/cooccur.species.rds")
           
##ITS##########################
###importing in the data set
psITS <- readRDS(file = "rds/ITS_crn.ps.deseq-new.Rds")

#Running the cooccur function for family:
pscooccurITS_pruned <- prune_samples(names(which(sample_sums(pscooccur1) >= 0)), pscooccur1)
pscooccurITS_pruned <-subset_samples(pscooccurITS_pruned, type != "DNA_control")
pscooccurITS_pruned_pa <- transform_sample_counts(pscooccurITS_pruned,function(x)1*(x>0))


#Creating the otu table from the pruned phyloseq object for family:
otutabITS <- otu_table(pscooccurITS_pruned_pa) 

#Running the cooccur function for family:
cooccur_ITS <- cooccur(mat = otutabITS,
                       type = "spp_site")
saveRDS(cooccur_ITS, file = file = "rds/cooccur.ITS.rds")

##both######################
###importing in the data set
psboth <- readRDS(file = "rds/combined.both_crn.ps.deseq-new.Rds")

#Running the cooccur function for family:
pscooccurboth_pruned <- prune_samples(names(which(sample_sums(pscooccur1) >= 0)), pscooccur1)
pscooccurboth_pruned <-subset_samples(pscooccurboth_pruned, type != "DNA_control")
pscooccurboth_pruned_pa <- transform_sample_counts(pscooccurboth_pruned,function(x)1*(x>0))


#Creating the otu table from the pruned phyloseq object for family:
otutab.both <- otu_table(pscooccurboth_pruned_pa) 

#Running the cooccur function for family:
cooccur_both <- cooccur(mat = otutab.both,
                       type = "spp_site")

saveRDS(cooccur_both, file = file = "rds/cooccur.both.rds")