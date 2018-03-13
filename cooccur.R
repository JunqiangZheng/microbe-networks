##16S############################
###importing in the data set
pscooccur1 <- readRDS(file = "rds/16S_crn.ps.deseq-new.Rds")
pscooccur1.10 <- readRDS(file = "rds/ITS_crn.ps.deseq-new10.Rds")

#Running the cooccur function for family:
pscooccur16S_pruned <- prune_samples(names(which(sample_sums(pscooccur1) >= 0)), pscooccur1)
pscooccur16S_pruned <-subset_samples(pscooccur16S_pruned, type != "DNA_control")
pscooccur16S_pruned_pa <- transform_sample_counts(pscooccur16S_pruned,function(x)1*(x>0))


#Creating the otu table from the pruned phyloseq object for family:
otutab16S.pa <- otu_table(pscooccur16S_pruned_pa)
otutab16S.x <- otu_table(pscooccur1)
otutab16S.10 <- otu_table(pscooccur1.10)

#cooccur.species.pa <- cooccur(mat = otutab16S.pa,
#                      type = "spp_site")
#saveRDS(cooccur.species.pa, file = "rds/cooccur.species.rds")

#cooccur.species.x <- cooccur(mat = otutab16S.x,
#                           type = "spp_site")
#saveRDS(cooccur.species, file = "rds/cooccur.species.rds")

cooccur.species.10 <- cooccur(mat = otutab16S.10,
                           type = "spp_site")
saveRDS(cooccur.species, file = "rds/cooccur.species.10.rds")
           
##ITS##########################
###importing in the data set
psITS <- readRDS(file = "rds/ITS_crn.ps.deseq-new.Rds")
psITS.10 <- readRDS(file = "rds/ITS_crn.ps.deseq-new10.Rds")

pscooccurITS_pruned <- prune_samples(names(which(sample_sums(pscooccur1) >= 0)), pscooccur1)
pscooccurITS_pruned_pa <- transform_sample_counts(pscooccurITS_pruned,function(x)1*(x>0))

otutabITS.pa <- otu_table(pscooccurITS_pruned_pa) 
otutabITS.10 <- otu_table(psITS.10)

cooccur_ITS.pa <- cooccur(mat = otutabITS.pa,
                       type = "spp_site")
saveRDS(cooccur_ITS.pa, file = "rds/cooccur.ITS.rds")


#cooccur_ITS.10 <- cooccur(mat = otutabITS.10,
#                           type = "spp_site")
#saveRDS(cooccur_ITS, file = "rds/cooccut.ITS.10.rds")

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

saveRDS(cooccur_both, file = "rds/cooccur.both.rds")
