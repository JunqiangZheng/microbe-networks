###importing in the data set
readRDS(file = "pscooccurold1.rds")

###Cutting down the information... Grouping it by family and putting in the taxonomic names:
pscooccurold.genus <- tax_glom(pscooccurold1, taxrank = "Genus")
otutab1old.genus <- otu_table(pscooccurold.genus)
taxtabold.genus <- tax_table(pscooccurold.genus)
rownames(otutab1old.genus) <- taxtabold.genus[,6]

pscooccurold.family <- tax_glom(pscooccurold.genus, taxrank = "Family")
otutab1old.family <- otu_table(pscooccurold.family)
taxtabold.family <- tax_table(pscooccurold.family)
rownames(otutab1old.family) <- taxtabold.family[,5]




#Running the cooccur function for family:
pscooccur_prunedold.family <- prune_samples(names(which(sample_sums(pscooccurold.family) >= 0)),pscooccurold.family)
pscooccur_prunedold.family <-subset_samples(pscooccur_prunedold.family, Type != "DNA_control")
pscooccur_prunedold.family <- prune_taxa(taxa_sums(pscooccur_prunedold.family) > 1000, pscooccur_prunedold.family)
pscooccur_pruned_paold.family <- transform_sample_counts(pscooccur_prunedold.family,function(x)1*(x>0))


#Creating the otu table from the pruned phyloseq object for family:
otutabold.family <- otu_table(pscooccur_pruned_paold.family) # this is a data frame
#otutab <- as.matrix(otu_table(pscooccur)) 
#this is a matrix



#Running the cooccur function for family:
cooccur_16Sold.family <- cooccur(mat = otutabold.family,
                              type = "spp_site",
                              thresh = TRUE,
                              spp_names = TRUE,
                              true_rand_classifier = .1, # default value
                              prob = "comb" #combinatory (comb) or hypergeometric (hyper)
                              #site_mask = "???"
                              #only_effects = TRUE,
                              #eff_standard = TRUE
)

summary(cooccur_16Sold.family)
plot(cooccur_16Sold.family)


#Doing the same process for Genus:
pscooccurold.genus_pruned_pa <- transform_sample_counts(pscooccurold.genus,function(x)1*(x>0))

otutabold.genus <- otu_table(pscooccurold.genus_pruned_pa)

cooccurold.genus_16S <- cooccur(mat = otutabold.genus,
                             type = "spp_site",
                             thresh = TRUE,
                             spp_names = TRUE,
                             true_rand_classifier = .1, # default value
                             prob = "comb" #combinatory (comb) or hypergeometric (hyper)
                             #site_mask = "???"
                             #only_effects = TRUE,
                             #eff_standard = TRUE
)

summary(cooccurold.genus_16S)
plot(cooccurold.genus_16S)