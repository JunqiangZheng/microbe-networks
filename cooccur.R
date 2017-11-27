###importing in the data set
readRDS(file = "pscooccur1.rds")


###Cutting down the information... Grouping it by family and putting in the taxonomic names:
pscooccur.genus <- tax_glom(pscooccur1_namedfus, taxrank = "Genus")
otutab1.genus <- otu_table(pscooccur.genus)
taxtab.genus <- tax_table(pscooccur.genus)
rownames(otutab1.genus) <- taxtab.genus[,6]

pscooccur.family <- tax_glom(pscooccur.genus, taxrank = "Family")
otutab1.family <- otu_table(pscooccur.family)
taxtab.family <- tax_table(pscooccur.family)
rownames(otutab1.family) <- taxtab.family[,5]




#Running the cooccur function for family:
pscooccur_pruned.family <- prune_samples(names(which(sample_sums(pscooccur.family) >= 0)),pscooccur.family)
pscooccur_pruned.family <-subset_samples(pscooccur_pruned.family, Type != "DNA_control")
pscooccur_pruned.family <- prune_taxa(taxa_sums(pscooccur_pruned.family) > 1000, pscooccur_pruned.family)
pscooccur_pruned_pa.family <- transform_sample_counts(pscooccur_pruned.family,function(x)1*(x>0))


#Creating the otu table from the pruned phyloseq object for family:
otutab.family <- otu_table(pscooccur_pruned_pa.family) # this is a data frame
#otutab <- as.matrix(otu_table(pscooccur)) 
#this is a matrix



#Running the cooccur function for family:
cooccur_16S.family <- cooccur(mat = otutab.family,
                      type = "spp_site"
                      #spp_names = TRUE,
                      #true_rand_classifier = .1, 
                      # default value prob = "comb" #combinatory (comb) or hypergeometric (hyper)
                      #site_mask = "???"
                      #only_effects = TRUE,
                      #eff_standard = TRUE
)

summary(cooccur_16S.family)
plot(cooccur_16S.family)

saveRDS(cooccur_16S.family, file = "cooccur_16S.family.rds")

#Doing the same process for Genus:
#pscooccur.genus_pruned <- prune_taxa(taxa_sums(pscooccur.genus) > 600, pscooccur.genus)
pscooccur.genus_pruned_pa <- transform_sample_counts(pscooccur.genus,function(x)1*(x>0))

otutab.genus <- otu_table(pscooccur.genus_pruned_pa)

cooccur.genus_16S <- cooccur(mat = otutab.genus,
                      type = "spp_site",
                      thresh = TRUE,
                      spp_names = TRUE,
                      true_rand_classifier = .1, # default value
                      prob = "comb" #combinatory (comb) or hypergeometric (hyper)
                      #site_mask = "???"
                      #only_effects = TRUE,
                      #eff_standard = TRUE
)

summary(cooccur.genus_16S)
plot.cooccur(cooccur.genus_16S)