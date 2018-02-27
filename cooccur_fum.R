###Import information
pscooccur1 <- readRDS(file = "rds/pscooccur1.rds")
sd <- sample_data(pscooccur1)




otutab1.species <- otu_table(pscooccur1)
taxtab.species <- tax_table(pscooccur1)


taxtab.species.named <- tax_table(pscooccur1)
Fus_taxon <- c("Fungi","Ascomycota","Sordariomycetes","Hypocreales","Nectriaceae","Fusarium","verticillioides")
taxtab2 <- rbind(taxtab.species.named,"Fusarium_vertillioides"=Fus_taxon)
tax_table(pscooccur1) <- taxtab2
otutab1.species <- otu_table(pscooccur1)
otutab1.speciesfum <- matrix(otutab1.species)

rownames(otutab1.speciesfum) == rownames(taxtab2)



rownames(otutab1.speciesfum) == rownames(taxtab2)

pscooccur1_namedfus <- phyloseq(otu_table(otutab1.species, taxa_are_rows = TRUE),
                                tax_table(taxtab2),
                                sample_data(sd)
                                )

saveRDS(pscooccur1_namedfus, file = "rds/pscooccur1_namedfus.rds")

###Change to presence absence
fus_pa <- transform_sample_counts(pscooccur1_namedfus,function(x)1*(x>0))

###Run the cooccur function (warning, can take like 40 minutes)
cooccur.species <- cooccur(mat = otu_table(fus_pa),
                           type = "spp_site",
                           thresh = TRUE,
                           spp_names = TRUE)

saveRDS(cooccur.species, file = "rds/cooccur.species.rds")

###check to see the interactions with fusarium
fum.cooccur <- pair(cooccur.species, "Fusarium_vertillioides", all = TRUE) 
summary(fum.cooccur, "Fusarium_vertillioides")
plot(fum.cooccur)



#####The only statistically significant one after adjusted p value - Limnobacter thioxidans


##Now looking to see if there are any that cooccur with fus amount
species_otu_pa <- transform_sample_counts(otutab1.species,function(x)1*(x>0))
otutab1.speciesfumamount <- rbind(otutab1.species, "Fusarium_vertillioides"= sd$log_avg_copies_fum)

pscooccur1_namedfusabundance <- phyloseq(otu_table(otutab1.speciesfumamount, taxa_are_rows = TRUE),
                                tax_table(taxtab2),
                                sample_data(sd)
)

cooccur.species.fusabundance <- cooccur(mat = otu_table(pscooccur1_namedfusabundance),
                           type = "spp_site",
                           thresh = TRUE,
                           spp_names = TRUE)
###############################

pscooccur1 <- readRDS(file = "rds/pscooccur1.rds")
sd <- sample_data(pscooccur1)




otutab1.species <- otu_table(pscooccur1)
taxtab.species <- tax_table(pscooccur1)


taxtab.species.named <- tax_table(pscooccur1)
Fus_taxon <- c("Fungi","Ascomycota","Sordariomycetes","Hypocreales","Nectriaceae","Fusarium","verticillioides")
taxtab2 <- rbind(taxtab.species.named,"Fusarium_vertillioides"=Fus_taxon)
tax_table(pscooccur1) <- taxtab2
otutab1.species <- otu_table(pscooccur1)
otutab1.speciesfum <- matrix(otutab1.species)

rownames(otutab1.speciesfum) == rownames(taxtab2)



rownames(otutab1.speciesfum) == rownames(taxtab2)

pscooccur1_namedfus <- phyloseq(otu_table(otutab1.species, taxa_are_rows = TRUE),
                                tax_table(taxtab2),
                                sample_data(sd)
)

saveRDS(pscooccur1_namedfus, file = "rds/pscooccur1_namedfus.rds")

###Change to presence absence
fus_pa <- transform_sample_counts(pscooccur1_namedfus,function(x)1*(x>0))assifier = .1, 
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
##



