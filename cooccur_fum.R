###Import information
pscooccur1 <- readRDS(file = "pscooccur1.rds")
sd <- import_qiime_sample_data("corn_rox.tsv")




otutab1.species <- otu_table(pscooccur1)
taxtab.species <- tax_table(pscooccur1)
pscooccur1_named <- rename_otus(pscooccur1)

sd <- sample_data(pscooccur1_named)
taxtab.species.named <- tax_table(pscooccur1_named)
Fus_taxon <- c("Fungi","Ascomycota","Sordariomycetes","Hypocreales","Nectriaceae","Fusarium","verticillioides")
taxtab2 <- rbind(taxtab.species.named,"Fusarium_vertillioides"=Fus_taxon)
tax_table(pscooccur1_named) <- taxtab2
otutab1.species <- otu_table(pscooccur1_named)
otutab1.speciesfum <- rbind(otutab1.species, "Fusarium_vertillioides"= sd$fum_presence_01)
rownames(otutab1.speciesfum) == rownames(taxtab2)

otutab1.speciesfum <- matrix(otutab1.speciesfum)

rownames(otutab1.speciesfum) == rownames(taxtab2)

pscooccur1_namedfus <- phyloseq(otu_table(otutab1.speciesfum, taxa_are_rows = TRUE),
                                tax_table(taxtab2),
                                sample_data(sd)
                                )

saveRDS(pscooccur1_namedfus, file = "pscooccur1_namedfus.rds")

###Change to presence absence
fus_pa <- transform_sample_counts(pscooccur1_namedfus,function(x)1*(x>0))

###Run the cooccur function (warning, can take like 40 minutes)
cooccur.species <- cooccur(mat = otu_table(fus_pa),
                           type = "spp_site",
                           thresh = TRUE,
                           spp_names = TRUE)

saveRDS(cooccur.species, file = "cooccur.species.rds")

###check to see the interactions with fusarium
fum.cooccur <- pair(cooccur.species, "Fusarium_vertillioides", all = TRUE) 
summary(fum.cooccur)
plot(fum.cooccur)




