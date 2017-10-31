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






##Check to see that Fusarium was correctly added into the species list (for some reason it adds it withough a label?? row 505)
"fum_presence_01" %in% rownames(otutab1.speciesfum)

taxtab.species.namedfum <- rbind(taxtab.species.named, "fum_presence_01")

rownames(otutab1.speciesfum)[505] <- "fum_presence_01"



###Script written by Lucas to rename the species
rename_otus <- function(ps, max.sp=2){
  taxtab <- tax_table(ps)
  taxa_list <- vector()
  for(row in 1:nrow(taxtab)){
    taxon <-taxtab[row,]
    NAs <- sum(is.na(taxon))
    # Genus_species nomenclature for OTUs IDed to species
    if(NAs == 0){
      name <- paste(taxon[,6],taxon[,7], sep="_")
      # Restrict ambiguous species calls
      if(length(unlist(strsplit(name,"/"))) > max.sp){
        NAs <- 1
      }
    }
    # If OTU is named to Genus, use "[Genus]_sp" convention
    if(NAs == 1){
      name <- paste(taxon[,length(taxon)-NAs],"sp", sep="_")
    }
    # If OTU identified to higher taxonomic level, just use name of that level
    if(NAs > 1){
      name <- paste(taxon[,length(taxon)-NAs])
    }
    taxa_list <- c(taxa_list,name)
  }
  # Rename duplicate taxa in taxa_list, by adding numbers to the end
  for(taxon in 1:length(taxa_list)){
    matches <- taxa_list[taxon] == taxa_list
    nmatches <- sum(matches)
    if(nmatches > 1){
      taxa_list[matches] <- paste(taxa_list[matches],1:(nmatches),sep = "_") 
    }
  }
  taxa_names(ps) <- taxa_list
  return(ps)
}
###End of scipt written by Lucas


