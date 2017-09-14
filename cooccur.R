trial, not working yet. Consult Lucas on settings that we should use.


Starting out, loading all libraries and setting the work domain, importing the otu table and using that function to label columns as "species", "genus", "family"...:
```{r}
library(cooccur)
library(rPython)
library(phyloseq)
library(ggplot2)
library("DESeq2")
library(vegan)
library(psych)
library(biom)
setwd("~/Documents")

pscooccur <- import_biom("16S_otu_table_newjson.biom",
                         parseFunction=parse_taxonomy_greengenes)
```


<<<<<<< HEAD
Cutting down the information... 
Grouping it by family:
```{r}
pscooccur.genus <- tax_glom(pscooccur, taxrank = "Genus")
pscooccur.family <- tax_glom(pscooccur.genus, taxrank = "Family")
```

Cutting down controls, (add remove samples with 1000 or less), and then making it presence absence for family:
```{r}
pscooccur_pruned <- prune_samples(names(which(sample_sums(pscooccur.family) >= 0)),pscooccur.family)
pscooccur_pruned <-subset_samples(pscooccur_pruned, Tissue != "DNA_control") # sam_data is empty, has it already been removed?
pscooccur_pruned <- prune_taxa(taxa_sums(pscooccur_pruned) > 1000, pscooccur_pruned)
pscooccur_pruned_pa <- transform_sample_counts(pscooccur_pruned,function(x)1*(x>0))
```

Creating the otu table from the pruned phyloseq object for family:
```{r}

otutab <- otu_table(pscooccur_pruned_pa) # this is a data frame
#otutab <- as.matrix(otu_table(pscooccur)) # this is a matrix
```

Running the cooccur function for just presence absence for family:
```{r}
cooccur_16S <- cooccur(mat = otutab,
=======
Cutting down the information -Grouping it by family, remove samples with 1000 or less (if its just looking at family abundance should I increase this number?), and then making it presence absence:
```{r}
pscooccur.genus <- tax_glom(pscooccur, taxrank = "Genus")
pscooccur.family <- tax_glom(pscooccur.genus, taxrank = "Family")
pscooccur_pruned <- prune_samples(names(which(sample_sums(pscooccur.family) >= 0)),pscooccur)
pscooccur_pruned <-subset_samples(pscooccur_pruned, Tissue != "DNA_control")
pscooccur_pruned_pa <- transform_sample_counts(pscooccur_pruned,function(x)1*(x>0))
```

Running the cooccur function - gives me error that the argument is length 0:
```{r}
cooccur_16S <- cooccur(mat = pscooccur_pruned_pa,
>>>>>>> 1ed463a2be6018250585ee260952ffbf979ee170
                      type = "spp_site",
                      thresh = TRUE,
                      spp_names = TRUE,
                      true_rand_classifier = .1, # default value
                      prob = "comb" #combinatory (comb) or hypergeometric (hyper)
                      #site_mask = "???"
                      #only_effects = TRUE,
                      #eff_standard = TRUE
)

summary(cooccur_16S)
plot(cooccur_16S)
```


Doing the same process for Genus:
```{r}
pscooccur.genus_pruned <- prune_samples(names(which(sample_sums(pscooccur.genus) >= 0)),pscooccur.genus)
pscooccur.genus_pruned <-subset_samples(pscooccur.genus_pruned, Tissue != "DNA_control") # sam_data is empty, has it already been removed?
pscooccur.genus_pruned <- prune_taxa(taxa_sums(pscooccur.genus_pruned) > 500, pscooccur.genus_pruned)
pscooccur.genus_pruned_pa <- transform_sample_counts(pscooccur.genus_pruned,function(x)1*(x>0))

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
plot(cooccur.genus_16S)
