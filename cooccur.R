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
```
                                                     
                                                     
Code wasnt syncing from rstudio correctly so I pasted the code I had into the bottom of the most recent thing uploaded to git...
                                                     
trial, not working yet. Consult Lucas on settings that we should use.


Starting out, loading all libraries and setting the work domain, importing the otu table and using that function to label columns as "species", "genus", "family"...:
```{r}
library(cooccur)
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


Cutting down the information... 
Grouping it by family and putting in the taxonomic names:
```{r}
pscooccur.genus <- tax_glom(pscooccur, taxrank = "Genus")
otutab1.genus <- otu_table(pscooccur.genus)
taxtab.genus <- tax_table(pscooccur.genus)
rownames(otutab1.genus) <- taxtab.genus[,6]

pscooccur.family <- tax_glom(pscooccur.genus, taxrank = "Family")
otutab1.family <- otu_table(pscooccur.family)
taxtab.family <- tax_table(pscooccur.family)
rownames(otutab1.family) <- taxtab.family[,5]
```


Cutting down controls, (add remove samples with 1000 or less), and then making it presence absence for family:
```{r}
pscooccur_pruned.family <- prune_samples(names(which(sample_sums(otutab1.family) >= 0)),otutab1.family)
pscooccur_pruned.family <-subset_samples(pscooccur_pruned.family, Tissue != "DNA_control") # sam_data is empty, has it already been removed?
pscooccur_pruned.family <- prune_taxa(taxa_sums(pscooccur_pruned.family) > 1000, pscooccur_pruned.family)
pscooccur_pruned_pa.family <- transform_sample_counts(pscooccur_pruned.family,function(x)1*(x>0))
```

Creating the otu table from the pruned phyloseq object for family:
```{r}
otutab.family <- otu_table(pscooccur_pruned_pa.family) # this is a data frame
#otutab <- as.matrix(otu_table(pscooccur)) # this is a matrix
```


Running the cooccur function for family:
```{r}
pscooccur_pruned.family <- prune_samples(names(which(sample_sums(pscooccur.family) >= 0)),pscooccur)
pscooccur_pruned.family <-subset_samples(pscooccur_pruned.family, Tissue != "DNA_control")
pscooccur_pruned_pa.family <- transform_sample_counts(pscooccur_pruned.family,function(x)1*(x>0))
```

Running the cooccur function - gives me error that the argument is length 0:
```{r}
cooccur_16S.family <- cooccur(mat = otutab.family,
                      type = "spp_site",
                      thresh = TRUE,
                      spp_names = TRUE,
                      true_rand_classifier = .1, # default value
                      prob = "comb" #combinatory (comb) or hypergeometric (hyper)
                      #site_mask = "???"
                      #only_effects = TRUE,
                      #eff_standard = TRUE
)

summary(cooccur_16S.family)
plot(cooccur_16S.family)
```


Doing the same process for Genus:
```{r}
pscooccur.genus_pruned <- prune_samples(names(which(sample_sums(pscooccur.genus) >= 0)),pscooccur.genus)
pscooccur.genus_pruned <-subset_samples(pscooccur.genus_pruned, Tissue != "DNA_control") # sam_data is empty, has it already been removed?
pscooccur.genus_pruned <- prune_taxa(taxa_sums(pscooccur.genus_pruned) > 600, pscooccur.genus_pruned)
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
 
