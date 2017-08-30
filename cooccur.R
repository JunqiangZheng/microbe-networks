trial, not working yet. Consult Lucas on settings that we should use.


Starting out:
```{r}
library(cooccur)
library(rPython)
setwd("~/Documents")

pscooccur <- import_biom("16S_otu_table_newjson.biom")
```


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
                      type = "spp_site",
                      thresh = TRUE,
                      spp_names = TRUE
                      #true_rand_classifier = "???",
                      #prob = "???",
                      #site_mask = "???"
                      #only_effects = TRUE,
                      #eff_standard = TRUE
)

summary(cooccur_16S)
plot(cooccur_16S)
```