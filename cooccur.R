trial, not working yet. Consult Lucas on settings that we should use.


Starting out:
```{r}
library(cooccur)
library(rPython)
setwd("~/Documents")

pscooccur <- import_biom("16S_otu_table_newjson.biom")
```


Cutting down the information (remove samples with 1000 or less and then making it presence absence):
```{r}
pscooccur_pruned <- prune_samples(names(which(sample_sums(pscooccur) >= 100)),pscooccur)
pscooccur_pruned[pscooccur_pruned > 0] <- 1
```

Running the cooccur function:
```{r}
##Using R python- python.call("filter_otus_from_otu_table.py -i 16S_otu_table_newjson.biom -o 16S_otu_table_newjson_no_singletons.biom -n 2")

cooccur_16S <- cooccur(mat = pscooccur_pruned,
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