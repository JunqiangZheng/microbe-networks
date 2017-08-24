trial, not working yet. Consult Lucas on settings that we should use.

```{r}
library(cooccur)
library(rPython)
setwd("~/Documents")

import_biom("16S_otu_table_newjson.biom")
python.call("filter_otus_from_otu_table.py -i 16S_otu_table_newjson.biom -o 16S_otu_table_newjson_no_singletons.biom -n 2")

cooccur_16S <- cooccur(mat = "16S_otu_table_newjson_no_singletons.biom",
                      type = "spp_site",
                      thresh = TRUE,
                      #spp_names = TRUE
                      #true_rand_classifier = "???",
                      #prob = "???",
                      #site_mask = "???"
                      only_effects = TRUE,
                      eff_standard = TRUE
)

summary(cooccur_16S)
plot(cooccur_16S)
```