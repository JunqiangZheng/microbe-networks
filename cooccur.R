##trial, not working yet. Consult Lucas on settings that we should use.
library(cooccur)
data(ps_16S)
cooccur_16S <- cooccur(mat = "ps_16S",
                      type = "???",
                      thresh = TRUE,
                      spp_names = TRUE,
                      true_rand_classifier = "???",
                      prob = "???",
                      site_mask = "???")

summary(cooccur_16S)
plot(cooccur_16S)