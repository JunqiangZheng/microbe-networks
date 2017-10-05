                                                   
Starting out, loading all libraries and setting the work domain, importing the otu table and using that function to label columns as "species", "genus", "family"...:
```{r}
library(cooccur)
library(phyloseq)
library(ggplot2)
library("DESeq2")
library(vegan)
library(psych)
library(biom)
library(wesanderson)
setwd("~/Documents")

pscooccur <- import_biom("16S_otu_table_newjson.biom",
                         parseFunction=parse_taxonomy_greengenes)
```

Adding in DESeq function:
```{r}
#Roo's DESeq script
DESeq_varstab <- function(phyloseq, design) {
  # phyloseq = the input phylose object that you want to get DESeq transformed counts for
  # design_variable = the design for the conversion to the DESeq object. must be in the form "as a function of", for example "~Host_Genus", must be a variable in the phyloseq object
  # Set variables to NULL
  deseq.vst = NULL
  geo_Means = NULL
  phyloseq.DESeq = NULL
  # Convert to a DESeq object
  deseq = phyloseq_to_deseq2(phyloseq, design)
  # calculate geometric means prior to estimate size factors
  gm_mean = function(x, na.rm=TRUE){
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  }
  geo_Means = apply(counts(deseq), 1, gm_mean)
  # Check to see if any columns (samples) don't have any OTUs in them:
  if(sum(colSums(counts(deseq)) == 0) == 0) { # if all samples have taxa, go on
    # Now we step through the size factors, dispersions, and varience stabilization:
    deseq = estimateSizeFactors(deseq, geoMeans = geo_Means)
    deseq = estimateDispersions(deseq) # long step
    deseq.vst = getVarianceStabilizedData(deseq)
    # replace negatives with zeros
    deseq.vst[deseq.vst <0] <- 0
    # add the varience stabilized otu numbers into the dataset:
    otu_table(phyloseq) <- otu_table(deseq.vst, taxa_are_rows = TRUE)
    # create a new object for the varience stabalized set
    phyloseq -> phyloseq.DESeq
    # And, filter any taxa that became 0s all the way across
    phyloseq.DESeq = filter_taxa(phyloseq.DESeq, function(x) sum(x) > 0.1, T)
    # return the new phyloseq object
    return(phyloseq.DESeq)
  } # end of IF loop 
  else {return("Error: your phyloseq object has samples with no taxa present.")}
} # end function
```

Accounting for error within the sample:
```{r}
pscooccur1 <- prune_samples(names(which(sample_sums(pscooccur) >= 0)),pscooccur)
pscooccur1 <- prune_taxa(taxa_sums(pscooccur1) > 11, pscooccur1)
sd <- import_qiime_sample_data("mapping_allsamples_simple.csv")
pscooccur1 <- merge_phyloseq(pscooccur1, sd)
pscooccur1 <- DESeq_varstab(pscooccur1, ~1)
```

Ordination:
```{r}
pscooccur1_ord <-phyloseq::ordinate(pscooccur1, method = "DCA", distance = "bray")
names(pscooccur1_ord)

col.pal <- wes_palette(n = 6, name = "Darjeeling2", type = "continuous")
palette(col.pal)
plot(pscooccur1_ord$rproj[,1:2], col=sample_data(pscooccur1)$Tissue, pch=16)
head(sample_data(pscooccur1))

ordihull(pscooccur1_ord$rproj[,1:2], sample_data(pscooccur1)$Tissue, col=sample_data(pscooccur1)$Tissue)
legend("topright", 
       legend = levels(sample_data(pscooccur1)$Tissue), 
       pch=16,
       col = as.factor(levels(sample_data(pscooccur1)$Tissue)))

# taxa by tissue
pscooccur1_bar <- merge_samples(pscooccur1, "Tissue")
pscooccur1_bar <- transform_sample_counts(pscooccur1_bar, function(x) 100 * x / sum(x))
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
pscooccur_pruned.family <- prune_taxa(taxa_sums(pscooccur_pruned.family) > 1000, pscooccur_pruned.family)
pscooccur_pruned_pa.family <- transform_sample_counts(pscooccur_pruned.family,function(x)1*(x>0))
```

Creating the otu table from the pruned phyloseq object for family:
```{r}
otutab.family <- otu_table(pscooccur_pruned_pa.family) # this is a data frame
#otutab <- as.matrix(otu_table(pscooccur)) # this is a matrix
```

<<<<<<< HEAD
=======

Running the cooccur function for family:
```{r}
pscooccur_pruned.family <- prune_samples(names(which(sample_sums(pscooccur.family) >= 0)),pscooccur)
pscooccur_pruned.family <-subset_samples(pscooccur_pruned.family, Tissue != "DNA_control")
pscooccur_pruned_pa.family <- transform_sample_counts(pscooccur_pruned.family,function(x)1*(x>0))
```
>>>>>>> a9e2b279086993839bbef1cc5b3273bed70c16dd

Running the cooccur function for family:
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
 
