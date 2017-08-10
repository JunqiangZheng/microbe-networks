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