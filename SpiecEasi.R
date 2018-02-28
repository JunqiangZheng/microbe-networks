setwd("C:/Users/DrPig/Desktop/microbe-networks")
library(phyloseq)
library(igraph)
library(SpiecEasi)
library(Matrix)
ps.ITS <- readRDS(file = "rds/ITS_crn.ps.deseq.Rds")
ps.16S <- readRDS(file = "rds/pscooccur1.rds")
ps.both <- readRDS(file = "rds/combined.both_crn.ps.deseq.Rds")

#16S############################

threshold.16S <- 0.1*nsamples(ps.16S)
ez.16S <-filter_taxa(ps.16S,function(x)sum(x>0)>=threshold.16S,TRUE)

spiec.out.16S <- spiec.easi(otu_table(ez.16S), method="mb",icov.select.params=list(rep.num=20, nc=1))
##Error in parallel::mclapply(1:rep.num, function(i) { : 'mc.cores' > 1 is not supported on Windows
spiec.graph.16S <- adj2igraph(spiec.out.16S$refit, vertex.attr=list(name=taxa_names(ez.16S)))
plot_network(spiec.graph.16S, ez.16S, type='taxa', color="Class", label=NULL)

# How many positive vs negative correlations are there?
betaMat.16S <- as.matrix(symBeta(as.matrix(getOptBeta(spiec.out.16S))))
positive.16S <- length(betaMat.16S[betaMat.16S>0])/2
negative.16S <- length(betaMat.16S[betaMat.16S<0])/2
total.16S <- length(betaMat.16S[betaMat.16S!=0])/2 

# For some reason, it's all positive correlations. Also, "lasso" method doesn't yield any correlations

# Try sparcc
sparcc.out.16S <- sparcc(otu_table(ez.16S))
library(Matrix)
# This converts all correlations to positive ones, using an arbitrary correlation cutoff value of 0.3.
# Try to use your current cooccur script to pull out both significant positive and negative correlations
sparcc.graph.16S <- abs(sparcc.out.16S$Cor) >= 0.3
diag(sparcc.graph.16S) <- 0
sparcc.graph.16S <- Matrix(sparcc.graph.16S, sparse=TRUE)
ig.sparcc.16S <- adj2igraph(sparcc.graph.16S, vertex.attr=list(name=taxa_names(ez.16S)))
plot(ig.sparcc.16S, vertex.label=NA, vertex.size = 1, main="sparcc")
plot_network(ig.sparcc.16S, ez.16S, type='taxa', color="Class", label=NULL)



#ITS###############################################################

threshold.ITS <- 0.1*nsamples(ps.ITS)
ez.ITS <-filter_taxa(ps.ITS,function(x)sum(x>0)>=threshold.ITS,TRUE)

spiec.out.ITS <- spiec.easi(otu_table(ez.ITS), method="mb",icov.select.params=list(rep.num=20, nc=1))
##Error in parallel::mclapply(1:rep.num, function(i) { : 'mc.cores' > 1 is not supported on Windows
spiec.graph.ITS <- adj2igraph(spiec.out.ITS$refit, vertex.attr=list(name=taxa_names(ez.ITS)))
plot_network(spiec.graph.ITS, ez.ITS, type='taxa', color="Class", label=NULL)

# How many positive vs negative correlations are there?
betaMat.ITS <- as.matrix(symBeta(as.matrix(getOptBeta(spiec.out.ITS))))
positive.ITS <- length(betaMat.ITS[betaMat.ITS>0])/2
negative.ITS <- length(betaMat.ITS[betaMat.ITS<0])/2
total.ITS <- length(betaMat.ITS[betaMat.ITS!=0])/2 

# For some reason, it's all positive correlations. Also, "lasso" method doesn't yield any correlations

# Try sparcc
sparcc.out.ITS <- sparcc(otu_table(ez.ITS))
# This converts all correlations to positive ones, using an arbitrary correlation cutoff value of 0.3.
# Try to use your current cooccur script to pull out both significant positive and negative correlations
sparcc.graph.ITS <- abs(sparcc.out.ITS$Cor) >= 0.3
diag(sparcc.graph.ITS) <- 0
sparcc.graph.ITS <- Matrix(sparcc.graph.ITS, sparse=TRUE)
ig.sparcc.ITS <- adj2igraph(sparcc.graph.ITS, vertex.attr=list(name=taxa_names(ez.ITS)))
plot(ig.sparcc.ITS, vertex.label=NA, vertex.size = 1, main="sparcc")
plot_network(ig.sparcc.ITS, ez.ITS, type='taxa', color="Class", label=NULL)





#BOTH###########################################
threshold.ITS <- 0.1*nsamples(ps.ITS)
ez.ITS <-filter_taxa(ps.ITS,function(x)sum(x>0)>=threshold.ITS,TRUE)

spiec.out.ITS <- spiec.easi(otu_table(ez.ITS), method="mb",icov.select.params=list(rep.num=20, nc=1))
spiec.graph.ITS <- adj2igraph(spiec.out.ITS$refit, vertex.attr=list(name=taxa_names(ez.ITS)))
plot_network(spiec.graph.ITS, ez.ITS, type='taxa', color="Class", label=NULL)

# How many positive vs negative correlations are there?
betaMat.ITS <- as.matrix(symBeta(as.matrix(getOptBeta(spiec.out.ITS))))
positive.ITS <- length(betaMat.ITS[betaMat.ITS>0])/2
negative.ITS <- length(betaMat.ITS[betaMat.ITS<0])/2
total.ITS <- length(betaMat.ITS[betaMat.ITS!=0])/2 

# For some reason, it's all positive correlations. Also, "lasso" method doesn't yield any correlations

# Try sparcc
sparcc.out.ITS <- sparcc(otu_table(ez.ITS))
# This converts all correlations to positive ones, using an arbitrary correlation cutoff value of 0.3.
# Try to use your current cooccur script to pull out both significant positive and negative correlations
sparcc.graph.ITS <- abs(sparcc.out.ITS$Cor) >= 0.3
diag(sparcc.graph.ITS) <- 0
sparcc.graph.ITS <- Matrix(sparcc.graph.ITS, sparse=TRUE)
ig.sparcc.ITS <- adj2igraph(sparcc.graph.ITS, vertex.attr=list(name=taxa_names(ez.ITS)))
plot(ig.sparcc.ITS, vertex.label=NA, vertex.size = 1, main="sparcc")
plot_network(ig.sparcc.ITS, ez.ITS, type='taxa', color="Class", label=NULL)

