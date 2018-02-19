library(SpiecEasi)

#16S
############################

threshold <- 0.1*nsamples(ps.both)
ps <-filter_taxa(ps.both,function(x)sum(x>0)>=threshold,TRUE)

spiec.out=spiec.easi(ps, method="mb",icov.select.params=list(rep.num=20, nc=4))
spiec.graph=adj2igraph(spiec.out$refit, vertex.attr=list(name=taxa_names(ps)))
plot_network(spiec.graph, ps, type='taxa', color="Class", label=NULL)

# How many positive vs negative correlations are there?
betaMat=as.matrix(symBeta(as.matrix(getOptBeta(spiec.out))))
positive=length(betaMat[betaMat>0])/2
negative=length(betaMat[betaMat<0])/2
total=length(betaMat[betaMat!=0])/2 

# For some reason, it's all positive correlations. Also, "lasso" method doesn't yield any correlations

# Try sparcc
sparcc.out <- sparcc(otu_table(ps))
library(Matrix)
# This converts all correlations to positive ones, using an arbitrary correlation cutoff value of 0.3.
# Try to use your current cooccur script to pull out both significant positive and negative correlations
sparcc.graph <- abs(sparcc.out$Cor) >= 0.3
diag(sparcc.graph) <- 0
sparcc.graph <- Matrix(sparcc.graph, sparse=TRUE)
ig.sparcc <- adj2igraph(sparcc.graph, vertex.attr=list(name=taxa_names(ps)))
plot(ig.sparcc, vertex.label=NA, vertex.size = 1, main="sparcc")
plot_network(ig.sparcc, ps, type='taxa', color="Class", label=NULL)





