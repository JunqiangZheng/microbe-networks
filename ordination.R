#Import data
readRDS(file = "pscooccur1.rds")

pscooccur1_ord <-phyloseq::ordinate(pscooccur1, method = "DCA", distance = "bray")
names(pscooccur1_ord)
names(sample_data(pscooccur1))

plot(1)
col.pal <- wes_palette(n = 6, name = "Moonrise1", type = "continuous")
palette(col.pal)
plot(pscooccur1_ord$rproj[,1:2], col=sample_data(pscooccur1)$type, pch=16)
head(sample_data(pscooccur1))

ordihull(pscooccur1_ord$rproj[,1:2], sample_data(pscooccur1)$type, col=sample_data(pscooccur1)$type)
legend("topright", 
       legend = levels(sample_data(pscooccur1)$type), 
       pch=16,
       col = as.factor(levels(sample_data(pscooccur1)$type)))

#determining statistical signficance of type on distance
tmp <- vegdist(otu_table(pscooccur), method = "bray")

# taxa by type
pscooccur1_bar <- merge_samples(pscooccur1, "type")
pscooccur1_bar <- transform_sample_counts(pscooccur1_bar, function(x) 100 * x / sum(x))