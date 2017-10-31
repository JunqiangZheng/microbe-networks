#Import data
readRDS(file = "pscooccur1.rds")

#What is the difference between methods... NMDS vs DCA


names(pscooccur1_ord)
names(sample_data(pscooccur1))

col.pal <- wes_palette(n = 6, name = "Darjeeling2", type = "continuous")
palette(col.pal)
#plot(pscooccur1_ord$rproj[,1:2], col=sample_data(pscooccur1)$type, pch=16)
#ordihull(pscooccur1_ord$rproj[,1:2], sample_data(pscooccur1)$type, col=sample_data(pscooccur1)$type)
#legend("topright", 
#      legend = levels(sample_data(pscooccur1)$type), 
#     pch=16,
#    col = as.factor(levels(sample_data(pscooccur1)$type)))

#pscooccur1_ord1<- ordinate(pscooccur1, method ="NMDS", pscooccur1_dis)

plot_ordination(pscooccur1, pscooccur1_ord, color = "type")
head(sample_data(pscooccur1))


#look into using vegan-adonis to determine statistical significance
pscooccur1_dis <- phyloseq::distance(pscooccur1, method = "bray")
df <- as(sample_data(pscooccur1),"data.frame")
adonis(pscooccur1_dis ~type, df)



#determining statistical signficance of type on distance
tmp <- vegdist(otu_table(pscooccur1), method = "bray")

# taxa by type
pscooccur1_bar <- merge_samples(pscooccur1, "type")
pscooccur1_bar <- transform_sample_counts(pscooccur1_bar, function(x) 100 * x / sum(x))