#Import data
ps.ITS <- readRDS(file = "rds/ITS_crn.ps.deseq.Rds")

#What kind of metadata can we look at to see if there is a significant difference in the community structure
names(sample_data(ps.ITS))

ps.ITS_dis <- phyloseq::distance(ps.ITS, method = "bray")
ps.ITS_ord1<- ordinate(ps.ITS, method ="NMDS", ps.ITS_dis)

#look into using vegan-adonis to determine statistical significance
df.ITS <- as(sample_data(ps.ITS),"data.frame")
adonis(ps.ITS_dis ~type, df.ITS)



#determining statistical signficance of type on distance
tmp.ITS <- vegdist(otu_table(ps.ITS), method = "bray")

# taxa by type
ps.ITS_bar <- merge_samples(ps.ITS, "type")
ps.ITS_bar <- transform_sample_counts(ps.ITS_bar, function(x) 100 * x / sum(x))



col.pal <- wes_palette(n = 6, name = "Darjeeling2", type = "continuous")
palette(col.pal)
#plot(ps.ITS_ord$rproj[,1:2], col=sample_data(ps.ITS)$type, pch=16)
#ordihull(ps.ITS_ord1$rproj[,1:2], sample_data(ps.ITS)$type, col=sample_data(ps.ITS)$type)
#legend("topright", 
#     legend = levels(sample_data(ps.ITS)$type), 
#   pch=16,
# col = as.factor(levels(sample_data(ps.ITS)$type)))



plot_ordination(ps.ITS, ps.ITS_ord1, color = "type")


