#Import data
ps.both<- readRDS(file = "rds/combined.both_crn.ps.Rds")

#What kind of metadata can we look at to see if there is a significant difference in the community structure
names(sample_data(ps.both))

ps.both_dis <- phyloseq::distance(ps.both, method = "bray")
ps.both_ord1<- ordinate(ps.both, method ="NMDS", ps.both_dis)

#look into using vegan-adonis to determine statistical significance
df.both <- as(sample_data(ps.both),"data.frame")
adonis(ps.both_dis ~type, df.both)



#determining statistical signficance of type on distance
tmp.both <- vegdist(otu_table(ps.both), method = "bray")

# taxa by type
ps.both_bar <- merge_samples(ps.both, "type")
ps.both_bar <- transform_sample_counts(ps.both_bar, function(x) 100 * x / sum(x))



col.pal <- wes_palette(n = 6, name = "Darjeeling2", type = "continuous")
palette(col.pal)
#plot(ps.both_ord$rproj[,1:2], col=sample_data(ps.both)$type, pch=16)
#ordihull(ps.both_ord1$rproj[,1:2], sample_data(ps.both)$type, col=sample_data(ps.both)$type)
#legend("topright", 
#     legend = levels(sample_data(ps.both)$type), 
#   pch=16,
# col = as.factor(levels(sample_data(ps.both)$type)))



plot_ordination(ps.both, ps.both_ord1, color = "type")


