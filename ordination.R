#Import data
pscooccur1 <- readRDS(file = "rds/pscooccur1.rds")

#What kind of metadata can we look at to see if there is a significant difference in the community structure
names(sample_data(pscooccur1))

popcorn <- sample_data(pscooccur1)[which(sample_data(pscooccur1)$type !='Popcorn') ] 

pscooccur1_dis <- phyloseq::distance(pscooccur1, method = "bray")
pscooccur1_ord1<- ordinate(pscooccur1, method ="NMDS", pscooccur1_dis)

#look into using vegan-adonis to determine statistical significance
df <- as(sample_data(pscooccur1),"data.frame")
adonis(pscooccur1_dis ~type, df)



#determining statistical signficance of type on distance
tmp <- vegdist(otu_table(pscooccur1), method = "bray")

# taxa by type
pscooccur1_bar <- merge_samples(pscooccur1, "type")
pscooccur1_bar <- transform_sample_counts(pscooccur1_bar, function(x) 100 * x / sum(x))



col.pal <- wes_palette(n = 6, name = "Darjeeling2", type = "continuous")
palette(col.pal)
#plot(pscooccur1_ord$rproj[,1:2], col=sample_data(pscooccur1)$type, pch=16)
#ordihull(pscooccur1_ord1$rproj[,1:2], sample_data(pscooccur1)$type, col=sample_data(pscooccur1)$type)
#legend("topright", 
 #     legend = levels(sample_data(pscooccur1)$type), 
  #   pch=16,
   # col = as.factor(levels(sample_data(pscooccur1)$type)))



plot_ordination(pscooccur1, pscooccur1_ord1, color = "type")


#####
####
#####
ps16S <- readRDS(file = "rds/16S_crn.ps.Rds")

#What kind of metadata can we look at to see if there is a significant difference in the community structure
names(sample_data(ps16S))

ps16S_dis <- phyloseq::distance(ps16S, method = "bray")
ps16S_ord1<- ordinate(ps16S, method ="NMDS", ps16S_dis)

#look into using vegan-adonis to determine statistical significance
df.ps16S <- as(sample_data(ps16S),"data.frame")
adonis(ps16S_dis ~type, df.ps16S)



#determining statistical signficance of type on distance
tmp.ps16S <- vegdist(otu_table(ps16S), method = "bray")


col.pal <- wes_palette(n = 6, name = "Darjeeling2", type = "continuous")
palette(col.pal)
#plot(pscooccur1_ord$rproj[,1:2], col=sample_data(pscooccur1)$type, pch=16)
#ordihull(pscooccur1_ord1$rproj[,1:2], sample_data(pscooccur1)$type, col=sample_data(pscooccur1)$type)
#legend("topright", 
#     legend = levels(sample_data(pscooccur1)$type), 
#   pch=16,
# col = as.factor(levels(sample_data(pscooccur1)$type)))



plot_ordination(ps16S, ps16S_ord1, color = "type")


