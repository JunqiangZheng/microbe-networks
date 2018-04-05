##16S####################
#Import data
ps.16S <- readRDS(file = "rds/16S_crn.ps.deseq-new.Rds")

#What kind of metadata can we look at to see if there is a significant difference in the community structure
ps.16S_dis <- phyloseq::distance(ps.16S, method = "bray")
ps.16S_ord<- ordinate(ps.16S, method ="NMDS", ps.16S_dis)

#look into using vegan-adonis to determine statistical significance
ps.16S_df <- as(sample_data(ps.16S),"data.frame")


for (i in c(ps.16S_dis~type,ps.16S_dis~year,ps.16S_dis~grower_ID,ps.16S_dis~fum_presence)){
    #z <- adonis(i, ps.16S_df)
    #x <- as.data.frame(z$aov.tab)[1,]
    #results_16S <- rbind(results, x)
    print(i)
}
print(results_16S)
write.csv(results_16S, "images/ordination_16S.csv")

col.pal <- wes_palette(n = 6, name = "Darjeeling2", type = "continuous")
palette(col.pal)
plot_ordination(ps.16S, ps.16S_ord, color = "type")

##ITS###########
#Import data
ps.ITS <- readRDS(file = "rds/ITS_crn.ps.deseq-new.Rds")

#What kind of metadata can we look at to see if there is a significant difference in the community structure
names(sample_data(ps.ITS))

ps.ITS_dis <- phyloseq::distance(ps.ITS, method = "bray")
ps.ITS_ord1<- ordinate(ps.ITS, method ="NMDS", ps.ITS_dis)

#look into using vegan-adonis to determine statistical significance
ps.ITS_df <- as(sample_data(ps.ITS),"data.frame")
adonis(ps.ITS_dis ~grower_ID*type, ps.ITS_df)

for (i in c(ps.ITS_dis ~type,ps.ITS_dis ~year,ps.ITS_dis ~grower_ID,ps.ITS_dis ~fum_presence)){
  z <- adonis(i, ps.ITS_df)
  x <- as.data.frame(z$aov.tab)[1,]
  results_ITS <- rbind(results, x)
}
print(results_ITS)
write.csv(results_ITS, "images/ordination_ITS.csv")

col.pal <- wes_palette(n = 6, name = "Darjeeling2", type = "continuous")
palette(col.pal)
plot_ordination(ps.ITS, ps.ITS_ord1, color = "type")
