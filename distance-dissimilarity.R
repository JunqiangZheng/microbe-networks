readRDS(file = "pscooccur1.rds")

##community dissimilarity
pscooccur1_dis <- phyloseq::distance(pscooccur1, method = "bray")
df <- as(sample_data(pscooccur1),"data.frame")

##https://www.rdocumentation.org/packages/ape/versions/4.1/topics/mantel.test
##this is not avaliable for my version of R

##https://stats.idre.ucla.edu/r/faq/how-can-i-perform-a-mantel-test-in-r/

##Creating the dist object for the lonlat
lonlat <- dist(cbind(df$longitude, df$latitude))


mantel.rtest(pscooccur1_dis, lonlat, nrepet = 9999)
#seem to be unrelated to geographical distances










#################################### Doing it for the OTU clustering set
pscooccurold1_dis <- phyloseq::distance(pscooccurold1, method = "bray")
dfold <- as(sample_data(pscooccurold1),"data.frame")

##https://www.rdocumentation.org/packages/ape/versions/4.1/topics/mantel.test
##this is not avaliable for my version of R

##https://stats.idre.ucla.edu/r/faq/how-can-i-perform-a-mantel-test-in-r/

##Creating the dist object for the lonlat
lonlatold <- dist(cbind(dfold$longitude, dfold$latitude))


mantel.rtest(pscooccurold1_dis, lonlatold, nrepet = 9999)
#seem to be unrelated to geographical distances
