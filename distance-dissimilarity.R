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
