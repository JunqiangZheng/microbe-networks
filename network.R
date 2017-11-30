pscooccur1_namedfus <- readRDS(file = "pscooccur1_namedfus.rds")
cooccur.species <- readRDS(file = "cooccur.species.rds")
cooccur_16S.family <- readRDS(file = "cooccur_16S.family.rds")

cooccur_16S.familyP <- prob.table(cooccur_16S.family)
head(cooccur_16S.familyP, n=2)

cooccur_16S.familyC <- effect.sizes(cooccur_16S.family)

cooccur_16S.familyPC <- cbind(cooccur_16S.familyC, cooccur_16S.familyP)

options(repr.plot.width = 10, repr.plot.height = 6)
plot(sort(cooccur_16S.familyPC$effect), 
     xlab='Family pairs', 
     ylab = 'Cooccurrence effect size',
     main='Wood endophyte cooccurrence effect sizes')

cooccur_16S.familyPC$p_gt_adj <- p.adjust(cooccur_16S.familyPC$p_gt, method = "BH")
strong_cooccur.family <- cooccur_16S.familyPC[cooccur_16S.familyPC$p_gt_adj <= 0.05,]
### does not yield any values 







Species





cooccur.speciesP <- prob.table(cooccur.species)
head(cooccur.speciesP, n=2)

cooccur.speciesC <- effect.sizes(cooccur.species)

cooccur.speciesPC <- cbind(cooccur.speciesC, cooccur.speciesP)

options(repr.plot.width = 10, repr.plot.height = 6)
plot(sort(cooccur.speciesPC$effect), 
     xlab='Family pairs', 
     ylab = 'Cooccurrence effect size',
     main='Wood endophyte cooccurrence effect sizes')

cooccur.speciesPC$p_gt_adj <- p.adjust(cooccur.speciesPC$p_gt, method = "BH")
cooccur.speciesPC$p_lt_adj <- p.adjust(cooccur.speciesPC$p_lt, method = "BH")
strong_cooccur.species <- cooccur.speciesPC[cooccur.speciesPC$p_gt_adj <= 0.05,]
negative_cooccur.species <- cooccur.speciesPC[cooccur.speciesPC$p_lt <= 0.05,] 
head(strong_cooccur.species)
#strong_cooccur.species <- strong_cooccur.species[,-c(12,13)]
#colnames(strong_cooccur.species)[10:11] <- c('OTU_A','OTU_B')

head(strong_cooccur.species, n=2)

options(repr.plot.width = 10,repr.plot.height = 6)
plot(sort(strong_cooccur.species$effect), 
     xlab='Family pairs',
     ylab = 'Cooccurrence effect size',
     main='Cooccurrence effect sizes')

dim(strong_cooccur.family)

pos <- strong_cooccur.species[,c(1, 2)]
neg <- negative_cooccur.species[,c(1,2)]
posneg <- rbind(pos, neg)
cooccur_graph.species <- graph_from_data_frame(posneg, directed=FALSE)
pos.graph <- graph_from_data_frame(pos, directed = FALSE)
neg.graph <- graph_from_data_frame(neg, directed = FALSE)





vcols <- vector(length = length(V(cooccur_graph.species)))
vcols[] <- 'black'
vcols[which(names(V(cooccur_graph.species)) == "Fusarium_vertillioides")] <- "pink"
ecol <- rep("red", ecount(cooccur_graph.species))
ecol[E(cooccur_graph.species) %in% E(pos.graph)] <- "blue"
#ecol[E(cooccur_graph.species) %in% E(neg.graph)] <- "red"
#ecol <- rep("red", ecount(cooccur_graph.species))
#ecol[(cooccur_graph.species) %in% pos.graph] <- "blue"
#ecol[(cooccur_graph.species) %in% neg.graph] <- "red"

plot(cooccur_graph.species, 
     vertex.color = vcols, 
     vertex.size = 10, 
     vertex.label = NA, 
     edge.color = ecol,
     edge.width = 5)
