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
strong_cooccur.species <- cooccur.speciesPC[cooccur.speciesPC$p_gt_adj <= 0.05,]
head(strong_cooccur.species)
strong_cooccur.species <- strong_cooccur.species[,-c(12,13)]
colnames(strong_cooccur.species)[10:11] <- c('OTU_A','OTU_B')

head(strong_cooccur.family, n=2)

options(repr.plot.width = 10,repr.plot.height = 6)
plot(sort(strong_cooccur.family$effect), 
     xlab='Family pairs',
     ylab = 'Cooccurrence effect size',
     main='Cooccurrence effect sizes')

dim(strong_cooccur.family)

aa <- strong_cooccur.family[,c(10,11,12,13)]
cooccur_graph.family <- graph_from_data_frame(aa, directed=FALSE)

vcols <- vector(length = length(V(cooccur_graph.family)))
vcols[] <- 'brown'
vcols[V(cooccur_graph.family) %in% aa] <- "purple"
>>>>>>> 708ffb163bbe7b9fd27f27ce21715f3dbae0c0fd
