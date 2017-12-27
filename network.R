cooccur.species <- readRDS(file = "cooccur.species.rds")

cooccur.speciesP <- prob.table(cooccur.species)
head(cooccur.speciesP, n=2)

cooccur.speciesC <- effect.sizes(cooccur.species)

cooccur.speciesPC <- cbind(cooccur.speciesC, cooccur.speciesP)

options(repr.plot.width = 10, repr.plot.height = 6)
plot(sort(cooccur.speciesPC$effect), 
     xlab='Species pairs', 
     ylab = 'Cooccurrence effect size',
     main='Endophyte cooccurrence effect sizes')

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
     xlab='Species pairs',
     ylab = 'Cooccurrence effect size',
     main='Cooccurrence effect sizes')

dim(strong_cooccur.species)

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

vertex <- names(V(cooccur_graph.species))
taxa_names(pscooccur1_namedfus)
pscooccur1_names_fus <- rename_otus(pscooccur1_namedfus)
ps_vertices <- prune_taxa(taxa_names(pscooccur1_namedfus) %in% vertex,pscooccur1_namedfus)
vsize <- as.data.frame(taxa_sums())
vsize <- count(1:nrow(vertex) %in% )))


#To make the same graph everytime
set.seed(486)

plot(cooccur_graph.species,
     layout = layout_with_fr,
     #label = TRUE,
     #label.font = 3,
     #label.cex = 0.75,
     #layout = cooccur.layout,
     vertex.color = vcols, 
     vertex.shape = "circle",
     vertex.label.cex = 0.75,
     vertex.size = 7, 
     vertex.label = NA, 
     edge.color = ecol,
     edge.width = 1.5)
#legends for the lines
legend(x = -0.75, y = -1.3, 
       #frame = FALSE,
       legend=c("Positive", "Negative"), 
       col = c("blue", "red"), 
       lty = (1)
       )
