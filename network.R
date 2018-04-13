#16S##############################################################
cooccur.species <- readRDS(file = "rds/cooccur.species.10.rds")

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
dim(strong_cooccur.species)
#strong_cooccur.species <- strong_cooccur.species[,-c(12,13)]
#colnames(strong_cooccur.species)[10:11] <- c('OTU_A','OTU_B')

head(strong_cooccur.species, n=2)

options(repr.plot.width = 10,repr.plot.height = 6)
plot(sort(strong_cooccur.species$effect), 
     xlab='Species pairs',
     ylab = 'Cooccurrence effect size',
     main='Cooccurrence effect sizes')

dim(strong_cooccur.species)

pos.effect <- strong_cooccur.species[,c(1, 2, 3)]
neg.effect <- negative_cooccur.species[,c(1, 2, 3)]
neg.effect$effects <- abs(neg.effect$effects)
posneg.effect <- rbind(pos.effect, neg.effect)
posneg.effect$effects <- posneg.effect$effects * 100

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
ps_vertices <- prune_taxa(taxa_names(ps.16S) %in% vertex,ps.16S)
nodes16S <- cbind(1:24,as.data.frame(tax_table(ps_vertices)))
#write.csv(nodes16S, "images/16Snodes.csv")


#To make the same graph everytime
set.seed(486)
network.16S <- plot(cooccur_graph.species,
     #layout = layout_with_fr,
     label = TRUE,
     #label.font = 3,
     label.cex = 0.75,
     #layout = cooccur.layout,
     vertex.color = vcols, 
     vertex.shape = "circle",
     vertex.label.cex = 0.75,
     vertex.label = vertex, 
     vertex.label.color = "blue",
     edge.color = ecol,
     edge.width = 1.5)
#legends for the lines
#legend(x = -0.75, y = -1.3, 
#       #frame = FALSE,
 #      legend=c("Positive", "Negative"), 
  #     col = c("blue", "red"), 
   #    lty = (1)
    #   )



#ITS###################################################################
cooccur.ITS <- readRDS(file = "rds/cooccur.ITS.10.rds")

cooccur.ITSP <- prob.table(cooccur.ITS)
head(cooccur.ITSP, n=2)

cooccur.ITSC <- effect.sizes(cooccur.ITS)

cooccur.ITSPC <- cbind(cooccur.ITSC, cooccur.ITSP)


cooccur.ITSPC$p_gt_adj <- p.adjust(cooccur.ITSPC$p_gt, method = "BH")
cooccur.ITSPC$p_lt_adj <- p.adjust(cooccur.ITSPC$p_lt, method = "BH")
strong_cooccur.ITS <- cooccur.ITSPC[cooccur.ITSPC$p_gt_adj <= 0.05,]
negative_cooccur.ITS <- cooccur.ITSPC[cooccur.ITSPC$p_lt <= 0.05,] 
dim(negative_cooccur.ITS)


pos.ITS.effect <- strong_cooccur.ITS[,c(1, 2, 3)]
neg.ITS.effect <- negative_cooccur.ITS[,c(1, 2, 3)]
neg.ITS.effect$effects <- abs(neg.ITS.effect$effects)
posneg.ITS.effect <- rbind(pos.ITS.effect, neg.ITS.effect)
posneg.ITS.effect$effects <- posneg.ITS.effect$effects * 100

pos.ITS <- strong_cooccur.ITS[,c(1, 2)]
neg.ITS <- negative_cooccur.ITS[,c(1,2)]
posneg.ITS <- rbind(pos.ITS, neg.ITS)
cooccur_graph.ITS <- graph_from_data_frame(posneg.ITS, directed=FALSE)
pos.graph.ITS <- graph_from_data_frame(pos.ITS, directed = FALSE)
neg.graph.ITS <- graph_from_data_frame(neg.ITS, directed = FALSE)


vcols.ITS <- vector(length = length(V(cooccur_graph.ITS)))
vcols.ITS[] <- 'black'
vcols.ITS[grep("^Fusarium", names(V(cooccur_graph.ITS)))] <- "pink"
ecol.ITS <- rep("red", ecount(cooccur_graph.ITS))
ecol.ITS[E(cooccur_graph.ITS) %in% E(pos.graph.ITS)] <- "blue"

vertex.ITS <- cbind(names(V(cooccur_graph.ITS)), 1:29)
ps_vertices.ITS <- prune_taxa(taxa_names(ps.ITS) %in% vertex.ITS,ps.ITS)
write.csv(as.data.frame(tax_table(ps_vertices.ITS)), "images/ITSnodes.csv")

vertex.attributes((cooccur_graph.ITS[-grep("^Fusarium", names(V(cooccur_graph.ITS)))])

vname = list()
vnames <- for(i in 1:29){
  if(vertex.ITS[i,1] == grep("^Fusarium", names(V(cooccur_graph.ITS)))){
    vnames[i] = 'F'
  }else{
    vnames[i] = vertex.ITS[i,2]
    
  }
}


set.seed(55)
plot(cooccur_graph.ITS,
     layout = layout_with_fr,
     label = TRUE,
     vertex.color = vcols.ITS, 
     vertex.shape = "circle",
     vertex.label = 0.75,
     vertex.label.cex = vertex.ITS[,2],
     edge.color = ecol.ITS,
     edge.width = 1.5)


#both##################################################################
cooccur.both <- readRDS(file = "rds/cooccur.both.rds")

cooccur.bothP <- prob.table(cooccur.both)
head(cooccur.bothP, n=2)

cooccur.bothC <- effect.sizes(cooccur.both)

cooccur.bothPC <- cbind(cooccur.bothC, cooccur.bothP)

cooccur.bothPC$p_gt_adj <- p.adjust(cooccur.bothPC$p_gt, method = "BH")
cooccur.bothPC$p_lt_adj <- p.adjust(cooccur.bothPC$p_lt, method = "BH")
strong_cooccur.both <- cooccur.bothPC[cooccur.bothPC$p_gt_adj <= 0.05,]
negative_cooccur.both <- cooccur.bothPC[cooccur.bothPC$p_lt_adj <= 0.05,] 
dim(negative_cooccur.both)

pos.both.effect <- strong_cooccur.both[,c(1, 2, 3)]
neg.both.effect <- negative_cooccur.both[,c(1, 2, 3)]
neg.both.effect$effects <- abs(neg.both.effect$effects)
posneg.both.effect <- rbind(pos.both.effect, neg.both.effect)
posneg.both.effect$effects <- posneg.both.effect$effects * 100

pos.both <- strong_cooccur.both[,c(1, 2)]
neg.both <- negative_cooccur.both[,c(1,2)]
posneg.both <- rbind(pos.both, neg.both)
cooccur_graph.both <- graph_from_data_frame(posneg.both, directed=FALSE)
pos.graph.both <- graph_from_data_frame(pos.both, directed = FALSE)
neg.graph.both <- graph_from_data_frame(neg.both, directed = FALSE)


vcols.both <- vector(length = length(V(cooccur_graph.both)))
vcols.both[] <- 'black'
vcols.both[which(names(V(cooccur_graph.both)) == "Fusarium_vertillioides")] <- "pink"
ecol.both <- rep("red", ecount(cooccur_graph.both))
ecol.both[E(cooccur_graph.both) %in% E(pos.graph.both)] <- "blue"
#ecol[E(cooccur_graph.species) %in% E(neg.graph)] <- "red"
#ecol <- rep("red", ecount(cooccur_graph.species))
#ecol[(cooccur_graph.species) %in% pos.graph] <- "blue"
#ecol[(cooccur_graph.species) %in% neg.graph] <- "red"

vertex.both <- names(V(cooccur_graph.both))
ps_vertices.both <- prune_taxa(taxa_names(ps.both.new10) %in% vertex,ps.both.new10)
write.csv(as.data.frame(tax_table(ps_vertices.both)), "images/bothnodes.csv")

set.seed(486)
plot(cooccur_graph.both,
     layout = layout_with_fr,
     #label = TRUE,
     #label.font = 3,
     #label.cex = 0.75,
     #layout = cooccur.layout,
     vertex.color = vcols.both, 
     vertex.shape = "circle",
     vertex.label.cex = 0.75,
     #vertex.size = vsize.taxasums, 
     vertex.label = NA, 
     edge.color = ecol.both,
     edge.width = 1.5)


