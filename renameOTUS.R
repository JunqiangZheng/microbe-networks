###Script written by Lucas to rename the species
<<<<<<< HEAD
rename_otus <- function(ps, max.sp=2, otu_map = FALSE, fasta_file = FALSE){
=======
<<<<<<< HEAD
rename_otus <- function(ps, max.sp=2, otu_map = FALSE, fasta_file = FALSE){
  taxtab <- tax_table(ps)
  taxa_list <- vector()
  seq_list <- vector()
  for(row in 1:nrow(taxtab)){
    taxon <-taxtab[row,]
    seq_list <- c(seq_list,row.names(taxon))
=======
rename_otus <- function(ps, max.sp=2){
>>>>>>> ac4cf4b3e97636c3bc4586826f6da23bd5d61ce8
  taxtab <- tax_table(ps)
  taxa_list <- vector()
  for(row in 1:nrow(taxtab)){
    taxon <-taxtab[row,]
<<<<<<< HEAD
    seq_list <- c(seq_list,row.names(taxon))
=======
>>>>>>> 4577483ec6c5e316b5bc82f681ee1bfd3cc06480
>>>>>>> ac4cf4b3e97636c3bc4586826f6da23bd5d61ce8
    NAs <- sum(is.na(taxon))
    # Genus_species nomenclature for OTUs IDed to species
    if(NAs == 0){
      name <- paste(taxon[,6],taxon[,7], sep="_")
      # Restrict ambiguous species calls
      if(length(unlist(strsplit(name,"/"))) > max.sp){
        NAs <- 1
      }
    }
<<<<<<< HEAD
    # If otu is named to Genus, use "[Genus]_sp" convention
    if(NAs == 1){
      name <- paste(taxon[,6],"sp", sep="_")
=======
<<<<<<< HEAD
    # If otu is named to Genus, use "[Genus]_sp" convention
    if(NAs == 1){
      name <- paste(taxon[,6],"sp", sep="_")
=======
    # If OTU is named to Genus, use "[Genus]_sp" convention
    if(NAs == 1){
      name <- paste(taxon[,length(taxon)-NAs],"sp", sep="_")
>>>>>>> 4577483ec6c5e316b5bc82f681ee1bfd3cc06480
>>>>>>> ac4cf4b3e97636c3bc4586826f6da23bd5d61ce8
    }
    # If OTU identified to higher taxonomic level, just use name of that level
    if(NAs > 1){
      name <- paste(taxon[,length(taxon)-NAs])
    }
    taxa_list <- c(taxa_list,name)
  }
  # Rename duplicate taxa in taxa_list, by adding numbers to the end
  for(taxon in 1:length(taxa_list)){
    matches <- taxa_list[taxon] == taxa_list
    nmatches <- sum(matches)
    if(nmatches > 1){
      taxa_list[matches] <- paste(taxa_list[matches],1:(nmatches),sep = "_") 
    }
  }
<<<<<<< HEAD
=======
<<<<<<< HEAD
>>>>>>> ac4cf4b3e97636c3bc4586826f6da23bd5d61ce8
  if(otu_map == TRUE){
    for(i in 1:length(taxa_list)){
      line <- paste0(taxa_list[i], "\t",seq_list[i])
      write(line, file="otu_name_mapfile.txt",append = TRUE)
    }
  }
  if(fasta_file == TRUE){
    for(i in 1:length(taxa_list)){
      line <- paste0(">",taxa_list[i],"\n",seq_list[i])
      write(line, file="16S_rep_set_renamed.fasta", append=TRUE)
    }
  }
  taxa_names(ps) <- taxa_list
  return(ps)
  closeAllConnections()
<<<<<<< HEAD
=======
=======
  taxa_names(ps) <- taxa_list
  return(ps)
>>>>>>> 4577483ec6c5e316b5bc82f681ee1bfd3cc06480
>>>>>>> ac4cf4b3e97636c3bc4586826f6da23bd5d61ce8
}
###End of scipt written by Lucas