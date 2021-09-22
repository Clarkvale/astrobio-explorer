
#GSE40648 Ecoli MICROARRAY
#R script by BENJAMIN CLARK
#Analysis by Hajar El Mouddene

library(GEOquery)
library(Biobase)
library(limma)
library(org.EcK12.eg.db)
source("microarray_functions.R")
library(GO.db)
library(dplyr)
#This pulls the all the samples from the microarry dataset. This returns a list containing a single expressionSet object. 
list.gse <- getGEO("GSE40648", GSEMatrix =TRUE, AnnotGPL=TRUE)

gse <- list.gse[[1]]

#take a look into the dataset
#head(gse)
#View(gse$title)
#dim(gse)

#here we format the feature names. fvarLabels belongs the biobase package and is used to extract features from ExpressionSet Objects
fvarLabels(gse) <- make.names(fvarLabels(gse))


#DE analysis


gse.rm <- gse[,-c(4,8)]
control <- c(1,2,3)
treatment <- c(4,5,6)

ecoli <- de.analysis(gse = gse.rm, microgravity_group = treatment, ground_group = control)

#print out the toptable 
ecoliname <- "datasets/GSE40648_Ecoli/GSE40648_Ecoli.csv"
out_tT <- pull.output.tT(ecoli$TopTable)






#adding new GO ids
symbols <- out_tT$Gene.symbol
symbols.1 <- sapply(symbols, FUN = function(x){return(strsplit(x, "///")[[1]][1])})
names(symbols.1) <- NULL
symbols.1 <- as.vector(symbols.1)

esym <- org.EcK12.egSYMBOL2EG
esym <- as.list(esym)
entrez_ids <- unlist(esym[symbols.1])

#for whatever reason the entrez ids supplied by GEO are mapping poorly to the 
#annotationDB so I need to use the above code for some extra mapping from gene symbols

#lets use all of them for good measure, i need to make sure i can still use gene symbols later for mapping 
geo_entrez_ids <-  out_tT$Entrez.ID 
geo_entrez_ids <- sapply(geo_entrez_ids, FUN = function(x){return(strsplit(x, "///")[[1]][1])})
names(geo_entrez_ids) <- out_tT$Gene.symbol



#combining the two
entrez_ids <- union(entrez_ids, geo_entrez_ids)

#mapping gos to our ids
gos <- org.EcK12.egGO
gos <- as.list(gos)

mapped_gos<- gos[entrez_ids]
mapped_gos <- mapped_gos[which(!is.na(names(mapped_gos)))]


add.go.tT <- function(tT, goTable, pos){
  
  if(any(grepl(dplyr::select(goTable, ontology)[,1], pattern = "MF"))){
    tT$GO.Function[pos] <- paste(dplyr::filter(goTable, ontology == "MF")$term , collapse = "///")
    tT$GO.Function.ID[pos] <- paste(dplyr::filter(goTable, ontology == "MF")$id, collapse = "///")
  }
  #GO ancestor ids are apparently not packaged with annoDB, these ones are found here:https://www.ebi.ac.uk/QuickGO/
  else{
    tT$GO.Function[pos]<- "molecular_function"
    tT$GO.Function.ID[pos] <- "GO:0003674"
  }
  
  if(any(grepl(dplyr::select(goTable, ontology)[,1], pattern = "BP"))){
    tT$GO.Process[pos] <- paste(dplyr::filter(goTable, ontology == "BP")$term , collapse = "///")
    tT$GO.Process.ID[pos] <- paste(dplyr::filter(goTable, ontology == "BP")$id, collapse = "///")
  }
  else{
    tT$GO.Process[pos] <- "regulation of biological process"
    tT$GO.Process.ID[pos] <- "GO:0050789"
  }
  
  if(any(grepl(dplyr::select(goTable, ontology)[,1], pattern = "CC"))){
    tT$GO.Component[pos] <- paste(dplyr::filter(goTable, ontology == "CC")$term , collapse = "///")
    tT$GO.Component.ID[pos] <- paste(dplyr::filter(goTable, ontology == "CC")$id, collapse = "///")
  }
  else{
    tT$GO.Component[pos] <- "cellular_component"
    tT$GO.Component.ID[pos]<- 'GO:0005575'
  }
  return(tT)
  
}
 
#getting go terms for my tables
go.db <- as.list(GO.db::GOTERM)
anno.go.genome <- function(go.ids){
  out <- list()
  for(i in 1:length(go.ids)){
    line_ids <- names(go.ids[[i]])
    line.go.items <- list(go.db[line_ids])
    terms <- sapply(names(line.go.items[[1]]),FUN = Term)
    names(terms) <- NULL
    terms <- unlist(terms)
    
    ontologies <- sapply(names(line.go.items[[1]]), FUN = Ontology)
    names(ontologies) <- NULL
    ontologies <- unlist(ontologies)
    
    inner_out <- data.frame(id = names(line.go.items[[1]]), term = terms, ontology = ontologies)
    out[names(go.ids)[i]]<- list(inner_out)
    
    
    
  }
  return(out)
}

eck12.names <- as.list(org.EcK12.egSYMBOL)
#mapping mined GOs to our toptable, this else if chain is kinda ugly but it do
#this is me just refusing to vectorize this for some reason, it very slow
go.genome <- anno.go.genome(mapped_gos)
valid_ind <- match(names(go.genome), out_tT$Entrez.ID)
for(i in 1:length(names(go.genome))){
  if(!is.na(valid_ind[i])){
    
    if((!is.null(go.genome[[i]]$ontology))){  
      go_entry <- go.genome[[i]]
      
      out_tT <- add.go.tT(out_tT, go_entry, i)
    }
    else{
      next
    }
  }
  else{
    if((!is.null(go.genome[[i]]$ontology))){
      go_entry <- go.genome[[i]]
      ided.rows <- grep(out_tT$Gene.symbol, pattern = eck12.names[names(go.genome[i])])
      out_tT <- add.go.tT(out_tT, goTable = go_entry, pos = ided.rows)
    }
    
  }
 
}




write.table(remove.controls(out_tT)$TopTable, ecoliname, row.names = FALSE, sep = ",")


#process and extract metadata for all datasets comparisons

metaName <- "datasets/GSE40648_Ecoli/GSE40648_meta"
strain <- "K12 MG1655"
gse_list <- list(ecoli)
labels <- c("ecoli")
#extractMetaData(filename = metaName, gse_groups = gse_list, microgravity_type = M.TYPE$RPM, metaLabels = labels, strain = strain)

