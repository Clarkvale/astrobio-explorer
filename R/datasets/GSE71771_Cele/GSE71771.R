library(GEOquery)
library(Biobase)
library(limma)
library(dplyr)
library(org.Ce.eg.db)
source("microarray_functions.R")

gse <- getGEO("GSE71771",GSEMatrix = T, AnnotGPL = T)

ex <- exprs(gse$GSE71771_series_matrix.txt.gz)
boxplot(ex)

ex[which(ex <= 0)] <- NaN

limma::plotDensities(log2(ex))
limma::plotMA(log2(ex))
n.ex <- limma::normalizeCyclicLoess(log2(ex)) #this makes things a bit more palatable

exprs(gse$GSE71771_series_matrix.txt.gz) <- n.ex



levs <- c(rep("NG", 3), rep("MG", 3))
targets <- data.frame((cbind(GSE = gse$geo_accession, Target = levs)))
f <- factor(targets$Target, levels = unique(levs))
design <- model.matrix(~0+f)
colnames(design) <- unique(levs)
fit <- lmFit(gse$GSE71771_series_matrix.txt.gz, design)

contrasts <- makeContrasts( Dif = MG - NG, levels = design)

fit2 <- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2)
tT <- topTable(fit2, confint = T, number = Inf, adjust.method = "fdr")

ses <- ci2se(tT$CI.R, tT$CI.L)
prob <- sapply(tT$B, plogis)

out_tT <- tT %>% dplyr::select(adj.P.Val, P.Value, t, B, logFC, GENE_SYMBOL, 
                               GENE_NAME, CHROMOSOMAL_LOCATION, GENE) %>% 
  
  dplyr::mutate(Standard.Error = ses, Probablity = prob) %>%
  dplyr::rename(Entrez.ID = GENE, Gene.title =  GENE_NAME , Chromosome.annotation = CHROMOSOMAL_LOCATION, Gene.Symbol = GENE_SYMBOL, Platform.ORF = GENE_SYMBOL) %>%
  dplyr::mutate(GO.Function = rep(NA, length(rownames(tT))), GO.Function.ID = rep(NA, length(rownames(tT))),
                                 GO.Component = rep(NA, length(rownames(tT))), GO.Component.ID = rep(NA, length(rownames(tT))),
                                 GO.Process = rep(NA, length(rownames(tT))), GO.Process.ID = rep(NA, length(rownames(tT))))


db <- org.Ce.eg.db
gos <- select(db, keys = as.character(out_tT$Entrez.ID), keytype= "ENTREZID", columns= c("GO", "ONTOLOGY"))





mapped.gos <- org.Ce.egGO
mapped.gos <- as.list(mapped.gos)

compress.ids <- function(entrezid, go.matrix){
  lines <- go.matrix %>% filter(ENTREZID == entrezid)
  cc <- lines %>% filter(ONTOLOGY == "CC")
  bp <- lines %>% filter(ONTOLOGY == "BP")
  mf <- lines %>% filter(ONTOLOGY == "MF")
  
  return(list(cc.go = paste(unique(cc$GO), collapse = "///"), 
              cc.term = paste(unique(cc$TERM), collapse = "///"),
              bp.go = paste(unique(bp$GO), collapse = "///"),
              bp.term = paste(unique(bp$TERM), collapse = "///"),
              mf.go = paste(unique(mf$GO), collapse = "///"),
              mf.term = paste(unique(mf$TERM), collapse = "///")))
}

make.vector <- function(index, list, name){
  return(list[[index]][[name]])
}

get.GOs <- function(TopTable, org.database, Entrez.name){
  gos <- select(db, keys = as.character(TopTable[[Entrez.name]]), keytype= "ENTREZID", columns= c("GO", "ONTOLOGY"))
  u.ids <- as.numeric(unique(gos$ENTREZID))
  TERM <- apply(FUN = Term, as.array(gos$GO), MARGIN = 1)
  gos <- cbind(gos, TERM)
  
  
  listed.gos <- lapply(u.ids, FUN = compress.ids, go.matrix = gos)
  vnames <- c("cc.go", "cc.term", "bp.go", "bp.term", "mf.go", "mf.term")
  
  df_out <- data.frame(u.ids)
  
  for(i in 1:length(vnames)){
    df_out <- cbind(df_out,  sapply(1:length(u.ids), FUN = make.vector, list = listed.gos, name = vnames[i]))
    
  }
  
  colnames(df_out) <- append("Entrez.id", vnames)
  return(df_out)
  
  ##adding null columns
  # TopTable <- TopTable %>% dplyr::mutate(GO.Function = rep(NA, length(rownames(TopTable))), GO.Function.ID = rep(NA, length(rownames(TopTable))),
  #                         GO.Component = rep(NA, length(rownames(TopTable))), GO.Component.ID = rep(NA, length(rownames(TopTable))),
  #                         GO.Process = rep(NA, length(rownames(TopTable))), GO.Process.ID = rep(NA, length(rownames(TopTable))))
  # 
  # #throwing stuff into lists
  # mapped.gos <- as.list(org.egGO)
  # go.db <- as.list(GO.db::GOTERM)
  # db.names <- as.list(org.egSYMBOL)
  # 
  # go.dataframes <- dataframe.go.genome(mapped.gos)
  # valid_ind <- match(names(go.dataframes), out_tT[Entrez.name])
  # 
  # for(i in 1:length(names(go.dataframes))){
  #   if(!is.na(valid_ind[i])){
  #     
  #     if((!is.null(go.dataframes[[i]]$ontology))){  
  #       go_entry <- go.dataframes[[i]]
  #       
  #       TopTable <- add.go.tT(TopTable, go_entry, i)
  #     }
  #     else{
  #       next
  #     }
  #   }
  #   else{
  #     if((!is.null(go.dataframes[[i]]$ontology))){
  #       go_entry <- go.dataframes[[i]]
  #       ided.rows <- grep(TopTable[Entrez.name], pattern = names(go.dataframes)[i])
  #       out_tT <- add.go.tT(TopTable, goTable = go_entry, pos = ided.rows)
  #     }
  #     
  #   }
  #   
  # }
  # return(TopTable)
  # 
  
  
  
}




dataframe.go.genome <- function(go.ids){
  go.db <- as.list(GO.db::GOTERM)
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


#go.list <- as.list(gos$ENTREZID)
#names(go.list) <- gos$GO

go.genome <- anno.go.genome(mapped.gos)
valid_ind <- match(names(go.genome), out_tT$Entrez.ID)

ce.names <- as.list(org.Ce.egSYMBOL)


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
      ided.rows <- grep(out_tT$Gene.symbol, pattern = ce.names[names(go.genome[i])])
      out_tT <- add.go.tT(out_tT, goTable = go_entry, pos = ided.rows)
    }
    
  }
  
}


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

copy.names <- out_tT$Platform.ORF
out_tT <- out_tT %>% mutate(Gene.Symbol = copy.names)
csv_file = "datasets/GSE71771_Cele/GSE71771.csv"
write.csv(out_tT, file = csv_file, sep = ",")


metaname <- "datasets/GSE71771_Cele/GSE71771_meta"
strain <- gse$GSE71771_series_matrix.txt.gz$`strain:ch1`[1]

#fixing the descriptions
desc <- gse$GSE71771_series_matrix.txt.gz$description
desc1 <- desc[1:3]
desc2 <- desc[4:6]

gse$GSE71771_series_matrix.txt.gz$description <- c(desc2,desc1)
extractMetaData(gse$GSE71771_series_matrix.txt.gz, design, list(contrasts), 
                microgravity_type = M.TYPE$SPACEFLOWN, filename = metaname, metaLabels = c(""), strain = strain, cellType = gse$GSE71771_series_matrix.txt.gz$`Stage:ch1`[1])

