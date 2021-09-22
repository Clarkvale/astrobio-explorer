###GSE50881 YEAST MICROARRAY TWO CHANNEL
#BENJAMIN CLARK 

library(GEOquery)
library(Biobase)
library(limma)
#library(arrayQualityMetrics)
#library(FoldGO)
library(GO.db)
library(dplyr)
source("microarray_functions.R")


list.gse <- getGEO("GSE50881", GSEMatrix =TRUE, AnnotGPL=TRUE)

gse <- list.gse[[1]]

#take a look into the dataset
#head(gse)
#View(gse$title)
#dim(gse)


fvarLabels(gse) <- make.names(fvarLabels(gse))




design <- cbind(DyeEffect = 1,SpaceVsGround = c(1,-1,1,-1,1,-1,1,-1))
#block <- fData(gse[,1])$Block


fit <- lmFit(gse,design)
fit <- eBayes(fit, 0.01)

topTable.dyes <- topTable(fit, coef = "DyeEffect")

anno.data <- fData(gse)


topTable.SvsG <- topTable(fit, adjust.method = "fdr",  coef = "SpaceVsGround", sort.by = "B", number = length(fit[[1]]), confint = TRUE)
f.tT <- remove.candida.controls(topTable.SvsG)

#f.tT <- subset(f.tT, select = c("CGD_Systematic_Name", "Description", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "SPOT_ID", "ID", "B", "CI.L", ))



rowid <- c() 
for(i in 1:length(f.tT$ID)){
  rowid <- append(rowid,((f.tT %>% filter(SPOT_ID == f.tT$SPOT_ID[i]) %>% arrange(P.Value))[1,])$ID)
  
}

best.rows.ids <- rowid[-which(duplicated(rowid))]
rm.ftT <- na.omit(f.tT[best.rows.ids,]) %>% arrange(adj.P.Val)

match.id.lowest.pval <- function(id){
  best.p <- as.data.frame(f.tT) %>% dplyr::filter(SPOT_ID == id) %>% dplyr::arrange(P.Value)[1,]
  return(best.p)
}



remove.candida.controls <- function(toptable){
  controls <- which(toptable$CGD_Systematic_Name == "")
  passed.probes <- which(toptable$CGD_Systematic_Name != "")
  return(topTable.SvsG[passed.probes,])
}






##getting GO annotations locally
setwd("C:/Users/Benja/repos/de-expression-git/de-expression-igem/datasets/GSE50881_Calb")
go_ids <- read.csv("Candida_albicans_GO.csv")
ordered.gos <- list()
for(annos in 1:length(rm.ftT$CGD_Systematic_Name)){
  indexes <- which(rm.ftT$CGD_Systematic_Name[annos] == go_ids$systemic_id)
  ids <- go_ids$go_id[indexes]  
  ordered.gos[rm.ftT$CGD_Systematic_Name[annos]] <- list(ids)
  
}
#expects go ids in a vector contained in a list 
make.GO.table <- function(GOIDs){
  GOIDs <- unname(unlist(GOIDs))
  #GOIDs <- as.vector(unname(sapply(full.BP.ID[1], FUN = function(x){return(strsplit(x, "///")[[1]])})))
  return(data.frame(ID = GOIDs, ontology = Ontology(GOIDs), term = Term(GOIDs)))
}



go_tables <- lapply(ordered.gos, FUN = make.GO.table)

#setting up null columns 
rm.ftT <- rm.ftT %>% mutate(GO.Function = rep(NA, length(rownames(rm.ftT))),
                            GO.Function.ID = rep(NA, length(rownames(rm.ftT))),
                            GO.Component = rep(NA, length(rownames(rm.ftT))),
                            GO.Component.ID = rep(NA, length(rownames(rm.ftT))),
                            GO.Process = rep(NA, length(rownames(rm.ftT))),
                            GO.Process.ID = rep(NA, length(rownames(rm.ftT))))

apply_gos <- function(tT, GO_tables){
  gene_ind <- which(names(GO_tables) == tT$CGD_Systematic_Name)
  
  if(length(GO_tables[1]) != 0){
    #Function
    f.rows <- filter(GO_tables[[1]], ontology == "MF")
    while(!is.null(f.rows)){
      tT$GO.Function[gene_ind] <- paste(f.rows$term, collapse = "///")
      tT$GO.Function.ID[gene_ind] <- paste(f.rows$ID, collapse = "///")
      #Process
      f.rows <- filter(GO_tables[[1]], ontology == "BP")
      tT$GO.Process[gene_ind] <-  paste(f.rows$term, collapse = "///")
      tT$GO.Process.ID[gene_ind] <- paste(f.rows$ID, collapse = "///")
      #Component
      f.rows <- filter(GO_tables[[1]], ontology == "CC")
      tT$GO.Component[gene_ind] <- paste(f.rows$term, collapse = "///")
      tT$GO.Component.ID[gene_ind] <- paste(f.rows$ID, collapse = "///")
      break
    }
  }
  return(tT)
  
}
for(i in 1:length(go_tables)){
  rm.ftT <- apply_gos(tT = rm.ftT, GO_tables = go_tables[i])
  
}

names <- sapply(rm.ftT$Description, USE.NAMES = FALSE, FUN = function(x){
  if(startsWith(x , prefix = "|")){
    return(NA)
  }
  else{
    return(stringr::str_extract(strsplit(x, split = "|", fixed = TRUE)[[1]][1], pattern = "\\w+"))
  }
})
names <- unlist(names)

Chromosome_Location <- stringr::str_extract_all(rm.ftT$Description, pattern = "(Contig\\d+:\\D+\\d+..\\d+\\)|Contig\\d+:\\d+..\\d+)")
Chromosome_Location <- sapply(Chromosome_Location, FUN = paste, collapse = "///", USE.NAMES = FALSE)


final.tT <- rm.ftT %>% mutate(
                               Gene.symbol = names, 
                               Chromosome.Location = Chromosome_Location, 
                               ID = NULL, SPOT_ID = NULL, Description = NULL,
                               SEQUENCE = NULL, Block = NULL, Column = NULL, Row = NULL,
                               Standard.Error = topTable.SE(rm.ftT), Probability = sapply(rm.ftT$B, plogis)) %>% 
          rename( Platform_ORF = CGD_Systematic_Name)

#final.tT <- pull.output.tT(final.tT)

write.csv(final.tT, file = "datasets/GSE50881_Calb/GSE50881.csv")


extractMetaData(gse_groups = list(list(GSE = gse[,c(1,3,5,7)]), list(GSE = gse[,c(2,4,6,8)])), microgravity_type = M.TYPE$SPACEFLOWN, 
                filename = "datasets/GSE50881_Calb/GSE50881_meta", metaLabels = c("dye_swap1", "dye_swap2"), strain = "SC5413")
