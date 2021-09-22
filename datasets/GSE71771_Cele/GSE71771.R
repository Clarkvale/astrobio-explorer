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
  dplyr::rename(Entrez.ID = GENE, Gene.title =  GENE_NAME , Chromosome.annotation = CHROMOSOMAL_LOCATION, Gene.symbol = GENE_SYMBOL, Platform.ORF = GENE_SYMBOL) %>%
  dplyr::mutate(GO.Function = rep(NA, length(rownames(tT))), GO.Function.ID = rep(NA, length(rownames(tT))),
                                 GO.Component = rep(NA, length(rownames(tT))), GO.Component.ID = rep(NA, length(rownames(tT))),
                                 GO.Process = rep(NA, length(rownames(tT))), GO.Process.ID = rep(NA, length(rownames(tT))))


db <- org.Ce.eg.db

gos <- get.GOs(out_tT, db, "Entrez.ID")


tT <-  out_tT[-which(duplicated(out_tT$Entrez.ID)),]
tT <- tT[-which(is.na(tT$Entrez.ID)),]
tT <- cbind(tT, gos)


tT <- tT[-which(is.na(tT$P.Value)),]



 Gene.symbol <- tT$Platform.ORF
 tT <- cbind(tT, Gene.symbol) 

 
csv_file = "datasets/GSE71771_Cele/GSE71771.csv"
write.csv(tT, file = csv_file, sep = ",")


metaname <- "datasets/GSE71771_Cele/GSE71771_meta"
strain <- gse$GSE71771_series_matrix.txt.gz$`strain:ch1`[1]

#fixing the descriptions
desc <- gse$GSE71771_series_matrix.txt.gz$description
desc1 <- desc[1:3]
desc2 <- desc[4:6]

gse$GSE71771_series_matrix.txt.gz$description <- c(desc2,desc1)
extractMetaData(gse$GSE71771_series_matrix.txt.gz, design, list(contrasts), 
                microgravity_type = M.TYPE$SPACEFLOWN, filename = metaname, metaLabels = c(""), strain = strain, cellType = gse$GSE71771_series_matrix.txt.gz$`Stage:ch1`[1])

