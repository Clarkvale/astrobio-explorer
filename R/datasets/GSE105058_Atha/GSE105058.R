library(GEOquery)
library(limma)
source("microarray_functions.R")
library(readr)
library(dplyr)
gse <- getGEO("GSE105058", AnnotGPL = T, GSEMatrix = T)[[1]]

#limma::plotDensities(gse)
#limma::plotMA(gse)

s <- c(1,2,3)
g <- c(4,5,6)


levs <-  c(rep("MG",3), rep("NG", 3))



targets <- data.frame(cbind(GSE = gse$geo_accession, Target = levs))
f <- factor(targets$Target, levels = unique(levs))
design <- model.matrix(~0+f)
colnames(design) <- unique(levs)
fit <- lmFit(gse, design)


cont.dif <- makeContrasts(
  Dif = MG - NG,
  levels = design
)
contrasts.fit <- contrasts.fit(fit, contrasts = cont.dif)

fit2 <- eBayes(contrasts.fit, 0.01)
tT <- topTable(fit2, adjust.method = "fdr", confint = TRUE, number = Inf)
tT <- remove.controls(tT)$TopTable

ses <- topTable.SE(tT)
prob <- sapply(tT$B, plogis)

#curating gene symbols
aliases <- read_tsv("datasets/GSE105058_Atha/gene_aliases_20191231.txt")
matched <-  match(toupper(tT$Platform_ORF), aliases$name )

matched_symbols <- aliases$symbol[matched]
tT$Gene.symbol <- matched_symbols

nas <- which(is.na(tT$Gene.symbol))
tT$Gene.symbol[nas] <- tT$Platform_ORF[nas]


#outputting final table 
out_tT <- tT %>% dplyr::select(adj.P.Val, P.Value, t, B, logFC, Gene.symbol, 
                               Gene.title, Platform_ORF, GO.Function, GO.Function.ID,
                               GO.Process, GO.Process.ID, GO.Component, GO.Component.ID,
                               Chromosome.annotation, ID, Gene.ID) %>% 
                               
                 dplyr::mutate(Standard.Error = ses, Probablity = prob) %>%
                 dplyr::rename(Entrez.ID = Gene.ID)

output.tT <- pull.output.tT(tT)


write.csv(out_tT, file = "datasets/GSE105058_Atha/GSE105058.csv")

extractMetaData(gse = gse, metaLabels = c("Atha"), 
                microgravity_type = M.TYPE$SPACEFLOWN, design = design, 
                contrasts = list(cont.dif), cellType = gse$`tissue:ch1`[1],
                filename = "datasets/GSE105058_Atha/GSE105058_meta", 
                description_label = gse$`factor value [microgravity]:ch1`)
