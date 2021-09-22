##Benjamin Clark and Marie-May B
library(limma)
library(edgeR)
source("microarray_functions.R")
library(org.Hs.eg.db)
library(GEOquery)
library(dplyr)

counts <- read.delim("datasets/GSE151411_Hsap/raw_counts.tabular", row.names = 1)
genes <- rownames(counts)
counts <- apply(counts, as.numeric, MARGIN = 2)

rownames(counts) <- genes

d0 <- DGEList(na.omit(counts))

d0 <- calcNormFactors(d0)



d <- filterByCPM(d0)

gse <- getGEO("GSE151411")[[1]]
levs <- make.names(sapply(gse$title, function(name){return(substr(name, 1, nchar(name)-2))}, USE.NAMES = F))
targets <- data.frame(cbind(GSE = gse$geo_accession, Target = levs))
f <- factor(targets$Target, levels = unique(levs))
design <- model.matrix(~0+f)
colnames(design) <- unique(levs)

y <- voom(d, design, plot = TRUE) #looks okay


fit <- lmFit(y, design)
Contrast <- makeContrasts( SMG.Infected - NG.Infected, levels = colnames(coef(fit)))
fit2 <- contrasts.fit(fit, Contrast)
fit2 <- eBayes(fit2)
tT <- topTable(fit2, n = Inf, adjust.method = "fdr", confint = T)


#pulling metadata
names <- select(org.Hs.eg.db, keys=rownames(tT), keytype="ENTREZID", columns="GENENAME")
symbols <- select(org.Hs.eg.db, keys=rownames(tT), keytype="ENTREZID", columns="SYMBOL")


tT <- tT %>% mutate(Gene.symbol = symbols$SYMBOL, Gene.title = names$GENENAME, 
                    Entrez.ID = rownames(tT),Platform.ORF = NA, 
                    Probability = sapply(tT$B, plogis),
                    Standard.Error = ci2se(tT$CI.R, tT$CI.L))



tT.out <- get.GOs(tT, org.Hs.eg.db, Entrez.name = "Entrez.ID") 

tT.out <- cbind(tT,tT.out)
write.csv(tT.out, file = "datasets/GSE151411_Hsap/GSE151411.csv", row.names = F)



metaName <- "datasets/GSE151411_Hsap/GSE151411_meta"
strain <- ""

extractMetaData(gse , filename = metaName, microgravity_type = M.TYPE$HARV, 
                strain = strain, contrasts = Contrast, design = design, 
                cellType = "EPEC infected U937 (ATCC CRL-1593.2)")

