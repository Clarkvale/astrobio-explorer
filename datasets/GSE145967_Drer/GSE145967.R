##Benjamin Clark
library(limma)
library(edgeR)
library(org.Dr.eg.db)
source("microarray_functions.R")
library(GEOquery)

counts <- read.delim("datasets/GSE145967_Drer/raw_counts.tabular", row.names = 1)
genes <- rownames(counts)
counts <- apply(counts, as.numeric, MARGIN = 2)

rownames(counts) <- genes

d0 <- DGEList(na.omit(counts))


d0 <- calcNormFactors(d0)

#filtering low expressed genes
cutoff <- 0.6
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,] 
dim(d) 

#factor names
gravity <- append(rep("NG", 3), rep("MG",3))
plotMDS(d, col = as.numeric(factor(gravity)))

mm <- model.matrix(~0 + gravity)
y <- voom(d, mm, plot = TRUE) #looks okay

colnames(mm) <- levels(factor(gravity))

fit <- lmFit(y, mm)
Contrast <- makeContrasts(MG - NG , levels = colnames(coef(fit)))
fit2 <- contrasts.fit(fit, Contrast)
fit2 <- eBayes(fit2)
tT <- topTable(fit2, n = Inf, adjust.method = "fdr", confint = T)


##adding annotations
annos <- read.columns("datasets/GSE145967_Drer/annos.tabular")
vmatch <- match(rownames(tT),annos$ENSEMBL)

vmatch_id <- annos[vmatch,]
tT <- cbind(tT, vmatch_id)
tT <- na.omit(tT)
tT <- tT[-which(duplicated(tT$ENTREZID)),]

#getting GOs 
gos <- get.GOs(tT, org.database = org.Dr.eg.db, Entrez.name = "ENTREZID")
tT <- cbind(tT, gos)

#formatting final table
tT_out <- tT %>% mutate(Platform.ORF = NA, 
                           Probability = sapply(tT$B, plogis), 
                           Standard.Error = ci2se(tT$CI.R, tT$CI.L)) %>%
                           rename(Gene.symbol = SYMBOL, Gene.title = GENENAME) %>%
                           select( -ENSEMBL, -ENTREZID)


##metaData
gse <- getGEO("GSE145967", GSEMatrix = TRUE, AnnotGPL=FALSE)[[1]]

gse <- gse[,c(1,2,3,7,8,9)]

extractMetaData(gse, design = mm, contrasts = list(Contrast), filename = "datasets/GSE145967_Drer/GSE145967_meta", 
                metaLabels = c(""), strain = gse$`strain:ch1`[1],cellType = "embryo", microgravity_type = M.TYPE$RCCS)

##printing out
write.table(tT_out, "datasets/GSE145967_Drer/GSE145967.csv", row.names = FALSE, sep = ",")



