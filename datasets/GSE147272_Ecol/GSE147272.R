library(edgeR)
library(dplyr)

#Reformatting output from featureCount

counts <- read.delim("raw_counts.tabular", row.names = 1)
genes <- rownames(counts)
counts <- na.omit(sapply(counts, as.numeric))
rownames(counts) <- genes[-1]
d0 <- DGEList(counts)


#Adding Normalizing factors
d0 <- calcNormFactors(d0)

#filtering low expressed genes
cutoff <- 1
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,] 
dim(d) 

#Setting up factor names
gravity <- append(rep("MG", 4), rep("NG",4))
growth.phase <- c("exp", "exp", "sta", "sta", "exp", "exp","sta", "sta")

group <- interaction(gravity, growth.phase)

#Checking clustering
plotMDS(d, col = as.numeric(group))


#Voom and modeling
mm <- model.matrix(~0 + group)
y <- voom(d, mm, plot = TRUE) #looks okay

fit <- lmFit(y, mm)
exp.Contrast <- makeContrasts(groupMG.exp - groupNG.exp , levels = colnames(coef(fit)))
fit2 <- contrasts.fit(fit, exp.Contrast)
fit2 <- eBayes(fit2)
tT.exp <- topTable(fit2, n = Inf, adjust.method = "fdr", confint = T)


#getting annotations 
anno <- rtracklayer::import("ecoli_nissle_1917.gtf")

anno  <-  anno %>% as.data.frame %>% filter(type == "exon") %>% select(-seqnames, -source, -score, -phase)
vmatch <- match(rownames(tT.exp),anno$gene_id)

vmatch_id <- anno$gene_name[vmatch]
tT.exp <- tT.exp %>% mutate(Gene.Symbol = vmatch_id, Platform.ORF = rownames(tT.exp), 
                            Probability = sapply(tT.exp$B, plogis), Standard.Error = ci2se(tT.exp$CI.R, tT.exp$CI.L)) %>% na.omit() 
write.csv(tT.exp, file = "datasets/GSE147272_Ecol/GSE147272_exponential.growth.csv")

#now for stationary phase
sta.Contrast <- makeContrasts(groupMG.sta - groupNG.sta , levels = colnames(coef(fit)))
fit3 <- contrasts.fit(fit, sta.Contrast)
fit3 <- eBayes(fit3)
tT.sta <- topTable(fit3, n = Inf, adjust.method = "fdr", confint = T)

vmatch <- match(rownames(tT.sta),anno$gene_id)


vmatch_id <- anno$gene_name[vmatch]
tT.sta <- tT.sta %>% mutate(Gene.Symbol = vmatch_id, Platform.ORF = rownames(tT.sta),
                            Probability = sapply(tT.sta$B, plogis), Standard.Error = ci2se(tT.sta$CI.R, tT.sta$CI.L)) %>% na.omit() 
write.csv(tT.sta, file = "datasets/GSE147272_Ecol/GSE147272_stationary.growth.csv")

#getting metadata
mdata <- GEOquery::getGEO("GSE147272")
mdata <- mdata$GSE147272_series_matrix.txt.gz

source("microarray_functions.R")
metaName <- "datasets/GSE147272_Ecol/GSE147272_meta"
strain <- "Nissle_1917"

gselist <- list(list(GSE = mdata[,c(1,2,5,6)]), list(GSE = mdata[,c(3,4,7,8)]))
metalabels <- c("exponential.growth", "stationary.growth")
extractMetaData(gse_groups = gselist, filename = metaName, microgravity_type = M.TYPE$RCCS, metaLabels = metalabels, strain = strain)

