#GSE90166
#Benjamiin Clark 

library(limma)
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
cutoff <- 0.008
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,] 
dim(d) 

#Setting up factor names
gravity <- append(rep("NG", 3), rep("MG",3))


#Checking clustering
plotMDS(d, col = as.numeric(factor(gravity)))

#Voom and modeling
mm <- model.matrix(~0 + gravity)
y <- voom(d, mm, plot = TRUE) #this isnt super

fit <- lmFit(y, mm)
exp.Contrast <- makeContrasts(gravityMG - gravityNG , levels = colnames(coef(fit)))
fit2 <- contrasts.fit(fit, exp.Contrast)
fit2 <- eBayes(fit2) 
tT.exp <- topTable(fit2, n = Inf, adjust.method = "fdr", confint = T)


#getting annotations 
anno <- rtracklayer::import("acc.gtf")

anno  <-  anno %>% as.data.frame %>% filter(type == "exon") %>% select(-seqnames, -source, -score, -phase)
vmatch <- match(rownames(tT.exp),anno$gene_id)

vmatch_id <- anno$gene_name[vmatch]
tT.exp <- tT.exp %>% mutate(Gene.Symbol = vmatch_id, Platform.ORF = rownames(tT.exp),
                            Probability = sapply(tT.exp$B, plogis), Standard.Error = ci2se(tT.exp$CI.R, tT.exp$CI.L) ) %>% na.omit() 

write.csv(tT.exp, file = "GSE90166.csv")

#getting metadata
mdata <- GEOquery::getGEO("GSE90166")
mdata <- mdata[[1]]

mdata$description <- c(rep("normal.gravity",3), rep("micro.gravity",3))

source("microarray_functions.R")
metaName <- "GSE90166_meta"
strain <- "UA159"

gselist <- list(list(GSE = mdata))
metalabels <- c("Smut")
extractMetaData(gse_groups = gselist, filename = metaName, microgravity_type = M.TYPE$HARV, metaLabels = metalabels, strain = strain)


