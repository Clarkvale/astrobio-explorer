source("microarray_functions.R")
library(limma)

gse <- getGEOFromRMA("GSE23880")
gse.slim <- gse[,c(c(1:6), c(10:15))]

levs  <- c(rep("MG", 6), rep("NG", 6))

targets <- data.frame(cbind(GSE = gse.slim$geo_accession, Target = levs))
f <- factor(targets$Target, levels = unique(levs))


design <- model.matrix(~0+f)
colnames(design) <- unique(levs)


fit <- lmFit(gse.slim, design)

cont <- makeContrasts(
  Dif =  MG - NG,
  levels = design
)


contrasts.fit <- contrasts.fit(fit, contrasts = cont)

fit2 <- eBayes(contrasts.fit, 0.01)

#Here we pull a toptable which will return all our logfc values and (hopefully) annotation data
tT<- topTable(fit2, adjust.method = "fdr", confint = TRUE, number = Inf)

tT <-  pull.output.tT(tT)

write.table(tT, "datasets/GSE23880_Dmel/GSE23880.csv", row.names = FALSE, sep = ",")

extractMetaData(gse = gse.slim, design = design, filename = "datasets/GSE23880_Dmel/GSE23880_meta",
                contrasts = cont, microgravity_type = M.TYPE$SPACEFLOWN, strain = gse.slim$`strain:ch1`[[1]],
                )
