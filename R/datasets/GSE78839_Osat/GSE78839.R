library(affy)
library(GEOquery)
library(limma)
library(dplyr)
source("microarray_functions.R")

celpath = "datasets/GSE78839_Osat/raw"
data = ReadAffy(celfile.path=celpath)

data.rma <- rma(data)
ex.norm <- exprs(data.rma)


gse <- getGEO("GSE78839", GSEMatrix = T, AnnotGPL = T)[[1]]


colnames(ex.norm) <- gse$geo_accession

exprs(gse) <- ex.norm

levs <- make.names(sapply(gse$title, function(x){return(substr(x, 1, nchar(x) - 17))}, USE.NAMES = F))
targets <- data.frame(cbind(GSE = gse$geo_accession, Target = levs))
f <- factor(targets$Target, levels = unique(levs))


design <- model.matrix(~0+f)
colnames(design) <- unique(levs)

#fit the data to a linear model
fit <- lmFit(gse, design)

cont <- makeContrasts(
  Dif = rice.calli.at.spaceflight.control - rice.calli.at.1g.ground.control,
  levels = design
)


contrasts.fit <- contrasts.fit(fit, contrasts = cont)

#limma can use naive bayes reduce p values even further by borrowing neighboring values.
fit2 <- eBayes(contrasts.fit, 0.01)

#Here we pull a toptable which will return all our logfc values and (hopefully) annotation data
tT<- topTable(fit2, adjust.method = "fdr", confint = TRUE, number = Inf)

tT <-  pull.output.tT(tT)
tT <- tT %>%  dplyr::filter(Entrez.ID != "")

write.table(tT, "datasets/GSE78839_Osat/GSE78839.csv", row.names = FALSE, sep = ",")
extractMetaData(gse, design, contrasts = cont, filename = "datasets/GSE78839_Osat/GSE78839_meta", microgravity_type = M.TYPE$SPACEFLOWN)

