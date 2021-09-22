library(GEOquery)
library(Biobase)
library(limma)
library(dplyr)
source("microarray_functions.R")

gse <- getGEO("GSE94381", GSEMatrix =TRUE, AnnotGPL=TRUE)[[1]]

logcheck(exprs(gse))
boxplot(exprs(gse))
plotDensities(exprs(gse))

fvarLabels(gse) <- make.names(fvarLabels(gse))

levs <- sapply(make.names(gse$title), FUN = function(x) {return(substr(x, 1, nchar(x)-1))}, USE.NAMES = F)
targets <- data.frame(cbind(GSE = gse$geo_accession, Target = levs))

#grab the factors from our targets dataframe
f <- factor(targets$Target, levels = unique(levs))
design <- model.matrix(~0+f)
colnames(design) <- unique(levs)

#fit the data to a linear model
fit <- lmFit(gse, design)


#Now we need to ask limma to draw contrasts between the groups we want to compare. It is important we get the direction right, with the treated #sample VS the null samples  

#5th gen, the names have to match our design names we made before
cont <- makeContrasts(
  Dif5 = BF.LD - BG.LD,
  levels = design
)

contrasts.fit <- contrasts.fit(fit, contrasts = cont)

#limma can use naive bayes reduce p values even further by borrowing neighboring values.
fit2 <- eBayes(contrasts.fit, 0.01)

#Here we pull a toptable which will return all our logfc values and (hopefully) annotation data
tT<- topTable(fit2, adjust.method = "fdr", confint = TRUE, number = Inf)
tT <-  pull.output.tT(tT)


filename<- "datasets/GSE94381_Mmus/GSE94381.csv"


write.table(tT, filename, row.names = FALSE, sep = ",")

metaName <- "datasets/GSE94381_Mmus/GSE94381_meta"


strain <- "C57/BL6"


#one label for each comparison we made
labels <- c("")

#microgravity_type needs to be derived from M.TYPE, design is the output from our model matrix, contrasts is a LIST of separate contrasts made from #make.Contrasts()

extractMetaData(gse = gse, metaLabels = labels, 
                microgravity_type = M.TYPE$SPACEFLOWN, design = design, 
                contrasts = list(cont), strain = strain,
                filename = metaName, cellType = gse$`tissue:ch1`[[1]])