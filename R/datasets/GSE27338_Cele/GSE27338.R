#Benjamin Clark 
library(GEOquery)
library(Biobase)
library(limma)
library(dplyr)
source("microarray_functions.R")
library(org.Ce.eg.db)

gse <- getGEO("GSE27338", AnnotGPL = T)
gse <- gse$GSE27338_series_matrix.txt.gz


exprs(gse)[which( exprs(gse) <= 0)] <- NaN
exprs(gse) <- log2(exprs(gse))

fvarLabels(gse) <- make.names(fvarLabels(gse))

#This is my levels array, just a short list of descriptors for each sample. Notice how I'm including both "NG" and "MG" for "normal gravity" and #vice versa and the generation number.
grav <- c(rep("ANG", 3), rep("MicroG",3),rep("NG",3), rep("ANG",3),rep("MicroG",3), rep("NG",3))
block <- c(rep(2, 9), rep(1,9))

#this is a way of recreating a very literal design matrix that we'll use later
targets <- data.frame(cbind(GSE = gse$geo_accession, Target = grav, Block = block))


#grab the factors from our targets dataframe
f <- factor(targets$Target, levels = unique(grav))
b <-  factor(targets$Block, levels = unique(block))
#build a proper design matrix that limma can use
design <- model.matrix(~b+f)



#fit the data to a linear model
fit <- lmFit(gse, design)


#Now we need to ask limma to draw contrasts between the groups we want to compare. It is important we get the direction right, with the treated #sample VS the null samples  

#5th gen, the names have to match our design names we made before
cont.dif <- makeContrasts(
  Dif = fMicroG - fNG,
  levels = design
)

#fit the contrasts now
contrasts.fit <- contrasts.fit(fit, contrasts = cont.dif)

#limma can use naive bayes reduce p values even further by borrowing neighboring values.
fit2 <- eBayes(contrasts.fit, 0.01)

#Here we pull a toptable which will return all our logfc values and (hopefully) annotation data
tT <- topTable(fit2, adjust.method = "fdr", confint = TRUE, number = Inf)

 

tT2 <- agilent2Affy(tT)

gos <- get.GOs(tT2, org.Ce.eg.db, "Entrez.ID")


tT.out <- inner_join(tT2, gos, by = "Entrez.ID")



 
#Save our data to a flat csv file. NOTE the function pull.output.tT() ONLY works on affymetrix microarrays. You'll have pull the info manually #otherwise. I suggest using dplyr
#Since we are still in the project directory we still need to walk down the repo to the sub directory the R file is in to save it.
name <- "datasets/GSE27338_Cele/GSE27338.csv"

write.csv(tT.out, file = name)

extractMetaData(gse = gse, design = design, contrasts = cont.dif, 
                microgravity_type =  M.TYPE$SPACEFLOWN, 
                filename = "datasets/GSE27338_Cele/GSE27338_meta"
                , strain = gse$`strain:ch1`[[1]])






