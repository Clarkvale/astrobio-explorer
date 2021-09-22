
#GSE4136 YEAST MICROARRAY
#BENJAMIN CLARK

library(GEOquery)
library(Biobase)
library(limma)
library(stats)
library(dendextend)
library(gplots)
library(dplyr)
source("microarray_functions.R")
#This pulls the all the samples from the microarry dataset. This returns a list containing a single expressionSet object. 
list.gse <- getGEO("GSE4136", GSEMatrix =TRUE, AnnotGPL=TRUE)

gse <- list.gse[[1]]

#take a look into the dataset
#head(gse)
#View(gse$title)
#dim(gse)

#here we format the feature names. fvarLabels belongs the biobase package and is used to extract features from ExpressionSet Objects
fvarLabels(gse) <- make.names(fvarLabels(gse))

#log2 transform
ex <- exprs(gse)
ex[which(ex <= 0)] <- NaN
exprs(gse) <- log2(ex)

levs <- c(rep("NG.5", 3), rep("NG.25",3),rep("MG.5",3), rep("MG.25",3))
targets <- data.frame(cbind(GSE = gse$geo_accession, Target = levs))
f <- factor(targets$Target, levels = unique(levs))
design <- model.matrix(~0+f)
colnames(design) <- unique(levs)
fit <- lmFit(gse, design)



#5th gen
cont.dif.5 <- makeContrasts(
  Dif5 = MG.5 - NG.5,
  
  
  levels = design
)
contrasts.fit <- contrasts.fit(fit, contrasts = cont.dif.5)

fit2 <- eBayes(contrasts.fit, 0.01)
g.tT5 <- topTable(fit2, adjust.method = "fdr", confint = TRUE, number = Inf)

gen5name <- "datasets/GSE4136_Scer/GSE4136_5thGen.csv"
write.table(pull.output.tT(g.tT5), gen5name, row.names = FALSE, sep = ",")



#25th gen
cont.dif.25 <- makeContrasts(
  Dif25 = MG.25 - NG.25,
  levels = design
)

contrasts.fit <- contrasts.fit(fit, contrasts = cont.dif.25)
fit2 <- eBayes(contrasts.fit, 0.01)
g.tT25 <- topTable(fit2, adjust.method = "fdr", confint = TRUE, number = Inf)

gen25name <- "datasets/GSE4136_Scer/GSE4136_25thGen.csv"
write.table(pull.output.tT(g.tT25), gen25name, row.names = FALSE, sep = ",")




metaName <- "datasets/GSE4136_Scer/GSE4136_meta"
strain <- "BY4743"

labels <- c("5thGen", "25thGen")
extractMetaData(gse = gse, metaLabels = labels, 
                microgravity_type = M.TYPE$SPACEFLOWN, design = design, 
                contrasts = list(cont.dif.5, cont.dif.25), strain = strain,
                filename = metaName)














#5th gen groups
control <- c(1,2,3)
treatment <- c(7,8,9)
gen5 <- de.analysis(gse = gse, microgravity_group = treatment, ground_group = control)

#print out the toptable 
gen5name <- "datasets/GSE4136_Scer/GSE4136_5thGen.csv"
write.table(gen5$TopTable, gen5name, row.names = FALSE, sep = ",")



#25th gen groups
control <- c(4,5,6)
treatment <- c(10,11,12)

gen25 <- de.analysis(gse = gse, microgravity_group = treatment, ground_group = control)
gen25name <- "datasets/GSE4136_Scer/GSE4136_25thGen.csv"
write.table(gen25$TopTable, gen25name, row.names = FALSE, sep = ",")


#process and extract metadata for all datasets comparisons

metaName <- "datasets/GSE4136_Scer/GSE4136_meta"
strain <- "BY4743"
gse_list <- list(gen5,gen25)
labels <- c("5thGen", "25thGen")
extractMetaData(filename = metaName, gse_groups = gse_list, microgravity_type = M.TYPE$HARV, metaLabels = labels, strain = strain)


#
#_____________________ESSENTIAL GENES FOR SPACEFLIGHT__________________________


homo14 <- readxl::read_xlsx("essential/TableS5_vsT1-hom1-dropsToBg.xlsx", sheet = "flight-ground 14G")
homo21 <- readxl::read_xlsx("essential/TableS5_vsT1-hom1-dropsToBg.xlsx", sheet = "flight-ground 21G")

homo_names <- union(homo14$Gene, homo21$Gene)



hetero <- readxl::read_xlsx("essential/TableS8_vsT1-het1-dropsToBg.xlsx", sheet = "flight-ground 21G")$Gene

full_names <- union(homo_names, hetero)


#getting genes excluded from study
ex_genes_g <- readxl::read_xlsx("essential/TableS3_excludedStrains-hom.xlsx", sheet = c("ground"))$Gene
ex_genes_f <- readxl::read_xlsx("essential/TableS3_excludedStrains-hom.xlsx", sheet = c("flight"))$Gene

ex_genes <- union(ex_genes_f, ex_genes_g)

#Adding to Gen5 Study 
essential_bool5 <- c()

essential_bool5[match(full_names, gen5$TopTable$Gene.symbol)] <- TRUE
essential_bool5[match(ex_genes, gen5$TopTable$Gene.symbol)] <- NA

listed_genes <- union(ex_genes, full_names)
falsed <- setdiff(gen5$TopTable$Gene.symbol, listed_genes)

essential_bool5[match(falsed, gen5$TopTable$Gene.symbol)] <- FALSE

gen5$TopTable <-  gen5$TopTable %>% mutate(Essential_For_SpaceFlight = essential_bool5)

#Adding to Gen25 Study
essential_bool25 <- c()

essential_bool25[match(full_names, gen25$TopTable$Gene.symbol)] <- TRUE
essential_bool25[match(ex_genes, gen25$TopTable$Gene.symbol)] <- NA

listed_genes <- union(ex_genes, full_names)
falsed <- setdiff(gen25$TopTable$Gene.symbol, listed_genes)

essential_bool25[match(falsed, gen25$TopTable$Gene.symbol)] <- FALSE

gen25$TopTable <-  gen25$TopTable %>% mutate(Essential_For_SpaceFlight = essential_bool25)


#________________________________________________________________________________

##DE analysis Between generations
micro.5thgen <- c(7,8,9)
micro.25thgen <- c(10,11,12)
gen5.vs.gen25 <- de.analysis(gse = gse, microgravity_group = micro.25thgen, ground_group = micro.5thgen)

failed.probes <- remove.controls(gen5.vs.gen25$TopTable)
##CLUSTERING ANALYSIS


#FILTERING SIGNIFICANT DE GENES
gen5.sig.df.genes <- which(gen5$TopTable$logFC >= 1 & gen5$TopTable$adj.P.Val <= 0.01)
gen5.orfs <- gen5$TopTable$ID[gen5.sig.df.genes]

gen25.sig.df.genes <- which(gen25$TopTable$logFC >= 1 & gen25$TopTable$adj.P.Val <= 0.01)
gen25.orfs <- gen25$TopTable$ID[gen25.sig.df.genes]

grouped.orfs <- union(gen5.orfs,gen25.orfs)

#pulling relevent expression data 
ex <- exprs(gse)
filtered.ex <- ex[rownames(ex) %in% grouped.orfs,]

#take out control columns
filtered.ex <- filtered.ex[,7:12]

#adding our factors
both <- intersect(gen5.orfs,gen25.orfs)
gen5.dif <- setdiff(gen5.orfs, gen25.orfs)
gen25.dif <- setdiff(gen25.orfs,gen5.orfs)


filtered.ex <- log2(filtered.ex)

#Making plots
colnames(filtered.ex) <- c("5gen1", "5gen2", "5gen3", "25gen1", "25gen2", "25gen3")

dev.new(width = 6, height = 256)
plot(heatmap.2(filtered.ex), distfun = function(x) as.dist(1-cor(t(x), method="pearson")))

hc <- as.dendrogram(hclust(as.dist(1-cor(t(filtered.ex), method="pearson")), method="complete"))


#c.index <- which(labels(hc) %in% both)

hc %>% set("labels_col", "white") %>% 
  set("by_labels_branches_col", value = both) %>% 
  plot(main = "Highlighted Union Genes Across Timescales")

#PLOTTING INTERSECTING Genes
filtered.ex.intersect <- ex[rownames(ex) %in% both,][,7:12]
colnames(filtered.ex.intersect) <- c("5gen1", "5gen2", "5gen3", "25gen1", "25gen2", "25gen3")
plot(heatmap(filtered.ex.intersect))


#build annotable
symbols <- gen5$TopTable$Gene.symbol[gen5$TopTable$ID %in% both]
logfc <- gen5$TopTable$logFC[gen5$TopTable$ID %in% both]
p.val <- gen5$TopTable$adj.P.Val[gen5$TopTable$ID %in% both]
