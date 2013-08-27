setwd("/home/comrade/Ubuntu One/DSEA/r-code")

library(grid)
library(gplots)
library(ggplot2)
library(reshape2)
require(Nozzle.R1)

 load("../datasets/sanger.RData")
 load("../datasets/oregon.RData")
# load("../datasets/merged_chembl.RData")

### Import Data Sets and Annotations ###
# ==================================== #
# merged_chembl_ATC_SCID -> chembl.ANNO  # Chembl Drug Annotations
# 
data.frame(exprs(sangerSet)) -> sanger.DATA
data.frame(exprs(xprOregonSet)) -> oregon.DATA
# 
# read.csv(file="../datasets/fimm_chembl_dict.csv", head=TRUE, sep=",") -> fimm.DICT
read.csv(file="../datasets/sanger_chembl_dict.csv", head=TRUE, sep=",") -> sanger.DICT
read.csv(file="../datasets/oregon_chembl_dict.csv", head=TRUE, sep=",") -> oregon.DICT

pData(sangerSet) -> sanger.sample.ANNO
pData(xprOregonSet) -> oregon.sample.ANNO

# cleanup
# remove(sangerSet, xprOregonSet, merged_chembl_ATC_SCID)


### Data Preparetion  ###
# ===================== #

# Snager Set:
sanger.DATA <- sanger.DATA[1:138,]  # Take only the first 138 drugs (those with IC_50 values)
rownames(sanger.DATA) <- gsub("_IC_50$", "", rownames(sanger.DATA))  # Trim "_IC_50$" from drug names 

# Oregon Set:
colnames(oregon.DATA) <- rownames(oregon.sample.ANNO)  # Trim "X" from sample name




################################
###  Below is very old code  ###
################################

load("../datasets/sanger.RData")

# arrange data
drug_names <- featureNames(sangerSet)
sample_names <- sampleNames(sangerSet)
sample_attrib <- colnames(pData(sangerSet))
pd <- pData(sangerSet)

# plot cell line correlations
cell_viability <- data.frame(exprs(sangerSet))
sset <- cell_viability[1:55,1:55]
qplot(x=Var1, y=Var2, data=melt(cor(sset, use="pairwise.complete.obs")), fill=value, geom="tile") + opts(axis.text.x = theme_text(angle = 90))

ind <- c()
for (i in seq(1:ncol(sset))) {
  if ( sum (is.na(sset[,i])) == nrow(sset) )
    ind <- c(ind, i)
}

sset <- sset[,-ind]  # del all-NAs columns
heatmap(as.matrix(sset, na.rm=FALSE, Colv=F, scale='none'))

# plot drug correlations
cell_viability <- data.frame(t(exprs(sangerSet)))
sset <- cell_viability[1:55,1:55]
qplot(x=Var1, y=Var2, data=melt(cor(cell_viability[1:50,1:50], use="na.or.complete")), fill=value, geom="tile")+ opts(axis.text.x = theme_text(angle = 90))

ind <- c()
for (i in seq(1:nrow(sset))) {
  if ( sum (is.na(sset[i,])) == ncol(sset) )
    ind <- c(ind, i)
}

sset <- sset[-ind,]  # del all-NAs columns
heatmap(as.matrix(sset, na.rm=FALSE, Colv=F, scale='none'))

# 1. Hierarchical for all the drugs together
cell_viability <- t(exprs(sangerSet))
cell_viability <- data.frame(sangerSet)

viability_matrix <- t(as.matrix(cell_viability))[1:55,] # create matrix where drugs are in rows
dist <- dist(as.matrix(viability_matrix))
hcl <- hclust(dist)
plot(hcl)