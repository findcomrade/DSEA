setwd("/home/comrade/Ubuntu One/DSEA/r-code")

library(grid)
library(ggplot2)
library(reshape2)

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