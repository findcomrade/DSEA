setwd("/home/comrade/Ubuntu One/DSEA/r-code")

library(grid)
library(ggplot2)

load("../datasets/sanger.RData")

# arrange data
drug_names <- featureNames(sangerSet)
sample_names <- sampleNames(sangerSet)
sample_attrib <- colnames(pData(sangerSet))
pd <- pData(sangerSet)

cell_viability <- t(exprs(sangerSet))
cell_viability <- data.frame(sangerSet)

heatmap(exprs(sangerSet))

# 1. Hierarchical for all the drugs together
viability_matrix <- t(as.matrix(cell_viability))[1:55,] # create matrix where drugs are in rows
dist <- dist(as.matrix(viability_matrix))
hcl <- hclust(dist)
plot(hcl)