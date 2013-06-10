setwd("/home/comrade/Ubuntu One/DSEA")
library(ggplot2)
library(grid)

source('testing-related.r')
load("oregon.RData")

# arrange data
drug_names <- featureNames(xprOregonSet)
sample_names <- sampleNames(xprOregonSet)
sample_attrib <- colnames(pData(xprOregonSet))
pd <- pData(xprOregonSet)

# extract cell viability data
cell_viability <- t(exprs(xprOregonSet))
cell_viability <- data.frame(cell_viability)

# ! check sensitibity bar plot - JP's script that shows drug sensitivuty ranged list ! 
# ! AND Drug Correlations also

### CLUSTERING ###


# 1. Hierarchical for all the drugs together
viability_matrix <- t(as.matrix(cell_viability))[1:55,] # create matrix where drugs are in rows
dist <- dist(as.matrix(viability_matrix))
hcl <- hclust(dist)
plot(hcl)