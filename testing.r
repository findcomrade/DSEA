setwd("/home/comrade/Ubuntu One/DSEA")
library(ggplot2)
library(grid)

source('testing-related.r')
# load and explore data
load("oregon.RData")
head(xprOregonSet)
summary(pData(xprOregonSet))

# arrange data
drug_names <- featureNames(xprOregonSet)
sample_names <- sampleNames(xprOregonSet)
sample_attrib <- colnames(pData(xprOregonSet))
pd <- pData(xprOregonSet)

# extract cell viability data from the data set
cell_viability <- t(exprs(xprOregonSet))
cell_viability <- data.frame(cell_viability)

# cell viability across all 66 drugs
data_frame <- collect_drugs(15,17)

# diagnosis VS viability
ggplot(data_frame, aes(diagnosis, viability)) + geom_jitter(size=2, alpha = I(1 / 1.5), aes(color=drug))
# diagnosis VS viability, trim viability
ggplot(data_frame, aes(diagnosis, viability)) + geom_jitter(size=2, alpha = I(1 / 1.5), aes(color=drug)) + ylim(0,10000)

# drug VS viabilty
ggplot(data_frame, aes(drug, viability)) + geom_jitter(size=2, alpha = I(1 / 1.5), aes(color=diagnosis)) + theme(legend.position = "none")
ggplot(data_frame, aes(drug, viability)) + 
  geom_boxplot(aes(fill = factor(diagnosis)))
# ggplot(data_frame, aes(diagnosis, viability)) + geom_point(aes(color = drug), size = 2)
# ggplot(data_frame, aes(diagnosis, viability)) + geom_point(aes(color = drug), size = 2) + ylim(0,10000)

# now we take one drug across 6 conditions
data_frame <- samples_by_drug("AP24534", cell_viability)
ggplot(data_frame, aes(drug, viability)) + 
  geom_boxplot(aes(fill = factor(diagnosis)))



### CLUSTERING ###
# create matrix where drugs are in rows
viability_matrix <- t(as.matrix(cell_viability))[1:55,]
dist <- dist(as.matrix(viability_matrix))   # find distance matrix 
hcl <- hclust(dist)
plot(hcl)

### small test
d1 <- samples_by_drug("Nilotinib", cell_viability)
d2 <- samples_by_drug("CYT387", cell_viability)
d4 <- samples_by_drug("EKB.569", cell_viability)
merged <- rbind(d1, d2, d3, d4)
ggplot(merged, aes(drug, viability)) + 
  geom_boxplot(aes(fill = factor(diagnosis)))

### take samples with certain diagnosis ###
condition_cell_viability <- filter_cell_viability("MATURE B-CELL NEOPLASMS")
# and clustering
viability_matrix <- t(as.matrix(condition_cell_viability))[1:66,]
dist <- dist(as.matrix(viability_matrix))   # find distance matrix 
hcl <- hclust(dist)
plot(hcl)

### small test

d1 <- samples_by_drug("Nilotinib", condition_cell_viability)
d2 <- samples_by_drug("Dasatinib", condition_cell_viability)
d3 <- samples_by_drug("Sunitinib", condition_cell_viability)
d4 <- samples_by_drug("AST.487", condition_cell_viability)
merged <- rbind(d1, d2, d3, d4)

ggplot(merged, aes(drug, viability)) + 
  geom_boxplot(aes(fill = factor(diagnosis)))