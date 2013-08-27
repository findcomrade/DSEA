setwd("/home/comrade/Ubuntu One/DSEA/r-code")
library(ggplot2)
library(grid)

source('testing-related.r')

read.csv(file="../datasets/updated_22_aug/AML_updated.csv", head=TRUE, sep="\t") -> aml.DATA
read.csv(file="../datasets/drug_dictionary_FIMM.csv", head=TRUE, sep="\t") -> fimm.DICT
read.csv(file="../datasets/updated_22_aug/fimm_collection.csv", head=TRUE, sep="\t") -> fimm.UPDATE


fimm.UPDATE <- unique(fimm.UPDATE)


aml.DATA[,"ID.Drug"] <- strtrim(aml.DATA[,"ID.Drug"],10)  # trim batch.id
fimm.DICT$ID.Drug <- strtrim(fimm.DICT[,"FIMM.batch.ID"],10)

collection.FIMM <- cbind(fimm.DICT[,5], fimm.DICT[,1:4])
colnames(collection.FIMM) <- c("FIMM.ID", "FIMM.Batch.ID", "Drug.Name", "ChEMBL.ID", "PubChem.CID")

collection.FIMM <- unique(collection.FIMM)

sum( unique(fimm.UPDATE[,"FIMM.ID"]) %in% unique(collection.FIMM[,"FIMM.ID"]) )  # 256
sum(collection.FIMM[,"FIMM.ID"] %in% fimm.UPDATE[,"FIMM.ID"] )  # 256
sum( unique(fimm.UPDATE[,"FIMM.BATCH.ID"]) %in% unique(collection.FIMM[,"FIMM.Batch.ID"]) )  # 201
sum( unique(fimm.UPDATE[,"Drug.Name"]) %in% unique(collection.FIMM[,"Alias"]) )  # 242

# bind Targets to collection
collection.FIMM$Target <- "Undefined"
wh <- which( collection.FIMM[,"FIMM.ID"] %in% unique(fimm.UPDATE[,"FIMM.ID"]) )
for(row in wh){
  id <- as.character( collection.FIMM[row,"FIMM.ID"] )
  collection.FIMM[row, "Target"] <- as.character( unique( fimm.UPDATE[ fimm.UPDATE[,"FIMM.ID"] == id, "Mechanism.Targets"] ))
}

save(collection.FIMM, file = "RData/FimmCollection.RData")

sum( aml.DATA[,"ID.Drug"] %in% fimm.UPDATE[,"FIMM.ID"] ) 
sum( aml.DATA[,"ID.Drug"] %in% fimm.DICT[,"FIMM.batch.ID"] ) 

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