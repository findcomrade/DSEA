setwd("/home/comrade/Ubuntu One/DSEA/r-code")

library(grid)
library(gplots)
library(ggplot2)
library(reshape2)

load("../datasets/sanger.RData")
load("../datasets/oregon.RData")
load("../datasets/merged_chembl.RData")

### Import Data Sets and Annotations ###
# ==================================== #
merged_chembl_ATC_SCID -> chembl.ANNO  # Chembl Drug Annotations

data.frame(exprs(sangerSet)) -> sanger.DATA
data.frame(exprs(xprOregonSet)) -> oregon.DATA

read.csv(file="../datasets/all_leukemia_cl_june_13_disha_astrid.csv", head=TRUE, sep=",") -> leukemia.DATA
read.csv(file="../datasets/original_dss_june_13_ovarian_akira_astrid.csv", head=TRUE, sep=",") -> ovarian.DATA
read.csv(file="../datasets/original_dss_breast_cancer_tamoxifen_project_sk_sh.csv", head=TRUE, sep=",") -> breast.DATA

read.csv(file="../datasets/fimm_chembl_dict.csv", head=TRUE, sep=",") -> fimm.DICT
read.csv(file="../datasets/sanger_chembl_dict.csv", head=TRUE, sep=",") -> sanger.DICT
read.csv(file="../datasets/oregon_chembl_dict.csv", head=TRUE, sep=",") -> oregon.DICT

pData(sangerSet) -> sanger.sample.ANNO
pData(xprOregonSet) -> oregon.sample.ANNO

#DB <- list(chembl.ANNO, oregon.DATA, sanger.DATA, leukemia.DATA, ovarian.DATA, breast.DATA,
#                            fimm.DICT, sanger.DICT, oregon.DICT, sanger.sample.ANNO,oregon.sample.ANNO)
# cleanup
remove(sangerSet, xprOregonSet, merged_chembl_ATC_SCID)


### Data Preparetion  ###
# ===================== #

# Snager Set:
sanger.DATA <- sanger.DATA[1:138,]  # Take only the first 138 drugs (those with IC_50 values)
rownames(sanger.DATA) <- gsub("_IC_50$", "", rownames(sanger.DATA))  # Trim "_IC_50$" from drug names 

# Oregon Set:
colnames(oregon.DATA) <- rownames(oregon.sample.ANNO)  # Trim "X" from sample name

# Ovarian Set:
# ovarian.DATA[duplicated(ovarian.DATA[,2]),2] =>
# The following drugs are duplicated in the initial table: "Vinorelbine" "Cytarabine" "Clofarabine" "Bortezomib"
# Corresponding rows can be merged, cause they complement each other:
# ovarian.DATA[ovarian.DATA[,2] == "Vinorelbine", 1:27]
ovarian.DATA[179,25:27] <- ovarian.DATA[28,25:27]  # Vinorelbine
ovarian.DATA[189,25:27] <- ovarian.DATA[38,25:27]  # Cytarabine
ovarian.DATA[198,25:27] <- ovarian.DATA[25,25:27]  # Clofarabine
ovarian.DATA[200,25:27] <- ovarian.DATA[19,25:27]  # Bortezomib
ovarian.DATA <- ovarian.DATA[-c(19,25,28,38),]  # del duolicated rows without IDs

# assign "NA" to empty FIMM_ID fields (there are 8 empty fields) - THIS IS DEPRICATED:
# Missing IDs refer to duplicated rows (drugs) that should be merged. See section above!
# ovarian.DATA[-grep( "FIMM", as.character(ovarian.DATA[,1]), ignore.case = FALSE ), 1] <- NA

### Produce FIMM Annotations  ###
# ============================= #
# Extract relevant information from 'ChEMBL db'
# to annotate drugs in FIMM data base.

# drop drugs without chembl_ID in FIMM Dctionary
drop <- which(is.na(fimm.DICT$CHEMBL_ID))
fimm.DICT <- fimm.DICT[-drop,]
fimm.drug.count <- length( unique(fimm.DICT$FIMM.batch.ID) )  # count drugs in FIMM Collection

# Build annotation table for FIMM Collection
# "Level_3_Description" is missing in "chembl.ANNO"
drop <-  chembl.ANNO[,"CHEMBL_ID"] %in% unique(fimm.DICT[,"CHEMBL_ID"]) 
fimm.ANNO <- chembl.ANNO[drop,c("CHEMBL_ID","Level1","Level2","Level3","Level4","Level5","Level_1_Description","Level_2_Description","Level_4_Description")]
drop <- which(is.na(fimm.ANNO$Level1))  # drop rows without annotations
fimm.ANNO <- fimm.ANNO[-drop,]
fimm.ANNO <- unique(fimm.ANNO)  # remove duplicated rows

remove(chembl.ANNO, drop)


#########################################
###          DATA EXPLORATION       ####
#########################################

### Breast DSS values:         ###
# Drugs: 239; Samples: 16
# Values: 3824; NAs: 0
# ============================== #

# Convert the data to 'double matrix' and: 
breast.matrix <- breast.DATA  # copy
breast.matrix <- data.matrix(breast.matrix[,-1])  # del 1st row and transform
breast.transponsed <- t(breast.matrix)
colnames(breast.transponsed) <- breast.DATA[,1]  # assign colnames with drug names
rownames(breast.matrix) <- breast.DATA[,1]  # assign colnames with drug names

# There are 58 Drugs with DSS == 0 for each cell line
# We should trim those in order to run proper clustering
# Drug List:  "  breast.DATA[drop,1] "
drop <- which(apply(breast.matrix,1,sum) == 0)
breast.matrix <- breast.matrix[-drop,]
breast.transponsed <- breast.transponsed[,-drop]


# Summary
drop <- fimm.DICT$Aliases %in% breast.DATA$Name.Drug  
breast.DICT <- fimm.DICT[drop,]
drop <- fimm.ANNO$CHEMBL_ID %in% breast.DICT$CHEMBL_ID
breast.ANNO <- fimm.ANNO[drop,]

breast.summary <- data.frame( Drugs.Count=length( unique(breast.DATA$Name.Drug) ), 
                              Drugs.Annotated=length( unique(breast.ANNO$CHEMBL_ID) ),
                              Annotations.Count=length( breast.ANNO$CHEMBL_ID ),
                              L1.Levels=length( unique(breast.ANNO$Level1) ),
                              L2.Levels=length( unique(breast.ANNO$Level2) ),
                              L3.Levels=length( unique(breast.ANNO$Level3) ),
                              L4.Levels=length( unique(breast.ANNO$Level4) ),
                              L5.Levels=length( unique(breast.ANNO$Level5) ) 
                                                )

x11()
par(mfrow=c(2,2))
for(k in 1:4){
  count <- table(breast.ANNO[,k+1])
  barplot(count, col=1:length(count), main=paste("Category Distribution on Level ", k, sep=""), xlab="\"Breast Data Set\"")
}
  
# Plot Drug Correlations:
#qplot(x=Var1, y=Var2, data=melt(cor(breast.transponsed[,1:50], use="na.or.complete")), fill=value, geom="tile") + opts(axis.text.x = theme_text(angle = 90))

# HeatMap of the Response: modify color.map to highlight Drug Classes
# color.map <- function(mol.biol) { if (mol.biol=="ALL1/AF4") "#FF0000" else "#0000FF" }
# patientcolors <- unlist(lapply(esetSel$mol.bio, color.map))
heatmap.2(as.matrix(breast.matrix, na.rm=FALSE, Colv=F, scale='none'), col=topo.colors(75), Rowv=TRUE, na.rm=TRUE,
                                    key=TRUE, symkey=FALSE, trace="none", cexRow=0.5, cexCol=0.5)  # ColSideColors=patientcolors

# Full Set Density
plot( density( breast.transponsed, na.rm=TRUE), main = "Full Set", xlab = "DSS" )

# Boxplot across drugs
boxplot(breast.transponsed, col=1:length(breast.transponsed), xlab = "A Drug", ylab = "Value Distribution", names=NA)
text(seq(1, length(rownames(breast.matrix)), by=1), -3, srt = 90, pos = 1, xpd = TRUE, labels=rownames(breast.matrix), cex=0.5)

# CLUSTERING:
# ==============

# (By "Pearson correlation as distance method")
# Alailable methods for HCLUST: "ward", "single", "complete", "average", "mcquitty", "median" or "centroid"
hr <- hclust(as.dist(1-cor(breast.transponsed, method="pearson")), method="complete")
hc <- hclust(as.dist(1-cor(breast.matrix, method="spearman")), method="complete")  

# Cuts the tree and creates color vector for clusters
mycl <- cutree(hr, h=max(hr$height)/1.5); 
clusters.frame <- as.data.frame(mycl)  # to use 'clusters.frame' in further analysis
mycolhc <- rainbow(length(unique(mycl)), start=0.1, end=0.9); mycolhc <- mycolhc[as.vector(mycl)] 
myheatcol <- topo.colors(75)

# Creates heatmap for entire data set where the obtained clusters are indicated in the color bar
heatmap.2(breast.matrix, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc), col=myheatcol, 
                  scale="row", trace="none", RowSideColors=mycolhc, cexRow=0.5, cexCol=0.5) 

# PROCESS OBTAINED CLUSTERS from "mycl"
clusters.frame$DrugName <- rownames(clusters.frame)
colnames(clusters.frame) <- c("Cluster", "DrugName")

# Prints color key for cluster assignments. The numbers next to the color boxes correspond to the cluster numbers in 'mycl'
# x11(height=6, width=2); names(mycolhc) <- names(mycl); barplot(rep(10, max(mycl)), col=unique(mycolhc[hr$labels[hr$order]]), horiz=T, names=unique(mycl[hr$order])) 

# Select sub-cluster number (here: cluster.id=c(1,2)) and generate corresponding dendrogram.
cluster.id <- c(1,2); ysub <- breast.matrix[names(mycl[mycl %in% cluster.id]),]
hrsub <- hclust(as.dist(1-cor(t(ysub), method="pearson")), method="complete") 

heatmap.2(ysub, Rowv=as.dendrogram(hrsub), Colv=as.dendrogram(hc), col=myheatcol, 
                  scale="row", trace="none", RowSideColors=mycolhc[mycl%in%clid]) 

# Print out row labels in same order as shown in the heatmap.
data.frame(Labels=rev(hrsub$labels[hrsub$order])) 


### Ovarian DSS values: !!!     ###
# Drugs: 276 (4 of those do not have a FIMM ID) 
# Samples: 26; Values: 7176;  NAs: 214
# ============================== #


### Leukemia DSS values:       ###
# Drugs: 240; Samples: 32
# Values: 7680;  NAs: 74
# ============================== #


