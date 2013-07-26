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


### Leukemia DSS values:       ###
# Drugs: 240; Samples: 32
# Values: 7680;  NAs: 74
# ============================== #

### Breast DSS values:         ###
# Drugs: 239; Samples: 16
# Values: 3824; NAs: 0
# ============================== #


### Data Exploration:  ###
# ====================== #

# Convert the data to 'double matrix' and: 
breast.matrix <- breast.DATA  # copy
breast.matrix <- data.matrix(breast.matrix[,-1])  # del 1st row and transform
breast.transponsed <- t(breast.matrix)
colnames(breast.transponsed) <- breast.DATA[,1]  # assign colnames with drug names
rownames(breast.matrix) <- breast.DATA[,1]  # assign colnames with drug names

# Plot Drug Correlations:
#qplot(x=Var1, y=Var2, data=melt(cor(breast.transponsed[,1:50], use="na.or.complete")), fill=value, geom="tile") + opts(axis.text.x = theme_text(angle = 90))

# HeatMap of the Response: modify color.map to highlight Drug Classes
# color.map <- function(mol.biol) { if (mol.biol=="ALL1/AF4") "#FF0000" else "#0000FF" }
# patientcolors <- unlist(lapply(esetSel$mol.bio, color.map))
heatmap.2(as.matrix(breast.transponsed, na.rm=FALSE, Colv=F, scale='none'), col=topo.colors(75), Rowv=FALSE, na.rm=TRUE,
                                    key=TRUE, symkey=FALSE, trace="none", cexRow=0.5)  # ColSideColors=patientcolors

# Full Set Density
plot( density( breast.transponsed, na.rm=TRUE), main = "Full Set", xlab = "DSS" )

# Boxplot across drugs
boxplot(breast.transponsed, col=1:length(breast.transponsed), xlab = "A Drug", ylab = "Value Distribution", names=NA)
text(seq(1, 239, by=1), -3, srt = 90, pos = 1, xpd = TRUE, labels=breast.DATA[,1], cex=0.5)

# CLUSTERING:
# Generates dendogrames
hr <- hclust(as.dist(1-cor(breast.transponsed, method="pearson")), method="complete")
hc <- hclust(as.dist(1-cor(breast.matrix, method="spearman")), method="complete")  

# Cuts the tree and creates color vector for clusters
mycl <- cutree(hr, h=max(hr$height)/1.5); mycolhc <- rainbow(length(unique(mycl)), start=0.1, end=0.9); mycolhc <- mycolhc[as.vector(mycl)] 


### Ovarian DSS values: !!!     ###
# Drugs: 276 (4 of those do not have a FIMM ID) 
# Samples: 26; Values: 7176;  NAs: 214
# ============================== #


y <- matrix(rnorm(500), 100, 5, dimnames=list(paste("g", 1:100, sep=""), paste("t", 1:5, sep=""))) # Creates a sample data set.
hr <- hclust(as.dist(1-cor(t(y), method="pearson")), method="complete"); hc <- hclust(as.dist(1-cor(y, method="spearman")), method="complete")  
# Generates row and column dendrograms.
mycl <- cutree(hr, h=max(hr$height)/1.5); mycolhc <- rainbow(length(unique(mycl)), start=0.1, end=0.9); mycolhc <- mycolhc[as.vector(mycl)] 
# Cuts the tree and creates color vector for clusters.
library(gplots); myheatcol <- redgreen(75) 
# Assign your favorite heatmap color scheme. Some useful examples: colorpanel(40, "darkblue", "yellow", "white"); 
# heat.colors(75); cm.colors(75); rainbow(75); redgreen(75); library(RColorBrewer); rev(brewer.pal(9,"Blues")[-1]). 
# Type demo.col(20) to see more color schemes.
heatmap.2(y, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc), col=myheatcol, scale="row", density.info="none", trace="none", RowSideColors=mycolhc) 
# Creates heatmap for entire data set where the obtained clusters are indicated in the color bar.
x11(height=6, width=2); names(mycolhc) <- names(mycl); barplot(rep(10, max(mycl)), col=unique(mycolhc[hr$labels[hr$order]]), horiz=T, names=unique(mycl[hr$order])) 
# Prints color key for cluster assignments. The numbers next to the color boxes correspond to the cluster numbers in 'mycl'.
clid <- c(1,2); ysub <- y[names(mycl[mycl%in%clid]),]; hrsub <- hclust(as.dist(1-cor(t(ysub), method="pearson")), method="complete") 
# Select sub-cluster number (here: clid=c(1,2)) and generate corresponding dendrogram.
x11(); heatmap.2(ysub, Rowv=as.dendrogram(hrsub), Colv=as.dendrogram(hc), col=myheatcol, scale="row", density.info="none", trace="none", RowSideColors=mycolhc[mycl%in%clid]) 
# Create heatmap for chosen sub-cluster.
data.frame(Labels=rev(hrsub$labels[hrsub$order])) 
# Print out row labels in same order as shown in the heatmap.

