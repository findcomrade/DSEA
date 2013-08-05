setwd("/home/comrade/Ubuntu One/DSEA/r-code")
load("../datasets/merged_chembl.RData")

library(vcd)
library(grid)
library(gplots)
library(lattice)
library(ggplot2)
library(reshape2)
require(Nozzle.R1)

REPORT <- newCustomReport( "DSEA:", asEmph(" Ovarian Project.") )
### Import Data Sets and Annotations ###
# ==================================== #
merged_chembl_ATC_SCID -> chembl.ANNO  # Chembl Drug Annotations
read.csv(file="../datasets/original_dss_june_13_ovarian_akira_astrid.csv", head=TRUE, sep=",") -> ovarian.DATA
read.csv(file="../datasets/fimm_chembl_dict.csv", head=TRUE, sep=",") -> fimm.DICT
remove(merged_chembl_ATC_SCID)  # cleanup

### Produce FIMM Annotations  ###
# ============================= #
fimmcollection.section.REPORT <- newSection( "Fimm Collection Annotations" )
atc.subsec.REPORT <- newSection( "Anatomical Therapeutic Chemical Classification System (ATC)" )
atc.subsec.REPORT <- addTo( atc.subsec.REPORT, newParagraph("This pharmaceutical coding system divides drugs into different groups according to the organ or system on which they act and/or their therapeutic and chemical characteristics. Each bottom-level ATC code stands for a pharmaceutically used substance, or a combination of substances, in a single indication (or use). This means that one drug can have more than one code. (Wikipedia)") )
atc.subsec.REPORT <- addTo( atc.subsec.REPORT, newList(
  newParagraph("The ", asStrong("first level"), " of the code indicates the anatomical main group and consists of one letter. There are 14 main groups;"),    
  newParagraph("The ", asStrong("second level"), " of the code indicates the therapeutic main group and consists of two digits;"),
  newParagraph("The ", asStrong("third level"), " of the code indicates the therapeutic/pharmacological subgroup and consists of one letter;"),    
  newParagraph("The ", asStrong("fourth level"), " of the code indicates the chemical/therapeutic/pharmacological subgroup and consists of one letter;"),
  newParagraph("The ", asStrong("fifth level"), " of the code indicates the chemical substance and consists of two digits.")
) 
)
# Extract relevant information from 'ChEMBL db'
# to annotate drugs in FIMM data base.
fimm.drug.count <- length( unique(fimm.DICT$FIMM.batch.ID) )
drop <- which(is.na(fimm.DICT$CHEMBL_ID))
fimm.DICT <- fimm.DICT[-drop,]  # drop drugs without chembl_ID in FIMM Dctionary
fimm.drug.chembl <- length( unique(fimm.DICT$FIMM.batch.ID) )
fimm.SUMMARY <- data.frame(Drugs.Count=fimm.drug.count, ChEMBL.IDs=fimm.drug.chembl)
remove(fimm.drug.chembl, fimm.drug.count)

# Build annotation table for FIMM Collection
# "Level_3_Description" is missing in "chembl.ANNO" !!!
drop <- chembl.ANNO[,"CHEMBL_ID"] %in% unique(fimm.DICT[,"CHEMBL_ID"]) 
fimm.ANNO <- chembl.ANNO[drop,c("CHEMBL_ID","Level1","Level2","Level3","Level4","Level5","Level_1_Description","Level_2_Description","Level_4_Description")]
drop <- which(is.na(fimm.ANNO$Level1))  # drop rows without annotations

fimm.ANNO <- fimm.ANNO[-drop,]
fimm.ANNO <- unique(fimm.ANNO)  # remove duplicated rows
fimm.SUMMARY$CheMBL.Annotated <- length(fimm.ANNO$CHEMBL_ID)

fimmcollection.section.REPORT <- addTo( fimmcollection.section.REPORT, atc.subsec.REPORT )
fimmcollection.section.REPORT <- addTo( fimmcollection.section.REPORT, newTable( fimm.SUMMARY, "FIMM Collection Summary" ) )
fimmcollection.section.REPORT <- addTo( fimmcollection.section.REPORT, newTable( fimm.DICT[1:7,], "FIMM Dictionary" ) )
fimmcollection.section.REPORT <- addTo( fimmcollection.section.REPORT, newTable( fimm.ANNO[1:3,], "FIMM Annotations" ) )
remove(chembl.ANNO, drop)


#########################################
###          DATA EXPLORATION       ####
#########################################

### Ovarian DSS values: !!!     ###
# Drugs: 276 (4 of those do not have a FIMM ID) 
# Samples: 26; Values: 7176;  NAs: 214
# ============================== #

dataset.section.REPORT <- newSection( "Data Set Exploration" )
# Convert the data to 'double matrix' and: 
ovarian.matrix <- ovarian.DATA  # copy
ovarian.matrix <- data.matrix(ovarian.matrix[,-c(1,2)])  # del 1st & 2nd rows 
ovarian.transponsed <- t(ovarian.matrix)  # and transform
colnames(ovarian.transponsed) <- ovarian.DATA[,2]  # assign colnames with drug names
rownames(ovarian.matrix) <- ovarian.DATA[,2]  # assign colnames with drug names

# substitute NAs with 0
nas <- is.na(ovarian.matrix); ovarian.matrix[nas] <- 0
nas <- is.na(ovarian.transponsed); ovarian.transponsed[nas] <- 0

# if DSS == 0 for each cell line: "Anastrozole" only
drop <- which(apply(ovarian.matrix,1,sum) == 0)

# We should trim those in order to run proper clustering
ovarian.matrix <- ovarian.matrix[-drop,]
ovarian.transponsed <- ovarian.transponsed[,-drop]

# SCREEN SUMMARY
drop <- fimm.DICT$FIMM.batch.ID %in% ovarian.DATA$ID.Drug
ovarian.DICT <- fimm.DICT[drop,]
ovarian.DICT <- unique(ovarian.DICT)
drop <- fimm.ANNO$CHEMBL_ID %in% ovarian.DICT$CHEMBL_ID
ovarian.ANNO <- fimm.ANNO[drop,]

ovarian.SUMMARY <- data.frame( Drugs.Count=length( unique(ovarian.DATA$Name.Drug) ), 
                                Drugs.Annotated=length( unique(ovarian.ANNO$CHEMBL_ID) ),
                                CL.Count=dim(ovarian.DATA)[2]-1, 
                                CL.Annotated=0, 
                                Anno.Count=length( ovarian.ANNO$CHEMBL_ID ),
                                ATC.L1.Count=length( unique(ovarian.ANNO$Level1) ),
                                ATC.L2.Count=length( unique(ovarian.ANNO$Level2) ),
                                ATC.L3.Count=length( unique(ovarian.ANNO$Level3) ),
                                ATC.L4.Count=length( unique(ovarian.ANNO$Level4) ),
                                ATC.L5.Count=length( unique(ovarian.ANNO$Level5) ) 
)

dataset.section.REPORT <- addTo( dataset.section.REPORT, newTable( ovarian.SUMMARY, "Data Set Summary" ) )
dataset.section.REPORT <- addTo( dataset.section.REPORT, newTable( ovarian.matrix[76:80,1:5], "Initial Data Set with DSS values" ) )

figure.file <- paste("class_distribution",".png",sep="")
png( paste( "../reports/Ovarian/", figure.file, sep="" ), width=2100, height=1100 )
#x11(); 
par(mfrow=c(2,2))
for(k in 1:4){  # "k" for levels
  count <- table(ovarian.ANNO[,k+1])
  #counts <- as.data.frame(count)
  #colnames(counts) <- c("Class", "Freq")
  #counts$Class <- reorder(counts$Class, order(counts$Freq))
  #ggplot(counts, aes(x=Class, y=Freq, fill=Class)) + geom_bar(stat="identity", position = "dodge")
  barplot(count, col=1:length(count), main=paste("Category Distribution on Level ", k, sep=""), xlab="\"Breast Data Set\"")
}
dev.off()

dataset.section.REPORT <- addTo( dataset.section.REPORT, newFigure(figure.file, "Category Distribution by ATC Levels") )

# HeatMap of the Response: modify color.map to highlight Drug Classes
# color.map <- function(mol.biol) { if (mol.biol=="ALL1/AF4") "#FF0000" else "#0000FF" }
# patientcolors <- unlist(lapply(esetSel$mol.bio, color.map))

#heatmap.2(as.matrix(ovarian.matrix, na.rm=FALSE, Colv=F, scale='none'), col=topo.colors(75), Rowv=TRUE, na.rm=TRUE,
#          key=TRUE, symkey=FALSE, trace="none", cexRow=0.5, cexCol=0.5)  # ColSideColors=patientcolors

# Full Set Density
# plot( density( ovarian.transponsed, na.rm=TRUE), main = "Full Set", xlab = "DSS" )

# Boxplot across drugs
figure.file <- paste("drugs_boxplot",".png",sep="")
png( paste( "../reports/Ovarian/", figure.file, sep="" ), width=1700, height=600 )
par(mfrow=c(1,1))
boxplot(ovarian.transponsed, col=1:length(ovarian.transponsed), xlab = "A Drug", ylab = "DSS Value Distribution", names=NA)
text(seq(1, length(rownames(ovarian.matrix)), by=1), -3, srt = 90, pos = 1, xpd = TRUE, labels=rownames(ovarian.matrix), cex=0.5)
dev.off()

dataset.section.REPORT <- addTo( dataset.section.REPORT, newFigure(figure.file, "DSS Value Distribution") )

# CLUSTERING:
# ==============
cluster.section.REPORT <- newSection( "Cluster Analysis" )
# (By "Pearson correlation as distance method")
# Alailable methods for HCLUST: "ward", "single", "complete", "average", "mcquitty", "median" or "centroid"
hr <- hclust(as.dist(1-cor(ovarian.transponsed, method="spearman")), method="complete")
hc <- hclust(as.dist(1-cor(ovarian.matrix, method="spearman")), method="complete")  

# Cuts the tree and creates color vector for clusters
mycl <- cutree(hr, h=max(hr$height)/1.5)
ovarian.ClUST <- as.data.frame(mycl)  # to use 'ovarian.ClUST' in further analysis
mycolhc <- rainbow(length(unique(mycl)), start=0.1, end=0.9); mycolhc <- mycolhc[as.vector(mycl)] 
myheatcol <- topo.colors(75)

# Creates heatmap for entire data set where the obtained clusters are indicated in the color bar
figure.file <- paste("heatmap",".png",sep="")
png( paste( "../reports/Ovarian/", figure.file, sep="" ), width=2500, height=2500 )
heatmap.2(ovarian.matrix, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc), col=myheatcol, 
          scale="row", trace="none", RowSideColors=mycolhc, cexRow=0.5, cexCol=0.5) 
dev.off()

cluster.section.REPORT <- addTo( cluster.section.REPORT, newFigure(figure.file, "Heatmap") )

# PROCESS OBTAINED CLUSTERS from "mycl"
ovarian.ClUST$DrugName <- names(mycl)
colnames(ovarian.ClUST) <- c("Cluster", "DrugName")

# Add "CHEMBL_ID" from "ovarian.DICT" -> "ovarian.ClUST"
for(alias in ovarian.ClUST$DrugName){
  if(sum(ovarian.DICT["Aliases"] == alias)){
    ovarian.ClUST[ovarian.ClUST[,"DrugName"] == alias,"CHEMBL_ID"] <- ovarian.DICT[ovarian.DICT["Aliases"] == alias,"CHEMBL_ID"][1]
  }
}

# Add "Cluster" from "ovarian.ClUST" -> "ovarian.ANNO"
not.na <- !is.na(ovarian.ClUST["CHEMBL_ID"])
for(id in ovarian.ANNO$CHEMBL_ID){
  if(sum(ovarian.ClUST[not.na,"CHEMBL_ID"] == id)){
    curr.id <- ovarian.ClUST[,"CHEMBL_ID"] == id
    curr.id <- curr.id & not.na
    ovarian.ANNO[ovarian.ANNO[,"CHEMBL_ID"] == id,"Cluster"] <- ovarian.ClUST[curr.id,"Cluster"]
  }
}

# Get Drug Lists by Clusters
ovarian.OUT <- data.frame()
for( clst in unique(ovarian.ClUST$Cluster) ){
  ovarian.OUT[clst,1] <- paste("Cluster ", clst, sep="")
  ovarian.OUT[clst,2] <- length(ovarian.ClUST[ovarian.ClUST[,"Cluster"] == clst,"DrugName"])
  ovarian.OUT[clst,3] <- gsub(",", "", toString( ovarian.ClUST[ovarian.ClUST[,"Cluster"] == clst,"DrugName"] )) 
}
colnames(ovarian.OUT) <- c("ClusterSym", "Cluster.Size", "Drug.List")

cluster.section.REPORT <- addTo( cluster.section.REPORT, newTable( ovarian.OUT, "Drug List by Clusters" ) )

figure.file <- paste("drugs_clust_boxplot",".png",sep="")
png( paste( "../reports/Ovarian/", figure.file, sep="" ), width=2500, height=2500 )
# x11(); 
par(mfrow=c(2,2))
for(k in 1:4){
  count <- table(ovarian.ANNO[,"Cluster"], ovarian.ANNO[,k+1])
  barplot(count, col=1:length(count), main=paste("Category Distribution on Level ", k, sep=""), xlab="\"Breast Data Set\"", legend=rownames(count))
}
par(mfrow=c(1,1)); dev.off() 
cluster.section.REPORT <- addTo( cluster.section.REPORT, newFigure(figure.file, "Box Plot") )


# Prints color key for cluster assignments. The numbers next to the color boxes correspond to the cluster numbers in 'mycl'
# x11(height=6, width=2); names(mycolhc) <- names(mycl); barplot(rep(10, max(mycl)), col=unique(mycolhc[hr$labels[hr$order]]), horiz=T, names=unique(mycl[hr$order])) 

# Save R Objects to a file
save(ovarian.ClUST, ovarian.ANNO, ovarian.DATA, ovarian.DICT, ovarian.SUMMARY,
     ovarian.matrix, ovarian.transponsed, ovarian.OUT, file = "leukemiaClust.RData")

REPORT <- addTo( REPORT, fimmcollection.section.REPORT, dataset.section.REPORT, cluster.section.REPORT  );
writeReport( REPORT, filename="../reports/Ovarian/REPORT" )