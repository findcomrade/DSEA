setwd("/home/comrade/Ubuntu One/DSEA/r-code")
load("../datasets/merged_chembl.RData")

library(vcd)
library(grid)
library(gplots)
library(lattice)
library(ggplot2)
library(reshape2)
require(Nozzle.R1)

REPORT <- newCustomReport( "DSEA:", asEmph(" Leukemia Project.") )
### Import Data Sets and Annotations ###
# ==================================== #
merged_chembl_ATC_SCID -> chembl.ANNO  # Chembl Drug Annotations
read.csv(file="../datasets/all_leukemia_cl_june_13_disha_astrid.csv", head=TRUE, sep=",") -> leukemia.DATA
# exclude a sample:
tmp <- which( colnames(leukemia.DATA) == "RPMI8226" )
leukemia.DATA <- leukemia.DATA[,-tmp]

read.csv(file="../datasets/fimm_chembl_dict.csv", head=TRUE, sep=",") -> fimm.DICT
remove(merged_chembl_ATC_SCID, tmp)  # cleanup

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

### Leukemia DSS values:       ###
# Drugs: 240; Samples: 32
# Values: 7680;  NAs: 74
# ============================== #
dataset.section.REPORT <- newSection( "Data Set Exploration" )
# Convert the data to 'double matrix' and: 
leukemia.matrix <- leukemia.DATA  # copy
leukemia.matrix <- data.matrix(leukemia.matrix[,-c(1,2)])  # del 1st & 2nd rows 
leukemia.transponsed <- t(leukemia.matrix)  # and transform
colnames(leukemia.transponsed) <- leukemia.DATA[,2]  # assign colnames with drug names
rownames(leukemia.matrix) <- leukemia.DATA[,2]  # assign colnames with drug names

# if DSS == 0 for each cell line: "Anastrozole" only
drop <- which(apply(leukemia.matrix,1,sum) == 0)

# We should trim those in order to run proper clustering
leukemia.matrix <- leukemia.matrix[-drop,]
leukemia.transponsed <- leukemia.transponsed[,-drop]

# substitute NAs with 0
nas <- is.na(leukemia.matrix); leukemia.matrix[nas] <- 0
nas <- is.na(leukemia.transponsed); leukemia.transponsed[nas] <- 0

# SCREEN SUMMARY
drop <- fimm.DICT$FIMM.batch.ID %in% leukemia.DATA$ID.Drug
leukemia.DICT <- fimm.DICT[drop,]
leukemia.DICT <- unique(leukemia.DICT)
drop <- fimm.ANNO$CHEMBL_ID %in% leukemia.DICT$CHEMBL_ID
leukemia.ANNO <- fimm.ANNO[drop,]

leukemia.SUMMARY <- data.frame( Drugs.Count=length( unique(leukemia.DATA$Name.Drug) ), 
                              Drugs.Annotated=length( unique(leukemia.ANNO$CHEMBL_ID) ),
                              CL.Count=dim(leukemia.DATA)[2]-1, 
                              CL.Annotated=0, 
                              Anno.Count=length( leukemia.ANNO$CHEMBL_ID ),
                              ATC.L1.Count=length( unique(leukemia.ANNO$Level1) ),
                              ATC.L2.Count=length( unique(leukemia.ANNO$Level2) ),
                              ATC.L3.Count=length( unique(leukemia.ANNO$Level3) ),
                              ATC.L4.Count=length( unique(leukemia.ANNO$Level4) ),
                              ATC.L5.Count=length( unique(leukemia.ANNO$Level5) ) 
)

dataset.section.REPORT <- addTo( dataset.section.REPORT, newTable( leukemia.SUMMARY, "Data Set Summary" ) )
dataset.section.REPORT <- addTo( dataset.section.REPORT, newTable( leukemia.matrix[76:80,1:5], "Initial Data Set with DSS values" ) )

figure.file <- paste("class_distribution",".png",sep="")
png( paste( "../reports/Leukemia/", figure.file, sep="" ), width=2100, height=1100 )
#x11(); 
par(mfrow=c(2,2))
for(k in 1:4){  # "k" for levels
  count <- table(leukemia.ANNO[,k+1])
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

#heatmap.2(as.matrix(leukemia.matrix, na.rm=FALSE, Colv=F, scale='none'), col=topo.colors(75), Rowv=TRUE, na.rm=TRUE,
#          key=TRUE, symkey=FALSE, trace="none", cexRow=0.5, cexCol=0.5)  # ColSideColors=patientcolors

# Full Set Density
# plot( density( leukemia.transponsed, na.rm=TRUE), main = "Full Set", xlab = "DSS" )

# Boxplot across drugs
figure.file <- paste("drugs_boxplot",".png",sep="")
png( paste( "../reports/Leukemia/", figure.file, sep="" ), width=1700, height=600 )
par(mfrow=c(1,1))
boxplot(leukemia.transponsed, col=1:length(leukemia.transponsed), xlab = "A Drug", ylab = "DSS Value Distribution", names=NA)
text(seq(1, length(rownames(leukemia.matrix)), by=1), -3, srt = 90, pos = 1, xpd = TRUE, labels=rownames(leukemia.matrix), cex=0.5)
dev.off()

dataset.section.REPORT <- addTo( dataset.section.REPORT, newFigure(figure.file, "DSS Value Distribution") )

# CLUSTERING:
# ==============
cluster.section.REPORT <- newSection( "Cluster Analysis" )
# (By "Pearson correlation as distance method")
# Alailable methods for HCLUST: "ward", "single", "complete", "average", "mcquitty", "median" or "centroid"
hr <- hclust(as.dist(1-cor(leukemia.transponsed, method="pearson")), method="complete")
hc <- hclust(as.dist(1-cor(leukemia.matrix, method="spearman")), method="complete")  

# Cuts the tree and creates color vector for clusters
mycl <- cutree(hr, h=max(hr$height)/1.5); 
leukemia.ClUST <- as.data.frame(mycl)  # to use 'leukemia.ClUST' in further analysis
mycolhc <- rainbow(length(unique(mycl)), start=0.1, end=0.9); mycolhc <- mycolhc[as.vector(mycl)] 
myheatcol <- topo.colors(75)

# Creates heatmap for entire data set where the obtained clusters are indicated in the color bar
figure.file <- paste("heatmap",".png",sep="")
png( paste( "../reports/Leukemia/", figure.file, sep="" ), width=2500, height=2500 )
heatmap.2(leukemia.matrix, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc), col=myheatcol, 
          scale="row", trace="none", RowSideColors=mycolhc, cexRow=0.5, cexCol=0.5) 
dev.off()

cluster.section.REPORT <- addTo( cluster.section.REPORT, newFigure(figure.file, "Heatmap") )

# PROCESS OBTAINED CLUSTERS from "mycl"
leukemia.ClUST$DrugName <- rownames(leukemia.ClUST)
colnames(leukemia.ClUST) <- c("Cluster", "DrugName")

# Add "CHEMBL_ID" from "leukemia.DICT" -> "leukemia.ClUST"
for(alias in leukemia.ClUST$DrugName){
  if(sum(leukemia.DICT["Aliases"] == alias)){
    leukemia.ClUST[leukemia.ClUST[,"DrugName"] == alias,"CHEMBL_ID"] <- leukemia.DICT[leukemia.DICT["Aliases"] == alias,"CHEMBL_ID"][1]
  }
}

# Add "Cluster" from "leukemia.ClUST" -> "breat.ANNO"
not.na <- !is.na(leukemia.ClUST["CHEMBL_ID"])
for(id in leukemia.ANNO$CHEMBL_ID){
  if(sum(leukemia.ClUST[not.na,"CHEMBL_ID"] == id)){
    curr.id <- leukemia.ClUST[,"CHEMBL_ID"] == id
    curr.id <- curr.id & not.na
    leukemia.ANNO[leukemia.ANNO[,"CHEMBL_ID"] == id,"Cluster"] <- leukemia.ClUST[curr.id,"Cluster"]
  }
}

# Get Drug Lists by Clusters
leukemia.OUT <- data.frame()
for( clst in unique(leukemia.ClUST$Cluster) ){
  leukemia.OUT[clst,1] <- paste("Cluster ", clst, sep="")
  leukemia.OUT[clst,2] <- length(leukemia.ClUST[leukemia.ClUST[,"Cluster"] == clst,"DrugName"])
  leukemia.OUT[clst,3] <- gsub(",", "", toString( leukemia.ClUST[leukemia.ClUST[,"Cluster"] == clst,"DrugName"] )) 
}
colnames(leukemia.OUT) <- c("ClusterSym", "Cluster.Size", "Drug.List")

cluster.section.REPORT <- addTo( cluster.section.REPORT, newTable( leukemia.OUT, "Drug List by Clusters" ) )

figure.file <- paste("drugs_clust_boxplot",".png",sep="")
png( paste( "../reports/Leukemia/", figure.file, sep="" ), width=2500, height=2500 )
# x11(); 
par(mfrow=c(2,2))
for(k in 1:4){
  count <- table(leukemia.ANNO[,"Cluster"], leukemia.ANNO[,k+1])
  barplot(count, col=1:length(count), main=paste("Category Distribution on Level ", k, sep=""), xlab="\"Breast Data Set\"", legend=rownames(count))
}
par(mfrow=c(1,1)); dev.off() 
cluster.section.REPORT <- addTo( cluster.section.REPORT, newFigure(figure.file, "Box Plot") )


# Prints color key for cluster assignments. The numbers next to the color boxes correspond to the cluster numbers in 'mycl'
# x11(height=6, width=2); names(mycolhc) <- names(mycl); barplot(rep(10, max(mycl)), col=unique(mycolhc[hr$labels[hr$order]]), horiz=T, names=unique(mycl[hr$order])) 

# Save R Objects to a file
save(leukemia.ClUST, leukemia.ANNO, leukemia.DATA, leukemia.DICT, leukemia.SUMMARY,
     leukemia.matrix, leukemia.transponsed, leukemia.OUT, file = "leukemiaClust.RData")

REPORT <- addTo( REPORT, fimmcollection.section.REPORT, dataset.section.REPORT, cluster.section.REPORT  );
writeReport( REPORT, filename="../reports/Leukemia/REPORT" )