setwd("/home/comrade/Ubuntu One/DSEA/r-code")
load("../datasets/merged_chembl.RData")

library(vcd)
library(grid)
library(gplots)
library(lattice)
library(ggplot2)
library(reshape2)
require(Nozzle.R1)

REPORT <- newCustomReport( "DSEA:", asEmph(" Breast Cancer Tamoxifen Project.") )
### Import Data Sets and Annotations ###
# ==================================== #
merged_chembl_ATC_SCID -> chembl.ANNO  # Chembl Drug Annotations
read.csv(file="../datasets/original_dss_breast_cancer_tamoxifen_project_sk_sh.csv", head=TRUE, sep=",") -> breast.DATA
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

### Breast DSS values:         ###
# Drugs: 239; Samples: 16
# Values: 3824; NAs: 0
# ============================== #
dataset.section.REPORT <- newSection( "Data Set Exploration" )
# Convert the data to 'double matrix' and: 
breast.matrix <- breast.DATA  # copy
breast.matrix <- data.matrix(breast.matrix[,-1])  # del 1st row and transform
breast.transponsed <- t(breast.matrix)
colnames(breast.transponsed) <- breast.DATA[,1]  # assign colnames with drug names
rownames(breast.matrix) <- breast.DATA[,1]  # assign colnames with drug names

# There are 58 Drugs with DSS == 0 for each cell line
drop <- which(apply(breast.matrix,1,sum) == 0)

# Not Effective Drugs
non.eff <- breast.matrix[drop,]
ind <- fimm.DICT[,"Aliases"] %in% rownames(non.eff)
ali <- fimm.DICT[ind,"CHEMBL_ID"]
ind <- fimm.ANNO[,"CHEMBL_ID"] %in% ali
non.ANNO <- fimm.ANNO[ind,]

for(k in 1:dim(non.ANNO)[1]){
  chembl.id <- non.ANNO[k,"CHEMBL_ID"]
  ind <- fimm.DICT[,"CHEMBL_ID"] == chembl.id
  non.ANNO[k,"Aliases"] <- gsub(",", "", toString(fimm.DICT[ind,"Aliases"]))
}


# We should trim those in order to run proper clustering
breast.matrix <- breast.matrix[-drop,]
breast.transponsed <- breast.transponsed[,-drop]


# SCREEN SUMMARY
drop <- fimm.DICT$Aliases %in% breast.DATA$Name.Drug  
breast.DICT <- fimm.DICT[drop,]
drop <- fimm.ANNO$CHEMBL_ID %in% breast.DICT$CHEMBL_ID
breast.ANNO <- fimm.ANNO[drop,]

breast.SUMMARY <- data.frame( Drugs.Count=length( unique(breast.DATA$Name.Drug) ), 
                              Drugs.Annotated=length( unique(breast.ANNO$CHEMBL_ID) ),
                              CL.Count=dim(breast.DATA)[2]-1, 
                              CL.Annotated=0, 
                              Anno.Count=length( breast.ANNO$CHEMBL_ID ),
                              ATC.L1.Count=length( unique(breast.ANNO$Level1) ),
                              ATC.L2.Count=length( unique(breast.ANNO$Level2) ),
                              ATC.L3.Count=length( unique(breast.ANNO$Level3) ),
                              ATC.L4.Count=length( unique(breast.ANNO$Level4) ),
                              ATC.L5.Count=length( unique(breast.ANNO$Level5) ) 
)

dataset.section.REPORT <- addTo( dataset.section.REPORT, newTable( breast.SUMMARY, "Data Set Summary" ) )
dataset.section.REPORT <- addTo( dataset.section.REPORT, newTable( breast.matrix[76:80,1:5], "Initial Data Set with DSS values" ) )

figure.file <- paste("class_distribution",".png",sep="")
png( paste( "../reports/BreastTamoxifen/", figure.file, sep="" ), width=2100, height=1100 )
#x11(); 
par(mfrow=c(2,2))
for(k in 1:4){  # "k" for levels
  count <- table(breast.ANNO[,k+1])
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

#heatmap.2(as.matrix(breast.matrix, na.rm=FALSE, Colv=F, scale='none'), col=topo.colors(75), Rowv=TRUE, na.rm=TRUE,
#          key=TRUE, symkey=FALSE, trace="none", cexRow=0.5, cexCol=0.5)  # ColSideColors=patientcolors

# Full Set Density
# plot( density( breast.transponsed, na.rm=TRUE), main = "Full Set", xlab = "DSS" )

# Boxplot across drugs
figure.file <- paste("drugs_boxplot",".png",sep="")
png( paste( "../reports/BreastTamoxifen/", figure.file, sep="" ), width=1700, height=600 )
par(mfrow=c(1,1))
boxplot(breast.transponsed, col=1:length(breast.transponsed), xlab = "A Drug", ylab = "DSS Value Distribution", names=NA)
text(seq(1, length(rownames(breast.matrix)), by=1), -3, srt = 90, pos = 1, xpd = TRUE, labels=rownames(breast.matrix), cex=0.5)
dev.off()

dataset.section.REPORT <- addTo( dataset.section.REPORT, newFigure(figure.file, "DSS Value Distribution") )

# CLUSTERING:
# ==============
cluster.section.REPORT <- newSection( "Cluster Analysis" )
# (By "Pearson correlation as distance method")
# Alailable methods for HCLUST: "ward", "single", "complete", "average", "mcquitty", "median" or "centroid"
hr <- hclust(as.dist(1-cor(breast.transponsed, method="spearman")), method="complete")
hc <- hclust(as.dist(1-cor(breast.matrix, method="spearman")), method="complete")  

# Cuts the tree and creates color vector for clusters
mycl <- cutree(hr, h=max(hr$height)/1.5); 
breast.ClUST <- as.data.frame(mycl)  # to use 'breast.ClUST' in further analysis
mycolhc <- rainbow(length(unique(mycl)), start=0.1, end=0.9); mycolhc <- mycolhc[as.vector(mycl)] 
myheatcol <- topo.colors(75)

# Creates heatmap for entire data set where the obtained clusters are indicated in the color bar
figure.file <- paste("heatmap",".png",sep="")
png( paste( "../reports/BreastTamoxifen/", figure.file, sep="" ), width=2500, height=2500 )
heatmap.2(breast.matrix, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc), col=myheatcol, 
          scale="row", trace="none", RowSideColors=mycolhc, cexRow=0.5, cexCol=0.5) 
dev.off()

cluster.section.REPORT <- addTo( cluster.section.REPORT, newFigure(figure.file, "Heatmap") )

# PROCESS OBTAINED CLUSTERS from "mycl"
breast.ClUST$DrugName <- rownames(breast.ClUST)
colnames(breast.ClUST) <- c("Cluster", "DrugName")

# Add "CHEMBL_ID" from "breast.DICT" -> "breast.ClUST"
for(alias in breast.ClUST$DrugName){
  if(sum(breast.DICT["Aliases"] == alias)){
    breast.ClUST[breast.ClUST[,"DrugName"] == alias,"CHEMBL_ID"] <- breast.DICT[breast.DICT["Aliases"] == alias,"CHEMBL_ID"][1]
  }
}

# Add "Cluster" from "breast.ClUST" -> "breat.ANNO"
not.na <- !is.na(breast.ClUST["CHEMBL_ID"])
for(id in breast.ANNO$CHEMBL_ID){
  if(sum(breast.ClUST[not.na,"CHEMBL_ID"] == id)){
    curr.id <- breast.ClUST[,"CHEMBL_ID"] == id
    curr.id <- curr.id & not.na
    breast.ANNO[breast.ANNO[,"CHEMBL_ID"] == id,"Cluster"] <- breast.ClUST[curr.id,"Cluster"]
  }
}

# Get Drug Lists by Clusters
breast.OUT <- data.frame()
for( clst in unique(breast.ClUST$Cluster) ){
  breast.OUT[clst,1] <- paste("Cluster ", clst, sep="")
  breast.OUT[clst,2] <- length(breast.ClUST[breast.ClUST[,"Cluster"] == clst,"DrugName"])
  breast.OUT[clst,3] <- gsub(",", "", toString( breast.ClUST[breast.ClUST[,"Cluster"] == clst,"DrugName"] )) 
}
colnames(breast.OUT) <- c("ClusterSym", "Cluster.Size", "Drug.List")

cluster.section.REPORT <- addTo( cluster.section.REPORT, newTable( breast.OUT, "Drug List by Clusters" ) )

figure.file <- paste("drugs_clust_boxplot",".png",sep="")
png( paste( "../reports/BreastTamoxifen/", figure.file, sep="" ), width=2500, height=2500 )
# x11(); 
par(mfrow=c(2,2))
for(k in 1:4){
  count <- table(breast.ANNO[,"Cluster"], breast.ANNO[,k+1])
  barplot(count, col=1:length(count), main=paste("Category Distribution on Level ", k, sep=""), xlab="\"Breast Data Set\"", legend=rownames(count))
}
par(mfrow=c(1,1)); dev.off() 
cluster.section.REPORT <- addTo( cluster.section.REPORT, newFigure(figure.file, "Box Plot") )


# Prints color key for cluster assignments. The numbers next to the color boxes correspond to the cluster numbers in 'mycl'
# x11(height=6, width=2); names(mycolhc) <- names(mycl); barplot(rep(10, max(mycl)), col=unique(mycolhc[hr$labels[hr$order]]), horiz=T, names=unique(mycl[hr$order])) 

REPORT <- addTo( REPORT, fimmcollection.section.REPORT, dataset.section.REPORT, cluster.section.REPORT  );
writeReport( REPORT, filename="../reports/BreastTamoxifen/REPORT" )