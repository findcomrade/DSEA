setwd("/home/comrade/Ubuntu One/DSEA/r-code")
source('pipeline_sup.R')

library(vcd)
library(grid)
library(Hmisc)
library(gplots)
library(stringr)
library(lattice)
library(reshape2)
library(RJSONIO)

### Import Data Sets and Annotations ###
# ==================================== #
read.csv(file="../datasets/merged_dss_new.csv", head=TRUE, sep="\t") -> data.MERGED
#read.csv(file="../datasets/leukemia_all_celllines_data_DSS.csv", head=TRUE, sep="\t") -> data.MERGED

# read.csv(file="~/Desktop/current/DSS2_all_SK_SH_4_LK.csv", head=TRUE, sep=",") -> data.MERGED
# read.csv(file="../datasets/merged_dss.csv", head=TRUE, sep=",") -> data.MERGED


# Exclude a Sample (to check clustering results)
tmp <- which( colnames(data.MERGED) == "SR" ) # X784_1 HUB_1613_270513_culture1_rock
data.MERGED <- data.MERGED[,-tmp]

load('RData/FimmDrugAnnotations.RData')

#########################################
###          DATA EXPLORATION       ####
#########################################

# Convert the data to 'double matrix' and: 
matrix.MERGED <- data.MERGED  # copy
matrix.MERGED <- data.matrix(matrix.MERGED[,-c(1,2)])  # del 1st & 2nd rows 
transponsed.MERGED <- t(matrix.MERGED)  # and transform
colnames(transponsed.MERGED) <- data.MERGED[,2]  # assign colnames with drug names
rownames(matrix.MERGED) <- data.MERGED[,2]  # assign colnames with drug names

# if DSS == 0 for each cell line: "Anastrozole" only
# drop <- which(apply(matrix.MERGED,1,sum) == 0)

# We should trim those in order to run proper clustering
#matrix.MERGED <- matrix.MERGED[-drop,]
#transponsed.MERGED <- transponsed.MERGED[,-drop]

# substitute NAs with 0
# nas <- is.na(matrix.MERGED); matrix.MERGED[nas] <- 0
# nas <- is.na(transponsed.MERGED); transponsed.MERGED[nas] <- 0

  
# Clean (trim) batch.ids 
data.MERGED$FIMM.Batch.ID <- data.MERGED[,"ID.Drug"]
data.MERGED[,"ID.Drug"] <- strtrim(data.MERGED[,"ID.Drug"],10)

# SCREEN SUMMARY
id <- data.MERGED[,"ID.Drug"] %in% dict.FIMM[,"FIMM.ID"]
al <- data.MERGED[,"Name.Drug"] %in% dict.FIMM[,"Drug.Name"]
drugs.MERGED <- data.frame( Drug.ID=data.MERGED[,"ID.Drug"], Drug.Name=data.MERGED[,"Name.Drug"], ID.in.Dict=id, Name.in.Dict=al )
remove(id, al)

drop <- dict.FIMM[,"FIMM.ID"] %in% data.MERGED[,"ID.Drug"]
dict.MERGED <- dict.FIMM[drop,]  # buid a dict only for leukemia drugs
dict.MERGED <- unique(dict.MERGED)
dict.MERGED$Drug.Name <- lapply(dict.MERGED[,"Drug.Name"],str_trim)
drop <- which(annotations.FIMM[,"Pubchem_CID"] %in% dict.MERGED[,"PubChem.CID"])
annotations.MERGED <- annotations.FIMM[drop,]  # buid annotation table only for leukemia drugs
remove(dict.FIMM, annotations.FIMM)

summary.MERGED <- data.frame( 
  Drugs.Count=length( unique(data.MERGED$Name.Drug) ), 
  Drugs.Annotated=length( unique(annotations.MERGED$Pubchem_CID) ),
  Anno.Count=length( annotations.MERGED$KEGG_id ),
  CL.Count=dim(data.MERGED)[2]-1, 
  CL.Annotated=0                            
)

# separate KEGG annotations:
# take only those rows drugs that are in AML project
# drop <- which( bio.class.FIMM[,"KEGG_id"] %in% annotations.MERGED[,"KEGG_id"]  )
# bio.class.AML <- bio.class.FIMM[drop,]  # Biological Classes
# plotDrugClassesDistibution(bio.class.AML, category.name='Biological')  
# 
drop <- which( target.class.FIMM[,"KEGG_id"] %in% annotations.MERGED[,"KEGG_id"]  )
target.class.AML <- target.class.FIMM[drop,]  # Targets
# plotDrugClassesDistibution(target.class.AML, category.name='Targets')  

# drop <- which( usp.class.FIMM[,"KEGG_id"] %in% annotations.MERGED[,"KEGG_id"]  )
# usp.class.AML <- usp.class.FIMM[drop,]  # Targets
# plotDrugClassesDistibution(usp.class.AML, category.name='USP')  
# 
drop <- which( neoplast.class.FIMM[,"KEGG_id"] %in% annotations.MERGED[,"KEGG_id"]  )
neoplast.class.AML <- neoplast.class.FIMM[drop,]  # Targets
# plotDrugClassesDistibution(neoplast.class.AML, category.name='Antineoplastic') 
# 
# drop <- which( atc.class.FIMM[,"KEGG_id"] %in% annotations.MERGED[,"KEGG_id"]  )
# atc.class.AML <- atc.class.FIMM[drop,]  # Targets
# plotDrugClassesDistibution(atc.class.AML, category.name='Antineoplastic') 


## Correlations & Pairs()
# qplot(x=Var1, y=Var2, data=melt(corr.sampl.MERGED), fill=value, geom="tile") + opts(axis.text.x = theme_text(angle = 90)) # <= unordered
# < cell lines >
#corr.sampl.MERGED <- cor(matrix.MERGED, use="complete.obs")
corr.sampl.MERGED <- rcorr(matrix.MERGED)
corr.sampl.MERGED <- corr.sampl.MERGED$r
#corr.sampl.MERGED <- corr.sampl.MERGED[nrow(corr.sampl.MERGED):1,]
ccl <- hclust(as.dist(1-corr.sampl.MERGED), method="complete")

par(mar=c(5.1, 4.1, 4.1, 2.1), oma=c(5,5,5,5))
heatmap.2(corr.sampl.MERGED, Rowv=TRUE, Colv=TRUE, trace="none", col=colorpanel(40, "white", "red"), dendrogram='row', revC=TRUE, cexCol=0.8, cexRow=0.7, key=TRUE)  

  # pairs(~THP.1+KG.1+SH.2+Kasumi.1+ME.1+SIG.M5+AP1060+MUTZ2+HL60.ATCC+MOLM.13_2+SKM.1_2+MONO.MAC.6_2+
#      AP1060_2+Kasumi6+HL60.ATCC_2+HL60.TB_2+OCI.AML2+OCI.AML3+MOLM.16+ML.2+NB.4+GDM.1+UT7+NOMO.1_2+AML.193_2 ,data=data.MERGED)  # +UT7+NOMO.1_2+AML.193_2+MOLT.4_2

# pairs(~CCRF.CEM+K562+RPMI8226+SR+MOLT.4_2 ,data=data.MERGED)

# < dugs >
#corr.drugs.MERGED <- cor(transponsed.MERGED, use="complete.obs")
corr.drugs.MERGED <- rcorr(transponsed.MERGED)
corr.drugs.MERGED <- corr.drugs.MERGED$r
nas <- is.na(corr.drugs.MERGED); corr.drugs.MERGED[nas] <- 0
ccl <- hclust(as.dist(1-corr.drugs.MERGED), method="complete")
heatmap.2(corr.drugs.MERGED, Rowv=as.dendrogram(ccl), Colv=as.dendrogram(ccl), trace="none", col=colorpanel(40, "blue", "white", "red"), dendrogram='none', revC=TRUE, cexCol=0.25, cexRow=0.25)


# HeatMap of the Response: modify color.map to highlight Drug Classes
color.map <- function(mol.biol) { if (mol.biol=="ALL1/AF4") "#FF0000" else "#0000FF" }
#patientcolors <- unlist(lapply(esetSel$mol.bio, color.map))
heatmap.2(as.matrix(matrix.MERGED, na.rm=FALSE, Colv=F, scale='none'), col=topo.colors(75), Rowv=TRUE, na.rm=TRUE,
          key=TRUE, symkey=FALSE, trace="none", cexRow=0.2, cexCol=0.5)  # ColSideColors=patientcolors


# CLUSTERING:
# ==============
# Alailable methods for HCLUST: "ward", "single", "complete", "average", "mcquitty", "median" or "centroid"
hr <- hclust(as.dist(1-corr.drugs.MERGED), method="complete")
hc <- hclust(as.dist(1-corr.sampl.MERGED), method="complete")  

# Cuts the tree and creates color vector for clusters
mycl <- cutree(hr, h=max(hr$height)/1.5) 
tree.DRUGS <- as.data.frame(mycl)  # to use 'tree.DRUGS' in further analysis
mysl <- cutree(hc, h=max(hr$height)/4) 
tree.SAMPLES <- as.data.frame( mysl )  # to use 'tree.SAMPLES' in further analysis

mycolhc <- rainbow(length(unique(mycl)), start=0.1, end=0.9); mycolhc <- mycolhc[as.vector(mycl)] 
mycolhs <- rainbow(length(unique(mysl)), start=0.1, end=0.9); mycolhs <- mycolhs[as.vector(mysl)] 

myheatcol <- topo.colors(75)

# Creates heatmap for entire data set where the obtained clusters are indicated in the color bar
heatmap.2(matrix.MERGED, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc), col=myheatcol, 
          scale="row", trace="none", RowSideColors=mycolhc, ColSideColors=mycolhs, cexRow=0.3, cexCol=0.5) 

# Generates JSON file to plot a circular dendogram 
halfway <- hclustToTree(hr, mycl)
jsonTree <- toJSON(halfway)
sink('/home/comrade/Projects/d3.v3/tree.json')
cat( substr(jsonTree, 4, nchar(jsonTree)-3) )
sink()


# PROCESS OBTAINED CLUSTERS from "mycl"
tree.DRUGS$DrugName <- rownames(tree.DRUGS)
colnames(tree.DRUGS) <- c("Cluster", "DrugName")
tree.DRUGS$DrugName <- lapply(tree.DRUGS[,"DrugName"], str_trim)

# The same for samples
tree.SAMPLES$SampleName <- rownames(tree.SAMPLES)
colnames(tree.SAMPLES) <- c("Cluster", "SampleName")
tree.SAMPLES$SampleName <- lapply(tree.SAMPLES[,"SampleName"], str_trim)
  
# Add "CHEMBL_ID" from "dict.MERGED" -> "tree.DRUGS"
for(alias in tree.DRUGS$DrugName){
  if(sum(dict.MERGED["Drug.Name"] == alias)){
    tree.DRUGS[tree.DRUGS[,"DrugName"] == alias,"ChEMBL.ID"] <- dict.MERGED[dict.MERGED["Drug.Name"] == alias,"ChEMBL.ID"][1]
    tree.DRUGS[tree.DRUGS[,"DrugName"] == alias,"PubChem.CID"] <- dict.MERGED[dict.MERGED["Drug.Name"] == alias,"PubChem.CID"][1]
    # tree.DRUGS[tree.DRUGS[,"DrugName"] == alias,"PubChem.CID"] <- dict.MERGED[dict.MERGED["Drug.Name"] == alias,"PubChem.CID"][1]
    tree.DRUGS[tree.DRUGS[,"DrugName"] == alias,"Class"] <- dict.MERGED[dict.MERGED["Drug.Name"] == alias,"Target"][1]
  }
}

# Set 'undefined' for NAs Classes
drop <- which(is.na(tree.DRUGS$Class))
tree.DRUGS[drop,"Class"] <- "undefined"

# Add "Cluster" from "tree.DRUGS" -> "annotations.MERGED"
not.na <- !is.na(tree.DRUGS["PubChem.CID"])
for(id in annotations.MERGED$Pubchem_CID){
  if(sum(tree.DRUGS[not.na,"PubChem.CID"] == id)){
    curr.id <- tree.DRUGS[,"PubChem.CID"] == id
    curr.id <- curr.id & not.na
    annotations.MERGED[annotations.MERGED[,"Pubchem_CID"] == id,"Cluster"] <- tree.DRUGS[curr.id,"Cluster"]
  }
}

# Get Drug Lists by Clusters
leukemia.OUT <- data.frame()
for( clst in unique(tree.DRUGS$Cluster) ){
  leukemia.OUT[clst,1] <- paste("Cluster ", clst, sep="")
  leukemia.OUT[clst,2] <- length(tree.DRUGS[tree.DRUGS[,"Cluster"] == clst,"DrugName"])
  leukemia.OUT[clst,3] <- gsub(",", "", toString( tree.DRUGS[tree.DRUGS[,"Cluster"] == clst,"DrugName"] )) 
}
colnames(leukemia.OUT) <- c("ClusterSym", "Cluster.Size", "Drug.List")

# add cluster info to AML annotations
drop <- which( target.class.AML$KEGG_id %in% annotations.MERGED$KEGG_id  )
target.class.AML[drop,"Cluster"] <- annotations.MERGED[drop, "Cluster"]


# layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE))
# for(k in 1:3){  # "k" for levels
#   count <- table(target.class.AML[,"Cluster"], target.class.AML[,k+1])    
#     
#   barplot(count, axisnames =TRUE, cex.axis = 0.8, cex.names = (0.99-k/10),
#           col=1:length(ordered.data$Freq), main=paste("Target", " : Level ", k, sep=""), las=2)
#   
#   #counts <- as.data.frame(count); #colnames(counts) <- c("Class", "Freq"); #counts$Class <- reorder(counts$Class, order(counts$Freq))
#   #ggplot(counts, aes(x=Class, y=Freq, fill=Class)) + geom_bar(stat="identity", position = "dodge")    
# }

# par(mfrow=c(2,2));
# for(k in 1:4){
#   count <- table(annotations.MERGED[,"Cluster"], annotations.MERGED[,k+1])
#   barplot(count, col=1:length(count), main=paste("Category Distribution on Level ", k, sep=""), xlab="\"Breast Data Set\"", legend=rownames(count))
# }
# par(mfrow=c(1,1)); dev.off() 


# Prints color key for cluster assignments. The numbers next to the color boxes correspond to the cluster numbers in 'mycl'
# x11(height=6, width=2); names(mycolhc) <- names(mycl); barplot(rep(10, max(mycl)), col=unique(mycolhc[hr$labels[hr$order]]), horiz=T, names=unique(mycl[hr$order])) 

# Save R Objects to a file
save(tree.DRUGS, tree.SAMPLES, annotations.MERGED, data.MERGED, dict.MERGED, summary.MERGED,
     matrix.MERGED, transponsed.MERGED, leukemia.OUT, target.class.AML, neoplast.class.AML, file = "RData/leukemiaClust.RData")
save(tree.DRUGS, tree.SAMPLES, data.MERGED, matrix.MERGED, transponsed.MERGED, leukemia.OUT, file = "RData/leukemiaClust.RData")

dev.off()
