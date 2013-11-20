setwd("/home/comrade/Ubuntu One/DSEA/r-code")
source('pipeline_sup.R')

library(vcd)
library(grid)
library(Hmisc)
library(gplots)
library(stringr)
library(lattice)
library(ggplot2)
library(reshape2)
library(RJSONIO)

### Import Data Sets and Annotations ###
# ==================================== #
read.csv(file="../datasets/all_leukemia_cl_june_13_disha_astrid.csv", head=TRUE, sep=",") -> data.AML

# Exclude a Sample (to check clustering results)
tmp <- which( colnames(data.AML) == "THP.1" )
data.AML <- data.AML[,-tmp]

load('RData/FimmDrugAnnotations.RData')

#########################################
###          DATA EXPLORATION       ####
#########################################

# Convert the data to 'double matrix' and: 
matrix.AML <- data.AML  # copy
matrix.AML <- data.matrix(matrix.AML[,-c(1,2)])  # del 1st & 2nd rows 
transponsed.AML <- t(matrix.AML)  # and transform
colnames(transponsed.AML) <- data.AML[,2]  # assign colnames with drug names
rownames(matrix.AML) <- data.AML[,2]  # assign colnames with drug names

# if DSS == 0 for each cell line: "Anastrozole" only
drop <- which(apply(matrix.AML,1,sum) == 0)

# We should trim those in order to run proper clustering
matrix.AML <- matrix.AML[-drop,]
transponsed.AML <- transponsed.AML[,-drop]

# substitute NAs with 0
nas <- is.na(matrix.AML); matrix.AML[nas] <- 0
nas <- is.na(transponsed.AML); transponsed.AML[nas] <- 0


# Clean (trim) batch.ids 
data.AML$FIMM.Batch.ID <- data.AML[,"ID.Drug"]
data.AML[,"ID.Drug"] <- strtrim(data.AML[,"ID.Drug"],10)

# SCREEN SUMMARY
id <- data.AML[,"ID.Drug"] %in% dict.FIMM[,"FIMM.ID"]
al <- data.AML[,"Name.Drug"] %in% dict.FIMM[,"Drug.Name"]
drugs.AML <- data.frame( Drug.ID=data.AML[,"ID.Drug"], Drug.Name=data.AML[,"Name.Drug"], ID.in.Dict=id, Name.in.Dict=al )
remove(id, al)

drop <- dict.FIMM[,"FIMM.ID"] %in% data.AML[,"ID.Drug"]
dict.AML <- dict.FIMM[drop,]  # buid a dict only for leukemia drugs
dict.AML <- unique(dict.AML)
dict.AML$Drug.Name <- lapply(dict.AML[,"Drug.Name"],str_trim)
drop <- which(annotations.FIMM[,"Pubchem_CID"] %in% dict.AML[,"PubChem.CID"])
annotations.AML <- annotations.FIMM[drop,]  # buid annotation table only for leukemia drugs
remove(dict.FIMM, annotations.FIMM)

summary.AML <- data.frame( 
            Drugs.Count=length( unique(data.AML$Name.Drug) ), 
            Drugs.Annotated=length( unique(annotations.AML$Pubchem_CID) ),
            Anno.Count=length( annotations.AML$KEGG_id ),
            CL.Count=dim(data.AML)[2]-1, 
            CL.Annotated=0                              
)

# separate KEGG annotations:
# take only those rows drugs that are in AML project
drop <- which( bio.class.FIMM[,"KEGG_id"] %in% annotations.AML[,"KEGG_id"]  )
bio.class.AML <- bio.class.FIMM[drop,]  # Biological Classes
plotDrugClassesDistibution(bio.class.AML, category.name='Biological')  

drop <- which( target.class.FIMM[,"KEGG_id"] %in% annotations.AML[,"KEGG_id"]  )
target.class.AML <- target.class.FIMM[drop,]  # Targets
plotDrugClassesDistibution(target.class.AML, category.name='Targets')  

drop <- which( usp.class.FIMM[,"KEGG_id"] %in% annotations.AML[,"KEGG_id"]  )
usp.class.AML <- usp.class.FIMM[drop,]  # Targets
plotDrugClassesDistibution(usp.class.AML, category.name='USP')  

drop <- which( neoplast.class.FIMM[,"KEGG_id"] %in% annotations.AML[,"KEGG_id"]  )
neoplast.class.AML <- neoplast.class.FIMM[drop,]  # Targets
plotDrugClassesDistibution(neoplast.class.AML, category.name='Antineoplastic') 

drop <- which( atc.class.FIMM[,"KEGG_id"] %in% annotations.AML[,"KEGG_id"]  )
atc.class.AML <- atc.class.FIMM[drop,]  # Targets
plotDrugClassesDistibution(atc.class.AML, category.name='Antineoplastic') 


## Correlations & Pairs()
# qplot(x=Var1, y=Var2, data=melt(corr.sampl.AML), fill=value, geom="tile") + opts(axis.text.x = theme_text(angle = 90)) # <= unordered
# < cell lines >
corr.sampl.AML <- cor(matrix.AML, use="complete.obs")
#corr.sampl.AML <- corr.sampl.AML[nrow(corr.sampl.AML):1,]
ccl <- hclust(as.dist(1-corr.sampl.AML), method="complete")
heatmap.2(corr.sampl.AML, Rowv=TRUE, Colv=TRUE, trace="none", col=colorpanel(40, "blue", "white", "red"), dendrogram='none', revC=TRUE)  

# pairs(~THP.1+KG.1+SH.2+Kasumi.1+ME.1+SIG.M5+AP1060+MUTZ2+HL60.ATCC+MOLM.13_2+SKM.1_2+MONO.MAC.6_2+
  #      AP1060_2+Kasumi6+HL60.ATCC_2+HL60.TB_2+OCI.AML2+OCI.AML3+MOLM.16+ML.2+NB.4+GDM.1+UT7+NOMO.1_2+AML.193_2 ,data=data.AML)  # +UT7+NOMO.1_2+AML.193_2+MOLT.4_2

# pairs(~CCRF.CEM+K562+RPMI8226+SR+MOLT.4_2 ,data=data.AML)

# < dugs >
corr.drugs.AML <- cor(transponsed.AML, use="complete.obs")
ccl <- hclust(as.dist(1-corr.drugs.AML), method="complete")
heatmap.2(corr.drugs.AML, Rowv=as.dendrogram(ccl), Colv=as.dendrogram(ccl), trace="none", col=colorpanel(40, "blue", "white", "red"), dendrogram='none', revC=TRUE, cexCol=0.5, cexRow=0.25)

# HeatMap of the Response: modify color.map to highlight Drug Classes
color.map <- function(mol.biol) { if (mol.biol=="ALL1/AF4") "#FF0000" else "#0000FF" }
#patientcolors <- unlist(lapply(esetSel$mol.bio, color.map))
heatmap.2(as.matrix(matrix.AML, na.rm=FALSE, Colv=F, scale='none'), col=topo.colors(75), Rowv=TRUE, na.rm=TRUE,
          key=TRUE, symkey=FALSE, trace="none", cexRow=0.2, cexCol=0.5)  # ColSideColors=patientcolors


# CLUSTERING:
# ==============
# Alailable methods for HCLUST: "ward", "single", "complete", "average", "mcquitty", "median" or "centroid"
hr <- hclust(as.dist(1-corr.drugs.AML), method="complete")
hc <- hclust(as.dist(1-corr.sampl.AML), method="complete")  

# Cuts the tree and creates color vector for clusters
mycl <- cutree(hr, h=max(hr$height)/1.5); 
leukemia.ClUST <- as.data.frame(mycl)  # to use 'leukemia.ClUST' in further analysis
mycolhc <- rainbow(length(unique(mycl)), start=0.1, end=0.9); mycolhc <- mycolhc[as.vector(mycl)] 
myheatcol <- topo.colors(75)

# Creates heatmap for entire data set where the obtained clusters are indicated in the color bar
heatmap.2(matrix.AML, Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc), col=myheatcol, 
          scale="row", trace="none", RowSideColors=mycolhc, cexRow=0.3, cexCol=0.5) 


halfway <- hclustToTree(hr, mycl)
jsonTree <- toJSON(halfway)
sink('/home/comrade/Projects/d3.v3/tree.json')
cat( substr(jsonTree, 4, nchar(jsonTree)-3) )
sink()


# PROCESS OBTAINED CLUSTERS from "mycl"
leukemia.ClUST$DrugName <- rownames(leukemia.ClUST)
colnames(leukemia.ClUST) <- c("Cluster", "DrugName")
leukemia.ClUST$DrugName <- lapply(leukemia.ClUST[,"DrugName"],str_trim)

# Add "CHEMBL_ID" from "dict.AML" -> "leukemia.ClUST"
for(alias in leukemia.ClUST$DrugName){
  if(sum(dict.AML["Drug.Name"] == alias)){
    leukemia.ClUST[leukemia.ClUST[,"DrugName"] == alias,"ChEMBL.ID"] <- dict.AML[dict.AML["Drug.Name"] == alias,"ChEMBL.ID"][1]
    leukemia.ClUST[leukemia.ClUST[,"DrugName"] == alias,"PubChem.CID"] <- dict.AML[dict.AML["Drug.Name"] == alias,"PubChem.CID"][1]
    leukemia.ClUST[leukemia.ClUST[,"DrugName"] == alias,"Class"] <- dict.AML[dict.AML["Drug.Name"] == alias,"Target"][1]
  }
}

# Set 'undefined' for NAs Classes
drop <- which(is.na(leukemia.ClUST$Class))
leukemia.ClUST[drop,"Class"] <- "undefined"

# Add "Cluster" from "leukemia.ClUST" -> "annotations.AML"
not.na <- !is.na(leukemia.ClUST["PubChem.CID"])
for(id in annotations.AML$Pubchem_CID){
  if(sum(leukemia.ClUST[not.na,"PubChem.CID"] == id)){
    curr.id <- leukemia.ClUST[,"PubChem.CID"] == id
    curr.id <- curr.id & not.na
    annotations.AML[annotations.AML[,"Pubchem_CID"] == id,"Cluster"] <- leukemia.ClUST[curr.id,"Cluster"]
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

# add cluster info to AML annotations
drop <- which( target.class.AML$KEGG_id %in% annotations.AML$KEGG_id  )
target.class.AML[drop,"Cluster"] <- annotations.AML[drop, "Cluster"]


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
#   count <- table(annotations.AML[,"Cluster"], annotations.AML[,k+1])
#   barplot(count, col=1:length(count), main=paste("Category Distribution on Level ", k, sep=""), xlab="\"Breast Data Set\"", legend=rownames(count))
# }
# par(mfrow=c(1,1)); dev.off() 


# Prints color key for cluster assignments. The numbers next to the color boxes correspond to the cluster numbers in 'mycl'
# x11(height=6, width=2); names(mycolhc) <- names(mycl); barplot(rep(10, max(mycl)), col=unique(mycolhc[hr$labels[hr$order]]), horiz=T, names=unique(mycl[hr$order])) 

# Save R Objects to a file
save(leukemia.ClUST, annotations.AML, data.AML, dict.AML, summary.AML,
     matrix.AML, transponsed.AML, leukemia.OUT, file = "RData/leukemiaClust.RData")

