source('source/dsea_aux.R')
source('pipeline_sup.R')
library(grid)
library(gplots)
library(ggplot2)
library(reshape2)
library(RJSONIO)

read.csv(file="../datasets/DSS2_merged_34-samples_FO3D-16-onwards_2014-04-15.csv", head=TRUE, sep="\t") -> initial.DSRT
read.csv(file="../datasets/Samples.csv", head=TRUE, sep="\t")        -> samples.Supp

# ==================================== #
### Data Preprocessing               ###
# ==================================== #
filt <- apply(data.matrix(initial.DSRT[,-c(1,2)]), 1, function(x) sum(is.na(x)) < 0.50 * length(x))
initial.DSRT <- initial.DSRT[filt,]                   # Drop Drugs with NAs more than 50%

matrix.DSRT  <- data.matrix(initial.DSRT[,-c(1,2)])   # del 1st (ID) & 2nd (Alias) cols 
tmatrix.DSRT <- t(matrix.DSRT)                        # and transform
rownames(matrix.DSRT)  <- initial.DSRT[,2]            # assign drug names
colnames(tmatrix.DSRT) <- initial.DSRT[,2]            # assign drug names

remove(filt)

### END OF Data Preprocessing        ###
# ==================================== #


### Upload corresponding data set wit clusters
load('RData/Clusters.RData')


### Prepare classes 
# ==================================== #
read.csv(file="../datasets/drug_annotations.csv", head=TRUE, sep="\t")   -> drug.dictionary

# Remove duplicated drugs (same drugs but from different vedors or in different concentrations)
dup.ind  <- which( duplicated(drug.dictionary$Drug.name) )
drug.dictionary  <- drug.dictionary[-dup.ind,]


# Drug Classes Data Frame
drug.classes  <- data.frame( Class=unique(drug.dictionary$Class), total.freq=NA )
rownames(drug.classes)  <- drug.classes$Class
# Count abundance in the entire data set
for(cl in drug.classes$Class){
  drug.classes[cl,"total.freq"]  <- sum(drug.dictionary[,"Class"] == cl)
}

# Add Class counts for each cluster
unique.clusters  <- unique(tree.Drugs$Cluster)
for(cluster in unique.clusters){
  row.ind       <- tree.Drugs[,"Cluster"] == cluster
  clst.drugs    <- data.frame( Drug=as.character(tree.Drugs[row.ind, "DrugName"]) )
  rownames(clst.drugs)  <- clst.drugs$Drug
  clst.drugs$class <- NA
  for(sub.drug in clst.drugs$Drug){
    ind  <- grep(sub.drug, drug.dictionary[,"Drug.name"], ignore.case=TRUE, value=FALSE)
    cat  <- as.character( drug.dictionary[ind[1],"Class"] ) 
    if(is.na(cat)) cat  <- "Other"
    clst.drugs[sub.drug,"class"]  <- cat
  }
  add.col  <- paste("Cluster",cluster,"freq",sep=".")
  drug.classes[,add.col]  <- NA
  for(item in drug.classes$Class){
    count  <- sum(clst.drugs$class == item)
    drug.classes[item,add.col]  <- count
  }
}
# END OF : Drug Classes Data Frame

# Add Class counts for each cluster
drug.partition.list[,"Aver.Class"] <- NA
rownames(drug.partition.list)  <- drug.partition.list$ClusterSym
col.ind  <- grep("Cluster", colnames(drug.classes))
if (length(drug.partition.list$ClusterSym) == length(col.ind)){
  for(clst in seq(1:length(drug.partition.list$ClusterSym))){
    tbl1.col  <- col.ind[clst]
    c.ind   <- which.max(drug.classes[,tbl1.col]) 
    c.name  <- as.character(drug.classes[c.ind,"Class"])
    drug.partition.list[clst,"Aver.Class"]  <- c.name
  }  
} else {print("Error: Clusters count does not match in drug.classes and drug.partition.list.")}
# END OF drug.partition


# Add Drug Class to tree representation (new col to 'tree.DRUGS' )
tree.Drugs[,"DrugClass"] <- NA
for(drug in tree.Drugs$DrugName){
  # find a drug in drug dictionary
  ind  <- grep(drug, drug.dictionary[,"Drug.name"], ignore.case=TRUE, value=FALSE)
  cat  <- as.character( drug.dictionary[ind[1],"Class"] )
  if(is.na(cat)) cat  <- "Other"
  tree.Drugs[drug,"DrugClass"]  <- cat
}

tree.Drugs[,"ClusterAver"] <- NA
for(drug in tree.Drugs$DrugName){
  # find a drug in drug dictionary
  clst.ind  <- tree.Drugs[drug,"Cluster"]
  cat  <- drug.partition.list[clst.ind,"Aver.Class"]
  if(is.na(cat)) cat  <- "Other"
  if(cat == "Other"){
    # when cluster aver is 'Other' we obtain annotation from dictionary
    ind  <- grep(drug, drug.dictionary[,"Drug.name"], ignore.case=TRUE, value=FALSE)
    cat  <- as.character( drug.dictionary[ind[1],"Class"] )
    if(is.na(cat)) cat  <- "Other"
    tree.Drugs[drug,"ClusterAver"]  <- cat    
  } else{
    # when class comes from cluster average
    tree.Drugs[drug,"ClusterAver"]  <- cat
  }  
}

dropJSONAnnotations(clusters.data=tree.Drugs, class.col="ClusterAver", path='Results/json/drug_annotations.json')
