drugSensitivity <- function(dss.matrix, cell.line, cut=25){
  # Looks for a col named as "cell.line", 
  # orders values and  generates a barplot
  
  ordered.data <- dss.matrix[ order(dss.matrix[,cell.line], decreasing = TRUE), cell.line ]  
  top.drugs <- ordered.data[ 1 : dim( topSensitive(dss.matrix,cell.line, cut) )[1] ] 
    
  col = c("mistyrose", "lightblue")
  myclorz <- ifelse(ordered.data %in% top.drugs, col[1], col[2])
  barplot( ordered.data, axisnames =TRUE, names.arg=names(ordered.data), horiz=FALSE, col= myclorz, las=2, cex.axis = 0.8, cex.names = 0.3,
                space=1, axis.lty=1, lwd=2, xaxs="i", yaxs="i", ylab="DSS", main =paste("Drug Sensitivity for \"", cell.line, "\" Cell Line", sep="") )
}

plotDrugClassesDistibution <- function(drug.classes, category.name=':'){
  
  # par(mfrow=c(2,2))
  layout(matrix(c(1,1,2,2), 2, 2, byrow = TRUE))
  for(k in 1:2){  # "k" for levels
    count <- table(drug.classes[,k+1])    
    ordered.data <- data.frame(count)
    ordered.data <- ordered.data[ order( ordered.data$Freq, decreasing = TRUE ), ]
    
    # get rid of zero-freq values
    # i have no idea why all the categories appear in the table!!!
    ordered.data <- ordered.data[ ordered.data$Freq > 0 ,]
    
    barplot(ordered.data$Freq, axisnames =TRUE, names.arg=ordered.data$Var1, cex.axis = 0.8, cex.names = 0.75,
            col=1:length(ordered.data$Freq), main=paste(category.name, " : Level ", k, sep=""), las=2)
    
    #counts <- as.data.frame(count); #colnames(counts) <- c("Class", "Freq"); #counts$Class <- reorder(counts$Class, order(counts$Freq))
    #ggplot(counts, aes(x=Class, y=Freq, fill=Class)) + geom_bar(stat="identity", position = "dodge")    
  }
  
}

topSensitive <- function(dss.matrix, cell.line, dss.threshold=25){
  # Returns a list of top sensitive drugs for a given cell line
  
  dss.values <- as.data.frame( dss.matrix[,cell.line] )
  rownames(dss.values) <- rownames(dss.matrix) 
  dss.values$DrugName <- rownames(dss.matrix) 
  colnames(dss.values) <- c(cell.line, "DrugName")
  
  # trim items with dss < threshold
  drugs.top <- as.data.frame(dss.values[dss.values[,cell.line] > dss.threshold,])
  drugs.top <- drugs.top[order(drugs.top[,cell.line], decreasing=TRUE),]  # order drugs 
  
  return (drugs.top)
}

topResistant <- function(dss.matrix, cell.line){
  # Returns a list of resistant drugs for a given cell line
  
  dss.values <- as.data.frame( dss.matrix[,cell.line] )
  rownames(dss.values) <- rownames(dss.matrix) 
  dss.values$DrugName <- rownames(dss.matrix) 
  colnames(dss.values) <- c(cell.line, "DrugName")
  
  # trim items with dss < threshold
  drugs.resist <- as.data.frame(dss.values[dss.values[,cell.line] == 0,])
  
  return (drugs.resist)
}

buildEnrichment <- function(cluster.data, drugs.sensitive, drugs.resistant){
  enrichment.table <- data.frame(Cluster=character(), Cluster.Size=numeric(), DataSet.Size=numeric(), 
                                 Sensit.Count=numeric(), Sensit.Total=numeric(), p.Value=numeric(),
                                 Resist.Count=numeric(), Resist.Total=numeric(), p.Value=numeric(), stringsAsFactors=FALSE)
  
  for(cluster in unique(cluster.data$Cluster)){
    cluster.set <- cluster.data[ cluster.data[,"Cluster"] == cluster, "DrugName" ]
    
    enrichment.table[cluster,"Cluster"] <- paste("Cluster ", cluster, sep="")
    enrichment.table[cluster,"Cluster.Size"] <- length(cluster.set)
    enrichment.table[cluster,"DataSet.Size"] <- length( unique(cluster.data$DrugName) ) 
    
    enrichment.table[cluster,"Sensit.Count"] <- sum(drugs.sensitive$DrugName %in% cluster.set)
    enrichment.table[cluster,"Sensit.Total"] <- length(drugs.sensitive$DrugName)
    
    enrichment.table[cluster,"Resist.Count"] <- sum(drugs.resistant$DrugName %in% cluster.set)
    enrichment.table[cluster,"Resist.Total"] <- length(drugs.resistant$DrugName)
  }
  
  return (enrichment.table)
}

dropJSON <- function(clusters.data, path="circular.json"){
  # Generates JSON file fot D3 library
  # 'clusters.data' - a data frame like e.g. 'leukemia.ClUST'
  
  sink(path)
  cat("{", "\n\t\"name\": \"Data Set\",", "\n\t\"children\": [\n")  # push root node
  
  for (cluster in unique(clusters.data$Cluster)){  # Over Clusters
    cat("\t{", "\n\t\"name\": ", paste("\"Cluster ", cluster, "\",", sep=""), "\n\t\"children\": [\n")
    
    last.index <- length(clusters.data[clusters.data[,"Cluster"] == cluster,"DrugName"])
    
    for(drug in clusters.data[clusters.data[,"Cluster"] == cluster,"DrugName"]){  # Over Drugs
      # extract drug respose code:
      # '1' - "sensitive"; '-1' - "resistant (DSS == 0)"; '0' - "intermediate".
      is.top <- clusters.data[clusters.data[,"DrugName"] == drug,"isTop"]  
      if( drug == clusters.data[clusters.data[,"Cluster"] == cluster,"DrugName"][last.index] ){
        cat("\t\t{", "\t\"name\": ", paste("\"",drug, "\", ",  "\"top\": ", is.top, sep=""), " } \n")    
      }
      else{ cat("\t\t{", "\t\"name\": ", paste("\"",drug, "\", ", "\"top\": ", is.top, sep=""), " }, \n") }
    }
    
    if( cluster == max(unique(clusters.data$Cluster)) ){ cat("\t]", "\n\t}") }
    else { cat("\t]", "\n\t},") }
  }
  cat("]", "\n}")
  sink()
}

cleanUpFactors <- function(my.frame){
  # This function iterates on frame's columns
  # If col is a factor we reduce the number of levels for the factor
  # to only those that currently appear in a list (in the col).
  # Also all NAs are substituted with 'undefined'
  
  for(col in 1:ncol(my.frame)){
    if(is.factor(my.frame[,col])){
      # clean levels
      my.frame[,col] <- factor(my.frame[,col])
      # add a new level
      my.frame[,col] <- factor(my.frame[,col], levels = c(levels(my.frame[,col]), "undefined"))
      my.frame[is.na( my.frame[,col] ),col] <- "undefined"      
    }
  }
  
  return(my.frame)
}

fulfillClassFrame <- function(class.frame, ids.list){
  # Takes a frame with drug classes assosiated with 
  # Kegg ids and accomplish the frame with missing ids
  # from annotations.FIMM. 
  len <- length(ids.list)
  new.frame <- data.frame( rep(NA,len), rep(NA,len), rep(NA,len),
                             rep(NA,len), rep(NA,len), ids.list   )
  colnames(new.frame) <- colnames(class.frame)
  new.frame <- rbind(class.frame, new.frame)  
  
  return(new.frame)
}