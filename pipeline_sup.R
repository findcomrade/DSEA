plotDrugSensitivity <- function(dss.matrix, cell.line){
  # Looks for a col named as "cell.line", 
  # orders values and  generates a barplot
  
  ordered.data <- dss.matrix[ order(dss.matrix[,cell.line], decreasing = TRUE), cell.line ]
  top.drugs <- ordered.data[ 1 : dim(topSensitive(dss.matrix,cell.line))[1] ]  # ordered.data[1:cutoff]
  
  col = c("mistyrose", "lightblue")
  myclorz <- ifelse(ordered.data %in% top.drugs, col[1], col[2])
  barplot( ordered.data, axisnames =TRUE, names.arg=names(ordered.data), horiz=FALSE, col= myclorz, las=2, cex.axis = 0.8, cex.names = 0.3,
                space=1, axis.lty=1, lwd=2, xaxs="i", yaxs="i", ylab="DSS", main =paste("Drug Sensitivity for \"", cell.line, "\" Cell Line", sep="") )
}

topSensitive <- function(dss.matrix, cell.line){
  # Returns a list of top sensitive drugs for a given cell line
  
  dss.threshold <- 25
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
  
  dss.threshold <- 25
  dss.values <- as.data.frame( dss.matrix[,cell.line] )
  rownames(dss.values) <- rownames(dss.matrix) 
  dss.values$DrugName <- rownames(dss.matrix) 
  colnames(dss.values) <- c(cell.line, "DrugName")
  
  # trim items with dss < threshold
  drugs.top <- as.data.frame(dss.values[dss.values[,cell.line] > dss.threshold,])
  drugs.top <- drugs.top[order(drugs.top[,cell.line], decreasing=TRUE),]  # order drugs 
  
  return (drugs.top)
}

dropJSON <- function(clusters.data, path="circular.json"){
  # Generates JSON file fot D3 library
  # 'clusters.data' - a data frame like e.g. 'leukemia.ClUST'
  
  sink(path)
  cat("{", "\n\t\"name\": \"Data Set\",", "\n\t\"children\": [\n")
  for (cluster in unique(clusters.data$Cluster)){
    cat("\t{", "\n\t\"name\": ", paste("\"Cluster ", cluster, "\",", sep=""), "\n\t\"children\": [\n")
    
    last.index <- length(clusters.data[clusters.data[,"Cluster"] == cluster,"DrugName"])
    for(drug in clusters.data[clusters.data[,"Cluster"] == cluster,"DrugName"]){
      if( drug == clusters.data[clusters.data[,"Cluster"] == cluster,"DrugName"][last.index] ){
        cat("\t\t{", "\t\"name\": ", paste("\"",drug,"\"", sep=""), " } \n")    
      }
      else{ cat("\t\t{", "\t\"name\": ", paste("\"",drug,"\"", sep=""), " }, \n") }
    }
    if( cluster == max(unique(clusters.data$Cluster)) ){ cat("\t]", "\n\t}") }
    else { cat("\t]", "\n\t},") }
  }
  cat("]", "\n}")
  sink()
}