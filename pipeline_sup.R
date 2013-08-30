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

pValue <- function(sensitive.target, sensitive.reference, overlap, total.count){
  # Performs Fisher's exact test for testing the null of independence of rows and columns;
  # sensitive.target <- 9;  sensitive.reference <- 20;  overlap <- 6;
  if (!overlap){
    return (0)
  }
  
  nonsense.target <- sensitive.target - overlap
  nonsense.reference <- total.count - sensitive.reference
  
  DrugTesting <-
    matrix( c( (overlap-1), nonsense.target, sensitive.reference, nonsense.reference),
            nrow = 2,
            dimnames = list( Guess = c("Sense", "NonSense"),
                             Truth = c("Target", "Ref")) )
  
  test <- fisher.test(DrugTesting, alternative = "greater")
  return ( as.numeric(test["p.value"]) )
}

# Enrichment: Drug-based path
buildEnrichmentD <- function(cluster.data, drugs.sensitive, drugs.resistant){
  enrichment.table <- data.frame(Cluster=character(), Cluster.Size=numeric(), DataSet.Size=numeric(), 
                                 Sensit.Overlap=numeric(), Sensit.Total=numeric(), pVal.S=numeric(),
                                 Resist.Overlap=numeric(), Resist.Total=numeric(), pVal.R=numeric(), stringsAsFactors=FALSE)
  
  for(cluster in unique(cluster.data$Cluster)){
    cluster.set <- cluster.data[ cluster.data[,"Cluster"] == cluster, "DrugName" ]
    
    cluster.size <- length(cluster.set)
    total.size <- length( unique(cluster.data$DrugName) ) 
    enrichment.table[cluster,"Cluster"] <- paste("Cluster ", cluster, sep="")
    enrichment.table[cluster,"Cluster.Size"] <- cluster.size
    enrichment.table[cluster,"DataSet.Size"] <- total.size
    
    sensitive.count <- length(drugs.sensitive$DrugName)
    sensitive.overlap <- sum(drugs.sensitive$DrugName %in% cluster.set)
    enrichment.table[cluster,"Sensit.Overlap"] <- sensitive.overlap
    enrichment.table[cluster,"Sensit.Total"] <- sensitive.count
    enrichment.table[cluster,"pVal.S"] <- pValue(sensitive.count, cluster.size, sensitive.overlap, total.size)
    
    resistant.count <- length(drugs.resistant$DrugName)
    resistant.overlap <- sum(drugs.resistant$DrugName %in% cluster.set)
    enrichment.table[cluster,"Resist.Overlap"] <- resistant.overlap
    enrichment.table[cluster,"Resist.Total"] <- length(drugs.resistant$DrugName)
    enrichment.table[cluster,"pVal.R"] <- pValue(resistant.count, cluster.size, resistant.overlap, total.size)
  }
  enrichment.table$pVal.Sensitive.Adj <- p.adjust(enrichment.table$pVal.S, "bonferroni" )
  enrichment.table$pVal.Resistant.Adj <- p.adjust(enrichment.table$pVal.R, "bonferroni" )
  return (enrichment.table)
}

# Enrichment: Cell-Line-based path
buildEnrichmentCL <- function(dss.db, dss.new, sample.name){
  # 'target.drugs' - top sensitive drugs (or resistant) of a sample for enrichment;
  # 'dss.db' - already processed and clustered data;
  # 'dss.new' - dss values for new samples;
  # 'sample.name' - target sample/cell line name - 
  # obviously, it have to exist in dss.new, but NOT in dss.db.
  
  drugs.in.db <- length(  unique( rownames(dss.db) )  )
  target.drugs <- topSensitive(dss.new, sample.name)
  target.sens.count <- length( unique(target.drugs$DrugName) )
  
  Cell.Line <- c();  top.count <- c();  screened.count <- c();  SCORE <- c(); p.value <- c();
  all.cell.lines <- colnames(dss.db)
  if ( !(sample.name %in% all.cell.lines) ){
    for (c.line in all.cell.lines){
      Cell.Line <- c( Cell.Line, paste(c.line,"",sep="") )
      ref.drugs <- topSensitive(dss.db, c.line)
      top.count <- c(  top.count, length( unique(ref.drugs$DrugName) )  )
      screened.count <- c(screened.count, sum( !is.na(dss.db[,c.line]) ))
      overlap <- sum( ref.drugs$DrugName %in% target.drugs$DrugName )
      SCORE <- c(SCORE, overlap)
      p.value <- c( p.value, pValue(target.sens.count, length(unique(ref.drugs$DrugName)), overlap, drugs.in.db) )
    }
  } else{
    return ('Specified sample name already exists in the DB.')
  }
  
  enrichment.table <- data.frame(Cell.Line, top.count, screened.count, SCORE, p.value)
  enrichment.table$bf.correction <- p.adjust(p.value, "bonferroni")
  enrichment.table <- enrichment.table[ enrichment.table$SCORE != 0, ]  # trim SCORE == 0
  enrichment.table <- enrichment.table[ order(enrichment.table[,6]), ]  # order by bonferroni corrected p-value
  enrichment.table$perc <- as.numeric(format(round(enrichment.table[,"SCORE"]/enrichment.table[,"top.count"] *100, 2), nsmall=2))  # calculate percentage
  enrichment.table$significant <- enrichment.table$bf.correction <= 0.05
  
  return (enrichment.table)
}

# Correlation Table
correlationTable <- function(dss.db, dss.sample, sample.name){
  # 'dss.sample' - dss values for a target sample
  # 'dss.db' - dss values in db.
  # We temporary use cbind() assuming that dss-value vectors
  # are of the same length and drugs are in the same order!
  
  merged <- cbind(dss.sample, dss.db)
  colnames(merged) <- c( sample.name, colnames(dss.db))
  corr.data <- rcorr(merged)
  rs <- corr.data$r[sample.name,]
  ps <- corr.data$P[sample.name,]
  corr.table <- data.frame(Corr=rs, pVal=ps)
  
  # drop NA row is correlation with itself
  drop <- which(is.na(corr.table$pVal))
  corr.table <- corr.table[-drop,]  
  
  # order according to correlation value
  corr.table <- corr.table[ order(corr.table$Corr, decreasing=TRUE), ]
  corr.table$Cell.Line <- as.factor(rownames(corr.table))
  return( corr.table )    
}

# Enriched Cell Lines
getEnrichedCellLines <- function(corr.table, enrich.table, r.min, p.max){
  corr.cl <- corr.table[ corr.table$Corr > r.min, c("Cell.Line", "Corr") ]
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

### JSON convertion for D3 library ###

createLeafNode <- function(hclust, i, clusters) {
  list(name = hclust$labels[[i]],
       order = hclust$order[[i]],
       cluster = as.numeric(clusters[hclust$labels[[i]]]),
       sensitive = 0)       
}

hclustToTree <- function(hclust, clusters) {
  # hcluste - takes hclust object returned by clustering
  # clusters - a data frame that assosiate each drug to a cluster
  #            where drug names are in rownames.
  if (length(hclust$merge) == 0)
    return(NULL)
  
  merges <- list()
  for (index in 1:nrow(hclust$merge)) {
    left <- hclust$merge[index, 1]
    right <- hclust$merge[index, 2]
    
    if (left < 0)
      left <- createLeafNode(hclust, -left, clusters)
    else
      left <- merges[[left]]
    if (right < 0)
      right <- createLeafNode(hclust, -right, clusters)
    else
      right <- merges[[right]]
    
    if (left$order > right$order) {
      tmp <- left
      left <- right
      right <- tmp
    }
    
    merges[[index]] <- list(
      children = list(
        left,
        right
      ),
      order = left$order
    )
  }
  
  return(merges[nrow(hclust$merge)])
}

render <- function(data) {
  matrix <- as.matrix(data)
  
  rng <- range(matrix)
  
  #   rowClust <- hclust(dist(matrix))
  #   matrix <- matrix[rowClust$order,]
  #   colClust <- hclust(dist(t(matrix)))
  #   matrix <- matrix[,colClust$order]
  tmp <- tempfile()
  png(tmp)
  hm <- heatmap(data, keep.dendro=TRUE)
  dev.off()
  unlink(tmp)
  rowClust <- as.hclust(hm$Rowv)
  colClust <- as.hclust(hm$Colv)
  matrix <- matrix[hm$rowInd, hm$colInd]
  
  rowDend <- toJSON(hclustToTree(rowClust)[[1]], pretty=TRUE)
  colDend <- toJSON(hclustToTree(colClust)[[1]], pretty=TRUE)
  
  matrix <- toJSON(list(data = as.numeric(t(matrix)),
                        dim = dim(matrix),
                        rows = row.names(matrix),
                        cols = names(matrix)))
  
  domain <- toJSON(seq.int(rng[1], rng[2], length.out = 100))
  
  colors <- topo.colors(100)
  colors <- toJSON(sub('FF$', '', colors))
  #colors <- "['yellow', 'green']"
  
  html <- paste(readLines('template.html', warn=FALSE), collapse='\n')
  html <- sub('{{domain}}', domain, html, fixed = TRUE)
  html <- sub('{{colors}}', colors, html, fixed = TRUE)
  html <- sub('{{rowDend}}', rowDend, html, fixed = TRUE)
  html <- sub('{{colDend}}', colDend, html, fixed = TRUE)
  html <- sub('{{matrix}}', matrix, html, fixed = TRUE)
  file <- tempfile('hclust', fileext='.html')
  writeLines(html, file)
  return(file)
}

### END OF JSON SECTION  ###