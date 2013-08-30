#
# Contains supplementary functions for DSEA module;
# 
# Author: Dmitrii Bychkov, FIMM 2013
# (dmitrii.bychkov@helsinki.fi)
#######################################################

buildEnrichmentTable <- function(cell.line){
  # each 
  target.drugs <- getTopDrugs(cell.line)
  target.senscount <- length(target.drugs)  # countScreenedDrugs(cell.line)
  
  Cell.Line <- c();  top.count <- c();  screened.count <- c();  SCORE <- c(); p.value <- c();
  # cell.lines  comes from loadDataSet()
  for (c.line in cell.lines){
    if (c.line != cell.line){
      Cell.Line <- c( Cell.Line, paste(c.line,"",sep="") )
      ref.drugs <- getTopDrugs(c.line)
      top.count <- c(top.count, length(ref.drugs))
      screened.count <- c(screened.count, countScreenedDrugs(paste(c.line,"",sep="")))
      overlap <- sum(ref.drugs %in% target.drugs)
      SCORE <- c(SCORE, overlap)
      p.value <- c( p.value, pValue(target.senscount, length(ref.drugs), overlap) )
    }
  }
  ENRICHMENT <- data.frame(Cell.Line, top.count, screened.count, SCORE, p.value)  # create df
  p <- p.value
  ENRICHMENT$bf.correction <- p.adjust(p, "bonferroni")
  ENRICHMENT <- ENRICHMENT[ ENRICHMENT$SCORE != 0, ]  # trim SCORE == 0
  ENRICHMENT <- ENRICHMENT[ order(ENRICHMENT[,6]), ]  # order by bonferroni corrected p-value
  ENRICHMENT$perc <- as.numeric(format(round(ENRICHMENT[,"SCORE"]/ENRICHMENT[,"top.count"] *100, 2), nsmall=2))  # calculate percentage
  ENRICHMENT$significant <- ENRICHMENT$bf.correction <= 0.05
  
  return (ENRICHMENT)
}

pValue <- function(sensitive.target, sensitive.reference, overlap){
  # Performs Fisher's exact test for testing the null of independence of rows and columns;
  # sensitive.target <- 9;  sensitive.reference <- 20;  overlap <- 6;
  if (!overlap){
    return (0)
  }
  
  total.count <- length(getDrugNames(ds.name))
  nonsense.target <- sensitive.target - overlap
  nonsense.reference <- total.count - sensitive.reference
    
  DrugTesting <-
    matrix( c((overlap-1), nonsense.target, sensitive.reference, nonsense.reference),
            nrow = 2,
            dimnames = list( Guess = c("Sense", "NonSense"),
                             Truth = c("Target", "Ref")) )
  
  test <- fisher.test(DrugTesting, alternative = "greater")
  return ( as.numeric(test["p.value"]) )
}

countScreenedDrugs <- function(cell.line){
  MATRIX <- as.matrix(FRAME)
  drugs.ordered <- FRAME$DRUG[as.numeric(MATRIX[,paste(cell.line,".order",sep="")])]
  values.ordered <- MATRIX[as.numeric(MATRIX[,paste(cell.line,".order",sep="")]),paste(cell.line,".value",sep="")]
  
  top.count <- 0
  for (x in values.ordered){
    if (!is.na(x)) {
      top.count <- top.count + 1    
    }
  }
  drugs.top <- drugs.ordered[1:top.count] 
  
  return (length(drugs.top))
}

plotDrugSensitivity <- function(cell.line){
  # TODO: implement choice of the cell.line
  cell.line <- "SW620" #  "T47D"
  cutoff <- length(getTopDrugs(cell.line))
  #MATRIX <- as.matrix(FRAME)
  #orderedData <- FRAME$DRUG[as.numeric(MATRIX[,paste(cell.line,".order",sep="")])]
  orderedData <- FRAME$SW620.value[FRAME$SW620.order, ]
  top.drugs <- orderedData[1:cutoff]
  
  col = c("mistyrose", "lightblue")
  myclorz <- ifelse(orderedData %in% top.drugs, col[1], col[2])
  barplot(orderedData, axisnames =TRUE,names.arg=names(orderedData),horiz=FALSE,col= myclorz,
          las=2,cex.axis = 0.8, cex.names = 0.3,space=1,axis.lty=1,lwd=2,xaxs="i",yaxs="i",ylab="Sanger Custom IC50", main ="Drug Sensitivity for T47D Cell Line across 139 Drugs")
}

loadDataSet <- function(dataset.name){
  # Loads and convert data from a given data set to
  # the data frame of a fixed layout.
  
  if (dataset.name == "sanger"){
    load("../datasets/sanger.RData")
    
    drug.names <- featureNames(sangerSet)
    cell.lines <- sampleNames(sangerSet)
    pd <- pData(sangerSet)
    
    drugs.span <- 1:139
    
    # trim "_IC_50$" in drug names
    DRUG <- gsub("_IC_50$","",drug.names[drugs.span])
    
    # FRAME: Drugs go in rows; 
    # Each cell line represented by two cols: 1) response value  2) ordering number
    # "2)" is used to arrange drugs in descending order of response (smaller value - better response)
    FRAME <- data.frame(DRUG)
    for (line in cell.lines){
      FRAME$newcol.value <- exprs(sangerSet[drugs.span,(as.character(sangerSet$Cell.Line) %in% line)])  
      FRAME$newcol.value[(FRAME$newcol.value > 10)] <- 10 # limit by 10 ?!
      FRAME$newcol.order <- order(FRAME$newcol.value,decreasing=FALSE)
      
      # fix col names according to the cell line
      names(FRAME)[length(FRAME)-1] <- paste(line, ".value", sep="")
      names(FRAME)[length(FRAME)] <- paste(line, ".order", sep="")
    }  # end for
  }  # endif
  
  return (FRAME)
}

getDrugNames <- function(dataset.name){
  drug.names <- c()
  
  if (dataset.name == "sanger"){
    load("../datasets/sanger.RData") 
    drug.names <- featureNames(sangerSet)
    drugs.span <- 1:139
    drug.names <- gsub("_IC_50$","",drug.names[drugs.span])  
  }
  
  return (drug.names)
}

getCellLines <- function(dataset.name){
  drug.names <- c()
  
  if (dataset.name == "sanger"){
    load("../datasets/sanger.RData") 
    cell.lines <- sampleNames(sangerSet)
  }
  
  return (cell.lines)
}

getTopDrugs <- function(cell.line){
  threshold <- -1
  MATRIX <- as.matrix(FRAME)
  drugs.ordered <- FRAME$DRUG[as.numeric(MATRIX[,paste(cell.line,".order",sep="")])]
  values.ordered <- MATRIX[as.numeric(MATRIX[,paste(cell.line,".order",sep="")]),paste(cell.line,".value",sep="")]
  top.count <- 0
  for (x in values.ordered){
    if ( (as.numeric(x) < threshold) & (!is.na(x)) ){
      top.count <- top.count + 1    
    }
  }
  drugs.top <- drugs.ordered[1:top.count] 
  
  return (drugs.top)
}

genOverlapPlot <- function(cell_line1, cell_line2, target_top, ref_top){
  MATRIX <- as.matrix(FRAME)
  target_CL <- cell_line1
  other_CL <- cell_line2
  # get top drugs for target CL
  ref_set <- FRAME$DRUG[as.numeric(MATRIX[,paste(target_CL,".order",sep="")])]
  ref_set <- ref_set[1:target_top]
  # get top values for target CL
  ref_val <- as.numeric(MATRIX[as.numeric(MATRIX[,paste(target_CL,".order",sep="")]),paste(target_CL,".value",sep="")])
  ref_val <- ref_val[1:target_top]
  
  # get top drugs for another CL
  tmp <- FRAME$DRUG[as.numeric(MATRIX[,paste(other_CL,".order",sep="")])]
  tmp <- tmp[1:ref_top]
  # get top values for another CL
  tmp_val <- as.numeric(MATRIX[as.numeric(MATRIX[,paste(other_CL,".order",sep="")]),paste(other_CL,".value",sep="")])
  tmp_val <- tmp_val[1:ref_top]
  
  dfTARGET <- data.frame(CL=rep(target_CL,length(ref_set)), Drug=ref_set, Sensitivity=ref_val, Common=ref_set %in% tmp)
  dfTMP <- data.frame(CL=rep(other_CL,length(tmp)), Drug=tmp, Sensitivity=tmp_val, Common=tmp %in% ref_set)
  df <- rbind(dfTARGET, dfTMP)
  
  #lst <- df$Drug
  #for (dr in lst){
  #  if (length(df[df[,"Drug"]==dr,"Common"]) == 1){
  #    cl <- ifelse(df[df[,"Drug"]==dr,"CL"] == target_CL, other_CL, target_CL)
  #    df <- rbind(data.frame(CL=cl, Drug=dr, Sensitivity=as.factor(0), Common=FALSE), df)
  #  }
  #}
  
  df$Drugs <- reorder(df$Drug, rank(df$Sensitivity))
  
  #png( paste( "reports/", figureFile1, sep="" ), width=750, height=500 );
  ggp <- ggplot(df, aes(x=Drugs, y=Sensitivity,fill=CL)) + geom_bar(stat="identity", position = "dodge") + 
    opts(axis.text.x = theme_text(angle = 90))
  #dev.off()
  
  return(ggp)
}

plotTissueFreq <- function(){
  unique.tissues <- table(pd$Tissue)
  unique.tissues <- as.data.frame(unique.tissues)
  colnames(unique.tissues) <- c("Tissue", "Freq")
  unique.tissues$Tissue <- reorder(unique.tissues$Tissue, rank(unique.tissues$Freq))
  
  ggp <- ggplot(unique.tissues, aes(x=Tissue, y=Freq)) + geom_bar(stat="identity", position = "dodge") + 
    opts(axis.text.x = theme_text(angle = 90))
  
  return (ggp)
}

plotEnrishedTissue <- function(){
  # get enriched cell lines
  enriched.celllines <- ENRICHMENT[as.logical(ENRICHMENT$significant),1]
  
  # collect tissues
  tissue.collection <- c()
  for (item in enriched.celllines){
    tissue.collection <- c( tissue.collection, as.character(pd[pd["Cell.Line"] == as.character(item),"Tissue"]) )
  }
  
  unique.tissues <- table(tissue.collection)
  unique.tissues <- as.data.frame(unique.tissues)
  colnames(unique.tissues) <- c("Tissue", "Freq")
  unique.tissues$Tissue <- reorder(unique.tissues$Tissue, rank(unique.tissues$Freq))
  
  ggp <- ggplot(unique.tissues, aes(x=Tissue, y=Freq)) + geom_bar(stat="identity", position = "dodge") + 
    opts(axis.text.x = theme_text(angle = 90))
  
  return (ggp)
}

plotDiagnosisFreq <- function(){
  unique.diagnos <- table(pd$Cancer.Type)
  unique.diagnos <- as.data.frame(unique.diagnos)
  colnames(unique.diagnos) <- c("Diagnosis", "Freq")
  unique.diagnos$Diagnosis <- reorder(unique.diagnos$Diagnosis, rank(unique.diagnos$Freq))
  
  ggp <- ggplot(unique.diagnos, aes(x=Diagnosis, y=Freq)) + geom_bar(stat="identity", position = "dodge") + 
    opts(axis.text.x = theme_text(angle = 90))
  
}

plotEnrichedDiagnos <- function(){
  # get enriched cell lines
  enriched.celllines <- ENRICHMENT[as.logical(ENRICHMENT$significant),1]
  
  # collect tissues
  diagnosis.collection <- c()
  for (item in enriched.celllines){
    diagnosis.collection <- c( diagnosis.collection, as.character(pd[pd["Cell.Line"] == as.character(item),"Cancer.Type"]) )
  } 
  
  unique.diagnos <- table(diagnosis.collection)
  unique.diagnos <- as.data.frame(unique.diagnos)
  colnames(unique.diagnos) <- c("Diagnosis", "Freq")
  unique.diagnos$Diagnosis <- reorder(unique.diagnos$Diagnosis, rank(unique.diagnos$Freq))
  
  ggp <- ggplot(unique.diagnos, aes(x=Diagnosis, y=Freq)) + geom_bar(stat="identity", position = "dodge") + 
    opts(axis.text.x = theme_text(angle = 90))
  
}