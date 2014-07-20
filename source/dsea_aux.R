auxDendroPlot <- function(hcl.out, cut.out, cluster.count, threshold){
  
  op <- par(bg = "#DDE3CA", cex=0.2, mar=c(5,5,5,3), oma=c(3,3,3,3)) 
  #, mfrow = c(1, 2))
  
  labelColors  <- rainbow( cluster.count)
  clusMember   <- cut.out
  
  colLab       <- function(n) {
    if (is.leaf(n)) {
      a       <- attributes(n)
      labCol  <- labelColors[clusMember[which(names(clusMember) == a$label)]]
      attr(n, "nodePar") <- c(a$nodePar, lab.col = labCol)
    }
    n
  }
  
  clusDendro = dendrapply( as.dendrogram(hcl.out), colLab )
  
  table(mycl, cutree(hclust(as.dist(1-corr.Drugs), method="ward"), k=15, h=3) )
  plot(clusDendro, col = "#487AA1", col.main = "#45ADA8", col.lab = "#7C8071", col.axis = "#F38630", 
       lwd = 3, lty = 3, sub = "", hang = -1, axes = FALSE)
  
  axis(side = 2, at = seq(0, max(hcl.out$height+1), 5), col = "#F38630", labels = FALSE, lwd = 2)
  # add text in margin
  mtext(seq(0, max(hcl.out$height+1), 5), side = 2, at = seq(0, max(hcl.out$height+1), 5), line = 1, col = "#A38630", las = 2, cex=0.6)
  abline(h=threshold, col="#487AA1")
  plot(hcl.out$height, col = "#487AA1", col.main = "#45ADA8", col.lab = "#7C8071", col.axis = "#F38630", 
       lwd = 3, lty = 3, sub = "", axes = FALSE, ylab="Ward's Merging Cost function", xlab="Iteration", cex.lab=2)
  axis(side = 4, at = seq(0, max(hcl.out$height+1), 5), col = "#F38630", labels = FALSE, lwd = 2)
  mtext(seq(0, max(hcl.out$height+1), 5), side = 4, at = seq(0, max(hcl.out$height+1), 5), line = 1, col = "#A38630", las = 2, cex=0.6)
  abline(h=threshold, col="#487AA1")
  
  par(op)  
}

auxDSSDensityEstimate <- function(vector, hh.cells=35, method="gaussian", title="Profile", xax.text="value", graph=TRUE){
  # Plots kernel density estimate of the data along with a histogram.
  # Red estimate adjusts for DSS in the range [0,inf) by flipping the left (negative)
  # side of the distribution bump around zero.
  #
  # Args:
  #   vector: the data 
  #   hh.cells: number of cells for the histogram
  #   method: a string giving the smoothing kernel to be used (see density() doc)
  # Returns:
  #   estimate: a data frame 
  raw   <- as.numeric(vector[!is.na(vector)])
    
  den <- tryCatch({
    density(raw, na.rm=TRUE, kernel=method, bw=bw.SJ(raw))  
  }, error = function(e){ 
    # print(e)
    print('Sample is too sparse => Estimate is not reliable!')
    return(-1) 
  })
  
  # The following will be exected when bw.SJ() fails due to sparse data:
  # We fix bandwidth == 0.31 (  based on intuition :)  )
  if (class(den) != "density") den <- density(raw, na.rm=TRUE, kernel=method, bw=0.31)
  
  # Plot
  if(graph){
    hh <- hist(vector, breaks=hh.cells, plot=graph, probability=TRUE, main = title, xlab = xax.text, ylim=c(0, max(den$y*2.1)))
    # den$y*diff(hh$mids[1:2])*length(vector)  # fit density probabilities to the Y axis in case probability=FALSE
    yfit <- den$y 
    lines(den$x, yfit, col="blue", lwd=2)
  } else{
    hh <- hist(vector, breaks=hh.cells, plot=graph)
  }
  
  # Flip probs (Y) from negative side of the bump around zero (neg X)
  # to correspondant positives (pos X)
  ind    <- which( den$x < 0 ) 
  if( length(ind) ){
    # if there is a negative tail of distribution
    trans  <- den$y[ind]                 # get Ys for negative Xs
    trans  <- rev(trans)                 # reverse the order of elements (probabilities)
    start  <- max(ind)+1                 # first element of positive Xs to be incremented
    end    <- max(ind)+length(trans)     # last element to be incremented
    insert <- den$y[start:end] + trans   # add reversed vector to the beginnning of pos Xs vector
    
    # Construct a new distribution - compensate for negative values.
    yfit1 <- c( rep(0,length(trans)),            # Zeros
                insert,                          # incremented interval
                den$y[(end+1):length(den$y)] )   # remaining values of the distribution 
  } else{
    # no need to adjust for negatives
    yfit1 <- den$y
  }
  
  if(graph) lines(den$x, yfit1, col="red", lwd=2)
  
  #print( sum(den$y*diff(den$x[1:2])) )
  control.sum <- sum(yfit1*diff(den$x[1:2]))
  if(control.sum > 0.95 & control.sum < 1.05){
    print( paste("AUC estimate (red): ", control.sum, sep="") )  
  } else{
    msg = paste('AUC: ', control.sum, ' => Estimate is not reliable!')
    print(msg)
    return(-1)
  }
  
  estimate <- data.frame(x=den$x, p=yfit1)
  
  return (estimate)
}

auxPolyCurve <- function(x, y, from, to, n = 50, miny, col = "red", border = col){
  drawPoly <- function(fun, from, to, n = 50, miny, col, border) {
    Sq <- seq(from = from, to = to, length = n)
    polygon(x = c(Sq[1], Sq, Sq[n]),y = c(miny, fun(Sq), miny), col = col, border = border)
  }
  
  lf <- length(from)
  stopifnot(identical(lf, length(to)))
  if(length(col) != lf)
    col <- rep(col, length.out = lf)
  if(length(border) != lf)
    border <- rep(border, length.out = lf)
  if(missing(miny))
    miny <- min(y)
  interp <- approxfun(x = x, y = y)
  mapply(drawPoly, from = from, to = to, col = col, border = border,
                          MoreArgs = list(fun = interp, n = n, miny = miny))
  invisible()
}

auxDSSSpecificityScore <- function(probdf, thresh, control=0, graph=TRUE, title="Profile"){
  # Estimates how a DSS is specific in comparison to the whole set.
  # It basically calculates area under curve restricted by thresh.
  #
  # Args:
  #   probdf:     a data frame - probability density function
  #            REQUIRE probdf$x, probdf$p - prob vector;
  #   thresh:     dss value for a sample.
  #   control:    some reference dss; that is MEAN value so far.
  # Returns:
  #   score: [ 0 , 100 ]; 100 - the most specific.
  score <- 0
  
  if (thresh < 1) return (0)
  
  # First calculate Score for the sample
  ind   <- which(probdf$x <= thresh)
  span  <- diff(probdf$x[1:2])
  score <- sum( probdf$p[ind] * span )
  
  if (score > 1) score <- 1
  score <- round(score*100, digits = 3)
  
  if (control){
    # Now calculate Score for the control
    ind   <- which(probdf$x <= control)
    span  <- diff(probdf$x[1:2])
    control.score <- sum( probdf$p[ind] * span )
    
    if (control.score > 1) control.score <- 1
    control.score <- round(control.score*100, digits = 3)
  }
  
  
  # Graphics
  if (graph){
    plot(probdf$x, probdf$p, main=title, type = "o", pch = 20, xlab="DSS", ylab="Prob density",
         panel.first = auxPolyCurve(probdf$x, probdf$p, n=1000, from = 0, to = thresh, col = c("orange"), border = "black"))  
    axis(side=1, at=thresh, labels=TRUE, tick=TRUE, col="blue", col.axis="blue",)
    mtext(paste("Sample SpS: ", score, "(", thresh, ")",sep=" "), side=3, line=0, outer=FALSE, col="blue") #, cex=0.3
    if (control) mtext(paste("Control SpS: ", control.score, "(", control, ")",sep=" "), side=3, line=-4, outer=FALSE, col="blue")
    
  }
  
  return (score)
}

auxSpecificityMatrix <- function(dss.table){
  # Args:
  #   dss.table: a dss matrix with drugs in rows. 
  # Returns:
  #   spec.table: a matrix of specificty scores.
  drug.names    <- rownames(dss.table)
  sample.names  <- colnames(dss.table)
  
  spec.table    <- matrix( data=NA, nrow=length(drug.names), ncol=length(sample.names) )
  rownames(spec.table)  <- drug.names
  colnames(spec.table)  <- sample.names
  
  for(drug in drug.names){  
    drug.profile <- dss.table[drug,]
    prob.df <- auxDSSDensityEstimate(vector=drug.profile, graph=FALSE)
    
    if(class(prob.df) == "data.frame"){
      for(sample in sample.names){
        dss <- dss.table[drug, sample]
        if( !is.na(dss) ){
          spec.table[drug, sample] <- auxDSSSpecificityScore(probdf=prob.df, thresh=dss, graph=FALSE)  
        }
      }  # over samples (cols)
    }  # end od IF
  }  # over drugs (rows)
  
  return(spec.table)  
}

auxTissueSpecificity <- function(dss.table, tissue.annotat){
  #
  # Args:
  #   dss.table: a dss matrix with drugs in rows;
  #   tissue.annotat: a data frame Nx3
  # Returns:
  #   spec.table: a matrix of specificty scores.
  
  drug.names    <- rownames(dss.table)
  sample.names  <- colnames(dss.table)
  
  tissue.names  <- unique(tissue.annotat$Sample.Origin)
  
  # First we create Drugs vs Tissue matrix with mean dss values
  tissue.mean.dss  <- matrix( data=NA, nrow=length(drug.names), ncol=length(tissue.names) )
  rownames(tissue.mean.dss)  <- drug.names
  colnames(tissue.mean.dss)  <- tissue.names
  
  for(tissue in tissue.names){
    # get sample names for a tissue
    ind  <- which(tissue.annotat[,"Sample.Origin"] == tissue)
    
    tissue.samples  <- tissue.annotat[ind,"Sample.Name"]
    tissue.samples.ind  <- which(sample.names %in% tissue.samples)
    # Now we iterate over drugs    
    for(drug in drug.names){
      vector  <- as.numeric( dss.table[drug, tissue.samples.ind] )      
      tissue.mean.dss[drug, tissue]  <- mean(vector)
    }
      
  }
  
  # Second part operates on tissue.mean.dss
  spec.table    <- matrix( data=NA, nrow=length(drug.names), ncol=length(tissue.names) )
  rownames(spec.table)  <- drug.names
  colnames(spec.table)  <- tissue.names
  
  for(drug in drug.names){  
    drug.profile <- dss.table[drug,]
    prob.df <- auxDSSDensityEstimate(vector=drug.profile, graph=FALSE)
    
    if(class(prob.df) == "data.frame"){
      for(tissue in tissue.names){
        dss <- tissue.mean.dss[drug, tissue]
        if( !is.na(dss) ){
          spec.table[drug, tissue] <- auxDSSSpecificityScore(probdf=prob.df, thresh=dss, graph=FALSE)  
        }
      }  # over samples (cols)
    }  # end od IF
  }  # over drugs (rows)
  
  return(spec.table)  
}

auxTissueDSSMeans <- function(dss.table, tissue.annotat){
  #
  # Args:
  #   dss.table: a dss matrix with drugs in rows;
  #   tissue.annotat: a data frame Nx3
  # Returns:
  #   spec.table: a matrix of specificty scores.
  
  drug.names    <- rownames(dss.table)
  sample.names  <- colnames(dss.table)
  
  tissue.names  <- unique(tissue.annotat$Sample.Origin)
  
  # First we create Drugs vs Tissue matrix with mean dss values
  tissue.mean.dss  <- matrix( data=NA, nrow=length(drug.names), ncol=length(tissue.names) )
  rownames(tissue.mean.dss)  <- drug.names
  colnames(tissue.mean.dss)  <- tissue.names
  
  for(tissue in tissue.names){
    # get sample names for a tissue
    ind  <- which(tissue.annotat[,"Sample.Origin"] == tissue)
    
    tissue.samples  <- tissue.annotat[ind,"Sample.Name"]
    tissue.samples.ind  <- which(sample.names %in% tissue.samples)
    # Now we iterate over drugs    
    for(drug in drug.names){
      vector  <- as.numeric( dss.table[drug, tissue.samples.ind] )      
      tissue.mean.dss[drug, tissue]  <- mean(vector)
    }
    
  }
  return(tissue.mean.dss)
}