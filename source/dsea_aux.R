auxDSSDensityEstimate <- function(vector, hh.cells=35, method="gaussian", title="Profile", xax.text="value"){
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
  
  # Plot
  den  <- density(raw, na.rm=TRUE, kernel=method, bw=bw.SJ(raw))
  hh   <- hist(vector, breaks=hh.cells, main = title, xlab = xax.text, probability=TRUE)
  # den$y*diff(hh$mids[1:2])*length(vector)  # fit density probabilities to the Y axis
  yfit <- den$y 
  lines(den$x, yfit, col="blue", lwd=2)
  
  # Flip probs (Y) from negative side of the bump around zero (neg X)
  # to correspondant positives (pos X)
  ind    <- which( den$x < 0 )         
  trans  <- den$y[ind]                 # get Ys for negative Xs
  trans  <- rev(trans)                 # reverse the order of elements (probabilities)
  start  <- max(ind)+1                 # first element of positive Xs to be incremented
  end    <- max(ind)+length(trans)     # last element to be incremented
  insert <- den$y[start:end] + trans   # add reversed vector to the beginnning of pos Xs vector
  
  # Construct a new distribution - compensate for negative values.
  yfit1 <- c( rep(0,length(trans)),            # Zeros
              insert,                          # incremented interval
              den$y[(end+1):length(den$y)] )   # remaining values of the distribution
  lines(den$x, yfit1, col="red", lwd=2)
  
  #print( sum(den$y*diff(den$x[1:2])) )
  print( paste("AUC (red): ", sum(yfit1*diff(den$x[1:2])), sep="") )
  estimate <- data.frame(x=den$x, p=yfit1)
  
  return (estimate)
}

auxAUC <- function(density, thresh){
  # Computes area under curve
  #
  # Args:
  #   density: A vectors 
  #   thresh:   
  # Returns:
  #   aux: [ 0 , 1 ]
  
  library("ROCR")
  aux <- 0
  
  return (aux)
}

