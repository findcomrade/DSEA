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
    
  tryCatch({
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
  trans  <- den$y[ind]                 # get Ys for negative Xs
  trans  <- rev(trans)                 # reverse the order of elements (probabilities)
  start  <- max(ind)+1                 # first element of positive Xs to be incremented
  end    <- max(ind)+length(trans)     # last element to be incremented
  insert <- den$y[start:end] + trans   # add reversed vector to the beginnning of pos Xs vector
  
  # Construct a new distribution - compensate for negative values.
  yfit1 <- c( rep(0,length(trans)),            # Zeros
              insert,                          # incremented interval
              den$y[(end+1):length(den$y)] )   # remaining values of the distribution
  if(graph) lines(den$x, yfit1, col="red", lwd=2)
  
  #print( sum(den$y*diff(den$x[1:2])) )
  print( paste("AUC estimate (red): ", sum(yfit1*diff(den$x[1:2])), sep="") )
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

auxDSSSpecificityScore <- function(probdf, thresh, graph=TRUE){
  # Estimates how a DSS is specific in comparison to the whole set.
  # It basically calculates area under curve restricted by thresh.
  #
  # Args:
  #   probdf:     a data frame - probability density function
  #            REQUIRE probdf$x, probdf$p - prob vector;
  #   thresh:    
  # Returns:
  #   score: [ 0 , 1 ]; 1 - the most specific.
  score <- 0
  
  if (thresh < 1) return (0)
  
  ind   <- which(probdf$x <= thresh)
  span  <- diff(probdf$x[1:2])
  score <- sum( probdf$p[ind] * span )
  
  if (score > 1) score <- 1
  
  score <- round(score*100, digits = 3)
  
  # Graphics
  if (graph){
    plot(probdf$x, probdf$p, type = "o", pch = 20, main="Specificity Score", xlab="DSS", ylab="Prob density",
         panel.first = auxPolyCurve(probdf$x, probdf$p, n=1000, from = 0, to = thresh, col = c("orange"), border = "black"))  
    axis(side=1, at=thresh, labels=TRUE, tick=TRUE, col="blue", col.axis="blue",)
    mtext(paste("Score:", score, sep="\t"), side=3, line=-3, outer=FALSE, col="blue")
  }
  
  return (score)
}