setwd("/home/comrade/Ubuntu One/DSEA")


###       Code Section            ###
barplot_fun <- function(xpr, cellline, ...){
  #xpr <- xpr[!is.na(c(xpr))]
  normed <- rep(0.5,length(xpr))
  xpr <- sweep(xpr, 1, normed , "-")
  xpr <- xpr[order(xpr,decreasing=FALSE),]
  xpr[(xpr > 10)] <- 10
  par(mar=c(10,5,4,5))
  #rainbow_hcl(1, c = 100, l = 50,start=240,end=280)
  col <- rainbow_hcl(2, c = 100, l = 10,start=210,end=225)
  myclorz <- ifelse(xpr>0, col[2], rainbow_hcl(1, c = 100, l = 50,start=240,end=280))
  main_text <- paste("Response profile of ",cellline, " cell line",sep="")
  barplot(xpr, axisnames =TRUE,names.arg=names(xpr),horiz=FALSE,las=2,cex.axis = 0.8, cex.names = 0.3,space=1,axis.lty=1,lwd=2,xaxs="i",yaxs="i",col= myclorz,main= main_text,ylab="Fraction of control")
  abline(h=0,lty=1,lwd=4, col = "black")
  
  legend("topleft",legend=c(" < 0.5"," >= 0.5"), cex=1,col="",fill=c(col[2], rainbow_hcl(1, c = 100, l = 50,start=240,end=280)),title="Fraction of control(IC50)",xpd=NA)
  
  #	plot(xpr,type = "h", lwd = 12, col = myclorz, xlab = "", ylab = "Fraction of control",main = main_text, xaxt="n")
  #	axis(side=1, at=1:length(names(xpr)),labels=names(xpr), lty=3, las=2,cex.axis=1.5, tck=-.01)
  #	abline(h=0,lty=1,lwd=4, col = "black")
  #	legend("topright",legend=c(" < 0.5"," >= 0.5"), inset=c(-0.3,-0.1),cex=1,col="",fill=c(col[2], rainbow_hcl(1, c = 100, l = 50,start=240,end=280)),title="Fraction of control(IC50)",xpd=NA)
  #       inset=c(-0.3,-0.1)
}


#####################################
### Parameters Definition Section ###
#####################################
Database <- "sanger.RData"
file_name <- "plot"
Output_type <- ".pdf"
Celllines <- c("ES3")


#####################################
###       Assembly Section        ###
#####################################
#Set the initialisation variables
library(vcd)
load(Database)
celllines <- Celllines
filename = paste(file_name,Output_type,sep="")
if(Output_type == ".pdf"){
  pdf(file = filename,width=11,title = "Visual overall view of drug inihibition")
}else{
  png(filename = filename,width=11,title = "Visual overall view of drug inihibition")
}

dataSet <-  sangerSet[,(as.character(sangerSet$Cell.Line) %in% Celllines)]
common_celllines <- as.character(dataSet$Cell.Line)
res_plots <- lapply(Celllines,function(cellline){
  if(cellline %in% as.character(dataSet$Cell.Line)){
    xpr1 <- exprs(dataSet[1:139,(as.character(dataSet$Cell.Line) %in% cellline)])
    #xpr1 <- xpr1[!is.na(c(xpr1)),]
    #Clean up the row names	 			
    #				dup_rmv <- grep("Camptothecin_IC_50.1",rownames(xpr1))
    #				xpr1 <- xpr1[-dup_rmv,]
    rownames(xpr1) <- gsub("_IC_50$","",rownames(xpr1))
    barplot_fun(xpr1, cellline)
  }else{
    plot.new(); 
    title(main = paste("Missing data for cell line ", cellline,sep=""))
  }			
  
})

dev.off()
# pdfopt(filename)
