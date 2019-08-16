library(HGNChelper)
library(pamr)
library(plyr)
library(pROC)
library(samr)

TCGA_GENERIC_CleanUpSampleNames <-function(GEN_Data,IDlength=12) {     
  SampleNames=colnames(GEN_Data)
  SampleNamesShort=as.character(apply(as.matrix(SampleNames),2,substr,1,IDlength))
  if (length(SampleNamesShort)!=length(unique(SampleNamesShort))) {
    # remove the doubles           
    Counts=table(SampleNamesShort)
    Doubles=rownames(Counts)[which(Counts>1)]
    
    cat("Removing doubles for",length(Doubles),"samples.\n")
    for(i in 1:length(Doubles)) {                         
      CurrentDouble=Doubles[i]          
      pos=grep(CurrentDouble,SampleNames)
      #GEN_Data[1:10,pos]
      #cor(GEN_Data[,pos])
      GEN_Data=GEN_Data[,-pos[2:length(pos)]]     
      SampleNames=colnames(GEN_Data) # need to update samplenames because pos is relative to this
    }
    SampleNames=colnames(GEN_Data)
    SampleNamesShort=as.character(apply(as.matrix(SampleNames),2,substr,1,IDlength))
    
    # now set the samplenames
    colnames(GEN_Data)=SampleNamesShort
  } else {
    colnames(GEN_Data)=SampleNamesShort     
  }     
  return(GEN_Data)
}

TCGA_GENERIC_GeneFiltering<- function(Type,MAdata,Percentage) {
  
  switch(Type,
         Variance={
           GeneVariances=rowVars(MAdata)            
           tmpResult=sort(GeneVariances,decreasing=TRUE)
           SortedGenes=names(tmpResult)
           tmpNrGenes=round(length(rownames(MAdata))*Percentage/100)
           MAdata_Filtered=MAdata[SortedGenes[1:tmpNrGenes],]
         },
         MAD={
           GeneVariances=rowMads(MAdata)            
           names(GeneVariances)=rownames(MAdata)
           tmpResult=sort(GeneVariances,decreasing=TRUE)
           SortedGenes=names(tmpResult)
           tmpNrGenes=round(length(rownames(MAdata))*Percentage/100)
           MAdata_Filtered=MAdata[SortedGenes[1:tmpNrGenes],]
         }
  )
  
  return(MAdata_Filtered)     
}


#Heatmap 3
heatmap.3 <- function(x,
                      Rowv = TRUE, Colv = if (symm) "Rowv" else TRUE,
                      distfun = dist,
                      hclustfun = hclust,
                      dendrogram = c("both","row", "column", "none"),
                      symm = FALSE,
                      scale = c("none","row", "column"),
                      na.rm = TRUE,
                      revC = identical(Colv,"Rowv"),
                      add.expr,
                      breaks,
                      symbreaks = max(x < 0, na.rm = TRUE) || scale != "none",
                      col = "heat.colors",
                      colsep,
                      rowsep,
                      sepcolor = "white",
                      sepwidth = c(0.05, 0.05),
                      cellnote,
                      notecex = 1,
                      notecol = "cyan",
                      na.color = par("bg"),
                      trace = c("none", "column","row", "both"),
                      tracecol = "cyan",
                      hline = median(breaks),
                      vline = median(breaks),
                      linecol = tracecol,
                      margins = c(5,5),
                      ColSideColors,
                      RowSideColors,
                      side.height.fraction=0.3,
                      cexRow = 0.2 + 1/log10(nr),
                      cexCol = 0.2 + 1/log10(nc),
                      labRow = NULL,
                      labCol = NULL,
                      key = TRUE,
                      keysize = 1.5,
                      density.info = c("none", "histogram", "density"),
                      denscol = tracecol,
                      symkey = max(x < 0, na.rm = TRUE) || symbreaks,
                      densadj = 0.25,
                      main = NULL,
                      xlab = NULL,
                      ylab = NULL,
                      lmat = NULL,
                      lhei = NULL,
                      lwid = NULL,
                      ColSideColorsSize = 1,
                      RowSideColorsSize = 1,
                      KeyValueName="Value",...){
  
  invalid <- function (x) {
    if (missing(x) || is.null(x) || length(x) == 0)
      return(TRUE)
    if (is.list(x))
      return(all(sapply(x, invalid)))
    else if (is.vector(x))
      return(all(is.na(x)))
    else return(FALSE)
  }
  
  x <- as.matrix(x)
  scale01 <- function(x, low = min(x), high = max(x)) {
    x <- (x - low)/(high - low)
    x
  }
  retval <- list()
  scale <- if (symm && missing(scale))
    "none"
  else match.arg(scale)
  dendrogram <- match.arg(dendrogram)
  trace <- match.arg(trace)
  density.info <- match.arg(density.info)
  if (length(col) == 1 && is.character(col))
    col <- get(col, mode = "function")
  if (!missing(breaks) && (scale != "none"))
    warning("Using scale=\"row\" or scale=\"column\" when breaks are",
            "specified can produce unpredictable results.", "Please consider using only one or the other.")
  if (is.null(Rowv) || is.na(Rowv))
    Rowv <- FALSE
  if (is.null(Colv) || is.na(Colv))
    Colv <- FALSE
  else if (Colv == "Rowv" && !isTRUE(Rowv))
    Colv <- FALSE
  if (length(di <- dim(x)) != 2 || !is.numeric(x))
    stop("`x' must be a numeric matrix")
  nr <- di[1]
  nc <- di[2]
  if (nr <= 1 || nc <= 1)
    stop("`x' must have at least 2 rows and 2 columns")
  if (!is.numeric(margins) || length(margins) != 2)
    stop("`margins' must be a numeric vector of length 2")
  if (missing(cellnote))
    cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
  if (!inherits(Rowv, "dendrogram")) {
    if (((!isTRUE(Rowv)) || (is.null(Rowv))) && (dendrogram %in%
                                                 c("both", "row"))) {
      if (is.logical(Colv) && (Colv))
        dendrogram <- "column"
      else dedrogram <- "none"
      warning("Discrepancy: Rowv is FALSE, while dendrogram is `",
              dendrogram, "'. Omitting row dendogram.")
    }
  }
  if (!inherits(Colv, "dendrogram")) {
    if (((!isTRUE(Colv)) || (is.null(Colv))) && (dendrogram %in%
                                                 c("both", "column"))) {
      if (is.logical(Rowv) && (Rowv))
        dendrogram <- "row"
      else dendrogram <- "none"
      warning("Discrepancy: Colv is FALSE, while dendrogram is `",
              dendrogram, "'. Omitting column dendogram.")
    }
  }
  if (inherits(Rowv, "dendrogram")) {
    ddr <- Rowv
    rowInd <- order.dendrogram(ddr)
  }
  else if (is.integer(Rowv)) {
    hcr <- hclustfun(distfun(x))
    ddr <- as.dendrogram(hcr)
    ddr <- reorder(ddr, Rowv)
    rowInd <- order.dendrogram(ddr)
    if (nr != length(rowInd))
      stop("row dendrogram ordering gave index of wrong length")
  }
  else if (isTRUE(Rowv)) {
    Rowv <- rowMeans(x, na.rm = na.rm)
    hcr <- hclustfun(distfun(x))
    ddr <- as.dendrogram(hcr)
    ddr <- reorder(ddr, Rowv)
    rowInd <- order.dendrogram(ddr)
    if (nr != length(rowInd))
      stop("row dendrogram ordering gave index of wrong length")
  }
  else {
    rowInd <- nr:1
  }
  if (inherits(Colv, "dendrogram")) {
    ddc <- Colv
    colInd <- order.dendrogram(ddc)
  }
  else if (identical(Colv, "Rowv")) {
    if (nr != nc)
      stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
    if (exists("ddr")) {
      ddc <- ddr
      colInd <- order.dendrogram(ddc)
    }
    else colInd <- rowInd
  }
  else if (is.integer(Colv)) {
    hcc <- hclustfun(distfun(if (symm)
      x
      else t(x)))
    ddc <- as.dendrogram(hcc)
    ddc <- reorder(ddc, Colv)
    colInd <- order.dendrogram(ddc)
    if (nc != length(colInd))
      stop("column dendrogram ordering gave index of wrong length")
  }
  else if (isTRUE(Colv)) {
    Colv <- colMeans(x, na.rm = na.rm)
    hcc <- hclustfun(distfun(if (symm)
      x
      else t(x)))
    ddc <- as.dendrogram(hcc)
    ddc <- reorder(ddc, Colv)
    colInd <- order.dendrogram(ddc)
    if (nc != length(colInd))
      stop("column dendrogram ordering gave index of wrong length")
  }
  else {
    colInd <- 1:nc
  }
  retval$rowInd <- rowInd
  retval$colInd <- colInd
  retval$call <- match.call()
  x <- x[rowInd, colInd]
  x.unscaled <- x
  cellnote <- cellnote[rowInd, colInd]
  if (is.null(labRow))
    labRow <- if (is.null(rownames(x)))
      (1:nr)[rowInd]
  else rownames(x)
  else labRow <- labRow[rowInd]
  if (is.null(labCol))
    labCol <- if (is.null(colnames(x)))
      (1:nc)[colInd]
  else colnames(x)
  else labCol <- labCol[colInd]
  if (scale == "row") {
    retval$rowMeans <- rm <- rowMeans(x, na.rm = na.rm)
    x <- sweep(x, 1, rm)
    retval$rowSDs <- sx <- apply(x, 1, sd, na.rm = na.rm)
    x <- sweep(x, 1, sx, "/")
  }
  else if (scale == "column") {
    retval$colMeans <- rm <- colMeans(x, na.rm = na.rm)
    x <- sweep(x, 2, rm)
    retval$colSDs <- sx <- apply(x, 2, sd, na.rm = na.rm)
    x <- sweep(x, 2, sx, "/")
  }
  if (missing(breaks) || is.null(breaks) || length(breaks) < 1) {
    if (missing(col) || is.function(col))
      breaks <- 16
    else breaks <- length(col) + 1
  }
  if (length(breaks) == 1) {
    if (!symbreaks)
      breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm),
                    length = breaks)
    else {
      extreme <- max(abs(x), na.rm = TRUE)
      breaks <- seq(-extreme, extreme, length = breaks)
    }
  }
  nbr <- length(breaks)
  ncol <- length(breaks) - 1
  if (class(col) == "function")
    col <- col(ncol)
  min.breaks <- min(breaks)
  max.breaks <- max(breaks)
  x[x < min.breaks] <- min.breaks
  x[x > max.breaks] <- max.breaks
  if (missing(lhei) || is.null(lhei))
    lhei <- c(keysize, 4)
  if (missing(lwid) || is.null(lwid))
    lwid <- c(keysize, 4)
  if (missing(lmat) || is.null(lmat)) {
    lmat <- rbind(4:3, 2:1)
    
    if (!missing(ColSideColors)) {
      #if (!is.matrix(ColSideColors))
      #stop("'ColSideColors' must be a matrix")
      if (!is.character(ColSideColors) || nrow(ColSideColors) != nc)
        stop("'ColSideColors' must be a matrix of nrow(x) rows")
      lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
      #lhei <- c(lhei[1], 0.2, lhei[2])
      lhei=c(lhei[1], side.height.fraction*ColSideColorsSize/2, lhei[2])
    }
    
    if (!missing(RowSideColors)) {
      #if (!is.matrix(RowSideColors))
      #stop("'RowSideColors' must be a matrix")
      if (!is.character(RowSideColors) || ncol(RowSideColors) != nr)
        stop("'RowSideColors' must be a matrix of ncol(x) columns")
      lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1), 1), lmat[,2] + 1)
      #lwid <- c(lwid[1], 0.2, lwid[2])
      lwid <- c(lwid[1], side.height.fraction*RowSideColorsSize/2, lwid[2])
    }
    lmat[is.na(lmat)] <- 0
  }
  
  if (length(lhei) != nrow(lmat))
    stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
  if (length(lwid) != ncol(lmat))
    stop("lwid must have length = ncol(lmat) =", ncol(lmat))
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  
  layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
  
  if (!missing(RowSideColors)) {
    if (!is.matrix(RowSideColors)){
      par(mar = c(margins[1], 0, 0, 0.5))
      image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
    } else {
      par(mar = c(margins[1], 0, 0, 0.5))
      rsc = t(RowSideColors[,rowInd, drop=F])
      rsc.colors = matrix()
      rsc.names = names(table(rsc))
      rsc.i = 1
      for (rsc.name in rsc.names) {
        rsc.colors[rsc.i] = rsc.name
        rsc[rsc == rsc.name] = rsc.i
        rsc.i = rsc.i + 1
      }
      rsc = matrix(as.numeric(rsc), nrow = dim(rsc)[1])
      image(t(rsc), col = as.vector(rsc.colors), axes = FALSE)
      if (length(rownames(RowSideColors)) > 0) {
        axis(1, 0:(dim(rsc)[2] - 1)/max(1,(dim(rsc)[2] - 1)), rownames(RowSideColors), las = 2, tick = FALSE)
      }
    }
  }
  
  if (!missing(ColSideColors)) {
    
    if (!is.matrix(ColSideColors)){
      par(mar = c(0.5, 0, 0, margins[2]))
      image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
    } else {
      par(mar = c(0.5, 0, 0, margins[2]))
      csc = ColSideColors[colInd, , drop=F]
      csc.colors = matrix()
      csc.names = names(table(csc))
      csc.i = 1
      for (csc.name in csc.names) {
        csc.colors[csc.i] = csc.name
        csc[csc == csc.name] = csc.i
        csc.i = csc.i + 1
      }
      csc = matrix(as.numeric(csc), nrow = dim(csc)[1])
      image(csc, col = as.vector(csc.colors), axes = FALSE)
      if (length(colnames(ColSideColors)) > 0) {
        axis(2, 0:(dim(csc)[2] - 1)/max(1,(dim(csc)[2] - 1)), colnames(ColSideColors), las = 2, tick = FALSE)
      }
    }
  }
  
  par(mar = c(margins[1], 0, 0, margins[2]))
  x <- t(x)
  cellnote <- t(cellnote)
  if (revC) {
    iy <- nr:1
    if (exists("ddr"))
      ddr <- rev(ddr)
    x <- x[, iy]
    cellnote <- cellnote[, iy]
  }
  else iy <- 1:nr
  image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col, breaks = breaks, ...)
  retval$carpet <- x
  if (exists("ddr"))
    retval$rowDendrogram <- ddr
  if (exists("ddc"))
    retval$colDendrogram <- ddc
  retval$breaks <- breaks
  retval$col <- col
  if (!invalid(na.color) & any(is.na(x))) { # load library(gplots)
    mmat <- ifelse(is.na(x), 1, NA)
    image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "",
          col = na.color, add = TRUE)
  }
  axis(1, 1:nc, labels = labCol, las = 2, line = -0.5, tick = 0,
       cex.axis = cexCol)
  if (!is.null(xlab))
    mtext(xlab, side = 1, line = margins[1] - 1.25)
  axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0,
       cex.axis = cexRow)
  if (!is.null(ylab))
    mtext(ylab, side = 4, line = margins[2] - 1.25)
  if (!missing(add.expr))
    eval(substitute(add.expr))
  if (!missing(colsep))
    for (csep in colsep) rect(xleft = csep + 0.5, ybottom = rep(0, length(csep)), xright = csep + 0.5 + sepwidth[1], ytop = rep(ncol(x) + 1, csep), lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
  if (!missing(rowsep))
    for (rsep in rowsep) rect(xleft = 0, ybottom = (ncol(x) + 1 - rsep) - 0.5, xright = nrow(x) + 1, ytop = (ncol(x) + 1 - rsep) - 0.5 - sepwidth[2], lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
  min.scale <- min(breaks)
  max.scale <- max(breaks)
  x.scaled <- scale01(t(x), min.scale, max.scale)
  if (trace %in% c("both", "column")) {
    retval$vline <- vline
    vline.vals <- scale01(vline, min.scale, max.scale)
    for (i in colInd) {
      if (!is.null(vline)) {
        abline(v = i - 0.5 + vline.vals, col = linecol,
               lty = 2)
      }
      xv <- rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
      xv <- c(xv[1], xv)
      yv <- 1:length(xv) - 0.5
      lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
    }
  }
  if (trace %in% c("both", "row")) {
    retval$hline <- hline
    hline.vals <- scale01(hline, min.scale, max.scale)
    for (i in rowInd) {
      if (!is.null(hline)) {
        abline(h = i + hline, col = linecol, lty = 2)
      }
      yv <- rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
      yv <- rev(c(yv[1], yv))
      xv <- length(yv):1 - 0.5
      lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
    }
  }
  if (!missing(cellnote))
    text(x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote),
         col = notecol, cex = notecex)
  par(mar = c(margins[1], 0, 0, 0))
  if (dendrogram %in% c("both", "row")) {
    plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
  }
  else plot.new()
  par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
  if (dendrogram %in% c("both", "column")) {
    plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
  }
  else plot.new()
  if (!is.null(main))
    title(main, cex.main = 1.5 * op[["cex.main"]])
  if (key) {
    par(mar = c(5, 4, 2, 1), cex = 0.75)
    tmpbreaks <- breaks
    if (symkey) {
      max.raw <- max(abs(c(x, breaks)), na.rm = TRUE)
      min.raw <- -max.raw
      tmpbreaks[1] <- -max(abs(x), na.rm = TRUE)
      tmpbreaks[length(tmpbreaks)] <- max(abs(x), na.rm = TRUE)
    }
    else {
      min.raw <- min(x, na.rm = TRUE)
      max.raw <- max(x, na.rm = TRUE)
    }
    
    z <- seq(min.raw, max.raw, length = length(col))
    image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks,
          xaxt = "n", yaxt = "n")
    par(usr = c(0, 1, 0, 1))
    lv <- pretty(breaks)
    xv <- scale01(as.numeric(lv), min.raw, max.raw)
    axis(1, at = xv, labels = lv)
    if (scale == "row")
      mtext(side = 1, "Row Z-Score", line = 2)
    else if (scale == "column")
      mtext(side = 1, "Column Z-Score", line = 2)
    else mtext(side = 1, KeyValueName, line = 2)
    if (density.info == "density") {
      dens <- density(x, adjust = densadj, na.rm = TRUE)
      omit <- dens$x < min(breaks) | dens$x > max(breaks)
      dens$x <- dens$x[-omit]
      dens$y <- dens$y[-omit]
      dens$x <- scale01(dens$x, min.raw, max.raw)
      lines(dens$x, dens$y/max(dens$y) * 0.95, col = denscol,
            lwd = 1)
      axis(2, at = pretty(dens$y)/max(dens$y) * 0.95, pretty(dens$y))
      title("Color Key\nand Density Plot")
      par(cex = 0.5)
      mtext(side = 2, "Density", line = 2)
    }
    else if (density.info == "histogram") {
      h <- hist(x, plot = FALSE, breaks = breaks)
      hx <- scale01(breaks, min.raw, max.raw)
      hy <- c(h$counts, h$counts[length(h$counts)])
      lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s",
            col = denscol)
      axis(2, at = pretty(hy)/max(hy) * 0.95, pretty(hy))
      title("Color Key\nand Histogram")
      par(cex = 0.5)
      mtext(side = 2, "Count", line = 2)
    }
    else title("Color Key")
  }
  else plot.new()
  retval$colorTable <- data.frame(low = retval$breaks[-length(retval$breaks)],
                                  high = retval$breaks[-1], color = retval$col)
  invisible(retval)
}

###########
#PAM functions
############


PAM_Analysis_MaxNrFeatures_CV_KB <- function (MA_Data,CL_Data,NrFolds=10,NrLambdaFolds=10,MaxNrFeatures) {
  
  # overlapping samples
  OverlapSamples=intersect(colnames(MA_Data),colnames(CL_Data))
  if (length(OverlapSamples)>0) {
    MA_Data=MA_Data[,OverlapSamples]
    CL_Data=CL_Data[,OverlapSamples,drop=FALSE]
    cat('There are',length(OverlapSamples),'overlapping samples between predictor data and outcome data.\n')     
  } else {
    cat('No samples are overlapping, returning.')     
    return()          
  }
  
  # creating a frame structure for calculating the folds. 
  CL_Data_frame=data.frame(t(CL_Data))     
  names(CL_Data_frame)[1]="outcome"
  
  # Check for Leave-One-Out CV
  if (NrFolds >ncol(CL_Data)) { 
    cat("There are not enough samples ???.\n")
    return()
  }
  if (NrFolds==ncol(CL_Data)) {
    cat('Doing leave-one-out cross validation.\n')
    CL_Data_frame$folds=1:ncol(CL_Data)
  } else {
    cat('Doing stratified nfold cross validation.\n')          
    # using a split
    CL_Data_frame <- ddply(CL_Data_frame,.(outcome),createFolds,k = NrFolds)
    rownames(CL_Data_frame)=colnames(CL_Data)
  }
  
  Predictions=matrix(0,ncol(CL_Data),1)
  rownames(Predictions)=colnames(CL_Data)
  colnames(Predictions)=rownames(CL_Data)
  NrGenesPerFold=matrix(0,NrFolds,1)
  Genes=list()
  GenesRanks=list()
  
  PosteriorProbabilities=matrix(0,ncol(CL_Data),length(unique(CL_Data[1,])))
  rownames(PosteriorProbabilities)=colnames(CL_Data)
  #colnames(PosteriorProbabilities)=rownames(CL_Data)
  
  for (i in 1:NrFolds) {
    TestSamples=rownames(CL_Data_frame)[CL_Data_frame$folds==i]
    TrainSamples=setdiff(rownames(CL_Data_frame),TestSamples)
    PAMresults=PAM_Analysis_MaxNrFeatures_KB(MA_Data[,TrainSamples],CL_Data[,TrainSamples,drop=FALSE],MA_Data[,TestSamples,drop=FALSE],NrLambdaFolds,MaxNrFeatures, HeteroNorm=NULL)          
    Predictions[TestSamples,1]=PAMresults$TestPredictions
    PosteriorProbabilities[TestSamples,1:ncol(PosteriorProbabilities)]=PAMresults$PosteriorProbabilites
    NrGenesPerFold[i]=PAMresults$OptimalSolution$minGenes
    #Need to work on getting the genes from this, make object to 
    Genes[[i]]=PAMresults$Nonzero
    GenesRanks[[i]]=PAMresults$PAM.Gene.Ranks
  }    
  #is there gene names in here???
  # Calculate the performance
  NrClasses=length(unique(CL_Data[1,]))     
  
  NrMisclassificationsPerClass=matrix(0,NrClasses,4)
  colnames(NrMisclassificationsPerClass)=c('NrMistakes','TotalNrCasesInClass','NrMisclassificationsPerc','Accuracy')
  tmpData=split(Predictions,CL_Data)
  ClassNames=c()
  for (i in 1:NrClasses) {
    NrMisclassificationsPerClass[i,1]=sum(tmpData[[i]]!=i)
    NrMisclassificationsPerClass[i,2]=length(tmpData[[i]])
    NrMisclassificationsPerClass[i,3]=NrMisclassificationsPerClass[i,1]/NrMisclassificationsPerClass[i,2]
    NrMisclassificationsPerClass[i,4]=1-NrMisclassificationsPerClass[i,3]
    ClassNames=c(ClassNames,paste('Class',i,sep=''))
  }
  rownames(NrMisclassificationsPerClass)=ClassNames
  
  OverallPerformance=matrix(0,1,3)
  colnames(OverallPerformance)=c('NrMisclassifications','NrMisclassificationsPerc','Accuracy')
  OverallPerformance[1]=sum(Predictions-t(CL_Data) !=0)
  OverallPerformance[2]=OverallPerformance[1]/ncol(CL_Data)
  OverallPerformance[3]=1-OverallPerformance[1]/ncol(CL_Data)
  
  return(list(Predictions=Predictions,PosteriorProbabilities=PosteriorProbabilities,OverallPerformance=OverallPerformance,NrMisclassificationsPerClass=NrMisclassificationsPerClass,NonZero.Genes=Genes,GenesRanks=GenesRanks, NrGenesPerFold=NrGenesPerFold,CL_Data_Overlap=CL_Data))
}

createFolds <- function(x,k){
  n <- nrow(x)
  x$folds <- rep(1:k,length.out = n)[sample(n,n)]
  x
}


PAM_Analysis_MaxNrFeatures_KB <- function (MA_Data_Train,CL_Data_Train,MA_Data_Test,NrLambdaFolds=5,MaxNrFeatures, HeteroNorm=NULL) {
  
  # Args[8]="/Users/ogevaert/Documents/WORK/Projects/Stanford/TCGAglioblastoma/data/2014-06-13/OutcomeDataTestPredictions_1.txt"
  # Args[9]="/Users/ogevaert/Documents/WORK/Projects/Stanford/TCGAglioblastoma/results/2014-06-13/CentroidFigure_1.tif"
  # Args[10]="/Users/ogevaert/Documents/WORK/Projects/Stanford/TCGAglioblastoma/results/2014-06-13/CentroidData_1.txt"
  
  # overlapping samples
  OverlapSamples=intersect(colnames(MA_Data_Train),colnames(CL_Data_Train))
  if (length(OverlapSamples)>0) {
    MA_Data_Train=MA_Data_Train[,OverlapSamples]
    CL_Data_Train=CL_Data_Train[,OverlapSamples,drop=FALSE]
    cat('There are',length(OverlapSamples),'overlapping samples between predictor data and outcome data.\n')     
  } else {
    cat('No samples are overlapping, returning.')     
    return()          
  }
  
  dim(MA_Data_Train)
  dim(CL_Data_Train)
  ModuleDataDimensions=dim(MA_Data_Train)
  
  dim(MA_Data_Test)
  ModuleDataDimensionsTest=dim(MA_Data_Test)
  
  # Creating the data object
  PAM.TrainData=list(x=data.matrix(MA_Data_Train),y=as.numeric(data.matrix(CL_Data_Train)),genenames=as.character(row.names(MA_Data_Train)),geneid=as.character(row.names(MA_Data_Train)))
  
  # Training the PAM model
  if(missing(HeteroNorm)){
    PAM.train <- pamr.train(PAM.TrainData)
  }else{
    PAM.train <- pamr.train(PAM.TrainData, hetero=HeteroNorm)
  }
  
  # cross validation
  PAM.cvresults<-pamr.cv(PAM.train,PAM.TrainData,NrLambdaFolds)
  
  #280,801
  
  # plotting the cross validation error
  #pamr.plotcv(PAM.cvresults)
  
  # confusion plot at lowest error
  #PAM.cvresults$error
  #PAM.cvresults$threshold
  
  # choosing threshold when balanced error rate is used
  BERresults=CalculateBalancedErrorRate(PAM.cvresults)
  BER=BERresults$BER
  #BER=1-BERresults$Specificity     
  
  minOneVariablePositions=PAM.cvresults$size>=2 & PAM.cvresults$size<=MaxNrFeatures # at least two genes in the PAM model
  minError=min(BER[minOneVariablePositions])
  tmpPos=which(BER==minError)
  minGenes=min(PAM.cvresults$size[tmpPos])
  Thresholds=PAM.cvresults$threshold[tmpPos] 
  #Threshold= Thresholds[1]
  #Changed treshold to one with min genes; previous version arbitrarily chose first thresholds in set of tresholds with min genes, which is always the lowest of 
  #those tresholds with less than MaxNrFeatures. So a less stringent treshold was actually being applied to predict the treshold than was suggested by MaxNrFeatures
  Threshold=PAM.cvresults$threshold[tmpPos][which.min(PAM.cvresults$size[tmpPos])]
  OptimalSolution=list(minGenes=minGenes,minError=minError,Threshold=Threshold)
  #PAM.trainOptimal <- pamr.train(PAM.TrainData,threshold=OptimalSolution$Threshold)
  
  ## Compute the confusion matrix for a particular model
  #pamr.confusion(PAM.cvresults, Threshold)
  
  ## Plot the cross-validated class probabilities by class
  #pamr.plotcvprob(PAM.cvresults, PAM.TrainData, Threshold)
  
  ## Plot the class centroids
  #pamr.plotcen(PAM.train, PAM.TrainData, Threshold)
  
  ## Make a gene plot of the most significant genes (less usefull)
  #pamr.geneplot(PAM.train, PAM.TrainData, Threshold)
  
  # Estimate false discovery rates and plot them
  #fdr.obj<- pamr.fdr(PAM.train, PAM.TrainData)
  
  #pamr.plotfdr(fdr.obj)
  
  #############################################################################################################
  ##### PAM predicting test set
  ##############################################################################################################
  
  # making the test data set
  PAM.TestData=list(x=data.matrix(MA_Data_Test),genenames=row.names(MA_Data_Test))
  
  # predicting the cohort 2 samples
  PAM.TestPredictions=pamr.predict(PAM.train,PAM.TestData$x,Threshold,type="class")
  PAM.PosteriorProbabilities=pamr.predict(PAM.train,PAM.TestData$x,Threshold,type="posterior")
  PAM.Centroid=pamr.predict(PAM.train,PAM.TestData$x,Threshold,type="centroid") # this gives the unshrunken centroids
  PAM.Nonzero=PAM.TestData$genenames[pamr.predict(PAM.train,PAM.TestData$x,Threshold,type="nonzero")]
  
  #get the ranking of the genes across earch round of cross validation# Nonzero genes (those that passed the threshold) should include the top ranking genes
  PAM.g=pamr.listgenes(PAM.train, PAM.TrainData,  Threshold, PAM.cvresults, genenames=FALSE)
  rownames(PAM.g)=PAM.g[,1]
  
  # Centroid outputs the unshrunken centroids, PAM.Nonzero has the genes, so save the shrunken ones here. 
  PAM.Centroid=PAM.Centroid[PAM.Nonzero,]
  
  PAM.TestPredictions=matrix(as.numeric(PAM.TestPredictions))
  rownames(PAM.TestPredictions)=colnames(MA_Data_Test)
  colnames(PAM.TestPredictions)="PAMprediction"
  
  #PAM.TestData$genenames[PAM.Nonzero]
  #PAM.AllTestPredictions<-array(dim=c(length(AllThresholds),1))
  #AllThresholds=PAM.cvresults$threshold
  #for (i in c(1:length(AllThresholds))) {
  #     results=pamr.predict(PAM.train,PAM.TestData$x, AllThresholds[i])
  #}
  
  #      # writing this to file
  
  
  #      write.table(AllResults,Args[8],sep="\t",row.names=TRUE)
  #      
  #      # writing also the PosteriorProbabilities to file
  #      Filename=paste(substr(Args[10], 1, nchar(Args[10])-4),"_PosteriorProbabilities.txt",sep='')
  #      write.table(PAM.PosteriorProbabilities,Filename,sep="\t",row.names=TRUE)
  
  return(list(TestPredictions=PAM.TestPredictions,PosteriorProbabilites=PAM.PosteriorProbabilities,Centroids=PAM.Centroid,Nonzero=PAM.Nonzero,PAM.Gene.Ranks=PAM.g,OptimalSolution=OptimalSolution))
  
}

CalculateBalancedErrorRate <- function(pamfit) {          
  BER=matrix(1,length(pamfit$threshold),1)
  Specificity=matrix(0,length(pamfit$threshold),1)
  Sensitivity=matrix(0,length(pamfit$threshold),1)
  for (i in 1:length(pamfit$threshold)) {
    tmpConfusionTable=pamr.confusion(pamfit,pamfit$threshold[i],extra=FALSE)
    Sensitivity[i,1]=tmpConfusionTable[2,2]/(sum(tmpConfusionTable[2,]))
    Specificity[i,1]=tmpConfusionTable[1,1]/(sum(tmpConfusionTable[1,]))
    
    ErrorsPerClass=matrix(0,nrow(tmpConfusionTable),1)
    for (j in 1:nrow(tmpConfusionTable)) {
      ErrorsPerClass[j,]=1-tmpConfusionTable[j,j]/sum(tmpConfusionTable[j,])
    }
    BER[i,]=mean(ErrorsPerClass)
  } 
  return(list(BER=BER,Sensitivity=Sensitivity,Specificity=Specificity))
}

##############
#SAM analysis
###################

SAM_Analysis_meth <-function(MA_Data,OutcomeData,OutcomeType,DeltaFDRthreshold,nrPermutations,StatisticalTest,OutputResultsFile,OutputFigureFile) {
  
  # overlapping samples
  OverlapSamples=intersect(colnames(MA_Data),colnames(OutcomeData))
  if (length(OverlapSamples)>0) {
    MA_Data=MA_Data[,OverlapSamples]
    OutcomeData=OutcomeData[,OverlapSamples,drop=FALSE]
    cat('There are',length(OverlapSamples),'overlapping samples between predictor data and outcome data.')     
  } else {
    cat('No samples are overlapping, returning.')     
    return()          
  }
  
  # Checking data dimensions
  dim(MA_Data)
  MA_DataDimensions=dim(MA_Data)
  dim(OutcomeData)
  OutcomeVariables=row.names(OutcomeData)
  
  if (ncol(MA_Data) != ncol(OutcomeData)) {
    cat('The number of samples in the data and the outcome file is not a match')     
    return()
  }
  
  # threshold to determine Delta
  #DeltaFDRthreshold=0.05
  
  deltaDefault=0.1
  AllResults=matrix(,nrow=0,ncol=4)
  if (OutcomeType=='Quantitative' ||  OutcomeType=='Multiclass') {
    AllResults=matrix(,nrow=0,ncol=4)
  }
  if (OutcomeType=='TwoClassUnpaired' ) {
    OutcomeType='Two class unpaired'
    AllResults=matrix(,nrow=0,ncol=4)
  }
  
  AllDeltaTables=list()
  for(i in 1:nrow(OutcomeData)) {
    # making the data structure
    data=list(x=MA_Data,y=OutcomeData[i,],geneid=row.names(MA_Data),genenames=row.names(MA_Data), logged2=TRUE, testStatistic=StatisticalTest)
    
    OutcomeVariables[i]	
    if (OutcomeType=='Quantitative' ) {
      data$x=matrix(data.matrix(data$x),nrow=MA_DataDimensions[1],ncol=MA_DataDimensions[2])
      data$y=as.numeric(data$y)
    }
    
    # sum((data$y==1))>1 en sum((data$y==2))>1) zijn voor binaire variabelen zodat er minstens 2 enen of 2 tweeen zijn
    # || OutcomeType=='Quantitative' is zodat als het quantitative output is, er met het vorige geen rekening wordt gehouden
    # sum(data$y)!=0 is voor quantitatieve varabielen die overal 0 zijn. 
    
    if (((sum((data$y==1))>1 && sum((data$y==2))>1) || OutcomeType=='Quantitative' || OutcomeType=='Multiclass') && sum(data$y)!=0 || OutcomeType=='Two class paired') {
      
      # Actual SAM analysis
      samr.obj=samr(data,resp.type=OutcomeType,nperms=nrPermutations, regression.method="ranks")
      samr.obj$foldchange=apply(as.data.frame(data$x),1,function(x) mean(x[data$y==2])-mean(x[data$y==1]))
      
      #need to change the foldchange to differential methylatio, 
      
      # Choosing the delta
      delta.table=samr.compute.delta.table(samr.obj)
      tmpPositions=which(delta.table[,5]<DeltaFDRthreshold)
      delta=delta.table[tmpPositions[1],1]
      
      # storing delta tables
      AllDeltaTables[[i]]=delta.table
      
      if (!is.na(delta)) {
        siggenes.table=samr.compute.siggenes.table(samr.obj,delta,data,delta.table)
      } else {
        siggenes.table=samr.compute.siggenes.table(samr.obj,deltaDefault,data,delta.table)
        deltaDefault
        deltaChoice='No results for computed delta, using default delta.'
        deltaChoice
        cat('No results for computed delta, using default delta. This should not happen. Try a more loose DeltaFDRthreshold.')
      }
      
      # getting this in a matrix format:
      UpMetaGenes=siggenes.table$genes.up
      DownMetaGenes=siggenes.table$genes.lo
      
      if (length(UpMetaGenes)==0 &&  length(DownMetaGenes)==0 ) {
        cat('There are no up and down regulated genes.\n')
      }
      
      if (OutcomeType=='Quantitative' ) {
        SamResults=rbind(UpMetaGenes[,c(3,7)],DownMetaGenes[,c(3,7)])
      }
      if (OutcomeType=='Two class unpaired' || OutcomeType=='Multiclass' || OutcomeType=='Two class paired') {
        SamResults=rbind(UpMetaGenes[,c(3,7,8)],DownMetaGenes[,c(3,7,8)])
      }
      if (!is.null(SamResults)) {
        Variable=rep(OutcomeVariables[i],times=dim(SamResults)[1])
        CurrentResult=cbind(Variable,SamResults)
        AllResults=rbind(AllResults,CurrentResult)
      }
    }	
  }
   ##### Writing all results to file
   if (OutputResultsFile != '') {
    write.table(AllResults,OutputResultsFile,sep="\t",row.names=FALSE)
  }     
  
  # plotting
  if (OutputFigureFile != '') {
    tiff(filename=OutputFigureFile,compression="none",pointsize=20,bg="white",width=1000,height=1000)
    samr.plot(samr.obj,delta)
    dev.off()
  }
  return(list(AllResults=AllResults,AllDeltaTables=AllDeltaTables))
}

###############
#rowwise wilcox rank sum test
#################################

#Need to change loops to apply, as this function takes ages to run
Rowwise_stat=function(dataframe=dataframe, var=var, selectfeatures=FALSE, mindiff=FALSE){
     if(!is.na(selectfeatures)){
          
          dat=array(NA,c(nrow(dataframe),2))
          rownames(dat)=rownames(dataframe)
          colnames(dat)=c("Diff","P.value")
          
          for( i in 1:nrow(dataframe)){
               dat[i,1]=mean(as.numeric(as.character(dataframe[i,which(var==levels(var)[2])])), na.rm=T)-mean(as.numeric(as.character(dataframe[i,which(var==levels(var)[1])])),  na.rm=T)
               dat[i,2]=wilcox.test(as.numeric(as.character(dataframe[i,]))~var)$p.value
          }
     } else {
          dat=array(NA,c(nrow(dataframe),2))
          rownames(dat)=rownames(dataframe)
          colnames(dat)=c("Diff","P.value")
          
          for( i in 1:nrow(dataframe)){
               dat[i,1]=mean(as.numeric(as.character(dataframe[i,which(var==levels(var)[2])])), na.rm=T)-mean(as.numeric(as.character(dataframe[i,which(var==levels(var)[1])])),  na.rm=T)
               dat[i,2]=wilcox.test(as.numeric(as.character(dataframe[i,]))~var)$p.value
          }
     }
     
     dat=as.data.frame(dat, drop=FALSE)
     dat$Q.value=p.adjust(dat$P.value, method="fdr")
     dat=dat[with(dat, order(abs(dat$Diff), decreasing=T)),]
     return(dat)
}

#
Plot_Beta_Univ_JustHist <-function(GeneName, MixtureModelResults, METdata, METnormal=0, FileName="") {
     
     Pos = which(rownames(MixtureModelResults$MethylationStates) %in% GeneName)
     Pos2 = which(rownames(METdata) %in% GeneName)
     if (length(Pos) > 0 & length(Pos2) > 0) {
          BetaModel = MixtureModelResults$Models[[Pos]]
          NrComponents = MixtureModelResults$NrComponents[Pos]
          METdataVector=METdata[Pos2,]
          
          if (FileName != "") {
               File = paste(FileName, "MixtureModel.tif", sep = '')
               tiff(filename = File, compression = "none", pointsize = 20, bg = "white", width = 600, height = 600)                     
          }
          histogram=hist(METdataVector,50,plot=FALSE)
          HistHeight=max(histogram$counts)*1.2
          histogram=hist(METdataVector,50,xlim=c(0,1),ylim=c(0,1.4*HistHeight),main=paste("Mixture model of",GeneName),xlab="DNA methylation",ylab="Frequency")
          
          # plotting the normal 95% confidence interval, if the normal data is present
          if (length(METnormal) > 1) {
               tmpTtest=t.test(METnormal[GeneName,])            
               NormalMean=mean(METnormal[GeneName,])
               lines(c(tmpTtest$conf.int[1],tmpTtest$conf.int[2]),c(1.2*HistHeight,1.2*HistHeight),ylim=c(0,1.2*HistHeight),lwd=5,type='l',col=1)
               points(NormalMean,1.2*HistHeight,ylim=c(0,1.2*HistHeight),lwd=5,type='l',col=1,pch=8)
          }   
          
          x = seq(0,1,0.001)
          MaxHeight = numeric(NrComponents)
          for (comp in 1:NrComponents) {
               MaxHeight[comp] = max(dbeta(x, BetaModel$a[comp], BetaModel$b[comp]))
          }
          MaxHeight = max(MaxHeight)
          Factor = max(HistHeight/MaxHeight)
          for (comp in 1:NrComponents) {
               lines(x, Factor * BetaModel$eta[comp] * dbeta(x, BetaModel$a[comp], BetaModel$b[comp]), ylim = c(0, 1.2 * HistHeight), lwd = 5, type = 'l', col = comp + 1, xlab = "DNA methylation", ylab = "Frequency")  
          }
     }
}


Plot_Beta_Univ_JustHist_with_oncocytoma <-function(GeneName, MixtureModelResults, METdata, METnormal=0, FileName="") {
     
     Pos = which(rownames(MixtureModelResults$MethylationStates) %in% GeneName)
     Pos2 = which(rownames(METdata) %in% GeneName)
     if (length(Pos) > 0 & length(Pos2) > 0) {
          BetaModel = MixtureModelResults$Models[[Pos]]
          NrComponents = MixtureModelResults$NrComponents[Pos]
          METdataVector=METdata[Pos2,]
          
          if (FileName != "") {
               File = paste(FileName, "MixtureModel.tif", sep = '')
               tiff(filename = File, compression = "none", pointsize = 20, bg = "white", width = 600, height = 600)                     
          }
          histogram=hist(METdataVector,50,plot=FALSE)
          HistHeight=max(histogram$counts)*1.2
          histogram=hist(METdataVector,50,xlim=c(0,1),ylim=c(0,1.4*HistHeight),main=paste("Mixture model of",GeneName),xlab="DNA methylation",ylab="Frequency")
          
          # plotting the normal 95% confidence interval, if the normal data is present
          if (length(METnormal) > 1) {
               tmpTtest=t.test(METnormal[GeneName,])            
               NormalMean=mean(METnormal[GeneName,])
               lines(c(tmpTtest$conf.int[1],tmpTtest$conf.int[2]),c(1.2*HistHeight,1.2*HistHeight),ylim=c(0,1.2*HistHeight),lwd=5,type='l',col=1)
               points(NormalMean,1.2*HistHeight,ylim=c(0,1.2*HistHeight),lwd=5,type='l',col=1,pch=8)
          }  
          
          # plotting the oncocytoma 95% confidence interval, if the oncocytoma data is present
          
          
          
          
          
          x = seq(0,1,0.001)
          MaxHeight = numeric(NrComponents)
          for (comp in 1:NrComponents) {
               MaxHeight[comp] = max(dbeta(x, BetaModel$a[comp], BetaModel$b[comp]))
          }
          MaxHeight = max(MaxHeight)
          Factor = max(HistHeight/MaxHeight)
          for (comp in 1:NrComponents) {
               lines(x, Factor * BetaModel$eta[comp] * dbeta(x, BetaModel$a[comp], BetaModel$b[comp]), ylim = c(0, 1.2 * HistHeight), lwd = 5, type = 'l', col = comp + 1, xlab = "DNA methylation", ylab = "Frequency")  
          }
     }
}


PAM_Analysis_CV_KB <- function (MA_Data,CL_Data,NrFolds=5,NrLambdaFolds=5, HeteroNorm=NULL) {
     #NrFolds tells the algorithm how many rounds of cross-validation to do
     # overlapping samples
     Genes=list()
     OverlapSamples=intersect(colnames(MA_Data),colnames(CL_Data))
     if (length(OverlapSamples)>0) {
          MA_Data=MA_Data[,OverlapSamples]
          CL_Data=CL_Data[,OverlapSamples,drop=FALSE]
          #drop=FALSE prevents the dataframe or matrix from being cooerced to a vector if there's only one column left
          cat('There are',length(OverlapSamples),'overlapping samples between predictor data and outcome data.\n')     
     } else {
          cat('No samples are overlapping, returning.')     
          return()          
     }
     
     # creating a frame structure for calculating the folds. 
     CL_Data_frame=data.frame(t(CL_Data))     
     names(CL_Data_frame)[1]="outcome"
     
     # Check for Leave-One-Out CV
     if (NrFolds >ncol(CL_Data)) { 
          cat("There are not enough samples ???.\n")
          return()
     }
     if (NrFolds==ncol(CL_Data)) {
          cat('Doing leave-one-out cross validation.\n')
          CL_Data_frame$folds=1:ncol(CL_Data)
     } else {
          cat('Doing stratified nfold cross validation.\n')          
          # using a split
          CL_Data_frame <- ddply(CL_Data_frame,.(outcome),createFolds,k = NrFolds)
          #ddply For each subset of a data frame, apply function then combine results into a data frame, like tapply
          #ddply(data,.outcome,.function)
          rownames(CL_Data_frame)=colnames(CL_Data)
     }
     
     Predictions=matrix(0,ncol(CL_Data),1)
     NrGenesPerFold=matrix(0,NrFolds,1)
     rownames(Predictions)=colnames(CL_Data)
     colnames(Predictions)=rownames(CL_Data)
     
     PosteriorProbabilities=matrix(0,ncol(CL_Data),length(unique(CL_Data[1,])))
     rownames(PosteriorProbabilities)=colnames(CL_Data)
     #colnames(PosteriorProbabilities)=rownames(CL_Data)
     for (i in 1:NrFolds) {
          TestSamples=rownames(CL_Data_frame)[CL_Data_frame$folds==i]
          TrainSamples=setdiff(rownames(CL_Data_frame),TestSamples)
          
          if(missing(HeteroNorm)){
               PAMresults=PAM_Analysis_KB(MA_Data[,TrainSamples],CL_Data[,TrainSamples,drop=FALSE],MA_Data[,TestSamples,drop=FALSE],NrLambdaFolds)          
          } else {
               PAMresults=PAM_Analysis_KB(MA_Data[,TrainSamples],CL_Data[,TrainSamples,drop=FALSE],MA_Data[,TestSamples,drop=FALSE],NrLambdaFolds, HeteroNorm=HeteroNorm)          
          }
          
          Predictions[TestSamples,1]=PAMresults$TestPredictions
          PosteriorProbabilities[TestSamples,1:ncol(PosteriorProbabilities)]=PAMresults$PosteriorProbabilites
          NrGenesPerFold[i]=PAMresults$OptimalSolution$minGenes
          Genes[[i]]=PAMresults$Nonzero   
     }    
     
     # Calculate the performance
     NrClasses=length(unique(CL_Data[1,]))     
     
     NrMisclassificationsPerClass=matrix(0,NrClasses,4)
     colnames(NrMisclassificationsPerClass)=c('NrMistakes','TotalNrCasesInClass','NrMisclassificationsPerc','Accuracy')
     tmpData=split(Predictions,CL_Data)
     ClassNames=c()
     for (i in 1:NrClasses) {
          NrMisclassificationsPerClass[i,1]=sum(tmpData[[i]]!=i)
          NrMisclassificationsPerClass[i,2]=length(tmpData[[i]])
          NrMisclassificationsPerClass[i,3]=NrMisclassificationsPerClass[i,1]/NrMisclassificationsPerClass[i,2]
          NrMisclassificationsPerClass[i,4]=1-NrMisclassificationsPerClass[i,3]
          ClassNames=c(ClassNames,paste('Class',i,sep=''))
     }
     rownames(NrMisclassificationsPerClass)=ClassNames
     
     OverallPerformance=matrix(0,1,3)
     colnames(OverallPerformance)=c('NrMisclassifications','NrMisclassificationsPerc','Accuracy')
     OverallPerformance[1]=sum(Predictions-t(CL_Data) !=0)
     OverallPerformance[2]=OverallPerformance[1]/ncol(CL_Data)
     OverallPerformance[3]=1-OverallPerformance[1]/ncol(CL_Data)
     
     return(list(Predictions=Predictions,PosteriorProbabilities=PosteriorProbabilities,OverallPerformance=OverallPerformance,NrMisclassificationsPerClass=NrMisclassificationsPerClass,NrGenesPerFold=NrGenesPerFold,CL_Data_Overlap=CL_Data, Genes=Genes))
}

#Made a version of PAM_Analysis that enables specification of hetero argument, i.e., allowing specification of a normal group
PAM_Analysis_KB <- function (MA_Data_Train,CL_Data_Train,MA_Data_Test,NrLambdaFolds=5, HeteroNorm=NULL) {
     
     # Args[8]="/Users/ogevaert/Documents/WORK/Projects/Stanford/TCGAglioblastoma/data/2014-06-13/OutcomeDataTestPredictions_1.txt"
     # Args[9]="/Users/ogevaert/Documents/WORK/Projects/Stanford/TCGAglioblastoma/results/2014-06-13/CentroidFigure_1.tif"
     # Args[10]="/Users/ogevaert/Documents/WORK/Projects/Stanford/TCGAglioblastoma/results/2014-06-13/CentroidData_1.txt"
     
     # overlapping samples
     OverlapSamples=intersect(colnames(MA_Data_Train),colnames(CL_Data_Train))
     if (length(OverlapSamples)>0) {
          MA_Data_Train=MA_Data_Train[,OverlapSamples]
          CL_Data_Train=CL_Data_Train[,OverlapSamples,drop=FALSE]
          cat('There are',length(OverlapSamples),'overlapping samples between predictor data and outcome data.\n')     
     } else {
          cat('No samples are overlapping, returning.')     
          return()          
     }
     
     dim(MA_Data_Train)
     dim(CL_Data_Train)
     ModuleDataDimensions=dim(MA_Data_Train)
     
     dim(MA_Data_Test)
     ModuleDataDimensionsTest=dim(MA_Data_Test)
     
     # Creating the data object
     PAM.TrainData=list(x=data.matrix(MA_Data_Train),y=as.numeric(data.matrix(CL_Data_Train)),genenames=as.character(row.names(MA_Data_Train)),geneid=as.character(row.names(MA_Data_Train)))
     
     # Training the PAM model
     if(missing(HeteroNorm)){
          PAM.train <- pamr.train(PAM.TrainData)
     }else{
          PAM.train <- pamr.train(PAM.TrainData, hetero=HeteroNorm)
     }
     
     # cross validation
     PAM.cvresults<-pamr.cv(PAM.train,PAM.TrainData,NrLambdaFolds)
     
     # plotting the cross validation error
     #pamr.plotcv(PAM.cvresults)
     
     # confusion plot at lowest error
     #PAM.cvresults$error
     #PAM.cvresults$threshold
     
     # choosing threshold when balanced error rate is used
     BERresults=CalculateBalancedErrorRate(PAM.cvresults)
     BER=BERresults$BER
     minOneVariablePositions=PAM.cvresults$size>=2 # at least two genes in the PAM model
     minError=min(BER[minOneVariablePositions])
     tmpPos=which(BER==minError)
     minGenes=min(PAM.cvresults$size[tmpPos])
     Thresholds=PAM.cvresults$threshold[tmpPos] 
     Threshold= Thresholds[1]
     minError
     minGenes
     Threshold
     OptimalSolution=list(minGenes=minGenes,minError=minError,Threshold=Threshold)
     
     ## Compute the confusion matrix for a particular model
     #pamr.confusion(PAM.cvresults, Threshold)
     
     ## Plot the cross-validated class probabilities by class
     #pamr.plotcvprob(PAM.cvresults, PAM.TrainData, Threshold)
     
     ## Plot the class centroids
     #pamr.plotcen(PAM.train, PAM.TrainData, Threshold)
     
     ## Make a gene plot of the most significant genes (less usefull)
     #pamr.geneplot(PAM.train, PAM.TrainData, Threshold)
     
     # Estimate false discovery rates and plot them
     #fdr.obj<- pamr.fdr(PAM.train, PAM.TrainData)
     
     #pamr.plotfdr(fdr.obj)
     
     #############################################################################################################
     ##### PAM predicting test set
     ##############################################################################################################
     
     # making the test data set
     PAM.TestData=list(x=data.matrix(MA_Data_Test),genenames=row.names(MA_Data_Test))
     
     # predicting the cohort 2 samples
     PAM.TestPredictions=pamr.predict(PAM.train,PAM.TestData$x,Threshold,type="class")
     PAM.PosteriorProbabilities=pamr.predict(PAM.train,PAM.TestData$x,Threshold,type="posterior")
     PAM.Centroid=pamr.predict(PAM.train,PAM.TestData$x,Threshold,type="centroid") # this gives the unshrunken centroids
     PAM.Nonzero=PAM.TestData$genenames[pamr.predict(PAM.train,PAM.TestData$x,Threshold,type="nonzero")]
     
     # Centroid outputs the unshrunken centroids, PAM.Nonzero has the genes, so save the shrunken ones here. 
     PAM.Centroid=PAM.Centroid[PAM.Nonzero,]
     
     PAM.TestPredictions=matrix(as.numeric(PAM.TestPredictions))
     rownames(PAM.TestPredictions)=colnames(MA_Data_Test)
     colnames(PAM.TestPredictions)="PAMprediction"
     
     #PAM.TestData$genenames[PAM.Nonzero]
     #PAM.AllTestPredictions<-array(dim=c(length(AllThresholds),1))
     #AllThresholds=PAM.cvresults$threshold
     #for (i in c(1:length(AllThresholds))) {
     #     results=pamr.predict(PAM.train,PAM.TestData$x, AllThresholds[i])
     #}
     
     #      # writing this to file
     
     
     #      write.table(AllResults,Args[8],sep="\t",row.names=TRUE)
     #      
     #      # writing also the PosteriorProbabilities to file
     #      Filename=paste(substr(Args[10], 1, nchar(Args[10])-4),"_PosteriorProbabilities.txt",sep='')
     #      write.table(PAM.PosteriorProbabilities,Filename,sep="\t",row.names=TRUE)
     
     return(list(TestPredictions=PAM.TestPredictions,PosteriorProbabilites=PAM.PosteriorProbabilities,Centroids=PAM.Centroid,Nonzero=PAM.Nonzero,OptimalSolution=OptimalSolution))
}

checkGeneSymbols.KB=function(x){
     x1=as.character(checkGeneSymbols(x)$Suggested.Symbol)
     ifelse(!is.na(x1), x1, x)
}

makedec=function(x){
     deciles = unique(quantile(x, seq(0, 10, 0.1), na.rm=T))
     ran_dec = as.numeric(cut(x, deciles, include.lowest = TRUE))
     return(ran_dec)
}

#Get nearby genes, input is 450k probes
GetnearbyGenes=function(probes1, upstreambp, downstreambp){
     library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
     ann450k=getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
     library(Homo.sapiens)
     org.Hs.egSYMBOL=as.data.frame(org.Hs.egSYMBOL)
     probes_annot=ann450k[probes1,c("chr","pos")]
     mycoords.gr=as.data.frame(cbind(probes_annot$chr, probes_annot$pos,probes_annot$pos))
     colnames(mycoords.gr)=c("chrom","start","end")
     mycoords.gr=as.data.table(mycoords.gr)
     mycoords.gr$start=as.numeric(as.character(mycoords.gr$start))-upstreambp
     mycoords.gr$end=as.numeric(as.character(mycoords.gr$end))+downstreambp
     mycoords.gr=makeGRangesFromDataFrame(mycoords.gr)
     OverlapGenes=subsetByOverlaps(genes(TxDb.Hsapiens.UCSC.hg19.knownGene), mycoords.gr)
     OverlappingGenes=org.Hs.egSYMBOL[org.Hs.egSYMBOL$gene_id %in% as.data.frame(OverlapGenes)$gene_id,"symbol"]
     OverlappingGenes=suppressWarnings(checkGeneSymbols(org.Hs.egSYMBOL[org.Hs.egSYMBOL$gene_id %in% as.data.frame(OverlapGenes)$gene_id,"symbol"])$Suggested.Symbol)
     OverlappingGenes=OverlappingGenes[which(OverlappingGenes %in% suppressWarnings(checkGeneSymbols(as.character(rownames(MAcancer)))$Suggested.Symbol))]
     OverlappingGenes=OverlappingGenes[!is.na(OverlappingGenes)]
     return(OverlappingGenes)
}

#
Process_Meth_data_KB <- function (MET_Data) { 
     library(impute)
     Genes=rownames(MET_Data)
     MET_Data=apply(MET_Data,2,as.numeric)
     rownames(MET_Data)=Genes
     
     SampleNames=colnames(MET_Data)
     
     # more than 10% = removal
     MissingValueThreshold=0.1
     
     # removing clones with too many missing values
     NrMissingsPerGene=apply(MET_Data,1,function(x) sum(is.na(x))) /ncol(MET_Data)
     cat("Removing",sum(NrMissingsPerGene>MissingValueThreshold),"genes with more than 10% missing values.\n")
     if (sum(NrMissingsPerGene>MissingValueThreshold)>0) MET_Data=MET_Data[NrMissingsPerGene<MissingValueThreshold,]
     
     # removing patients with too many missings values     
     NrMissingsPerSample=apply(MET_Data,2,function(x) sum(is.na(x))) /nrow(MET_Data)
     cat("Removing",sum(NrMissingsPerSample>MissingValueThreshold),"patients with more than 10% missing values.\n")
     if (sum(NrMissingsPerSample>MissingValueThreshold)>0) MET_Data=MET_Data[,NrMissingsPerSample<MissingValueThreshold]
     
     # knn impute using Tibshirani's method     
     k=15
     KNNresults=impute.knn(as.matrix(MET_Data),k) # this does not work, gives storage error. ???
     MET_Data_KNN=KNNresults$data
     
     # cleaning up sample names
     MET_Data_KNN_Clean=TCGA_GENERIC_CleanUpSampleNames(MET_Data_KNN,15)
     return(MET_Data_KNN_Clean)
}




#Make data frame indicating the genes that were used by the PAM_CV model to predict class in each round of cross validation
make.gene.ranks.array=function(PAMres){
     library(abind)
     library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
     ann450k = getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
     
     GenesRanks=PAMres$GenesRanks
     genes=unique(abind(GenesRanks, along=1)[,1])
     
     GenesRankArray=array(NA,c(length(genes),length(GenesRanks)))
     rownames(GenesRankArray)=genes
     colnames(GenesRankArray)=c(1:length(GenesRanks))
     for(i in 1:length(GenesRanks)){
          GenesRankArray[rownames(GenesRanks[[i]]),i]=as.numeric(as.character(GenesRanks[[i]][,"prop-selected-in-CV"]))
     }
     GenesRankArray2=as.matrix(apply(GenesRankArray,2,function(x) as.numeric(as.character(x))))
     rownames(GenesRankArray2)=rownames(GenesRankArray)
     GenesRankArray=GenesRankArray2
     
     #Add annotation for probes
     GenesRankArray=cbind(as.data.frame(GenesRankArray),ann450k[rownames(GenesRankArray),c("chr","pos","UCSC_RefGene_Name")])
     GenesRankArray$Gene=lapply(GenesRankArray$UCSC_RefGene_Name, function(x) paste(unique(unlist(strsplit(x,";"))), collapse=","))
     GenesRankArray$Gene=as.character(GenesRankArray$Gene)
     
     return(GenesRankArray)
}

checkGeneSymbols.KB=function(x){
     x1=as.character(checkGeneSymbols(x)$Suggested.Symbol)
     ifelse(!is.na(x1), x1, x)
}



PAM_Analysis_MaxNrFeatures_CV_KB2 <- function (MA_Data,CL_Data,NrFolds=10,NrLambdaFolds=10,MaxNrFeatures) {
  
  # overlapping samples
  OverlapSamples=intersect(colnames(MA_Data),colnames(CL_Data))
  if (length(OverlapSamples)>0) {
    MA_Data=MA_Data[,OverlapSamples]
    CL_Data=CL_Data[,OverlapSamples,drop=FALSE]
    cat('There are',length(OverlapSamples),'overlapping samples between predictor data and outcome data.\n')     
  } else {
    cat('No samples are overlapping, returning.')     
    return()          
  }
  
  # creating a frame structure for calculating the folds. 
  CL_Data_frame=data.frame(t(CL_Data))     
  names(CL_Data_frame)[1]="outcome"
  
  # Check for Leave-One-Out CV
  if (NrFolds >ncol(CL_Data)) { 
    cat("There are not enough samples ???.\n")
    return()
  }
  if (NrFolds==ncol(CL_Data)) {
    cat('Doing leave-one-out cross validation.\n')
    CL_Data_frame$folds=1:ncol(CL_Data)
  } else {
    cat('Doing stratified nfold cross validation.\n')          
    # using a split
    CL_Data_frame <- ddply(CL_Data_frame,.(outcome),createFolds,k = NrFolds)
    rownames(CL_Data_frame)=colnames(CL_Data)
  }
  
  Predictions=matrix(0,ncol(CL_Data),1)
  rownames(Predictions)=colnames(CL_Data)
  colnames(Predictions)=rownames(CL_Data)
  NrGenesPerFold=matrix(0,NrFolds,1)
  Genes=list()
  GenesRanks=list()
  
  PosteriorProbabilities=matrix(0,ncol(CL_Data),length(unique(CL_Data[1,])))
  rownames(PosteriorProbabilities)=colnames(CL_Data)
  #colnames(PosteriorProbabilities)=rownames(CL_Data)
  
  for (i in 1:NrFolds) {
    print(i)
    TestSamples=rownames(CL_Data_frame)[CL_Data_frame$folds==i]
    TrainSamples=setdiff(rownames(CL_Data_frame),TestSamples)
    PAMresults=PAM_Analysis_MaxNrFeatures_KB2(MA_Data[,TrainSamples],CL_Data[,TrainSamples,drop=FALSE],MA_Data[,TestSamples,drop=FALSE],NrLambdaFolds,MaxNrFeatures, HeteroNorm=NULL)          
    Predictions[TestSamples,1]=PAMresults$TestPredictions
    PosteriorProbabilities[TestSamples,1:ncol(PosteriorProbabilities)]=PAMresults$PosteriorProbabilites
    NrGenesPerFold[i]=PAMresults$OptimalSolution$minGenes
    #Need to work on getting the genes from this, make object to 
    Genes[[i]]=PAMresults$Nonzero
    GenesRanks[[i]]=as.data.frame(as.matrix(PAMresults$PAM.Gene.Ranks))
  }    
   NrClasses=length(unique(CL_Data[1,]))     
  
  NrMisclassificationsPerClass=matrix(0,NrClasses,4)
  colnames(NrMisclassificationsPerClass)=c('NrMistakes','TotalNrCasesInClass','NrMisclassificationsPerc','Accuracy')
  tmpData=split(Predictions,CL_Data)
  ClassNames=c()
  for (i in 1:NrClasses) {
    NrMisclassificationsPerClass[i,1]=sum(tmpData[[i]]!=i)
    NrMisclassificationsPerClass[i,2]=length(tmpData[[i]])
    NrMisclassificationsPerClass[i,3]=NrMisclassificationsPerClass[i,1]/NrMisclassificationsPerClass[i,2]
    NrMisclassificationsPerClass[i,4]=1-NrMisclassificationsPerClass[i,3]
    ClassNames=c(ClassNames,paste('Class',i,sep=''))
  }
  rownames(NrMisclassificationsPerClass)=ClassNames
  
  OverallPerformance=matrix(0,1,3)
  colnames(OverallPerformance)=c('NrMisclassifications','NrMisclassificationsPerc','Accuracy')
  OverallPerformance[1]=sum(Predictions-t(CL_Data) !=0)
  OverallPerformance[2]=OverallPerformance[1]/ncol(CL_Data)
  OverallPerformance[3]=1-OverallPerformance[1]/ncol(CL_Data)
  
  return(list(Predictions=Predictions,PosteriorProbabilities=PosteriorProbabilities,OverallPerformance=OverallPerformance,NrMisclassificationsPerClass=NrMisclassificationsPerClass,NonZero.Genes=Genes,GenesRanks=GenesRanks, NrGenesPerFold=NrGenesPerFold,CL_Data_Overlap=CL_Data))
}

createFolds <- function(x,k){
  n <- nrow(x)
  x$folds <- rep(1:k,length.out = n)[sample(n,n)]
  x
}

PAM_Analysis_MaxNrFeatures_KB2 <- function (MA_Data_Train,CL_Data_Train,MA_Data_Test,NrLambdaFolds=5,MaxNrFeatures, HeteroNorm=NULL) {
  
  # Args[8]="/Users/ogevaert/Documents/WORK/Projects/Stanford/TCGAglioblastoma/data/2014-06-13/OutcomeDataTestPredictions_1.txt"
  # Args[9]="/Users/ogevaert/Documents/WORK/Projects/Stanford/TCGAglioblastoma/results/2014-06-13/CentroidFigure_1.tif"
  # Args[10]="/Users/ogevaert/Documents/WORK/Projects/Stanford/TCGAglioblastoma/results/2014-06-13/CentroidData_1.txt"
  
  # overlapping samples
  OverlapSamples=intersect(colnames(MA_Data_Train),colnames(CL_Data_Train))
  if (length(OverlapSamples)>0) {
    MA_Data_Train=MA_Data_Train[,OverlapSamples]
    CL_Data_Train=CL_Data_Train[,OverlapSamples,drop=FALSE]
    cat('There are',length(OverlapSamples),'overlapping samples between predictor data and outcome data.\n')     
  } else {
    cat('No samples are overlapping, returning.')     
    return()          
  }
  
  dim(MA_Data_Train)
  dim(CL_Data_Train)
  ModuleDataDimensions=dim(MA_Data_Train)
  
  dim(MA_Data_Test)
  ModuleDataDimensionsTest=dim(MA_Data_Test)
  
  # Creating the data object
  PAM.TrainData=list(x=data.matrix(MA_Data_Train),y=as.numeric(data.matrix(CL_Data_Train)),genenames=as.character(row.names(MA_Data_Train)),geneid=as.character(row.names(MA_Data_Train)))
  
  # Training the PAM model
  if(missing(HeteroNorm)){
    PAM.train <- pamr.train(PAM.TrainData)
  }else{
    PAM.train <- pamr.train(PAM.TrainData, hetero=HeteroNorm)
  }
  
  # cross validation
  PAM.cvresults<-pamr.cv(PAM.train,PAM.TrainData,NrLambdaFolds)
  #There were 30 thresholds
  
  # choosing threshold when balanced error rate is used. Get minimal BER from 30 thresholds
  BERresults=CalculateBalancedErrorRate(PAM.cvresults)
  BER=BERresults$BER
  #BER=1-BERresults$Specificity     
  
  minOneVariablePositions=PAM.cvresults$size>=2 & PAM.cvresults$size<=MaxNrFeatures # at least two genes in the PAM model
  minError=min(BER[minOneVariablePositions])
  tmpPos=which(BER==minError & PAM.cvresults$size>=2 & PAM.cvresults$size<=MaxNrFeatures)
  
  #When you have a small number of samples, the same error often shows up for several thresholds, so you can't just select the row with that error
  minGenes=PAM.cvresults$size[tmpPos][1]
  #It's cheating because PAM reports the lowest number of genes with that error rate, but then uses a lower thresholds for cross validation. Need to use the treshold with reported number of genes
  #this is how many genes were used for the model with the minimal BER rate
  Thresholds=PAM.cvresults$threshold[tmpPos] 
  Threshold=Thresholds[1]
  #This is the least stringent treshold with a number of genes equal to or lesss than the MaxNrFeature, but with more than 2 genes
  #!bug in code here, arbitraily chooses the first threshold here if there are multiple. This is not typically the thrshold at which mingene was selected, so you get a larger number of genes than specified
  OptimalSolution=list(minGenes=minGenes,minError=minError,Threshold=Threshold)
  #Changed treshold to one with min genes; previous version arbitrarily chose first thresholds in set of tresholds with min genes, which is always the lowest of 
  #those tresholds with less than MaxNrFeatures. So a less stringent treshold was actually being applied to predict the treshold than was suggested by MaxNrFeatures
  #Threshold=PAM.cvresults$threshold[tmpPos][which.min(PAM.cvresults$size[tmpPos])]
  
  #PAM.trainOptimal <- pamr.train(PAM.TrainData,threshold=OptimalSolution$Threshold)
  
  ## Compute the confusion matrix for a particular model
  #pamr.confusion(PAM.cvresults, Threshold)
  
  ## Plot the cross-validated class probabilities by class
  #pamr.plotcvprob(PAM.cvresults, PAM.TrainData, Threshold)
  
  ## Plot the class centroids
  #pamr.plotcen(PAM.train, PAM.TrainData, Threshold)
  
  ## Make a gene plot of the most significant genes (less usefull)
  #pamr.geneplot(PAM.train, PAM.TrainData, Threshold)
  
  # Estimate false discovery rates and plot them
  #fdr.obj<- pamr.fdr(PAM.train, PAM.TrainData)
  
  #pamr.plotfdr(fdr.obj)
  
  #############################################################################################################
  ##### PAM predicting test set
  ##############################################################################################################
  
  # making the test data set
  PAM.TestData=list(x=data.matrix(MA_Data_Test),genenames=row.names(MA_Data_Test))
  
  # predicting the cohort 2 samples
  PAM.TestPredictions=pamr.predict(PAM.train,PAM.TestData$x,Threshold,type="class")
  PAM.PosteriorProbabilities=pamr.predict(PAM.train,PAM.TestData$x,Threshold,type="posterior")
  PAM.Centroid=pamr.predict(PAM.train,PAM.TestData$x,Threshold,type="centroid") # this gives the unshrunken centroids
  PAM.Nonzero=PAM.TestData$genenames[pamr.predict(PAM.train,PAM.TestData$x,Threshold,type="nonzero")]
  
  #get the ranking of the genes across earch round of cross validation# Nonzero genes (those that passed the threshold) should include the top ranking genes
  PAM.g=pamr.listgenes(PAM.train, PAM.TrainData,  Threshold, PAM.cvresults, genenames=FALSE)
  rownames(PAM.g)=PAM.g[,1]
  #Throwing error:  length of 'dimnames' [2] not equal to array extent. Think its where there are thresholds that don' pass
 
  # Centroid outputs the unshrunken centroids, PAM.Nonzero has the genes, so save the shrunken ones here. 
  PAM.Centroid=PAM.Centroid[PAM.Nonzero,]
  
  PAM.TestPredictions=matrix(as.numeric(PAM.TestPredictions))
  rownames(PAM.TestPredictions)=colnames(MA_Data_Test)
  colnames(PAM.TestPredictions)="PAMprediction"
  
  #PAM.TestData$genenames[PAM.Nonzero]
  #PAM.AllTestPredictions<-array(dim=c(length(AllThresholds),1))
  #AllThresholds=PAM.cvresults$threshold
  #for (i in c(1:length(AllThresholds))) {
  #     results=pamr.predict(PAM.train,PAM.TestData$x, AllThresholds[i])
  #}
  
  #      # writing this to file
  
  
  #      write.table(AllResults,Args[8],sep="\t",row.names=TRUE)
  #      
  #      # writing also the PosteriorProbabilities to file
  #      Filename=paste(substr(Args[10], 1, nchar(Args[10])-4),"_PosteriorProbabilities.txt",sep='')
  #      write.table(PAM.PosteriorProbabilities,Filename,sep="\t",row.names=TRUE)
  
  return(list(TestPredictions=PAM.TestPredictions,PosteriorProbabilites=PAM.PosteriorProbabilities,Centroids=PAM.Centroid,Nonzero=PAM.Nonzero,PAM.Gene.Ranks=PAM.g,OptimalSolution=OptimalSolution))
  
}

##################################################################
#MethylMix no script
#################################################################

MethylMix_No_filter=function(cancer, METpath, EXPpath){
dat=load(paste(METpath, cancer, "/gdac_20160715/ProcessedData_",cancer,".RData",sep=""))
METcancer=ProcessedData$MET_Data_Cancer
colnames(METcancer)=substr(colnames(METcancer),1,12)

METnormal=ProcessedData$MET_Data_Normal
colnames(METnormal)=substr(colnames(METnormal),1,12)

#Removing SNP probes, as they seem to stop MethylMix running, not sure why
METcancer=METcancer[-c(grep("rs", rownames(METcancer))),]
METnormal=METnormal[-c(grep("rs", rownames(METnormal))),]

OverlapProbes=intersect(rownames(METcancer),rownames(METnormal))
METcancer=METcancer[OverlapProbes,]
METnormal=METnormal[OverlapProbes,]

filelist=paste(EXPpath, cancer,"/gdac_20160715/", sep="")
file=list.files(filelist, pattern=".Rdata|.RData")[1]
dat=load(paste(filelist, file, sep=""))

MAcancer=ProcessedData$MA_Data_Cancer
colnames(MAcancer)=substr(colnames(MAcancer),1,12)

library(parallel)
library(doParallel)
library(RPMM)

MethResults=MethylMix("univ.beta", METcancer, METnormal, MAcancer, Parallel = T, filter="none2")
save(MethResults, file=paste("/srv/gevaertlab/data/TCGA/Kevin_Data/Full_MethylMix_No_filter_",cancer,".RData", sep=""))
return(MethResults=MethResults)
}


####
# Function: MethylMix
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# This function implements...
#
# Arguments:
#   method: "univ.beta", "univ.m", "biv.m", version of MethylMix to run
#   METcancer, METnormal, MAcancer: required data sets
#   MAnormal: optional data set, gene expression in normal samples
#   NoNormalMode: do not compare to normal samples to define states (default: false)
#   test: how to do comparison to normal: "wilcoxon" (univariate, only with MET
#           data, default if MAnormal not provided), "ttest" (univariate, only
#           with MET data), "manova" (bivariate, MA and MET data, default if
#           MA is provided, if MA is not provided is changed to wilcoxon)
#   maxComp: up to which number of components should we keep testing if adding a
#               new component improves the fit (default = 3)
#   OutputRoot, Parallel: taken from original package, so far only available
#                            for the original version, univariate with beta values
#   filter: "lin.reg", "none", "none2"
#
# Returns:
#   A list containing:
#       - 
####

MethylMix <- function(method, METcancer, METnormal, MAcancer, MAnormal = NULL, 
                      NoNormalMode = FALSE, test = NULL, maxComp = 3,
                      OutputRoot='', Parallel = FALSE, filter = "lin.reg", listOfGenes = NULL, 
                      UseChisquare = TRUE) { 
    
    # Package needed for the original version
    if (method == "univ.beta") {
        if(require("RPMM")){
            cat("RPMM is loaded correctly\n")
        } else {
            cat("trying to install mclust\n")
            install.packages("RPMM")
            if(require("RPMM")){
                cat("RPMM installed and loaded\n")
            } else {
                stop("could not install RPMM")
            }
        }
    }
    
    # Things necessary to do only for the m-values versions
    if (method == "univ.m" || method == "biv.m") {
        
        # Required packages
        # Package mclust for fitting mixture model
        if(require("mclust")){
            cat("mclust is loaded correctly\n")
        } else {
            cat("trying to install mclust\n")
            install.packages("mclust")
            if(require(mclust)){
                cat("mclust installed and loaded\n")
            } else {
                stop("could not install mclust")
            }
        }
        
        # Take care of beta values out of range and transform into M values
        METcancer = getMvalues(METcancer)
        METnormal = getMvalues(METnormal)
        
        # Things necessary to do only for the bivariate versions
        if (method == "biv.m") {       
            # Keep only those samples with both METcancer and MAcancer data
            OverlapSamples = intersect(colnames(METcancer), colnames(MAcancer))
            MAcancer = MAcancer[, OverlapSamples, drop = FALSE]
            METcancer = METcancer[, OverlapSamples, drop = FALSE]
            # Same for normal samples if MAnormal is provided
            if (!is.null(MAnormal)) {
                OverlapSamples = intersect(colnames(METnormal), colnames(MAnormal))      
                # If the intersection is empty, we won't be able to use the bivariate normal data,
                # use only the methylation normal data, as if MAnormal wasn't provided
                if (length(OverlapSamples) == 0) {
                    MAnormal = NULL
                    cat("Normal methylation and microarray data come from different samples. Microarray data for normal samples will not be used.\n")
                } else {
                    MAnormal = MAnormal[, OverlapSamples, drop = FALSE]
                    METnormal = METnormal[, OverlapSamples, drop = FALSE]   
                }
            }
        }
    }
    
    # Step 1: modeling the gene expression using cancer methylation data
    if (!is.null(listOfGenes)) {
        # ModelGeneExpression iterates through each gene in METcancer and then looks it up in MAcancer,
        # then I only need to modify METcancer to have the genes I'm interested in.
        METcancer = METcancer[listOfGenes, , drop = F]
    }
    if (filter == "lin.reg") {
        # This is the default, the same as in original MethylMix
        FunctionalGenes = MethylMix_ModelGeneExpression(METcancer, MAcancer, Method = "Regression", NULL)
    }
    else if (filter == "lin.reg.both.dir") {
        # This is the default, the same as in original MethylMix
        FunctionalGenes = MethylMix_ModelGeneExpression_BothDir(METcancer, MAcancer, Method = "Regression", NULL)
    }
    else if (filter == "none") {
        # No filtering, run mixture model in all genes.
        # In the option above, when we run MethylMix_ModelGeneExpression, the 
        # resulting FunctionalGenes appear in both METcancer and MAcancer. So
        # for this non-filter option I will put in FunctionalGenes all the genes
        # shared between both data sets
        METcancer.split.names  <- sapply(strsplit(rownames(METcancer),  '---'), function(x) x[1])
        genes.to.keep.MET = METcancer.split.names %in% rownames(MAcancer)
        FunctionalGenes = rownames(METcancer)[genes.to.keep.MET]
    }
    else if (filter == "none2") {
        # No filtering, run mixture model in all genes.
        # In this one I don't do the intersection between METcancer and MAcancer
        FunctionalGenes = rownames(METcancer)
    }
    #     else if (filter == "nat.spl") {
    #         FunctionalGenes = SplineFilter(METcancer, MAcancer)
    #     }
    else stop("Invalid argument for filter.")
    
    
    # Step 2: modeling the methylation data as a mixture of beta/normal distributions 
    if (length(FunctionalGenes) > 0) {
        if (method == "univ.beta") {
            MixtureModelResults = MethylMix_MixtureModel(METcancer, METnormal, FunctionalGenes, Parallel, 1, length(FunctionalGenes), NoNormalMode, maxComp, UseChisquare)
        }
        else if (method == "univ.m") {
            MixtureModelResults = Mixt.Model.Univ.Mvalues(METcancer, METnormal, FunctionalGenes, Parallel, NoNormalMode, test, maxComp)
        }
        else if (method == "biv.m") {
            MixtureModelResults = Mixt.Model.Biv.Mvalues(METcancer, METnormal, MAcancer, MAnormal, FunctionalGenes, test, NoNormalMode, maxComp)
        }
        else stop("Invalid method argument.")
    } else {
        cat("No functional genes found.\n")
    }
    
    # Detach packages
    if (method == "univ.m" || method == "biv.m") {
        detach("package:mclust", unload = TRUE)
    }
    
    # Writing to file
    if (length(OutputRoot) > 1) {
        MethylMix_WriteToFile(OutputRoot, MixtureModelResults)
    }
    
    #KB addition so that it doesn't break if there are no MethylMix genes
    if (!exists("MixtureModelResults"))
    {
        cat("No MethylMix genes")
    } else {
        return(MixtureModelResults)
    }
}

####
# Function: MethylMix_ModelGeneExpression
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Remains the same as in the original package
####
MethylMix_ModelGeneExpression <- function(METcancer, MAcancer, Method = c("Regression", "Pearson", "Spearman"), CovariateData = NULL) {      
    
    # overlapping samples     
    OverlapSamples = intersect(colnames(METcancer), colnames(MAcancer))
    cat("Found", length(OverlapSamples), "samples with both methylation and expression data.\n")
    MAcancer = MAcancer[, OverlapSamples, drop = FALSE]
    METcancer = METcancer[, OverlapSamples, drop = FALSE]
    if (!is.null(CovariateData)) CovariateData = as.matrix(CovariateData[OverlapSamples, ])
    
    Rsquares = matrix(0, nrow = length(rownames(METcancer)), ncol = 1)    
    Genes = rownames(METcancer)  
    PvalueThreshold = 0.001  
    RsquareThreshold = 0.1
    
    pb = txtProgressBar(1, length(rownames(METcancer)), style = 3, width = 100)
    cat("Correlating methylation data with gene expression.\n")
    for(i in 1:length(rownames(METcancer))) {
        setTxtProgressBar(pb, i)
        tmpGene = unlist(strsplit(Genes[i], '---'))[1]
        pos = which(rownames(MAcancer) == tmpGene)
        
        if (length(pos) > 0) {
            if (!is.null(CovariateData)) {
                res = lm(MAcancer[pos, ] ~ METcancer[Genes[i], ] + CovariateData)
            } else {           
                res = lm(MAcancer[pos, ] ~ METcancer[Genes[i], ])
            }
            res.summary = summary(res)
            if (!is.null(CovariateData)) {
                if (res$coefficients[2] < 0 & res.summary$coefficients[2, 4] < PvalueThreshold & abs(res.summary$coefficients[2, 1]) > abs(res.summary$coefficients[3, 1])) { # methylation effect bigger than tissue
                    Rsquares[i] = res.summary$r.squared
                }
            } else if (res$coefficients[2] < 0 & res.summary$coefficients[2, 4] < PvalueThreshold) {
                Rsquares[i] = res.summary$r.squared    
            }
        }
    }
    
    # Rsquare threshold
    FunctionalGenes = Genes[Rsquares > RsquareThreshold]
    cat("\nFound", length(FunctionalGenes), "functional genes.\n")
    return(FunctionalGenes)
}

MethylMix_ModelGeneExpression_BothDir <- function(METcancer, MAcancer, Method = c("Regression", "Pearson", "Spearman"), CovariateData = NULL) {      
    
    # overlapping samples     
    OverlapSamples = intersect(colnames(METcancer), colnames(MAcancer))
    cat("Found", length(OverlapSamples), "samples with both methylation and expression data.\n")
    MAcancer = MAcancer[, OverlapSamples, drop = FALSE]
    METcancer = METcancer[, OverlapSamples, drop = FALSE]
    if (!is.null(CovariateData)) CovariateData = as.matrix(CovariateData[OverlapSamples, ])
    
    Rsquares = matrix(0, nrow = length(rownames(METcancer)), ncol = 1)    
    Genes = rownames(METcancer)  
    PvalueThreshold = 0.001  
    RsquareThreshold = 0.1
    
    pb = txtProgressBar(1, length(rownames(METcancer)), style = 3, width = 100)
    cat("Correlating methylation data with gene expression.\n")
    for(i in 1:length(rownames(METcancer))) {
        setTxtProgressBar(pb, i)
        tmpGene = unlist(strsplit(Genes[i], '---'))[1]
        pos = which(rownames(MAcancer) == tmpGene)
        
        if (length(pos) > 0) {
            if (!is.null(CovariateData)) {
                res = lm(MAcancer[pos, ] ~ METcancer[Genes[i], ] + CovariateData)
            } else {           
                res = lm(MAcancer[pos, ] ~ METcancer[Genes[i], ])
            }
            res.summary = summary(res)
            if (!is.null(CovariateData)) {
                if (res.summary$coefficients[2, 4] < PvalueThreshold & abs(res.summary$coefficients[2, 1]) > abs(res.summary$coefficients[3, 1])) { # methylation effect bigger than tissue
                    Rsquares[i] = res.summary$r.squared
                }
            } else if (res.summary$coefficients[2, 4] < PvalueThreshold) {
                Rsquares[i] = res.summary$r.squared    
            }
        }
    }
    
    # Rsquare threshold
    FunctionalGenes = Genes[Rsquares > RsquareThreshold]
    cat("\nFound", length(FunctionalGenes), "functional genes.\n")
    return(FunctionalGenes)
}

####
# Function: Mixt.Model.Univ.Mvalues
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Arguments:
#   METcancer, METnormal: data sets
#   FunctionalGenes: vector with the names of the functional genes found
#   Parallel: logical for parallel computing, not implemented yet
#   NoNormalMode: do not compare to normal sample to define states
#   test: "wilcoxon" (default) or "ttest", test to use to compare to normal
#   maxComp: up to which number of components should we keep testing if adding a
#               new component improves the fit (default = 3)
#
# Returns:
#   MethylationStates: matrix with DM values, rows are driven genes, cols are samples
#   NrComponents: matrix with the number of components identified for each driver gene
#   Models: list of the output objects from the mixture model fit
#   MethylationDrivers: vector of strings with the names of the driver genes
#   MixtureStates: list, DM values for every component in each driver gene
#   AllFlipOverStates: should be removed, I used this when working with the flip over function
#   Parameters: list, for each gene, a list with the vectors mean, sd (estimated
#               parameters from the mixture model) and prop (proportion of samples
#               in each component)
#   Classification: list, for each gene, a vector indicating to which component
#                   each sample was assigned
####
Mixt.Model.Univ.Mvalues <- function(METcancer, METnormal, FunctionalGenes, Parallel = FALSE, NoNormalMode = FALSE, test = NULL, maxComp = 3) {
    
    # overlap of genes
    GeneOverlap = intersect(rownames(METcancer), rownames(METnormal))
    METcancer = METcancer[GeneOverlap, , drop = FALSE]
    METnormal = METnormal[GeneOverlap, , drop = FALSE]
    
    cat("\nStarting Gaussian mixture modeling.\n")
    METcancer = METcancer[which(rownames(METcancer) %in% FunctionalGenes), , drop = FALSE]  
    METnormal = METnormal[which(rownames(METnormal) %in% FunctionalGenes), , drop = FALSE]
    
    MethylationStates = matrix(0, length(rownames(METcancer)), length(colnames(METcancer)))    
    rownames(MethylationStates) = rownames(METcancer)
    colnames(MethylationStates) = colnames(METcancer)
    
    AllNrComponents = matrix(1, length(rownames(METcancer)), 1)
    AllFlipOverStates = matrix(1, length(rownames(METcancer)), 1)
    AllMixtureStates = array(list(), length(rownames(METcancer)))
    Models = array(list(), length(rownames(METcancer)))
    Parameters = array(list(), length(rownames(METcancer)))
    Classification = array(list(), length(rownames(METcancer)))
    GeneNames = rownames(MethylationStates)
    
    cat("Running Gaussian mixture model on", length(rownames(METcancer)), "functional genes and on", length(colnames(METcancer)), "samples.\n")
    
    options(warn = -1) 
    for(i in 1:length(rownames(METcancer))) {
        MixtureModelResults_SingleGene = ModelSingleGene.Univ.Mvalues(GeneNames[i], METcancer[i,], METnormal[i,], NoNormalMode, test, maxComp)
        MethylationStates[i, ] = MixtureModelResults_SingleGene$MethylationState
        Models[[i]] = MixtureModelResults_SingleGene$Model
        AllNrComponents[i, 1] = MixtureModelResults_SingleGene$NrComponents
        AllFlipOverStates[i, 1] = MixtureModelResults_SingleGene$FlipOverState
        AllMixtureStates[[i]] = MixtureModelResults_SingleGene$MixtureStates
        Parameters[[i]] = MixtureModelResults_SingleGene$parameters
        Classification[[i]] = MixtureModelResults_SingleGene$classification
    } 
    options(warn = 0)
    
    # Removing the genes without any differential methylation. 
    FunctionalGenes = rownames(METcancer)
    if (!NoNormalMode) {  # I use NoNormalMode = T when I run for all normals, in that case I don't want to remove any gene
        NonZeroPositions = rowSums(MethylationStates)!= 0
        FunctionalGenes = rownames(METcancer)[NonZeroPositions]
        AllNrComponents = as.matrix(AllNrComponents[NonZeroPositions, ])
        AllFlipOverStates = as.matrix(AllFlipOverStates[NonZeroPositions, ])
        AllMixtureStates = AllMixtureStates[NonZeroPositions]
        Models = Models[NonZeroPositions]
        Parameters = Parameters[NonZeroPositions]
        Classification = Classification[NonZeroPositions]
        MethylationStates = MethylationStates[NonZeroPositions, , drop = FALSE]
    }
    
    rownames(AllNrComponents) = FunctionalGenes
    rownames(AllFlipOverStates) = FunctionalGenes
    names(AllMixtureStates) = FunctionalGenes
    names(Models) = FunctionalGenes
    names(Parameters) = FunctionalGenes
    names(Classification) = FunctionalGenes
    MethylationDrivers = rownames(MethylationStates)
    
    return(list(MethylationStates = MethylationStates,
                NrComponents = AllNrComponents,
                Models = Models,
                MethylationDrivers = MethylationDrivers,
                MixtureStates = AllMixtureStates,
                AllFlipOverStates = AllFlipOverStates,  # I added AllFlipOverStates here to test it, it's not returned in original MethylMix
                Parameters = Parameters,                # I added this to have easy access to estimates, access that don't depend on the package and model used
                Classification = Classification))       # The resulting classification
}

####
# Function: ModelSingleGene.Univ.Mvalues
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# For each gene, this function fits the mixture model, selects the number of 
# components and defines the corresponging methylation states.
#
# Arguments:
#   GeneName: string with the name of the gene
#   METdataVector, METdataNormalVector: vectors witht the data for GeneName
#   NoNormalMode: do not compare to normal sample to define states
#   test: "wilcoxon" (default) or "ttest", test to use to compare to normal
#   maxComp: up to which number of components should we keep testing if adding a
#               new component improves the fit (default = 3)
#
# Returns:
#   MethylationState: DM value for each sample
#   NrComponents: number of components identified
#   Model: output object from the mixture model fit
#   MixtureState: DM values for every component
#   FlipOverState: FlipOverState
#   parameters: list with the vectors mean, sd (estimated parameters from the 
#               mixture model) and prop (proportion of samples in each component)
#   classification: a vector indicating to which component each sample was assigned
####
# i=1;GeneName=GeneNames[i]; METdataVector=METcancer[i,] ;METdataNormalVector=METnormal[i,]
ModelSingleGene.Univ.Mvalues <- function(GeneName, METdataVector, METdataNormalVector, NoNormalMode = FALSE, test = NULL, maxComp = 3) {
    
    OrigOrder = order(METdataVector)
    PvalueThreshold = 0.01
    MeanDifferenceTreshold = 0.10  # This is in Beta-value scale
    
    # Package "mclust"
    # Advantage: by default selects automatically the number of components using BIC (not used here)
    # and the default initialization is to classify the data by quantile.
    # Disadvantage:  I wasn't able to provide in the "initialization" argument other type of 
    # initial values or initial estimates of the parameters of the model.
    # It produces a warning when we force it to use a fixed number of components that is not appropriate
    # saying that one of the components has no observation ("In map(out$z) : no assignment to ...")
    
    
    # 1 component model
    mods = vector("list", maxComp + 1)
    mods[[1]] = Mclust(METdataVector, G = 1, model = "V")
    bic = numeric(maxComp + 1)
    bic[1] = - mods[[1]]$bic # mclust provides BIC with opposite sign
    # 2- to maxComp+1- components model
    for (comp in 2:(maxComp + 1)) {
        res = Mclust(METdataVector, G = comp, model = "V")
        if (is.null(res)) mods[[comp]] = NA  else mods[[comp]] = res
        if (is.na(mods[[comp]])) {
            # Try again with a conjugate prior on the variances
            # cat("Trying again.\n")
            res = Mclust(METdataVector, G = comp, model = "V", prior=priorControl(scale = rep(1, comp)))
            if (is.null(res)) mods[[comp]] = NA  else mods[[comp]] = res
            # If it didn't converge again:
            if (is.na(mods[[comp]])) {    
                cat("Not able to fit ", comp, " components. ", comp-1, " component", ifelse(comp-1 == 1, "", "s"), " selected.\n")
                NrComponents = comp - 1
                break
            }
        }
        bic[comp] = - mods[[comp]]$bic
        model.means = summary(mods[[comp]])$mean
        model.means.in.beta = sort((2^model.means) / (2^model.means + 1))
        not.different.means = ifelse(all(abs(diff(model.means.in.beta)) > MeanDifferenceTreshold), F, T)
        
        # As in original MethylMix, we try adding another component if the following 2 conditions are satisfied:
        #   A: Adding one component reduces BIC
        #   B: All absolute differences between means in model with one extra component are above the MeanDifferenceThreshold
        # Then continue if A & B, else not continue, which is the same as saying not continue if !A OR !B                
        if (bic[comp] >= bic[comp - 1] | not.different.means) {
            NrComponents = comp - 1
            break
        }
        if (comp == maxComp + 1) NrComponents = maxComp + 1 # If we get here it means that maxComp + 1 is better than maxComp components, and we should try one more, but we won't
    }                
    
    # If more than maxComp components are better, provide a message but use only maxComp
    if (NrComponents == maxComp + 1) {
        cat(paste(GeneName, ": more than ", maxComp, " components provide optimal BIC, but ", maxComp, " are applied.\n"))
        NrComponents = maxComp
    }
    mod = mods[[NrComponents]]
    
    # Results of the classification
    classification = mod$classification
    # If there's only one component mclust doesn't provide the names of the samples, I add  it here
    if (NrComponents == 1) names(classification) = rownames(mod$data)
    
    # List with means, sd and mixing proportions
    parameters = list(mean = summary(mod)$mean, sd = sqrt(summary(mod)$var), prop = summary(mod)$pro)
    
    MethylationState = matrix(0, 1, length(METdataVector))
    MixtureStates = matrix(0, NrComponents, 1)
    
    for (comp in 1:NrComponents) {
        METdataVector_comp = METdataVector[classification == comp]
        
        res=list(p.value = 1)
        if (length(METdataVector_comp) > 0) {   
            # Choose between doing the comparison to normal with non-parametric Wilcoxon or T-test
            # If user didn't provide argument test, test is NULL and will be set to Wilcoxon (default)
            if (is.null(test)) test = "wilcoxon"
            if (test == "wilcoxon") res = wilcox.test(METdataVector_comp, METdataNormalVector)
            else if (test == "ttest") {
                res = try(t.test(METdataVector_comp, METdataNormalVector))
                if (class(res) == "try-error") {
                    cat("Not enough observations for t-test. Wilcoxon performed.\n")
                    res = wilcox.test(METdataVector_comp, METdataNormalVector)
                }
            }
            else stop("Invalid test argument")
        }
        
        # Methylation means in M-value and Beta-value
        M.means = c(mean(METdataVector_comp), mean(METdataNormalVector))
        Beta.means = (2^M.means) / (2^M.means + 1)
        # Difference in Beta-scale should be greater that threshold (0.1)
        # Difference between means in M and Beta values
        Difference = M.means[1] - M.means[2]
        Beta.Difference = Beta.means[1] - Beta.means[2]
        if ((res$p.value < PvalueThreshold & abs(Beta.Difference) > MeanDifferenceTreshold) | NoNormalMode) {
            MethylationState[1, classification == comp] = Difference            
            MixtureStates[comp, 1] = Difference
        }
    }
    
    FlipOverState = 0
    if (length(unique(MixtureStates)) >= 2) {
        # correcting the flipover effect
        FlipOverResults = MethylMix_RemoveFlipOver2(OrigOrder, MethylationState, METdataVector, NrComponents)
        MethylationState = FlipOverResults$MethylationState
        FlipOverState = FlipOverResults$LearnedState     
    }
    
    message = ifelse(NrComponents == 1, " component is best.\n", " components are best.\n")
    cat(c(GeneName, ": ", NrComponents, message), sep = "")
    
    return(list(MethylationState = MethylationState,
                NrComponents = NrComponents,
                Model = mod,
                MixtureStates = MixtureStates,
                FlipOverState = FlipOverState,
                parameters = parameters,
                classification = classification))
}

####
# Function: MethylMix_RemoveFlipOver2
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# This function deals with the classification of the samples when the overlap
# of the estimated densities overlap.
# It has the number 2 in the name to differentiate it from the function that 
# already exists in the original package.
####
MethylMix_RemoveFlipOver2 <- function(OrigOrder, MethylationState, METdataVector, NrComponents, UseTrainedFlipOver = FALSE, FlipOverState = 0) {
    
    Differences = diff(MethylationState[1, OrigOrder]) 
    DifferencesSel = Differences[Differences != 0]  
    LearnedState = 0
    
    # If we have 2 components in the model, we perform the original flip over modification in Methylmix
    # This handles the escenario where one of the components has observations at both sides of the other
    if (NrComponents == 2) {
        # If (length(DifferencesSel) < 2, there's no mixing
        if (length(DifferencesSel) == 2) {
            if (DifferencesSel[1] * -1 == DifferencesSel[2]) {
                
                posDiff1 = which(Differences == DifferencesSel[1])
                stateSize1 = posDiff1
                posDiff2 = which(Differences == DifferencesSel[2])
                stateSize2 = length(Differences) - posDiff2
                
                if (UseTrainedFlipOver == TRUE) {
                    if (FlipOverState == 1) {
                        # I changed OrigOrder[posDiff2 + 1:length(Differences)] into OrigOrder[(posDiff2 + 1) : (length(Differences) + 1)] so that index is not out of range and NAs are not introduced (it didn't affect anything anyway)
                        MethylationState[1, OrigOrder[(posDiff2 + 1) : (length(Differences) + 1)]] = MethylationState[1, OrigOrder[posDiff2]]
                    } else if (FlipOverState == 2) {
                        MethylationState[1, OrigOrder[1:posDiff1]] = MethylationState[1, OrigOrder[posDiff1 + 1]]
                    }                    
                } else {
                    if (stateSize2 > stateSize1) {                    
                        MethylationState[1, OrigOrder[1:posDiff1]] = MethylationState[1, OrigOrder[posDiff1 + 1]]
                        LearnedState = 2
                    } else if (stateSize1 > stateSize2) {                    
                        MethylationState[1, OrigOrder[(posDiff2 + 1) : (length(Differences) + 1)]] = MethylationState[1, OrigOrder[posDiff2]]
                        LearnedState = 1
                    }  
                }
            }
        }
    } 
    
    # For models with 3 components, we have this new implementation that handles escenarios where only one of the
    # components is divided and separated between the other 2.
    if (NrComponents == 3) {
        # If length(DifferencesSel) == 2, there is no overlapping and mixing
        if (length(DifferencesSel) > 2) {
            
            # seq.of.states: Shows the order in which the different states are being mixing
            # posDiff: for each sequence of the same state, which is the position of the last element
            seq.of.states = numeric(0)
            posDiff = numeric(0)
            for (i in 1:length(DifferencesSel)) {
                posDiff[i] = which(Differences == DifferencesSel[i])
                seq.of.states[i] = MethylationState[1, OrigOrder][posDiff[i]]
            }
            # Catch the last group:
            seq.of.states = c(seq.of.states, tail(MethylationState[1, OrigOrder], 1))
            posDiff = c(0, posDiff, length(MethylationState))
            
            # Size of each group:
            size = diff(posDiff)
            
            # A table of seq.of.states will show which state is separated
            tab = table(seq.of.states)
            
            # I will deal only with the case where only one state is separated
            if (sum(tab > 1) == 1) {
                
                # Mean of each subgroup
                means = tapply(METdataVector[OrigOrder], rep(1:length(size), size), mean)
                
                # Identify which is the state that is divided into subgroups
                separated.state = round(as.numeric(names(tab)[tab > 1]), 4)
                
                # Subgroups that correspond to this separated state
                subgr = which(round(seq.of.states, 4) == separated.state)
                
                # Subgroups well defined and separated
                subgr.ok = which(round(seq.of.states, 4) != separated.state)
                
                # In the separated state, the largest subgroup will remain, and the others
                # will be allocated to one of the ok subgroups, to the one of closest mean
                subgr.remains = subgr[which.max(size[subgr])]
                subgr.allocate = subgr[!subgr %in% subgr.remains]
                for (gr in subgr.allocate) {
                    allocate.in.group = subgr.ok[which.min(abs(means[gr] - means[subgr.ok]))]
                    pos.to.change = (posDiff[gr] + 1) : posDiff[gr + 1]
                    MethylationState[1, OrigOrder[pos.to.change]] = seq.of.states[allocate.in.group]
                }
                LearnedState = 3 # I chose a 3 just to distiguish from the 1 or 2 that the original flipover function returns
            }
        }
    }
    return(list(MethylationState = MethylationState, LearnedState = LearnedState))
}

####
# Function: Mixt.Model.Biv.Mvalues
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Arguments:
#   METcancer, METnormal, MAcancer: required data sets
#   MAnormal: optional data set, gene expression in normal samples
#   FunctionalGenes: vector with the names of the functional genes found
#   test: how to do comparison to normal: "wilcoxon" (univariate, only with MET
#           data, default if MAnormal not provided), "ttest" (univariate, only
#           with MET data), "manova" (bivariate, MA and MET data, default if
#           MA is provided, if MA is not provided is changed to wilcoxon)
#   NoNormalMode: do not compare to normal samples to define states (default: false)
#   maxComp: up to which number of components should we keep testing if adding a
#               new component improves the fit (default = 3)
#
# Returns:
#   MethylationStates: matrix with DM values, rows are driven genes, cols are samples
#   NrComponents: matrix with the number of components identified for each driver gene
#   Models: list of the output objects from the mixture model fit
#   MethylationDrivers: vector of strings with the names of the driver genes
#   MixtureStates: list, DM values for every component in each driver gene
#   Parameters: list, for each gene, a list with the vectors mean, sd (estimated
#               parameters from the mixture model) and prop (proportion of samples
#               in each component)
#   Classification: list, for each gene, a vector indicating to which component
#                   each sample was assigned
####
Mixt.Model.Biv.Mvalues <- function(METcancer, METnormal, MAcancer, MAnormal = NULL, FunctionalGenes, test = NULL, NoNormalMode = FALSE, maxComp = 3) {
    
    # overlap of genes between all data sets
    overlapMET = intersect(rownames(METcancer), rownames(METnormal))
    overlapMA = rownames(MAcancer)
    if (!is.null(MAnormal)) overlapMA = intersect(rownames(MAcancer), rownames(MAnormal))
    
    # Intersect genes in MET and MA data sets without loosing the "---Cluster" in MET data sets
    overlapMETsplit  <- sapply(strsplit(overlapMET,  '---'), function(x) x[1])
    genes.to.keep.MET = overlapMET[overlapMETsplit %in% overlapMA]
    genes.to.keep.MA = overlapMA[overlapMA %in% overlapMETsplit]
    METcancer = METcancer[genes.to.keep.MET, , drop = FALSE]
    METnormal = METnormal[genes.to.keep.MET, , drop = FALSE]
    MAcancer = MAcancer[genes.to.keep.MA, , drop = FALSE]
    if (!is.null(MAnormal)) MAnormal = MAnormal[genes.to.keep.MA, , drop = FALSE]
    
    # Keep only those genes that were identified as functional
    METcancer = METcancer[rownames(METcancer) %in% FunctionalGenes, , drop = FALSE]  
    METnormal = METnormal[rownames(METnormal) %in% FunctionalGenes, , drop = FALSE]
    FunctionalGenesMA  <- sapply(strsplit(FunctionalGenes,  '---'), function(x) x[1])
    MAcancer = MAcancer[rownames(MAcancer) %in% FunctionalGenesMA, , drop = FALSE]
    if (!is.null(MAnormal)) MAnormal = MAnormal[rownames(MAnormal) %in% FunctionalGenesMA, , drop = FALSE]
    # METcancer and METnormal will have more rows (genes) than MAcancer and MAnormal
    # if there are cases where there are more than one in MET for a gene in MA, 
    # like "name---Cluster1", "name---Cluster2" in MET and "name" in MA
    
    MethylationStates = matrix(0, length(rownames(METcancer)), length(colnames(METcancer)))    
    rownames(MethylationStates) = rownames(METcancer)
    colnames(MethylationStates) = colnames(METcancer)
    
    AllNrComponents = matrix(1, length(rownames(METcancer)), 1)
    AllMixtureStates = array(list(), length(rownames(METcancer)))
    Models = array(list(), length(rownames(METcancer)))
    Parameters = array(list(), length(rownames(METcancer)))
    Classification = array(list(), length(rownames(METcancer)))
    GeneNamesMET = rownames(METcancer)
    GeneNamesMA = sapply(strsplit(GeneNamesMET,  '---'), function(x) x[1])
    
    cat("Running Gaussian mixture model on", length(rownames(METcancer)), "functional genes and on", length(colnames(METcancer)), "samples.\n")
    options(warn = -1)
    for(i in 1:length(GeneNamesMET)) {
        geneNameMET = GeneNamesMET[i]
        geneNameMA = GeneNamesMA[i]
        
        if (!is.null(MAnormal)) { 
            MAdataNormalVector = MAnormal[geneNameMA, ] 
        } else MAdataNormalVector = NULL
        Biv.Mix.Mod.Res.SingleGene = ModelSingleGene.Biv.Mvalues(geneNameMET, METcancer[geneNameMET, ], MAcancer[geneNameMA, ], METnormal[geneNameMET, ], MAdataNormalVector, test, NoNormalMode, maxComp)
        MethylationStates[i, ] = Biv.Mix.Mod.Res.SingleGene$MethylationState
        Models[[i]] = Biv.Mix.Mod.Res.SingleGene$Model
        AllNrComponents[i, 1] = Biv.Mix.Mod.Res.SingleGene$NrComponents
        AllMixtureStates[[i]] = Biv.Mix.Mod.Res.SingleGene$MixtureStates
        Parameters[[i]] = Biv.Mix.Mod.Res.SingleGene$parameters
        Classification[[i]] = Biv.Mix.Mod.Res.SingleGene$classification
    }         
    options(warn = 0)
    
    # Removing the genes without any differential methylation. 
    FunctionalGenes = rownames(METcancer)
    if (!NoNormalMode) {  # I use NoNormalMode = T when I run for all normals, in that case I don't want to remove any gene
        NonZeroPositions = rowSums(MethylationStates)!= 0
        FunctionalGenes = rownames(METcancer)[NonZeroPositions] # No FunctionalGenes[NonZeroPositions] because FunctionalGenes can be longer than the genes evaluated by SingleGene, which are the ones in METcancer, METnormal and MAcancer
        AllNrComponents = as.matrix(AllNrComponents[NonZeroPositions, ])
        AllMixtureStates = AllMixtureStates[NonZeroPositions]
        Models = Models[NonZeroPositions]
        Parameters = Parameters[NonZeroPositions]
        Classification = Classification[NonZeroPositions]
        MethylationStates = MethylationStates[NonZeroPositions, , drop = FALSE]
    }
    
    rownames(AllNrComponents) = FunctionalGenes
    names(AllMixtureStates) = FunctionalGenes
    names(Models) = FunctionalGenes
    names(Parameters) = FunctionalGenes
    names(Classification) = FunctionalGenes
    MethylationDrivers = rownames(MethylationStates)
    
    return(list(MethylationStates = MethylationStates,
                NrComponents = AllNrComponents,
                Models = Models,
                MethylationDrivers = MethylationDrivers,
                MixtureStates = AllMixtureStates,
                Parameters = Parameters,                # I added this to have easy access to estimates, access that don't depend on the package and model used
                Classification = Classification))       # The resulting classification
}

####
# Function: ModelSingleGene.Biv.Mvalues
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# For each gene, this function fits the mixture model, selects the number of 
# components and defines the respective methylation states. 
#
# Arguments:
#   GeneName: string with the name of the gene
#   METdataVector, MAdataVector, METdataNormalVector: required vectors with the data for GeneName
#   MAdataNormalVector: optional, gene expression in normal samples
#   test: how to do comparison to normal: "wilcoxon" (univariate, only with MET
#           data, default if MAnormal not provided), "ttest" (univariate, only
#           with MET data), "manova" (bivariate, MA and MET data, default if
#           MA is provided, if MA is not provided is changed to wilcoxon)
#   NoNormalMode: do not compare to normal samples to define states (default: false)
#   maxComp: up to which number of components should we keep testing if adding a
#               new component improves the fit (default = 3)
#
# Returns:
#   MethylationState: DM value for each sample
#   NrComponents: number of components identified
#   Model: output object from the mixture model fit
#   MixtureState: DM values for every component
#   FlipOverState: FlipOverState
#   parameters: list with the vectors mean, sd (estimated parameters from the 
#               mixture model) and prop (proportion of samples in each component)
#   classification: a vector indicating to which component each sample was assigned
####

# GeneName=GeneNames[i]; METdataVector=METcancer[i,] ; MAdataVector = MAcancer[i, ]; METdataNormalVector=METnormal[i,]
ModelSingleGene.Biv.Mvalues <- function(GeneName, METdataVector, MAdataVector, METdataNormalVector, MAdataNormalVector = NULL, test = NULL, NoNormalMode = FALSE, maxComp = 3) {
    
    PvalueThreshold = 0.01
    MeanDifferenceTreshold = 0.10  # This is in Beta-value scale
    
    # 1 component model
    mods = vector("list", maxComp + 1)
    mods[[1]] = Mclust(data.frame(METdataVector, MAdataVector), G = 1, modelNames = "VVV")
    bic = numeric(maxComp + 1)
    bic[1] = - mods[[1]]$bic # mclust provides BIC with opposite sign                     
    
    # 2- to maxComp+1- components model
    for (comp in 2:(maxComp + 1)) {
        res = Mclust(data.frame(METdataVector, MAdataVector), G = comp, modelNames = "VVV")
        if (is.null(res)) mods[[comp]] = NA  else mods[[comp]] = res
        if (is.na(mods[[comp]])) {
            # Try again with a conjugate prior on the variances
            # cat("Trying again.\n")
            res = Mclust(data.frame(METdataVector, MAdataVector), G = comp, modelNames = "VVV", prior=priorControl(scale = diag(1, comp)))
            if (is.null(res)) mods[[comp]] = NA  else mods[[comp]] = res
            # If it didn't converge again:
            if (is.na(mods[[comp]])) {    
                cat("Not able to fit ", comp, " components. ", comp-1, " component", ifelse(comp-1 == 1, "", "s"), " selected.\n")
                NrComponents = comp - 1
                break
            }                
        }
        bic[comp] = - mods[[comp]]$bic
        
        model.means = summary(mods[[comp]])$mean
        met.means.in.beta = sort((2^model.means[1, ]) / (2^model.means[1, ] + 1))
        not.different.met.means = ifelse(all(abs(diff(met.means.in.beta)) > MeanDifferenceTreshold), F, T)
        
        # As in original MethylMix, we try adding another component if the following 2 conditions are satisfied:
        #   A: Adding one component reduces BIC
        #   B: All absolute differences between methylation means in model with one extra component are above the MeanDifferenceThreshold
        # Then continue if A & B, else not continue, which is the same as saying not continue if !A OR !B
        if (bic[comp] >= bic[comp - 1] | not.different.met.means) {
            NrComponents = comp - 1
            break
        }
        if (comp == maxComp + 1) NrComponents = maxComp + 1 # If we get here it means that maxComp + 1 is better than maxComp components, and we should try one more, but we won't
    }                
    
    # If more than maxComp components are better, provide a message but use only maxComp
    if (NrComponents == maxComp + 1) {
        cat(paste(GeneName, ": more than ", maxComp, " components provide optimal BIC, but ", maxComp, " are applied.\n"))
        NrComponents = maxComp
    }
    mod = mods[[NrComponents]]
    
    # Results of the classification
    classification = mod$classification
    if (NrComponents == 1) names(classification) = names(METdataVector)
    # List with means, sd and mixing proportions
    sigma = list()
    for (l in 1:NrComponents) sigma[[l]] = summary(mod)$var[,,l]
    parameters = list(mean = summary(mod)$mean, sigma = sigma, prop = summary(mod)$pro)
    
    MethylationState = matrix(0, 1, length(METdataVector))
    MixtureStates = matrix(0, NrComponents, 1)
    
    # If I don't have MAnormal, I do like in the univariate case, take MET data and
    # compare to normal with Wilcoxon, but here I added the option to perform a t-test.
    # If I have MAnormal, instead of Wilc/t-test I do a MANOVA to compare to normal both MET and MA data
    # In any case, to decide if I have a new methylation state, I use the pvalue and 
    # a threshold for the difference in the met means as we did in the univariate case.
    # We also have the NoNormalMode, if true nothing of this is considered
    for (comp in 1:NrComponents) {
        res=list(p.value = 1)
        METdataVector_comp = METdataVector[classification == comp]
        if (length(METdataVector_comp) > 0) {
            
            # Argument test not provided and MA data available, do manova
            if (is.null(test) & !is.null(MAdataNormalVector)) test = "manova"
            # Argument test not provided and MA data not available, do wilcoxon only with MET data
            if (is.null(test) & is.null(MAdataNormalVector)) test = "wilcoxon"
            # User said "manova" but not MA data available, then do it only with MET data, with Wilcoxon
            if (test == "manova" & is.null(MAdataNormalVector)) test = "wilcoxon"
            
            gaveError = FALSE
            if (test == "manova") {
                MAdataVector_comp = MAdataVector[classification == comp]
                met = c(METdataVector_comp, METdataNormalVector)
                ma = c(MAdataVector_comp, MAdataNormalVector)
                type = rep(c("cancer", "normal"), c(length(METdataVector_comp), length(METdataNormalVector)))
                manovaRes = try(summary(manova(cbind(met, ma) ~ type)), silent=T)
                if (class(manovaRes) == "try-error"){
                    gaveError = TRUE # in once case there was a component with only one data point and this produce error
                } else {
                    res = list(p.value = manovaRes$stats[1, 6])
                    if (is.null(res$p.value)) gaveError = TRUE # Sometimes when there's only one data point manova does not through error and the pvalue is NULL
                }
            } else if (test == "wilcoxon") {
                res = wilcox.test(METdataVector_comp, METdataNormalVector)
            } else if (test == "ttest") {
                ttest = try(t.test(METdataVector_comp, METdataNormalVector), silent=TRUE)
                # this produces error if there's only one data point in the component
                if (class(ttest) == "try-error") {
                    gaveError = TRUE
                } else {
                    res = ttest
                }
            } else stop("Invalid test argument")
            
            if (gaveError) {
                # do wilcoxon
                res = wilcox.test(METdataVector_comp, METdataNormalVector)
            }
        }
        
        # Methylation means in M-value and Beta-value
        M.means = c(mean(METdataVector_comp), mean(METdataNormalVector))
        Beta.means = (2^M.means) / (2^M.means + 1)
        
        # Difference in Beta-scale should be greater that threshold (0.1)
        # Difference between means in M and Beta values
        Difference = M.means[1] - M.means[2]
        Beta.Difference = Beta.means[1] - Beta.means[2]
        
        if ((res$p.value < PvalueThreshold & abs(Beta.Difference) > MeanDifferenceTreshold) | NoNormalMode) {
            MethylationState[1, classification == comp] = Difference            
            MixtureStates[comp, 1] = Difference
        }
    }  
    
    message = ifelse(NrComponents == 1, " component is best.\n", " components are best.\n")
    cat(c(GeneName, ": ", NrComponents, message), sep = "")
    
    return(list(MethylationState = MethylationState,
                NrComponents = NrComponents,
                Model = mod,
                MixtureStates = MixtureStates,
                parameters = parameters,
                classification = classification))
}


#Start=1;Stop=length(FunctionalGenes)
#MethylMix_MixtureModel(METcancer, METnormal, FunctionalGenes, Parallel, 1, length(FunctionalGenes), NoNormalMode, maxComp)
MethylMix_MixtureModel <- function(METcancer,METnormal,FunctionalGenes,Parallel=FALSE,Start,Stop,NoNormalMode=FALSE,maxComp=3,UseChisquare=TRUE) {
    
    # overlap of samples
    GeneOverlap=intersect(rownames(METcancer),rownames(METnormal))
    METcancer=METcancer[GeneOverlap,,drop=FALSE]
    METnormal=METnormal[GeneOverlap,,drop=FALSE]
    
    if (Stop>length(FunctionalGenes)) Stop=length(FunctionalGenes)
    cat("\nStarting Beta mixture modeling.\n")
    METcancer=METcancer[which(rownames(METcancer) %in% FunctionalGenes),,drop=FALSE]  
    METnormal=METnormal[which(rownames(METnormal) %in% FunctionalGenes),,drop=FALSE]
    
    MethylationStates=matrix(0,length(rownames(METcancer)),length(colnames(METcancer)))    
    rownames(MethylationStates)=rownames(METcancer)
    colnames(MethylationStates)=colnames(METcancer)
    
    AllNrComponents=matrix(1,length(rownames(METcancer)),1)
    AllFlipOverStates=matrix(1,length(rownames(METcancer)),1)
    AllMixtureStates=array(list(),length(rownames(METcancer)))
    Models=array(list(),length(rownames(METcancer)))
    GeneNames=rownames(MethylationStates)
    
    cat("Running Beta mixture model on",length(rownames(METcancer)),"functional genes and on",length(colnames(METcancer)),"samples.\n")
    #cat("Running genes from position",Start,"to",Stop,".\n")
    if (Parallel & .Platform$OS.type!="windows") {
        #           tmp=tryCatch(
        #                {
        # registering an appropriate nr of cores
        NrCores=detectCores()
        NrCores=NrCores-2
        if (NrCores<2) NrCores=1
        registerDoParallel(cores=NrCores) # this will register nr of cores/threads
        
        # alternative to do it in parallel
        tmpResults=foreach::foreach(i=Start:Stop) %dopar% {      #length(rownames(METcancer))               
            MixtureModelResults_SingleGene=MethylMix_ModelSingleGene(GeneNames[i],METcancer[i,],METnormal[i,],NoNormalMode,maxComp,UseChisquare=UseChisquare)
        }
        for (i in 1:(Stop-Start+1)) {   #length(rownames(METcancer))
            MethylationStates[Start+i-1,]=tmpResults[[i]]$MethylationState
            Models[[Start+i-1]]=tmpResults[[i]]$Model
            AllNrComponents[Start+i-1,1]=tmpResults[[i]]$NrComponents
            AllFlipOverStates[Start+i-1,1]=tmpResults[[i]]$FlipOverState
            AllMixtureStates[[Start+i-1]]=tmpResults[[i]]$MixtureStates
            #                     }  
            #                }, finally = {
            #                     if (.Platform$OS.type=="windows") {
            #                          stopImplicitCluster()
            #                     }
        }
        #           )
    } else {
        #pb=txtProgressBar(1,length(rownames(METcancer)),style=3,width=100)
        for(i in Start:Stop) {      #length(rownames(METcancer))               
            #setTxtProgressBar(pb,i)
            MixtureModelResults_SingleGene=MethylMix_ModelSingleGene(GeneNames[i],METcancer[i,],METnormal[i,],NoNormalMode,maxComp,UseChisquare=UseChisquare)
            MethylationStates[i,]=MixtureModelResults_SingleGene$MethylationState
            Models[[i]]=MixtureModelResults_SingleGene$Model
            AllNrComponents[i,1]=MixtureModelResults_SingleGene$NrComponents
            AllFlipOverStates[i,1]=MixtureModelResults_SingleGene$FlipOverState
            AllMixtureStates[[i]]=MixtureModelResults_SingleGene$MixtureStates
        }   
    }   
    
    # removing the genes without any differential methylation. 
    if (!NoNormalMode) {  # I use NoNormalMode = T when I run for all normals, in that case I don't want to remove any gene
        NonZeroPositions=rowSums(MethylationStates)!=0
        AllNrComponents=as.matrix(AllNrComponents[NonZeroPositions,])
        AllFlipOverStates=as.matrix(AllFlipOverStates[NonZeroPositions,])
        FunctionalGenes=FunctionalGenes[NonZeroPositions]
        AllMixtureStates=AllMixtureStates[NonZeroPositions]
        Models=Models[NonZeroPositions]
        MethylationStates=MethylationStates[NonZeroPositions,,drop=FALSE]
    }
    rownames(AllNrComponents)=FunctionalGenes
    MethylationDrivers=rownames(MethylationStates)
    
    return(list(MethylationStates=MethylationStates,NrComponents=AllNrComponents,Models=Models,MethylationDrivers=MethylationDrivers,MixtureStates=AllMixtureStates))
}

#i=1;GeneName=GeneNames[i];METdataVector=METcancer[i,];METdataNormalVector=METnormal[i,]
MethylMix_ModelSingleGene <-function(GeneName,METdataVector,METdataNormalVector,NoNormalMode=FALSE,maxComp=3,UseChisquare=TRUE) {
    
    OrigOrder = order(METdataVector)
    ChiSquareThreshold = ifelse(UseChisquare, qchisq(0.95, df = 3), 0)
    PvalueThreshold = 0.01
    MeanDifferenceTreshold = 0.10
    
    # 1 component model
    w0.m = matrix(1, length(METdataVector), 1)
    mods = vector("list", maxComp + 1)
    mods[[1]] = blc_2(matrix(METdataVector, ncol = 1), w = w0.m, maxiter = 100, tol = 1e-06, verbose = FALSE)
    bic = numeric(maxComp + 1)
    bic[1] = -2 * mods[[1]]$llike + 2 * log(length(METdataVector))  
    
    # 2 to maxComp components model
    for (comp in 2:maxComp) {
        
        # Divide initial groups using quantiles
        prob = seq(1/comp, (comp-1)/comp, 1/comp)
        tmpQuantiles = quantile(METdataVector, prob)
        w0.m = matrix(0, nrow = length(METdataVector), ncol = comp)
        tmpQuantiles = c(tmpQuantiles, Inf)
        w0.m[METdataVector < tmpQuantiles[1], 1] <- 1
        for (i in 2:comp) {
            w0.m[METdataVector >= tmpQuantiles[i - 1] & METdataVector < tmpQuantiles[i], i] <- 1
        }
        
        # Fit beta mixture model
        mods[[comp]] <- blc_2(matrix(METdataVector, ncol=1), w = w0.m, maxiter = 100, tol = 1e-06, verbose = FALSE)
        if (sum(is.na(mods[[comp]]$mu)) > 0) mods[[comp]]$llike=0
        df = comp * 3 - 1
        bic[[comp]] = -2 * mods[[comp]]$llike + df * log(length(METdataVector))
        
        # See differences between model's means to compare them to threshold
        model.means = sort(mods[[comp]]$mu)
        different.means = ifelse(all(abs(diff(model.means)) > MeanDifferenceTreshold), T, F)
        
        # If the model was improved try one more, else break
        if ((bic[[comp - 1]] - bic[[comp]]) <= ChiSquareThreshold | !different.means) {
            NrComponents = comp - 1
            break
        }
        NrComponents = comp
    }
    
    Model = mods[[NrComponents]]
    MethylationState = matrix(0, 1, length(METdataVector))
    FlipOverState = 0
    res = list(p.value = 1)
    MixtureStates = matrix(0, NrComponents, 1)
    
    if (NrComponents == 1) {
        res = wilcox.test(METdataVector, METdataNormalVector)
        Difference = mean(METdataVector) - mean(METdataNormalVector)
        if ((res$p.value < PvalueThreshold & abs(Difference) > MeanDifferenceTreshold) | NoNormalMode) {
            MethylationState[1, ] = Difference
            MixtureStates[1, 1] = Difference
        }
        cat(c(GeneName,": 1 component is best.\n"))
    } else {
        for (comp in 1:NrComponents) {
            METdataVector_comp = METdataVector[Model$w[, comp] > (1/NrComponents)]
            if (length(METdataVector_comp) > 0) res = wilcox.test(METdataVector_comp, METdataNormalVector) else res$p.value = 1
            Difference = mean(METdataVector_comp) - mean(METdataNormalVector)
            if ((res$p.value < PvalueThreshold & abs(Difference) > MeanDifferenceTreshold) | NoNormalMode) {
                MethylationState[1, Model$w[, comp] > (1/NrComponents)] = Difference            
                MixtureStates[comp, 1] = Difference
            }  
        }
        # Flipover in original package only when there are two components
        if (NrComponents == 2) {
            FlipOverResults=MethylMix_RemoveFlipOver(OrigOrder,MethylationState)
            MethylationState=FlipOverResults$MethylationState
            FlipOverState=FlipOverResults$LearnedState             
        }
        cat(c(GeneName,": ", NrComponents, " components are best.\n"))
    }                 
    return(list(MethylationState = MethylationState,
                NrComponents = NrComponents,
                Model = Model,
                MixtureStates = MixtureStates,
                FlipOverState = FlipOverState))
}

MethylMix_RemoveFlipOver <-function(OrigOrder,MethylationState,UseTrainedFlipOver=FALSE,FlipOverState=0) {
    
    Differences=diff(MethylationState[1,OrigOrder])
    DifferencesSel=Differences[Differences!=0]
    LearnedState=0
    if (length(DifferencesSel)==2) {
        if (DifferencesSel[1]*-1==DifferencesSel[2]) {     
            
            posDiff1=which(Differences==DifferencesSel[1])
            stateSize1=posDiff1
            posDiff2=which(Differences==DifferencesSel[2])
            stateSize2=length(Differences)-posDiff2
            
            if (UseTrainedFlipOver==TRUE) {
                if (FlipOverState==1) {
                    MethylationState[1,OrigOrder[posDiff2+1:length(Differences)]]=MethylationState[1,OrigOrder[posDiff2]]
                } else if (FlipOverState==2) {
                    MethylationState[1,OrigOrder[1:posDiff1]]=MethylationState[1,OrigOrder[posDiff1+1]]
                }                    
            } else {
                if (stateSize2>stateSize1) {                    
                    MethylationState[1,OrigOrder[1:posDiff1]]=MethylationState[1,OrigOrder[posDiff1+1]]
                    LearnedState=2
                } else if (stateSize1>stateSize2) {                    
                    MethylationState[1,OrigOrder[posDiff2+1:length(Differences)]]=MethylationState[1,OrigOrder[posDiff2]]
                    LearnedState=1
                }  
            }
        }
    }
    return(list(MethylationState=MethylationState,LearnedState=LearnedState))
}  


blc_2 <- function (Y, w, maxiter = 25, tol = 1e-06, weights = NULL, verbose = TRUE) 
{
    Ymn <- min(Y[Y > 0], na.rm = TRUE)
    Ymx <- max(Y[Y < 1], na.rm = TRUE)
    Y <- pmax(Y, Ymn/2)
    Y <- pmin(Y, 1 - (1 - Ymx)/2) 
    Yobs = !is.na(Y)
    J <- dim(Y)[2] 
    K <- dim(w)[2]
    n <- dim(w)[1]
    if (n != dim(Y)[1]) 
        stop("Dimensions of w and Y do not agree")
    if (is.null(weights)) 
        weights <- rep(1, n)
    mu <- a <- b <- matrix(Inf, K, J)
    crit <- Inf
    for (i in 1:maxiter) {
        warn0 <- options()$warn
        options(warn = -1)
        eta <- apply(weights * w, 2, sum)/sum(weights)
        mu0 <- mu
        for (k in 1:K) {
            for (j in 1:J) {  
                ab <- betaEst_2(Y[, j], w[, k], weights)
                a[k, j] <- ab[1]
                b[k, j] <- ab[2]
                mu[k, j] <- ab[1]/sum(ab)
            }
        }     
        ww <- array(0, dim = c(n, J, K))
        for (k in 1:K) {
            for (j in 1:J) {
                ww[Yobs[, j], j, k] <- dbeta(Y[Yobs[, j], j], 
                                             a[k, j], b[k, j], log = TRUE)
            }
        }          
        options(warn = warn0)
        w <- apply(ww, c(1, 3), sum, na.rm = TRUE)
        wmax <- apply(w, 1, max)
        for (k in 1:K) w[, k] <- w[, k] - wmax
        w <- t(eta * t(exp(w)))
        like <- apply(w, 1, sum)
        w <- (1/like) * w
        llike <- weights * (log(like) + wmax)
        crit <- max(abs(mu - mu0))
        if (verbose) 
            print(crit)          
        if (is.na(crit) || crit < tol) 
            break
    }
    return(list(a = a, b = b, eta = eta, mu = mu, w = w, llike = sum(llike)))
}


betaEst_2 <-function (Y, w, weights) 
{
    y=Y
    yobs = !is.na(y)
    if (sum(yobs) <= 1) 
        return(c(1, 1))
    y = y[yobs]
    w = w[yobs]
    weights = weights[yobs]
    N <- sum(weights * w)
    p <- sum(weights * w * y)/N
    v <- sum(weights * w * y * y)/N - p * p
    logab <- log(c(p, 1 - p)) + log(pmax(1e-06, p * (1 - p)/v - 
                                             1))
    if (sum(yobs) == 2) 
        return(exp(logab))
    opt <- try(optim(logab, betaObjf, ydata = y, wdata = w, weights = weights, 
                     method = "BFGS"), silent = TRUE)
    if (inherits(opt, "try-error")) 
        return(c(1, 1)) 
    exp(opt$par) # if using optimx, exp(as.numeric(opt$par))
}


MethylMix_WriteToFile <-function(OutputRoot,MethylMixResults) {
    OutputFile=paste(OutputRoot,'_MethylationStates.txt',sep="")
    write.table(MethylMixResults$MethylationStates,OutputFile,sep="\t")
    
    OutputFile=paste(OutputRoot,"_NrComponents.txt",sep="")
    write.table(MethylMixResults$NrComponents,OutputFile,sep="\t")
    
    # also saving an Robject
    OutputFile=paste(OutputRoot,"_MethylMixResults.RData",sep="")
    save(MethylMixResults,file=OutputFile)
}


OneGene.Biv.Mvalues <- function(gene, METcancer, METnormal, MAcancer, MAnormal = NULL, 
                                NoNormalMode = FALSE, test = NULL, maxComp = 3,
                                OutputRoot='') {
    
    if(require("mclust")){
        cat("mclust is loaded correctly\n")
    } else {
        cat("trying to install mclust\n")
        install.packages("mclust")
        if(require(mclust)){
            cat("mclust installed and loaded\n")
        } else {
            stop("could not install mclust")
        }
    }
    
    index = c(grep(paste0("^", gene, "$") , rownames(METcancer)), grep(paste0("^", gene, "---") , rownames(METcancer)))
    if (length(index) == 0) stop(paste(gene, "not found in METcancer data."))
    METcancer = METcancer[index, , drop = F]
    index = c(grep(paste0("^", gene, "$") , rownames(METnormal)), grep(paste0("^", gene, "---") , rownames(METnormal)))
    if (length(index) == 0) stop(paste(gene, "not found in METnormal data."))
    METnormal = METnormal[index, , drop = F]
    index = c(grep(paste0("^", gene, "$") , rownames(MAcancer)), grep(paste0("^", gene, "---") , rownames(MAcancer)))
    if (length(index) == 0) stop(paste(gene, "not found in MAcancer data."))
    MAcancer = MAcancer[index, , drop = F]
    if (!is.null(MAnormal)) {
        index = c(grep(paste0("^", gene, "$") , rownames(MAnormal)), grep(paste0("^", gene, "---") , rownames(MAnormal)))
        if (length(index) == 0) stop(paste(gene, "not found in MAnormal data."))
        MAnormal = MAnormal[index, , drop = F]
    }
    
    # Take care of beta values out of range and transform into M values
    METcancer = getMvalues(METcancer)
    METnormal = getMvalues(METnormal)
    
    # Keep only those samples with both METcancer and MAcancer data
    OverlapSamples = intersect(colnames(METcancer), colnames(MAcancer))
    MAcancer = MAcancer[, OverlapSamples, drop = FALSE]
    METcancer = METcancer[, OverlapSamples, drop = FALSE]
    # Same for normal samples if MAnormal is provided
    if (!is.null(MAnormal)) {
        OverlapSamples = intersect(colnames(METnormal), colnames(MAnormal))      
        # If the intersection is empty, we won't be able to use the bivariate normal data,
        # use only the methylation normal data, as if MAnormal wasn't provided
        if (length(OverlapSamples) == 0) {
            MAnormal = NULL
            cat("Normal methylation and microarray data come from different samples. Microarray data for normal samples will not be used.\n")
        } else {
            MAnormal = MAnormal[, OverlapSamples, drop = FALSE]
            METnormal = METnormal[, OverlapSamples, drop = FALSE]   
        }
    }
    
    MethylationStates = matrix(0, length(rownames(METcancer)), length(colnames(METcancer)))    
    rownames(MethylationStates) = rownames(METcancer)
    colnames(MethylationStates) = colnames(METcancer)
    
    AllNrComponents = matrix(1, length(rownames(METcancer)), 1)
    rownames(AllNrComponents) = rownames(METcancer)
    AllMixtureStates = array(list(), length(rownames(METcancer)))
    names(AllMixtureStates) = rownames(METcancer)
    Models = array(list(), length(rownames(METcancer)))
    names(Models) = rownames(METcancer)
    Parameters = array(list(), length(rownames(METcancer)))
    names(Parameters) = rownames(METcancer)
    Classification = array(list(), length(rownames(METcancer)))
    names(Classification) = rownames(METcancer)
    GeneNamesMET = rownames(METcancer)
    GeneNamesMA = sapply(strsplit(GeneNamesMET,  '---'), function(x) x[1])
    MethylationDrivers = rownames(METcancer)    
    
    options(warn = -1)
    for(i in 1:length(GeneNamesMET)) {
        geneNameMET = GeneNamesMET[i]
        geneNameMA = GeneNamesMA[i]
        
        if (!is.null(MAnormal)) { 
            MAdataNormalVector = MAnormal[geneNameMA, ] 
        } else MAdataNormalVector = NULL
        Biv.Mix.Mod.Res.SingleGene = ModelSingleGene.Biv.Mvalues(geneNameMET, METcancer[geneNameMET, ], MAcancer[geneNameMA, ], METnormal[geneNameMET, ], MAdataNormalVector, test, NoNormalMode, maxComp)
        MethylationStates[i, ] = Biv.Mix.Mod.Res.SingleGene$MethylationState
        Models[[i]] = Biv.Mix.Mod.Res.SingleGene$Model
        AllNrComponents[i, 1] = Biv.Mix.Mod.Res.SingleGene$NrComponents
        AllMixtureStates[[i]] = Biv.Mix.Mod.Res.SingleGene$MixtureStates
        Parameters[[i]] = Biv.Mix.Mod.Res.SingleGene$parameters
        Classification[[i]] = Biv.Mix.Mod.Res.SingleGene$classification
    }         
    options(warn = 0)
    
    return(list(MethylationStates = MethylationStates,
                NrComponents = AllNrComponents,
                Models = Models,
                MethylationDrivers = MethylationDrivers,
                MixtureStates = AllMixtureStates,
                Parameters = Parameters,
                Classification = Classification))
}

#gene="HOXB4"; NoNormalMode = FALSE; test = NULL; maxComp = 3;
OneGene.Univ.Mvalues <- function(gene, METcancer, METnormal, # MAcancer, MAnormal = NULL, WE DON'T USE MA DATA SETS HERE
                                 NoNormalMode = FALSE, test = NULL, maxComp = 3,
                                 OutputRoot='') {
    
    if(require("mclust")){
        cat("mclust is loaded correctly\n")
    } else {
        cat("trying to install mclust\n")
        install.packages("mclust")
        if(require(mclust)){
            cat("mclust installed and loaded\n")
        } else {
            stop("could not install mclust")
        }
    }
    
    index = c(grep(paste0("^", gene, "$") , rownames(METcancer)), grep(paste0("^", gene, "---") , rownames(METcancer)))
    if (length(index) == 0) stop(paste(gene, "not found in METcancer data."))
    METcancer = METcancer[index, , drop = F]
    index = c(grep(paste0("^", gene, "$") , rownames(METnormal)), grep(paste0("^", gene, "---") , rownames(METnormal)))
    if (length(index) == 0) stop(paste(gene, "not found in METnormal data."))
    METnormal = METnormal[index, , drop = F]
#     index = c(grep(paste0("^", gene, "$") , rownames(MAcancer)), grep(paste0("^", gene, "---") , rownames(MAcancer)))
#     if (length(index) == 0) stop(paste(gene, "not found in MAcancer data."))
#     MAcancer = MAcancer[index, , drop = F]
#     if (!is.null(MAnormal)) {
#         index = c(grep(paste0("^", gene, "$") , rownames(MAnormal)), grep(paste0("^", gene, "---") , rownames(MAnormal)))
#         if (length(index) == 0) stop(paste(gene, "not found in MAnormal data."))
#         MAnormal = MAnormal[index, , drop = F]
#     }
    
    # Take care of beta values out of range and transform into M values
    METcancer = getMvalues(METcancer)
    METnormal = getMvalues(METnormal)
    
    MethylationStates = matrix(0, length(rownames(METcancer)), length(colnames(METcancer)))    
    rownames(MethylationStates) = rownames(METcancer)
    colnames(MethylationStates) = colnames(METcancer)
    AllNrComponents = matrix(1, length(rownames(METcancer)), 1)
    rownames(AllNrComponents) = rownames(METcancer)
    AllFlipOverStates = matrix(1, length(rownames(METcancer)), 1)
    rownames(AllFlipOverStates) = rownames(METcancer)
    AllMixtureStates = array(list(), length(rownames(METcancer)))
    names(AllMixtureStates) = rownames(METcancer)
    Models = array(list(), length(rownames(METcancer)))
    names(Models) = rownames(METcancer)
    Parameters = array(list(), length(rownames(METcancer)))
    names(Parameters) = rownames(METcancer)
    Classification = array(list(), length(rownames(METcancer)))
    names(Classification) = rownames(METcancer)
    GeneNames = rownames(MethylationStates)
    MethylationDrivers = GeneNames
    
    options(warn = -1) 
    for(i in 1:length(rownames(METcancer))) {
        MixtureModelResults_SingleGene = ModelSingleGene.Univ.Mvalues(GeneNames[i], METcancer[i,], METnormal[i,], NoNormalMode, test, maxComp)
        MethylationStates[i, ] = MixtureModelResults_SingleGene$MethylationState
        Models[[i]] = MixtureModelResults_SingleGene$Model
        AllNrComponents[i, 1] = MixtureModelResults_SingleGene$NrComponents
        AllFlipOverStates[i, 1] = MixtureModelResults_SingleGene$FlipOverState
        AllMixtureStates[[i]] = MixtureModelResults_SingleGene$MixtureStates
        Parameters[[i]] = MixtureModelResults_SingleGene$parameters
        Classification[[i]] = MixtureModelResults_SingleGene$classification
    } 
    options(warn = 0)
    
    return(list(MethylationStates = MethylationStates,
                NrComponents = AllNrComponents,
                Models = Models,
                MethylationDrivers = MethylationDrivers,
                MixtureStates = AllMixtureStates,
                AllFlipOverStates = AllFlipOverStates,  
                Parameters = Parameters,                
                Classification = Classification))
}

OneGene.Univ.Beta <- function(gene, METcancer, METnormal, #MAcancer, WE DON'T NEED MAcancer HERE
                              NoNormalMode = FALSE, test = NULL, maxComp = 3,
                              OutputRoot='', UseChisquare=FALSE) {
    
    if(require("RPMM")){
        cat("RPMM is loaded correctly\n")
    } else {
        cat("trying to install mclust\n")
        install.packages("RPMM")
        if(require("RPMM")){
            cat("RPMM installed and loaded\n")
        } else {
            stop("could not install RPMM")
        }
    }
    
    index = c(grep(paste0("^", gene, "$") , rownames(METcancer)), grep(paste0("^", gene, "---") , rownames(METcancer)))
    if (length(index) == 0) stop(paste(gene, "not found in METcancer data."))
    METcancer = METcancer[index, , drop = F]
    index = c(grep(paste0("^", gene, "$") , rownames(METnormal)), grep(paste0("^", gene, "---") , rownames(METnormal)))
    if (length(index) == 0) stop(paste(gene, "not found in METnormal data."))
    METnormal = METnormal[index, , drop = F]
#     index = c(grep(paste0("^", gene, "$") , rownames(MAcancer)), grep(paste0("^", gene, "---") , rownames(MAcancer)))
#     if (length(index) == 0) stop(paste(gene, "not found in MAcancer data."))
#     MAcancer = MAcancer[index, , drop = F]
    
    # overlap of samples
    GeneOverlap=intersect(rownames(METcancer),rownames(METnormal))
    METcancer=METcancer[GeneOverlap,,drop=FALSE]
    METnormal=METnormal[GeneOverlap,,drop=FALSE]
    
    MethylationStates=matrix(0,length(rownames(METcancer)),length(colnames(METcancer)))    
    rownames(MethylationStates)=rownames(METcancer)
    colnames(MethylationStates)=colnames(METcancer)
    
    AllNrComponents=matrix(1,length(rownames(METcancer)),1)
    rownames(AllNrComponents)=rownames(METcancer)
    AllFlipOverStates=matrix(1,length(rownames(METcancer)),1)
    AllMixtureStates=array(list(),length(rownames(METcancer)))
    Models=array(list(),length(rownames(METcancer)))
    GeneNames=rownames(MethylationStates)
    
    for(i in 1:length(rownames(METcancer))) {
        MixtureModelResults_SingleGene=MethylMix_ModelSingleGene(GeneNames[i],METcancer[i,],METnormal[i,],NoNormalMode,maxComp,UseChisquare)
        MethylationStates[i,]=MixtureModelResults_SingleGene$MethylationState
        Models[[i]]=MixtureModelResults_SingleGene$Model
        AllNrComponents[i,1]=MixtureModelResults_SingleGene$NrComponents
        AllFlipOverStates[i,1]=MixtureModelResults_SingleGene$FlipOverState
        AllMixtureStates[[i]]=MixtureModelResults_SingleGene$MixtureStates
    } 
    MethylationDrivers = rownames(METcancer)
    
    return(list(MethylationStates=MethylationStates,NrComponents=AllNrComponents,Models=Models,MethylationDrivers=MethylationDrivers,MixtureStates=AllMixtureStates))
    
}

# SplineFilter <- function(METcancer, MAcancer) {
#     
#     # overlapping samples     
#     OverlapSamples = intersect(colnames(METcancer), colnames(MAcancer))
#     cat("Found", length(OverlapSamples), "samples with both methylation and expression data.\n")
#     MAcancer = MAcancer[, OverlapSamples, drop = FALSE]
#     METcancer = METcancer[, OverlapSamples, drop = FALSE]
#     
#     pvalues = matrix(1, nrow = length(rownames(METcancer)), ncol = 1)    
#     Genes = rownames(METcancer)  
#     PvalueThreshold = 0.001  
#     
#     pb = txtProgressBar(1, length(rownames(METcancer)), style = 3, width = 100)
#     cat("Analyzing methylation and gene expression data.\n")
#     for(i in 1:length(rownames(METcancer))) {
#         setTxtProgressBar(pb, i)
#         tmpGene = unlist(strsplit(Genes[i], '---'))[1]
#         pos = which(rownames(MAcancer) == tmpGene)
#         
#         if (length(pos) > 0) {
#             res = glm(MAcancer[pos, ] ~ ns(METcancer[Genes[i], ], df = 4))
#             pvalues[i] = anova(res,test = "Chisq")$Pr[2]
#         }
#     }
#     
#     # Rsquare threshold
#     FunctionalGenes = Genes[pvalues <= PvalueThreshold]
#     cat("\nFound", length(FunctionalGenes), "functional genes.\n")
#     return(FunctionalGenes)
#     
# }

getMvalues <- function(matr) {   
    # Handle values out of range
    for (i in 1:nrow(matr)) {
        xmax = max(matr[i, ], na.rm = TRUE)
        xmin = min(matr[i, ], na.rm = TRUE)
        if (xmax >= 1) { # check this so we only run this where it's needed
            row = matr[i, ]
            max2 = max(row[row < 1], na.rm = TRUE)
            matr[i, ] = pmin(row, 1 - (1 - max2)/2)
        }
        if (xmin <= 0) {
            row = matr[i, ]
            min2 = min(row[row > 0], na.rm = TRUE)
            matr[i, ] = pmax(matr[i, ], min2/2)
        }
    }
    # Transform inot M values
    matr <- log2(matr/(1 - matr))
    return (matr)
}

Plot_Mvalues_Univ <-function(GeneName, MixtureModelResults, METdata, MAdata=0, METnormal=0, FileName="") {
    
    GeneNameMA = unlist(strsplit(GeneName, '---'))[1]
    
    Pos = which(rownames(MixtureModelResults$MethylationStates) %in% GeneName)
    Pos2 = which(rownames(METdata) %in% GeneName)
    if (length(Pos) > 0 & length(Pos2) > 0) {
        NormalModel = MixtureModelResults$Parameters[[Pos]]
        NrComponents = MixtureModelResults$NrComponents[Pos]
        METdataVector = getMvalues(METdata[Pos2, , drop = F])
        
        if (FileName != "") {
            File = paste(FileName, "MixtureModel.tif", sep = '')
            tiff(filename = File, compression = "none", pointsize = 20, bg = "white", width = 600, height = 600)                     
        }
        histogram = hist(METdataVector, 50, plot=FALSE)
        HistHeight = max(histogram$counts) * 1.2
        histogram = hist(METdataVector, 50, ylim = c(0, 1.4 * HistHeight), main = paste("Gaussian Mixture Model of", GeneName), xlab = "DNA methylation (M-value)", ylab = "Frequency")
        
        # plotting the normal 95% confidence interval, if the normal data is present
        if (length(METnormal) > 1) {
            METnormalVector = getMvalues(METnormal[GeneName, , drop = F])
            tmpTtest = t.test(METnormalVector)            
            NormalMean = mean(METnormalVector)
            lines(c(tmpTtest$conf.int[1], tmpTtest$conf.int[2]), c(1.2 * HistHeight, 1.2 * HistHeight), ylim = c(0, 1.2 * HistHeight), lwd=5, type='l', col=1)
            points(NormalMean, 1.2 * HistHeight, ylim = c(0, 1.2 * HistHeight), lwd=5, type = 'l', col = 1, pch = 8)
        }   
        
        x = seq(min(METdataVector, na.rm = T) - 0.5, max(METdataVector, na.rm = T) + 0.5, 0.01) 
        MaxHeight = numeric(NrComponents)
        for (comp in 1:NrComponents) {
            MaxHeight[comp] = max(dnorm(x, NormalModel$mean[comp], NormalModel$sd[comp]))
        }
        MaxHeight = max(MaxHeight)
        Factor = max(HistHeight/MaxHeight)
        for (comp in 1:NrComponents) {
            lines(x, Factor * NormalModel$prop[comp] * dnorm(x, NormalModel$mean[comp], NormalModel$sd[comp]), lwd = 5, col = comp + 1)
        }
        
        #plotting the Gene expression values
        if (length(MAdata) > 1) {
            OverlapSamples = intersect(colnames(METdata),colnames(MAdata))
            METdataVector = METdataVector[, OverlapSamples, drop=FALSE]
            MAdataVector = MAdata[GeneNameMA, OverlapSamples, drop=FALSE]
            if (FileName != "") {
                dev.off()
                File = paste(FileName, "Expression.tif", sep = '')
                tiff(filename = File, compression="none", pointsize = 20, bg = "white", width = 600, height = 600)                     
            }
            
            # The number of components and the mixture states may not be the same
            # For ex, we can have 3 components, but two of them were identified with a 0 DM value
            # If so, we plot different colors for the components and different symbols for the states
            # Vector with the classification into components of each sample
            ClassVector = MixtureModelResults$Classification[[GeneName]][OverlapSamples]        
            
            plot(METdataVector, MAdataVector, type = 'n', xlab = "DNA methylation (M-values)", ylab = "Gene expression", main = paste("Expression correlation for", GeneName))
            
            for (comp in 1:NrComponents) {
                if (length(unique(MixtureModelResults$MixtureStates[[GeneName]])) == MixtureModelResults$NrComponents[GeneName,]) {
                    points(METdataVector[, ClassVector == comp], MAdataVector[, ClassVector == comp], col = comp + 1, pch = 19)
                } else {
                    # Vector of unique DM values (unique Methylation States)
                    uniqueStates = unique(MixtureModelResults$MixtureStates[[GeneName]])
                    # Vector indicating to which mixture states belong each sample
                    stateIndicator = sapply(MixtureModelResults$MethylationStates[GeneName, ], function(x) which(uniqueStates == x))[OverlapSamples]
                    points(METdataVector[, ClassVector == comp], MAdataVector[, ClassVector == comp], col = comp + 1, pch = stateIndicator[ClassVector == comp])                    
                }
            }
        }
    } else {
        cat("This gene does not exist.\n") 
    }
    if (FileName!="") {
        dev.off()
    }
}

Plot_Mvalues_Biv <-function(GeneName, MixtureModelResults, METdata, MAdata, METnormal=0, MAnormal=0, FileName="") {
    
    # Load package ellipse
    if(require("ellipse")){
        cat("ellipse is loaded correctly\n")
    } else {
        cat("trying to install ellipse\n")
        install.packages("ellipse")
        if(require(ellipse)){
            cat("ellipse installed and loaded\n")
        } else {
            stop("could not install ellipse")
        }
    }
    
    GeneNameMA = unlist(strsplit(GeneName, '---'))[1]
    
    Pos = which(rownames(MixtureModelResults$MethylationStates) %in% GeneName)
    Pos2 = which(rownames(METdata) %in% GeneName)
    Pos3 = which(rownames(MAdata) %in% GeneNameMA)
    if (length(Pos) > 0 & length(Pos2) > 0 & length(Pos3) > 0) {
        
        NormalModel = MixtureModelResults$Parameters[[Pos]]
        NrComponents = MixtureModelResults$NrComponents[Pos]
        
        OverlapSamples = intersect(colnames(METdata),colnames(MAdata))
        METdataVector = getMvalues(METdata[Pos2, OverlapSamples, drop = F])
        MAdataVector = MAdata[GeneNameMA, OverlapSamples, drop = F]
        
        if (FileName != "") {
            File = paste(FileName, "MixtureModel.tif", sep = '')
            tiff(filename = File, compression = "none", pointsize = 20, bg = "white", width = 600, height = 600)                     
        }
        
        plot(METdataVector, MAdataVector, type = 'n', xlab = "DNA methylation (M-values)", ylab = "Gene expression", main = paste("Gene: ", GeneName))
        
        ClassVector = MixtureModelResults$Classification[[GeneName]][OverlapSamples]        
        
        for (comp in 1:NrComponents) {
            points(METdataVector[,  ClassVector == comp], MAdataVector[, ClassVector == comp], col = comp + 1, pch = 19, cex = 0.8)            
            sigma = MixtureModelResults$Parameters[[GeneName]]$sigma[[comp]]
            means = MixtureModelResults$Parameters[[GeneName]]$mean[, comp]
            points(ellipse(sigma, centre = means, level = 0.90), type = 'l', col = comp + 1, lwd = 3)
        }
        
        # Add ellipse for normal data
        if (length(METnormal) > 1 & length(MAnormal) > 1) {
            OverlapSamples = intersect(colnames(METnormal),colnames(MAnormal))
            # only add ellipse for samples with both MET and MA data
            if (length(OverlapSamples) > 1) {
                METnormal2 = getMvalues(METnormal[GeneName, OverlapSamples, drop=FALSE])
                MAnormal2 = MAnormal[GeneNameMA, OverlapSamples, drop=FALSE]
                data <- t(rbind(METnormal2, MAnormal2))
                points(ellipse(cov(data), centre = colMeans(data, na.rm = T), level = 0.9), type = 'l', lwd = 3)
            }
        }
        # Add line with confidence interval for methylation with normal data
        if (length(METnormal) > 1) {
            tmpTtest = t.test(getMvalues(METnormal[GeneName, , drop = F]))            
            NormalMean = mean(getMvalues(METnormal[GeneName, , drop = F]))
            height = par('usr')[3] + (par('usr')[4] - par('usr')[3]) * 0.08
            lines(c(tmpTtest$conf.int[1], tmpTtest$conf.int[2]), c(height, height), lwd=5)
        }
        # Add line with confidence interval for gene expression with normal data
        if (length(MAnormal) > 1) {
            tmpTtest = t.test(MAnormal[GeneNameMA, ])            
            NormalMean = mean(MAnormal[GeneNameMA, ])
            height = par('usr')[1] + (par('usr')[2] - par('usr')[1]) * 0.05
            lines(c(height, height), c(tmpTtest$conf.int[1], tmpTtest$conf.int[2]), lwd=5)
        }
        
    } else {
        cat("This gene does not exist.\n") 
    }
    if (FileName!="") {
        dev.off()
    }
}

#GeneName = "A2ML1---Cluster1"; MixtureModelResults = luad.orig; METdata = METcancer; MAdata = MAcancer
Plot_Beta_Univ <-function(GeneName, MixtureModelResults, METdata, MAdata=0, METnormal=0, FileName="") {
    
    GeneNameMA = unlist(strsplit(GeneName, '---'))[1]
    
    Pos = which(rownames(MixtureModelResults$MethylationStates) %in% GeneName)
    Pos2 = which(rownames(METdata) %in% GeneName)
    if (length(Pos) > 0 & length(Pos2) > 0) {
        BetaModel = MixtureModelResults$Models[[Pos]]
        NrComponents = MixtureModelResults$NrComponents[Pos]
        METdataVector=METdata[Pos2,]
        
        if (FileName != "") {
            File = paste(FileName, "MixtureModel.tif", sep = '')
            tiff(filename = File, compression = "none", pointsize = 20, bg = "white", width = 600, height = 600)                     
        }
        histogram=hist(METdataVector,50,plot=FALSE)
        HistHeight=max(histogram$counts)*1.2
        histogram=hist(METdataVector,50,xlim=c(0,1),ylim=c(0,1.4*HistHeight),main=paste("Mixture model of",GeneName),xlab="DNA methylation",ylab="Frequency")
        
        # plotting the normal 95% confidence interval, if the normal data is present
        if (length(METnormal) > 1) {
            tmpTtest=t.test(METnormal[GeneName,])            
            NormalMean=mean(METnormal[GeneName,])
            lines(c(tmpTtest$conf.int[1],tmpTtest$conf.int[2]),c(1.2*HistHeight,1.2*HistHeight),ylim=c(0,1.2*HistHeight),lwd=5,type='l',col=1)
            points(NormalMean,1.2*HistHeight,ylim=c(0,1.2*HistHeight),lwd=5,type='l',col=1,pch=8)
        }   
        
        x = seq(0,1,0.001)
        MaxHeight = numeric(NrComponents)
        for (comp in 1:NrComponents) {
            MaxHeight[comp] = max(dbeta(x, BetaModel$a[comp], BetaModel$b[comp]))
        }
        MaxHeight = max(MaxHeight)
        Factor = max(HistHeight/MaxHeight)
        for (comp in 1:NrComponents) {
            lines(x, Factor * BetaModel$eta[comp] * dbeta(x, BetaModel$a[comp], BetaModel$b[comp]), ylim = c(0, 1.2 * HistHeight), lwd = 5, type = 'l', col = comp + 1, xlab = "DNA methylation", ylab = "Frequency")  
        }
        
        #plotting the Gene expression values
        if (length(MAdata) > 1) {
            OverlapSamples = intersect(colnames(METdata),colnames(MAdata))
            METdataVector = METdata[GeneName, OverlapSamples, drop=FALSE]
            MAdataVector = MAdata[GeneNameMA, OverlapSamples, drop=FALSE]
            if (FileName != "") {
                dev.off()
                File = paste(FileName, "Expression.tif", sep = '')
                tiff(filename = File, compression="none", pointsize = 20, bg = "white", width = 600, height = 600)                     
            }
            
            # The number of components and the mixture states may not be the same
            # For ex, we can have 3 components, but two of them were identified with a 0 DM value
            # If so, we plot different colors for the components and different symbols for the states
            # Vector with the classification into components of each sample
            # ClassVector = MixtureModelResults$Classification[[GeneName]][OverlapSamples]
            names(MixtureModelResults$Models) = MixtureModelResults$MethylationDrivers
            ClassVector <- apply(MixtureModelResults$Models[[GeneName]]$w, 1, which.max)
            names(ClassVector) = colnames(MixtureModelResults$MethylationStates)
            ClassVector <- ClassVector[OverlapSamples]
            # Vector of unique DM values (unique Methylation States)
            names(MixtureModelResults$MixtureStates) = MixtureModelResults$MethylationDrivers
            uniqueStates = unique(MixtureModelResults$MixtureStates[[GeneName]])
            # Vector indicating to which mixture states belong each sample
            stateIndicator = sapply(MixtureModelResults$MethylationStates[GeneName, ], function(x) which(uniqueStates == x))[OverlapSamples]
            
            plot(METdataVector, MAdataVector, type = 'n', xlim=c(0,1),col=2,lty=1, xlab = "DNA methylation", ylab = "Gene expression", main = paste("Expression correlation for", GeneName))
            
            for (comp in 1:NrComponents) {
                if (length(unique(MixtureModelResults$MixtureStates[[GeneName]])) == MixtureModelResults$NrComponents[GeneName,]) {
                    points(METdataVector[, ClassVector == comp], MAdataVector[, ClassVector == comp], col = comp + 1, pch = 19)
                } else {
                    # cat(stateIndicator[ClassVector == comp], "\n")        
                    points(METdataVector[, ClassVector == comp], MAdataVector[, ClassVector == comp], col = comp + 1, pch = stateIndicator[ClassVector == comp])                    
                }
            }
        }
    } else {
        cat("This gene does not exist.\n") 
    }
    if (FileName!="") {
        dev.off()
    }
}





MethylMix_ModelGeneExpression_keep_Rsquares <- function(METcancer, MAcancer, Method = c("Regression", "Pearson", "Spearman"), CovariateData = NULL) {      
    
    # overlapping samples     
    OverlapSamples = intersect(colnames(METcancer), colnames(MAcancer))
    cat("Found", length(OverlapSamples), "samples with both methylation and expression data.\n")
    MAcancer = MAcancer[, OverlapSamples, drop = FALSE]
    METcancer = METcancer[, OverlapSamples, drop = FALSE]
    if (!is.null(CovariateData)) CovariateData = as.matrix(CovariateData[OverlapSamples, ])
    
    Rsquares = matrix(0, nrow = length(rownames(METcancer)), ncol = 1)    
    Genes = rownames(METcancer)  
    PvalueThreshold = 0.001  
    RsquareThreshold = 0.1
    
    pb = txtProgressBar(1, length(rownames(METcancer)), style = 3, width = 100)
    cat("Correlating methylation data with gene expression.\n")
    for(i in 1:length(rownames(METcancer))) {
        setTxtProgressBar(pb, i)
        tmpGene = unlist(strsplit(Genes[i], '---'))[1]
        pos = which(rownames(MAcancer) == tmpGene)
        
        if (length(pos) > 0) {
            if (!is.null(CovariateData)) {
                res = lm(MAcancer[pos, ] ~ METcancer[Genes[i], ] + CovariateData)
            } else {           
                res = lm(MAcancer[pos, ] ~ METcancer[Genes[i], ])
            }
            res.summary = summary(res)
            if (!is.null(CovariateData)) {
                if (res$coefficients[2] < 0 & res.summary$coefficients[2, 4] < PvalueThreshold & abs(res.summary$coefficients[2, 1]) > abs(res.summary$coefficients[3, 1])) { # methylation effect bigger than tissue
                    Rsquares[i] = res.summary$r.squared
                }
            } else if (res$coefficients[2] < 0 & res.summary$coefficients[2, 4] < PvalueThreshold) {
                Rsquares[i] = res.summary$r.squared    
            }
        }
    }
    
    # Rsquare threshold
    FunctionalGenes = Genes[Rsquares > RsquareThreshold]
    Rsquares=Rsquares[Rsquares > RsquareThreshold]
    names(Rsquares)=FunctionalGenes
    Rsquares=as.data.frame(Rsquares,drop=FALSE)
    cat("\nFound", length(FunctionalGenes), "functional genes.\n")
    return(list(FunctionalGenes,Rsquares))
}

##########################################################################################
#########################################################################################
#Scripts to download and preprocess TCGA DNA methylatio data 
#######################################################################################
#############################################################################################

library(RCurl)
library(limma)

Download_CancerSite_MethylMix <- function(CancerSite,TargetDirectory,downloadData=TRUE) {    
	cat('Downloading MA and MET data for:',CancerSite,'\n')
     
	# downloading the RNAseq data data
	MAdirectories=c("NA")
     MAdirectories=Download_CancerSite_GeneExpression(CancerSite,TargetDirectory,downloadData)    
		
     # downloading the methylation data
	METdirectories=c("NA","NA")
	METdirectories=Download_CancerSite_DNAmethylation(Cancer,TargetDirectory,downloadData)
	
	return(list(MAdirectories=MAdirectories,METdirectories=METdirectories))     
}

Download_CancerSite_GeneExpression <- function(CancerSite,TargetDirectory,downloadData=TRUE) {    
     #system(command)
     dir.create(TargetDirectory,showWarnings=FALSE)
     
     # Settings
     TCGA_acronym_uppercase=toupper(CancerSite)
     
     # get RNA seq data (GBM does not have much RNAseq data.)
     dataType='stddata'	
     dataFileTag='mRNAseq_Preprocess.Level_3'	 
     
     #special case for GBM and OV, not enough RNAseq data, so using the microarray data instead
     if (CancerSite=="GBM") { 	             
          dataFileTag=c('Merge_transcriptome__agilentg4502a_07_1__unc_edu__Level_3__unc_lowess_normalization_gene_level__data','Merge_transcriptome__agilentg4502a_07_2__unc_edu__Level_3__unc_lowess_normalization_gene_level__data')        	         
     } else if(CancerSite=="OV") {	               
          dataFileTag='Merge_transcriptome__agilentg4502a_07_3__unc_edu__Level_3__unc_lowess_normalization_gene_level__data'        
     }	
     cat('Searching MA data for:',CancerSite,"\n")
     if (length(dataFileTag)==1) {	  
          MAdirectories=get_firehoseData(downloadData,saveDir=TargetDirectory,TCGA_acronym_uppercase=TCGA_acronym_uppercase,dataFileTag=dataFileTag)    	
     } else {	    # a few data sets have multiple gene expression data sets. 
          MAdirectories=c()	  
          for (i in 1:length(dataFileTag)) {
               MAdirectories=c(MAdirectories,get_firehoseData(downloadData,saveDir=TargetDirectory,TCGA_acronym_uppercase=TCGA_acronym_uppercase,dataFileTag=dataFileTag[i]))	 
          }        
     }
     
     return(MAdirectories=MAdirectories)
}


Download_CancerSite_DNAmethylation_KB <- function(CancerSite,TargetDirectory,downloadData=TRUE, tag27kor450k) {    
     
     # need to redo this, because calling this function also independently of Download_CancerSite_MethylMix
     dir.create(TargetDirectory,showWarnings=FALSE)
     
     # download the 27k data
     #dataType='stddata'
     #dataFileTag='Merge_methylation__humanmethylation27'
     #cat('Searching 27k MET data for:',CancerSite,'\n')
     #METdirectory27k=get_firehoseData(downloadData,TargetDirectory,CancerSite,dataType,dataFileTag)
    
     # download the 450k data
     dataFileTag=paste('Merge_methylation__humanmethylation', tag27kor450k, sep="")
     cat('Searching 450k MET data for:',CancerSite,'\n')
     METdirectory=get_firehoseData(downloadData,TargetDirectory,CancerSite,dataType,dataFileTag)

     return(METdirectories=list(METdirectory=METdirectory))
}


#Marcos' version

Download_CancerSite_DNAmethylation <- function(CancerSite,TargetDirectory,downloadData=TRUE) {    
  
  # need to redo this, because calling this function also independently of Download_CancerSite_MethylMix
  dir.create(TargetDirectory,showWarnings=FALSE)
  
  # download the 27k data
  dataType='stddata'
  dataFileTag='Merge_methylation__humanmethylation27'
  cat('Searching 27k MET data for:',CancerSite,'\n')
  METdirectory27k=get_firehoseData(downloadData,TargetDirectory,CancerSite,dataType,dataFileTag)
  
  # download the 450k data
  dataFileTag='Merge_methylation__humanmethylation450'
  cat('Searching 450k MET data for:',CancerSite,'\n')
  METdirectory450k=get_firehoseData(downloadData,TargetDirectory,CancerSite,dataType,dataFileTag)
  
  
  
  return(METdirectories=list(METdirectory27k=METdirectory27k,METdirectory450k=METdirectory450k))
}



get_firehoseData <- function(downloadData=TRUE,saveDir = "./",TCGA_acronym_uppercase = "LUAD",dataType="stddata",dataFileTag = "mRNAseq_Preprocess.Level_3",
                             FFPE=FALSE,fileType= "tar.gz",  gdacURL= "http://gdac.broadinstitute.org/runs/",untarUngzip=TRUE,printDisease_abbr=FALSE){  
     
     # Cases Shipped by BCR  # Cases with Data*  Date Last Updated (mm/dd/yy)
     cancers <- c("Acute Myeloid Leukemia [LAML] \n","Adrenocortical carcinoma [ACC]	\n",
                  "Bladder Urothelial Carcinoma [BLCA] \n",	"Brain Lower Grade Glioma [LGG] \n",
                  "Breast invasive carcinoma [BRCA] \n","Cervical squamous cell carcinoma and endocervical adenocarcinoma [CESC] \n",
                  "Cholangiocarcinoma [CHOL] \n",	"Colon adenocarcinoma [COAD] \n",	"Esophageal carcinoma [ESCA] \n",
                  "Glioblastoma multiforme [GBM] \n",	"Head and Neck squamous cell carcinoma [HNSC]	\n",
                  "Kidney Chromophobe [KICH]	\n","Kidney renal clear cell carcinoma [KIRC]	\n",
                  "Kidney renal papillary cell carcinoma [KIRP]	\n","Liver hepatocellular carcinoma [LIHC]	\n",
                  "Lung adenocarcinoma [LUAD]	\n", "Lung squamous cell carcinoma [LUSC] \n",
                  "Lymphoid Neoplasm Diffuse Large B-cell Lymphoma [DLBC]	\n","Mesothelioma [MESO] \n",
                  "Ovarian serous cystadenocarcinoma [OV]	\n","Pancreatic adenocarcinoma [PAAD]	\n",
                  "Pheochromocytoma and Paraganglioma [PCPG] \n","Prostate adenocarcinoma [PRAD] \n",
                  "Rectum adenocarcinoma [READ]	\n","Sarcoma [SARC]	\n","Skin Cutaneous Melanoma [SKCM]	\n",
                  "Stomach adenocarcinoma [STAD] \n","Testicular Germ Cell Tumors [TGCT] \n","Thymoma [THYM] \n",
                  "Thyroid carcinoma [THCA]	\n","Uterine Carcinosarcoma [UCS]	 \n",
                  "Uterine Corpus Endometrial Carcinoma [UCEC]	\n","Uveal Melanoma [UVM] \n");
     
     if(printDisease_abbr){      
          return(cat("here are the possible TCGA database disease acronyms. \nRe-run this function with printDisease_abbr=FALSE to then run an actual query.\n\n",cancers));      
     }
     gdacURL_orig <- gdacURL
     urlData <- getURL(gdacURL)
     urlData <- strsplit2(urlData,paste(dataType,"__",sep=""))
     #This is where it screws up
     #The 'latest' data is not the latest, according ot the dates
     #All of the text is in the first columns, the dates are 13 characters long, so I just choose these to get into the next folder
     
     urlData=urlData[nchar(urlData)==13]
     urlData=gsub("/.*","", urlData)
     
     #first column is junk
     #urlData <- urlData[,2:dim(urlData)[2]]
     #urlData <- strsplit2(urlData,"/")
     #urlData <- urlData[,1]
     #This is failing because the structure of the page has changed I think
     
     #see the strptime codes here: http://stat.ethz.ch/R-manual/R-devel/library/base/html/strptime.html
     #I am transferring this text to dates: that way, we can programmatically find the latest one.
     #the POSIXct class in R is the date-time class
     urlData <- as.POSIXct(strptime(urlData, "%Y_%m_%d"))
     #want to remove time zone: do as.Date
     dateData <- as.Date(urlData[which(!is.na(urlData))])
     lastDate <- dateData[match( summary(dateData)[which(names(summary(dateData))=="Max.")], dateData)]
     
     #sub back _ symbols
     lastDate <- gsub("-","_",as.character(lastDate))
     
     lastDateCompress <- gsub("_","",lastDate)
     #need last "/" for it to find this page.
     
      gdacURL <- paste(gdacURL,dataType,"__",lastDate,"/data/",TCGA_acronym_uppercase,"/",lastDateCompress,"/",sep="")
     
     #now get full dataset name. we just have the http link to the last page we select this dataset from right now.
     urlData <- getURL(gdacURL)
     
     
     #regular expressions: need \ to have R recognize any " or \ that's actually in our text
     urlData <- strsplit2(urlData,"href=\\\"")
     
     while(length(grep("was not found",urlData))>0) {
          warning(paste0("\nThe TCGA run dated ",lastDate, "for ", dataType,
                         " for disease ",TCGA_acronym_uppercase," isn't available for download yet.
                         Taking the run dated just before this one.\n"))
          dateData <-  dateData[-which(dateData==(summary(dateData)[which(names(summary(dateData))=="Max.")]))]
          lastDate <- dateData[match( summary(dateData)[which(names(summary(dateData))=="Max.")], dateData)]
          #sub back _ symbols
          lastDate <- gsub("-","_",as.character(lastDate))
          lastDateCompress <- gsub("_","",lastDate)
          #need last "/" for it to find this page.
          gdacURL <- paste(gdacURL_orig,dataType,"__",lastDate,"/data/",TCGA_acronym_uppercase,"/",lastDateCompress,"/",sep="")
          
          #now get full dataset name. we just have the http link to the last page we select this dataset from right now.
          urlData <- getURL(gdacURL)
          #regular expressions: need \ to have R recognize any " or \ that's actually in our text
          urlData <- strsplit2(urlData,"href=\\\"")
          
          #did we reach the end of the dates - ie only one left? leave the loop then.
          if(length(dateData)<=1){		
               break	
          }  
     } 
     
     #cat("Using data from date ",lastDate,"\n")
     #should be OK now!
     if(length(grep("was not found",urlData))>0){  
          #this disease may not even be in the analyses directory yet.
          stop( paste0("\nNot returning any viable url data paths after searching by date for disease ",TCGA_acronym_uppercase,". No data was downloaded.\n"))
     }
     
     #remove any FFPE datasets, or only keep those depending on user inputs.
     if (FFPE) { 
          urlData <- urlData[grep("FFPE",urlData)]	
          if(length(urlData)==0){		
               stop("\nNo FFPE data found for this query. Try FFPE=FALSE.\n")		
          }	  
     } else {	
          #we DON'T want FFPE data.
          #but if no FFPE data to begin with: don't subset on this.
          if(length(grep("FFPE",urlData))>0){		
               urlData <- urlData[-grep("FFPE",urlData)]		
          }
          if(length(urlData)==0){		
               stop("\nNo non-FFPE data found for this query. Try FFPE=TRUE.\n")		
          }
     }
     #now get full dataset name.
     fileName <- urlData[grep(dataFileTag,urlData)]
     
     if(length(fileName)==0){	  
          #warnMessage <- paste0("\nNot returning any viable url data paths after searching by date for disease ",TCGA_acronym_uppercase," for data type ",dataFileTag ,".No data was downloaded.\n")
          #warning(warnMessage)
          cat("\tThere is no",dataFileTag,"data for",TCGA_acronym_uppercase,"\n")
          return(NA)	  
     }
     #some redundancy..but that' OK because we'll add back on the unique tar.gz file tag.
     #first file is one we want - not md5 file.
     fileName <- strsplit2(fileName,"tar.gz")[1,1]
     fileName <- paste(fileName,fileType,sep="")
     
     #final download url
     gdacURL <- paste(gdacURL,fileName,sep="")
     
     # Keep the savedir when we don't download !!!!!!!!
     saveDir <- paste(saveDir,"gdac_",lastDateCompress,'/',sep="")
     
     if(downloadData){		
          cat("\tDownloading",dataFileTag,"data, version:",lastDate,"\n")				
          cat("\tThis may take 10-60 minutes depending on the size of the data set.\n")
          
          # create dirs
          dir.create(saveDir,showWarnings=FALSE)
          
          # download file		
          setwd(saveDir)				
          download.file(gdacURL,fileName,quiet=FALSE,mode="wb")
          
          #this assumes a tar.gz file.
          if(fileType=="tar.gz" && untarUngzip) {					
               cat("\tUnpacking data.\n")
               tarfile=paste0(saveDir,fileName)
               untar(tarfile)
               
               #remove tarred file
               fileToRemove <- strsplit2(gdacURL,"/")[ ,ncol(strsplit2(gdacURL,"/"))]
               file.remove(paste0(saveDir,fileToRemove))
               
          } else if(untarUngzip) {		
               warning("File expansion/opening only built in for tar.gz files at the moment.\n")		
          }
          
          finalDir <- strsplit2(gdacURL,"/")[ ,ncol(strsplit2(gdacURL,"/"))]
          finalDir <- strsplit2(finalDir,fileType)		
          #must remove LAST period (ie laster character) only. 
          finalDir <- substr(finalDir,start=0,stop=(nchar(finalDir)-1))
          finalDir <- paste0(saveDir,finalDir)		
          cat("\tFinished downloading",dataFileTag,"data to",finalDir,"\n")
          
     } else {
          #just spit out the command you need
          cat("\tdownload data url is :\n ",gdacURL,'\n')
          finalDir <- strsplit2(gdacURL,"/")[ ,ncol(strsplit2(gdacURL,"/"))]
          finalDir <- strsplit2(finalDir,fileType)
          
          #must remove LAST period (ie laster character) only. 
          finalDir <- substr(finalDir,start=0,stop=(nchar(finalDir)-1))
          finalDir <- paste0(saveDir,finalDir)	
     }
     DownloadedFile=paste0(finalDir,'/')
     return(DownloadedFile)
     
}

SampleCount_CancerSite <- function(CancerSite) {
          
     METdirectories=Download_CancerSite_DNAmethylation(CancerSite,TargetDirectory,FALSE)
     cat("\tGetting sample count for",CancerSite,"\n")
   
     SamplegroupNumbers27k=c()
     SamplegroupNumbers450k=c()
     if (!is.na(METdirectories$METdirectory27k)) {
          # Load data
          METfiles=dir(METdirectories$METdirectory27k)
          MatchedFilePosition=grep('methylation__humanmethylation27',METfiles)             
          Filename=paste0(METdirectories$METdirectory27k,METfiles[MatchedFilePosition])          
          MET_Data=read.csv(Filename,sep="\t",nrows=10)          
          MET_Data=as.matrix(MET_Data)
          Probes=MET_Data[,1]
          rownames(MET_Data)=Probes
          MET_Data=MET_Data[,-1]
          MET_Data=MET_Data[-1,]
          MET_Data=MET_Data[,seq(1,ncol(MET_Data),4)]     
          class(MET_Data)='numeric'
          Samplegroups27k=TCGA_GENERIC_GetSampleGroups(colnames(MET_Data))
          SamplegroupNumbers27k=list()
          SamplegroupNumbers27k$Primary=length(Samplegroups27k$Primary)
          SamplegroupNumbers27k$SolidNormal=length(Samplegroups27k$SolidNormal)
          cat('The 27k data set has',SamplegroupNumbers27k$Primary,'primary cancer samples and',SamplegroupNumbers27k$SolidNormal,'normal samples.\n')
     }
     
     if (!is.na(METdirectories$METdirectory450k)) {
          # Load data
          METfiles=dir(METdirectories$METdirectory450k)
          MatchedFilePosition=grep('methylation__humanmethylation450',METfiles)             
          Filename=paste0(METdirectories$METdirectory450k,METfiles[MatchedFilePosition])          
          
          MET_Data=read.csv(Filename,sep="\t",nrows=10)
          MET_Data=as.matrix(MET_Data)
          Probes=MET_Data[,1]
          rownames(MET_Data)=Probes
          MET_Data=MET_Data[,-1]
          MET_Data=MET_Data[-1,]
          MET_Data=MET_Data[,seq(1,ncol(MET_Data),4)]     
          class(MET_Data)='numeric'
          
          Samplegroups450k=TCGA_GENERIC_GetSampleGroups(colnames(MET_Data))
          SamplegroupNumbers450k=list()
          SamplegroupNumbers450k$Primary=length(Samplegroups450k$Primary)
          SamplegroupNumbers450k$SolidNormal=length(Samplegroups450k$SolidNormal)
          #SamplegroupNumbers450k$BloodNormal=length(Samplegroups450k$Primary)
          cat('The 450k data set has',SamplegroupNumbers450k$Primary,'primary cancer samples and',SamplegroupNumbers450k$SolidNormal,'normal samples.\n')
     }
     
     return(list(SamplegroupNumbers27k=SamplegroupNumbers27k,SamplegroupNumbers450k=SamplegroupNumbers450k))
}

Preprocess_CancerSite_GeneExpression <- function(CancerSite,MAdirectories) {    
     #below dataset name is: BatchData once loaded.
  dat=load("/home/kbren/rscripts/BatchData.rda")
  MinPerBatchCancer=5    
     MinPerBatchNormal=2    
     
     # Processing MA data, special case for OV and GBM where no RNA seq data is available
     if (CancerSite=="OV" || CancerSite=="GBM") { 
          MAstring='transcriptome__agilent'
     } else if (CancerSite=="STAD" || CancerSite=="ESCA") { # for these cancers RSEM data does not exist. 
          MAstring='mRNAseq_RPKM_log2.txt'
     } else {
          MAstring='mRNAseq_RSEM_normalized_log2.txt'
     }
     
     #cat("Loading mRNA data.\n")
     if (length(MAdirectories)>1) {      		
          cat("\tFound multiple MA data sets.\n")
          DataSetsCancer=list()
          GeneListsCancer=list()
          SampleListsCancer=list()                
          MetaBatchDataCancer=data.frame()
          
          DataSetsNormal=list()
          GeneListsNormal=list()
          SampleListsNormal=list()                
          MetaBatchDataNormal=data.frame()
          for (i in 1:length(MAdirectories)) {        
               cat("\tProcessing data set",i,"\n")
               MAfiles=dir(MAdirectories[i])
               MatchedFile=grep(MAstring,MAfiles)        
               if (length(MatchedFile)>0) {        
                    # Getting the cancer data first
                    DataSetsCancer[[i]]=Preprocess_MAdata_Cancer(CancerSite,MAdirectories[i],MAfiles[MatchedFile])
                    GeneListsCancer[[i]]=rownames(DataSetsCancer[[i]])
                    SampleListsCancer[[i]]=colnames(DataSetsCancer[[i]])                    
                    currentBatchCancer=matrix(i,length(colnames(DataSetsCancer[[i]])),1) # growing a batch data object
                    currentBatchDataCancer=data.frame(ArrayName=colnames(DataSetsCancer[[i]]),SampleName=colnames(DataSetsCancer[[i]]),Batch=currentBatchCancer)
                    MetaBatchDataCancer=rbind(MetaBatchDataCancer,currentBatchDataCancer)
                    
                    # Getting the normal data as well.
                    DataSetsNormal[[i]]=Preprocess_MAdata_Normal(CancerSite,MAdirectories[i],MAfiles[MatchedFile])
                    GeneListsNormal[[i]]=rownames(DataSetsNormal[[i]])
                    SampleListsNormal[[i]]=colnames(DataSetsNormal[[i]])                    
                    currentBatchNormal=matrix(i,length(colnames(DataSetsNormal[[i]])),1) # growing a batch data object
                    currentBatchDataNormal=data.frame(ArrayName=colnames(DataSetsNormal[[i]]),SampleName=colnames(DataSetsNormal[[i]]),Batch=currentBatchNormal)
                    MetaBatchDataNormal=rbind(MetaBatchDataNormal,currentBatchDataNormal)
                    
               } else {
                    cat("MA file not found for this cancer.\n")
               }           
          }
          # combine data sets with Combat. 
          cat("Combining data sets.\n")
          OverlapProbesCancer=Reduce(intersect,GeneListsCancer)
          OverlapProbesNormal=Reduce(intersect,GeneListsNormal)
          OverlapSamplesCancer=Reduce(intersect,SampleListsCancer)    
          OverlapSamplesNormal=Reduce(intersect,SampleListsNormal)    
          if (length(OverlapSamplesCancer)>0 | length(OverlapSamplesNormal)>0) {
               cat('This should not happen. There is overlap between cancer or normal samples. No solution yet.\n')           
          }        
          
          for (i in 1:length(MAdirectories)) {
               DataSetsCancer[[i]]=DataSetsCancer[[i]][OverlapProbesCancer,]
               DataSetsNormal[[i]]=DataSetsNormal[[i]][OverlapProbesNormal,]
          }
          # combat on cancer data sets. 
          MA_TCGA_Cancer=Reduce(cbind,DataSetsCancer)
          MA_TCGA_Cancer=TCGA_BatchCorrection_MolecularData(MA_TCGA_Cancer,MetaBatchDataCancer,MinPerBatchCancer)    
          
          # combat on normal data sets. 
          MA_TCGA_Normal=Reduce(cbind,DataSetsNormal)
          MA_TCGA_Normal=TCGA_BatchCorrection_MolecularData(MA_TCGA_Normal,MetaBatchDataNormal,MinPerBatchNormal)    
          
     } else {
          
          MAfiles=dir(MAdirectories)
          MatchedFile=grep(MAstring,MAfiles)        
          if (length(MatchedFile)>0) {                  
               MA_TCGA_Cancer=Preprocess_MAdata_Cancer(CancerSite,MAdirectories,MAfiles[MatchedFile])
               MA_TCGA_Normal=Preprocess_MAdata_Normal(CancerSite,MAdirectories,MAfiles[MatchedFile])               
               cat("There are",length(colnames(MA_TCGA_Cancer)),"cancer samples and",length(colnames(MA_TCGA_Normal)),"normal samples in",CancerSite,"\n")
          } else {               
               stop("MA file not found for this cancer.\n")
          }           
     }         
     return(list(MA_Data_Cancer=MA_TCGA_Cancer,MA_Data_Normal=MA_TCGA_Normal))
}

Preprocess_MAdata_Cancer <- function(CancerSite,Directory,File) {    
    dat=load("/home/kbren/rscripts/BatchData.rda")
    MinPerBatch=5   
     cat("Loading cancer mRNA data.\n")
     cat("\tMissing value estimation.\n")
     MA_TCGA=TCGA_Load_MolecularData(paste(Directory,File,sep=''))        
     Samplegroups=TCGA_GENERIC_GetSampleGroups(colnames(MA_TCGA))        
     if (CancerSite =='LAML') {
          MA_TCGA=MA_TCGA[,Samplegroups$PeripheralBloodCancer,drop=F]
     } else {
          MA_TCGA=MA_TCGA[,Samplegroups$Primary,drop=F]
     }          
     cat("\tBatch correction.\n")
     #Here's where it breaks
     MA_TCGA=TCGA_BatchCorrection_MolecularData(MA_TCGA,BatchData,MinPerBatch)
     
     cat("\tProcessing gene ids and merging.\n")
     Genes=rownames(MA_TCGA)
     SplitGenes=strsplit2(Genes,'\\|')
     rownames(MA_TCGA)=SplitGenes[,1]        
     MA_TCGA=MA_TCGA[!rownames(MA_TCGA) %in% '?',,drop=F]        
     MA_TCGA=TCGA_GENERIC_MergeData(unique(rownames(MA_TCGA)),MA_TCGA)  
     
     return(MA_TCGA=MA_TCGA)
}

Preprocess_MAdata_Normal <- function(CancerSite,Directory,File) {    
     dat=load("/home/kbren/rscripts/BatchData.rda")
    MinPerBatch=2 # less samples in one batch when dealing with normals.
     
     cat("Loading normal mRNA data.\n")
     cat("\tMissing value estimation.\n")
     MA_TCGA=TCGA_Load_MolecularData(paste(Directory,File,sep=''))        
     Samplegroups=TCGA_GENERIC_GetSampleGroups(colnames(MA_TCGA))        
     if (CancerSite =='LAML') {
          MA_TCGA=MA_TCGA[,Samplegroups$BloodNormal,drop=F]
     } else {
          # MA_TCGA=MA_TCGA[,Samplegroups$SolidNormal]
          MA_TCGA=MA_TCGA[,Samplegroups$SolidNormal,drop=F]
     }          
     cat("\tBatch correction.\n")
     MA_TCGA=TCGA_BatchCorrection_MolecularData(MA_TCGA,BatchData,MinPerBatch)
     
     cat("\tProcessing gene ids and merging.\n")
     Genes=rownames(MA_TCGA)
     SplitGenes=strsplit2(Genes,'\\|')
     rownames(MA_TCGA)=SplitGenes[,1]        
     # MA_TCGA=MA_TCGA[!rownames(MA_TCGA) %in% '?',]        
     MA_TCGA=MA_TCGA[!rownames(MA_TCGA) %in% '?', , drop=F]
     MA_TCGA=TCGA_GENERIC_MergeData(unique(rownames(MA_TCGA)),MA_TCGA)  
     
     return(MA_TCGA=MA_TCGA)
}

Preprocess_CancerSite_DNAmethylation <- function(CancerSite,METdirectories) {    
     
     cat("\tProcessing data for",CancerSite,"\n")
     ProcessedData27k=c()
     ProcessedData450k=c()
     if (!is.na(METdirectories$METdirectory27k)) {
          cat('\tLoading data for 27k.\n')
          ProcessedData27k=Preprocess_CancerSite_Methylation27k(CancerSite,METdirectories$METdirectory27k)
     }
     
     if (!is.na(METdirectories$METdirectory450k)) {
          cat('\tLoading data for 450k.\n')
          # No need to use this specific Preprocess function for BRCA
#           if (CancerSite == "BRCA") {
#               ProcessedData450k=Preprocess_CancerSite_Methylation450k_BRCA(CancerSite,METdirectories$METdirectory450k)
#           } else {
              ProcessedData450k=Preprocess_CancerSite_Methylation450k(CancerSite,METdirectories$METdirectory450k)
#           }
     }
     
     # check if we want to combine 27k and 450k
     if (length(ProcessedData27k)!=0 & length(ProcessedData450k)!=0 && CancerSite !="LAML") {
          # only do it when enough samples are in both
          if (ncol(ProcessedData27k$MET_Data_Cancer)>50 & ncol(ProcessedData450k$MET_Data_Cancer)>50) {               
               Mode='450kon27k'               
               cat("\tCombining 450k and 27k by mapping 450k probes to 27k array.\n")
               # Check if there are duplicate samples, remove the 27k ones. 
               OverlapSamplesCancer=intersect(colnames(ProcessedData27k$MET_Data_Cancer),colnames(ProcessedData450k$MET_Data_Cancer))
               if (length(OverlapSamplesCancer)>0) {
                    cat("\tCancer sample overlap is not empty: ",length(OverlapSamplesCancer))
                    ProcessedData27k$MET_Data_Cancer=ProcessedData27k$MET_Data_Cancer[,-OverlapSamplesCancer,drop=FALSE]
               }
               OverlapSamplesNormal=intersect(colnames(ProcessedData27k$MET_Data_Normal),colnames(ProcessedData450k$MET_Data_Normal))
               if (length(OverlapSamplesNormal)>0) {
                    cat("\tNormal sample overlap is not empty: ",length(OverlapSamplesNormal))
                    ProcessedData27k$MET_Data_Normal=ProcessedData27k$MET_Data_Normal[,-OverlapSamplesNormal,drop=FALSE]
               }
               
               # Overlap the probes
               ProcessedData=list(MET_Data_Cancer=c(),MET_Data_Normal=c())
               OverlapProbesCancer=intersect(rownames(ProcessedData27k$MET_Data_Cancer),rownames(ProcessedData450k$MET_Data_Cancer))               
               ProcessedData$MET_Data_Cancer=cbind(ProcessedData27k$MET_Data_Cancer[OverlapProbesCancer,,drop=FALSE],ProcessedData450k$MET_Data_Cancer[OverlapProbesCancer,,drop=FALSE])
               if ( length(colnames(ProcessedData27k$MET_Data_Normal))>0 & length(colnames(ProcessedData450k$MET_Data_Normal))>0 ) {
                    OverlapProbesNormal=intersect(rownames(ProcessedData27k$MET_Data_Normal),rownames(ProcessedData450k$MET_Data_Normal))               
                    ProcessedData$MET_Data_Normal=cbind(ProcessedData27k$MET_Data_Normal[OverlapProbesNormal,,drop=FALSE],ProcessedData450k$MET_Data_Normal[OverlapProbesNormal,,drop=FALSE])
               } else if ( length(colnames(ProcessedData27k$MET_Data_Normal))>0 ) {
                    ProcessedData$MET_Data_Normal=ProcessedData27k$MET_Data_Normal                    
               } else if ( length(colnames(ProcessedData450k$MET_Data_Normal))>0 ) {
                    ProcessedData$MET_Data_Normal=ProcessedData450k$MET_Data_Normal                    
               } 
               
               # Batch correction on combined Tumor data.
               Batch=matrix(1,length(colnames(ProcessedData$MET_Data_Cancer)),1)
               Batch[1:length(colnames(ProcessedData27k$MET_Data_Cancer)),1]=2
               BatchData=data.frame(ArrayName=colnames(ProcessedData$MET_Data_Cancer),SampleName=colnames(ProcessedData$MET_Data_Cancer),Batch=Batch)
               ProcessedData$MET_Data_Cancer=TCGA_BatchCorrection_MolecularData(ProcessedData$MET_Data_Cancer,BatchData,0)
               
               if (length(colnames(ProcessedData$MET_Data_Normal))>0) {
                    # Batch correction on combined Normal data.
                    Batch=matrix(1,length(colnames(ProcessedData$MET_Data_Normal)),1)
                    Batch[1:length(colnames(ProcessedData27k$MET_Data_Normal)),1]=2
                    BatchData=data.frame(ArrayName=colnames(ProcessedData$MET_Data_Normal),SampleName=colnames(ProcessedData$MET_Data_Normal),Batch=Batch)
                    ProcessedData$MET_Data_Normal=TCGA_BatchCorrection_MolecularData(ProcessedData$MET_Data_Normal,BatchData,5)
               }               
               
          } else if (ncol(ProcessedData27k$MET_Data_Cancer)>ncol(ProcessedData450k$MET_Data_Cancer)) { 
          	cat("\tNot enough 450k samples, only using the 27k (need min 50 samples).\n")
               Mode='27k'
               ProcessedData=ProcessedData27k          
          } else {
          	cat("\tNot enough 27k samples, only using the 450k (need min 50 samples).\n")
               Mode='450k'
               ProcessedData=ProcessedData450k     
          }
     } else if (CancerSite == "LAML") {
          cat("\tLAML is a special case, only using 450k data.\n")
          OverlapSamplesCancer=intersect(colnames(ProcessedData27k$MET_Data_Cancer),colnames(ProcessedData450k$MET_Data_Cancer))
          cat("\tOverlap length is:",length(OverlapSamplesCancer),".\n")
          Mode='450k'          
          ProcessedData=ProcessedData450k          
     } else if (length(ProcessedData27k)!=0) {
     	  cat("\tOnly 27k samples.\n")
          Mode='27k'
          ProcessedData=ProcessedData27k          
     } else {       
     	  cat("\tOnly 450k samples.\n")
          Mode='450k'
          ProcessedData=ProcessedData450k          
     }     
     return(ProcessedData=ProcessedData)
}



Preprocess_CancerSite_Methylation27k <- function(CancerSite,METdirectory) {
     
    # Settings
    data("BatchData")
    MinPerBatch=5
    MissingValueThreshold=0.2
    
    # Load data
    METfiles=dir(METdirectory)
    MatchedFilePosition=grep('methylation__humanmethylation27',METfiles)             
    Filename=paste0(METdirectory,METfiles[MatchedFilePosition])          
    MET_Data=TCGA_GENERIC_LoadIlluminaMethylationData(Filename)
    
    # Split up normal and cancer data
    Samplegroups=TCGA_GENERIC_GetSampleGroups(colnames(MET_Data))
    if (CancerSite=='LAML') {
          MET_Data_Cancer=MET_Data[,Samplegroups$PeripheralBloodCancer,drop=FALSE]
    } else {
     	MET_Data_Cancer=MET_Data[,Samplegroups$Primary]
    }
    MET_Data_Cancer=TCGA_GENERIC_CleanUpSampleNames(MET_Data_Cancer,15)
    if (CancerSite=='LAML') {
     	MET_Data_Normal=MET_Data[,Samplegroups$BloodNormal,drop=FALSE]
    } else {
     	MET_Data_Normal=MET_Data[,Samplegroups$SolidNormal,drop=FALSE]
    }
    if (length(MET_Data_Normal)>0) {
     	MET_Data_Normal=TCGA_GENERIC_CleanUpSampleNames(MET_Data_Normal,15)
    }
    cat("There are",length(colnames(MET_Data_Cancer)),"cancer samples and",length(colnames(MET_Data_Normal)),"normal samples in",CancerSite,"\n")
    
    # Missing value estimation
    cat("\tMissing value estimation for the cancer samples.\n")
    MET_Data_Cancer=TCGA_Process_EstimateMissingValues(MET_Data_Cancer,MissingValueThreshold)
    if (length(MET_Data_Normal)>0) {
         	cat("\tMissing value estimation for the normal samples.\n")
         	MET_Data_Normal=TCGA_Process_EstimateMissingValues(MET_Data_Normal,MissingValueThreshold)
    }
    
    # Batch correction for cancer and normal. 
    cat("\tBatch correction for the cancer samples.\n")
    BatchEffectCheck=TCGA_GENERIC_CheckBatchEffect(MET_Data_Cancer,BatchData)
    MET_Data_Cancer=TCGA_BatchCorrection_MolecularData(MET_Data_Cancer,BatchData,MinPerBatch)
    
    if (length(MET_Data_Normal)>0) {
         	cat("\tBatch correction for the normal samples.\n")
         	BatchEffectCheck=TCGA_GENERIC_CheckBatchEffect(MET_Data_Normal,BatchData)
         	MET_Data_Normal=TCGA_BatchCorrection_MolecularData(MET_Data_Normal,BatchData,MinPerBatch)
    } else {
    	     MET_Data_Normal=c()
    }
    
    # Reducing to 12 ids. 
    #MET_Data_Cancer=TCGA_GENERIC_CleanUpSampleNames(MET_Data_Cancer,12)
    #MET_Data_Normal=TCGA_GENERIC_CleanUpSampleNames(MET_Data_Normal,12)     
    
    MET_Data_Cancer[MET_Data_Cancer<0]=0
    MET_Data_Cancer[MET_Data_Cancer>1]=1
    if (length(MET_Data_Normal)>0) {
         	MET_Data_Normal[MET_Data_Normal<0]=0
         	MET_Data_Normal[MET_Data_Normal>1]=1
    }
    
    return(list(MET_Data_Cancer=MET_Data_Cancer,MET_Data_Normal=MET_Data_Normal))
}

Preprocess_CancerSite_Methylation450k_KB <- function(CancerSite,METdirectory) {
     
     #data("BatchData")
     MinPerBatch=5
     MissingValueThreshold=0.2
     
     METfiles=dir(METdirectory)
     MatchedFilePosition=grep('methylation__humanmethylation450',METfiles)             
     Filename=paste0(METdirectory,METfiles[MatchedFilePosition])     
     MET_Data=TCGA_GENERIC_LoadIlluminaMethylationData(Filename)
          
     # Split up normal and cancer data
     Samplegroups=TCGA_GENERIC_GetSampleGroups(colnames(MET_Data))
     if (CancerSite=='LAML') {
         MET_Data_Cancer=MET_Data[,Samplegroups$PeripheralBloodCancer]
     } else {
         MET_Data_Cancer=MET_Data[,Samplegroups$Primary]
     }
     MET_Data_Cancer=TCGA_GENERIC_CleanUpSampleNames(MET_Data_Cancer,15)
     if (CancerSite=='LAML') {
         MET_Data_Normal=MET_Data[,Samplegroups$BloodNormal]
     } else {
         MET_Data_Normal=MET_Data[,Samplegroups$SolidNormal]
     }     
     if (length(MET_Data_Normal)>0) {
          MET_Data_Normal=TCGA_GENERIC_CleanUpSampleNames(MET_Data_Normal,15)
     }
     cat("There are",length(colnames(MET_Data_Cancer)),"cancer samples and",length(colnames(MET_Data_Normal)),"normal samples in",CancerSite,"\n")
     
     # Clear space
     rm(MET_Data); gc()
     
     # Missing value estimation
     cat("\tMissing value estimation for the cancer samples.\n")
     MET_Data_Cancer=TCGA_Process_EstimateMissingValues(MET_Data_Cancer,MissingValueThreshold)     
     if (length(MET_Data_Normal)>0) {
          cat("\tMissing value estimation for the normal samples.\n")
          MET_Data_Normal=TCGA_Process_EstimateMissingValues(MET_Data_Normal,MissingValueThreshold)
     }     
     
     # Split up cancer data set
     MET_Data_Cancer1=MET_Data_Cancer[1:50000,]
     MET_Data_Cancer2=MET_Data_Cancer[50001:100000,]
     MET_Data_Cancer3=MET_Data_Cancer[100001:150000,]
     MET_Data_Cancer4=MET_Data_Cancer[150001:200000,]
     MET_Data_Cancer5=MET_Data_Cancer[200001:250000,]
     MET_Data_Cancer6=MET_Data_Cancer[250001:300000,]
     MET_Data_Cancer7=MET_Data_Cancer[300001:350000,]
     MET_Data_Cancer8=MET_Data_Cancer[350001:nrow(MET_Data_Cancer),]
     
     # clearing some memory
     rm(MET_Data_Cancer); gc()
     
     # Split up normal data set
     if (length(MET_Data_Normal)>0) {
          MET_Data_Normal1=MET_Data_Normal[1:100000,,drop=FALSE]
          MET_Data_Normal2=MET_Data_Normal[100001:200000,,drop=FALSE]
          MET_Data_Normal3=MET_Data_Normal[200001:300000,,drop=FALSE]
          MET_Data_Normal4=MET_Data_Normal[300001:nrow(MET_Data_Normal),,drop=FALSE]
          # clearing some memory
          rm(MET_Data_Normal); gc()
     } else {
          MET_Data_Normal1=c()
     }

     cat("\tBatch correction for the cancer samples.\n")
     MET_Data_Cancer1=TCGA_BatchCorrection_MolecularData(MET_Data_Cancer1,BatchData,MinPerBatch)
     MET_Data_Cancer2=TCGA_BatchCorrection_MolecularData(MET_Data_Cancer2,BatchData,MinPerBatch)
     MET_Data_Cancer3=TCGA_BatchCorrection_MolecularData(MET_Data_Cancer3,BatchData,MinPerBatch)
     MET_Data_Cancer4=TCGA_BatchCorrection_MolecularData(MET_Data_Cancer4,BatchData,MinPerBatch)
     MET_Data_Cancer5=TCGA_BatchCorrection_MolecularData(MET_Data_Cancer5,BatchData,MinPerBatch)
     MET_Data_Cancer6=TCGA_BatchCorrection_MolecularData(MET_Data_Cancer6,BatchData,MinPerBatch)
     MET_Data_Cancer7=TCGA_BatchCorrection_MolecularData(MET_Data_Cancer7,BatchData,MinPerBatch)
     MET_Data_Cancer8=TCGA_BatchCorrection_MolecularData(MET_Data_Cancer8,BatchData,MinPerBatch)

    if (length(MET_Data_Normal1)>0) {
         cat("\tBatch correction for the normal samples.\n")
          MET_Data_Normal1=TCGA_BatchCorrection_MolecularData(MET_Data_Normal1,BatchData,2)
          MET_Data_Normal2=TCGA_BatchCorrection_MolecularData(MET_Data_Normal2,BatchData,2)
          MET_Data_Normal3=TCGA_BatchCorrection_MolecularData(MET_Data_Normal3,BatchData,2)
          MET_Data_Normal4=TCGA_BatchCorrection_MolecularData(MET_Data_Normal4,BatchData,2)
     }
        
     # Combine batch corrected data
     # It's possible that samples get deleted due to missings in one part and not another.
     OverlapSamples=Reduce(intersect,list(colnames(MET_Data_Cancer1),colnames(MET_Data_Cancer2),colnames(MET_Data_Cancer3),colnames(MET_Data_Cancer4)
                              ,colnames(MET_Data_Cancer5),colnames(MET_Data_Cancer6),colnames(MET_Data_Cancer7),colnames(MET_Data_Cancer8)))
     MET_Data_Cancer=rbind(MET_Data_Cancer1[,OverlapSamples],MET_Data_Cancer2[,OverlapSamples],MET_Data_Cancer3[,OverlapSamples]
                          ,MET_Data_Cancer4[,OverlapSamples],MET_Data_Cancer5[,OverlapSamples],MET_Data_Cancer6[,OverlapSamples]
                          ,MET_Data_Cancer7[,OverlapSamples],MET_Data_Cancer8[,OverlapSamples])
    
    rm(MET_Data_Cancer1);rm(MET_Data_Cancer2);rm(MET_Data_Cancer3);rm(MET_Data_Cancer4);rm(MET_Data_Cancer5);rm(MET_Data_Cancer6);rm(MET_Data_Cancer7);rm(MET_Data_Cancer8);gc()
    
     if (length(MET_Data_Normal1)>0) {
          OverlapSamples=Reduce(intersect,list(colnames(MET_Data_Normal1),colnames(MET_Data_Normal2),colnames(MET_Data_Normal3),colnames(MET_Data_Normal4)))
          MET_Data_Normal=rbind(MET_Data_Normal1[,OverlapSamples,drop=FALSE],MET_Data_Normal2[,OverlapSamples,drop=FALSE],
                                MET_Data_Normal3[,OverlapSamples,drop=FALSE],MET_Data_Normal4[,OverlapSamples,drop=FALSE])
     } else {
          MET_Data_Normal=c()
     }

    rm(MET_Data_Normal1);rm(MET_Data_Normal2); rm(MET_Data_Normal3);rm(MET_Data_Normal4);gc()
    
     # Set values <0 to 0 and >1 to 1, because of batch correction
     MET_Data_Cancer[MET_Data_Cancer<0]=0
     MET_Data_Cancer[MET_Data_Cancer>1]=1
     if (length(MET_Data_Normal)>0) {
          MET_Data_Normal[MET_Data_Normal<0]=0
          MET_Data_Normal[MET_Data_Normal>1]=1
     }
     
     return(list(MET_Data_Cancer=MET_Data_Cancer,MET_Data_Normal=MET_Data_Normal))
}

Preprocess_CancerSite_Methylation450k_BRCA <- function(CancerSite,METdirectory) {
    
    cat("\tStarting Preprocess of BRCA")
    
    data("BatchData")
    MinPerBatch=5
    MissingValueThreshold=0.2
    
    METfiles=dir(METdirectory)
    MatchedFilePosition=grep('methylation__humanmethylation450',METfiles)             
    Filename=paste0(METdirectory,METfiles[MatchedFilePosition])     
    MET_Data=TCGA_GENERIC_LoadIlluminaMethylationData(Filename)
    
    # Split up normal and cancer data
    Samplegroups=TCGA_GENERIC_GetSampleGroups(colnames(MET_Data))
    MET_Data_Cancer=MET_Data[,Samplegroups$Primary]
    MET_Data_Cancer=TCGA_GENERIC_CleanUpSampleNames(MET_Data_Cancer,15)
    MET_Data_Normal=MET_Data[,Samplegroups$SolidNormal]   
    if (length(MET_Data_Normal)>0) {
        MET_Data_Normal=TCGA_GENERIC_CleanUpSampleNames(MET_Data_Normal,15)
    }
    cat("There are",length(colnames(MET_Data_Cancer)),"cancer samples and",length(colnames(MET_Data_Normal)),"normal samples in",CancerSite,"\n")
    
    # Clear space
    rm(MET_Data); gc()
    
    # Missing value estimation
    cat("\tMissing value estimation for the cancer samples.\n")
    MET_Data_Cancer=TCGA_Process_EstimateMissingValues(MET_Data_Cancer,MissingValueThreshold)     
    if (length(MET_Data_Normal)>0) {
        cat("\tMissing value estimation for the normal samples.\n")
        MET_Data_Normal=TCGA_Process_EstimateMissingValues(MET_Data_Normal,MissingValueThreshold)
    }     
    
    
    # Split up cancer data set in parts of nprobes number of probes
    nprobes = 10000
    nrows = nrow(MET_Data_Cancer)
    nparts_Cancer = ceiling(nrows/nprobes)
    cat("Splitting cancer data in", nparts_Cancer, "parts of", nprobes, "\n")
    for (i in 1:nparts_Cancer) {
        first = nprobes*(i-1)+1
        last = ifelse(i==nparts_Cancer, nrows, nprobes*i)
        name = paste0("MET_Data_Cancer", i)
        assign(name, MET_Data_Cancer[first:last, , drop = F])
        # Save and delete this piece
        save(list = name, file = paste0("/home/marcosp/temp_brca_pieces/", name, ".RData"))
        rm(list = name); gc()
    }
    cat("Splitting cancer data finished.\n")
    
    # clearing some memory
    MET_Data_Cancer=c()
    
    # Split up normal data set in parts of nprobes number of probes
    if (length(MET_Data_Normal)>0) {
        has.normal = TRUE
        nprobes = 10000
        nrows = nrow(MET_Data_Normal)
        nparts_Normal = ceiling(nrows/nprobes)
        cat("Splitting normal data in", nparts_Normal, "parts of", nprobes, "\n")
        for (i in 1:nparts_Normal) {
            first = nprobes*(i-1)+1
            last = ifelse(i==nparts_Normal, nrows, nprobes*i)
            name = paste0("MET_Data_Normal", i)
            assign(name, MET_Data_Normal[first:last, , drop = F])
            # Save and delete this piece
            save(list = name, file = paste0("/home/marcosp/temp_brca_pieces/", name, ".RData"))
            rm(list = name); gc()
        }
        # clearing some memory
        MET_Data_Normal=c()
        cat("Splitting normal data finished.\n")
    } else {
        cat("No normal data to split.\n")
        has.normal = FALSE
    }
    
    cat("\tBatch correction for the cancer samples.\n")
    # list for keeping sample names
    list_tmp_cancer = list()
    for (i in 1:nparts_Cancer) {
        name = paste0("MET_Data_Cancer", i)
        cat("Doing Batch Correction", i, "/", nparts_Cancer, "\n")
        # Load piece of data
        load(paste0("/home/marcosp/temp_brca_pieces/", name, ".RData"))
        # Batch correction
        assign(name, TCGA_BatchCorrection_MolecularData(get(name),BatchData,MinPerBatch))
        # Save samples names
        list_tmp_cancer[[i]] = colnames(get(name))
        # Save and delete
        save(list = name, file = paste0("/home/marcosp/temp_brca_pieces/", name, ".RData"))
        rm(list = name); gc()
    }
    
    if (has.normal) {
        # list for keeping sample names
        list_tmp_normal = list()
        cat("\tBatch correction for the normal samples.\n")
        for (i in 1:nparts_Normal) {
            name = paste0("MET_Data_Normal", i)
            cat("Doing Batch Correction", i, "/", nparts_Normal, "\n")
            # Load piece of data
            load(paste0("/home/marcosp/temp_brca_pieces/", name, ".RData"))
            # Batch correction
            assign(name, TCGA_BatchCorrection_MolecularData(get(name),BatchData,2))
            # Save samples names
            list_tmp_normal[[i]] = colnames(get(name))
            # Save and delete
            save(list = name, file = paste0("/home/marcosp/temp_brca_pieces/", name, ".RData"))
            rm(list = name); gc()
        } 
    }
    
    # Combine batch corrected data - Cancer
    # It's possible that samples get deleted due to missings in one part and not another.
    cat("Combining batch corrected data - Cancer\n")
    OverlapSamples = Reduce(intersect, list_tmp_cancer)
    rm(list_tmp_cancer); gc()
    load(paste0("/home/marcosp/temp_brca_pieces/MET_Data_Cancer1.RData"))
    MET_Data_Cancer = MET_Data_Cancer1[, OverlapSamples]
    rm(MET_Data_Cancer1); gc()
    for (i in 2:nparts_Cancer) {
        name = paste0("MET_Data_Cancer", i)
        cat("Loading part", i, "/", nparts_Cancer, "\n")
        load(paste0("/home/marcosp/temp_brca_pieces/", name, ".RData"))
        MET_Data_Cancer = rbind(MET_Data_Cancer, get(name)[,OverlapSamples])
        rm(list = name); gc()
    }
    
    # Combine batch corrected data - Normal
    if (has.normal) {
        cat("Combining batch corrected data - Normal\n")
        OverlapSamples=Reduce(intersect,list_tmp_normal)
        rm(list_tmp_normal); gc()
        load(paste0("/home/marcosp/temp_brca_pieces/MET_Data_Normal1.RData"))
        MET_Data_Normal = MET_Data_Normal1[, OverlapSamples, drop = FALSE]
        rm(MET_Data_Normal1); gc()
        for (i in 2:nparts_Normal) {
            name = paste0("MET_Data_Normal", i)
            cat("Loading part", i, "/", nparts_Normal, "\n")
            load(paste0("/home/marcosp/temp_brca_pieces/", name, ".RData"))
            MET_Data_Normal = rbind(MET_Data_Normal, get(name)[, OverlapSamples, drop = FALSE])
            rm(list = name); gc()
        }
    } else {
        cat("No normal data to combine\n")
        MET_Data_Normal=c()
    }
    
    # Set values <0 to 0 and >1 to 1, because of batch correction
    cat("Setting values to be in range.\n")
    MET_Data_Cancer[MET_Data_Cancer<0]=0
    MET_Data_Cancer[MET_Data_Cancer>1]=1
    if (length(MET_Data_Normal)>0) {
        MET_Data_Normal[MET_Data_Normal<0]=0
        MET_Data_Normal[MET_Data_Normal>1]=1
    }
    cat("All finished.\n")
    return(list(MET_Data_Cancer=MET_Data_Cancer,MET_Data_Normal=MET_Data_Normal))
    
}


TCGA_GENERIC_LoadIlluminaMethylationData <- function(Filename) {
     
     # read in an illumina methylation file with the following format: 
     # header row with sample labels
     # 2nd header row with 4 columns per sample: beta-value, geneSymbol, chromosome and GenomicCoordinate
     # The first column has the probe names. 
     MET_Data<-fread(Filename)
     MET_Data=as.matrix(MET_Data)
     Probes=MET_Data[,1]
     rownames(MET_Data)=Probes
     MET_Data=MET_Data[,-1]
     MET_Data=MET_Data[-1,]
     MET_Data=MET_Data[,seq(1,ncol(MET_Data),4)]     
     class(MET_Data)='numeric'
     
     return(MET_Data)
}


















