pv.peaks2DataType = function(peaks,datatype=DBA_DATA_DEFAULT) {

   res = peaks
   
   if(is.null(peaks)) {
      return(NULL)
   }   
   
   if(datatype==DBA_DATA_FRAME) {
      return(peaks)	
   }
   
   if(datatype==DBA_DATA_RANGEDDATA) {
      if(class(peaks)!="RangedData") {
         #require(IRanges)
         cnames = colnames(peaks)
         cnames[1:3] = c("space","start","end")
         colnames(peaks) = cnames
      }
      res = as(peaks,"RangedData")
   }

   if(datatype==DBA_DATA_GRANGES) {
      #require(GenomicRanges)
      if(class(peaks)=="RangedData") {
         res = suppressWarnings(as(peaks,"GRanges"))
      } else if (class(peaks) != "GRanges") {
         res = GRanges(Rle(peaks[,1]),IRanges(peaks[,2],width=peaks[,3]-peaks[,2]+1,names=rownames(peaks)),
                       strand = Rle("*", length(seqnames)))
         if(ncol(peaks)>3) {
         	mdata = data.frame(peaks[,4:ncol(peaks)])
         	colnames(mdata)  = colnames(peaks)[4:ncol(peaks)]
         	elementMetadata(res) = mdata
         }   	
      }
   }   
   
   return(res)
   
}

pv.DataType2Peaks = function(RDpeaks){
   if(class(RDpeaks)=='logical') {
     return(RDpeaks)
   }
   if(class(RDpeaks)=='integer') {
     return(RDpeaks)
   }
   if(class(RDpeaks)=='numeric') {
     return(RDpeaks)
   }
   if(class(RDpeaks)=='character') {
     return(RDpeaks)
   }
   if(class(RDpeaks)!="data.frame") {
      #require(IRanges)
      res = as.data.frame(RDpeaks)[,-4]
      cnames = colnames(res)
      cnames[1:3] = c("CHR","START","END")
      colnames(res) = cnames
      if(class(RDpeaks)=="GRanges") {
         res = res[,-4]	
      }
   } else {
      res = RDpeaks
   }
   return(res)
}

pv.getPlotData = function(pv,attributes=PV_GROUP,contrast=1,method=DBA_EDGER,th=.1,bUsePval=FALSE,bNormalized=T,report,
                          bPCA=F,bLog=T,minval,maxval,mask) {
                          	
   if(contrast > length(pv$contrasts)) {
      stop('Specified contrast number is greater than number of contrasts')
      return(NULL)
   }
  
   con = pv$contrasts[[contrast]]
  
   if(missing(report)) {
     report = pv.DBAreport(pv,contrast=contrast,method=method,th=th,bUsePval=bUsePval,
                           bNormalized=bNormalized,bCounts=T,bSupressWarning=T)
     if(is.null(report)) {
        stop('Unable to plot -- no sites within threshold')	
     }
   }
   
   if(!missing(mask)){
      if (!is.logical(mask)) {
         if (max(mask) > length(pv$peaks)) {
            stop("Invalid sample number in mask.",call.=F)
         }
         temp = rep(F, length(pv$peaks))
         temp[mask] = T
         mask = temp
      }
      sites = as.numeric(rownames(report))
      group1 = con$group1 & mask
      group2 = con$group2 & mask
      extra = mask & !(group1 | group2)
      allsamps = c(which(group1), which(group2), which(extra))
      numsamps = length(allsamps)
      domap = matrix(0,length(sites),0)
      if(sum(group1)) {     
         domap = cbind(domap,pv$vectors[sites,3+which(group1)])
      }
      if(sum(group2)) {
         domap = cbind(domap,pv$vectors[sites,3+which(group2)])
      }
      if(sum(extra)) {
         domap = cbind(domap,pv$vectors[sites,3+which(extra)])
      }
      rownames(domap) = rownames(report)
      colnames(domap) = pv$class[PV_ID,allsamps] 
      con$group1 = group1
      con$group2 = group2
      peaks = pv$peaks[allsamps]
      for(i in 1:length(peaks)){
         peaks[[i]] = peaks[[i]][sites,]
      }  
   } else {   
      allsamps = c(which(con$group1),which(con$group2))
      extra = rep(F,ncol(pv$class)) 
      repcols = colnames(report)
      numsamps = sum(con$group1)+sum(con$group2)
      if(length(repcols) < (numsamps+9)) {
         stop('Report does not have count data, re-run dba.report with bCounts=T')
      }
      first = 10
      if(repcols[10]=="Called1") {
         if(length(repcols) < (numsamps+11)) {
            stop('Report does not have count data, re-run dba.report with bCounts=T')
         }
         first = 12	 
      }
      domap = report[,first:(first+numsamps-1)]
      group1 = rep(F,numsamps)
      group2 = rep(T,numsamps)
      group1[1:sum(con$group1)] = T
      group2[1:sum(con$group1)] = F
      con$group1 = group1
      con$group2 = group2

      sites = as.numeric(rownames(report))
      peaks = pv$peaks[allsamps]
      for(i in 1:length(peaks)){
         peaks[[i]] = peaks[[i]][sites,]
      }  
   }
   
   if(bLog) {
      domap[domap<=0]=1
      domap = log2(domap)
      if(missing(minval)) {
         minval = 0
      } else {
      	minval = max(0,minval)
      }
   }
   
   if(!missing(minval)) {
      domap[domap< minval]= minval
   }
   if(!missing(maxval)) {
      domap[domap>maxval] = maxval
   }
     
   pv$vectors = cbind(report[,1:3],domap)
   pv$allvectors = pv$vectors
   pv$class = pv$class[,allsamps]
   pv$peaks = peaks
   pv$sites = sites
   pv$contrasts = list(pv$contrasts[[contrast]])
   pv$contrasts[[1]]$group1 = rep(F,ncol(pv$class))
   pv$contrasts[[1]]$group2 = rep(F,ncol(pv$class))
   if(sum(con$group1)) {
      pv$contrasts[[1]]$group1[1:sum(con$group1)]=T
   }
   if(sum(con$group2)) {
      pv$contrasts[[1]]$group2[(sum(con$group1)+1):(sum(con$group1)+sum(con$group2))]=T   
   }
      
   if(bPCA)  {
   	
      #pv$pc = princomp(domap,cor=bPCAcor)

      if(attributes[1] == PV_GROUP) {
         pv$class[PV_ID,] = c( rep(con$name1,sum(con$group1)), rep(con$name2,sum(con$group2)), rep("none",sum(extra)) )
      }
   } 
   
   return(pv)     
}


pv.get_reads = function(pv,peaksets,bSubControl=T){
   if(is.null(bSubControl)) {
      bSubControl = T
   }
   reads = NULL
   for(peakset in peaksets) {
      reads = cbind(reads,pv$peaks[[peakset]]$Reads)
      if(bSubControl) {
         reads[,ncol(reads)] = reads[,ncol(reads)] - pv$peaks[[peakset]]$cReads
      }
   }

   reads[reads<1]=1

   return(reads)
}

pv.getoneorNA = function(vec) {
 res = unique(vec)
 if (length(res) == 1) {
    return(res)
 } else {
   return(NA)
 }
}

pv.add_consensus = function(model,psets){
   
   if(missing(psets)) {
      toadd = unique(model$class[PV_ID,duplicated(model$class[PV_ID,])])
      for(pset in toadd) {
         tomerge = which(model$class[PV_ID,] %in% pset)
         model = pv.consensus(model,tomerge,minOverlap=2)
      }   
      return(model) 
   }
   
   if(is.null(psets)) {
      return(model)
   }
   
   model$psets = psets
   for(pset in psets){
      if(length(pset)>1){
         model = pv.consensus(model,pset,minOverlap=2)
      }   
   }
   return(model)	
}

pv.removeComp = function (data,numRemove=0) {

   if(numRemove == 0) {
      return(data)
   }
   
   res = svd(data)
   U   = res$u
   d   = res$d
   Vt  = t(res$v)
   
   tot = sum(d)
   run = 0
   for(i in 1:numRemove) {
   	  run = run + d[i]
      #cat(i,d[i],d[i]/tot,run/tot,'\n')
      d[i] = 0
   }	
   M =  U %*% diag(d) %*% Vt
   return(M)
}

pv.occ2matrix = function(occ,col=PV_COR,div){
   #occ[,1] = as.character(occ[,1])
   #occ[,2] = as.character(occ[,2])
   els = unique(c(occ[,1],occ[,2]))
   
   nels = length(els)
   els  = els[order(els,decreasing=F)]
   res = matrix(NA,nels,nels)

   for(i in 1:nels) {
      for(j in 1:nels) {
         if(i==j) {
            res[i,j] = 1
         } else {
           entry = (occ[,1] == els[i]) & (occ[,2] == els[j])
           entry = entry | (occ[,1] == els[j]) & (occ[,2] == els[i])
           if(sum(entry)>0) {
              res[i,j] = as.numeric(occ[entry,col])
              if(!missing(div)) {
              	 if(div != PV_TOTAL) {
              	    todiv = as.numeric(occ[entry,div])
              	 } else {
              	    todiv = as.numeric(occ[entry,PV_ONLYA]) +
              	            as.numeric(occ[entry,PV_ONLYB]) +
              	            as.numeric(occ[entry,PV_INALL])
              	 }
                 res[i,j] = res[i,j] / todiv
              }
           }
         }
      }
   } 
  
  colnames(res) = els
  rownames(res) = els
  return(res)
}

pv.nums2labels = function(pv,els,atts) {
   if(is.null(atts)) {
      return(els)
   }
   for(i in 1:length(els)) {
       els[i] = pv.namestrings(pv$class[atts,as.numeric(els[i])])$tstring
   }
   return(els)
}

pv.getOverlapData = function(pv,contrast,report) {
   con = pv$contrasts[[contrast]]                       	
   repcols = colnames(report)
   numsamps = sum(con$group1)+sum(con$group2)
   if(length(repcols) < (numsamps+9)) {
      return(pv)
   }
   first = 10
   if(repcols[10]=="Called1") {
      if(length(repcols) < (numsamps+11)) {
         return(pv)
      }
      first = 12	 
   }
   
   domap = report[,first:(first+numsamps-1)]
   
   pv$vectors    = cbind(report[,1:3],domap)
   pv$allvectors = pv$vectors
   pv$class      = cbind(pv$class[,con$group1],pv$class[,con$group2])
   
   return(pv)     
}

#PV_CHIP_RPKM      = 6
#PV_CHIP_READS     = 7
#PV_CONTROL_RPKM   = 8
#PV_CONTROL_READS  = 9


pv.get_scores = function(pv,peaksets){
   control = rep(0,nrow(pv$peaks[[peaksets[1]]]))
   scores = NULL
   for(peakset in peaksets) {
      control = control + pv$peaks[[peakset]]$cRPKM
      scores = cbind(scores,pv$peaks[[peakset]]$RPKM)
   }
   control = control/length(peaksets)
   for(i in 1:ncol(scores)) {
      scores[,i] = log2(scores[,i] / control)	
   }
   return(scores)
}

### Thanks to Obi Griffith 
### Obtained at http://www.biostars.org/post/show/18211/how-do-i-draw-a-heatmap-in-r-with-both-a-color-key-and-multiple-color-side-bars/
### Downloaded from http://dl.dropbox.com/u/16769159/Biostar/heatmap.3.R
checkinvalid = function (x) 
{
    if (missing(x) || is.null(x) || length(x) == 0) 
        return(TRUE)
    if (is.list(x)) 
        return(all(sapply(x, checkinvalid)))
    else if (is.vector(x)) 
        return(all(is.na(x)))
    else return(FALSE)
}


heatmap.3=function (x, Rowv = TRUE, Colv = if (symm) "Rowv" else TRUE, 
    distfun = dist, hclustfun = hclust, dendrogram = c("both", 
        "row", "column", "none"), symm = FALSE, scale = c("none", 
        "row", "column"), na.rm = TRUE, revC = identical(Colv, 
        "Rowv"), add.expr, breaks, symbreaks = min(x < 0, na.rm = TRUE) || 
        scale != "none", col = "heat.colors", colsep, rowsep, 
    sepcolor = "white", sepwidth = c(0.05, 0.05), cellnote, notecex = 1, 
    notecol = "cyan", na.color = par("bg"), trace = c("column", 
        "row", "both", "none"), tracecol = "cyan", hline = median(breaks), 
    vline = median(breaks), linecol = tracecol, margins = c(5, 
        5), ColSideColors, RowSideColors, cexRow = 0.2 + 1/log10(nr), 
    cexCol = 0.2 + 1/log10(nc), labRow = NULL, labCol = NULL, 
    key = TRUE, keysize = 1.5, density.info = c("histogram", 
        "density", "none"), denscol = tracecol, symkey = min(x < 
        0, na.rm = TRUE) || symbreaks, densadj = 0.25, main = NULL, 
    xlab = NULL, ylab = NULL, lmat = NULL, lhei = NULL, lwid = NULL,
    NumColSideColors = 1, NumRowSideColors = 1, KeyValueName="Value",
    ...) 
{
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
    if (missing(breaks) || is.null(breaks) || length(breaks) < 
        1) {
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
      if(!is.null(ColSideColors)) {	
       #if (!is.matrix(ColSideColors)) 
       #stop("'ColSideColors' must be a matrix")
        if (!is.character(ColSideColors) || dim(ColSideColors)[1] != 
            nc) 
            stop("'ColSideColors' dim()[2] must be of length ncol(x)")
        lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
        #lhei <- c(lhei[1], 0.2, lhei[2])
         lhei=c(lhei[1], 0.1*NumColSideColors, lhei[2]) 
      }
    }
    if (!missing(RowSideColors)) {
      if(!is.null(RowSideColors)) { 
       #if (!is.matrix(RowSideColors)) 
       #stop("'RowSideColors' must be a matrix")
        if (!is.character(RowSideColors) || dim(RowSideColors)[1] != 
            nr) 
            stop("'RowSideColors' must be a character vector of length nrow(x)")
        lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1), 1), lmat[,2] + 1)
        #lwid <- c(lwid[1], 0.2, lwid[2])
        lwid <- c(lwid[1], 0.1*NumRowSideColors, lwid[2])
      }
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
     if (!is.null(RowSideColors)) {
	    #if (!is.matrix(RowSideColors)){
	    if(revC) {
	      xrowInd = rev(rowInd)	
	    } else {
	      xrowInd = rowInd	
	    }	
	    if(ncol(RowSideColors)==1) {	
       		par(mar = c(margins[1], 0, 0, 0.5))
        	image(rbind(1:nr), col = RowSideColors[,1][xrowInd], axes = FALSE)
        	axis(1,0,colnames(RowSideColors)[1],las=2,tick=F)
        } 
        else{
        par(mar = c(margins[1], 0, 0, 0.5))
        rsc = RowSideColors[xrowInd, ]
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
        if (length(colnames(RowSideColors)) > 0) {
            axis(1, 0:(dim(rsc)[2] - 1)/(dim(rsc)[2] - 1),
                colnames(RowSideColors), 
                las = 2, tick = FALSE)
        }
     }
    }
 }

    if (!missing(ColSideColors)) {
     if (!is.null(ColSideColors)) {
        #if (!is.matrix(ColSideColors)){
         if(ncol(ColSideColors)==1) {
        	par(mar = c(0.5, 0, 0, margins[2]))
        	image(cbind(1:nc), col = ColSideColors[,1][colInd], axes = FALSE)
        	axis(2,0,colnames(ColSideColors)[1],las=2,tick=F)
        }
        else {    
    	 par(mar = c(0.5, 0, 0, margins[2]))
         csc = ColSideColors[colInd, ]
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
            axis(2, 0:(dim(csc)[2] - 1)/(dim(csc)[2] - 1),
                 colnames(ColSideColors), 
                 las = 2, tick = FALSE)
         }
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
    image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + 
        c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col, 
        breaks = breaks, ...)
    retval$carpet <- x
    if (exists("ddr")) 
        retval$rowDendrogram <- ddr
    if (exists("ddc")) 
        retval$colDendrogram <- ddc
    retval$breaks <- breaks
    retval$col <- col
    if (!checkinvalid(na.color) & any(is.na(x))) {
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
        for (csep in colsep) rect(xleft = csep + 0.5, ybottom = rep(0, 
            length(csep)), xright = csep + 0.5 + sepwidth[1], 
            ytop = rep(ncol(x) + 1, csep), lty = 1, lwd = 1, 
            col = sepcolor, border = sepcolor)
    if (!missing(rowsep)) 
        for (rsep in rowsep) rect(xleft = 0, ybottom = (ncol(x) + 
            1 - rsep) - 0.5, xright = nrow(x) + 1, ytop = (ncol(x) + 
            1 - rsep) - 0.5 - sepwidth[2], lty = 1, lwd = 1, 
            col = sepcolor, border = sepcolor)
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
            tmpbreaks[1] <- -max(abs(x))
            tmpbreaks[length(tmpbreaks)] <- max(abs(x))
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


pv.attname = function(attribute) {
   if(attribute == PV_GROUP)     return("Group")
   if(attribute == PV_TISSUE)    return("Tissue")
   if(attribute == PV_FACTOR)    return("Factor")
   if(attribute == PV_CONDITION) return("Condition")
   if(attribute == PV_TREATMENT) return("Treatment")
   if(attribute == PV_REPLICATE) return("Replicate")
   if(attribute == PV_CALLER)    return("Caller")
   
   return("Group")
      	
}
pv.attributematrix = function(pv,mask,contrast,attributes,cols,bReverse=F,bAddGroup=F) {

   if(is.null(attributes)){
      attributes=c(PV_TISSUE,PV_FACTOR,PV_CONDITION,PV_TREATMENT,PV_REPLICATE,PV_CALLER)
      if(bAddGroup) {
         attributes=c(PV_GROUP,attributes)	
      }	
   }
   
   classdb = pv$class[,mask]
   numsamps = ncol(classdb)
   atts = NULL
   for(num in length(attributes):1) {
      attribute = attributes[num]
      if(attribute==PV_GROUP) {
         if(!missing(contrast)) {
            gps = rep(3,length(pv$contrasts[[contrast]]$group1))
            gps[pv$contrasts[[contrast]]$group1]=1
            gps[pv$contrasts[[contrast]]$group2]=2
            classdb = rbind(classdb,gps)
            attribute = nrow(classdb)
            vals = 1:max(gps)
         } else {
            vals = NA	
         }
      } else {
         vals = unique(classdb[attribute,])
      }
      if ( (sum(!is.na(vals))>1) && (sum(!is.na(vals))<numsamps) ) {
         addcol = matrix(classdb[attribute,],numsamps,1)
         for(i in 1:length(vals)) {
            addcol[addcol[,1]==vals[i],1]=cols[i]
         }
         colnames(addcol) = pv.attname(attribute)
         if(bReverse) {
            atts = cbind(addcol,atts)  
         } else { 	
            atts = cbind(atts,addcol)  	
         }
      }         	
   }
   return(atts)	
}

pv.attributePCA = function(DBA) {

   if(pv.morethanone(DBA,DBA_TISSUE))    return(DBA_TISSUE)	
   if(pv.morethanone(DBA,DBA_FACTOR))    return(DBA_FACTOR)	
   if(pv.morethanone(DBA,DBA_CONDITION)) return(DBA_CONDITION)	
   if(pv.morethanone(DBA,DBA_TREATMENT)) return(DBA_TREATMENT)	
   if(pv.morethanone(DBA,DBA_REPLICATE)) return(DBA_REPLICATE)	
   if(pv.morethanone(DBA,DBA_CALLER))    return(DBA_CALLER)

   return(DBA_ID)
}

pv.morethanone = function(DBA,att){

   vals = unique(DBA$class[att,])
   
   if(sum(!is.na(vals))>=2) {
      return(TRUE)
   } else {
      return(FALSE)	
   }	
}

pv.checkValue = function(val,check) {
   if(is.null(val)) {
      return(FALSE)
   } 
   if (val != check) {
      return(FALSE)
   }
   return(TRUE)
}



