pv.peaks2DataType = function(peaks,datatype=DBA_DATA_DEFAULT) {
   
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
                          bPCA=F,bLog=T,minval,maxval) {
                          	
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
   pv$class = cbind(pv$class[,con$group1],pv$class[,con$group2])
   
   if(bPCA)  {
   	
      #pv$pc = princomp(domap,cor=bPCAcor)

      if(attributes[1] == PV_GROUP) {
         pv$class[PV_ID,] = c(rep(con$name1,sum(con$group1)),rep(con$name2,sum(con$group2)))
      }
   } 
   
   return(pv)     
}


pv.get_reads = function(pv,peaksets,bSubControl=T){
   reads = NULL
   for(peakset in peaksets) {
      reads = cbind(reads,pv$peaks[[peakset]][,PV_CHIP_READS])
      if(bSubControl) {
         reads[,ncol(reads)] = reads[,ncol(reads)] - pv$peaks[[peakset]][,PV_CONTROL_READS]
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

PV_CHIP_RPKM      = 5
PV_CHIP_READS     = 6
PV_CONTROL_RPKM   = 7
PV_CONTROL_READS  = 8


pv.get_scores = function(pv,peaksets){
   control = rep(0,nrow(pv$peaks[[peaksets[1]]]))
   scores = NULL
   for(peakset in peaksets) {
      control = control + pv$peaks[[peakset]][,PV_CONTROL_RPKM]
      scores = cbind(scores,pv$peaks[[peakset]][,PV_CHIP_RPKM])
   }
   control = control/length(peaksets)
   for(i in 1:ncol(scores)) {
      scores[,i] = log2(scores[,i] / control)	
   }
   return(scores)
}



