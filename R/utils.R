
pv.change_ext = function(fname,ext){
  split = strsplit(as.character(fname),".bed")
  fn    = split[[1]][1]
  res = sprintf("%s.%s",fn,ext)
  return(res)
}

pv.find_sample = function(samples,sname) {

   if(is.na(sname)) {
      return(NULL)
   }
   
   matching = which(samples$GenomicsID %in% sname)
   return(matching[1])
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


pv.plot_contrasts = function(model,maxPlot=12){

   warning("pv.plot_contrasts unsupported")
   return(NULL)
   
   pnums = length(model$psets)
   if(pnums > maxPlot) {
      pnums = maxPlot
   }
   for(pnum in 1:pnums){
   	  pset = model$psets[[pnum]]
      if(length(pset)==2) {
         pv.contrast(model,pset[1],pset[2])
      } else if(length(pset)==3) {
         pv.contrast(model,pset[1],pset[2],pset[3])
      }
   }
   
   if(pnums < maxPlot){
      cons = which(pv.mask(model,PV_CONSENSUS,T))
      if(length(cons) == 2) {
         pv.contrast(model,cons[1],cons[2])
      }	else  if(length(cons) == 3) {
         pv.contrast(model,cons[1],cons[2],cons[3])
      }
   }
   return(NULL)
}

rmExt <- function(fileName,extension=c(".bed",".bam"))  {
  for (x in extension){
    fileName <- sub(x,"",fileName)
  }
  fileName
}



pv.bamReads = function(bamFile){
  res = pv.BAMstats(bamFile)
  return(res[1])	
}

pv.BAMstats <- function(bamFile){
   stats <-  system(paste("/usr/local/bin/samtools flagstat",bamFile),intern=T)
   stats <- as.integer(sapply(strsplit(stats[c(1,4,3)]," "),"[",1))
   stats[4] <- stats[1]-stats[3]
   stats[5] <- round((stats[3]/stats[1])*100,2)
   names(stats) <- c("Total","Mapped","Duplicates","Unique","DuplicationRate")
   return(stats)
}



pv.getCounts = function(bamfile,intervals,insertLength=0) {

   fdebug(sprintf('pv.getCounts: ENTER %s',bamfile))
   
   ## write peaks to a file 
#   bedf = tempfile(as.character(Sys.getpid()))#tmpdir='.')
#   pv.do_peaks2bed(intervals[,1:3],fn=bedf)

#   tmpFile = tempfile(as.character(Sys.getpid()))
#   fdebug(sprintf('pv.getCounts: call countReadsOnIntervals %s',bamfile))
#   system(paste("/home/brown22/bin/countReadsOnIntervals",bamfile,bedf,tmpFile),intern=T)
#   fdebug(sprintf('pv.getCounts: return %s',tmpFile))
   
#   counts = read.table(tmpFile)[,1]
#   unlink(tmpFile)
#   unlink(bedf)
   
#   counts[counts==0]=1
#   fdebug(sprintf('pv.getCounts: number of counts = %d',length(counts)))
#   libsize = pv.bamReads(bamfile)
#   fdebug(sprintf('pv.getCounts: libsize = %d',libsize))

   fdebug("Starting croi_load_reads...")
   bamtree <- .Call("croi_load_reads",as.character(bamfile),as.integer(insertLength))
   fdebug("Loaded...")
   libsize.croi <- .Call("croi_tree_size",bamtree)
   fdebug("Starting croi_count_reads...")
   counts.croi <- .Call("croi_count_reads",bamtree,
                                           as.character(intervals[[1]]),
                                           as.integer(intervals[[2]]),
                                           as.integer(intervals[[3]]),
                                           as.integer(length(intervals[[1]])))
   fdebug("Done croi_count_reads...")
   counts.croi[counts.croi==0]=1
   fdebug(sprintf("Counted %d reads...",libsize.croi))

#   print(length(counts))
#   print(length(counts.croi))
#   print(sum(counts==counts.croi))

#   d = counts != counts.croi
#   print(length(d))
#   print(sum(d))
#   counts_diff = counts[d]
#   counts_croi_diff = counts.croi[d]
#   stuff = cbind(intervals[d,],counts_diff,counts_croi_diff)
#   print(stuff[1:100,])

   counts = counts.croi
   libsize = libsize.croi

   widths = intervals[,3] - intervals[,2]
   rpkm = (counts/(widths/1000))/(libsize/1E6)
   
   return(list(counts=counts,rpkm=rpkm))
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

pv.addRowNamesToCounts = function(pv) {
   for(i in 1:length(pv$contrasts)) {
      rownames(pv$contrasts[[i]]$edgeR$counts) = 1:nrow(pv$contrasts[[i]]$edgeR$counts)       
   }
   return(pv)
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
   occ[,1] = as.character(occ[,1])
   occ[,2] = as.character(occ[,2])
   els = unique(c(occ[,1],occ[,2]))
   
   nels = length(els)
   res = matrix(NA,nels,nels)

   for(i in 1:nels) {
      for(j in 1:nels) {
         if(i==j) {
            res[i,j] = 1
         } else {
           entry = (occ[,1] == els[i]) & (occ[,2] == els[j])
           entry = entry | (occ[,1] == els[j]) & (occ[,2] == els[i])
           if(sum(entry)>0) {
              res[i,j] = as.numeric(levels(occ[,col])[as.numeric(occ[entry,col])])
              if(!missing(div)) {
              	 if(div != PV_TOTAL) {
              	    todiv = as.numeric(levels(occ[,div])[as.numeric(occ[entry,div])])
              	 } else {
              	    todiv = as.numeric(levels(occ[,PV_ONLYA])[as.numeric(occ[entry,PV_ONLYA])]) +
              	            as.numeric(levels(occ[,PV_ONLYB])[as.numeric(occ[entry,PV_ONLYB])]) +
              	            as.numeric(levels(occ[,PV_INALL])[as.numeric(occ[entry,PV_INALL])])
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


pv.getoneorNA = function(vec) {
 res = unique(vec)
 if (length(res) == 1) {
    return(res)
 } else {
   return(NA)
 }
}

pv.peaks2RangedData = function(peaks) {

   if(class(peaks)!="RangedData") {
      #require(IRanges)
      cnames = colnames(peaks)
      cnames[1:3] = c("space","start","end")
      colnames(peaks) = cnames
   }
   
   res = as(peaks,"RangedData")
   
   return(res)
   
}

pv.RangedData2Peaks = function(RDpeaks){

   if(class(RDpeaks)=="RangedData") {
      #require(IRanges)
      res = as.data.frame(RDpeaks)[,-4]
      cnames = colnames(res)
      cnames[1:3] = c("CHR","START","END")
      colnames(res) = cnames
   } else {
      res = RDpeaks
   }
   return(res)
}



