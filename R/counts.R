##########€€€€€######################
## pv_counts.R -- count-dependant  ##
## 20 October 2009                 ##
## 3 February 2011 -- packaged     ##
## Rory Stark                      ##
## Cancer Research UK              ##
#####################################
PV_DEBUG = FALSE
## pv.model -- build model, e.g. from sample sheet
pv.model = function(model,mask,minOverlap=2,
                    samplesheet='sampleSheet.csv',config=data.frame(RunParallel=FALSE),
                    caller="raw",skipLines=0,bAddCallerConsensus=T,bRemoveM=T, bRemoveRandom=T,
                    bKeepAll=T,bAnalysis=T,attributes) {

   
   
   if(!missing(model)) {
   	  if(missing(attributes)) {
         if(is.null(model$attributes)) {   
            attributes = PV_ID
         } else {
            attributes = model$attributes
         }
      }
      config = model$config
      model = pv.vectors(model,mask=mask,minOverlap=minOverlap,
                          bKeepAll=bKeepAll,bAnalysis=bAnalysis,attributes=attributes)
      model$config = config
      return(model)
   }
   
   if(missing(attributes)) {   
      attributes = PV_ID
   }

   if(is.character(samplesheet)) {
      samples = read.table(samplesheet,sep=',',stringsAsFactors=F,header=T)
      if(is.null(samples$SampleID)){
         samples$SampleID = 1:nrow(samples)
      }
      if(is.null(samples$Tissue)){
         samples$Tissue = ""
      } 
      if(is.null(samples$Factor)){
         samples$Factor = ""
      }
      if(is.null(samples$Condition)){
         samples$Condition = ""
      }
      if(is.null(samples$Treatment)){
         samples$Treatment = ""
      }
      if(is.null(samples$Replicate)){
         samples$Replicate = ""
      }

      samples$SampleID[is.na(samples$SampleID)]=""
      samples$Tissue[is.na(samples$Tissue)]=""
      samples$Factor[is.na(samples$Factor)]=""
      samples$Condition[is.na(samples$Condition)]=""
      samples$Treatment[is.na(samples$Treatment)]=""
      samples$Replicate[is.na(samples$Replicate)]=""
   } else samples = samplesheet
   
   model = NULL
   if(is.character(config)) {
      if(!is.null(config)) {
         config  = data.frame(t(read.table(config,comment.char="#",row.names=1)),
                              stringsAsFactors=F)
      }
   }
   if(is.null(config$parallelPackage)){
      config$parallelPackage=DBA_PARALLEL_MULTICORE
   }
   model$config = config
   
   for(i in 1:nrow(samples)) {
   	  
   	  if(is.null(samples$PeakCaller[i])) {
   	     peakcaller  = caller
   	     normCol = 0
   	  } else if(is.na(samples$PeakCaller[i])) {
   	     peakcaller  = caller
   	     normCol = 0
      } else {
   	     peakcaller = as.character(samples$PeakCaller[i])
   	  }

   	  if(is.null(samples$ControlID[i])) {
   	     controlid  = ''
   	  } else if(is.na(samples$ControlID[i])) {
   	     controlid  = ''
      } else {
   	     controlid = as.character(samples$ControlID[i])
   	  }
   	  
     message(as.character(samples$SampleID[i]),' ',
         as.character(samples$Tissue[i]),' ',
         as.character(samples$Factor[i]),' ',
         as.character(samples$Condition[i]),' ',
         as.character(samples$Treatment[i]),' ',
         as.integer(samples$Replicate[i]),' ',peakcaller)
         
     model = pv.peakset(model,
                         peaks       = as.character(samples$Peaks[i]),
                         sampID      = as.character(samples$SampleID[i]),
                         tissue      = as.character(samples$Tissue[i]),
                         factor      = as.character(samples$Factor[i]),
                         condition   = as.character(samples$Condition[i]),
                         treatment   = as.character(samples$Treatment[i]),
                         consensus   = F,
                         peak.caller = peakcaller,
                         control     = controlid,
                         reads       = NA,
                         replicate   = as.integer(samples$Replicate[i]),
                         readBam     = as.character(samples$bamReads[i]),
                         controlBam  = as.character(samples$bamControl[i]),
                         bRemoveM=bRemoveM, bRemoveRandom=bRemoveRandom,skipLines=skipLines)
      }

   model$samples = samples

   if(bAddCallerConsensus){
      model = pv.add_consensus(model)
   }
   
    model = pv.vectors(model,mask=mask,minOverlap=minOverlap,
                       bKeepAll=bKeepAll,bAnalysis=bAnalysis,attributes=attributes) 

     
   return(model)
}


## pv.counts -- add peaksets with scores based on read counts
PV_RES_RPKM             = 1
PV_RES_RPKM_FOLD        = 2
PV_RES_READS            = 3
PV_RES_READS_FOLD       = 4
PV_RES_READS_MINUS      = 5
PV_SCORE_RPKM           = PV_RES_RPKM
PV_SCORE_RPKM_FOLD      = PV_RES_RPKM_FOLD
PV_SCORE_READS          = PV_RES_READS
PV_SCORE_READS_FOLD     = PV_RES_READS_FOLD
PV_SCORE_READS_MINUS    = PV_RES_READS_MINUS
PV_SCORE_TMM_MINUS_FULL       = 6
PV_SCORE_TMM_MINUS_EFFECTIVE  = 7
PV_SCORE_TMM_READS_FULL       = 8
PV_SCORE_TMM_READS_EFFECTIVE  = 9

pv.counts = function(pv,peaks,minOverlap=2,defaultScore=PV_SCORE_RPKM_FOLD,bLog=T,insertLength=0,
                     bOnlyCounts=T,bCalledMasks=T,minMaxval,
                     bParallel=F,bUseLast=F,bWithoutDupes=F) {
   
   pv = pv.check(pv)
   
   if(!missing(peaks)) {
      if(is.character(peaks[1,1])){
         tmp = pv.peakset(NULL,peaks)
         pv$chrmap = tmp$chrmap
         peaks = tmp$peaks[[1]]
      } else {
         pv$chrmap = unique(as.character(peaks[,1]))
      }
      colnames(peaks)[1:3] = c("CHR","START","END")
      bed = pv.dovectors(peaks[,1:3],bKeepAll=T)
   } else {
     if(minOverlap == pv$minOverlap) {
        bed = pv$vectors[,1:3]
     } else if (minOverlap == 1) {
        bed = pv$allvectors[,1:3]
     } else {
        bed = pv.consensus(pv,1:length(pv$peaks),
                           minOverlap=minOverlap,bFast=T)$peaks[[length(pv$peaks)+1]][,1:3]
     }
   }
   
   bed[,1] = pv$chrmap[bed[,1]]
   bed = pv.peaksort(bed)

   numChips = ncol(pv$class)
   chips  = unique(pv$class[PV_BAMREADS,])
   inputs = pv$class[PV_BAMCONTROL,]
   inputs = unique(inputs[!is.na(inputs)])
   todo   = unique(c(chips,inputs))
  
   if(!pv.checkExists(todo)) {
      stop('Some read files could not be accessed. See warnings for details.')
   }
   
   if(!bUseLast) {
      pv = dba.parallel(pv)
      if((pv$config$parallelPackage>0) && bParallel) {   	     
   	     params  = dba.parallel.params(pv$config,c("pv.getCounts","pv.bamReads","pv.BAMstats","fdebug"))            
         results = dba.parallel.lapply(pv$config,params,todo,
                                       pv.getCounts,bed,insertLength,bWithoutDupes=bWithoutDupes)
      } else {
         results = NULL
         for(job in todo) {
      	    message('Sample: ',job)
            results = pv.listadd(results,pv.getCounts(job,bed,insertLength,bWithoutDupes=bWithoutDupes))
         }	
      }
      if(PV_DEBUG){
         #save(results,file='dba_last_result.RData')
      }
   } else {
   	  if(PV_DEBUG) {
         load('dba_last_result.RData')
      } else {
         warning("Can't load last result: debug off")
      }
   }
   
   if ((defaultScore >= DBA_SCORE_TMM_MINUS_FULL) || (defaultScore <= DBA_SCORE_TMM_READS_EFFECTIVE) ) {
      redoScore = defaultScore
      defaultScore = PV_SCORE_READS_MINUS	
   } else redoScore = 0
      
   allchips = unique(pv$class[c(PV_BAMREADS,PV_BAMCONTROL),])
   numAdded = 0
   for(chipnum in 1:numChips) {
   	  if (pv.nodup(pv,chipnum)) {
      	 jnum = which(todo %in% pv$class[PV_BAMREADS,chipnum])
   	     cond = results[[jnum]]
   	     if(length(cond$counts)==0){
   	        warning('ERROR IN PROCESSING ',todo[jnum])
   	     }
 	  
   	     if(!is.na(pv$class[PV_BAMCONTROL,chipnum])) {
   	        cnum = which(todo %in% pv$class[PV_BAMCONTROL,chipnum])
   	        cont = results[[cnum]]
   	        if(length(cont$counts)==0){
   	           warning('ERROR IN PROCESSING ',todo[cnum])
   	        }
   	     } else {
   	  	    cont = NULL
   	        cont$counts = rep(1,length(cond$counts))
   	        cont$rpkm   = rep(1,length(cond$rpkm))   
   	     }
   	     	  
   	     rpkm_fold   = cond$rpkm   / cont$rpkm
   	     reads_fold  = cond$counts / cont$counts
   	     reads_minus = cond$counts - cont$counts
  	  
         if(bLog) {
            rpkm_fold  = log2(rpkm_fold)
            reads_fold = log2(reads_fold)
         }
         if(defaultScore == PV_RES_RPKM) {
            scores = cond$rpkm
         } else if (defaultScore == PV_RES_RPKM_FOLD ) {
            scores = rpkm_fold
         } else if (defaultScore == PV_RES_READS) {
            scores = cond$counts    
         } else if (defaultScore == PV_RES_READS_FOLD) {
            scores = reads_fold
         } else if (defaultScore == PV_RES_READS_MINUS) {
            scores = reads_minus
         }
      
         res = cbind(bed,scores,cond$rpkm,cond$counts,cont$rpkm,cont$counts)
         colnames(res) = c("Chr","Start","End","Score","RPKM","Reads","cRPKM","cReads")
         pv = pv.peakset(pv,
                         peaks       = res,
                         sampID      = pv$class[PV_ID,chipnum],
                         tissue      = pv$class[PV_TISSUE,chipnum],
                         factor      = pv$class[PV_FACTOR,chipnum],
                         condition   = pv$class[PV_CONDITION,chipnum],
                         treatment   = pv$class[PV_TREATMENT,chipnum],
                         consensus   = T,
                         peak.caller = 'counts',
                         control     = pv$class[PV_CONTROL,chipnum],
                         reads       = cond$libsize, #pv$class[PV_READS,chipnum],
                         replicate   = pv$class[PV_REPLICATE,chipnum],
                         readBam     = pv$class[PV_BAMREADS,chipnum],
                         controlBam  = pv$class[PV_BAMCONTROL,chipnum],
                         bNormCol    = 0,
                         bRemoveM = F, bRemoveRandom=F,bMakeMasks=F)
         numAdded = numAdded + 1
      }                  
   }
      
   if(bOnlyCounts) {
   	  numpeaks = length(pv$peaks)
      res = pv.vectors(pv,(numpeaks-numAdded+1):numpeaks,minOverlap=1,bAnalysis=F,bAllSame=T)
      if(redoScore > 0) {
         res = pv.setScore(res,redoScore)	
      }   
      if(!missing(minMaxval)) {
         data = res$allvectors[,4:ncol(res$allvectors)]
         maxs = apply(res$allvectors[,4:ncol(res$allvectors)],1,max)
         tokeep = maxs>=minMaxval
         if(sum(tokeep)>1) {
            res$allvectors = res$allvectors[tokeep,]
            rownames(res$allvectors) = 1:sum(tokeep)
            res$vectors    = res$allvectors
            for(i in 1:length(res$peaks)) {
               res$peaks[[i]] = res$peaks[[i]][tokeep,]
               rownames(res$peaks[[i]]) = 1:sum(tokeep)
            } 
            res = pv.vectors(res,minOverlap=1,bAnalysis=F,bAllSame=T)
         } else {
            stop('No sites have activity greater than minMaxval')
         }
      }
      if(bCalledMasks && missing(peaks)) {
         res$sites = pv.CalledMasks(pv,res,bed)
      }
   } else {
      if(redoScore > 0) {
         res = pv.setScore(res,redoScore,minMaxval=minMaxval)	
      } 
      res = pv.vectors(pv)   
   }
   
   return(res)	
}

pv.nodup = function(pv,chipnum) {

   if(chipnum == 1) {
      return(TRUE)
   }
   
   chips = pv$class[PV_BAMREADS,1:(chipnum-1)] == pv$class[PV_BAMREADS,chipnum]
   conts = pv$class[PV_BAMCONTROL,1:(chipnum-1)] == pv$class[PV_BAMCONTROL,chipnum]
   
   if(sum(chips&conts)>0) {
      return(FALSE)
   } else {
      return(TRUE)
   }

}

pv.checkExists = function(filelist){
   res = file.access(filelist,mode=4)
   for(i in 1:length(filelist)) {
      if(res[i]==-1) {
         warning(filelist[i]," not accessible",call.=FALSE)	
      }	
   }
   return(sum(res)==0)
}

pv.getCounts = function(bamfile,intervals,insertLength=0,bWithoutDupes=F) {

   fdebug(sprintf('pv.getCounts: ENTER %s',bamfile))
   
   fdebug("Starting croi_load_reads...")
   bamtree <- .Call("croi_load_reads",as.character(bamfile),as.integer(insertLength))
   fdebug("Loaded...")
   libsize.croi <- .Call("croi_tree_size",bamtree)
   fdebug("Starting croi_count_reads...")
   counts.croi <- .Call("croi_count_reads",bamtree,
                                           as.character(intervals[[1]]),
                                           as.integer(intervals[[2]]),
                                           as.integer(intervals[[3]]),
                                           as.integer(length(intervals[[1]])),
                                           as.logical(bWithoutDupes))
   fdebug("Done croi_count_reads...")
   counts.croi[counts.croi==0]=1
   fdebug(sprintf("Counted %d reads...",libsize.croi))

   counts = counts.croi
   libsize = libsize.croi

   widths = intervals[,3] - intervals[,2]
   rpkm = (counts/(widths/1000))/(libsize/1E6)
   
   return(list(counts=counts,rpkm=rpkm,libsize=libsize))
}


