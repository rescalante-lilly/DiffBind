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
      model = pv.vectors(model,mask=mask,minOverlap=minOverlap,
                          bKeepAll=bKeepAll,bAnalysis=bAnalysis,attributes=attributes)
      return(model)
   }
   
   if(missing(attributes)) {   
      attributes = PV_ID
   }

   if(is.character(samplesheet)) {
      samples = read.table(samplesheet,sep=',',stringsAsFactors=F,header=T)
      samples$SampleID[is.na(samples$SampleID)]=""
      samples$Tissue[is.na(samples$Tissue)]=""
      samples$Factor[is.na(samples$Factor)]=""
      samples$Condition[is.na(samples$Condition)]=""
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
         as.integer(samples$Replicate[i]),' ',peakcaller)
         
     model = pv.peakset(model,
                         peaks       = as.character(samples$Peaks[i]),
                         sampID      = as.character(samples$SampleID[i]),
                         tissue      = as.character(samples$Tissue[i]),
                         factor      = as.character(samples$Factor[i]),
                         condition   = as.character(samples$Condition[i]),
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

pv.counts = function(pv,peaks,minOverlap=2,defaultScore=PV_RES_READS_MINUS,bLog=T,insertLength=0,
                     bOnlyCounts=T,bCalledMasks=T,minMaxval=0,
                     bParallel=F,bUseLast=F) {
   
   pv = pv.check(pv)
   
   if(!missing(peaks)) {
      if(is.character(peaks[1,1])){
         tmp = pv.peakset(NULL,peaks)
         pv$chrmap = tmp$chrmap
         peaks = tmp$peaks[[1]]
      }
      colnames(peaks)[1:3] = c("CHR","START","END")
      bed = pv.dovectors(peaks[,1:3],bKeepAll=T)
   } else {
     if(minOverlap == 2) {
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
                                       pv.getCounts,bed,insertLength)
      } else {
         results = NULL
         for(job in todo) {
      	    message('Sample: ',job)
            results = pv.listadd(results,pv.getCounts(job,bed,insertLength))
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
      if(minMaxval>0) {
         data = res$allvectors[,4:ncol(res$allvectors)]
         maxs = apply(res$allvectors[,4:ncol(res$allvectors)],1,max)
         tokeep = maxs>=minMaxval
         if(sum(tokeep)>1) {
            res$allvectors = res$allvectors[tokeep,]
            res$vectors    = res$allvectors
            for(i in 1:length(res$peaks)) {
               res$peaks[[i]] = res$peaks[[i]][tokeep,]
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
      res = pv.vectors(pv)   
   }   
   
   return(res)	
}


## pv.plotScatter -- scatter plots of scores for XY pairs
pv.plotScatter = function(pv,peaksets,overlaps,olmask,sites,attributes=pv$attributes,
                            bSmooth=T,CorMethod="pearson",...) {
  
   if(!missing(overlaps)) cres  = overlaps
   if(!missing(olmask))   crecs = olmask
   
   plotfun = plot
   if(bSmooth) {
      plotfun=smoothScatter
   }
   
   if(missing(peaksets)) {
      mask = rep(T,ncol(pv$class))
   } else {
      mask = rep(T,ncol(pv$class))
      mask[peaksets] = T
   }
   
   vals = pv.FixMin(pv$vectors)
   
   if(!missing(overlaps)){
   	  if(!missing(olmask)) {
   	     cres = cres[crecs,]
   	  }
   	  if(is.null(nrow(cres))){
   	     cres = matrix(cres,1,length(cres))	
   	  }
      for(i in 1:nrow(cres)){
         s1 = cres[i,1] + 3
         first = pv$class[attributes,s1-3]
         s2 = cres[i,2] + 3
         second = pv$class[attributes,s2-3]
   	     res = pv.namestrings(first,second)
         corval = cres[i,6]
         plotfun(vals[sites,s1],vals[sites,s2],pch=20,cex=.1,
                 xlab=sprintf("%i:[%s]",s1-3,res$n1),
                 ylab=sprintf("%i:[%s]",s2-3,res$n2),
                 main=sprintf("CONTRAST:%s",res$tstring),
                 sub=sprintf('R=%1.2f',corval),...)
      }
   } else {
      numSets = sum(mask)
      wmask = which(c(F,F,F,mask))
      for(i in 1:(numSets-1)) {
   	     s1 = wmask[i]
   	     first = pv$class[attributes,s1-3]
         for(j in (i+1):numSets){
      	    s2=wmask[j]
      	    second = pv$class[attributes,s2-3]
      	    res = pv.namestrings(first,second)
      	    corval = cor(vals[sites,s1],vals[sites,s2],method=CorMethod)
            plotfun(vals[sites,s1],vals[sites,s2],pch=20,cex=.1,
                 xlab=sprintf("%i:[%s]",s1,res$n1),
                 ylab=sprintf("%i:[%s]",s2,res$n2),
                 main=sprintf("CONTRAST:%s",res$tstring),
                 sub=sprintf('R=%1.2f',corval),...)
         }           
      }
   }   	
} 


## pv.compare -- differential binding affinity table
PV_GENOME_HUMAN = "H_sapiens_Mar_2006"
PV_GENOME_MOUSE = "M_musculus_Jul_2007"
pv.compare = function(pv,X,Y,pX=X,pY=Y,Xlab="X",Ylab="Y",sites,minval=-1,bPeaks=F,bCounts=F,
                      bNoControl=F,GenomeString=PV_GENOME_HUMAN,names,bPlot=F,
                      IntervalList,Colors=c(1,3,2,4:10),...) {
   	
   	if(!bPeaks) {
   	   if(missing(sites)) {
   	      mask = rep(T,nrow(pv$vectors))
       } else mask = sites
   	
   	   DB  = pv$vectors[mask,X+3] - pv$vectors[mask,Y+3]
   	   DBX = pv$vectors[mask,pX+3] > minval
   	   DBY = pv$vectors[mask,pY+3] > minval
       DBmean = (pv$vectors[mask,X+3] + pv$vectors[mask,Y+3]) / 2

       if(missing(GenomeString)) {
   	      sites = apply(pv$vectors[mask,1:3],1,pv.dositename_map,pv$chrmap)
   	   } else {
   	      sites = apply(pv$vectors[mask,1:3],1,pv.dositename_map,pv$chrmap,GenomeString)    
       } 
       
       o = order(abs(DB),decreasing=T)
   	   res = cbind(sites,DB,abs(DB),DBmean,pv$vectors[mask,X+3],pv$vectors[mask,Y+3])    
   	   ecor = cor(pv$vectors[mask,X+3],pv$vectors[mask,Y+3])
   	   
    } else { # bPeaks==1
       
       if(missing(mask)) {
          mask = rep(T,nrow(pv$peaks[[X]]))
       }
   	   
       if(missing(GenomeString)) {
   	      sites = apply(pv$peaks[[X]][mask,1:3],1,pv.dositename)
   	   } else {
   	      sites = apply(pv$peaks[[X]][mask,1:3],1,pv.dositename,GenomeString)    
       } 
       
       if(bCounts) {	
   	      Xdata = pv$peaks[[X]][mask,]	
   	      Ydata = pv$peaks[[Y]][mask,]
   	      if(bNoControl) {
   	         control = 1
   	      } else {
   	      	 control = (Xdata[,6] + Ydata[,6])/2
   	         control[control==0]=1
   	      }
   	      Xreads = Xdata[,5]
   	      Yreads = Ydata[,5]   	   
   	      Xreads[Xreads==0]=1
   	      Yreads[Yreads==0]=1
   	      Xfold = log((Xreads / control),2)
   	      Yfold = log((Yreads / control),2)
   	      Xfold[Xfold<0]=0
   	      Yfold[Yfold<0]=0
   	      #cat("Correlations - ChiP:",cor(Xreads,Yreads),"Control:",cor(Xdata[,6], Ydata[,6]),"\n")
   	      #cat(Xlab,"Reads in peaks/control:",sum(Xdata[,5]),"/",sum(Xdata[,6]),'\n')                  
   	      #cat(Ylab,"Reads in peaks/control:",sum(Ydata[,5]),"/",sum(Ydata[,6]),'\n') 
   	   } else {
   	      Xfold = as.numeric(pv$peaks[[X]][mask,4])
   	      Yfold = as.numeric(pv$peaks[[Y]][mask,4])
   	   }
   	   
   	   if(pX==X) {
   	      DBX = Xfold > minval
   	   } else {
   	      DBX = pv$vectors[mask,pX+3] > minval
   	   }
   	   if(pY==Y) {
   	      DBY = Yfold > minval
   	   } else {
   	      DBY = pv$vectors[mask,pY+3] > minval	
   	   }
   	   Cdb = Xfold - Yfold
   	   Cmean = (Xfold+Yfold)/2
   	   res = cbind(sites,Cdb,abs(Cdb),Cmean,Xfold,Yfold)
   	   o = order(abs(Cdb),decreasing=T)
   	   ecor = cor(Xfold,Yfold)
   	} 
   	
   	res = cbind(res,DBX,DBY,DBX & DBY) 	                             
   	colnames(res) = c("Site_loc","log_fold_diff","abs(log_fold_diff)","Mean_log_fold",
   	                  sprintf("%s_log_fold",Xlab),sprintf("%s_log_fold",Ylab),
   	                  sprintf("%s_peak?",Xlab),sprintf("%s_peak?",Ylab),"Both_Peak?")
   	
   	#cat('Correlation:',ecor,'\n')
   	   	   
   	if(bPlot) {
       if(missing(IntervalList)) {
          IntervalList = list(DBX & DBY,DBX & (!DBY),DBY & (!DBX))
       }
       maxval = ceiling(as.numeric(max(res[,5:6])))
       plot(res[IntervalList[[1]],5],res[IntervalList[[1]],6],xlim=c(0,maxval),ylim=c(0,maxval),
            pch=20,cex=.1,col=Colors[1],sub=sprintf("R=%1.2f",ecor),xlab=Xlab,ylab=Ylab,...)
       for(i in 2:length(IntervalList)) {
          points(res[IntervalList[[i]],5],res[IntervalList[[i]],6],pch=20,cex=.1,col=Colors[i])
       }
    }
   	
   	if(!missing(names)) {
   	   cnames= colnames(res)
   	   res = cbind(res,names)
   	   colnames(res) = c(cnames,"NAME")	
   	}
   	return(data.frame(res[o,]))   	 
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
