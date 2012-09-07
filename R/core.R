#####################################
## pv_core.R -- Peak Vectorization ##
## 20 October 2009                 ##
## 3 February 2011 -- packaged     ##
## Rory Stark                      ##
## Cancer Research UK              ##
#####################################

############################
## Top Level Entry Points ##
############################
## pv.peakset       -- add a peakset w/ attributes to the model
## pv.vectors       -- build the binding expression vectors and do clustering/PCA analysis
## pv.list          -- list peaksets w/ attributes in model
## pv.consensus     -- add a consensus peakset based on peaksets already in model

## pv.mask          -- create a mask to define a subset of peaksets in a model
## pv.whichSites    -- create a mask of sites belonging to specific peakset(s)

## pv.plotClust     -- hierarchical cluster plot
## pv.plotPCA       -- interactive 3D PCA plot
## pv.plotHeatmap   -- plot a binding site heatmap w/ hierarchical clustering 
## pv.sort          -- sort binding sites (e.g. for heatmap)

## pv.overlap       -- generate overlapping/unique peaksets
## pv.plotVenn      -- draw venn diagrams

## pv.occupancy     -- generate occupancy statistics for peaksets in a model 
## pv.plotScatter   -- scatter plots of contrasts


################
## Constants  ##
################
PV_GROUP      = 0
PV_ID         = 1
PV_TISSUE     = 2
PV_FACTOR     = 3
PV_CONDITION  = 4
PV_CONSENSUS  = 5
PV_CALLER     = 6
PV_CONTROL    = 7
PV_READS      = 8
PV_REPLICATE  = 9
PV_BAMREADS   = 10
PV_BAMCONTROL = 11
PV_TREATMENT  = 12

PV_DEBUG = FALSE

###########################################

## pv.peakset -- add a peakset to the model
pv.peakset = function(pv=NULL,peaks, sampID, tissue, factor,condition, treatment, replicate,
                      control, peak.caller, peak.format, reads=0, consensus=F, readBam, controlBam,
                      scoreCol=NULL, bLowerScoreBetter=NULL, bRemoveM=T, bRemoveRandom=T,
                      minOverlap=2,bFast=F,bMakeMasks=T,skipLines=1){
	
   zeroVal = -1
   bLog=F
     
   if(missing(peaks)) {
     peaks = 1:length(pv$peaks)
   }
     
   if(missing(peak.format))       peak.format=NULL
   if(missing(scoreCol))          scoreCol=NULL
   if(missing(bLowerScoreBetter)) bLowerScoreBetter=NULL   
   
   bConsensus = F
   if(is.numeric(consensus)) { ## Add a set of consensus peaksets
   	  bConsensus = T
      pv = pv.consensusSets(pv,peaks=peaks,minOverlap=minOverlap,attributes=consensus,
                            tissue,factor,condition,treatment,replicate,control,peak.caller,
                            readBam, controlBam)
                            
   } else { ## add a specific consensus peakset
	   if(is.vector(peaks) && length(peaks) > 1) { # consensus
	   	  bConsensus = T
	      pv = pv.consensus(pv,peaks,minOverlap=minOverlap,bFast=bFast)
	      if(!is.null(minOverlap)) {
	
		      nset = length(pv$peaks)
		      if(!missing(sampID)){
		         pv$class[PV_ID,nset]=sampID
		         colnames(pv$class)[nset]=sampID
		      }
		  }
	
	      if(!missing(tissue))      pv$class[PV_TISSUE,nset]=tissue
	      if(!missing(factor))      pv$class[PV_FACTOR,nset]=factor
	      if(!missing(condition))   pv$class[PV_CONDITION,nset]=condition
	      if(!missing(treatment))   pv$class[PV_TREATMENT,nset]=treatment
	      if(!missing(replicate))   pv$class[PV_REPLICATE,nset]=replicate
	      if(!missing(control))     pv$class[PV_CONTROL,nset]=control
	      if(!missing(peak.caller)) pv$class[PV_CALLER,nset]=peak.caller
	      if(!missing(readBam))     pv$class[PV_BAMREADS,nset]=readBam
	      if(!missing(controlBam))  pv$class[PV_BAMCONTROL,nset]=controlBam
	   }
   }
   if(bConsensus) {      
      if(bMakeMasks) {
         pv$masks = pv.mask(pv)
      }
      return(pv)
   }
   
   if(missing(tissue))       tissue=''
   if(missing(factor))       factor=''
   if(missing(condition))    condition=''
   if(missing(treatment))    treatment=''
   if(missing(replicate))    replicate=''
   if(missing(control))      control=''
   if(missing(peak.caller))  peak.caller=''
   if(missing(readBam))      readBam=NA
   if(length(readBam)==0)    readBam=NA
   if(missing(controlBam))   controlBam=NA
   if(length(controlBam)==0) controlBam=NA
   
                  
   if(is.character(peaks)){ # Read in peaks from a file
     
      pcaller = strtrim(peak.caller,6)
      if(is.null(peak.format)) {
          peak.format = pcaller
      }	
      if(is.null(scoreCol)) {
         scoreCol = pv.defaultScoreCol(peak.format)	
      }
      if(is.null(bLowerScoreBetter)) bLowerScoreBetter = FALSE
      
      peaks = pv.readPeaks(peaks,peak.format,skipLines)
   }
   
   if(ncol(peaks) < scoreCol) {
      peaks = cbind(peaks=1)
      scoreCol = 0
   }
      
   if(scoreCol > 0) {
      peaks[,scoreCol] = pv.normalize(peaks,scoreCol,zeroVal=zeroVal,bLog=bLog)
      if(bLowerScoreBetter) {
         peaks[,scoreCol] = 1 - peaks[,scoreCol]   	
      }
   }      
      
   if(bRemoveM) {
      idx = peaks[,1] != "chrM"
      peaks = peaks[idx,]
      if(sum(!idx)>0) {
         peaks[,1] = as.factor(as.character(peaks[,1]))
      }
   }

   if(bRemoveRandom) {
   	  for(i in c(1:22,"X","Y")) {
   	     ch = sprintf("chr%s_random",i)
         idx = peaks[,1] != ch
         peaks = peaks[idx,]
         if(sum(!idx)>0) {
            peaks[,1] = as.factor(as.character(peaks[,1]))
         }
      }
   }
   
   newchrs = as.character(peaks[,1])
   pv$chrmap  = unique(c(pv$chrmap,newchrs))
   peaks[,1] = factor(peaks[,1],pv$chrmap)

   pv$peaks = pv.listadd(pv$peaks,peaks)
   
   if(missing(sampID)) {
      if(is.null(pv)) {
         sampID=1
      } else if (is.null(pv$peaks)) {
         sampID=1
      } else {
         sampID = length(pv$peaks)
      }
   }
   clascol = cbind(NULL,c(sampID,tissue,factor,condition,consensus,peak.caller,control,
                          reads,replicate,readBam,controlBam,treatment))
   colnames(clascol) = sampID
   pv$class = cbind(pv$class,clascol)
   rownames(pv$class) = c("ID","Tissue","Factor","Condition", "Consensus",
                          "Peak caller","Control","Reads","Replicate","bamRead",
                          "bamControl","Treatment")
   pv$vectors = pv$allvectors = NULL                            
   if(bMakeMasks) {
      pv$masks = pv.mask(pv)
   }
   return(pv)	
}



pv.peakset_all = function(pv, addpv, minOverlap) {

   for(i in 1:length(addpv$peaks)) {
   	
   	  message(addpv$class[PV_ID,i],' ',
              addpv$class[PV_TISSUE,i],' ',
              addpv$class[PV_FACTOR,i],' ',
              addpv$class[PV_CONDITION,i],' ',
              addpv$class[PV_REPLICATE,i],' ',
              addpv$class[PV_CALLER,i])
    
      pv = pv.peakset(pv,peaks=addpv$peaks[[i]],
                      sampID      = addpv$class[PV_ID,i],
                      tissue      = addpv$class[PV_TISSUE,i],
                      factor      = addpv$class[PV_FACTOR,i],
                      condition   = addpv$class[PV_CONDITION,i],
                      replicate   = addpv$class[PV_REPLICATE,i],
                      control     = addpv$class[PV_CONTROL,i],
                      peak.caller = addpv$class[PV_CALLER,i],
                      reads       = addpv$class[PV_READS,i],
                      consensus   = addpv$class[PV_CONSENSUS,i],
                      readBam     = addpv$class[PV_BAMREADS,i],
                      controlBam  = addpv$class[PV_BAMCONTROL,i]
                     )
   }
   
   pv = dba(pv, minOverlap=minOverlap)

   return(pv)

}

## pv.vectors -- build the binding expression vectors and do clustering/PCA
pv.vectors = function(pv,mask,minOverlap=2,bKeepAll=T,bAnalysis=T,attributes,bAllSame=F){


   if(missing(attributes)) {
      if(is.null(pv$attributes)) {   
         attributes = PV_ID
      } else {
         attributes = pv$attributes
      }
   }
        
   if(!missing(mask)){
      if(is.logical(mask)) {
         mask = which(mask)
      }
   	  peaks = NULL
      for(i in mask){
         peaks = pv.listadd(peaks,pv$peaks[[i]])
      }
      class     = pv$class[,mask]
      chrmap    = pv$chrmap
      config    = pv$config
      samples   = pv$samples
      contrasts = pv$contrasts
      
      pv = NULL
      pv$peaks     = peaks
      pv$class     = class
      pv$chrmap    = chrmap 
      pv$config    = config
      pv$samples   = samples 
      pv$contrasts = contrasts  
   } 
   
   if(is.vector(pv$class)) {
      pv$class = matrix(pv$class,length(pv$class),1)
      colnames(pv$class) = pv$class[1,]
   }
   
   peaks = pv$peaks
   
   numvecs = length(peaks)
   ncols = numvecs+3
   
   npeaks=0
   defval = -1

   if(!bAllSame) {
      for(pnum in 1:numvecs){
         npeaks = npeaks + nrow(peaks[[pnum]])
      }
      allpeaks = matrix(defval,npeaks,numvecs+3)
      st = 1
      end = 0
      for(pnum in 1:numvecs){
 
         peak = peaks[[pnum]]
         end = end + nrow(peak)
      
         allpeaks[st:end,(pnum+3)] = peak[,4]

         allpeaks[st:end,1] = match(as.character(peak[,1]),pv$chrmap)
         allpeaks[st:end,2] = peak[,2]
         allpeaks[st:end,3] = peak[,3]
             
         st = end + 1
      }
      colnames(allpeaks) = c("CHR","START","END",1:numvecs)

      if(bKeepAll) {
   	     #result = pv.dovectors(allpeaks,pv$class,bKeepAll=F,bDel=F,bUseLast=F)
         pv$allvectors  = pv.dovectors(allpeaks,pv$class,bKeepAll=T)
         rownames(pv$allvectors) = 1:nrow(pv$allvectors)
         if((ncol(pv$allvectors)>4) && (minOverlap>1)) {
            pv$overlapping = apply(pv$allvectors[,4:ncol(pv$allvectors)],1,pv.minOverlap,minOverlap)
            result = pv$allvectors[pv$overlapping,]
         } else {
            result = pv$allvectors
         }
      } else {
         result = pv.dovectors(allpeaks,pv$class,bKeepAll=F)
      }
   
      pv$vectors = result
   
   } else { ## ALL SAME
      pv$allvectors = pv$peaks[[1]][,1:4]

      for(i in 2:numvecs){
         pv$allvectors = cbind(pv$allvectors,pv$peaks[[i]][,4])
      }	
      colnames(pv$allvectors) = c("CHR","START","END",pv$class[PV_ID,1:numvecs])
      pv$vectors = pv$allvectors
   }
   
   pv$attributes = attributes
   pv$minOverlap = minOverlap
   
   if(nrow(pv$vectors)>0) {
      rownames(pv$vectors) = 1:nrow(pv$vectors)
   }
   
   if(bAnalysis && ncol(pv$class)>1) {
      #pv = pv.analysis(pv)
   }
  pv$hc = NULL
  pv$pc = NULL
  pv$masks = pv.mask(pv)
  return(pv)   	
}

## pv.list -- list attributes of samples in model
pv.deflist = c(PV_ID,PV_TISSUE,PV_FACTOR,PV_CONDITION,PV_TREATMENT,PV_REPLICATE,PV_CALLER)
pv.list = function(pv,mask,bContrasts=F,attributes=pv.deflist,th=0.1,bUsePval=F){
 
   if(!missing(mask)){
      if(!is.logical(mask)) {
         tmp  = rep (F,length(pv$peaks))
         tmp[mask] = T
         mask = tmp
      }
   }
   
   if(bContrasts) {
      return(pv.listContrasts(pv,th=th,bUsePval=bUsePval))
   }
   
   if(missing(attributes)) {
      attributes = pv.deflist
   }
      
   if(missing(mask)){
      mask = rep(T,ncol(pv$class))
   }
   
   if(!is.logical(mask)) {
      newm = rep(F,length(pv$peaks))
      for(ps in mask) {
         newm[ps]=T
      }
      mask = newm
   }
   
   res = t(pv$class[attributes,mask])
   rownames(res) = which(mask)
   
   intervals = NULL
   for(i in 1:length(mask)) {
      if(mask[i]) {
         intervals = c(intervals,nrow(pv$peaks[[i]]))
      }	
   }
   res = cbind(res,intervals)
   colnames(res)[ncol(res)]='Intervals'
   
   j = ncol(res)
   for(i in j:1) {
      x = unique(res[,i])
      if(colnames(res)[i]=='Peak caller') {
         if(all.equal(attributes,pv.deflist) == TRUE) {
            if(length(x)==1) {
               res = res[,-i]	
            }
         }
      } else {
         if(length(x)==1 && x[1]=="") {
            res = res[,-i]	
         }    	
      }
   }
   
   return(data.frame(res))
}

## pv.consensus -- add a consensus peakset based on peaksets already in model
pv.consensus = function(pv,sampvec,minOverlap=2,bFast=F,sampID){
   
   
   if(missing(sampvec)) {
      sampvec = 1:length(pv$peaks)
   }
   if(is.null(sampvec)) {
      sampvec = 1:length(pv$peaks)
   }
   if(class(sampvec)=='logical') {
      sampvec = which(sampvec)
   }

   tmp = NULL
   if((bFast | is.null(minOverlap)) & (max(sampvec)<=ncol(pv$class)))  {
   	  if(length(sampvec)<length(pv$peaks)) {
   	     pv = pv.vectors(pv,sampvec,minOverlap=1,bAnalysis=F)
   	     sampvec = 1:length(pv$peaks)
   	  } else {
         pv = pv.check(pv)
      }
      sites = pv.whichSites(pv,sampvec,bUseAllvecs=T)
      tmp$vectors = pv$allvectors[sites,c(1:3,(sampvec+3))]
      if(is.null(minOverlap)) {
         res = rep(0,length(sampvec))
         counts = apply(tmp$vectors[,4:ncol(tmp$vectors)],1,pv.countOverlap,-1)
         for(i in 1:length(sampvec)){
            res[i] = sum(counts>=i)
         }
         return(res)
      }      
   } else {
      peaklist  = NULL
      classlist = NULL
      for(samp in sampvec){
         peaklist  = pv.listadd(peaklist,pv$peaks[[samp]])
      }
      tmp$peaks  = peaklist
      tmp$class  = pv$class[, sampvec]
      tmp$chrmap = pv$chrmap
   
      if(is.null(minOverlap)) { ## SPECIAL CASE: OVERLAP RATES
   	     tmp = pv.vectors(tmp,bKeepAll=T,bAnalysis=F)
         res = rep(0,length(sampvec))
         for(i in 1:length(sampvec)){
            res[i] = sum(apply(tmp$allvectors[,4:ncol(tmp$vectors)],1,pv.minOverlap,i))
         }
         return(res)
      } else if(minOverlap == 1) {
         tmp = pv.vectors(tmp,bKeepAll=T,bAnalysis=F)
         tmp$vectors = tmp$allvectors
      } else {
         tmp = pv.vectors(tmp,bKeepAll=F,bAnalysis=F)
      }
   }
   
   if((minOverlap != 2) | bFast) { 
      goodvecs = apply(tmp$vectors[,4:ncol(tmp$vectors)],1,pv.minOverlap,minOverlap)
   } else {
      goodvecs = rep(T,nrow(tmp$vectors))
   }
   tmp$vectors  = tmp$vectors[goodvecs,]
   
   mean.density = apply(tmp$vectors[,4:ncol(tmp$vectors)],1,pv.domean)
   tmp$vectors = cbind(tmp$vectors[,1:3],mean.density)
   
   #kludge to get peakset in correct format
   tmpf = tempfile(as.character(Sys.getpid())) #tmpdir='.')
   pv.do_peaks2bed(tmp$vectors, pv$chrmap,tmpf)
   tmp$vectors = pv.readbed(tmpf)
   unlink(tmpf)
   
   if (length(unique(pv$class[PV_REPLICATE, sampvec]))==1) {
     replicate = unique(pv$class[PV_REPLICATE, sampvec])
   } else {
      replicate = pv.catstr(pv$class[PV_REPLICATE, sampvec])
   }
   if(missing(sampID)) {
   	 sampID = pv.catstr(pv$class[PV_ID, sampvec])
   }
   pv = pv.peakset(pv,peaks=tmp$vectors,
                   sampID      = sampID,
                   tissue      = pv.catstr(pv$class[PV_TISSUE, sampvec]),
                   factor      = pv.catstr(pv$class[PV_FACTOR, sampvec]),
                   condition   = pv.catstr(pv$class[PV_CONDITION, sampvec]),
                   treatment   = pv.catstr(pv$class[PV_TREATMENT, sampvec]),
                   peak.caller = pv.catstr(pv$class[PV_CALLER, sampvec]),
                   control     = pv.catstr(pv$class[PV_CONTROL, sampvec]),
                   reads       = mean(as.numeric(pv$class[PV_READS, sampvec])),
                   replicate   = replicate,
                   consensus   = T,
                   readBam     = pv.getoneorNA(pv$class[PV_BAMREADS, sampvec]),
                   controlBam  = pv.getoneorNA(pv$class[PV_BAMCONTROL, sampvec]),
                   scoreCol    = 0)
                        
  pv$vectors    = NULL
  pv$allvectors = NULL
  pv$hc         = NULL
  pv$pc         = NULL

  return(pv)
}  

pv.consensusSets = function(pv,peaks=NULL,minOverlap,attributes,
                            tissue,factor,condition,treatment,replicate,control,peak.caller,
                            readBam, controlBam)	{

   if(is.character(peaks)) {
      stop("\"peaks\" parameter can not be a filename when \"consensus\" specifies attributes",call.=FALSE)	
   }

   if(is.null(peaks)) {
     peaks = rep(T,ncol(pv$class))
   }
   
   include = F
   exclude = F
   if(sum(attributes<0)) exclude = T
   if(sum(attributes>0)) include = T       
   
   if(include & exclude) {
      stop('Consensus attributes must be all inclusive (positive) or all exclusive (negative)',call.=F)	
   }
   
   if(exclude) {
      atts = NULL
      if(!(-PV_TISSUE    %in% attributes)) atts = c(atts,PV_TISSUE)  	
      if(!(-PV_FACTOR    %in% attributes)) atts = c(atts,PV_FACTOR)  	
      if(!(-PV_CONDITION %in% attributes)) atts = c(atts,PV_CONDITION)  	
      if(!(-PV_TREATMENT %in% attributes)) atts = c(atts,PV_TREATMENT) 
      if(!(-PV_REPLICATE %in% attributes)) atts = c(atts,PV_REPLICATE)
      if(!(-PV_CALLER    %in% attributes)) atts = c(atts,PV_CALLER) 
      attributes = atts
   }
   
   numatts = length(attributes)
   class = pv$class[attributes,]
   if(is.vector(class)) {
      class = matrix(class,1,length(class))	
   }
   
   specs = unique(class[,peaks],MARGIN=2)
   if(is.vector(specs)) {
      specs = matrix(specs,1,length(specs))	
   }
   
   if(ncol(specs) == ncol(class)) {
      warning('All peaksets unique for specified attributes; no consensus peaksets added.',call.=F)	
      return(pv)
   }
   for(i in 1:ncol(specs)) {
   	  cand = class %in% specs[,i]
   	  if(is.vector(cand)) {
   	     cand = matrix(cand,numatts,ncol(class))
   	  }
      samples = apply(cand,MARGIN=2,function(x){sum(x) == numatts}) & peaks
      diffatts = apply(class,MARGIN=1,function(x){length(unique(x))>1})
      if(sum(samples)>1) {
      	 message('Add consensus: ',paste(specs[diffatts,i],collapse=" "))
         pv = pv.consensus(pv,samples,sampID=paste(specs[diffatts,i],collapse=":"),minOverlap=minOverlap)
         sampnum = ncol(pv$class)
         if(pv$class[PV_ID,sampnum]=="") pv$class[PV_ID,sampnum]="ALL"
         if(!missing(tissue))      pv$class[PV_TISSUE,sampnum]    = tissue
         if(!missing(factor))      pv$class[PV_FACTOR,sampnum]    = factor
         if(!missing(condition))   pv$class[PV_CONDITION,sampnum] = condition
         if(!missing(treatment))   pv$class[PV_TREATMENT,sampnum] = treatment
         if(!missing(replicate))   pv$class[PV_REPLICATE,sampnum] = replicate
         if(!missing(control))     pv$class[PV_CONTROL,sampnum]   = control
         if(!missing(peak.caller)) pv$class[PV_CALLER,sampnum]    = peak.caller
         if(!missing(readBam))     pv$class[PV_BAMREADS,sampnum]  = readBam
         if(!missing(controlBam))  pv$class[PV_BAMCONTROL,sampnum]= controlBam
      }   	
   }
  
   return(pv)   
}


## pv.mask -- create a mask to define a subset of peaksets in a model
pv.mask = function(pv,attribute,value,combine='or',mask,merge='or',bApply=F){
   
   numsamps = ncol(pv$class)
   
   if(missing(mask)) {
      if((merge =='or') | (merge=='and')) {
         mask = rep(F, numsamps)
      } else {
      	 mask = rep(T, numsamps)
      }
   }
   
   if(missing(attribute)) {
      masks = NULL
      for(att in c(PV_TISSUE,PV_FACTOR,PV_CONDITION,PV_TREATMENT,PV_CALLER)) {
         vals = unique(pv$class[att,])
         for (v in vals) {
         	res = list(x=pv.mask(pv,att,v))
         	names(res) = v
            masks = c(masks,res)
         }
      }
      res = list(x=pv.mask(pv,PV_CONSENSUS,T))
      if(sum(res[[1]])) {
      	 names(res) = "Consensus"
         masks = c(masks,res)       
      }
      reps = unique(pv$class[PV_REPLICATE,])
      reps = reps[!is.na(reps)]
      if(length(reps>1)) {
         for(rep in reps) {
            res = list(x=pv.mask(pv,PV_REPLICATE,rep))
            names(res) = sprintf("Replicate.%s",rep)
            masks = c(masks,res)
         }
      }
      
      masks$All  = rep(T,ncol(pv$class))
      masks$None = rep(F,ncol(pv$class))
     
      return(masks)  
   }
   
   curmask = NULL
   for(v in value) {
      newmask = pv$class[attribute,] == v
      if(is.null(curmask)) curmask = newmask
      if((combine =='or') | (combine=='nor')) {
         curmask = curmask | newmask
      }
      if((combine =='and') | (combine=='nand')) {
         curmask = curmask & newmask
      }
   }
   
   if((combine =='nor') | (combine=='nand')) {
      curmask = !curmask
   }
   
   if(merge == 'or') {
      mask = mask | curmask	
   }
   
   if(merge == 'and') {
      mask = mask & curmask	
   }

   if(merge == 'nor') {
      mask = !(mask | curmask)	
   }

   if(merge == 'nand') {
      mask = !(mask & curmask)	
   }           
   if(bApply) {
      pv$vectors = NULL
      pv$allvectors = NULL
      pv$class = pv$class[,mask]
      newpeaks = NULL
      for(i in 1:length(mask)) {
         if(mask[i]) {
            newpeaks = pv.listadd(newpeaks,pv$peaks[[i]])
         }
      }
      pv$peaks = newpeaks
      return(pv)
   } else {
      return(mask)
   }
}

## pv.whichSites -- return index vector of sites belonging to a specific peakset
pv.whichSites = function(pv,pnum,combine="or",minVal=-1,bUseAllvecs=F){
   
   if(bUseAllvecs) {
      vecs = pv$allvectors	
   } else {
      vecs = pv$vectors
   }
   if(length(pnum)==1) {
      res = vecs[,pnum+3]>minVal
   } else {
   	  res = vecs[,pnum[1]+3]>minVal
      for(p in 2:length(pnum)) {
         newvec = vecs[,pnum[p]+3]>minVal
         if(combine=='or') {
            res = res | newvec
         }
         if(combine=='and') {
            res = res & newvec
         }
         if(combine=='nor') {
            res = !(res | newvec)
         }
         if(combine=='nand') {
            res = !(res & newvec)
         }
      }	
   }
   return(res)
}

## pv.plotClust  -- hierarchical cluster plot
pv.plotClust = function(pv,mask,numSites,sites,attributes=pv$attributes,distMeth="pearson") {
   if(missing(mask)) {
      mask = rep(T,ncol(pv$class))
   }
   if(missing(sites)) {
      sites = rep(T,nrow(pv$vectors))
   }
   if(missing(numSites)) {
      numSites = length(sites)
   }
   pv$vectors = pv$vectors[sites,mask][1:numSites,]
   pv$class   = pv$class[,mask]
   pv         = pv.analysis(pv,attributes,bPCA=F,distMeth=distMeth)
   plot(pv$hc)
   return(pv$hc)
}

## pv.plotPCA -- 3D plot of PCA
pv.plotPCA = function(pv,attributes=PV_ID,second,third,fourth,size,mask,
                      numSites,sites,cor=F,startComp=1,b3D=T,vColors,...){
   
   pv = pv.check(pv)
   
   class  = attributes[1]   
   if(length(attributes)>1) {
      second = attributes[2]	
      if(length(attributes)>2) {
         third = attributes[3]
      }
      if(length(attributes)>3) {
         fourth = attributes[4]
      }
   }
   
   if(missing(sites)) sites = NULL

    if(missing(numSites)){
       numSites = nrow(pv$vectors)
    }
         	  
   if(!missing(mask) || !missing(numSites)){
   	  if(missing(mask)) {
   	     mask = rep(T,ncol(pv$class))
   	  }
      pv = pv.pcmask(pv,numSites,mask,sites,cor=cor)
   } 
   
   pc = pv$pc
   
   #if(!is.null(pv$mask)) {
   #   classes = pv$class[,which(pv$mask)]
   #} else {
      classes = pv$class
   #}
   
   if(max(class) > nrow(classes)){
      return(F)
   }
   
   vr = rep(0,length(pc$sdev))
   for(i in 1:length(vr)) { vr[i] = pc$sdev[i]^2 }
   
   if(b3D){
   	  startComp=1
      pvar = sum(vr[startComp:(startComp+2)])/sum(vr)*100
   } else {
      pvar = sum(vr[startComp:(startComp+1)])/sum(vr)*100
   }

   if(!missing(second)){
   	  if(!missing(third)) {
   	     if(!missing(fourth)) {
	     	classvec = sprintf("%s:%s:%s:%s",classes[class,],classes[second,],classes[third,],classes[fourth,])      
            thetitle = sprintf("PCA: %s:%s:%s:%s [%2.0f%% of total variance]",rownames(classes)[class],
                                               rownames(classes)[second],
                                               rownames(classes)[third],
                                               rownames(classes)[fourth],pvar)
   	     } else {
            classvec = sprintf("%s:%s:%s",classes[class,],classes[second,],classes[third,])      
            thetitle = sprintf("PCA: %s:%s:%s [%2.0f%% of total variance]",rownames(classes)[class],
                                               rownames(classes)[second],
                                               rownames(classes)[third],pvar)
            }
         } else {
         classvec = sprintf("%s:%s",classes[class,],classes[second,])      
         thetitle = sprintf("PCA: %s:%s [%2.0f%% of total variance]",rownames(classes)[class],
                                            rownames(classes)[second],pvar)
      }
   } else {
      classvec = classes[class,]	
      thetitle = sprintf("PCA: %s [%2.0f%% of total variance]",rownames(classes)[class],pvar)
  }
   
   numsamps = ncol(classes)
   if(numsamps <=10) {
      sval = 1.5
   } else if (numsamps <= 25) {
      sval = 1.25
   } else {
      sval = 1
   }   
   if(!missing(size)){
      if(!is.null(size)) {
         sval = size
      }
   }
   
   if(missing(vColors)) {
      vColors = pv.colsv
   }

   if(b3D) {
    if (length(find.package(package='rgl',quiet=T))>0) {
       library(rgl)
       plot3d(pc$loadings[,startComp:(startComp+2)],col=pv.colorv(classvec,vColors),type='s',size=sval,
              aspect=c(1,1,1),main=thetitle,...)
    } else {
       warning("Package rgl not installed")
       plot(pc$loadings[,startComp:(startComp+1)],col=pv.colorv(classvec,vColors),type='p',pch=19,cex=sval,
            xlab=sprintf('Principal Component #%d',startComp),ylab=sprintf('Principal Component #%d',startComp+1),
            main = thetitle,...)
    }
   } else {
   	  if(!missing(size)) {
   	     #sval = sval + .5
   	  }
      plot(pc$loadings[,startComp:(startComp+1)],col=pv.colorv(classvec,vColors),type='p',pch=19,cex=sval,
           xlab=sprintf('Principal Component #%d',startComp),ylab=sprintf('Principal Component #%d',startComp+1),
           main = thetitle,...)
   }      
   uclass = unique(classvec)
   res = matrix(vColors[1:length(uclass)],length(uclass),1)
   rownames(res) = uclass
   colnames(res) = "Legend"
   return(res)
   #return(cbind(uclass,vColors[1:length(uclass)])) 
}


## pv.plotHeatmap -- draw a heatmap using a subset of binding sites
PV_ONLYA = 3
PV_ONLYB = 4
PV_INALL = 5
PV_COR   = 6
PV_OLAP  = 7
PV_TOTAL = 0
pv.plotHeatmap = function(pv,numSites=1000,attributes=pv$attributes,mask,sites,contrast,
                          overlaps, olmask, olPlot=PV_COR,divVal,RowAttributes,ColAttributes,rowSideCols,colSideCols,
                          bTop=T,minval,maxval,bReorder=F,ColScheme="Greens",distMeth="pearson",...) {
#require(gplots)
#require(RColorBrewer)
#require(amap)
  
   pv = pv.check(pv)

   if(missing(mask)){
      mask = rep(T,ncol(pv$class))
    } else if(is.null(mask)) {
      	 mask = rep(T,ncol(pv$class))
    }
           
   ocm = NULL
   if(!missing(overlaps)) {
      cres  = overlaps
      if(!missing(olmask)) {
         cres = cres[olmask,]
   	  }
   	  if(is.null(nrow(cres))){
   	     cres = matrix(cres,1,length(cres))	
   	  }
   	  #cres = pv.overlapToLabels(pv,cres,attributes)
   	  if(missing(divVal)) {
   	     ocm = pv.occ2matrix(cres,olPlot)
   	  } else {
   	     ocm = pv.occ2matrix(cres,olPlot,divVal)
   	  }
   	  labels = pv.nums2labels(pv,rownames(ocm),attributes)
   	  rowlab = collab = labels
   	  rownames(ocm) = labels
   	  colnames(ocm) = labels
   	  domap = ocm

   } else {


      if(missing(sites)){
         sites = 1:nrow(pv$vectors)
         numSites = min(length(sites),numSites)
      } else {
   	     if(sum(sites)<=length(sites)){
            numSites = min(sum(sites),numSites)
            sites = which(sites)
         } else {
            numSites = length(sites)
         }
      }
   
      if(bTop==T) {
         sites = sites[1:numSites]
      } else {
         tsites = length(sites)
         sites = sites[(tsites-numSites+1):tsites]
      }
      colnames=NULL
      for(i in 1:ncol(pv$class)) {
         colnames = c(colnames,pv.namestrings(pv$class[attributes,i])$tstring)
      }
      domap = matrix(0.1,length(sites),sum(mask))
      for(i in 1:ncol(domap)) {
         domap[,i] = as.numeric(pv$vectors[sites,c(F,F,F,mask)][,i])
      }
      rowlab = ""
      collab = colnames[mask]
   }
   
   if(!missing(minval)) {
      domap[domap< minval]= minval
   }
   if(!missing(maxval)) {
      domap[domap>maxval]=maxval
   }
   cols = colorRampPalette(brewer.pal(9,ColScheme))(256)
   
   if(missing(rowSideCols)) {
      rowSideCols = pv.colsv	
   }
   rowatts = NULL
   rowcols = 0
   if(missing(RowAttributes)){
      if(!missing(contrast)) {
         rowatts = pv.attributematrix(pv,mask,contrast=contrast,PV_GROUP,rowSideCols)	
         rowcols = length(unique(as.vector(rowatts)))
      }
   } else if(!is.null(RowAttributes)) {
      rowatts = pv.attributematrix(pv,mask,contrast=contrast,RowAttributes,rowSideCols,bReverse=T)
      rowcols = length(unique(as.vector(rowatts)))
   } 

   if(missing(colSideCols)) {
      colSideCols = pv.colsv	
   }
   if(missing(ColAttributes)){
      colatts = pv.attributematrix(pv,mask,contrast=contrast,NULL,colSideCols,bAddGroup=is.null(ocm))	
      colcols = length(unique(as.vector(colatts)))
   } else if(!is.null(ColAttributes)) {
      colatts = pv.attributematrix(pv,mask,contrast=contrast,ColAttributes,colSideCols)
      colcols = length(unique(as.vector(colatts)))
   } else {
      colatts = NULL
      colcols = 0
   }
   
   if(is.null(ocm)){
   	  if(!missing(RowAttributes)) {
   	     warning("Row color bars not permitted for peak score heatmaps.",call.=F)	
   	  }
      heatmap.3(domap,labCol=collab,col=cols,trace="none",labRow=rowlab,
                distfun=function(x) Dist(x,method=distMeth),
                ColSideColors=colatts,NumColSideColors=colcols,
                ...)
   } else {

      res = heatmap.3(domap,labCol=collab,col=cols,trace="none",labRow=rowlab,
                      distfun=function(x) Dist(x,method=distMeth),symm=T,revC=T,Colv=T,
                      RowSideColors=rowatts,ColSideColors=colatts,NumRowSideColors=rowcols,NumColSideColors=colcols,
                      ...)         
      if(bReorder) {
      	 if(length(unique(rownames(ocm)))==nrow(ocm)) {
            ocm = pv.reorderM(ocm,res$rowDendrogram)
         } else {
            warning("Unable to re-order returned correlation matrix as labels are non-unique")	
         }
      }

      return(ocm)
   }

}

## pv.sort  - sort binding sites (e.g. for heatmap)
pv.sort = function(pv,fun=sd,mask,...) {
  
  if(missing(mask)){
     mask = rep(T,ncol(pv$class))
  }

  scores = apply(pv$vectors[,c(F,F,F,mask)],1,fun,...)
  ranked = order(scores,decreasing=T)
  
  pv$vectors   = pv$vectors[ranked,]
  
  return(pv)	
}
## pv.overlap -- generate overlapping/unique peaksets
pv.overlap = function(pv,mask,bFast=F,minVal=0) {
   
   if(!missing(mask)){
   	  if(!is.logical(mask)) {
         peaksets = mask
      }	else {
         peaksets = which(mask)
      }
      if (length(peaksets) <= 4) {
         A = peaksets[1]
         B = peaksets[2]
         if(length(peaksets) >= 3) {
            C = peaksets[3]
         } 
         if(length(peaksets) == 4) {
            D = peaksets[4]
         }
      } else {
         warning('Too many peaksets in mask.')
         return(NULL)
      }	
   } else {
   	  stop('Must specify mask for peaksets to overlap.',call.=F)
   }
   
   numcounts = sum(pv$class[PV_CALLER,peaksets] %in% "counts")
   if(numcounts == length(peaksets)) {
      if(is.null(pv$sites)) {
         stop("Called masks not present; re-run dba.count with bCalledMasks=TRUE",call.=F)	
      }
      maskA = pv$sites[[A]]
      maskB = pv$sites[[B]]   
      if(length(peaksets)>=3) maskC = pv$sites[[C]]
      if(length(peaksets)==4) maskD = pv$sites[[D]]         
   } else if (numcounts > 0) {
      stop("Mixed counted and uncounted peaksets, can't compute overlap",call.=F)	
   } else {
   	  #scores = pv$allvectors[,peaksets+3]
      pv = pv.vectors(pv,mask=peaksets,bAnalysis=F)
      #pv$allvectors[,4:ncol(pv$allvectors)] = scores
      pv$allvectors[,1] = pv$chrmap[pv$allvectors[,1]]
      if(length(peaksets) == 2) {
         A = 1
         B = 2
      } else if (length(peaksets) ==3) {
         A = 1
         B = 2
         C = 3       
      } else {
         A = 1
         B = 2
         C = 3
         D = 4      	
      }
      maskA = pv$allvectors[,3+A] > minVal   	
      maskB = pv$allvectors[,3+B] > minVal
      if(length(peaksets)>=3) maskC = pv$allvectors[,3+C] > minVal 
      if(length(peaksets)==4) maskD = pv$allvectors[,3+D] > minVal            
   }
   
   if(length(peaksets)<4){
      if(length(peaksets)<3) {
         res = pv.contrast2(pv$allvectors,A,B,minVal=minVal,v1=maskA,v2=maskB)        
      } else {
         res = pv.contrast3(pv$allvectors,A,B,C,minVal=minVal,v1=maskA,v2=maskB,v3=maskC)
      }
   } else {
      res = pv.contrast4(pv$allvectors,A,B,C,D,minVal=minVal,v1=maskA,v2=maskB,v3=maskC,v4=maskD) 
   }

   return(res)
}

## pv.plotVenn -- draw venn diagrams 
pv.plotVenn = function(ovrec,label1="A",label2="B",label3="C",label4="D",main="",sub="") {

   if(length(ovrec)==3) {
      pv.venn2(ovrec,label1,label2,main,sub)
   }

   if(length(ovrec)==7) {
      pv.venn3(ovrec,label1,label2,label3,main,sub)
   }
   
   if(length(ovrec)==15) {
      pv.venn4(ovrec,label1,label2,label3,label4,main,sub)
   }
}

## pv.occupancy-- generate co-occupancy stats from peaksets in a model
pv.occupancy = function(pv,mask,sites,byAttribute,Sort='inall',CorMethod="pearson",
                           labelAtts=pv$attributes,bPlot=F,minVal=0,bCorOnly=F,bNonZeroCors=F,chrmask) {
  
   pv = pv.check(pv)
   
   vecs=pv$allvectors
  
   if(missing(mask)) {
      mask = rep(T,ncol(pv$class))
   } else if(is.null(mask)) {
      mask = rep(T,ncol(pv$class))	
   } else {
      if(!is.logical(mask)) {
         tmp  = rep (F,length(pv$peaks))
         tmp[mask] = T
         mask = tmp
      }
   }
   if(missing(sites)) {
      sites = 1:nrow(vecs)
   }
  
   res=NULL     
   if(missing(byAttribute)){
   	  if(length(sites) < nrow(vecs)) {
   	     pv$allvectors = vecs[sites,]
   	     pv$vectors = pv$allvectors
   	  }
   	  if(!missing(chrmask)) {
   	     chrindex = match(chrmask,pv$chrmap)
   	     vecindex = pv$allvectors[,1] == chrmask
   	     pv$allvectors = pv$allvectors[vecindex,]
   	     vecindex = pv$vectors[,1] == chrmask
   	     pv$vectors = pv$vectors[vecindex,]   	     
   	  }
      res = pv.pairs(pv,mask=mask,CorMethod=CorMethod,bPlot=bPlot,minVal=minVal,bCorOnly=bCorOnly,bNonZeroCors=bNonZeroCors)
      if(!is.null(nrow(res))) {
         if(!missing(labelAtts)) {
            res = pv.overlapToLabels(pv,res,labelAtts)
         }
      }
   } else { ## by attribute
      vals = unique(pv$class[byAttribute,mask])
      for(i in 1:length(vals)) {
         comps = which(pv$class[byAttribute,] %in% vals[i])
         vmask = rep(F,length(mask))
         vmask[comps]=T
         if(sum(vmask)>1) {
            res = rbind(res,pv.occupancy(pv,mask=vmask,sites=sites,
                                        Sort=Sort,CorMethod=CorMethod,minVal=minVal,bCorOnly=bCorOnly))
         }
      }
   } 
   
   if(!is.null(nrow(res))) {
      if(Sort == 'cor') {
         res = res[pv.orderfacs(res[,6],decreasing=T),]
      } else if(Sort == 'percent') {
         res = res[pv.orderfacs(res[,7],decreasing=T),]
      } else {
         res = res[pv.orderfacs(res[,5],decreasing=T),]
      }   
   }
   
   return(res)	
}


## pv.plotBoxplot -- Boxplots
pv.plotBoxplot = function(DBA, contrast, method = DBA_EDGER, th=0.1, bUsePval=F, bNormalized=T, attribute=DBA_GROUP, 
                          bAll, bAllIncreased, bAllDecreased, bDB, bDBIncreased, bDBDecreased,
                          pvalMethod=wilcox.test,  bReversePos=FALSE, attribOrder, vColors, varwidth=T, notch=T, ...) {

   if(missing(bAll) && missing(bAllIncreased) && missing(bAllDecreased)) {
     bMissingAll = T	
   } else bMissingAll = F

   if(missing(bDB) && missing(bDBIncreased) && missing(bDBDecreased)) {
     bMissingDB = T	
   } else bMissingDB = F
   
   if(missing(contrast)) {
      if(bMissingAll) {
         bAll          = T
         bAllIncreased = F
         bAllDecreased = F
         bDB           = F
         bDBIncreased  = F
         bDBDecreased  = F         	
      }	
   } else {
      if(bMissingAll && bMissingDB) {
         bAll          = F
         bAllIncreased = F
         bAllDecreased = F
         bDB           = T
         bDBIncreased  = F
         bDBDecreased  = F   
      }	
   }

   if(missing(vColors)) {
      vColors = pv.colsv
      vColors = vColors[2:length(vColors)]
   }
   
   if(attribute==DBA_GROUP) {
      numPlots = 2
      cols = vColors[1:2]
      groups = list(DBA$class[PV_ID,DBA$contrasts[[contrast]]$group1],DBA$class[PV_ID,DBA$contrasts[[contrast]]$group2])   
      names = c(DBA$contrasts[[contrast]]$name1,DBA$contrasts[[contrast]]$name2) 
      if(!missing(attribOrder)) {
         if(attribOrder[1]==2 & attribOrder[2] ==1) {
            groups = list(DBA$class[PV_ID,DBA$contrasts[[contrast]]$group2],DBA$class[PV_ID,DBA$contrasts[[contrast]]$group1])   
            names = c(DBA$contrasts[[contrast]]$name2,DBA$contrasts[[contrast]]$name1) 
         }
      } 
   } else {
      samples = which(DBA$contrasts[[contrast]]$group1 | DBA$contrasts[[contrast]]$group2)
      classes = DBA$class[, samples]
      names = unique(classes[attribute,])
      numPlots = length(names)
      if(!missing(attribOrder)) {
      	  if(length(attribOrder < numPlots)) {
      	     neworder = 1:numPlots
      	     neworder[1:length(attribOrder)] = attribOrder
      	     attribOrder = neworder  
      	  }
         names  = names[attribOrder[1:numPlots]]
      }
      cols = vColors[1:numPlots]
      groups = NULL
      for(name in names) {
         groups = pv.listadd(groups,classes[PV_ID, classes[attribute,]==name])   
      }
   }
   
   subtitle = FALSE
   
   if(bAll | bAllIncreased | bAllDecreased) {
      report = pv.DBAreport(DBA, contrast=contrast, method = method,th=1,bNormalized=bNormalized,bCounts=T)
      if(bReversePos) {
         increase = report$Fold > 0   
         posgroup = DBA$contrasts[[contrast]]$name1
         neggroup = DBA$contrasts[[contrast]]$name2
      } else {
         increase = report$Fold < 0
         posgroup = DBA$contrasts[[contrast]]$name2
         neggroup = DBA$contrasts[[contrast]]$name1
      }
      if(bUsePval) {
         DB = report$p <= th
      } else {
         DB = report$FDR <= th
      }
   } else {
      report = NULL
   }
   
   toplot = NULL
   vcols    = NULL
   vnames = NULL

   if(bAll) {
      res       = lapply(groups,pv.box,report)
      for(i in 1:length(res)) {
         names(res)[i] = sprintf("%s",names[i])
      }
      toplot   = c(toplot,res)
      vnames = c(vnames,names)
      vcols    = c(vcols,cols)   
   }
   
   if(bAllIncreased) {
      res       = lapply(groups,pv.box,report[increase,])
      for(i in 1:length(res)) {
         names(res)[i] = sprintf("%s+",names[i])
      }      
      toplot   = c(toplot,res)
      vnames = c(vnames,rep("+",numPlots))
      vcols    = c(vcols,cols)
      subtitle = TRUE   
   }
   
   if(bAllDecreased) {
      res       = lapply(groups,pv.box,report[!increase,])
      for(i in 1:length(res)) {
         names(res)[i] = sprintf("%s-",names[i])
      }          
      toplot   = c(toplot,res)
      vnames = c(vnames,rep("-",numPlots))
      vcols    = c(vcols,cols)
      subtitle = TRUE      
   }
   
   
   if(is.null(report)) {
      report = pv.DBAreport(DBA, contrast=contrast, method=method,th=th,bUsePval=bUsePval,bNormalized=bNormalized,bCounts=T)
   } else {
      report = report[DB,]
   }
   
   if(nrow(report)==1) {
      stop('Need more than one DB site for boxplot')
   }
    
   if(bReversePos) {
      increase = report$Fold > 0   
      posgroup = DBA$contrasts[[contrast]]$name1
      neggroup = DBA$contrasts[[contrast]]$name2
   } else {
      increase = report$Fold < 0
      posgroup = DBA$contrasts[[contrast]]$name2
      neggroup = DBA$contrasts[[contrast]]$name1
   }
           
   if(bDB) {
      res       = lapply(groups,pv.box,report)
      for(i in 1:length(res)) {
         names(res)[i] = sprintf("%s.DB",names[i])
      }
      toplot   = c(toplot,res)
      vnames = c(vnames,names)
      vcols    = c(vcols,cols)   
   }   

   if(bDBIncreased) {
      res = lapply(groups,pv.box,report[increase,])
      for(i in 1:length(res)) {
         names(res)[i] = sprintf("%s.DB+",names[i])
      }
      toplot = c(toplot,res)
      vnames = c(vnames,rep("+",numPlots))
      vcols = c(vcols,cols) 
      subtitle = TRUE     
   }
   
   if(bDBDecreased) {
      res = lapply(groups,pv.box,report[!increase,])
       for(i in 1:length(res)) {
         names(res)[i] = sprintf("%s.DB-",names[i])
      }     
      toplot = c(toplot,res)
      vnames = c(vnames,rep("-",numPlots))
      vcols = c(vcols,cols) 
      subtitle = TRUE     
   }
   
   if(bNormalized==T) {
     ystr = "log2 normalized reads in binding sites"
   } else {
     ystr = "log2 reads in binding sites"   
   }
   
   if(subtitle == T) {
     subt = sprintf("+ indicates sites with increased affinity in %s\n- indicates sites with increased affinity in %s",
                    posgroup,neggroup)	
   } else {
     subt = ""	
   }
   boxplot(toplot,notch=notch, varwidth=varwidth,
           col=vcols,names=vnames,main="Binding affinity",        
           sub=subt,ylab=ystr)

   if(!is.null(pvalMethod)){
      pvals = matrix(1,length(toplot),length(toplot))
      rownames(pvals) = names(toplot)
      colnames(pvals) = names(toplot)
      for(i in 1:(length(toplot)-1)) {
         for(j in (i+1):length(toplot)) {
            pvals[i,j] = pvalMethod(toplot[[i]],toplot[[j]])$p.value
            pvals[j,i] = pvals[i,j]
         }
      }
   } else {
      pvals=NULL
   }     
   return(pvals)   
}

pv.box = function(ids,report) {
   idx = match(ids,colnames(report))
   res = log2(apply(report[,idx],1,mean))
   return(res)
}

