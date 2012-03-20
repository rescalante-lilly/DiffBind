pv.DBA = function(pv,method='edgeR',bSubControl=T,bFullLibrarySize=F,bTagwise=T,
                  minMembers=3,bParallel=F, block) {
  
   if(bParallel) {
      setParallel = TRUE
      bParallel = FALSE
   } else {
      setParallel = FALSE	
   }
   
   if(is.null(pv$contrasts)) {
   	  if(missing(block)) {
   	  	 pv = pv.contrast(pv,minMembers=minMembers)
   	  } else {
         pv = pv.contrast(pv,minMembers=minMembers,block=block)
      }
   }
   
   if(is.null(pv$contrasts)) {
      stop('Unable to perform analysis: no contrasts specified.')	
   }
   
   noreps = FALSE
   for(contrast in pv$contrasts) {
      if(sum(contrast$group1)<2) {
         noreps = TRUE
      }
      if(sum(contrast$group2)<2) {
         noreps = TRUE
      }
   }
   if(noreps) {
      warning("Some groups have no replicates. Results may be unreliable.",call.=FALSE)	
   }
     
  if(bParallel) {
  	 pv = dba.parallel(pv)
     jobs = NULL
     numjobs = 0
  }
  
  
  results = NULL
  if('edgeR' %in% method) {
     #require(edgeR)
     if(bParallel && (pv$config$parallelPackage > 0)) {
     	numjobs = numjobs + 1
        params = dba.parallel.params(pv$config,c('pv.allDEedgeR','pv.DEedgeR','pv.contrast','pv.listadd'))
        fdebug('submit job: pv.all')
     	jobs = pv.listadd(jobs,dba.parallel.addjob(pv$config,params,
     	                                          pv.allDEedgeR,pv,
     	                                          bFullLibrarySize=bFullLibrarySize,bParallel=T,
     	                                          bSubControl=bSubControl,bTagwise=bTagwise,bGLM=F))
     } else {
        results = pv.listadd(results, pv.allDEedgeR(pv,block=block,
                                                    bSubControl=bSubControl,bFullLibrarySize=bFullLibrarySize,
                                                    bParallel=setParallel,bTagwise=bTagwise,bGLM=F))
     }
  }
  if('edgeRGLM' %in% method) {
     #require(edgeR)
     if(bParallel && (pv$config$parallelPackage > 0)) {
     	numjobs = numjobs + 1
        params = dba.parallel.params(pv$config,c('pv.allDEedgeR','pv.DEedgeR','pv.contrast','pv.listadd'))
        fdebug('submit job: pv.all')
     	jobs = pv.listadd(jobs,dba.parallel.addjob(pv$config,params,
     	                                          pv.allDEedgeR,pv,
     	                                          bFullLibrarySize=bFullLibrarySize,bParallel=T,
     	                                          bSubControl=bSubControl,bTagwise=bTagwise,bGLM=T))
     } else {
        results = pv.listadd(results, pv.allDEedgeR(pv,block=block,
                                                    bSubControl=bSubControl,bFullLibrarySize=bFullLibrarySize,
                                                    bParallel=setParallel,bTagwise=bTagwise,bGLM=T))
     }
  }
  
   if('DESeq' %in% method) {
      if (length(find.package(package='DESeq',quiet=T))>0) {
        require(DESeq)
      } else {
        stop("Package DESeq not installed")
      }
      if(bParallel && (pv$config$parallelPackage > 0)) {
         numjobs = numjobs + 1
         params = dba.parallel.params(pv$config,c('pv.DESeq'))
         jobs = pv.listadd(jobs,dba.parallel.addjob(pv$config,params,pv.allDESeq,pv,
                                                    bSubControl=bSubControl,bFullLibrarySize=bFullLibrarySize,
                                                    bTagwise=bTagwise,bGLM=F,
                                                    bParallel=T))
      } else {
         results = pv.listadd(results,pv.allDESeq(pv,bSubControl=bSubControl,bFullLibrarySize=bFullLibrarySize,
                              bTagwise=bTagwise,bGLM=F,bParallel=setParallel))
      }
   }

   if('DESeqGLM' %in% method) {
      if (length(find.package(package='DESeq',quiet=T))>0) {
        require(DESeq)
      } else {
        stop("Package DESeq not installed")
      }
      if(bParallel && (pv$config$parallelPackage > 0)) {
         numjobs = numjobs + 1
         params = dba.parallel.params(pv$config,c('pv.DESeq'))
         jobs = pv.listadd(jobs,dba.parallel.addjob(pv$config,params,pv.allDESeq,pv,
                                                    bSubControl=bSubControl,bFullLibrarySize=bFullLibrarySize,
                                                    bTagwise=bTagwise,bGLM=T,
                                                    bParallel=T))
      } else {
         results = pv.listadd(results,pv.allDESeq(pv,bSubControl=bSubControl,bFullLibrarySize=bFullLibrarySize,
                              bTagwise=bTagwise,bGLM=T,bParallel=setParallel))
      }
   }
        
  
  if(bParallel && (pv$config$parallelPackage > 0)) {
     results = dba.parallel.wait4jobs(pv$config,jobs)
  }
  
  jnum = 1
  if( ('edgeR' %in% method) || ('edgeRGLM' %in% method) ) {
     edger = results[[jnum]]
     for(i in 1:length(edger)){
        pv$contrasts[[i]]$edgeR = edger[[i]]
     }
     jnum = jnum+1
  }
  if( ('DESeq' %in% method) || ('DESeqGLM' %in% method) ) {
     deseq = results[[jnum]]
     for(i in 1:length(deseq)){
        pv$contrasts[[i]]$DESeq = deseq[[i]]
     }
     jnum = jnum+1
  }

  fdebug(sprintf('Exit pv.DBA: %f',pv$contrasts[[1]]$edgeR$counts[7,1]))
  return(pv)
}


	
pv.DEinit = function(pv,mask1,mask2,group1=1,group2=2,method='edgeR',meanTH=0,
                     bSubControl=F,bFullLibrarySize=F,removeComps=0,bRawCounts=F,targets=NULL) {
  
  fdebug('enter pv.DEinit')
  
  edgeR = F
  DESeq = F
  if(method == 'edgeR') {
     #require('edgeR')
     edgeR = T   
  } else if (method == 'DESeq'|| method=='DESeqGLM') {
    if (length(find.package(package='DESeq',quiet=T))>0) {
       require(DESeq)
    } else {
       stop("Package DESeq not installed")
    }
     DESeq = T 
  } else {
    warning('Invalid method: ',method,call.=FALSE)
    return(NULL)
  }
	
  srcmask = pv.mask(pv,PV_CALLER,"source") | pv.mask(pv,PV_CALLER,"counts")
  
  fdebug(sprintf('pv.DEinit: %s %2.0f %s %2.0f srcmask %2.0f',
                 group1,sum(mask1),group2,sum(mask2),sum(srcmask)))
  g1 = which(mask1 & srcmask)
  g2 = which(mask2 & srcmask)
  numfirst = length(g1)
  
  s1 = pv.get_reads(pv,g1,bSubControl=bSubControl)
  s2 = pv.get_reads(pv,g2,bSubControl=bSubControl)
  
  if(meanTH > 0){
  	 mean1 = apply(s1,1,mean)
  	 mean2 = apply(s2,1,mean)
     idx1  = mean1>meanTH
     idx2  = mean2>meanTH
     keep  = idx1 | idx2
     s1    = s1[keep,]
     s2    = s2[keep,]
     mean1 = mean1[keep]
     mean2 = mean2[keep]    	
  } else {
     keep = rep(T,nrow(pv$peaks[[g1[1]]]))
  }
  counts = cbind(s1,s2)
  
  fdebug(sprintf('Generic counts: %f',counts[7,1]))

  rownames(counts) = as.character(1:nrow(counts))
  colnames(counts) = c(pv$class[PV_ID,mask1],pv$class[PV_ID,mask2])
  
  if(bRawCounts) {
     return(counts)	
  }
  
  groups = factor(c(rep(group1,length(g1)),rep(group2,length(g2))))
  if(bFullLibrarySize) {
     libsize = as.numeric(pv$class[PV_READS,c(g1,g2)])
  } else {
     libsize = NULL
  }
  if(!bRawCounts) {
     if(edgeR) {
        res = DGEList(counts,lib.size=libsize,
                      group=groups,genes=as.character(1:nrow(counts)))
        rownames(res$counts) = 1:nrow(res$counts)
        fdebug(sprintf('DGEList counts: %f',res$counts[7,1]))
     }
     if(DESeq) {
  	    colnames(counts) = NULL
  	    if(is.null(targets)) {
           res = newCountDataSet(counts,groups)
        } else {
           res = newCountDataSet(counts,targets)	
        }
        if(bFullLibrarySize) {
           sizeFactors(res) = libsize/min(libsize)
        }
     }
  }                
  return(res)
}

pv.DEedgeR = function(pv,group1,group2,label1="Group 1",label2="Group 2",blockList=NULL,
                      bSubControl=F,bFullLibrarySize=F,bTagwise=T,bGLM=T,bNormOnly=F) {

  fdebug('Enter pv.DEedgeR')
  
  #require(edgeR)

  res = pv.DEinit(pv,group1,group2,label1,label2,method='edgeR',
                  bSubControl=bSubControl,bFullLibrarySize=bFullLibrarySize)
  res = calcNormFactors(res,method="TMM")
  fdebug(sprintf('calcNormFactors: %f',res$counts[7,1]))
  
  if(bNormOnly) {
     return(res)	
  }
  
  if(is.null(blockList)) {
     fdebug('blockList is NULL')
  } else {
     fdebug('blockList is not NULL')
  }
  fdebug('pv.DEedgeR: check for blocking factor')
  if(is.null(blockList)) {
     fdebug('pv.DEedgeR: NO blocking factor')
     res = estimateCommonDisp(res)
     fdebug(sprintf('estimateCommonDisp: %f',res$counts[7,1]))
     if(bGLM){
        res$design = model.matrix(~res$samples$group)
        if(bTagwise) {
           res = estimateGLMCommonDisp(res,res$design)
           res = estimateGLMTagwiseDisp(res,res$design)
        } else {
           res = estimateGLMCommonDisp(res,res$design)
        }
        res$GLM = glmFit(res,res$design)
        res$LRT = glmLRT(res,res$GLM,2)
     } else {
        if(bTagwise){
  	       res = estimateTagwiseDisp(res,prior.n=50/(ncol(res$counts)-2),trend="none")
  	       #res = estimateTagwiseDisp(res,prior.n=getPriorN(res),trend="movingave")
  	       res$db     = exactTest(res,dispersion='tagwise')
        } else {
           res$db     = exactTest(res,dispersion='common')	
        }
     }
     fdebug(sprintf('Fit and test: %f',res$counts[7,1]))
     fdebug('pv.DEedgeR: estimateTagwiseDisp complete')
     
     fdebug(sprintf('pv.DEedgeR: exactTest complete:%s-%s',res$db$comparison[1],res$db$comparison[2]))
     #res$db$fdr = topTags(res$db,nrow(res$db$counts))
  } else {
  	 fdebug('pv.DEedgeR: BLOCKING FACTOR')
  	 
  	 targets = pv.blockFactors(pv,group1,group2,label1,label2,blockList)
  	 if(is.null(targets)){
  	    return(res)	
  	 }
  	 res$samples = data.frame(cbind(targets,res$samples[,2:3]))
  	 
     attr =  blockList[[1]]$attribute
     if(attr=='Replicate') {
        res$designmatrix = model.matrix(~ Replicate + group,data = targets)
     } else if(attr=='Tissue') {
        res$designmatrix = model.matrix(~ Tissue + group,data = targets)
     } else if(attr=='Factor') {
        res$designmatrix = model.matrix(~ Factor + group,data = targets)
     } else if(attr=='Condition') {
        res$designmatrix = model.matrix(~ Condition + group,data = targets)
     } else if(attr=='Caller') {
        res$designmatrix = model.matrix(~ Caller + group,data = targets)
     } else if(attr=='Treatment') {
        res$designmatrix = model.matrix(~ Treatment + group,data = targets)
     } else if(attr=='Block') {
        res$designmatrix = model.matrix(~ Block + group,data = targets)
     } else {
       warning('Unsupported blocking attribute: ',attr,call.=FALSE)
       return(NULL)	
     }
     message('edgeR multi-factor analysis.')
     res = calcNormFactors(res)
     res = estimateGLMCommonDisp(res,res$designmatrix)
     if(bTagwise) {
      res = estimateGLMTagwiseDisp(res,res$designmatrix)
     }
     res$GLM = glmFit(res,res$designmatrix)
     res$LRT = glmLRT(res,res$GLM,ncol(res$designmatrix))
     res$counts=NULL	 
     #res$fdr = topTags(res$LRT,nrow(res$counts))
  }
  
  res$bSubControl      = bSubControl
  res$bFullLibrarySize = bFullLibrarySize
  
  fdebug(sprintf('Exit pv.DEedgeR: %f',res$counts[7,1]))
  return(res)	

}

pv.DESeq = function(pv,group1,group2,label1="Group 1",label2="Group 2",
                    bSubControl=T,bFullLibrarySize=F,bTagwise=T,bGLM=T,blockList=NULL){
    if (length(find.package(package='DESeq',quiet=T))>0) {
       require(DESeq)
    } else {
       stop("Package DESeq not installed")
    }
   res = NULL
   
   targets=NULL
   if(!is.null(blockList)) {
      targets = pv.blockFactors(pv,group1,group2,label1,label2,blockList)
  	  if(is.null(targets)){
  	     return(res)	
  	  }
   } 
   res$DEdata = pv.DEinit(pv,group1,group2,label1,label2,method='DESeq',
	                      bSubControl=bSubControl,bFullLibrarySize=bFullLibrarySize,targets=targets)
   res$counts = counts(res$DEdata)
   if(!bFullLibrarySize) {
      res$DEdata = estimateSizeFactors(res$DEdata)
   }
   res$facs = sizeFactors(res$DEdata)
   if(sum(group1)+sum(group2)==2){
      res$DEdata = estimateDispersions(res$DEdata,fitType='local',method='blind',sharingMode='fit-only')
	} else {
	   if(bTagwise && !bGLM && is.null(blockList)) {
	      res$DEdata = estimateDispersions(res$DEdata,fitType='local',method='per-condition')
	   } else {
	      res$DEdata = estimateDispersions(res$DEdata,fitType='local',method='pooled')
	      if(bTagwise && is.null(blockList)) {
	         warning('Unable to use tagwise dispersion estimates with GLM',call.=FALSE)	
	      }
	   }
   }
   if(bGLM && is.null(blockList)){
      res$fullFit      = fitNbinomGLMs(res$DEdata,count ~ condition)
	  res$reducedFit   = fitNbinomGLMs(res$DEdata,count ~ 1)
	  res$de           = nbinomGLMTest(res$fullFit,res$reducedFit)
	  res$de           = cbind(1:length(res$de),res$de,p.adjust(res$de,method="BH"))
	  colnames(res$de) = c('id','pval','padj')
	  res$de           = data.frame(res$de)
	} else if (!is.null(blockList)) {
	   message('DESeq multi-factor analysis')
       attr =  blockList[[1]]$attribute
       if(attr=='Replicate') {
          res$fullFit    = fitNbinomGLMs(res$DEdata,count ~ group + Replicate)
          res$reducedFit = fitNbinomGLMs(res$DEdata,count ~ Replicate)
       } else if(attr=='Tissue') {
          res$fullFit    = fitNbinomGLMs(res$DEdata,count ~ group + Tissue)
          res$reducedFit = fitNbinomGLMs(res$DEdata,count ~ Tissue)
       } else if(attr=='Factor') {
          res$fullFit    = fitNbinomGLMs(res$DEdata,count ~ group + Factor)
          res$reducedFit = fitNbinomGLMs(res$DEdata,count ~ Factor)
       } else if(attr=='Condition') {
          res$fullFit    = fitNbinomGLMs(res$DEdata,count ~ group + Condition)
          res$reducedFit = fitNbinomGLMs(res$DEdata,count ~ Condition)
       } else if(attr=='Caller') {
          res$fullFit    = fitNbinomGLMs(res$DEdata,count ~ group + Caller)
          res$reducedFit = fitNbinomGLMs(res$DEdata,count ~ Caller)
       } else if(attr=='Treatment') {
          res$fullFit    = fitNbinomGLMs(res$DEdata,count ~ group + Treatment)
          res$reducedFit = fitNbinomGLMs(res$DEdata,count ~ Treatment)
       } else if(attr=='Block') {
          res$fullFit    = fitNbinomGLMs(res$DEdata,count ~ group + Block)
          res$reducedFit = fitNbinomGLMs(res$DEdata,count ~ Block)
       } else {
          warning('Unsupported blocking attribute: ',attr,call.=FALSE)
          return(NULL)	
       }
      res$de           = nbinomGLMTest(res$fullFit,res$reducedFit)
	  res$de           = cbind(1:length(res$de),res$de,p.adjust(res$de,method="BH"))
	  colnames(res$de) = c('id','pval','padj')
	  res$de           = data.frame(res$de)
    } else {
	   res$de = nbinomTest(res$DEdata,label1,label2)[,c(1,7:8)]
	}
	  
   res$de = res$de[order(res$de$padj),]

   res$bSubControl      = bSubControl
   res$bFullLibrarySize = bFullLibrarySize

   return(res)

}

pv.DEedgeR_parallel = function(contrast,pv,blockList,bSubControl,bFullLibrarySize,bTagwise,bGLM) {
   crec = pv$contrasts[[contrast]]
   if(!is.null(blockList)) {
      blockList = crec$blocklist
   }
   res = pv.DEedgeR(pv,crec$group1,crec$group2,crec$name1,crec$name2,blockList=blockList,
                    bSubControl=bSubControl,bFullLibrarySize=bFullLibrarySize,bTagwise=bTagwise,bGLM=bGLM)
   
   fdebug(sprintf('Exit pv.DEedgeR_parallel: %f',res$counts[7,1]))
   return(res)
}

                      	
pv.allDEedgeR = function(pv,block,bFullLibrarySize=F,bParallel=F,bSubControl=F,bTagwise=T,bGLM=F) {
 
fdebug('ENTER pv.allDEedgeR')
#require(edgeR)
   
   if(is.null(pv$contrasts)) {
   	  if(missing(block)) {
   	  	 pv$contrasts = pv.contrast(pv)
   	  } else {
         pv$contrasts = pv.contrast(pv,block=block)
      }
   }
   
   if(bParallel) {
      pv = dba.parallel(pv)
   	  jobs = NULL
   	  blocks = NULL
   }
   
   fdebug('pv.allDEedgeR: for each contrast')
   reslist = NULL
  
   if(bParallel && (pv$config$parallelPackage > 0)) {
      params = dba.parallel.params(pv$config,c('pv.DEedgeR_parallel','pv.DEedgeR','pv.DEinit','calcNormFactors',
                                              'estimateCommonDisp','estimateTagwiseDisp',
                                              'estimateGLMCommonDisp','estimateGLMTagwiseDisp',
                                              'exactTest','topTags',
   	                                          'glmFit','glmLRT'))
      reslist = dba.parallel.lapply(pv$config,params,1:length(pv$contrasts),pv.DEedgeR_parallel,pv,
                                    NULL,bSubControl,bFullLibrarySize,bTagwise,bGLM=bGLM)
                                    
      fdebug(sprintf('Return from parallel call to pv.DEedgeR_parallel: %f',reslist[[1]]$counts[7,1]))
      
      blist = NULL
      for(i in 1:length(pv$contrasts)) {
         if(!is.null(pv$contrasts[[i]]$blocklist)) {
            blist = c(blist,i)
         }
      }
      if(length(blist > 0)) {
         bres =  dba.parallel.lapply(pv$config,params,blist,pv.DEedgeR_parallel,pv,
                                    TRUE,bSubControl,bFullLibrarySize,bTagwise)
         for(i in 1:length(blist)) {
            reslist[[blist[i]]]$block = bres[[i]]
         }    
      }     
   } else { #SERIAL
      for(i in 1:length(pv$contrast)) { 	
         res = pv.DEedgeR(pv,pv$contrasts[[i]]$group1,pv$contrasts[[i]]$group2,
                          pv$contrasts[[i]]$name1,pv$contrasts[[i]]$name2,
                          bSubControl=bSubControl,bFullLibrarySize=bFullLibrarySize,bTagwise=bTagwise,bGLM=bGLM)
      
         if(!is.null(pv$contrasts[[i]]$blocklist)) {
            res$block = pv.DEedgeR(pv,pv$contrasts[[i]]$group1,pv$contrasts[[i]]$group2,
                                   pv$contrasts[[i]]$name1,pv$contrasts[[i]]$name2,
                                   pv$contrasts[[i]]$blocklist,
                                   bSubControl=bSubControl,bFullLibrarySize=bFullLibrarySize,bTagwise=bTagwise)   
         }
         reslist = pv.listadd(reslist,res)   
      }
   }
   
   fdebug(sprintf('Exit pv.allDEedgeR: %f',reslist[[1]]$counts[7,1]))
   return(reslist)
}


pv.DESeq_parallel = function(contrast,pv,blockList,bSubControl,bFullLibrarySize,bTagwise=T,bGLM=F) {
   crec = pv$contrasts[[contrast]]
   if(!is.null(blockList)) {
      blockList = crec$blocklist
   }
   res = pv.DESeq(pv,crec$group1,crec$group2,crec$name1,crec$name2,
                    bSubControl=bSubControl,bFullLibrarySize=bFullLibrarySize,bTagwise=bTagwise,bGLM=bGLM,blockList=blockList)
   return(res)
}

pv.allDESeq = function(pv,block,bSubControl=F,bFullLibrarySize=F,bTagwise=T,bGLM=F,bParallel=F) {

  	if (length(find.package(package='DESeq',quiet=T))>0) {
       require(DESeq)
    } else {
       stop("Package DESeq not installed")
    }

   if(is.null(pv$contrasts)) {
   	  if(missing(block)) {
   	  	 pv = pv.contrast(pv)
   	  } else {
         pv = pv.contrast(pv,block=block)
      }
   }
   
   if(bParallel) {
      pv = dba.parallel(pv)
      jobs = NULL
      blocks = NULL
   }
   
   fdebug('pv.allDESeq: for each contrast')
   reslist = NULL
   
   if(bParallel && (pv$config$parallelPackage > 0)) {   
   	  params = dba.parallel.params(pv$config,c('pv.DESeq_parallel','pv.DESeq'))
   	  
      reslist  = dba.parallel.lapply(pv$config,params,1:length(pv$contrasts),pv.DESeq_parallel,pv,NULL, 
                                     bSubControl,bFullLibrarySize,bTagwise=bTagwise,bGLM=bGLM)
      blist = NULL
      for(i in 1:length(pv$contrasts)) {
         if(!is.null(pv$contrasts[[i]]$blocklist)) {
            blist = c(blist,i)
         }
      }
      if(length(blist > 0)) {
         bres =  dba.parallel.lapply(pv$config,params,1:length(pv$contrasts),pv.DESeq_parallel,pv,TRUE, 
                                     bSubControl,bFullLibrarySize,bTagwise=bTagwise,bGLM=bGLM)
         for(i in 1:length(blist)) {
            reslist[[blist[i]]]$block = bres[[i]]
         }    
      }     
   } else { # SERIAL
      for(i in 1:length(pv$contrast)) { 	
         res = pv.DESeq(pv,pv$contrasts[[i]]$group1,pv$contrasts[[i]]$group2,
                          pv$contrasts[[i]]$name1,pv$contrasts[[i]]$name2,
                          bSubControl=bSubControl,bFullLibrarySize=bFullLibrarySize,bTagwise=bTagwise,bGLM=bGLM)
      
         if(!is.null(pv$contrasts[[i]]$blocklist)) {
            res$block = pv.DESeq(pv,pv$contrasts[[i]]$group1,pv$contrasts[[i]]$group2,
                                   pv$contrasts[[i]]$name1,pv$contrasts[[i]]$name2,
                                   bSubControl=bSubControl,bFullLibrarySize=bFullLibrarySize,bTagwise=bTagwise,
                                   blockList=pv$contrasts[[i]]$blocklist)   
         }
         reslist = pv.listadd(reslist,res)   
      }
   }  
   
   return(reslist)
}


pv.blockFactors = function(pv,group1,group2,label1,label2,blockList) {
	 samples = group1 | group2
  	 targets = NULL
  	 for(block in blockList) {
  	    samps = block$samples & samples &  
  	            (pv.mask(pv,PV_CALLER,"source") | pv.mask(pv,PV_CALLER,"counts"))
  	    IDs   = pv$class[PV_ID,samps]
  	    groups = NULL
  	    wsamps = which(samps)
  	    for(sw in wsamps) {
  	       if(group1[sw]) {
  	          groups = c(groups,label1)
  	       } else {
  	          groups = c(groups,label2)
  	       }   
  	    }
  	    block = rep(block$label, sum(samps))
  	    block = cbind(groups,block)
  	    rownames(block) = IDs
  	    targets = rbind(targets,block)
  	 }
  	 
  	 snames = c(colnames(pv$class[,group1]),colnames(pv$class[,group2]))
  	 if(length(unique(snames))!=sum(group1|group2)){
  	    warning('Error: all samples must have unique IDs for blocking analysis',call.=FALSE)
  	    return(NULL)	
  	 }
  	 tnames = rownames(targets)
  	 #if(length(snames)!=length(tnames)){
  	 #   warning('Error: all samples must be matched for blocking analysis')
  	 #   return(res)	
  	 #}  	 
  	 
  	 newt = targets
  	 for(i in 1:nrow(targets)) {
  	    idx = match(tnames[i],snames)
  	    newt[idx,] = targets[i,]	
  	 }
  	 targets = newt
  	 rownames(targets) = snames
  	 
  	 colnames(targets) = c("group",blockList[[1]]$attribute)

     targets = data.frame(targets)
     targets[[1]] = factor(targets[[1]])
     targets[[2]] = factor(targets[[2]])
	
	 return(targets)
	
}



pv.DBAreport = function(pv,contrast=1,method='edgeR',th=.1,bUsePval=F,bCalled=F,
                        bCounts=F,bCalledDetail=F,
                        file,initString='reports/DBA',bNormalized=T,ext="csv",minFold=0,bSupressWarning=F) {
   
   if(contrast > length(pv$contrasts)) {
      stop('Specified contrast number is greater than number of contrasts')
      return(NULL)
   }
   con = pv$contrasts[[contrast]]
   if(method=='edgeR' || method=='edgeRGLM'){
      if(is.null(con$edgeR) || class(con$edgeR)=="try-error") {
         stop('edgeR analysis has not been run for this contrast')
         return(NULL)
      }
      if(is.null(con$edgeR$counts)) {
      	 counts = pv.DEinit(pv,con$group1,con$group2,con$label1,con$label2,method='edgeR',
                            bSubControl=con$edgeR$bSubControl,bFullLibrarySize=con$edgeR$bFullLibrarySize,bRawCounts=TRUE)
      } else {
         counts = con$edgeR$counts	
      }
      if(!is.null(con$edgeR$LRT)) {
         siteCol = 1
         pvCol   = 5
         fdrCol  = 6
         data = topTags(con$edgeR$LRT,nrow(counts))$table   	  	
   	  } else {
         siteCol = 1
         pvCol   = 4
         fdrCol  = 5
         data = topTags(con$edgeR$db,nrow(counts))$table
      }

      if(bNormalized){
      	 sizes = con$edgeR$samples$lib.size * con$edgeR$samples$norm.factors
      	 counts = t(t(counts)/sizes)
      	 counts = counts * con$edgeR$common.lib.size
      } 
   } else if (method=='DESeq' || method=='DESeqGLM' || method=='DESeqBlock') {
      if (length(find.package(package='DESeq',quiet=T))>0) {
         require(DESeq)
      } else {
         stop("Package DESeq not installed")
      }
   	  if(is.null(con$DESeq) || class(con$DESeq)=="try-error") {
         stop('DESeq analysis has not been run for this contrast') 
         return(NULL) 
      }
      siteCol = 1
      pvCol   = 2
      fdrCol =  3
      counts = pv.DEinit(pv,con$group1,con$group2,con$label1,con$label2,method='DESeq',
	                     bSubControl=con$DESeq$bSubControl,bFullLibrarySize=con$DESeq$bFullLibrarySize,bRawCounts=T)
      if(method=='DESeqBlock') {
         data = con$DESeq$block$de
         if(bNormalized){
      	    counts = t(t(counts)/con$DESeq$block$facs)
         }      	
      } else {
         data = con$DESeq$de
         if(bNormalized){
      	    counts = t(t(counts)/con$DESeq$facs)
         }
      }   
   } else if(method=='edgeRlm'){
   	  if(is.null(con$edgeR$counts)) {
         counts = pv.DEinit(pv,con$group1,con$group2,con$label1,con$label2,method='edgeR',
                            bSubControl=con$edgeR$block$bSubControl,bFullLibrarySize=con$edgeR$block$bFullLibrarySize,
                            bRawCounts=TRUE)
      } else {
         counts = con$edgeR$counts	
      }

      siteCol = 1
      pvCol   = 5  
      fdrCol  = 6
      data = topTags(con$edgeR$block$LRT,nrow(counts))$table

      if(bNormalized){
      	 sizes = con$edgeR$samples$lib.size * con$edgeR$samples$norm.factors
      	 counts = t(t(counts)/sizes)
      	 counts = counts * con$edgeR$common.lib.size
      } 
   } else {
      stop('Unknown DE method: ',method)
      return(NULL)
   }
   if(bUsePval) {
     thCol = pvCol	
   } else {
     thCol = fdrCol	
   }

   x = which(is.na(data[,pvCol]))
   if(length(x)>0){
      data[x,pvCol] = 1
   }
   x = which(is.na(data[,fdrCol]))
   if(length(x)>0){
      data[x,fdrCol] = 1
   }
     
   keep =  data[,thCol]<=th
   sites = as.numeric(data[keep,siteCol])
   if(sum(keep)==0) {
   	  if(!bSupressWarning) {
         warning('No sites above threshold',call.=FALSE)
      }
      return(NULL)
   } else if(sum(keep)==1) {
      cnames = colnames(counts)
      counts = matrix(counts[sites,],nrow=1,ncol=ncol(counts))
      colnames(counts) = cnames
   } else {
      counts = counts[sites,]
   }
   
   if(length(sites)==1) {
      conc = log2(mean(counts))
      if(sum(con$group1)>1) {
         con1 = log2(mean(counts[,1:sum(con$group1)]))
      } else {
         con1 = log2(counts[,1])
      }
      if(sum(con$group2)>1) {
         con2 = log2(mean(counts[,(sum(con$group1)+1):ncol(counts)]))
      } else {
         con2 = log2(counts[,sum(con$group1)+1])
      }        
   } else {
      conc = log2(apply(counts,1,mean))
      if(sum(con$group1)>1) {
         con1 = log2(apply(counts[,1:sum(con$group1)],1,mean))
      } else {
         con1 = log2(counts[,1])
      }
      if(sum(con$group2)>1) {
         con2 = log2(apply(counts[,(sum(con$group1)+1):ncol(counts)],1,mean))
      } else {
         con2 = log2(counts[,sum(con$group1)+1])
      }
   }
   fold = con1 - con2
   
   #sites = apply(pv$peaks[[which(con$group1)[1]]][data[keep,siteCol],1:3],1,pv.dositename)
   
   data = cbind(pv.getsites(pv,sites),conc,con1,con2,fold,data[keep,c(pvCol,fdrCol)])
   
   conc1 = sprintf('Conc_%s',con$name1)
   conc2 = sprintf('Conc_%s',con$name2)
   
   colnames(data) = c('Chr','Start','End','Conc',conc1,conc2,'Fold','p-value','FDR')
   
   if(bCalled & !is.null(pv$sites)) {
   	  called1 = rep(F,length(pv$sites[[1]]))
      toadd = which(con$group1)
      for(i in toadd) {
         called1 = called1 + pv$sites[[i]]
      }
      Called1 = called1[sites]
      called2 = rep(F,length(pv$sites[[1]]))
      toadd = which(con$group2)
      for(i in toadd) {
         called2 = called2 + pv$sites[[i]]
      }
      Called2 = called2[sites]
      data = cbind(data,Called1,Called2)
   }
   
   
   if(bCounts) {
      colnames(counts) = c(pv$class[PV_ID,con$group1],pv$class[PV_ID,con$group2])
      if(length(sites)>1){
         data = cbind(data,counts)
      } else {
         dnames = colnames(data)
         cnames = colnames(counts)
         data = cbind(data,matrix(counts,1,ncol(counts)))
         colnames(data) = c(dnames,cnames)
      }
   }

   if(bCalledDetail & !is.null(pv$sites)) {
      toadd = c(which(con$group1),which(con$group2))
      newd = NULL
      for(i in toadd) {
         newd = cbind(newd,pv$sites[[i]][sites])
      }
      newd[newd==T] = '+'
      newd[newd==F] = '-'
      colnames(newd) = c(pv$class[PV_ID,con$group1],pv$class[PV_ID,con$group2])
      data = cbind(data,newd)
   }   

   if(minFold>0) {
      data = data[abs(data$Fold)>=minFold,]
   }
   
   data = data[order(data$'p-value'),]
   
   if(!missing(file)) {
   	 if(is.null(file)) {
        file=sprintf("%s_%s_vs_%s_%s.%s",initString,con$name1,con$name2,method,ext)
     } else {
        file=sprintf("%s_%s.%s",initString,file,ext)
     }
     write.csv(data,row.names=F,file=file)
   }
      
   return(data)
   
}

pv.getsites = function(pv,sites){
   if(is.logical(sites)) {
      sites = which(sites)
   }
   idx   = match(sites,rownames(pv$allvectors))
   sites = pv$allvectors[idx,1:3]
   sites[,1] = pv$chrmap[sites[,1]]
   return(sites)
}

pv.DBAplotMA = function(pv,contrast,method='edgeR',bMA=T,bXY=F,th=0.1,bUsePval=F,facname="",bNormalized=T,
                        cex=.15,...) {

   if(missing(contrast)){
      contrast=1:length(pv$contrasts)
   } else {
     if(contrast > length(pv$contrasts)) {
       stop('Specified contrast number is greater than number of contrasts')
       return(NULL)
     }
   }
   
   plotfun = plot
#   if (bSmooth) {
#      plotfun = smoothScatter
#   }

   
   numSites = nrow(pv$vectors)
   
   for(con in 1:length(contrast)) {
   	  conrec = pv$contrasts[[contrast[con]]]
   	  for(meth in method) {
         res = pv.DBAreport(pv,contrast=contrast[con],method=meth,bUsePval=T,th=100,bNormalized=bNormalized)
         if(!is.null(res)) {
            if(bUsePval) {
              idx = res$"p-value" <= th
              tstr = "p"
            } else {
              idx = res$FDR <= th
              tstr = "FDR"
           }
           if(bMA){
              xmin  = floor(min(res$Conc))
              xmax  = ceiling(max(res$Conc))
              ymin  = floor(min(res$Fold))
              ymax  = ceiling(max(res$Fold))
              plotfun(res$Conc[!idx],res$Fold[!idx],pch=20,cex=cex,
                      xaxp=c(xmin,xmax,xmax-xmin),xlim=c(xmin,xmax),
                      xlab='log concentration',
                      yaxp=c(ymin,ymax,(ymax-ymin)),ylim=c(ymin,ymax),
                      ylab=sprintf('log fold change: %s - %s',conrec$name1,conrec$name2),
                      main=sprintf('%s Binding Affinity: %s vs. %s (%s %s < %1.3f)',
                                   facname, conrec$name1,conrec$name2,sum(idx),tstr,th),...)
              points(res$Conc[idx],res$Fold[idx],pch=20,cex=cex,col=2)
              abline(h=0,col='dodgerblue')
           }
           if(bXY){
              xmin  = floor(min(res[,5]))
              xmax  = ceiling(max(res[,5]))
              ymin  = floor(min(res[,6]))
              ymax  = ceiling(max(res[,6]))
              xymin = min(xmin,ymin)
              xymin = max(xymin,0)
              xymax = max(xmax,ymax)
              plotfun(res[!idx,6],res[!idx,5],pch=20,cex=cex,col=1,
                      xaxp=c(xymin,xymax,xymax-xymin),xlim=c(xymin,xymax),
                      xlab=sprintf('log concentration :%s',conrec$name2),
                      yaxp=c(xymin,xymax,(xymax-xymin)),ylim=c(xymin,xymax),
                      ylab=sprintf('log concentration :%s',conrec$name1),
                      main=sprintf('%s Binding Affinity: %s vs. %s (%s %s < %1.3f)',
                                   facname, conrec$name1,conrec$name2,sum(idx),tstr,th),...)
              points(res[idx,6],res[idx,5],pch=20,cex=cex,col=2)
              abline(0,1,col='dodgerblue')
            }
         }
      }      	
   }	
}

pv.normTMM = function(pv,bMinus=TRUE,bFullLib=FALSE){
   
   if(length(pv$peaks)<2) {
      warning('Unable to TMM normalize -- not enough peaksets',call.=FALSE)
      return(pv)	
   }
   
   g1     = rep(F,length(pv$peaks))
   g1[1]  = T
   res    = pv.DEedgeR(pv,g1,!g1,"1","2",bSubControl=bMinus,bFullLibrarySize=bFullLib,bNormOnly=T)
   #res    = estimateCommonDisp(res)
   counts = res$counts
   sizes  = res$samples$lib.size * res$samples$norm.factors
   counts = t(t(counts)/sizes)
   counts = counts * mean(res$samples$lib.size)

   return(counts)
                       
}

pv.stripDBA = function(conrec) {
   conrec$edgeR = pv.stripEdgeR(conrec$edgeR)
   conrec$DESeq = pv.stripDESeq(conrec$DESeq)
   return(conrec)
}

pv.stripEdgeR = function(erec) {
   if(!is.null(erec)) {
   	  erec$counts     = NULL
      erec$pseudo.alt = NULL
      erec$conc       = NULL
      #erec$genes      = NULL
      erec$all.zeros  = NULL
      erec$tagwise.dispersion = NULL

      if(!is.null(erec$GLM)) {
         erec$GLM = pv.stripEdgeRGLM(erec$GLM)
      }
      if(!is.null(erec$LRT)) {
         erec$LRT = pv.stripEdgeRLRT(erec$LRT)
      }
      if(!is.null(erec$block)) {
         erec$block = pv.stripEdgeR(erec$block)
      }
      if(!is.null(erec$db)) {
         erec$db = pv.stripEdgeR(erec$db)
      }      
   }
   return(erec)
}

pv.stripEdgeRGLM = function(grec) {
   if(!is.null(grec)) {
   	  grec$counts       = NULL
      grec$fitted.values= NULL
      grec$offset       = NULL
      grec$coefficients = NULL
      grec$deviance     = NULL
      grec$df.residual  = NULL
      grec$abundance    = NULL
      #grec$genes        = NULL
   }
   return(grec)
}

pv.stripEdgeRLRT = function(lrec) {
   if(!is.null(lrec)) {
      lrec = pv.stripEdgeR(lrec)
      lrec = pv.stripEdgeRGLM(lrec)
      lrec$GLM = pv.stripEdgeRGLM(lrec$GLM)
   
      lrec$dispersion.used   = NULL
   }
   return(lrec)
}      
   
pv.stripDESeq = function(drec) {
	if(!is.null(drec)) {
	   drec$counts     = NULL
	   drec$DEdata     = NULL
	   drec$fullFit    = NULL
	   drec$reducedFit = NULL
	   
	   if(!is.null(drec$block)) {
	      drec$block = pv.stripDESeq(drec$block)	
	   }	
	}
	return(drec)
}



