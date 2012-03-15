pv.DBA = function(pv,method='edgeR',bSubControl=T,bFullLibrarySize=F,bTagwise=T,
                  minMembers=3,bParallel=F, block) {
  
   if(bParallel) {
      setParallel = TRUE
      bParallel   = FALSE
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
      warning("Some groups have no replicates. Results may be unreliable.")	
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
     	                                          bSubControl=bSubControl,bTagwise=bTagwise))
     } else {
        results = pv.listadd(results, pv.allDEedgeR(pv,block=block,
                                                    bSubControl=bSubControl,bFullLibrarySize=bFullLibrarySize,
                                                    bParallel=setParallel,bTagwise=bTagwise))
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
                                                    bTagwise=bTagwise,
                                                    bParallel=T))
      } else {
         results = pv.listadd(results,pv.allDESeq(pv,bSubControl=bSubControl,bFullLibrarySize=bFullLibrarySize,
                              bTagwise=bTagwise,bParallel=setParallel))
      }
   }     
  if('t-test' %in% method) {
     if(bParallel && (pv$config$parallelPackage > 0)) {
       numjobs = numjobs + 1
       params = dba.parallel.params(pv$config,c('pv.Ttests'))
       jobs =  pv.listadd(jobs,dba.parallel.addjob(pv$config,params,pv.Ttests,pv,bParallel=T))
     } else {
       results = pv.listadd(results,pv.Ttests(pv,bParallel=setParallel))
     }
  }
  
  if(bParallel && (pv$config$parallelPackage > 0)) {
     results = dba.parallel.wait4jobs(pv$config,jobs)
  }
  
  jnum = 1
  if('edgeR' %in% method) {
     edger = results[[jnum]]
     for(i in 1:length(edger)){
        pv$contrasts[[i]]$edgeR = edger[[i]]
     }
     jnum = jnum+1
  }
  if('DESeq' %in% method) {
     deseq = results[[jnum]]
     for(i in 1:length(deseq)){
        pv$contrasts[[i]]$DESeq = deseq[[i]]
     }
     jnum = jnum+1
  }
  if('t-test' %in% method) {
     ttestr = results[[jnum]]
     for(i in 1:length(ttestr)){
        pv$contrasts[[i]]$Ttest = ttestr[[i]]
     }
     jnun = jnum+1
  }  
  return(pv)
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

pv.DBAplotHeatmap = function(pv,contrast,method='edgeR',th=.1, bUsePval=F,
                             maxSites=1000,PCA=PV_GROUP,minval=0,maxval,removeComps=0,bLog=T,
                             bPCAonly=F,b3D=T,vColors,...) {

   if(missing(contrast)){
      contrast=1:length(pv$contrasts)
   }
   
   if(missing(minval)) {
      minval = -10000000
   }
   
   for(conrec in contrast) {
      for(meth in method) {
         if((meth == 'edgeR') || (meth == 'edgeRtw')) {
            res = pv.DEedgeRplots(pv,contrast=conrec,bMA=F,bHeatmap=!bPCAonly,bUsePval=bUsePval,th=th,
                             maxSites=maxSites,PCA=PCA,minval=minval,maxval=maxval,
                             removeComps=removeComps,bLog=bLog,b3D=b3D,vColors=vColors,...)
         }
         if(meth == 'DESeq') {
            res = pv.DESeqplots(pv,contrast=conrec,bMA=F,bHeatmap=!bPCAonly,bUsePval=bUsePval,th=th,
                             maxSites=maxSites,PCA=PCA,minval=minval,maxval=maxval,
                             removeComps=removeComps,bLog=bLog,b3D=b3D,vColors=vColors,...)
         }         
         if(meth == 't-test') {
            warning('t-test Heatmap plots not implemented.')
         }   
      }	
   }
   return(res)	
}

	
pv.DEinit = function(pv,mask1,mask2,group1=1,group2=2,method='edgeR',meanTH=0,
                     bSubControl=F,bFullLibrarySize=F,removeComps=0,bRawCounts=F) {
  
  fdebug('enter pv.DEinit')
  
  edgeR = F
  DESeq = F
  if(method == 'edgeR') {
     #require('edgeR')
     edgeR = T   
  } else if (method == 'DESeq') {
    if (length(find.package(package='DESeq',quiet=T))>0) {
       require(DESeq)
    } else {
       stop("Package DESeq not installed")
    }
     DESeq = T 
  } else {
    warning('Invalid method: ',method)
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

  rownames(counts) = as.character(1:nrow(counts))
  colnames(counts) = c(pv$class[PV_ID,mask1],pv$class[PV_ID,mask2])
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
     }
     if(DESeq) {
  	    colnames(counts) = NULL
        res = newCountDataSet(counts,groups)
        if(bFullLibrarySize) {
           sizeFactors(res) = libsize/min(libsize)
        }
     }
  }                
  return(res)
}

pv.DEedgeR = function(pv,group1,group2,label1="Group 1",label2="Group 2",blockList=NULL,
                      bSubControl=F,bFullLibrarySize=F,bTagwise=T) {

  fdebug('Enter pv.DEedgeR')
  
  #require(edgeR)

  res = pv.DEinit(pv,group1,group2,label1,label2,method='edgeR',
                  bSubControl=bSubControl,bFullLibrarySize=bFullLibrarySize)
  res = calcNormFactors(res,method="TMM")
  if(is.null(blockList)) {
     fdebug('blockList is NULL')
  } else {
     fdebug('blockList is not NULL')
  }
  fdebug('pv.DEedgeR: check for blocking factor')
  if(is.null(blockList)) {
     fdebug('pv.DEedgeR: NO blocking factor')
     res = estimateCommonDisp(res)
     if(bTagwise){
  	    res = estimateTagwiseDisp(res,prior.n=50/(ncol(res$counts)-2),trend="none")
  	    res$db     = exactTest(res,dispersion='tagwise')
     } else {
        res$db     = exactTest(res,dispersion='common')	
     }
     fdebug('pv.DEedgeR: estimateTagwiseDisp complete')
     
     fdebug(sprintf('pv.DEedgeR: exactTest complete:%s-%s',res$db$comparison[1],res$db$comparison[2]))
     #res$db$fdr = topTags(res$db,nrow(res$db$counts))
  } else {
  	 fdebug('pv.DEedgeR: BLOCKING FACTOR')
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
  	 
  	 snames = rownames(res$samples)
  	 if(length(unique(snames))!=nrow(res$samples)){
  	    warning('Error: all samples must have unique IDs for blocking analysis')
  	    return(res)	
  	 }
  	 tnames = rownames(targets)
  	 if(length(snames)!=length(tnames)){
  	    warning('Error: all samples must be matched for blocking analysis')
  	    return(res)	
  	 }  	 
  	 
  	 newt = targets
  	 for(i in 1:nrow(targets)) {
  	    idx = match(tnames[i],snames)
  	    newt[idx,] = targets[i,]	
  	 }
  	 targets = newt
  	 rownames(targets) = snames
  	 
  	 colnames(targets) = c("group",blockList[[1]]$attribute)
     res$samples = data.frame(cbind(targets,res$samples[,2:3]))
     targets = data.frame(targets)
     targets[[1]] = factor(targets[[1]])
     targets[[2]] = factor(targets[[2]])
     attr =  blockList[[1]]$attribute
     if(attr=='Replicate') {
        res$designmatrix = model.matrix(~ Replicate + group,data = targets)
     } else {
       warning('Unsupported blocking attribute: ',attr)
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
  
  return(res)	

}


pv.DESeq = function(pv,group1,group2,label1="Group 1",label2="Group 2",
                    bSubControl=T,bFullLibrarySize=F,bTagwise=T){
    if (length(find.package(package='DESeq',quiet=T))>0) {
       require(DESeq)
    } else {
       stop("Package DESeq not installed")
    }
   res = NULL
   res$DEdata = pv.DEinit(pv,group1,group2,label1,label2,method='DESeq',
                          bSubControl=bSubControl,bFullLibrarySize=bFullLibrarySize)
   res$counts = counts(res$DEdata)
   if(!bFullLibrarySize) {
      res$DEdata = estimateSizeFactors(res$DEdata)
   }
   res$facs = sizeFactors(res$DEdata)
   #res$counts = t(t(res$counts)/res$facs) #scale(res$counts,FALSE,1/facs)
   if(sum(group1)+sum(group2)==2){
      res$DEdata = estimateDispersions(res$DEdata,fitType='local',method='blind',sharingMode='fit-only')
   } else {
   	  if(bTagwise) {
   	  	 res$DEdata = estimateDispersions(res$DEdata,fitType='local',method='per-condition')
   	  } else {
         res$DEdata = estimateDispersions(res$DEdata,fitType='local',method='pooled')
      }
   }
   res$de = nbinomTest(res$DEdata,label1,label2)
  
   res$de = res$de[order(res$de$padj),]

   return(res)

}

pv.DEedgeR_parallel = function(contrast,pv,blockList,bSubControl,bFullLibrarySize,bTagwise) {
   crec = pv$contrasts[[contrast]]
   if(!is.null(blockList)) {
      blockList = crec$blocklist
   }
   res = pv.DEedgeR(pv,crec$group1,crec$group2,crec$name1,crec$name2,blockList=blockList,
                    bSubControl=bSubControl,bFullLibrarySize=bFullLibrarySize,bTagwise=bTagwise)
   return(res)
}

                      	
pv.allDEedgeR = function(pv,block,bFullLibrarySize=F,bParallel=F,bSubControl=F,bTagwise=T) {
 
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
                                    NULL,bSubControl,bFullLibrarySize,bTagwise)
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
                          bSubControl=bSubControl,bFullLibrarySize=bFullLibrarySize,bTagwise=bTagwise)
      
         if(!is.null(pv$contrasts[[i]]$blocklist)) {
            res$block = pv.DEedgeR(pv,pv$contrasts[[i]]$group1,pv$contrasts[[i]]$group2,
                                   pv$contrasts[[i]]$name1,pv$contrasts[[i]]$name2,
                                   pv$contrasts[[i]]$blocklist,
                                   bSubControl=bSubControl,bFullLibrarySize=bFullLibrarySize,bTagwise=bTagwise)   
         }
         reslist = pv.listadd(reslist,res)   
      }
   }
   
   fdebug('Exit pv.allDEedgeR')
   return(reslist)
}


pv.DESeq_parallel = function(contrast,pv,bSubControl,bFullLibrarySize,bTagwise=bTagwise) {
   crec = pv$contrasts[[contrast]]
   res = pv.DESeq(pv,crec$group1,crec$group2,crec$name1,crec$name2,
                    bSubControl=bSubControl,bFullLibrarySize=bFullLibrarySize,bTagwise=bTagwise)
   return(res)
}

pv.allDESeq = function(pv,block,bSubControl=F,bFullLibrarySize=F,bParallel=F,bTagwise=T) {

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
   	  params = dba.parallel.params(pv$config,c('pv.DESeq_parallel','pv.DESeq'))
      reslist  = dba.parallel.lapply(pv$config,params,1:length(pv$contrasts),pv.DESeq_parallel,pv, 
                                             bSubControl,bFullLibrarySize,bTagwise=bTagwise)

   } else {
   	  reslist = lapply(1:length(pv$contrasts),pv.DESeq_parallel,pv,
                         bSubControl=bSubControl,bFullLibrarySize=bFullLibrarySize,bTagwise=bTagwise)
   }  
   
   return(reslist)
}


pv.Ttests = function(pv,block,scores,bParallel=F) {
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
   
   reslist = NULL
   
   for(i in 1:length(pv$contrasts)) {
   	  contr = pv$contrasts[[i]]
      numG1 = sum(contr$group1)
      if(missing(scores)) {
         sc1 = NULL
         sc2 = NULL
      } else if(scores=='edgeR') {
         sc1 = contr$edgeR$pseudo.alt[,1:numG1]
         sc2 = contr$edgeR$pseudo.alt[,(numG1+1):ncol(contr$edgeR$pseudo.alt)]
      } else if(scores=='DESeq') {
      	 sc1 = contr$DESeq$counts[,1:numG1]
      	 sc2 = contr$DESeq$counts[,(numG1+1):ncol(contr$DESeq$counts)] 
      }
      
      if(bParallel && (pv$config$parallelPackage > 0)) {
   	     params = dba.parallel.params(pv$config,c('pv.test_groups'))
         jobs = pv.listadd(jobs,dba.parallel.addjob(pv$config,params,pv.test_groups,pv,
                                                    contr$group1,contr$group2,
                                                    scores1=sc1,scores2=sc2))
     } else {
         reslist = pv.listadd(reslist,pv.test_groups(pv,contr$group1,contr$group2,
                                                     scores1=sc1,scores2=sc2)) 
      }   
   }
   
   
   if(bParallel && (pv$config$parallelPackage > 0)) {
      reslist = dba.parallel.wait4jobs(pv$config,jobs)
      #if(length(pv$contrasts)==1) {
      #   reslist = list(reslist)
      #}
   }

   return(reslist)
}


pv.test_groups = function(pv,mask1,mask2,meanTH=0,scores1=NULL,scores2=NULL) {
	
  srcmask = pv.mask(pv,PV_CALLER,"source") | pv.mask(pv,PV_CALLER,"counts")
  g1 = which(mask1 & srcmask)
  g2 = which(mask2 & srcmask)
  numfirst = length(g1)
  
  if(is.null(scores1)) {
     s1 = pv.get_scores(pv,g1)
  } else s1 = scores1
  mean1 = apply(s1,1,mean)
  
  if(is.null(scores2)) {
     s2 = pv.get_scores(pv,g2)
  } else s2 = scores2
  mean2 = apply(s2,1,mean)
  
  if(meanTH > 0){
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
  
  scores = cbind(s1,s2)
  pvals = apply(scores,1,pv_dotest,numfirst)
  pvals[is.na(pvals)]=1
  
  adjust = p.adjust(pvals,"BH")
  
  res = cbind(pv$peaks[[g1[1]]][keep,1:3],pvals,adjust,mean1,mean2,scores)
  colnames(res) = c("chr","start","end","p-value","adjusted","mean 1","mean 2",pv$class[PV_ID,g1],pv$class[PV_ID,g2])
  
  return(res)

}


pv_dotest = function(tvec,numfirst) {
  ret = tryCatch(t.test(tvec[1:numfirst],tvec[(numfirst+1:length(tvec))])$p.value,error=function(e) 1)
  return(ret)
}

pv.compareDE = function(pv,contrast=1,numTop=100) {
 con = pv$contrasts[[contrast]]
 ttest   = order(con$Ttest[,4],decreasing=F)
 deseq   = as.character(con$DESeq$de$id)
 edger   = as.character(con$edgeR$de$fdr$table$genes)
 edgertw = as.character(con$edgeR$db$fdr$table$genes)
 
 if(is.null(con$edgeR$block)) {
    res = matrix(0,4,4)
    rownames(res) = c('TTest','edgeR','edgeRtw','DESeq')
 } else {
 	edgerlm = as.character(con$edgeR$block$fdr$table$genes)
    res = matrix(0,5,5)
    rownames(res) = c('TTest','edgeR','edgeRtw','edgeRlm','DESeq')
 }
 colnames(res) = rownames(res)
 
 res[1,1] = mean(match(ttest[1:numTop],ttest))
 res[1,2] = mean(match(ttest[1:numTop],edger))
 res[1,3] = mean(match(ttest[1:numTop],edgertw))
 res[1,4] = mean(match(ttest[1:numTop],deseq))
 res[2,1] = mean(match(edger[1:numTop],ttest))
 res[2,2] = mean(match(edger[1:numTop],edger))
 res[2,3] = mean(match(edger[1:numTop],edgertw)) 
 res[2,4] = mean(match(edger[1:numTop],deseq))
 res[3,1] = mean(match(edgertw[1:numTop],ttest))
 res[3,2] = mean(match(edgertw[1:numTop],edger))
 res[3,3] = mean(match(edgertw[1:numTop],edgertw))
 res[3,4] = mean(match(edgertw[1:numTop],deseq))
 res[4,1] = mean(match(deseq[1:numTop],ttest))
 res[4,2] = mean(match(deseq[1:numTop],edger)) 
 res[4,3] = mean(match(deseq[1:numTop],edgertw))
 res[4,4] = mean(match(deseq[1:numTop],deseq))
   
 if(nrow(res)==5){
    res[,5]	= res[,4]
    res[5,] = res[4,]
    res[1,4] = mean(match(ttest[1:numTop],edgerlm))
    res[2,4] = mean(match(edger[1:numTop],edgerlm))
    res[3,4] = mean(match(edgertw[1:numTop],edgerlm))
    res[4,4] = mean(match(edgertw[1:numTop],edgertw))
    res[5,4] = mean(match(deseq[1:numTop],edgerlm))
    res[4,1] = mean(match(edgerlm[1:numTop],ttest))
    res[4,2] = mean(match(edgerlm[1:numTop],edger))
    res[4,3] = mean(match(edgerlm[1:numTop],edgertw))
    res[4,5] = mean(match(edgerlm[1:numTop],deseq))  
 }
 
 return(res-.5)
 
}

pv.compareDE.wilcox = function(pv,contrast=1,method1='edgeRtw',method2='DESeq') {
   con = pv$contrasts[[contrast]]
   c1 = pv.getDE(con,method1)
   c2 = pv.getDE(con,method2)
   
   o1 = order(c1$IDs,decreasing=F)
   o2 = order(c2$IDs,decreasing=F)
   
   res = wilcox.test(c1$fdrs[o1],c2$fdrs[o2],paired=T)
   
   return(res)	

}

pv.getDE = function(con,method) {
   res = NULL
   if(method=='edgeR') {
      res$IDs  = con$edgeR$de$fdr$table[,1]
      res$fdrs = con$edgeR$de$fdr$table[,5]
   }
   if(method=='edgeRtw') {
      res$IDs  = con$edgeR$db$fdr$table[,1]
      res$fdrs = con$edgeR$db$fdr$table[,5]
   }
   if(method=='edgeRblock') {
      res$IDs  = con$edgeR$block$fdr$table[,1]
      res$fdrs = con$edgeR$block$fdr$table[,5]
   }
    if(method=='DESeq') {
      res$IDs  = con$DESeq$de$id
      res$fdrs = con$DESeq$de$padj
   }
   if(method=='t-test') {
      res$IDs  = as.character(1:nrow(con$Ttest))
      res$fdrs = con$TTest[,5]
   }   
   res$IDs = as.character(res$IDs)
   return(res)      	
}




pv.DEedgeRplots = function(pv,contrast=1,bBlock=F,maxSites=1000,bMA=T,bHeatmap=T,bLog=T,
                           PCA=NULL,bUsePval=F,th=.1,minval=0,maxval,
                           addLines=2,removeComps=0,b3D=T,vColors,...) {
   #require('edgeR')
   con = pv$contrasts[[contrast]]
   
   if(is.null(con$edgeR) || class(con$edgeR)=="try-error") {
      stop("edgeR analysis has not been run for this contrast")
      return	
   }

   if(is.null(maxSites)) {
      maxSites=nrow(pv$vectors)
   }
   if(bBlock) {
      erec = con$edgeR$block
      etable = topTags(erec,nrow(erec$counts))$table
      method='edgeRlm'
   } else {
      erec = con$edgeR$db
      etable = topTags(erec,nrow(erec$counts))$table
      method='edgeRtw'
   } 
   if(bUsePval) {
     #numsig = sum(erec$fdr$table[,4] <= th)
     numsig = sum(etable[,4] <= th)
     nstat = "Pval"
   } else {
     #numsig = sum(erec$fdr$table[,5] <= th)
     numsig = sum(etable[,5] <= th)
     nstat = "FDR"
   }
   numuse = min(numsig,maxSites)
   if(numuse <1) {
      warning('No binding sites in plot.')
      return(NULL)
   }
   desites = etable[1:numsig,1]
   if(bMA){
      #plotSmear(con$edgeR,de.tags=desites,pch=20,cex=.1,main=sprintf("%s: %s vs %s (%d %s <= %1.3f)",
      #                                                 method,con$name1,con$name2,numsig,nstat,th),...)
      # if(addLines>0) {
      #    abline(h=c(-addLines,addLines),col='dodgerblue')
      # }
   }
   
   pv1 = NULL
   pv1$attributes = pv$attributes
   
   
   data = con$edgeR$pseudo.alt
   data[which(data<1)]=1
   if(bLog) {
      data = log2(data)
   }
   if(removeComps!=0) {
     data = pv.removeComp(data,removeComps)
   }        
   pv1$vectors = cbind(data[,1],data[,1],data[,1],data)
   pv1$class   = cbind(pv$class[,con$group1],pv$class[,con$group2])
   fctab = etable[1:numuse,]
   pos   = fctab[fctab[,3]>0,]
   poso  = order(pos[,3],decreasing=F)
   neg   = fctab[fctab[,3]<=0,]
   nego  = order(neg[,3],decreasing=F)
   fctab = rbind(pos[poso,],neg[nego,])
   pv1$vectors = pv1$vectors[fctab[,1],]
   #pv1$vectors = pv1$vectors[1:numuse,]

   pv1$vectors[which(pv1$vectors<minval)]=minval
   if(!missing(maxval)) {
      pv1$vectors[which(pv1$vectors>maxval)]=maxval
   }
   #sums = apply(pv1$vectors[1:numuse,],2,sum)
   sums = apply(pv1$vectors,2,sum)
   pv1$vectors[1,sums==0]=.1   
   
   if(bHeatmap){
      pv.plotHeatmap(pv1,numSites=numuse,main=sprintf("%s vs %s (%d %s <= %1.3f)", 
                                                       con$name1,con$name2,numsig,nstat,th),...) 
   }
   
   if(!is.null(PCA)) { 
      pv1 = pv.analysis(pv1)
      if(PCA[1] == PV_GROUP) {
         pv1$class[PV_ID,] = c(rep(con$name1,sum(con$group1)),rep(con$name2,sum(con$group2)))
         PCA = PV_ID
      }
      if(missing(vColors)) {
         vColors = pv.colsv
      }
      pv.plotPCA(pv1,PCA,b3D=b3D,vColors=vColors,...)	
   }                                                       
}

pv.DESeqMA = function(DErec,th=0.05,bUsePval=F,bFlip=F,...) {
   if(bFlip) {
      flip = -1
   } else {
      flip = 1 
   }
   if(bUsePval) {
      cols =ifelse(DErec$de$pval <= th, "red","black")
   } else {
      cols =ifelse(DErec$de$padj <= th, "red","black")
   }
   plot(DErec$de$baseMean, flip * DErec$de$log2FoldChange,log="x",pch=20,cex=.1,
        col = cols,...) 
   abline(h=c(-2,2),col='dodgerblue')
           
}

pv.DESeqplots = function(pv,contrast=1,maxSites=1000,bMA=T,bHeatmap=T,PCA=NULL,
                         bUsePval=F,th=.1,minval=0,maxval,removeComps=0,bLog=T,b3D=T,vColors,...) {
    if (length(find.package(package='DESeq',quiet=T))>0) {
       require(DESeq)
    } else {
       stop("Package DESeq not installed")
    }
   con = pv$contrasts[[contrast]]

   if(is.null(con$DESeq) || class(con$DESeq)=="try-error") {
      stop("DESeq analysis has not been run for this contrast")
      return(NULL) 
   }
   
   if(bUsePval) {
      numsig = sum(con$DESeq$de$pval <= th)
      nstat = "Pval"
   } else {
      numsig = sum(con$DESeq$de$padj <= th)
      nstat = "FDR"
   }
   numuse = min(numsig,maxSites)
   if(numuse <1) {
      warning('No binding sites in plot.')
      return(NULL)
   }
   
   desites = con$DESeq$de$id[1:numsig]
   if(bMA){
   	  pv.DESeqMA(con$DESeq,th,bUsePval=bUsePval,main=sprintf("DESeq: %s vs %s (%d %s <= %1.3f)",
                                                       con$name1,con$name2,numsig,nstat,th))
   }
   
   pv1 = NULL
   pv1$attributes = pv$attributes
   
   #con$DESeq$counts[which(con$DESeq$counts<1)]=1
   data = con$DESeq$counts
   data[which(data<1)]=1
   if(bLog){
      data = log2(data)
   }  
   if(removeComps!=0) {
     data = pv.removeComp(data,removeComps)
   }  
   pv1$vectors = cbind(data[,1],data[,1],data[,1],data)
   
   pv1$class   = cbind(pv$class[,con$group1],pv$class[,con$group2])
   fctab = con$DESeq$de[1:numuse,]
   pos   = fctab[fctab[,6]>0,]
   poso  = order(pos[,6],decreasing=F)
   neg   = fctab[fctab[,6]<=0,]
   nego  = order(neg[,6],decreasing=F)
   fctab = rbind(pos[poso,],neg[nego,])
   pv1$vectors = pv1$vectors[fctab[,1],]
   #pv1$vectors = pv1$vectors[1:numuse,]
      
   pv1$vectors[which(pv1$vectors<minval)]=minval
   if(!missing(maxval)) {
      pv1$vectors[which(pv1$vectors>maxval)]=maxval
   }
   #sums = apply(pv1$vectors[1:numuse,],2,sum)
   sums = apply(pv1$vectors,2,sum)
   pv1$vectors[1,sums==0]=.1   
      
   if(bHeatmap){
      pv.plotHeatmap(pv1,numSites=numuse,main=sprintf("%s vs %s (%d FDR <= %1.3f)", 
                                                                 con$name1,con$name2,numsig,th),...) 
   }
   
   if(!is.null(PCA)) { 
      pv1 = pv.analysis(pv1)
      if(PCA == PV_GROUP) {
         pv1$class[PV_ID,] = c(rep(con$name1,sum(con$group1)),rep(con$name2,sum(con$group2)))
         PCA = PV_ID
      }     
      if(missing(vColors)) {
         vColors = pv.colsv
      } 
      pv.plotPCA(pv1,PCA,b3D=b3D,vColors=vColors,...)	
   }                                                       
}

pv.DBAreport = function(pv,contrast=1,method='edgeR',th=.1,bUsePval=F,bCalled=F,
                        bCounts=F,bCalledDetail=F,
                        file,initString='reports/DBA',bNormalized=T,ext="csv",minFold=0) {
   
   if(contrast > length(pv$contrasts)) {
      stop('Specified contrast number is greater than number of contrasts')
      return(NULL)
   }
   con = pv$contrasts[[contrast]]
   if(method=='edgeR'){
      siteCol = 1
      pvCol   = 4
      fdrCol  = 5
      #data = con$edgeR$db$fdr$table
      if(is.null(con$edgeR) || class(con$edgeR)=="try-error") {
         stop('edgeR analysis has not been run for this contrast')
         return(NULL)
      }
      data = topTags(con$edgeR$db,nrow(con$edgeR$db$counts))$table
      counts = con$edgeR$counts
      if(bNormalized){
      	 sizes = con$edgeR$samples$lib.size * con$edgeR$samples$norm.factors
      	 counts = t(t(counts)/sizes)
      	 counts = counts * con$edgeR$common.lib.size
      } 
   } else if (method=='DESeq') {
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
      pvCol   = 7
      fdrCol =  8
      data = con$DESeq$de
      #negs = data[,4] < data[,3]
      #data[negs,5] = -(data[negs,3]/data[negs,4])
      #data[,5] = -data[,5]
      counts = con$DESeq$counts

      if(bNormalized){
      	 counts = t(t(counts)/con$DESeq$facs)
      	 #facs = sizeFactors(con$DESeq$DEdata)
         #for(i in length(facs))  {
         #   counts[,i] = counts[,i] / facs[i]
         #}
      }   
   } else if(method=='edgeRlm'){
      siteCol = 1
      pvCol   = 5  
      fdrCol  = 6
      #data = con$edgeR$block$fdr$table
      data = topTags(con$edgeR$block$LRT,nrow(con$edgeR$counts))$table
      counts = con$edgeR$counts
      if(bNormalized){
      	 sizes = con$edgeR$samples$lib.size * con$edgeR$samples$norm.factors
      	 counts = t(t(counts)/sizes)
      	 counts = counts * con$edgeR$common.lib.size
      } 
   } else {
      warning('Unknown DE method: ',method)
      return(NULL)
   }
   if(bUsePval) {
     thCol = pvCol	
   } else {
     thCol = fdrCol	
   }

   keep =  data[,thCol]<=th
   sites = as.numeric(data[keep,siteCol])
   if(sum(keep)==0) {
      stop('No sites above threshold')
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
