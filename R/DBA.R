#############################################
## DBA.R -- Differential Binding Analysis  ##
## 20 April 2011                           ##
## Rory Stark                              ##
## Cancer Research UK                      ##
#############################################

## dba	            Construct a dba object
## dba.peakset	    Add a peakset to a dba object
	
## dba.overlap	    Compute binding site overlaps
## dba.count  	    Count reads in binding sites
	
## dba.contrast	    Establish contrast(s) for analysis
## dba.analyze  	Execute affinity analysis
## dba.report	    Generate report for a contrast analysis
	
## dba.plotClust	Cluster dendrogram plo
## dba.plotHeatmap	Heatmap plot
## dba.plotPCA	    Principal Components plot
## dba.plotBox	    Boxplots
## dba.plotMA	    MA/scatter plot
## dba.plotVenn	    2, 3, or 4-way Venn diagram plot
	
## dba.show	        List dba metadata
## dba.mask	        Mask samples or sites 
	
## dba.save	        Save dba object
## dba.load	        Load dba object

### NOTE: DBA is a front-end to a package formerly known as pv, with most DBA
### functions being simple pass-thoughs to pv functions

#########################################################
## dba -- construct DBA object, e.g. from sample sheet ##
#########################################################
DBA_VERSION1  = 1
DBA_VERSION2  = 1
DBA_VERSION3  = 4

DBA_GROUP     = PV_GROUP
DBA_ID        = PV_ID 
DBA_TISSUE    = PV_TISSUE 
DBA_FACTOR    = PV_FACTOR
DBA_CONDITION = PV_CONDITION
DBA_TREATMENT = PV_TREATMENT
DBA_CONSENSUS = PV_CONSENSUS
DBA_CALLER    = PV_CALLER
DBA_CONTROL   = PV_CONTROL
DBA_READS     = PV_READS
DBA_REPLICATE = PV_REPLICATE



DBA_EDGER_CLASSIC = 'edgeR'
DBA_EDGER_BLOCK   = 'edgeRlm'
DBA_EDGER_GLM     = 'edgeRGLM'
DBA_EDGER         = DBA_EDGER_GLM

DBA_DESEQ_CLASSIC = 'DESeq'
DBA_DESEQ_BLOCK   = 'DESeqBlock'
DBA_DESEQ_GLM     = 'DESeqGLM'
DBA_DESEQ         = DBA_DESEQ_GLM

DBA_DATA_FRAME      = 0
DBA_DATA_RANGEDDATA = 1
DBA_DATA_GRANGES    = 2
DBA_DATA_DEFAULT    = DBA_DATA_GRANGES

dba = function(DBA,mask, minOverlap=2,
               sampleSheet="dba_samples.csv", 
               config=data.frame(RunParallel=TRUE, reportInit="DBA"),
               caller="raw", skipLines=0, bAddCallerConsensus=FALSE, 
               bRemoveM=TRUE, bRemoveRandom=TRUE, 
               bCorPlot=FALSE, attributes) 
{
   if(!missing(DBA)){
      DBA = pv.check(DBA)	
   }
   
   res = pv.model(DBA, mask=mask, minOverlap=minOverlap, samplesheet=sampleSheet, config=config, 
                    caller=caller, skipLines=skipLines,bAddCallerConsensus=bAddCallerConsensus, 
                    bRemoveM=bRemoveM, bRemoveRandom=bRemoveRandom,
                    bKeepAll=TRUE, bAnalysis=TRUE, 
                    attributes=attributes)
                    
   res$contrasts=NULL
   
   if(is.null(res$config$DataType)) {
      res$config$DataType=DBA_DATA_DEFAULT
   }

   if(is.null(res$config$parallelPackage)){
      res$config$parallelPackage=DBA_PARALLEL_MULTICORE
   }
   if(is.null(res$config$RunParallel)){
      res$config$RunParallel=TRUE
   }
   if(is.null(res$config$AnalysisMethod)){
      res$config$AnalysisMethod=DBA_EDGER
   }

   if(missing(DBA)){
      DBA=NULL
   } 
   if(is.null(DBA$config$reportInit)){
      res$config$reportInit="DBA"
   } else {
      res$config$reportInit=DBA$config$reportInit
   }
   if(class(res)!="DBA") {
      class(res) = "DBA"
   }
   
   if(bCorPlot) {
      dba.plotHeatmap(res)
   }
    
   return(res)                 
}

###############################################
## dba.peakset -- add a peakset to the model ##
###############################################

dba.peakset = function(DBA=NULL, peaks, sampID, tissue, factor, condition, treatment, replicate,
                       control, peak.caller, reads=0, consensus=FALSE, bamReads, bamControl,
                       normCol=4, bRemoveM=TRUE, bRemoveRandom=TRUE,
                       minOverlap=2, bMerge=TRUE,
                       bRetrieve=FALSE, writeFile, numCols=4,
                       DataType=DBA$config$DataType)
{
   res = NULL
   
   if(!missing(peaks)){
   	 if(class(peaks) != "DBA") {
        peaks = pv.DataType2Peaks(peaks)
     } else {
        peaks = pv.check(peaks)	
     }
   }
   
   if(bRetrieve==TRUE || !missing(writeFile)) {
   
      if(missing(writeFile)) {
         writeFile = NULL
      }
     
      if(missing(peaks)) {
      	 if(!is.null(DBA)) {
           DBA = pv.check(DBA)
         } else {
            stop("DBA object is NULL; can't retrieve peakset.")	
         }	
      } else {
         if(is.vector(peaks)) {
      	    if(is.null(DBA)) {
              stop("DBA object is NULL; can't retrieve peakset.")	
            }	
      	    if(length(peaks) > 1) {
      	 	   if(minOverlap > length(peaks)) {
      	 	      stop('minOverlap is greater than number of specified peaks.')	
      	 	   } else {
                  DBA = dba(DBA,mask=peaks,minOverlap=minOverlap,bCorPlot=FALSE)
                  peaks = NULL
               }
            }      
         }	
      }
      
      res = pv.writePeakset(DBA, fname=writeFile, peaks=peaks, numCols=numCols)     
      
      if(DataType!=DBA_DATA_FRAME) {
         res = pv.peaks2DataType(res,DataType)
      }    
   
   } else {
   	  if(!is.null(DBA)) {
   	     DBA = pv.check(DBA)
   	  }
   	  if(!missing(peaks)) {
   	     if(class(peaks)=="DBA") {
   	        res = pv.peakset_all(DBA, addpv=peaks, minOverlap=minOverlap)
   	     }
   	   }
   	   if(is.null(res)) {
   
         res = pv.peakset(DBA, peaks=peaks, 
                          sampID=sampID, tissue=tissue, factor=factor,condition=condition,treatment=treatment,
                          replicate=replicate,control=control,
                          peak.caller=peak.caller,reads=reads, consensus=consensus, 
                          readBam=bamReads, controlBam=bamControl,
                          bNormCol=normCol, bRemoveM=bRemoveM, bRemoveRandom=bRemoveRandom,
                          minOverlap=minOverlap)
      }
      
      if(class(res)!="DBA") {
         class(res) = "DBA"
      }
    
      if(is.null(res$config$DataType)) {
         res$config$DataType=DBA_DATA_DEFAULT
      }
      if(is.null(res$config$parallelPackage)){
         res$config$parallelPackage=DBA_PARALLEL_MULTICORE
      }
      if(is.null(res$config$RunParallel)){
         res$config$RunParallel=TRUE
      }
      if(is.null(DBA$config$reportInit)){
         res$config$reportInit="DBA"
      } else {
         res$config$reportInit=DBA$config$reportInit
      }
      if(is.null(res$config$AnalysisMethod)){
         res$config$AnalysisMethod=DBA_EDGER
      }
            
      if(bMerge) {
         res = pv.check(res)
      }
   }   
                       
   return(res)                       

}                      

##################################################
## dba.overlap -- compute binding site overlaps ####################################################

DBA_OLAP_PEAKS = 1 # Return list of peaksets (common/unique peaks) 
DBA_OLAP_ALL   = 2 # Return overlap report with statstics for peakset pairs
DBA_OLAP_RATE  = 3 # Return vector of number of peaks in overlap for all values of minOverlap

DBA_OLAP  = PV_OLAP
DBA_COR   = PV_COR
DBA_INALL = PV_INALL

dba.overlap = function(DBA, mask, mode=DBA_OLAP_PEAKS, minVal=0,
                       contrast, method=DBA$config$AnalysisMethod, th=.1, bUsePval=FALSE, report,
                       byAttribute, bCorOnly=TRUE, CorMethod="pearson", 
                       DataType=DBA$config$DataType)
{                      
   DBA = pv.check(DBA)
   
   if( (mode == DBA_OLAP_ALL) | (!missing(contrast)) | (!missing(report)) ) {
   	
      if( (!missing(contrast)) | (!missing(report)) ) {
         
         if(missing(report)) {
            report   = dba.report(DBA,method=method, contrast=contrast,th=th,bUsePval=bUsePval,DataType=DBA_DATA_FRAME)
         } else {
         	if(class(report)!="data.frame") {
         	   stop('Class not supported for report parameter. Call dba.report with DataType=DBA_DATA_FRAME.')	
         	}
         	if(!missing(contrast)) {
         	   DBA = pv.getOverlapData(DBA,contrast,report)
         	}
         }
         
         sites = as.numeric(rownames(report))
         
         if(missing(mask)) {
         	if(missing(contrast)) {
         	   mask = 1:length(DBA$peaks)
         	} else {
         	   mask  = DBA$contrasts[[contrast]]$group1 | DBA$contrasts[[contrast]]$group2
            }   
         }  else if (length(mask)==1) {
            mask = 1:length(DBA$peaks)         
         }
         
         res = pv.occupancy(DBA, mask=mask, sites = sites, byAttribute=byAttribute, 
                            Sort='cor', CorMethod=CorMethod,
                            minVal=minVal, bCorOnly=bCorOnly)         
      
      } else {
   	
         res = pv.occupancy(DBA, mask=mask, byAttribute=byAttribute, 
                            Sort='cor', CorMethod=CorMethod,
                            minVal=minVal, bCorOnly=bCorOnly)
      } 
      
   }  else if(mode == DBA_OLAP_RATE) {
   
      res = pv.consensus(DBA,sampvec=mask,minOverlap=NULL)
   
   }  else {
   
      res = pv.overlap(DBA,mask=mask,minVal=minVal)
      
      if(DataType!=DBA_DATA_FRAME) {
         for(i in 1:length(res)) {
            res[[i]] = pv.peaks2DataType(res[[i]],DataType)
         }
      }   
   }                       
   
   return(res)
}
  
###############################################  
## dba.count -- count reads in binding sites ##
###############################################   

DBA_SCORE_RPKM                = PV_SCORE_RPKM
DBA_SCORE_RPKM_FOLD           = PV_SCORE_RPKM_FOLD
DBA_SCORE_READS               = PV_SCORE_READS
DBA_SCORE_READS_FOLD          = PV_SCORE_READS_FOLD
DBA_SCORE_READS_MINUS         = PV_SCORE_READS_MINUS
DBA_SCORE_TMM_MINUS_FULL      = PV_SCORE_TMM_MINUS_FULL
DBA_SCORE_TMM_MINUS_EFFECTIVE = PV_SCORE_TMM_MINUS_EFFECTIVE
DBA_SCORE_TMM_READS_FULL      = PV_SCORE_TMM_READS_FULL
DBA_SCORE_TMM_READS_EFFECTIVE = PV_SCORE_TMM_READS_EFFECTIVE

dba.count = function(DBA, peaks, minOverlap=2, score=DBA_SCORE_TMM_MINUS_EFFECTIVE, bLog=FALSE,
                     insertLength, maxFilter, bRemoveDuplicates=FALSE,
                     bCalledMasks=TRUE, bCorPlot=TRUE, bParallel=DBA$config$RunParallel) 
{
   DBA = pv.check(DBA)            
   
   bUseLast = F
  
   if(!missing(peaks)) {
      if(is.null(peaks)) {
         callers = unique(DBA$class[DBA_CALLER,])
         if((length(callers)==1) & (callers=='counts')) {
            DBA = pv.check(DBA)
            res = pv.setScore(DBA,score=score,bLog=bLog,minMaxval=maxFilter)
            return(res)	
         }	
      } else {
         peaks = pv.DataType2Peaks(peaks)
      }
   }
   
   if(missing(insertLength)) {
      insertLength=0
   }
  
   res = pv.counts(DBA, peaks=peaks, minOverlap=minOverlap, 
                   defaultScore=score, bLog=bLog, insertLength=insertLength, bOnlyCounts=T,
                   bCalledMasks=bCalledMasks, minMaxval=maxFilter, bParallel=bParallel, bUseLast=bUseLast, bWithoutDupes=bRemoveDuplicates)
   
   if(bCorPlot){
      x = dba.plotHeatmap(res,correlations=T)
   }

   if(class(res)!="DBA") {
      class(res) = "DBA"
   }
      
   return(res)
}

########################################################
## dba.contrast -- establish contrast(s) for analysis ##
########################################################

dba.contrast = function(DBA, group1, group2=!group1, name1="group1", name2="group2",
                        minMembers=3, block,
                        categories = c(DBA_TISSUE,DBA_FACTOR,DBA_CONDITION,DBA_TREATMENT))
{
   if(minMembers < 2) {
      stop('minMembers must be at least 2. Use of replicates strongly advised.')	
   }
   
   DBA = pv.check(DBA)
      
   res = pv.contrast(DBA, group1=group1, group2=group2, name1=name1, name2=name2,
                     minMembers=minMembers, categories=categories,block=block)

   if(class(res)!="DBA") {
      class(res) = "DBA"
   }
                        
   return(res)                       	
}

###################################################################
## dba.analyze -- perform differential binding affinity analysis ##
###################################################################

dba.analyze = function(DBA, method=DBA$config$AnalysisMethod, 
                       bSubControl=TRUE, bFullLibrarySize=FALSE, bTagwise=TRUE,
                       bCorPlot=TRUE, bReduceObjects=T, bParallel=DBA$config$RunParallel)
{
   
   #if(bParallel && DBA$config$parallelPackage==DBA_PARALLEL_MULTICORE) {
   #   warning('Parallel operation currently unreliable. Executing serially.')#,immediate.=TRUE)
   #   bParallel=F	
   #}
   	
   DBA = pv.check(DBA)
      
   res = pv.DBA(DBA, method ,bSubControl,bFullLibrarySize,bTagwise=bTagwise,minMembers=3,bParallel)
   
   if(bReduceObjects) {
      if(!is.null(res$contrasts)) {
         for(i in 1:length(res$contrasts)) {
            res$contrasts[[i]] = pv.stripDBA(res$contrasts[[i]])   	
         }
      }      	
   }
    
   if(bCorPlot){
   	  warn = T
   	  rep = pv.DBAreport(res,contrast=1,method=method[1],th=.1,bSupressWarning=T)
   	  if(!is.null(rep)) {
   	     if(!is.null(dim(rep))) {
   	        if(nrow(rep)>1) {
   	           warn=F
   	           x = dba.plotHeatmap(res,contrast=1,method=method[1],correlations=T)
   	        }	
   	     }
   	  }
   	  if(warn) {
   	     warning('No correlation heatmap plotted -- contrast 1 has no differentially bound sites.')	
   	  }
   }



   if(class(res)!="DBA") {
      class(res) = "DBA"
   }
      	
   return(res)
}

###########################################################
## dba.report -- generate report for a contrast analysis ##
###########################################################

dba.report = function(DBA, contrast=1, method=DBA$config$AnalysisMethod, th=.1, bUsePval=FALSE, 
                      fold=0, bNormalized=TRUE,
                      bCalled=FALSE, bCounts=FALSE, bCalledDetail=FALSE,
                      file,initString=DBA$config$reportInit,ext='csv',DataType=DBA$config$DataType) 
                     
{

   DBA = pv.check(DBA) 

   res = pv.DBAreport(pv=DBA,contrast=contrast,method=method,th=th,bUsePval=bUsePval,
                      bCalled=bCalled,bCounts=bCounts,bCalledDetail=bCalledDetail,
                      file=file,initString=initString,bNormalized=bNormalized,ext=ext,minFold=fold) 

   if(DataType!=DBA_DATA_FRAME) {
      res = pv.peaks2DataType(res,DataType)
   }
   
   return(res)	

}                      

################################################
## dba.plotHeatmap -- Heatmap with clustering ##
################################################

dba.plotHeatmap = function(DBA, attributes=DBA$attributes, maxSites=1000, minval, maxval,
                           contrast, method=DBA$config$AnalysisMethod, th=.1, bUsePval=FALSE, report, score,
                           mask, sites, sortFun,
                           correlations=TRUE, olPlot=DBA_COR, ColAttributes,RowAttributes, colSideCols, rowSideCols=colSideCols,
                           margin=10, colScheme="Greens", distMethod="pearson",
                           ...)
{
   DBA = pv.check(DBA)
   
   if(missing(contrast) && !missing(score)) {
      DBA = dba.count(DBA,peaks=NULL,score=score)	
   }
   
   if(!missing(contrast)) {
   	 if(!missing(report)) {
         report = pv.DataType2Peaks(report)
      }   	
      DBA = pv.getPlotData(DBA,attributes=attributes,contrast=contrast,report=report,
   	                       method=method,th=th,bUsePval=bUsePval,bNormalized=T,
   	                       bPCA=F,bLog=T,minval=minval,maxval=maxval)
   	  contrast = 1                     
   }
   	                          	  
   if(length(correlations)==1 & ((correlations[1] == DBA_OLAP_ALL) | (correlations[1] == TRUE)))  { 	
   	  correlations = pv.occupancy(DBA, mask=mask, sites=sites, Sort='cor', bCorOnly=T,CorMethod=distMethod) 
   }
   	  
   if(correlations[1]!=FALSE) {
      res = pv.plotHeatmap(DBA, attributes=attributes, overlaps=correlations, olPlot=olPlot,
                           ColScheme=colScheme, distMeth=distMethod, bReorder=TRUE, contrast=contrast,
                           RowAttributes=RowAttributes,ColAttributes=ColAttributes,rowSideCols=rowSideCols,colSideCols=colSideCols,
                           minval=minval, maxval=maxval, margins=c(margin,margin),
                           ...)
      return(res)
   }
      
   if(!missing(contrast)) {
      res = pv.plotHeatmap(DBA, numSites=maxSites, attributes=attributes, contrast=contrast,
                           RowAttributes=RowAttributes,ColAttributes=ColAttributes,rowSideCols=rowSideCols,colSideCols=colSideCols,
                           ColScheme=colScheme, distMeth=distMethod, 
                           margins=c(margin,margin),...)
      res = NULL
                      
   } else {
     
	  if(!missing(sortFun)) {
	     savevecs = DBA$vectors
		 DBA = pv.sort(DBA, sortFun, mask=mask)
      }
      
      res = pv.plotHeatmap(DBA, numSites=maxSites, attributes=attributes, mask=mask, sites=sites,
                           RowAttributes=RowAttributes,ColAttributes=ColAttributes,rowSideCols=rowSideCols,colSideCols=colSideCols,
                           minval=minval, maxval=maxval, ColScheme=colScheme, distMeth=distMethod, 
                           margins=c(margin,margin),...)
      res = NULL
         
	  if(!missing(sortFun)) {
         DBA$vectors = savevecs
	  }
   }
     
   return(res)	
}

#######################################################
## dba.plotPCA -- Principal Components Analysis plot ##
#######################################################

dba.plotPCA = function(DBA, attributes, minval, maxval,
                       contrast, method=DBA$config$AnalysisMethod, th=.1, bUsePval=FALSE, report, score,
                       mask, sites, cor=FALSE,
                       b3D=FALSE, vColors, dotSize, ...)
                       
{
   DBA = pv.check(DBA)
   
   if(missing(contrast) && !missing(score)) {
      DBA = dba.count(DBA,peaks=NULL,score=score)	
   }   
   
   if(!missing(contrast)){
   	  if(missing(attributes)) {
   	     attributes = DBA_GROUP
   	  }
   	  if(missing(dotSize)) {
   	     dotSize=NULL
   	  }
      if(!missing(report)) {
         report = pv.DataType2Peaks(report)
      }  	  
   	  DBA = pv.getPlotData(DBA,attributes=attributes,contrast=contrast,method=method,th=th,
   	                       bUsePval=bUsePval,report=report,bPCA=T,minval=minval,maxval=maxval)
   	  if(attributes[1] == PV_GROUP) {
   	     attributes = PV_ID
   	  }
   	  res = pv.plotPCA(DBA,attributes=attributes,size=dotSize,cor=cor,b3D=b3D,vColors=vColors,...)
   } else {
   	  if(missing(attributes)) {
   	     attributes = pv.attributePCA(DBA)
   	  }
   
      res = pv.plotPCA(DBA, attributes=attributes, size=dotSize, mask=mask, 
                       sites=sites, b3D=b3D, vColors=vColors, ...)  
   }

   return(res)	
}

#############################
## dba.plotBox --Boxplots  ##
#############################
                                      
dba.plotBox = function(DBA, contrast=1, method=DBA$config$AnalysisMethod, th=0.1, bUsePval=FALSE, bNormalized=TRUE,
                       attribute=DBA_GROUP, 
                       bAll=FALSE, bAllIncreased=FALSE, bAllDecreased=FALSE, 
                       bDB=TRUE, bDBIncreased=TRUE, bDBDecreased=TRUE,
                       pvalMethod=wilcox.test,  bReversePos=FALSE, attribOrder, 
                       vColors, varwidth=TRUE, notch=TRUE, ...) 

{
   DBA = pv.check(DBA)
   
   if(contrast > length(DBA$contrasts)) {
      stop('Supplied contrast greater than number of contrasts')	
   }
   
   res = pv.plotBoxplot(DBA, contrast=contrast, method=method, th=th, bUsePval=bUsePval, bNormalized=bNormalized,
                        attribute=attribute,bAll=bAll, bAllIncreased=bAllIncreased, bAllDecreased=bAllDecreased, 
                        bDB=bDB, bDBIncreased=bDBIncreased, bDBDecreased=bDBDecreased,
                        pvalMethod=pvalMethod,  bReversePos=bReversePos, attribOrder=attribOrder, vColors=vColors, 
                        varwidth=varwidth, notch=notch, ...)
 
   return(res)	
}

#########################################
## dba.plotMA -- MA or XY scatter plot ##
#########################################
                                      
dba.plotMA = function(DBA, contrast=1, method=DBA$config$AnalysisMethod, th=.1, bUsePval=FALSE, fold=0, bNormalized=TRUE,
                      factor="", bXY=FALSE, dotSize=.33, ...)

{
   DBA = pv.check(DBA)
   res = pv.DBAplotMA(DBA, contrast=contrast, method=method, bMA=!bXY, bXY=bXY, th=th, bUsePval=bUsePval, fold=fold,
                      facname=factor, bNormalized=bNormalized, cex=dotSize, ...)
 
   return(res)	
}
                                                                                                           
###########################################################
## dba.plotClust -- Hierarchical clister dengrogram plot ##
###########################################################

dba.plotClust = function(DBA, mask, sites, attributes=DBA$attributes, distMethod="pearson",
                         contrast, method=DBA$config$AnalysisMethod, th=.1, bUsePval=FALSE)
{
   
   if(!missing(contrast)) {
      
      message('dba.plotClust: Contrast clustering not implemented.')
      res = NULL
      
   } else {
   
      res = pv.plotClust(DBA, mask=mask, sites=sites, 
                         attributes = attributes, distMeth=distMethod)
   }

   return(res)
}                         

####################################################
## dba.plotVenn -- Venn diagram plots of overlaps ##
####################################################

dba.plotVenn = function(DBA, mask, overlaps, label1, label2, label3, ...)
{

   if(!missing(mask)) {
      overlaps = dba.overlap(DBA,mask,mode=DBA_OLAP_PEAKS,DataType=DBA_DATA_FRAME)
      
      res = pv.whichPeaksets(DBA,mask)
      if(missing(label1)) {
         label1 = DBA$class[PV_ID,res$A]
      }
      if(missing(label2)) {
         label2 = DBA$class[PV_ID,res$B]
      }
      if(missing(label3)) {
         label3 = DBA$class[PV_ID,res$C]
      }

   } else {
   
      if(missing(label1)) {
         label1 = "A"
      }
      if(missing(label2)) {
         label2 = "B"
      }
      if(missing(label3)) {
         label3 = "C"
      }
      for(i in 1:length(overlaps)) {
         overlaps[[i]] = pv.DataType2Peaks(overlaps[[i]])
      }
   }
         
   pv.plotVenn(overlaps,label1=label1,label2=label2,label3=label3,...)
  
}

###################################
## dba.show -- List DBA metadata ##
###################################

dba.show = function(DBA, mask, attributes, bContrasts=FALSE, th=0.1, bUsePval=FALSE) 
{
   DBA = pv.check(DBA)
   
   res = pv.list(DBA, mask=mask, bContrasts=bContrasts, attributes=attributes, th=th, bUsePval=bUsePval)
   
   return(res)
}

#######################################
## dba.mask -- mask samples or sites ##
#######################################

dba.mask = function(DBA, attribute, value, combine='or', mask, merge='or', bApply=FALSE,
                    peakset, minValue=-1)
                    
{
   DBA = pv.check(DBA)
      
   if(missing(peakset)) {
      
      res = pv.mask(DBA, attribute=attribute, value=value, combine=combine,
                         mask=mask, merge=merge, bApply=bApply)
      
   } else {
   
      res = pv.whichSites(DBA, pnum=peakset, combine=combine, minVal=minValue)
   
   }
   
   return(res)
 
 }
 
 ###################################
 ## dba.save -- save a DBA object ##
 ###################################

dba.save = function(DBA, file='DBA', dir='.', pre='dba_', ext='RData', bMinimize=FALSE)
{

   if(bMinimize) {
      DBA$allvectors = NULL
   }
   
   if(nrow(DBA$class)<DBA_TREATMENT) {
     DBA$class = rbind(DBA$class,'')
     rownames(DBA$class)[DBA_TREATMENT]='Treatment'	
   }
   
   DBA$vectors = NULL	
   DBA$values  = NULL
   DBA$hc      = NULL
   DBA$pc      = NULL

   DBA$config$lsfInit      = NULL
   DBA$config$parallelInit = NULL
   DBA$config$initFun      = NULL
   DBA$config$paramFun     = NULL
   DBA$config$addjobFun    = NULL
   DBA$config$lapplyFun    = NULL
   DBA$config$wait4jobsFun = NULL
   DBA$config$parallelInit = NULL

   DBA$config$Version1 = DBA_VERSION1
   DBA$config$Version2 = DBA_VERSION2
   DBA$config$Version3 = DBA_VERSION3
         
   res = pv.save(DBA,file=file ,
                 dir=dir, pre=pre, ext=ext,
                 compress=TRUE)
   
   return(res)
} 

 ###################################
 ## dba.load -- load a DBA object ##
 ###################################

dba.load = function(file='DBA', dir='.', pre='dba_', ext='RData')
{
	
   res = pv.load(file=file, dir=dir, pre=pre, ext=ext)
   
   if(is.null(res$allvectors)) {
      if(is.null(res$minOverlap)) {
         minOverlap=2
      } else {
         minOverlap = res$minOverlap
      }
      contrasts = res$contrasts
      res = dba(res,minOverlap=minOverlap)
      res$contrasts = contrasts
   }
   
   if(is.null(res$vectors)) {
      if(!is.null(res$overlapping)) {
         res$vectors = res$allvectors[res$overlapping,]
      } else {
         res$vectors = res$allvectors
      }
      if(!is.null(res$vectors)) {
         rownames(res$vectors) = 1:nrow(res$vectors)
      }
   }
   
   if(is.null(res$config$DataType)) {
   	  res$config$DataType=DBA_DATA_DEFAULT
      if(!is.null(res$config$RangedData)) {
         if(res$config$RangedData==F) {
            res$config$DataType=DBA_DATA_FRAME   	
         } 
      }
   }
   
   if(is.null(res$config$saveDir)) {
      res$config$saveDir = "Robjects"
   }

   if(is.null(res$config$savePrefix)) {
      res$config$saveDir = "dba_"
   }
   
   if(is.null(res$config$saveExt)) {
      res$config$saveDir = "RData"
   }

   if(is.null(res$config$reportInit)) {
      res$config$saveDir = "reports/DBA"
   }
   if(is.null(res$AnalysisMethod)){
      res$config$AnalysisMethod=DBA_EDGER
   }
   
   res$config$lsfInit      = NULL
   res$config$parallelInit = NULL
   res$config$initFun      = NULL
   res$config$paramFun     = NULL
   res$config$addjobFun    = NULL
   res$config$lapplyFun    = NULL
   res$config$wait4jobsFun = NULL
   res$config$parallelInit = NULL

   res = pv.version(res,DBA_VERSION1,DBA_VERSION2, DBA_VERSION3)
            
   if(class(res)!="DBA") {
      class(res) = "DBA"
   }
   
   return(res)
} 

##########################
## DBA object functions ##
##########################

print.DBA = function(x,...){
   #if(is.null(x$allvectors)) {
   #   print(dba.show(x))
   #} else {
   	  x = pv.check(x)
   	  cat(sprintf("%s:\n",summary(x)))
      print(dba.show(x))
      if(!is.null(x$contrasts)){
   	     cat("\n")

         if(length(x$contrasts) == 1) {
           cat(sprintf("1 Contrast:\n") )
   	     } else {
           cat(sprintf("%d Contrasts:\n",length(x$contrasts)))
         }
         print(dba.show(x,bContrasts=T))
      }
   #}
}

summary.DBA = function(object,...) {
   if(is.null(object$allvectors)) {
      cat('Run dba first\n')
      return
   }
   res = sprintf("%d Samples, %d sites in matrix",
          length(object$peaks),nrow(object$vectors))
   if(nrow(object$vectors) != nrow(object$allvectors)) {
      res = sprintf("%s (%d total)",res,nrow(object $allvectors))
   }
   return(res)
}

plot.DBA = function(x,...){
   DBA = pv.check(x)
   res = dba.plotHeatmap(x,...)
}
