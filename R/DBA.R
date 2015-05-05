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
DBA_VERSION2  = 14
DBA_VERSION3  = 2

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
DBA_INTERVALS = PV_INTERVALS
DBA_FRIP      = PV_SN_RATIO
DBA_ALL_ATTRIBUTES = c(DBA_ID,DBA_TISSUE,DBA_FACTOR,
                       DBA_CONDITION,DBA_TREATMENT,
                       DBA_REPLICATE,DBA_CALLER)

DBA_EDGER_CLASSIC = 'edgeR'
DBA_EDGER_BLOCK   = 'edgeRlm'
DBA_EDGER_GLM     = 'edgeRGLM'
DBA_EDGER         = DBA_EDGER_GLM

DBA_DESEQ_CLASSIC = 'DESeq1'
DBA_DESEQ_BLOCK   = 'DESeq1Block'
DBA_DESEQ_GLM     = 'DESeq1GLM'
DBA_DESEQ2        = 'DESeq2'
DBA_DESEQ2_BLOCK  = 'DESeq2Block'
DBA_DESEQ         = DBA_DESEQ_GLM

DBA_ALL_METHODS = c(DBA_EDGER,DBA_DESEQ,DBA_DESEQ2)
DBA_ALL_BLOCK   = c(DBA_EDGER_BLOCK,DBA_DESEQ_BLOCK,DBA_DESEQ2_BLOCK)
DBA_ALL_METHODS_BLOCK = c(DBA_ALL_METHODS, DBA_ALL_BLOCK)

DBA_DATA_FRAME                    = 0
DBA_DATA_RANGEDDATA               = 1
DBA_DATA_GRANGES                  = 2
DBA_DATA_SUMMARIZED_EXPERIMENT    = 3
DBA_DATA_DBAOBJECT                = 4
DBA_DATA_DEFAULT                  = DBA_DATA_GRANGES

dba = function(DBA,mask, minOverlap=2,
               sampleSheet="dba_samples.csv", 
               config=data.frame(RunParallel=TRUE, reportInit="DBA", DataType=DBA_DATA_GRANGES, 
                                 AnalysisMethod=DBA_EDGER, minQCth=15, fragmentSize=125, 
                                 bCorPlot=TRUE, th=.1, bUsePval=FALSE),
               peakCaller="raw", peakFormat, scoreCol, bLowerScoreBetter, filter, skipLines=0, bAddCallerConsensus=FALSE, 
               bRemoveM=TRUE, bRemoveRandom=TRUE, bSummarizedExperiment=FALSE,
               bCorPlot, attributes) 
{
    if(!missing(DBA)){
        if(class(DBA)=="character") {
            stop("DBA object is a character string; perhaps meant to be argument \'sampleSheet\'?")	
        }
        DBA = pv.check(DBA,bCheckSort=F)	
    } 
    
    res = pv.model(DBA, mask=mask, minOverlap=minOverlap, samplesheet=sampleSheet, config=config, 
                   caller=peakCaller, format=peakFormat, scorecol=scoreCol,bLowerBetter=bLowerScoreBetter,
                   skipLines=skipLines,bAddCallerConsensus=bAddCallerConsensus, 
                   bRemoveM=bRemoveM, bRemoveRandom=bRemoveRandom,
                   bKeepAll=TRUE, bAnalysis=TRUE, filter=filter,
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
    if(is.null(res$config$bCorPlot)){
        if(missing(bCorPlot)){
            res$config$bCorPlot=TRUE
        } else {
            res$config$bCorPlot=bCorPlot   
        }
    }
    if(is.null(res$config$th)){
        res$config$th=0.1
    }
    if(is.null(res$config$bUsePval)){
        res$config$bUsePval=FALSE
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
    
    if(missing(bCorPlot)) {
        try(dba.plotHeatmap(res),silent=TRUE)
    } else if(bCorPlot) {
        try(dba.plotHeatmap(res),silent=TRUE)      
    }
    
    if(bSummarizedExperiment) {
        res = pv.DBA2SummarizedExperiment(res)
    } else {
        if(!is.null(DBA$ChIPQCobj)) {
            #          resQC = DBA$ChIPQCobj
            #          resQC@DBA = res
            #          res = resQC
            warning('Returning new DBA object (not ChIPQCexperiment object)')
        }      
    }
    
    return(res)                 
}

###############################################
## dba.peakset -- add a peakset to the model ##
###############################################

dba.peakset = function(DBA=NULL, peaks, sampID, tissue, factor, condition, treatment, replicate,
                       control, peak.caller, peak.format, reads=0, consensus=FALSE, bamReads, bamControl,
                       scoreCol, bLowerScoreBetter, filter, counts, bRemoveM=TRUE, bRemoveRandom=TRUE,
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
    
    if(missing(DataType)) {
        DataType = DBA_DATA_DEFAULT
    }
    
    if(bRetrieve==TRUE || !missing(writeFile)) { ## RETRIEVE/WRITE PEAKSETS
        
        if(missing(writeFile)) {
            writeFile = NULL
        }
        
        if(missing(peaks)) {
            if(!is.null(DBA)) {
                DBA = pv.check(DBA,T)
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
        
    } else { ## ADD PEAKSET(S)
        
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
                             peak.caller=peak.caller, peak.format=peak.format, reads=reads, consensus=consensus, 
                             readBam=bamReads, controlBam=bamControl,
                             scoreCol=scoreCol, bLowerScoreBetter=bLowerScoreBetter, 
                             bRemoveM=bRemoveM, bRemoveRandom=bRemoveRandom,
                             minOverlap=minOverlap, filter=filter, counts=counts)
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
        if(is.null(res$config$th)){
            res$config$th=0.1
        }
        if(is.null(res$config$bUsePval)){
            res$config$bUsePval=FALSE
        }
        
        if(is.null(DBA$config$reportInit)){
            res$config$reportInit="DBA"
        } else {
            res$config$reportInit=DBA$config$reportInit
        }
        if(is.null(res$config$AnalysisMethod)){
            res$config$AnalysisMethod=DBA_EDGER
        }
        if(is.null(res$config$bCorPlot)){
            res$config$bCorPlot=TRUE
        } 
        
        if(bMerge) {
            res = pv.check(res)
        }
        
        if(!is.null(DBA$ChIPQCobj)) {
            #          resQC = DBA$ChIPQCobj
            #          resQC@DBA = res
            #          res = resQC
            warning('Returning new DBA object (not ChIPQCexperiment object)')
        }   
    }   
    
    return(res)                       
    
}                      

##################################################
## dba.overlap -- compute binding site overlaps ##
##################################################

DBA_OLAP_PEAKS = 1 # Return list of peaksets (common/unique peaks) 
DBA_OLAP_ALL   = 2 # Return overlap report with statstics for peakset pairs
DBA_OLAP_RATE  = 3 # Return vector of number of peaks in overlap for all values of minOverlap

DBA_OLAP  = PV_OLAP
DBA_COR   = PV_COR
DBA_INALL = PV_INALL

dba.overlap = function(DBA, mask, mode=DBA_OLAP_PEAKS, minVal=0,
                       contrast, method=DBA$config$AnalysisMethod, 
                       th=DBA$config$th, bUsePval=DBA$config$bUsePval, report,
                       byAttribute, bCorOnly=TRUE, CorMethod="pearson", 
                       DataType=DBA$config$DataType)
{                      
    DBA = pv.check(DBA,T)
    
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
            res = lapply(res,pv.peaks2DataType,DataType)
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
DBA_SCORE_TMM_MINUS_FULL_CPM      = PV_SCORE_TMM_MINUS_FULL_CPM
DBA_SCORE_TMM_MINUS_EFFECTIVE_CPM = PV_SCORE_TMM_MINUS_EFFECTIVE_CPM
DBA_SCORE_TMM_READS_FULL_CPM      = PV_SCORE_TMM_READS_FULL_CPM
DBA_SCORE_TMM_READS_EFFECTIVE_CPM = PV_SCORE_TMM_READS_EFFECTIVE_CPM
DBA_SCORE_SUMMIT              = PV_SCORE_SUMMIT
DBA_SCORE_SUMMIT_ADJ          = PV_SCORE_SUMMIT_ADJ
DBA_SCORE_SUMMIT_POS          = PV_SCORE_SUMMIT_POS

DBA_READS_DEFAULT = PV_READS_DEFAULT
DBA_READS_BAM     = PV_READS_BAM
DBA_READS_BED     = PV_READS_BED

dba.count = function(DBA, peaks, minOverlap=2, score=DBA_SCORE_TMM_MINUS_FULL, bLog=FALSE,
                     fragmentSize=DBA$config$fragmentSize, summits, filter=0, bRemoveDuplicates=FALSE, bScaleControl=TRUE,
                     mapQCth=DBA$config$mapQCth, filterFun=max, bCorPlot=DBA$config$bCorPlot, 
                     bUseSummarizeOverlaps=FALSE, readFormat=DBA_READS_DEFAULT,
                     bParallel=DBA$config$RunParallel) 
{
    DBA = pv.check(DBA,missing(peaks))
    
    #if(!missing(bLowMem)) {
    #   stop("parameter bLowMem deprecated, replaced with bUseSummarizeOverlaps")   
    #}
    
    if(minOverlap > length(DBA$peaks)) {
        stop(sprintf("minOverlap can not be greater than the number of peaksets [%s]",length(DBA$peaks)))	
    }           
    
    bUseLast = F
    
    res=NULL
    resetContrasts = TRUE
    if(missing(peaks) && !missing(summits) && !is.null(DBA$summits)) {
        if(DBA$summits == 0 ) {
            peaks=NULL
        } else {
            if(summits != DBA$summits) {
                stop('Can not change value of summits. Re-run from peaks.')
            } else {
                warning('No action taken, returning passed object...')
                res = DBA
            }
        }
    }
    if(!missing(peaks) || length(filter)>1) {
        if(is.null(peaks) || length(filter)>1) {
            callers = unique(DBA$class[DBA_CALLER,])
            if((length(callers)==1) & (callers=='counts')) {
                DBA = pv.check(DBA)
                if(!missing(summits)) {
                    if(summits>0) {
                        newpeaks = pv.Recenter(DBA,summits,DBA$sites)
                        res = pv.counts(DBA,peaks=newpeaks,
                                        defaultScore=score, bLog=bLog, insertLength=fragmentSize, bOnlyCounts=T,
                                        bCalledMasks=TRUE, minMaxval=filter, bParallel=bParallel, bUseLast=bUseLast,
                                        bWithoutDupes=bRemoveDuplicates,bScaleControl=bScaleControl,filterFun=filterFun,
                                        bLowMem=bUseSummarizeOverlaps,readFormat=readFormat,summits=0,
                                        minMappingQuality=mapQCth)
                        res$summits = summits
                    } else {
                        stop('Error: summits=0')
                    }
                } else {
                    if(length(filter)>1) {
                        DBA = pv.setScore(DBA,score=score,bLog=bLog,minMaxval=0,filterFun=filterFun)
                        res = pv.filterRate(DBA,filter,filterFun=filterFun)	
                    } else {
                        res = pv.setScore(DBA,score=score,bLog=bLog,minMaxval=filter,filterFun=filterFun)
                    }
                    resetContrasts=FALSE	
                }
            } else {
                stop('DBA object must contains only counts')	
            }	
        } else {
            peaks = pv.DataType2Peaks(peaks)
        }
    }
    
    if(is.null(res)) {
        res = pv.counts(DBA, peaks=peaks, minOverlap=minOverlap, 
                        defaultScore=score, bLog=bLog, insertLength=fragmentSize, bOnlyCounts=T,
                        bCalledMasks=TRUE, minMaxval=filter, bParallel=bParallel, bUseLast=bUseLast,
                        bWithoutDupes=bRemoveDuplicates,bScaleControl=bScaleControl,filterFun=filterFun,
                        bLowMem=bUseSummarizeOverlaps,readFormat=readFormat,summits=summits,
                        minMappingQuality=mapQCth)
        if(!missing(summits)) {
            res$summits = summits
        }
    }
    if(resetContrasts && length(res$contrasts)>0) {
        for(i in 1:length(res$contrasts)) {
            res$contrasts[[i]]$edgeR = NULL
            res$contrasts[[i]]$DESeq = NULL         	
        }
    }
    
    if(bCorPlot){
        try(dba.plotHeatmap(res,correlations=T),silent=TRUE)
    }
    
    if(class(res)!="DBA") {
        class(res) = "DBA"
    }
    
    if(!is.null(DBA$ChIPQCobj)) {
        res = checkQCobj(DBA$ChIPQCobj,res)
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
    
    DBA = pv.check(DBA,TRUE)
    
    res = pv.contrast(DBA, group1=group1, group2=group2, name1=name1, name2=name2,
                      minMembers=minMembers, categories=categories,block=block)
    
    if(class(res)!="DBA") {
        class(res) = "DBA"
    }
    
    if(!is.null(DBA$ChIPQCobj)) {
        res = checkQCobj(DBA$ChIPQCobj,res)
    }
    
    return(res)                       	
}

###################################################################
## dba.analyze -- perform differential binding affinity analysis ##
###################################################################

dba.analyze = function(DBA, method=DBA$config$AnalysisMethod, 
                       bSubControl=TRUE, bFullLibrarySize=TRUE, bTagwise=TRUE,
                       bCorPlot=DBA$config$bCorPlot, bReduceObjects=T, bParallel=DBA$config$RunParallel)
{
    
    #if(bParallel && DBA$config$parallelPackage==DBA_PARALLEL_MULTICORE) {
    #   warning('Parallel operation currently unreliable. Executing serially.')#,immediate.=TRUE)
    #   bParallel=F	
    #}
    
    DBA = pv.check(DBA,TRUE)
    
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
                    x = try(dba.plotHeatmap(res,contrast=1,method=method[1],correlations=T),silent=TRUE)
                }	
            }
        }
        if(warn) {
            warning('No correlation heatmap plotted -- contrast 1 does not have enough differentially bound sites.')	
        }
    }
    
    if(class(res)!="DBA") {
        class(res) = "DBA"
    }
    
    if(!is.null(DBA$ChIPQCobj)) {
        res = checkQCobj(DBA$ChIPQCobj,res)
    }
    
    return(res)
}

###########################################################
## dba.report -- generate report for a contrast analysis ##
###########################################################

dba.report = function(DBA, contrast, method=DBA$config$AnalysisMethod, th=DBA$config$th, bUsePval=DBA$config$bUsePval, 
                      fold=0, bNormalized=TRUE,
                      bCalled=FALSE, bCounts=FALSE, bCalledDetail=FALSE,
                      bDB, bNotDB, bAll=TRUE, bGain=FALSE, bLoss=FALSE,
                      file,initString=DBA$config$reportInit,ext='csv',DataType=DBA$config$DataType) 
    
{
    
    DBA = pv.check(DBA,TRUE) 
    
    if(DataType==DBA_DATA_SUMMARIZED_EXPERIMENT) {
        bCounts=T
    }
    
    if(!missing(bDB)|!missing(bNotDB)) {
        if(missing(bDB)) {
            bDB=FALSE
        }
        if(missing(bNotDB)) {
            bNotDB=FALSE
        }
        res = pv.resultsDBA(DBA,contrasts=contrast,methods=method,th=th,bUsePval=bUsePval,fold=fold,
                            bDB=bDB,bNotDB=bNotDB,bUp=bGain,bDown=bLoss,bAll=bAll)
        
        res$resultObject = TRUE
        return(res)                    
    }
    
    if(missing(contrast)) {
        contrast=1
    }
    
    res = pv.DBAreport(pv=DBA,contrast=contrast,method=method,th=th,bUsePval=bUsePval,
                       bCalled=bCalled,bCounts=bCounts,bCalledDetail=bCalledDetail,
                       file=file,initString=initString,bNormalized=bNormalized,ext=ext,minFold=fold) 
    
    if(DataType==DBA_DATA_SUMMARIZED_EXPERIMENT) {
        DBA = pv.getPlotData(DBA,contrast=contrast,report=res,
                             method=method,th=th,bUsePval=bUsePval,bNormalized=T)
        res = pv.DBA2SummarizedExperiment(DBA,report=res)
        return(res)
    }
    
    if(DataType!=DBA_DATA_FRAME) {
        res = pv.peaks2DataType(res,DataType)
    }
    
    return(res)	
    
}                      

################################################
## dba.plotHeatmap -- Heatmap with clustering ##
################################################

dba.plotHeatmap = function(DBA, attributes=DBA$attributes, maxSites=1000, minval, maxval,
                           contrast, method=DBA$config$AnalysisMethod, 
                           th=DBA$config$th, bUsePval=DBA$config$bUsePval, 
                           report, score, bLog=TRUE, mask, sites, sortFun, 
                           correlations=TRUE, olPlot=DBA_COR, ColAttributes,RowAttributes, colSideCols, rowSideCols=colSideCols,
                           margin=10, colScheme="Greens", distMethod="pearson",
                           ...)
{
    DBA = pv.check(DBA,TRUE)
    
    if( (missing(contrast) || !missing(mask)) && !missing(score) ) {
        DBA = dba.count(DBA,peaks=NULL,score=score,bCorPlot=FALSE)	
    }
    
    mask = pv.setMask(DBA,mask,contrast)
    
    if(!missing(contrast)) {
        if(!missing(report)) {
            report = pv.DataType2Peaks(report)
        }   	
        DBA = pv.getPlotData(DBA,attributes=attributes,contrast=contrast,report=report,
                             method=method,th=th,bUsePval=bUsePval,bNormalized=T,
                             bPCA=F,bLog=F,minval=minval,maxval=maxval,mask=mask)
        contrast = 1
        mask = NULL                     
    }
    
    if(bLog) {
        vectors = DBA$vectors[,4:ncol(DBA$vectors)]
        vectors[vectors<=0]=1
        vectors = log2(vectors)
        DBA$vectors[,4:ncol(DBA$vectors)] = vectors
        if(missing(minval)) {
            minval = 0
        } else {
            minval = max(0,minval)
        }
    }
    DBA$allvectors = DBA$vectors
    
    if(length(correlations)==1 & ((correlations[1] == DBA_OLAP_ALL) | (correlations[1] == TRUE)))  {
        if(nrow(DBA$allvectors)>1) {
            if(!missing(sites)) {
                if(is.logical(sites)) {
                    sites = which(sites)	
                }	   
            } 	
            correlations = pv.occupancy(DBA, mask=mask, sites=sites, Sort='cor', bCorOnly=T,CorMethod=distMethod)
        } else {
            warning('No correlation heatmap plotted -- contrast does not have enough differentially bound sites.')	
            return(NULL)   	     	
        }
    }
    
    if(correlations[1]!=FALSE) {
        res = pv.plotHeatmap(DBA, attributes=attributes, overlaps=correlations, olPlot=olPlot, mask=mask,
                             ColScheme=colScheme, distMeth=distMethod, bReorder=TRUE, contrast=contrast,
                             RowAttributes=RowAttributes,ColAttributes=ColAttributes,rowSideCols=rowSideCols,colSideCols=colSideCols,
                             minval=minval, maxval=maxval, margins=c(margin,margin),
                             ...)
    } else {
        
        if(!missing(contrast)) {
            if(nrow(DBA$allvectors)<2) { 	
                warning('No heatmap plotted -- contrast does not have enough differentially bound sites.')	
                return(NULL)   	     	
            }
            res = pv.plotHeatmap(DBA, numSites=maxSites, attributes=attributes, contrast=contrast,
                                 RowAttributes=RowAttributes,ColAttributes=ColAttributes,rowSideCols=rowSideCols,colSideCols=colSideCols,
                                 ColScheme=colScheme, distMeth=distMethod, 
                                 margins=c(margin,margin), ...)
            res = DBA$vectors[1:maxSites,][res$rowInd,c(1:3,3+res$colInd)]
            if(!is.character(res[1,1])) {
                res[,1] = DBA$chrmap[res[,1]]
            }
            res = as(res,"GRanges")
        } else {
            
            if(!missing(sortFun)) {
                savevecs = DBA$vectors
                DBA = pv.sort(DBA, sortFun, mask=mask)
            }
            
            res = pv.plotHeatmap(DBA, numSites=maxSites, attributes=attributes, mask=mask, sites=sites,
                                 RowAttributes=RowAttributes,ColAttributes=ColAttributes,rowSideCols=rowSideCols,colSideCols=colSideCols,
                                 minval=minval, maxval=maxval, ColScheme=colScheme, distMeth=distMethod, 
                                 margins=c(margin,margin),...)
            res = DBA$vectors[1:maxSites,][res$rowInd,c(1:3,3+res$colInd)]
            if(!is.character(res[1,1])) {
                res[,1] = DBA$chrmap[res[,1]]
            }
            res = as(res,"GRanges")
            if(!missing(sortFun)) {
                DBA$vectors = savevecs
            }
        }
    }
    
    invisible(res)	
}

#######################################################
## dba.plotPCA -- Principal Components Analysis plot ##
#######################################################

dba.plotPCA = function(DBA, attributes, minval, maxval,
                       contrast, method=DBA$config$AnalysisMethod, 
                       th=DBA$config$th, bUsePval=DBA$config$bUsePval, 
                       report, score, bLog=TRUE, mask, sites, label, cor=FALSE,
                       b3D=FALSE, vColors, dotSize, labelSize, labelCols, ...)
    
{
    DBA = pv.check(DBA,TRUE)
    
    mask = pv.setMask(DBA,mask,contrast)
    
    if(missing(contrast) && !missing(score)) {
        DBA = dba.count(DBA,peaks=NULL,score=score,bCorPlot=FALSE)	
    } else if (!missing(score)) {
        warning('score parameter ignored when contrast is specified')	
    }  
    
    if(missing(label)) {
        label=NULL
    } else {
        if(missing(labelSize)) {
            labelSize=.8
        }
        if(missing(labelCols)) {
            labelCols="black"
        }
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
        DBA = pv.getPlotData(DBA,attributes=attributes,contrast=contrast,report=report,
                             method=method,th=th,bUsePval=bUsePval,bNormalized=T,
                             bPCA=T,minval=minval,maxval=maxval,mask=mask,bLog=F)                     
        if(attributes[1] == PV_GROUP) {
            attributes = PV_ID
        }
        res = pv.plotPCA(DBA,attributes=attributes,size=dotSize,cor=cor,
                         b3D=b3D,vColors=vColors,label=label,bLog=bLog,
                         labelSize=labelSize,labelCols=labelCols,...)
    } else {
        if(missing(attributes)) {
            attributes = pv.attributePCA(DBA)
        }
        
        res = pv.plotPCA(DBA, attributes=attributes, size=dotSize, mask=mask, 
                         sites=sites, b3D=b3D, vColors=vColors, label=label,
                         bLog=bLog,labelSize=labelSize,labelCols=labelCols,...)  
    }
    
    invisible(res)	
}

#############################
## dba.plotBox --Boxplots  ##
#############################

dba.plotBox = function(DBA, contrast=1, method=DBA$config$AnalysisMethod, 
                       th=DBA$config$th, bUsePval=DBA$config$bUsePval, 
                       bNormalized=TRUE, attribute=DBA_GROUP, 
                       bAll=FALSE, bAllIncreased=FALSE, bAllDecreased=FALSE, 
                       bDB=TRUE, bDBIncreased=TRUE, bDBDecreased=TRUE,
                       pvalMethod=wilcox.test,  bReversePos=FALSE, attribOrder, 
                       vColors, varwidth=TRUE, notch=TRUE, ...) 
    
{
    DBA = pv.check(DBA,TRUE)
    
    if(contrast > length(DBA$contrasts)) {
        stop('Supplied contrast greater than number of contrasts')	
    }
    
    res = pv.plotBoxplot(DBA, contrast=contrast, method=method, th=th, bUsePval=bUsePval, bNormalized=bNormalized,
                         attribute=attribute,bAll=bAll, bAllIncreased=bAllIncreased, bAllDecreased=bAllDecreased, 
                         bDB=bDB, bDBIncreased=bDBIncreased, bDBDecreased=bDBDecreased,
                         pvalMethod=pvalMethod,  bReversePos=bReversePos, attribOrder=attribOrder, vColors=vColors, 
                         varwidth=varwidth, notch=notch, ...)
    
    invisible(res)	
}

#########################################
## dba.plotMA -- MA or XY scatter plot ##
#########################################

dba.plotMA = function(DBA, contrast=1, method=DBA$config$AnalysisMethod, 
                      th=DBA$config$th, bUsePval=DBA$config$bUsePval, 
                      fold=0, bNormalized=TRUE,
                      factor="", bXY=FALSE, dotSize=.45, bSignificant=TRUE, bSmooth=TRUE, ...)
    
{
    DBA = pv.check(DBA,TRUE)
    
    res = pv.DBAplotMA(DBA, contrast=contrast, method=method, bMA=!bXY, bXY=bXY, th=th, bUsePval=bUsePval, fold=fold,
                       facname=factor, bNormalized=bNormalized, cex=dotSize, 
                       bSignificant = bSignificant, bSmooth=bSmooth,  ...)
    
    invisible(res)
}

###########################################################
## dba.plotClust -- Hierarchical clister dengrogram plot ##
###########################################################

dba.plotClust = function(DBA, mask, sites, attributes=DBA$attributes, distMethod="pearson",
                         contrast, method=DBA$config$AnalysisMethod, th=DBA$config$th, bUsePval=DBA$config$bUsePval)
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

dba.plotVenn = function(DBA, mask, overlaps, label1, label2, label3, label4, main, sub, 
                        contrast, method=DBA$config$AnalysisMethod, 
                        th=DBA$config$th, bUsePval=DBA$config$bUsePval,
                        bDB=TRUE, bNotDB, bAll=TRUE, bGain=FALSE, bLoss=FALSE,
                        labelAttributes, bReturnPeaksets=FALSE, DataType=DBA$config$DataType)
{
    DBA = pv.check(DBA,TRUE)
    
    newDBA = NULL
    
    if (!missing(overlaps)){
        
        if(missing(label1)) {
            label1 = "A"
        }
        if(missing(label2)) {
            label2 = "B"
        }
        if(missing(label3)) {
            label3 = "C"
        }
        if(missing(label4)) {
            label4 = "D"
        }
        for(i in 1:length(overlaps)) {
            overlaps[[i]] = pv.DataType2Peaks(overlaps[[i]])
        }
        
    } else if(!missing(contrast)){
        if(max(contrast)>length(DBA$contrasts)) {
            stop('Contrast greater than number of contrasts.')   
        }
        newDBA = dba.report(DBA,contrast=contrast,method=method,th=th,bUsePval=bUsePval,
                            bDB=bDB, bNotDB=bNotDB,bAll=bAll, bGain=bGain, bLoss=bLoss)
        
        if(is.null(newDBA$peaks)){
            stop('No peaksets meet specified criteria.')   
        }
        if(missing(mask)) {
            mask = 1:length(newDBA$peaks)
            if(length(mask)>4){
                stop('Too many peaksets meet specified criteria.')
            }
            if(length(mask)==1) {
                stop('Only one peakset meets specified criteria.')
            }
        } else {
            if(is.logical(mask)) {
                if(length(mask)!=length(newDBA$peaks)) {
                    stop('Logical mask doe not have same number of elements as there are peaksets.')
                }
                if(sum(mask)>4) {
                    stop('Too many peaksets in mask.')
                } else if(length(mask)>4){
                    stop('Too many peaksets in mask.')
                }
                mask = which(mask)
            }
            if(length(mask)>length(newDBA$peaks) | max(mask)>length(newDBA$peaks) ) {
                stop('Peakset specified in mask is out of range.')
            }
        }
        
        overlaps = dba.overlap(newDBA,mask,mode=DBA_OLAP_PEAKS,DataType=DBA_DATA_FRAME)
        
        res = pv.whichPeaksets(newDBA,mask)
        
        if(missing(labelAttributes)) {
            labelAttributes=c(DBA_ID,DBA_FACTOR,DBA_TISSUE,DBA_CONDITION,DBA_TREATMENT)  
        }
        crec = matrix(newDBA$class[labelAttributes,mask],length(labelAttributes),length(mask))
        if(length(mask)==2) {
            labels = pv.namestrings(crec[,1],crec[,2])
        } else if (length(mask)==3) {
            labels = pv.namestrings(crec[,1],crec[,2],crec[,3])         
        } else {
            labels = pv.namestrings(crec[,1],crec[,2],crec[,3],crec[,4])         
        }
        
        if(missing(label1)) {
            label1 = labels$n1
        }
        if(missing(label2)) {
            label2 = labels$n2
        }
        if(missing(label3)) {
            label3 = labels$n3
        }
        if(missing(label4)) {
            label4 = labels$n4
        }
        
        if (missing(sub)) {
            sub = labels$tstring
        }
        
    } else if(!missing(mask)) {
        
        if(is.logical(mask)) {
            if(length(mask)!=length(DBA$peaks)) {
                stop('Logical mask doe not have same number of elements as there are peaksets.')
            }
            if(sum(mask)>4) {
                stop('Too many peaksets in mask.')
            }
            mask = which(mask)
        } else if(length(mask)>4){
            stop('Too many peaksets in mask.')
        }      
        
        overlaps = dba.overlap(DBA,mask,mode=DBA_OLAP_PEAKS,DataType=DBA_DATA_FRAME)
        
        res = pv.whichPeaksets(DBA,mask)
        if(missing(labelAttributes)) {
            if(pv.checkValue(DBA$resultObject,TRUE)) {
                labelAttributes=c(DBA_ID,DBA_FACTOR,DBA_TISSUE,DBA_CONDITION,DBA_TREATMENT)  
            } else {
                labelAttributes=DBA_ID
            }
        }
        crec = matrix(DBA$class[labelAttributes,mask],length(labelAttributes),length(mask))
        if(length(mask)==2) {
            labels = pv.namestrings(crec[,1],crec[,2])
        } else if (length(mask)==3) {
            labels = pv.namestrings(crec[,1],crec[,2],crec[,3])         
        } else {
            labels = pv.namestrings(crec[,1],crec[,2],crec[,3],crec[,4])         
        }
        if(missing(label1)) {
            label1 = labels$n1
        }
        if(missing(label2)) {
            label2 = labels$n2
        }
        if(missing(label3)) {
            label3 = labels$n3
        }
        if(missing(label4)) {
            label4 = labels$n4
        }
        if (missing(sub)) {
            sub = labels$tstring
        }
    } else {
        stop("Must specify one of mask, overlaps, or contrast.")
    }
    
    if(missing(main)) {
        main = "Binding Site Overlaps"
    }
    if (missing(sub)) {
        sub = ""
    }
    
    pv.plotVenn(overlaps,label1=label1,label2=label2,label3=label3,label4=label4,main,sub)
    
    if(bReturnPeaksets) {
        if(DataType == DBA_DATA_DBAOBJECT) {
            if(!is.null(newDBA)) {
                return(newDBA)
            } else {
                warning('No DBA object to return.')
            }
        } else {
            for(i in 1:length(overlaps)) {
                overlaps[[i]] = pv.peaks2DataType(overlaps[[i]],DataType)
            }
            return(overlaps)
        }   
    }   
}

###################################
## dba.show -- List DBA metadata ##
###################################

dba.show = function(DBA, mask, attributes, bContrasts=FALSE, 
                    th=DBA$config$th, bUsePval=DBA$config$bUsePval) 
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
    
    if(class(DBA)=="ChIPQCexperiment") {
        saveChIPQC = DBA
        DBA = DBA@DBA
    } else saveChIPQC = NULL
    
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
    
    if(!is.null(saveChIPQC)) {
        saveChIPQC@DBA = DBA
        DBA = saveChIPQC
    }
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
    
    if(class(res)=="ChIPQCexperiment") {
        saveChIPQC = res
        res = res@DBA
    } else saveChIPQC = NULL
    
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
    
    if(is.null(res$bCorPlot)){
        res$config$bCorPlot=TRUE
    }
    
    if(is.null(res$config$th)){
        res$config$th=0.1
    }
    if(is.null(res$config$bUsePval)){
        res$config$bUsePval=FALSE
    }   
    
    res$config$lsfInit      = NULL
    res$config$parallelInit = NULL
    res$config$initFun      = NULL
    res$config$paramFun     = NULL
    res$config$addjobFun    = NULL
    res$config$lapplyFun    = NULL
    res$config$wait4jobsFun = NULL
    res$config$parallelInit = NULL
    
    res$config = as.list(res$config)
    
    if (is.null(res$config$mapQCth)) {
        res$config$mapQCth=15   
    }
    
    if (is.null(res$config$fragmentSize)) {
        res$config$fragmentSize=125
    }   
    
    res = pv.version(res,DBA_VERSION1,DBA_VERSION2, DBA_VERSION3)
    
    if(class(res)!="DBA") {
        class(res) = "DBA"
    }
    
    if(!is.null(saveChIPQC)) {
        saveChIPQC@DBA = res
        res = saveChIPQC
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
        print(dba.show(x,bContrasts=T,th=x$config$th,bUsePval=x$config$bUsePval))
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
