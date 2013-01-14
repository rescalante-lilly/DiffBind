##########€€€€€######################
## pv_helper.R -- support for pv   ##
## 20 October 2009                 ##
## 3 February 2011 -- packaged     ##
## Rory Stark                      ##
## Cancer Research UK              ##
#####################################

#########################
## pv HELPER FUNCTIONS ##
#########################

pv.check = function(pv) {
   
   if(missing(pv)) {
      stop('DBA object missing!',call.=F)
   }
   
   if(is.null(pv)) {
      return(NULL)	
   }
   if(is.null(pv$vectors)) {
   	  if(is.null(pv$minOverlap)) {
   	     minOverlap=2
   	  } else {
   	     minOverlap = pv$minOverlap
   	  }
      pv = pv.vectors(pv,minOverlap=minOverlap)
   }
   if(is.null(pv$config$DataType)) {
   	  pv$config$DataType=DBA_DATA_DEFAULT
      if(!is.null(pv$config$RangedData)) {
         if(pv$config$RangedData==F) {
            pv$config$DataType=DBA_DATA_FRAME   	
         } 
      }
   }
   if(nrow(pv$class)<PV_TREATMENT) {
     pv$class = rbind(pv$class,'')
     rownames(pv$class)[PV_TREATMENT]='Treatment'	
   }
   
   return(pv)
}

pv.version = function(pv,v1,v2,v3){

   warn = FALSE
   if(is.null(pv$config$Version1)) {
      warn = T   	
   } else {
      if(pv$config$Version1 != v1) {
        warn = T   	
      }	
   }
   if(is.null(pv$config$Version2)) {
      warn = T   	
   } else {
      if(pv$config$Version2 != v2) {
        warn = T   	
      }	
   }
   
   if(warn) {
      warning('Loading DBA object from a previous version -- updating...',call.=F)
   }
   
   if(nrow(pv$class)<PV_TREATMENT) {
     pv$class = rbind(pv$class,'')
     rownames(pv$class)[PV_TREATMENT]='Treatment'	
   }
   pv$config$Version1 = v1
   pv$config$Version2 = v2
   pv$config$Version3 = v3   	
   
   return(pv)
}

pv.setScore = function(pv,score,bLog=F,minMaxval) {

   if(!is.null(pv$score)) {
      if(pv$score == score) {
         if(!is.null(pv$maxFilter)) {
            if(pv$maxFilter == minMaxval) {
               return(pv)	
            }	
         }
      }	
   }
   	
   if ((score >= DBA_SCORE_TMM_MINUS_FULL) && (score <= DBA_SCORE_TMM_READS_EFFECTIVE) ) {

      if(score == DBA_SCORE_TMM_MINUS_FULL) {
         bMinus   = TRUE
         bFullLib = TRUE	
      }
      if(score == DBA_SCORE_TMM_MINUS_EFFECTIVE) {
         bMinus   = TRUE
         bFullLib = FALSE
      }
      if(score == DBA_SCORE_TMM_READS_FULL) {
         bMinus   = FALSE
         bFullLib = TRUE	
      }
      if(score == DBA_SCORE_TMM_READS_EFFECTIVE) {
         bMinus   = FALSE
         bFullLib = FALSE	
      }
      
      pv$allvectors[,4:ncol(pv$allvectors)] = pv.normTMM(pv,bMinus=bMinus,bFullLib=bFullLib)
   
   } else {	
   	
  	   for(i in 1:length(pv$peaks)) {
	      colnum = 3+i
	      if(score == PV_SCORE_RPKM) {
	         pv$allvectors[,colnum] = pv$peaks[[i]]$RPKM	
	      }   		
	      if(score == PV_SCORE_RPKM) {
	         pv$allvectors[,colnum] = pv$peaks[[i]]$RPKM	
	      } else if(score == PV_SCORE_RPKM_FOLD) {
	         pv$allvectors[,colnum] = pv$peaks[[i]]$RPKM/pv$peaks[[i]]$cRPKM
	         if(bLog) {
	           pv$allvectors[,colnum] = log2(pv$allvectors[,colnum])	
	         }
	      } else if(score == PV_SCORE_READS) {
	         pv$allvectors[,colnum] = pv$peaks[[i]]$Reads	
	      } else if(score == PV_SCORE_READS_FOLD) {
	         pv$allvectors[,colnum] = pv$peaks[[i]]$Reads/pv$peaks[[i]]$cReads	
	         if(bLog) {
	           pv$allvectors[,colnum] = log2(pv$allvectors[,colnum])	
	         }
	      } else if(score == PV_SCORE_READS_MINUS) {
	         pv$allvectors[,colnum] = pv$peaks[[i]]$Reads-pv$peaks[[i]]$cReads	
	      }  
	   }
   }
      
   pv$vectors = pv$allvectors
       
   if(!missing(minMaxval)) {
      data = pv$allvectors[,4:ncol(pv$allvectors)]
      maxs = apply(pv$allvectors[,4:ncol(pv$allvectors)],1,max)
      tokeep = maxs>=minMaxval
      if(sum(tokeep)>1) {
         pv$allvectors = pv$allvectors[tokeep,]
         rownames(pv$allvectors) = 1:sum(tokeep)
         pv$vectors = pv$allvectors
         for(i in 1:length(pv$peaks)) {
               pv$peaks[[i]] = pv$peaks[[i]][tokeep,]
               rownames(pv$peaks[[i]]) = 1:sum(tokeep)
         }
         pv$overlapping = pv$overlapping[tokeep]
         pv = pv.vectors(pv,minOverlap=1,bAnalysis=F,bAllSame=T)
      } else {
         stop('No sites have activity greater than maxFilter')
      }
      if(!is.null(pv$contrasts)) {
         for(i in 1:length(pv$contrasts)) {
            pv$contrasts[[i]]$edgeR=NULL
            pv$contrasts[[i]]$DESeq=NULL
         }
      }
      pv$maxFilter = minMaxval
   }
   
   pv$score = score
      
   return(pv)
}

pv.whichPeaksets = function(pv,mask) {

   if(missing(mask)) {
      warning('mask required',call.=F)
      return(NULL)
   }
   if(is.null(mask)) {
      warning('mask required',call.=F)
      return(NULL)
   }
   if(class(mask)=='logical') {
      mask = which(mask)
   }
   
   A = mask[1]
   B = mask[2]
   if(length(mask) >= 3) {
      C = mask[3]
   } else {
      C = NULL
   }
   if(length(mask) >= 4) {
      D = mask[4]
   } else {
      D = NULL
   }   

   return(list(A=A,B=B,C=C,D=D))
}

pv.listadd = function(a,b){
   b = list(b)
   if (is.null(a)) return(b)
   return(c(a,b))
}

pv.listaddto = function(a,b){
   if (is.null(a)) return(b)
   return(c(a,b))
}

fdebug = function(str,file='debug.txt'){

   PV_DEBUG=FALSE
   
   if(PV_DEBUG == FALSE){
      return
   }
   
   #write(sprintf('%s\n',str),file=file,append=T)

}

pv.peaksort = function(peaks){
   p2 = peaks[,2]
   o1 = order(p2,decreasing=F)
   peaks = peaks[o1,]
   p1 = peaks[,1]
   o2 = order(p1,decreasing=F)
   peaks = peaks[o2,]
   return(peaks)
}

pv.contrast2 = function(vectors,n1,n2,minVal=0,v1,v2){
   
   if(missing(v1)) v1 = vectors[,n1+3] > minVal
   if(missing(v2)) v2 = vectors[,n2+3] > minVal
   
   allpeaks = v1 | v2
   inAll    = v1 & v2
   onlyA    = v1 & !v2
   onlyB    = !v1 & v2
   
   res = list(onlyA = vectors[onlyA,c(1:3,(n1+3))],
   	          onlyB = vectors[onlyB,c(1:3,(n2+3))],
   	          inAll = vectors[inAll,c(1:3,(n1+3),(n2+3))])
   
   for(i in 1:3){
      if(is.null(nrow(res[[i]]))){
         res[[i]] = as.matrix(t(res[[i]]))
      }
   }
   
   colnames(res[[1]]) = c("chr","start","end","score")	
   colnames(res[[2]]) = c("chr","start","end","score")	
   colnames(res[[3]]) = c("chr","start","end","scoreA","scoreB")	
   
   return(res)   	
}

pv.contrast3 = function(vectors,n1,n2,n3,minVal=0,v1,v2,v3){
	
   if(missing(v1)) v1 = vectors[,n1+3] > minVal
   if(missing(v2)) v2 = vectors[,n2+3] > minVal
   if(missing(v3)) v3 = vectors[,n3+3] > minVal
   allpeaks = v1 | v2 | v3
   
   
   inAll  =  v1 &  v2 &  v3
   onlyA  =  v1 & !v2 & !v3
   onlyB  = !v1 &  v2 & !v3 
   onlyC  = !v1 & !v2 &  v3
   
   notA   = !v1 &  (v2 & v3)
   notB   = !v2 &  (v1 & v3)
   notC   = !v3 &  (v1 & v2)    
	    
   res = list(onlyA = vectors[onlyA,c(1:3,(n1+3))],
   	          onlyB = vectors[onlyB,c(1:3,(n2+3))],
   	          onlyC = vectors[onlyC,c(1:3,(n3+3))],
   	          notA  = vectors[notA,c(1:3,(n2+3),(n3+3))],
   	          notB  = vectors[notB,c(1:3,(n1+3),(n3+3))],
   	          notC  = vectors[notC,c(1:3,(n1+3),(n2+3))],
   	          inAll = vectors[inAll,c(1:3,(n1+3),(n2+3),(n3+3))])      
   
   for(i in 1:7){
      if(is.null(nrow(res[[i]]))){
         res[[i]] = as.matrix(t(res[[i]]))
      }
   }	

   colnames(res[[1]]) = c("chr","start","end","score")	
   colnames(res[[2]]) = c("chr","start","end","score")	
   colnames(res[[3]]) = c("chr","start","end","score")	
   colnames(res[[4]]) = c("chr","start","end","scoreB","scoreC")
   colnames(res[[5]]) = c("chr","start","end","scoreA","scoreC")
   colnames(res[[6]]) = c("chr","start","end","scoreA","scoreB")
   colnames(res[[7]]) = c("chr","start","end","scoreA","scoreB","scoreC")
           
   return(res)
}

pv.contrast4 = function(vectors,n1,n2,n3,n4,minVal=0,v1,v2,v3,v4){
	
   if(missing(v1)) v1 = vectors[,n1+3] > minVal
   if(missing(v2)) v2 = vectors[,n2+3] > minVal
   if(missing(v3)) v3 = vectors[,n3+3] > minVal
   if(missing(v4)) v4 = vectors[,n4+3] > minVal
   allpeaks = v1 | v2 | v3 | v4
   
   
   inAll  =  v1 &  v2 &  v3 & v4
   onlyA  =  v1 & !v2 & !v3 & !v4
   onlyB  = !v1 &  v2 & !v3 & !v4
   onlyC  = !v1 & !v2 &  v3 & !v4
   onlyD  = !v1 & !v2 & !v3 & v4
   
   notA   = !v1 &  (v2 & v3 & v4)
   notB   = !v2 &  (v1 & v3 & v4)
   notC   = !v3 &  (v1 & v2 & v4)
   notD   = !v4 &  (v1 & v2 & v3)      
   
   AandB  = (v1 & v2) & !(v3 | v4)
   AandC  = (v1 & v3) & !(v2 | v4)
   AandD  = (v1 & v4) & !(v2 | v3)
   BandC  = (v2 & v3) & !(v1 | v4)
   BandD  = (v2 & v4) & !(v1 | v3)
   CandD  = (v3 & v4) & !(v1 | v2)             
	    
   res = list(onlyA = vectors[onlyA,c(1:3,(n1+3))],
   	          onlyB = vectors[onlyB,c(1:3,(n2+3))],
   	          onlyC = vectors[onlyC,c(1:3,(n3+3))],
   	          onlyD = vectors[onlyD,c(1:3,(n4+3))],   	          
   	          AandB = vectors[AandB,c(1:3,(n1+3),(n2+3))],
   	          AandC = vectors[AandC,c(1:3,(n1+3),(n3+3))],
   	          AandD = vectors[AandD,c(1:3,(n1+3),(n4+3))],
   	          BandC = vectors[BandC,c(1:3,(n2+3),(n3+3))],
   	          BandD = vectors[BandD,c(1:3,(n3+3),(n4+3))],
   	          CandD = vectors[CandD,c(1:3,(n3+3),(n4+3))],
   	          notA  = vectors[notA,c(1:3,(n2+3),(n3+3),(n4+3))],
   	          notB  = vectors[notB,c(1:3,(n1+3),(n3+3),(n4+3))],
   	          notC  = vectors[notC,c(1:3,(n1+3),(n2+3),(n4+3))],
   	          notD  = vectors[notD,c(1:3,(n1+3),(n2+3),(n3+3))],
   	          inAll = vectors[inAll,c(1:3,(n1+3),(n2+3),(n3+3),(n4+3))])      
   
   for(i in 1:15){
      if(is.null(nrow(res[[i]]))){
         res[[i]] = as.matrix(t(res[[i]]))
      }
   }	

   colnames(res[[1]])  = c("chr","start","end","score")	
   colnames(res[[2]])  = c("chr","start","end","score")	
   colnames(res[[3]])  = c("chr","start","end","score")	
   colnames(res[[4]])  = c("chr","start","end","score")
   colnames(res[[5]])  = c("chr","start","end","scoreA","scoreB")
   colnames(res[[6]])  = c("chr","start","end","scoreA","scoreC")
   colnames(res[[7]])  = c("chr","start","end","scoreA","scoreD")
   colnames(res[[8]])  = c("chr","start","end","scoreB","scoreC")   
   colnames(res[[9]])  = c("chr","start","end","scoreB","scoreD") 
   colnames(res[[10]]) = c("chr","start","end","scoreC","scoreD")
   colnames(res[[11]]) = c("chr","start","end","scoreB","scoreC","scoreD")    
   colnames(res[[12]]) = c("chr","start","end","scoreA","scoreC","scoreD")       
   colnames(res[[13]]) = c("chr","start","end","scoreA","scoreB","scoreD")    
   colnames(res[[14]]) = c("chr","start","end","scoreA","scoreB","scoreC")      
   colnames(res[[15]]) = c("chr","start","end","scoreA","scoreB","scoreC","scoreD")
           
   return(res)
}



pv.analysis = function(pv,attributes=pv$attributes,bPCA=T,distMeth="pearson") {

   #require(amap)
   
   peaks  = pv$vectors 
   values = matrix(as.numeric(as.matrix(peaks[,4:ncol(peaks)])),nrow(peaks),ncol(peaks)-3)
    
   if(sum(pv$vectors[,4:ncol(pv$vectors)] != 1) == 0) {
      return(pv)
   }
   
   if(sum(pv$vectors[,4:ncol(pv$vectors)] != -1) == 0) {
      return(pv)
   }
      
   if(nrow(peaks) <= length(pv$peaks) ) {
      bPCA=F
   }
   #values = unique(values)
   cnames=NULL
   for(i in 1:ncol(pv$class)) {
      cnames = c(cnames,pv.namestrings(pv$class[attributes,i])$tstring)
   }
   colnames(values) = cnames
   x = apply(values,1,pv.howmany)
   values = values[x>1,]
   #pv$hc = hclust(Dist(t(values),method=distMeth))
   if(bPCA) pv$pc = princomp(values)
   #pv$values = values

   return(pv)  
}

pv.howmany = function(vals){
   return(sum(vals>0))	
}


pv.readPeaks = function(peaks,peak.format,skipLines=0){
   if(peak.format == "macs") {
      peaks = pv.macs(peaks)
   } 
   else if(peak.format == "bayes") {
      peaks = pv.bayes(peaks)
   }
   else if(peak.format == "swembl") {
      peaks = pv.swembl(peaks)
   } 
   else if(peak.format == "raw") {
      peaks = pv.readbed(peaks)
   } 
   else if(peak.format == "fp4") {
      peaks =  pv.readbed(peaks,1)
   } 
   else if(peak.format == "bed") {
      peaks =  pv.readbed(peaks,skipLines)
   } 
   else if(peak.format == "tpic") {
      peaks =  pv.tpic(peaks)
   } 
   else if(peak.format == "sicer") {
      peaks =  pv.sicer(peaks)
   }    
   else if(peak.format == "narrow") {
      peaks =  pv.readbed(peaks,skipLines)
   } 
   else if(peak.format == "raw") {
      peaks =  pv.readbed(peaks,skipLines)
   } else {
      peaks =  pv.readbed(peaks,skipLines)   	
   }     
}


pv.defaultScoreCol = function(peak.format){
   if(peak.format == "macs") {
      val = 7
   } 
   else if(peak.format == "bayes") {
      val = 0
   }
   else if(peak.format == "swembl") {
      val = 4
   } 
   else if(peak.format == "raw") {
      val = 4
   } 
   else if(peak.format == "fp4") {
      val = 5
   } 
   else if(peak.format == "bed") {
      val = 5
   } 
   else if(peak.format == "tpic") {
      val = 0
   } 
   else if(peak.format == "sicer") {
      val = 7
   }    
   else if(peak.format == "narrow") {
      val = 8
   } 
   else if(peak.format == "raw") {
      val = 4
   } else {
      val = 4 	
   } 
   return(val)    
}


FDRth=100
pv.macs = function(fn){
 data = read.table(fn,blank.lines.skip=T,header=T)
 res  = pv.peaksort(data)
 return(res)
}

pv.swembl = function(fn){
 data = read.table(fn,skip=14)
 res  = pv.peaksort(data)
 return(res) 
}

pv.readbed = function(fn,skipnum=0){
 data = read.table(fn,skip=skipnum)
 res  = pv.peaksort(data)
 if(ncol(res)==3) {
    res=cbind(res,1)	
 }
 return(res)
}

pv.bayes = function(fn){
 data = read.table(fn)
 idx = data[,4]>0.5
 data = data[idx,]
 res  = pv.peaksort(data)
 return(res)
}

pv.tpic = function(fn){
 data = read.table(fn)
 res  = pv.peaksort(data)
 return(cbind(res,1))
}

pv.sicer = function(fn){
 data = read.table(fn)
 res  = pv.peaksort(data)
 return(res[,c(1:3,7)]) 
}

pv.sourcedata = function(fn,maxval){
 data = read.table(fn)
 vals = data[,6]
 if(!missing(maxval)) {
    vals[vals>maxval]=maxval
 }
 #vals = vals/100
 vals = log2(vals)
 vals[vals<0]=0
 data = cbind(data[,1:3],vals,data[,4:5])
 data = pv.peaksort(data)
 return(data)
}

pv.peakset_all = function(pv, addpv, minOverlap) {

   for(i in 1:length(addpv$peaks)) {
   	
   	  message(addpv$class[PV_ID,i],' ',
              addpv$class[PV_TISSUE,i],' ',
              addpv$class[PV_FACTOR,i],' ',
              addpv$class[PV_CONDITION,i],' ',
              addpv$class[PV_TREATMENT,i],' ',
              addpv$class[PV_REPLICATE,i],' ',
              addpv$class[PV_CALLER,i])
    
      pv = pv.peakset(pv,peaks=addpv$peaks[[i]],
                      sampID      = addpv$class[PV_ID,i],
                      tissue      = addpv$class[PV_TISSUE,i],
                      factor      = addpv$class[PV_FACTOR,i],
                      condition   = addpv$class[PV_CONDITION,i],
                      treatment   = addpv$class[PV_TREATMENT,i],
                      replicate   = addpv$class[PV_REPLICATE,i],
                      control     = addpv$class[PV_CONTROL,i],
                      peak.caller = addpv$class[PV_CALLER,i],
                      reads       = addpv$class[PV_READS,i],
                      consensus   = addpv$class[PV_CONSENSUS,i],
                      readBam     = addpv$class[PV_BAMREADS,i],
                      controlBam  = addpv$class[PV_BAMCONTROL,i]
                     )
   }
   
   if(minOverlap>0 && minOverlap<1) {
      minOverlap = ceiling(length(pv$peaks) * minOverlap)	
   }
   pv = dba(pv, minOverlap=minOverlap)

   return(pv)

}

pv.minOverlap = function(vec,minval){
   res = sum(vec != -1)
   if(minval==0) {
     if (length(vec) == res) {
       return(T)
     } else {
       return(F)
     }
   } else {
      if(res >= minval) {
         return(T)
      } else {
         return(F)
      }
   }	
}

pv.countOverlap = function(vec,minval= -1){
   res = sum(vec > minval)
   return(res)
}

pv.domean = function(vals){
   return(mean(vals[vals>0]))	
}

pv.do_peaks2bed = function(peaks,chrmap=NULL,fn,numCols=4) {
	    if(!is.null(chrmap)) {
	       peaks[,1] = chrmap[peaks[,1]]
	    }
	    if(numCols>0) {
	       numCols = max(3,numCols)
           mcols = min(numCols, ncol(peaks))
        } else {
           mcols = ncol(peaks)
        }
        
        if(!is.null(fn)) {
           ds = options("scipen")
           options(scipen=8)
           write.table(peaks[,1:mcols],file=fn,#sprintf("%s.bed",fn), 
                       quote=F,row.names=F,col.names=F,sep="\t")
           options(scipen=ds$scipen)
        }
   return(peaks[,1:mcols])
}

pv.catstr = function(strvec){
   unq = unique(strvec)
   if(length(unq) == 1){
      return(unq)
   }
   str = unq[1]
   for(i in 2:length(unq)){
      str = sprintf("%s-%s",str,unq[i])	
   }
   return(str)	
}

pv.namestrings = function(crec1,crec2,crec3) {
   s1 = NULL
   s2 = NULL
   s3 = NULL
   t1 = NULL
   if(missing(crec2)) {
      crec2=crec1
   }
   if(missing(crec3)) {
      crec3=crec1   
   }
   for(i in 1:length(crec1)) {
   	  if( (crec1[i]==crec2[i]) && (crec1[i]==crec3[i])) {
   	     t1 = pv.addstr(t1,crec1[i])
   	  } else {
   	     s1 = pv.addstr(s1,crec1[i])
   	     s2 = pv.addstr(s2,crec2[i])
   	     s3 = pv.addstr(s3,crec3[i])
   	  }
   }
   if(is.null(s1)) { s1 = ""}
   if(is.null(s2)) { s2 = ""}
   if(is.null(s3)) { s3 = ""}
   if(is.null(t1)) { t1 = ""}
   
   return(list(n1=s1,n2=s2,n3=s3,tstring=t1))	
}

pv.addstr = function(s1,a1) {
   if(is.null(a1)) {
      return(s1)
   } else if (a1 == "") {
      return(s1)
   }
   
   if(a1=="TRUE") {
      a1 = "T"
   }
   if(a1=="FALSE") {
      a1 = "F"
   }
      
   if(is.null(s1)){
      s1 = a1
    } else {
   	  s1 = sprintf("%s:%s",s1,a1)
   	}
   	return(s1)	
}

pv.dovectors = function(allpeaks,classes,bKeepAll=F,maxgap=0,useExternal=TRUE,useC=TRUE){

   if(!useC) {
      stop('pv.dovectors called with useC set to FALSE.')	
   }

   if(is.character(allpeaks[1,1])){
      warning('chromosome names are strings in pv.dovectors',call.=F)
   }  

   allpeaks = pv.peaksort(allpeaks)
   allpeaks_copy <- as.data.frame(allpeaks)
   allpeaks_copy$CHR <- as.integer(allpeaks_copy$CHR)
   allpeaks_copy$START <- as.integer(allpeaks_copy$START)
   allpeaks_copy$END <- as.integer(allpeaks_copy$END)
   maxGap_copy <- as.integer(maxgap)
   result <- .Call("mo_mergeOne",allpeaks_copy,bKeepAll,maxGap_copy)
   
   colnames(result)[1:3] = c("CHR","START","END")
   if(ncol(result)>3) {
      colnames(result)[4:ncol(result)] = colnames(classes)
   }
   
   return(result)
}

pv.CalledMasks = function(pv,newpv,master) {
   master = cbind(master[,1:3],1)
   spare = pv.peakset(pv,master,peak.caller='raw',scoreCol=4,bLowerScoreBetter=F)
   spare = pv.model(spare)
   res   = pv.list(spare,spare$masks$counts)
   masternum = length(spare$peaks)
   resl = NULL
   for(i in 1:nrow(res)) {
      sampl = pv.matchingSamples(res[i,1],spare$class)
      sampvec = rep(F,nrow(master))
      for(samp in sampl) {
         sampvec = sampvec | pv.whichCalled(spare,samp,masternum)
      }
      resl = pv.listadd(resl,sampvec)      
   }
   names(resl) = res$ID
   return(resl)
}

pv.matchingSamples = function(id,classes) {
   res = which(!classes[PV_CALLER,] %in% "counts" & classes[PV_ID,] %in% id)
   return(res)  
}

pv.whichCalled = function(pv,called,master,minVal=-1) {
   called = pv$vectors[,called+3] > minVal
   master = pv$vectors[,master+3] > minVal
   res = called[master]     	
   return(res)
}


## pv.pairs -- compare all pairs of peaksets, rank by % overlap
pv.pairs = function(pv,mask,bPlot=F,attributes=pv$attributes,bAllVecs=T,
                    CorMethod="pearson",bCorOnly=F,bNonZeroCors=F,minVal=0,bFixConstantVecs=T) {

   if(missing(mask)) {
      mask=rep(T,ncol(pv$class))
      peaks = pv$peaks
   } else {
      peaks = NULL
      for(i in 1:length(mask)){
         if(mask[i]) {
            peaks = pv.listadd(peaks,pv$peaks[[i]])
         }
      }
   }
   
   tmp = NULL
   if(bAllVecs==T) {
      tmp$allvectors = pv$allvectors[,c(T,T,T,mask)]
   } else {
   	  tmp$allvectors = pv$vectors[,c(T,T,T,mask)]
   }
   tmp$vectors    = pv$vectors[,c(T,T,T,mask)]
   tmp$class      = pv$class[,mask]
   tmp$peaks      = peaks
   tmp$chrmap     = pv$chrmap
   
   numSets = sum(mask)
   resm = NULL
   resl = NULL
   cvecs = NULL
   for(first in 1:(numSets-1)) {
   	  if(bCorOnly==F){
   	     #cat(".")
   	  }
      for(second in (first+1):numSets) {
      	 if(!bCorOnly) {
            res  = pv.overlap(tmp,mask = c(first,second),minVal=minVal)
            resl = pv.listadd(resl,res)
            inall = nrow(res$inAll)
            onlya = nrow(res$onlyA)
            onlyb = nrow(res$onlyB)
            prop = inall / (inall + onlya + onlyb)
         } else {
            inall = 0
            onlya = 0
            onlyb = 0
            prop  = 0
         }
         v1 = tmp$vectors[,first+3]
         v2 = tmp$vectors[,second+3]
         if(bNonZeroCors) {
           zeros = v1==0 & v2==0
           v1 = v1[!zeros]
           v2 = v2[!zeros]
         }
         if(bFixConstantVecs) {
            if(sd(v1)==0) {
               #fval = v1[1]
               #v1[1] = fval+.000001
               #v1[2] = fval-.000001
               cvecs = c(cvecs,colnames(tmp$vectors)[first+3])	
            }
            if(sd(v2)==0) {
               #fval = v2[1]
               #v2[1] = fval+.000001
               #v2[2] = fval-.000001
               cvecs = c(cvecs,colnames(tmp$vectors)[second+3])	
            }		
         }
         if(sd(v1) && sd(v2)) {
            corval = cor(v1,v2,method=CorMethod)
         } else corval = 0
         resm = rbind(resm,c(which(mask)[first],which(mask)[second],
                             onlya,onlyb,inall,corval,prop))
      }
   }
   if(bCorOnly==F) {
      #cat("\n")
   }
   cvecs = unique(cvecs)
   if(!is.null(cvecs)) {
      cvecs = sort(cvecs)
      for(cv in cvecs) {
         warning(sprintf('Scores for peakset %s are all the same -- correlations set to zero.',cv),call.=F)	
      }	
   }
   
   o = order(resm[,7],decreasing=T)
   colnames(resm) = c("A","B","onlyA","onlyB","inAll","Cor","Overlap")
   
   if(bPlot & !bCorOnly){
      warning('Plotting in pv.occupancy unsupported',call.=F)
      #for(i in 1:length(o)) {
      #	 recnum = o[i]
      #   pv.PlotContrast(pv,resl[[recnum]],resm[recnum,1],resm[recnum,2],attributes=attributes)
      #}
   }
   
   return(resm[o,])	
}

pv.overlapToLabels = function(pv,overlap,labelatts=PV_ID) {
   for(i in 1:nrow(overlap)) {
      overlap[i,1] = pv.namestrings(pv$class[labelatts,as.numeric(overlap[i,1])])$tstring
      overlap[i,2] = pv.namestrings(pv$class[labelatts,as.numeric(overlap[i,2])])$tstring
   }
   #overlap = overlap[order(overlap[,1],overlap[,2]),]
   return(data.frame(overlap))
}

pv.orderfacs = function(facvec,decreasing=F) {
   res = order(as.numeric(as.character(facvec)),decreasing=decreasing)
   return(res)	
}

pv.normalize = function(peaks,pCol,zeroVal=-1,bLog=F){
   width   = peaks[,3] - peaks[,2]
   density = peaks[,pCol]/width 
   res = density/max(density)
   if(bLog) {
      res = log2(res)
      x = res == -Inf
      res[x] = 1
      res = res - min(res)
      res[x] = zeroVal
   } else {
      res[res == 0] = zeroVal
   }
   return(res)
}

pv.colsv = c("black","red","dodgerblue","darkgreen",
             "cyan","yellow","grey50","purple3",
             "sienna","limegreen","deeppink","lightblue",
             "rosybrown1","violet","seagreen1","slategrey",
             "lavender","orange","lightgrey","olivedrab")

pv.colorv = function(classes,cols=pv.colsv){
             
   colv = rep(0,length(classes))
   uv = unique(classes)
   for(i in 1:length(uv)){
      colv[classes==uv[i]] = cols[i]
   }
   return(colv)	
} 

pv.activefun = function(x){
   if(sum(x>0)>0){
      return(TRUE)
   } else {
      return(FALSE)
   }	
}

pv.addrow = function(x,a,y){
   if(is.null(y)) return(x)
   nm = nrow(y)
   if(is.null(nm)) return(rbind(x,a))
   if(nm == 0)     return(x)
   ncl = length(a)
   return(rbind(x,matrix(a,nm,ncl,byrow=T)))
}

pv.reorderM = function(ocm,dgram) {
   lab = rownames(ocm)	
   newlab = rev(rapply(dgram,pv.dval))
   neword = match(newlab,lab)
   newocm = ocm[neword,]
   newocm = newocm[,neword]
   return(newocm)
}
pv.dval = function(dgram) {
   att = attributes(dgram)
   if(!is.null(att$leaf)) {
   	  if(!is.null(att$label)) {
         if(att$leaf==TRUE) {
            return(att$label)
         }
      }
   }
}

pv.pcmask = function(pv,numSites, mask, sites,removeComps,cor=F){

   if(missing(numSites)) numSites = nrow(pv$vectors)
   if(is.null(numSites)) numSites = nrow(pv$vectors)  
   numSites = min(numSites,nrow(pv$vectors))
  
   if(missing(sites)) sites = 1:numSites
   if(is.null(sites)) sites = 1:numSites

   if(missing(mask)) mask = rep(T,ncol(pv$class))
   
   res = NULL   
   res$class = pv$class
   pv$values = pv$vectors[sites,c(F,F,F,mask)]
   active   = apply(pv$values,1,pv.activefun)
   numSites = min(numSites,sum(active))
   
   pv$values = pv$values[active,][1:numSites,]
   
   if(!missing(removeComps)) {
      pv$values = pv.removeComp(pv$values,numRemove=removeComps)
   }

   if(nrow(pv$values) >= sum(mask)) {
      res$pc   = princomp(pv$values,cor=cor)
   }
   res$mask = mask
   
   return(res)
}

pv.venn2 = function(mrec,n1,n2,...){
#require(limma)
   
   res = NULL
   res = pv.addrow(res,c(1,0),mrec$onlyA)
   res = pv.addrow(res,c(0,1),mrec$onlyB)
   res = pv.addrow(res,c(1,1),mrec$inAll)   

   vennDiagram(res,names=c(n1,n2), circle.col=2:3,counts.col=1,...)
}

pv.venn2 = function(olaps,l1,l2,main="",sub="") {

   counts = c(nrow(olaps$onlyA),
              nrow(olaps$onlyB),
              nrow(olaps$inAll))
   names(counts) = c("A","B","A_B")		
   counts = list(counts)
   vennPlot(counts,setlabels=c(l1,l2),mysub=sub,mymain=main)
   		
}

pv.venn3 = function(m3way,n1,n2,n3,...){
#require(limma)
   
   res = NULL
   res = pv.addrow(res,c(1,0,0),m3way$onlyA)
   res = pv.addrow(res,c(0,1,0),m3way$onlyB)
   res = pv.addrow(res,c(0,0,1),m3way$onlyC)

   res = pv.addrow(res,c(1,1,0),m3way$notC)
   res = pv.addrow(res,c(1,0,1),m3way$notB)
   res = pv.addrow(res,c(0,1,1),m3way$notA)
   
   res = pv.addrow(res,c(1,1,1),m3way$inAll)   
   
   vennDiagram(res,names=c(n1,n2,n3), circle.col=2:4,counts.col=1,...)
}

pv.venn3 = function(olaps,l1,l2,l3,main="",sub="") {

   counts = c(nrow(olaps$onlyA),
              nrow(olaps$onlyB),
              nrow(olaps$onlyC),	
              nrow(olaps$notC),
              nrow(olaps$notB),
              nrow(olaps$notA),
              nrow(olaps$inAll))
   names(counts) = c("A","B","C","A_B","A_C","B_C","A_B_C")		
   counts = list(counts)
   vennPlot(counts,setlabels=c(l1,l2,l3),mysub=sub,mymain=main)
   		
}


pv.venn4 = function(olaps,l1,l2,l3,l4,main="",sub="") {

   counts = c(nrow(olaps$onlyA),
              nrow(olaps$onlyB),
              nrow(olaps$onlyC),	
              nrow(olaps$onlyD),
              nrow(olaps$AandB),
              nrow(olaps$AandC),
              nrow(olaps$AandD),
              nrow(olaps$BandC),
              nrow(olaps$BandD),
              nrow(olaps$CandD),
              nrow(olaps$notD),
              nrow(olaps$notC),
              nrow(olaps$notB),
              nrow(olaps$notA),
              nrow(olaps$inAll))
   names(counts) = c("A","B","C","D","A_B","A_C","A_D","B_C","B_D","C_D",
                     "A_B_C","A_B_D","A_C_D","B_C_D","A_B_C_D")		
   counts = list(counts)
   
   vennPlot(counts,setlabels=c(l1,l2,l3,l4),mysub=sub,mymain=main)
   		
}

