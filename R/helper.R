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
FDRth=100
pv.macs = function(fn){
 data = read.table(fn,blank.lines.skip=T,header=T)
 #if(data[1,1] =='chr'){
 #   data = read.table(fn,skip=18)
 #}
 #data = data[data[,9]<=FDRth,]
 res  = pv.peaksort(data)
 return(res[,c(1:3,7)]) 
}

pv.wold = function(fn){
 data = read.table(fn)[,2:5]	
 res  = pv.peaksort(data)
 return(res[,1:4]) 
}

pv.swembl = function(fn){
 data = read.table(fn,skip=14)
 res  = pv.peaksort(data)
 return(res[,1:4]) 
}

pv.readbed = function(fn,skipnum=0){
 data = read.table(fn,skip=skipnum)
 res  = pv.peaksort(data)
 return(res)
}

pv.bayes = function(fn){
 data = read.table(fn)
 idx = data[,4]>0.5
 data = data[idx,]
 res  = pv.peaksort(data)
 return(res)
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

pv.domean = function(vals){
   return(mean(vals[vals>0]))	
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

pv.dovectors = function(allpeaks,classes,bKeepAll=F,maxgap=0,useExternal=TRUE,useC=TRUE){

   if(is.character(allpeaks[1,1])){
      warning('chromosome names are strings in pv.dovectors')
   }  
   
   if(useC) { ## Use embedded C routine

      allpeaks = pv.peaksort(allpeaks)
      allpeaks_copy <- as.data.frame(allpeaks)
      allpeaks_copy$CHR <- as.integer(allpeaks_copy$CHR)
      allpeaks_copy$START <- as.integer(allpeaks_copy$START)
      allpeaks_copy$END <- as.integer(allpeaks_copy$END)
      maxGap_copy <- as.integer(maxgap)
      result <- .Call("mo_mergeOne",allpeaks_copy,bKeepAll,maxGap_copy)

   } else if(useExternal) { ## USE PERL SCRIPT

      dovectors = tempfile(as.character(Sys.getpid()))
      dovectors.bed = sprintf("%s.bed",dovectors)
      dovectors.merge = sprintf("%s_merge.bed",dovectors)
      
      allpeaks = pv.peaksort(allpeaks)
      pv.do_peaks2bed(allpeaks,fn=dovectors.bed,numCols=0)
      
      if(bKeepAll) {
   	     cmdstr = sprintf("/home/brown22/bin/pv.vectors %s %d T",dovectors,maxgap)
      } else {
   	    cmdstr = sprintf("/home/brown22/bin/pv.vectors %s %d F",dovectors,maxgap)
      }
      system(cmdstr)
      result = read.table(dovectors.merge)

      system(sprintf("rm %s %s", dovectors.bed, dovectors.merge))

     
   } else {  ## R VERSION
	
      allpeaks = pv.peaksort(allpeaks) 
      chrs=unique(allpeaks[,1])
      result=NULL
      for(i in 1:length(chrs)){
         cpeaks = allpeaks[allpeaks[,1] == chrs[i],]
         if(is.null(dim(cpeaks))) {
            cpeaks = rbind(NULL,cpeaks)
         }
         result = rbind(result,pv.dovectors.chr(cpeaks,bKeepAll=bKeepAll,maxgap=maxgap))
      }
      tmpf = tempfile(as.character(Sys.getpid()))#tmpdir='.')
      pv.do_peaks2bed(result, NULL,tmpf)
      result = read.table(tmpf)
      unlink(tmpf)
   }

   colnames(result)[1:3] = c("CHR","START","END")
   if(ncol(result)>3) {
      colnames(result)[4:ncol(result)] = colnames(classes)
   }
   
   return(result)
}

pv.dovectors.chr = function(peaks,bKeepAll=F,maxgap=0){

   #return(dovectors.chr(peaks,bKeepAll,maxgap))
  
   chr=peaks[1,1]
   nvecs = ncol(peaks)-3
   vref  = nvecs + 3
   merged=NULL
   g  = pv.gaps(peaks)
   firsttime = T
   while((sum(g<0)>0) || firsttime) {
      #cat('-')
      firsttime = F
      gr = g > maxgap
      r  = pv.regions(gr)
      
      if(bKeepAll) {
   	     #insert non-overlaps back in: 
         inmerge = c(apply(r,1,pv.listnums),recursive=T)
         oneoffs = which(!((1:nrow(peaks)) %in% (inmerge)))
         toadd = cbind(oneoffs,oneoffs)
         r = rbind(r,toadd)
         o = order(r[,1])
         r = r[o,]
      } else {
         r = r[r[,1] != r[,2],]
      }
      
      if(is.null(dim(r))){
         r = rbind(NULL,r)
      }
      if(nrow(r) == 0) {
         return(NULL)
      }

      merged = NULL
      epeaks = peaks[,3]
      merged = t(apply(r,1,pv.dovec,peaks,chr))
      if(ncol(merged)>0){
      	 if(nvecs>0) {
            colnames(merged) = c("CHR","START","END",1:nvecs)
         } else {
         	colnames(merged) = c("CHR","START","END")
         }
         merged[,2] = peaks[merged[,2],2]
         merged[,3] = peaks[merged[,3],3]
      } else {
         return(NULL)
      }
      g     = pv.gaps(merged)
      peaks = merged 
      bKeepAll = T
   }

   return(merged)
}

PV_IDX_START = 2
PV_IDX_END   = 3

pv.gaps = function(data){
   numEls = nrow(data)
   if(is.null(numEls)) return(0)
   if(numEls < 2) return(0)
   gaps   = c(0,data[2:numEls,PV_IDX_START]-data[1:(numEls-1),PV_IDX_END])
   return(gaps)
}

pv.regions = function(gapsRange){
   numGaps = length(gapsRange)
   before  = c(T,gapsRange[1:(numGaps-1)])
   after   = c(gapsRange[2:(numGaps)],T) 
   start   = (gapsRange & !after)
   if(!gapsRange[1]) start[1] = T
   end     = (!before & gapsRange)
   starti  = which(start)
   endi    = which(end)-1
   if(length(endi) < length(starti)) endi = c(endi,numGaps)
   res     = cbind(starti,endi)
   return(res)
}

pv.listnums = function(x){
   return(x[1]:x[2])	
}


pv.dovec = function(r,peaks,chr){ 
   start = r[1]
   end   = r[2]
   vecs = ncol(peaks)
   if (vecs > 3) {
      if (start != end) {
   	     newm = peaks[start:end,4:vecs]
         nrev = apply(newm,2,max)
         end  = r[1] + which.max(peaks[start:end,3]) - 1
      } else {
         nrev = peaks[r[1],4:vecs]  
      }
      res = c(chr,start,end,nrev)
   } else {
      res = c(chr,start,end)
   }

   return(res)
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
   
   if(PV_DEBUG == FALSE){
      return
   }
   
   #write(sprintf('%s\n',str),file=file,append=T)

}


pv.contrast2 = function(vectors,n1,n2,minVal=0){
   v1 = vectors[,n1+3] > minVal
   v2 = vectors[,n2+3] > minVal
   
   allpeaks = v1 | v2
   inAll   = v1 & v2
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

   return(res)   	
}

pv.contrast3 = function(vectors,n1,n2,n3,minVal=0){
	
   v1 = vectors[,n1+3] > minVal
   v2 = vectors[,n2+3] > minVal
   v3 = vectors[,n3+3] > minVal
   allpeaks = v1 | v2 | v3
   
   nv1 = sum(v1); nv2 = sum(v2); nv3 = sum(v3)
   
   inAll  =  v1 &  v2 &  v3
   onlyA  =  v1 & !v2 & !v3
   onlyB  = !v1 &  v2 & !v3 
   onlyC  = !v1 & !v2 &  v3
   
   sofar = inAll | onlyA | onlyB | onlyC
   v1 = v1 & !sofar
   v2 = v2 & !sofar
   v3 = v3 & !sofar
   
   notA   = !v1 &  (v2 | v3)
   v2 = v2 & !notA
   v3 = v3 & !notA
   notB   = !v2 &  (v1 | v3)
   v1 = v1 & !notB
   v3 = v3 & !notB
   notC   = !v3 &  (v1 | v3)       
	    
   res = list(onlyA = vectors[onlyA,c(1:3,(n1+3),(n2+3),(n3+3))],
   	          onlyB = vectors[onlyB,c(1:3,(n1+3),(n2+3),(n3+3))],
   	          onlyC = vectors[onlyC,c(1:3,(n1+3),(n2+3),(n3+3))],
   	          notA  = vectors[notA,c(1:3,(n1+3),(n2+3),(n3+3))],
   	          notB  = vectors[notB,c(1:3,(n1+3),(n2+3),(n3+3))],
   	          notC  = vectors[notC,c(1:3,(n1+3),(n2+3),(n3+3))],
   	          inAll = vectors[inAll,c(1:3,(n1+3),(n2+3),(n3+3))])      
   
   for(i in 1:7){
      if(is.null(nrow(res[[i]]))){
         res[[i]] = as.matrix(t(res[[i]]))
      }
   }	
           
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

pv.addrow = function(x,a,y){
   if(is.null(y)) return(x)
   nm = nrow(y)
   if(is.null(nm)) return(rbind(x,a))
   if(nm == 0)     return(x)
   ncl = length(a)
   return(rbind(x,matrix(a,nm,ncl,byrow=T)))
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

## pv.pairs -- compare all pairs of peaksets, rank by % overlap
pv.pairs = function(pv,mask,bPlot=T,attributes=pv$attributes,bAllVecs=T,
                    CorMethod="pearson",bCorOnly=F,bNonZeroCors=F,minVal=0) {

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
   for(first in 1:(numSets-1)) {
   	  if(bCorOnly==F){
   	     #cat(".")
   	  }
      for(second in (first+1):numSets) {
      	 if(!bCorOnly) {
            res  = pv.overlap(tmp,first,second,minVal=minVal)
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
         corval = cor(v1,v2,method=CorMethod)
         resm = rbind(resm,c(which(mask)[first],which(mask)[second],
                             onlya,onlyb,inall,corval,prop))
      }
   }
   if(bCorOnly==F) {
      #cat("\n")
   }
   
   o = order(resm[,7],decreasing=T)
   colnames(resm) = c("A","B","onlyA","onlyB","inAll","Cor","Overlap")
   
   if(bPlot & !bCorOnly){
      for(i in 1:length(o)) {
      	 recnum = o[i]
         pv.PlotContrast(pv,resl[[recnum]],resm[recnum,1],resm[recnum,2],attributes=attributes)
      }
   }
   
   return(resm[o,])	
}


pv.PlotContrast = function(pv,mrec,A,B,C=NULL,attributes=pv$attributes) {
   first  = pv$class[attributes,A]
   second = pv$class[attributes,B]
   if(!is.null(C)) { 
      third = pv$class[attributes,C]
      res = pv.namestrings(first,second,third)
      pv.venn3(mrec,res$n1,res$n2,res$n3,main=res$tstring)
   } else {
      res = pv.namestrings(first,second)
      pv.venn2(mrec,res$n1,res$n2,main=res$tstring)
   }    	
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

pv.whichPeaksets = function(pv,mask) {

   if(missing(mask)) {
      warning('mask required')
      return(NULL)
   }
   if(is.null(mask)) {
      warning('mask required')
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

   return(list(A=A,B=B,C=C))
}

## pv.triads -- compare all triads of peaksets, rank by % overlap
pv.triads = function(pv,mask) {

   if(missing(mask)) {
      mask=rep(T,ncol(pv$class))
   }
   
   tmp = NULL
   tmp$allvectors = pv$allvectors[,c(T,T,T,mask)]
   tmp$class      = pv$class[,mask]
   
   numSets = sum(mask)
   resm = NULL
   for(first in 1:(numSets-2)) {
      for(second in (first+1):(numSets-1)) {
         for(third in (second+1):numSets) {
            res  = pv.overlap(tmp,first,second,third)
            inAll = nrow(res$inAll)
            onlyA = nrow(res$onlyA)
            onlyB = nrow(res$onlyB)
            onlyC = nrow(res$onlyC)
            notA  = nrow(res$notA)
            notB  = nrow(res$notB)
            notC  = nrow(res$notC)
            A = onlyA + inAll + notB + notC
            B = onlyB + inAll + notA + notC
            C = onlyC + inAll + notA + notB
            prop = inAll / min(A,B,C)
            resm = rbind(resm,c(first,second,third,onlyA,onlyB,onlyC,notA,notB,notC,inAll,prop))
         }
      }
   }
   
   o = order(resm[,11],decreasing=T)
   colnames(resm) = c("A","B","C","onlyA","onlyB","onlyC","notA","notB","notC","inBoth","Overlap")
   return(resm[o,])	
}

pv.fixvals = function(vec,oldval,newval) {
   vec[vec == oldval] = newval
   return(vec)	
} 

pv.FixMin = function(vecs,oldval=-1,newval=0) {
   res= apply(vecs[,4:ncol(vecs)],1,pv.fixvals,oldval,newval)
   return(cbind(vecs[,1:3],t(res)))	
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

pv.dositename = function(svec,GenomeString) {
   if(missing(GenomeString)) {
      return(sprintf("%s:%s-%s",svec[1],as.numeric(svec[2]),as.numeric(svec[3])))
   } else {
      return(sprintf('=HYPERLINK("http://localhost:7085/UnibrowControl?version=%s&seqid=%s&start=%s&end=%s","%s:%s-%s")',
                      GenomeString,svec[1],as.numeric(svec[2])-250,as.numeric(svec[3])+250,
                      svec[1],svec[2],svec[3]))   
   }	
}

pv.dositename_map = function(svec,chrmap,GenomeString) {
   if(missing(GenomeString)) {
      return(sprintf("%s:%s-%s",chrmap[svec[1]],as.numeric(svec[2]),as.numeric(svec[3])))
   } else {
      return(sprintf('=HYPERLINK("http://localhost:7085/UnibrowControl?version=%s&seqid=%s&start=%s&end=%s","%s:%s-%s")',
                      GenomeString,chrmap[svec[1]],as.numeric(svec[2])-250,as.numeric(svec[3])+250,
                      chrmap[svec[1]],svec[2],svec[3]))   
   }	
}

## pv.peaks2bed -- write out peaks
pv.peaks2bed = function(pv,pnum,fn) {
   pv.do_peaks2bed(pv$peaks[[pnum]],pv$chrmap,fn)	
}

pv.overlap_util = function(g1,g2) {
  if(nrow(g1) >= nrow(g2)) {
     longer = g1
     shorter = g2
  } else {
     longer = g2
     shorter = g1
  }
#  if (pv.peaksort(longer) != longer) {
#    cat('Longer set must be peaksorted first!\n')
#    return(NULL)
#  }
  shorter = pv.peaksort(shorter)
  chrs=unique(as.character(longer[,1]))
  res=NULL
  for(i in 1:length(chrs)){
     res = c(res,pv.ol_chr(shorter[shorter[,1]==chrs[i],2:3],
                     longer[longer[,1]==chrs[i],2:3]))
  }
  return(res)
}

pv.ol_chr = function(g1,g2) {

  if(is.null(g1)){
     return(rep(F,nrow(g2)))	
  }
  if(is.null(nrow(g1))) {
     g1 = matrix(g1,1,2)
  }

  z = apply(g1,1,pv_dool,g2)
  #x = g1[,1] <= g2[,2]
  #y = g1[,2] >= g2[,1]
  #z = x & y
  
  return(z)  	
}

pv_dool = function(val,arr) {
   x =  arr[,2] > val[1]
   y =  arr[,1] < val[2]
   z = x & y
   return(which(z))	
}

 
dovectors.chr = function(peaks,bKeepAll=F,maxgap=0){
   chr     = peaks[1,1]
   start   = 0
   end     = 0
   npeaks  = nrow(peaks)
   lastcol = ncol(peaks)
   cols    = lastcol-3
   i       = 1
   res     = matrix(0,npeaks,lastcol)
   pnum    = 1
   while(i < npeaks){
   	  single = T
      start = peaks[i,2]
      end   = peaks[i,3]
      vals  = peaks[i,4:lastcol] 
      expand = T
      while(expand & (i<npeaks)) {
         if (end - peaks[i+1,2] >= maxgap) {
            single = F
            i = i+1
            end = max(end,peaks[i,3])
            # get maximum value for each column
            vals = sapply(1:cols,function(x) max(vals[x],peaks[i,x+3]))
         } else {
            expand = F
         }
       }
      if (!single || bKeepAll) {
         res[pnum,] = c(chr,start,end,vals)
         pnum = pnum + 1
      }
      i=i+1  	
      single=T 
   }
   if((i==npeaks) & single & bKeepAll){
      res[pnum,]  = peaks[i,]
      pnum=pnum+1
   }
   if(pnum ==1) return(NULL)
   return(res[1:(pnum-1),])	
}

pv.attstring = function(pv,atcode) {
	
	res = rownames(pv$class[atcode])
    return(res)
}


pv.overlapToLabels = function(pv,overlap,labelatts=PV_ID) {
   for(i in 1:nrow(overlap)) {
      overlap[i,1] = pv.namestrings(pv$class[labelatts,as.numeric(overlap[i,1])])$tstring
      overlap[i,2] = pv.namestrings(pv$class[labelatts,as.numeric(overlap[i,2])])$tstring
   }
   #overlap = overlap[order(overlap[,1],overlap[,2]),]
   return(data.frame(overlap))
}


pv.sortfacs = function(facvec,decreasing=F) {
   res = sort(as.numeric(as.character(facvec)),decreasing=decreasing)
   return(res)	
}

pv.orderfacs = function(facvec,decreasing=F) {
   res = order(as.numeric(as.character(facvec)),decreasing=decreasing)
   return(res)	
}


pv.CalledMasks = function(pv,newpv,master) {
   master = cbind(master[,1:3],1)
   spare = pv.peakset(pv,master,peak.caller='raw')
   spare = pv.model(spare)
   res = pv.list(spare,spare$masks$source)
   masternum = length(spare$peaks)
   resl = NULL
   for(i in 1:nrow(res)) {
      sampl = pv.matchingSamples(res[i,],spare$class)
      sampvec = rep(F,nrow(master))
      for(samp in sampl) {
         sampvec = sampvec | pv.whichCalled(spare,samp,masternum)
      }
      resl = pv.listadd(resl,sampvec)      
   }
   names(resl) = res$ID
   return(resl)
}

pv.matchingSamples = function(sprops,classes) {
   res = which(!classes[PV_CALLER,] %in% "source" & classes[PV_ID,] %in% sprops$ID)
   return(res)  
}

pv.whichCalled = function(pv,called,master,minVal=-1) {
   called = pv$vectors[,called+3] > minVal
   master = pv$vectors[,master+3] > minVal
   res = called[master]     	
   return(res)
}

pv.getMatching = function(one,many) {
  z = apply(many,2,function(x,y){x==y},one)
  res = which(z[1,] & z[2,])
  return(res) 
}

pv.CalledMasks = function(pv,newpv,master){
   numNew = length(newpv$peaks)
   numOld = length(pv$peaks) - numNew
   master = cbind(master[,1:3],1)
   pv = pv.peakset(pv,master,readBam="",controlBam="")
   pv = pv.vectors(pv,c(1:numOld,length(pv$peaks)),minOverlap=1)
   masternum = length(pv$peaks)
   classes = pv$class[c(PV_BAMREADS,PV_BAMCONTROL),]
   resl = NULL
   for(i in 1:numNew) {
     peaksets = pv.getMatching(newpv$class[c(PV_BAMREADS,PV_BAMCONTROL),i],classes)
     sampvec = rep(F,nrow(master))
     for(samp in peaksets) {
         sampvec = sampvec | pv.whichCalled(pv,samp,masternum)
     }
     resl = pv.listadd(resl,sampvec)
   }
   names(resl) = newpv$class[PV_ID,]
   return(resl)
}

pv.check = function(pv) {
   if(is.null(pv$vectors)) {
   	  if(is.null(pv$minOverlap)) {
   	     minOverlap=2
   	  } else {
   	     minOverlap = pv$minOverlap
   	  }
      pv = pv.vectors(pv,minOverlap=minOverlap)
   }
   return(pv)
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

pv.reorderM = function(ocm,dgram) {
   lab = rownames(ocm)	
   newlab = rev(rapply(dgram,pv.dval))
   neword = match(newlab,lab)
   newocm = ocm[neword,]
   newocm = newocm[,neword]
   return(newocm)
}
   
pv.setScore = function(pv,score,bLog=F,minMaxval) {
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
         pv = pv.vectors(pv,minOverlap=1,bAnalysis=F,bAllSame=T)
      } else {
         stop('No sites have activity greater than minMaxval')
      }
   }   
   return(pv)
}

