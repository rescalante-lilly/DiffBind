pv.save = function(DBAobject,file='model',dir='Robjects',pre='pv_',ext='RData',
                   compress=TRUE,compression_level=9,ascii=FALSE) {
   fn = sprintf('%s/%s%s.%s',dir,pre,file,ext)
   save(DBAobject,file=fn,compress=compress,
        compression_level=compression_level,ascii=ascii)
   return(fn)
}

pv.load = function(file='model',dir='Robjects',pre='pv_',ext='RData') {
   	DBAobject = NULL
        pv = NULL
   	load(sprintf('%s/%s%s.%s',dir,pre,file,ext))
   	if(is.null(DBAobject)) {
   	   DBAobject = pv
   	}
   	return(DBAobject)
}




## pv.writePeakset --- write out vectorized peaks as a bed file for external 
pv.writePeakset = function(pv,fname,peaks,numCols=4){
   
   if(missing(peaks)) {
      peaks = rep(T,nrow(pv$vectors))
   } else {
      if(class(peaks)=='logical') {
         peaks = which(peaks)[1]	
      }	
   }
   
   if(missing(fname)) {
      fname = NULL
   }
   
   if((class(peaks)=='numeric') || (class(peaks)=='integer')) {
      bed = pv.do_peaks2bed(pv$peaks[[peaks]],pv$chrmap,fname,numCols=numCols)
      return(bed)
   }           
   
   if(!is.null(dim(peaks))) {
      if(class(peaks[1,1])=="character") {
         bed = pv.do_peaks2bed(peaks,NULL,fname,numCols=numCols)
      } else {
         bed = pv.do_peaks2bed(peaks,pv$chrmap,fname,numCols=numCols)
      }
   } else {
      bed = pv.do_peaks2bed(pv$vectors,pv$chrmap,fname,numCols=ncol(pv$vectors))
   }

   return(bed)
}
