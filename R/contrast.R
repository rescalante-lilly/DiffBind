

pv.contrast = function(pv,group1,group2=!group1,name1="group1",name2="group2",
                       minMembers=3,categories,bMulti=T,block) {
   
   numStart = length(pv$contrasts)
    
   if(missing(group1)) {    
   	  
      if( (sum(pv.mask(pv,PV_CALLER,'source'))==0) &
          (sum(pv.mask(pv,PV_CALLER,'counts'))==0) ) {
         warning('Model must include count data for contrasts.',call.=F)
         return(pv)
      }
      
      pv$contrasts = NULL # clear existing contrasts
      
   	  if(missing(categories)) {
   	     attributes = c(PV_TISSUE,PV_FACTOR,PV_CONDITION,PV_TREATMENT)
   	  } else {
   	     attributes=categories
   	  }
   	  if(missing(block)) {
         res = pv.getContrasts(pv,minMembers=minMembers,attributes=attributes)
         if(bMulti) {
            res = pv.listaddto(res,pv.contrastPairs(pv,minMembers=minMembers,attributes=attributes,conlist=res))  
         }
      } else {
         res = pv.getContrasts(pv,minMembers=minMembers,attributes=attributes,block=block)
         if(bMulti) {
            res = pv.listaddto(res,pv.contrastPairs(pv,minMembers=minMembers,attributes=attributes,block=block))   
         }
      }
      if(!is.null(res)) {
         res = pv.contrastDups(res)
      } else {
         warning("No contrasts added. Perhaps try more categories, or lower value for minMembers.",call.=F)	
      }
      if(!missing(block)){
         problem = NULL
      	 issues = 0
         for(i in 1:length(res)){
            if(!pv.checkBlock(res[[i]],bWarning=F)){
               #warning("Blocking factor has unmatched sample(s).")
               issues = issues+1
               problem = res[[i]]
               res[[i]]$blocklist = NULL   		
            }
         }
         if(issues > 0) {
            if(issues < length(res)) {
               warning('Blocking factor not used for some contrasts:',call.=F)	
            } else {
               warning('Blocking factor invalid for all contrasts:',call.=F)
            }
            x = pv.checkBlock(problem)		
         }      	
      }
   } else {
      res = pv.addContrast(pv,group1,group2,name1,name2)
      if(!missing(block)) {
         res$contrasts[[length(res$contrasts)]]$blocklist = pv.BlockList(pv,block)
         if(!pv.checkBlock(res$contrasts[[length(res$contrasts)]])) {
            #warning("Blocking factor has unmatched sample(s).")
            res$contrasts[[length(res$contrasts)]]$blocklist = NULL   	
         }
      }
      if(!is.null(res$contrasts)) {
         res$contrasts = pv.contrastDups(res$contrasts)	
      }
      if(length(res$contrasts) == numStart) {
         warning('Unable to add redundant contrast.',call.=F)  	
      }
      return(res)
   }

   pv$contrasts = pv.listaddto(pv$contrasts,res)
   
   return(pv)

}

pv.getContrasts = function(pv,minMembers=3,attributes=c(PV_TISSUE,PV_FACTOR,PV_CONDITION,PV_TREATMENT),block){

   srcidx =  pv.mask(pv,PV_CALLER,"source") | pv.mask(pv,PV_CALLER,"counts")
   mdata = pv$class[,srcidx]
   
   if(!missing(block)) {
   	  #if(block != PV_REPLICATE) {
   	  #   warning('Unsupported blocking attribute')
   	  #}
      block = pv.BlockList(pv,block)
   } else block = NULL
   
   mcats = attributes
   jobs = NULL
   for(mcat in mcats) {
      vals = unique(mdata[mcat,])
      if (length(vals)>1) {
      	 members = NULL
         for(i in 1:length(vals)) {
            members = c(members,sum(mdata[mcat,]== vals[i]))            
         }
         for(i in 1:length(vals)) {
            if(members[i] >= minMembers) {
               if(i < length(vals)){
                  for(j in (i+1):length(vals)) {
                     if(members[j] >= minMembers) {
                        job = NULL
                        job$group1 = pv.mask(pv,mcat,vals[i]) & srcidx
                        job$group2 = pv.mask(pv,mcat,vals[j]) & srcidx
                        job$name1  = vals[i]
                        job$name2  = vals[j]
                        job$blocklist = block
                        jobs = pv.listadd(jobs,job)  
                     }
                  }      
               }
               if((sum(members)-members[i]) >= minMembers) {
                  job = NULL
                  job$group1 = pv.mask(pv,mcat,vals[i]) & srcidx
                  job$group2 = pv.mask(pv,mcat,vals[i],merge='nand') & srcidx
                  job$name1  = vals[i]
                  job$name2  = sprintf("!%s",vals[i])
                  job$blocklist = block
                  jobs = pv.listadd(jobs,job)                 
               }   
            }
         }     
      }
   }
   return(jobs)   
}

pv.contrastPairs = function(pv,minMembers=3,attributes=c(PV_TISSUE,PV_FACTOR,PV_CONDITION,PV_TREATMENT),block=NULL,conlist=NULL) {
   
   if(length(attributes)==1) {
      return(pv$contrasts)
   }
   
   clist = conlist
   srcmask = pv.mask(pv,PV_CALLER,"source") | pv.mask(pv,PV_CALLER,"counts")
   numatts = length(attributes)
   res = NULL
   for(a1 in 1:(numatts-1)) {
   	  att1 = attributes[a1]
      val1 = unique(pv$class[att1,])
      if(length(val1)>1) {
         for(i in 1:length(val1)) {
            for (a2 in (a1+1):numatts) {
         	   att2 = attributes[a2]
               val2 = unique(pv$class[att2,])
               if(length(val2)>1) {
                  for(j in 1:length(val2)) {
                     res = pv.listadd(res,list(a1=att1,v1=val1[i],a2=att2,v2=val2[j]))
                  }
               }
            }
         }
      }
   }
   
   if(!missing(block)) {
      block = pv.BlockList(pv,block)
   } else block = NULL
   
   if(!is.null(res)) {
      for(i in 1:length(res)) {
   	     m1 = pv.mask(pv,res[[i]]$a1,res[[i]]$v1,mask=srcmask,merge='and')
   	     m1 = pv.mask(pv,res[[i]]$a2,res[[i]]$v2,mask=m1,merge='and')
   	     if(sum(m1)>=minMembers) {
   	  	    if(i < length(res)) {
               for(j in (i+1):length(res)) {
   	              m2 = pv.mask(pv,res[[j]]$a1,res[[j]]$v1,mask=srcmask,merge='and')
   	              m2 = pv.mask(pv,res[[j]]$a2,res[[j]]$v2,mask=m2,merge='and') 
   	              if( (sum(m2)>=minMembers) && (sum(m1 & m2) == 0) ) {
   	                 crec = NULL
   	                 crec$group1 = m1
   	                 crec$group2 = m2
   	                 crec$name1  = sprintf("%s:%s",res[[i]]$v1,res[[i]]$v2)
   	                 crec$name2  = sprintf("%s:%s",res[[j]]$v1,res[[j]]$v2)
   	                 crec$blocklist = block
   	                 clist = pv.listadd(clist,crec)
   	              }
   	           }                      
            }
            if(sum((!m1)&srcmask) >=minMembers) {
   	           crec = NULL
   	           crec$group1 = m1
   	           crec$group2 = (!m1) & srcmask
   	           crec$name1  = sprintf("%s:%s",res[[i]]$v1,res[[i]]$v2)
   	           crec$name2  = sprintf("!%s:%s",res[[i]]$v1,res[[i]]$v2)
   	           crec$blocklist = block
   	           clist = pv.listadd(clist,crec)	
   	         }     
         }
      }	
   }
   #if(!is.null(block)){
   #   for(i in 1:length(clist)){
   #      if(!pv.checkBlock(clist[[i]])){
   #         warning("Unable to add blocking factor: unmatched samples.")
   #         clist[[i]]$blocklist = NULL   		
   #      }
   #   }      	
   #} 
   return(clist)
}


pv.addContrast = function(pv,group1,group2=!group1,name1="group1",name2="group2") {
  
  if(is.null(group1) | is.null(group2)) {
     stop('Null group, can not add contrast')	
  }
  
  if(!is.logical(group1)) {
  	 if(max(group1) > length(pv$peaks)) {
  	    stop('Invalid sample number in first group.')	
  	 }
     temp = rep(F,length(pv$peaks))
     temp[group1] = T
     group1 = temp
  }
  
  if(!is.logical(group2)) {
  	 if(max(group2) > length(pv$peaks)) {
  	    stop('Invalid sample number in second group.')	
  	 }  	
     temp = rep(F,length(pv$peaks))
     temp[group2] = T
     group2 = temp
  }
    
  if(sum(group1)==0) {
     return(pv)
  }
  if(sum(group2)==0) {
     return(pv)
  }
  
  if( length(group1) != length(pv$peaks) || length(group2) != length(pv$peaks) ) {
     stop('Length of vector specifying groups greater than number of samples.')
  }
      
  crec = NULL
  crec$name1  = name1
  crec$name2  = name2
  crec$group1 = group1
  crec$group2 = group2
  pv$contrasts = pv.listadd(pv$contrasts,crec)
  return(pv)	
}

pv.contrastDups = function(clist) {
   numc = length(clist)
   if(numc <= 1) {
      return(clist)
   }
   res = rep(T,numc)
   for(i  in 1:(numc-1)) {
      crec = clist[[i]]
      for(j in (i+1):numc) {
         rec2 = clist[[j]]
         if ( (identical(crec$group1,rec2$group1) && identical(crec$group2,rec2$group2)) ||
              (identical(crec$group1,rec2$group2) && identical(crec$group2,rec2$group1)) ){
             res[j] = F   
          }
      }      	
   }
   newc = NULL
   for(i in 1:numc) {
      if (res[i]) {
         newc = pv.listadd(newc,clist[[i]])
      }	
   }
   return(newc)	
}


pv.BlockList = function(pv,attribute=PV_REPLICATE) {

   if(is.numeric(attribute)) {
      if(length(attribute)>1) {
         vec = rep(F,length(pv$peaks))
         vec[attribute]=T
         attribute=vec	
      }
   }
   if(class(attribute)=='numeric') {   
      vals    = sort(unique(pv$class[attribute,]))
      attname = rownames(pv$class)[attribute]
      if(attname == 'Peak caller') {
         attname = 'Caller'
      }
      res = NULL
      for(val in vals) {
   	     newmask = pv.mask(pv,attribute,val)
         res = pv.listadd(res,list(attribute=attname,label=val,samples=newmask))   
      }
   } else { #logical vector(s)
      if(class(attribute)=='logical') {
      	 if(length(attribute)!=length(pv$peaks)) {
      	    stop('Length of attribute vector must equal total number of peaksets.')	
      	 }
         attribute = list(true=attribute,false=!attribute)	
      }
      if(class(attribute)!='list') {
         stop('attribute must be a DBA_ attribute, a logical vector, or a list of logical vectors.')	
      }
      attname = 'Block'
      if(is.null(names(attribute))) {
         names(attribute) = 1:length(attribute)
      }
      res = NULL
      hasatt = rep(F,length(pv$peaks))
      for (i in 1:length(attribute)) {
      	 att = attribute[[i]]
      	 if(is.numeric(att)) {
      	 	att = rep(F,length(pv$peaks))
      	 	att[attribute[[i]]]=T
      	 }
      	 hasatt = hasatt | att
         res = pv.listadd(res,list(attribute=attname,label=names(attribute)[i],samples=att))	
      }
      if(sum(!hasatt)>0) {
         res = pv.listadd(res,list(attribute=attname,label="other",samples=!hasatt))	
      }
   }
   
   return(res)	
}

pv.checkBlock = function(contrast,bCheckBalanced=F,bCheckMultiple=T,bCheckCross=T,bCheckUnique=T,bCheckAllInOne=T,bWarning=T) {
   
   if(bCheckBalanced){
      if(sum(contrast$group1)!=sum(contrast$group2)) {
         if(bWarning) warning("Blocking factor has unmatched sample(s).",call.=F)
         return(FALSE)	
      }
   
      for(att in contrast$blocklist) {
         if(sum(contrast$group1 & att$samples) != sum(contrast$group2 & att$samples)) {
            if(bWarning) warning("Blocking factor has unmatched sample(s).",call.=F)
            return(FALSE)
         }
      }
   }
   
   if(bCheckMultiple) {
      if(length(contrast$blocklist)<2) {
         if(bWarning) warning('Blocking factor has only one value',call.=F)
         return(FALSE)	
      }	
   }
   
   if(bCheckCross || bCheckAllInOne) {
   	  cross  = FALSE
   	  allone = FALSE
      for(att in contrast$blocklist) {
         if(sum(contrast$group1 & att$samples) & sum(contrast$group2 & att$samples)) {
            cross = TRUE
         }
         if(sum(att$samples & (contrast$group1 | contrast$group2)) ==
            sum(contrast$group1 | contrast$group2))
            allone=TRUE
      }
      if(!cross && bCheckCross) {
         if(bWarning) warning('No blocking values are present in both groups',call.=F)	
         return(FALSE)
      }
      if(allone && bCheckAllInOne) {
         if(bWarning) warning('A blocking value applies to all samples in both groups',call.=F)	
         return(FALSE)
      }
   }
   
   if(bCheckUnique) {
      unique = rep(0,sum(contrast$group1)+sum(contrast$group2)) 	
      for(att in contrast$blocklist) {
         unique = unique + (contrast$group1 & att$samples)
         unique = unique + (contrast$group2 & att$samples)	
      }
      if(sum(unique>1)>0) {
         if(bWarning) warning('Some sample(s) have more than one value for blocking factor',call.=F)	
         return(FALSE)
      }
   }
   
   return(TRUE)
}

EDGER_COL_PVAL = 4
EDGER_COL_FDR  = 5
pv.listContrasts = function(pv,th=0.1,bUsePval=F) {
   if(is.null(pv$contrasts)) {
      return(NULL)
      clist = pv.contrast(pv)   
   } else {
      clist = pv$contrasts
   }
   if(is.null(clist)) {
      return(NULL)
   }
   res = NULL
   edger    = F
   deseq1   = F
   deseq2   = F
   edgerlm  = F
   deseq1lm = F
   deseq2lm = F
   maxbvals = 0
   for(crec in clist) {
      newrec = c(crec$name1,sum(crec$group1),crec$name2,sum(crec$group2))
      bvals = 0
      if(!is.null(crec$blocklist)) {
         for(brec in crec$blocklist) {
         	bvals = bvals+1
            newrec = c(newrec,brec$label,sum(brec$samples&(crec$group1|crec$group2)))
         }
      }
      maxbvals = max(bvals,maxbvals)
      if(!is.null(crec$edgeR)) {
         edger = T
      }
      if(!is.null(crec$edgeR$block)) {
         edgerlm = T
      }
      if(!is.null(crec$DESeq1)) {
         deseq1 = T
      }
      if(!is.null(crec$DESeq1$block)) {
         deseq1lm = T
      }
      if(!is.null(crec$DESeq2)) {
        deseq2 = T
      }
      if(!is.null(crec$DESeq2$block)) {
        deseq2lm = T
      }
      if(!is.null(res)) {
         if(length(newrec)>ncol(res)) {
            res = cbind(res,matrix("-",nrow(res),length(newrec)-ncol(res)))	
         } else if (length(newrec)<ncol(res)) {
            newrec = c(newrec,rep("-",ncol(res)-length(newrec)))	
         }
      }
      res = rbind(res,newrec)
   }
   cvec = c("Group1","Members1","Group2","Members2")
   if(maxbvals > 0) {
     for(i in 1:maxbvals) {
        cvec = c(cvec,sprintf("Block%sVal",i),sprintf("InBlock%i",i))
     }
   }
   eres = NULL
   if(edger) {
      for(crec in clist) {
         if(!is.null(names(crec$edgeR))){
         	if(is.null(crec$edgeR$LRT)) {  
               if(bUsePval) {
                  eres = c(eres,sum(topTags(crec$edgeR$db,
                                            nrow(crec$edgeR$db$counts))$table[,EDGER_COL_PVAL]<=th,na.rm=T))
               } else {
                  eres = c(eres,sum(topTags(crec$edgeR$db,
                                            nrow(crec$edgeR$db$counts))$table[,EDGER_COL_FDR]<=th,na.rm=T))
               } 
            } else {
               if(bUsePval) {
                  eres = c(eres,sum(topTags(crec$edgeR$LRT,
                                            nrow(crec$edgeR$db$counts))$table[,EDGER_COL_PVAL+1]<=th,na.rm=T))
               } else {
                  eres = c(eres,sum(topTags(crec$edgeR$LRT,
                                            nrow(crec$edgeR$db$counts))$table[,EDGER_COL_FDR+1]<=th,na.rm=T))
               } 
            }
         } else {
            eres = c(eres,"-")   
         }
      }
      res = cbind(res,eres)
      cvec = c(cvec,'DB edgeR')
   }
   
   eres = NULL
   if(edgerlm) {
      for(crec in clist) {
         if(!is.null(names(crec$edgeR$block))){
            if(bUsePval) {
               eres = c(eres,sum(topTags(crec$edgeR$block$LRT,
                                 nrow(crec$edgeR$block$counts))$table[,EDGER_COL_PVAL+1]<=th,na.rm=T))
            } else {
               eres = c(eres,sum(topTags(crec$edgeR$block$LRT,
                                 nrow(crec$edgeR$block$counts))$table[,EDGER_COL_FDR+1]<=th,na.rm=T))
            }
         } else {
            eres = c(eres,"-")   
         }
      }
      res = cbind(res,eres)
      cvec = c(cvec,'DB edgeR-block')
   }
   
   eres = NULL
   if(deseq1) {
      for(crec in clist) {
         if(!is.null(names(crec$DESeq1)) && (class(crec$DESeq1) != "try-error") ){
            if(bUsePval) {
               eres = c(eres,sum(crec$DESeq1$de$pval<=th,na.rm=T))
            } else {
               eres = c(eres,sum(crec$DESeq1$de$padj<=th,na.rm=T))
            }
         } else {
            eres = c(eres,"-")   
         }
      }
      res = cbind(res,eres)
      cvec = c(cvec,'DB DESeq')
   }
   
   eres = NULL
   if(deseq1lm) {
      for(crec in clist) {
         if(!is.null(names(crec$DESeq1$block)) && (class(crec$DESeq1) != "try-error") ){
            if(bUsePval) {
               eres = c(eres,sum(crec$DESeq1$block$de$pval<=th,na.rm=T))
            } else {
               eres = c(eres,sum(crec$DESeq1$block$de$padj<=th,na.rm=T))
            }
         } else {
            eres = c(eres,"-")   
         }
      }
      res = cbind(res,eres)
      cvec = c(cvec,'DB DESeq-block')
   }
   
   eres = NULL
   if(deseq2) {
     for(crec in clist) {
       if(!is.null(names(crec$DESeq2)) && (class(crec$DESeq2) != "try-error") ){
         if(bUsePval) {
           eres = c(eres,sum(crec$DESeq2$de$pval<=th,na.rm=T))
         } else {
           eres = c(eres,sum(crec$DESeq2$de$padj<=th,na.rm=T))
         }
       } else {
         eres = c(eres,"-")   
       }
     }
     res = cbind(res,eres)
     cvec = c(cvec,'DB DESeq2')
   }
   
   eres = NULL
   if(deseq2lm) {
     for(crec in clist) {
       if(!is.null(names(crec$DESeq2$block)) && (class(crec$DESeq2) != "try-error") ){
         if(bUsePval) {
           eres = c(eres,sum(crec$DESeq2$block$de$pval<=th,na.rm=T))
         } else {
           eres = c(eres,sum(crec$DESeq2$block$de$padj<=th,na.rm=T))
         }
       } else {
         eres = c(eres,"-")   
       }
     }
     res = cbind(res,eres)
     cvec = c(cvec,'DB DESeq2-block')
   }  
   
   colnames(res) = cvec
   rownames(res) = 1:nrow(res)
   return(data.frame(res))
}

pv.design = function(DBA,categories=c(DBA_CONDITION,DBA_TREATMENT,DBA_TISSUE,DBA_FACTOR,DBA_REPLICATE)) {

   facs = NULL
   for(cond in categories) {
      if(length(unique(DBA$class[cond,]))>1) {
         facs = c(facs,cond)
      } else {
      	 catname = "UNKNOWN"
         if(cond == DBA_CONDITION) {
         	catname = "DBA_CONDITION"
         }	
         if(cond == DBA_TREATMENT) {
         	catname = "DBA_TREATMENT"
         }	
         if(cond == DBA_TISSUE) {
         	catname = "DBA_TISSUE"
         }	
         if(cond == DBA_FACTOR) {
         	catname = "DBA_FACTOR"
         }	
         if(cond == DBA_REPLICATE) {
         	catname = "DBA_REPLICATE"
         }
         warning(sprintf("Category %s lacks multiple unique values, ignored in design",catname),call.=F)	
      }	
   }
   
   red = DBA$class[facs,]
   if(is.null(nrow(red))) {
      red = matrix(red,1,length(red))
      rownames(red) = rownames(DBA$class)[facs]
      colnames(red) = colnames(DBA$class)
   }
   dmatrix = data.frame(t(red))
   terms = colnames(dmatrix)
   form = "~"
   for(term in terms) {
      form = sprintf("%s+%s",form,term)	
   }
   
   design = model.matrix(as.formula(form),data=dmatrix)
   
   return(design)   	
}

pv.resultsDBA = function(DBA,contrasts,methods=DBA$config$AnalysisMethod,th=0.1,bUsePval=F,fold=0,bDB=T,bNotDB=F,bUp=F,bDown=F,bAll=T) {

   if(missing(contrasts)) {
      contrasts = 1:length(DBA$contrasts)
   }
   
   res = NULL
   for(contrast in contrasts) {
      for(method in methods) {
         res = pv.doResults(res,DBA,contrast,method,th,bUsePval,fold,bDB=bDB,bNotDB=bNotDB,bUp=bUp,bDown=bDown,bAll=bAll)		
      }
   }
   
   if(is.null(res)) {
      warning('No valid contrasts/methods specified.',call.=F)	
   }
   
   return(res)	
	
}

pv.doResults = function(res,DBA,contrast,method,th,bUsePval,fold=0,bDB=T,bNotDB=F,bAll=F,bUp=F,bDown=F) {

    if(method=='edgeR' || method=="edgeRGLM") {
       if(is.null(DBA$contrasts[[contrast]]$edgeR)) {
          return(res)	
       }
       methname = "edgeR"
       block = ""
    }
    if (method=='edgeRlm') {
       if(is.null(DBA$contrasts[[contrast]]$edgeR$block)) {
          return(res)	
       }
       methname="edgeR"
       block="block"	
    }
    if(method=='DESeq1' || method=='DESeq1GLM') {
       if(is.null(DBA$contrasts[[contrast]]$DESeq1)) {
          return(res)	
       }      
       methname = "DESeq"
       block = ""
    }
    if (method=='DESeq1Block') {
       if(is.null(DBA$contrasts[[contrast]]$DESeq1$block)) {
          return(res)	
       }           
       methname="DESeq"
       block="block"	
    }
    if(method=='DESeq2' || method=='DESeq2GLM') {
      if(is.null(DBA$contrasts[[contrast]]$DESeq2)) {
        return(res)  
      }      
      methname = "DESeq2"
      block = ""
    }
    if (method=='DESeq2Block') {
      if(is.null(DBA$contrasts[[contrast]]$DESeq2$block)) {
        return(res)	
      }           
      methname="DESeq2"
      block="block"	
    }
    
    if(bNotDB) {
     useth = 1
    } else {
     useth = th
    }
    
    rep = suppressWarnings(dba.report(DBA,contrast=contrast,method=method,th=useth,bUsePval=bUsePval,fold=0,DataType=DBA_DATA_FRAME))
    
    if(is.null(rep)) {
       return(res)	
    }
    
    id = sprintf("%s_vs_%s",DBA$contrast[[contrast]]$name1,DBA$contrast[[contrast]]$name2)
    
    if(bUsePval) {
    	scores = rep$"p-value"
    } else {
    	scores = rep$FDR
    }    
    
    db = (scores <= th) & (abs(rep$Fold) >= fold)
    up = rep$Fold >= 0
    peaks = cbind(rep[,1:3],abs(rep$Fold))
           
    if(bDB) {
       if(bAll) {
          if(sum(db)) {
             res = dba.peakset(res,peaks=peaks[db,],sampID=id,factor="DB",tissue="All",condition=methname,treatment=block)	
          }
       }
       if(bUp) {
          if(sum(db&up)) {	
             res = dba.peakset(res,peaks=peaks[db&up,],sampID=id, factor ="DB", tissue ="Gain",condition=methname,treatment=block)	
          }
       }
       if(bDown) {
       	  if(sum(db&!up)) {
             res = dba.peakset(res,peaks=peaks[db&!up,],sampID=id, factor ="DB", tissue ="Loss",condition=methname,treatment=block)	
          }
       }		
    }
    if(bNotDB) {
       if(bAll) {
       	  if(sum(!db)) {
             res = dba.peakset(res,peaks=peaks[!db,],sampID=id, factor ="!DB", tissue ="All",condition=methname,treatment=block)	
          }
       }
       if(bUp) {
       	  if(sum(!db&up)) {
             res = dba.peakset(res,peaks=peaks[!db&up,],sampID=id, factor ="!DB", tissue ="Gain",condition=methname,treatment=block)	
          }
       }
       if(bDown) {
       	  if(sum(!db&!up)) {
             res = dba.peakset(res,peaks=peaks[!db&!up,],sampID=id, factor ="!DB", tissue ="Loss",condition=methname,treatment=block)	
          }
       }		
    }    
    
    return(res)	
}

