
pv.contrast = function(pv,group1,group2=!group1,name1="group1",name2="group2",
                       minMembers=3,categories,bMulti=T,block) {
   
   numStart = length(pv$contrasts)
    
   if(missing(group1)) {    
   	  
      if( (sum(pv.mask(pv,PV_CALLER,'source'))==0) &
          (sum(pv.mask(pv,PV_CALLER,'counts'))==0) ) {
         warning('Model must include count data for contrasts.')
         return(pv)
      }
      
      pv$contrasts = NULL # clear existing contrasts
      
   	  if(missing(categories)) {
   	     attributes = c(PV_TISSUE,PV_FACTOR,PV_CONDITION)
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
         warning("No contrasts added. Perhaps try with lower value for minMembers?")	
      }
      if(!missing(block)){
         for(i in 1:length(res)){
            if(!pv.checkBlock(res[[i]])){
               warning("Blocking factor has unmatched sample(s).")
               #res[[i]]$blocklist = NULL   		
            }
         }      	
      }
   } else {
      res = pv.addContrast(pv,group1,group2,name1,name2)
      if(!missing(block)) {
         res$contrasts[[length(res$contrasts)]]$blocklist = pv.BlockList(pv,block)
         if(!pv.checkBlock(res$contrasts[[length(res$contrasts)]])) {
            warning("Blocking factor has unmatched sample(s).")
            #res$contrasts[[length(res$contrasts)]]$blocklist = NULL   	
         }
      }
      if(!is.null(res$contrasts)) {
         res$contrasts = pv.contrastDups(res$contrasts)	
      }
      if(length(res$contrasts) == numStart) {
         warning('Unable to add redundant contrast.')  	
      }
      return(res)
   }

   pv$contrasts = pv.listaddto(pv$contrasts,res)
   
   return(pv)

}

pv.getContrasts = function(pv,minMembers=3,attributes=c(PV_TISSUE,PV_FACTOR,PV_CONDITION),block){

   srcidx =  pv.mask(pv,PV_CALLER,"source") | pv.mask(pv,PV_CALLER,"counts")
   mdata = pv$class[,srcidx]
   
   if(!missing(block)) {
   	  if(block != PV_REPLICATE) {
   	     warning('Unsupported blocking attribute')
   	  }
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

pv.contrastPairs = function(pv,minMembers=3,attributes=c(PV_TISSUE,PV_FACTOR,PV_CONDITION),block=NULL,conlist=NULL) {
   
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

   vals    = sort(unique(pv$class[attribute,]))
   attname = rownames(pv$class)[attribute]
   res = NULL
   for(val in vals) {
   	  newmask = pv.mask(pv,attribute,val)
      res = pv.listadd(res,list(attribute=attname,label=val,samples=newmask))   
   }
   return(res)	
}

pv.checkBlock = function(contrast) {
   
   if(sum(contrast$group1)!=sum(contrast$group2)) {
      return(FALSE)	
   }
   
   for(att in contrast$blocklist) {
      if(sum(contrast$group1 & att$samples) != sum(contrast$group2 & att$samples)) {
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
   edger   = F
   deseq   = F
   edgerlm = F
   for(crec in clist) {
      newrec = c(crec$name1,sum(crec$group1),crec$name2,sum(crec$group2))
      bvals = 0
      if(!is.null(crec$blocklist)) {
         for(brec in crec$blocklist) {
         	bvals = bvals+1
            newrec = c(newrec,brec$label,sum(brec$samples&(crec$group1|crec$group2)))
         }
      }
      if(!is.null(crec$edgeR)) {
         edger = T
      }
      if(!is.null(crec$edgeR$block)) {
         edgerlm = T
      }
      if(!is.null(crec$DESeq)) {
         deseq = T
      }
      
      res = rbind(res,newrec)
   }
   cvec = c("Group1","Members1","Group2","Members2")
   if(bvals > 0) {
     for(i in 1:bvals) {
        cvec = c(cvec,sprintf("Block%s",i),sprintf("InBlock%i",i))
     }
   }
   eres = NULL
   if(edger) {
      for(crec in clist) {
         if(!is.null(names(crec$edgeR))){
            if(bUsePval) {
               #eres = c(eres,sum(crec$edgeR$db$fdr$table[,EDGER_COL_PVAL]<=th,na.rm=T))
               eres = c(eres,sum(topTags(crec$edgeR$db,
                                 nrow(crec$edgeR$db$counts))$table[,EDGER_COL_PVAL]<=th,na.rm=T))
            } else {
               #eres = c(eres,sum(crec$edgeR$db$fdr$table[,EDGER_COL_FDR]<=th,na.rm=T))
               eres = c(eres,sum(topTags(crec$edgeR$db,
                             nrow(crec$edgeR$db$counts))$table[,EDGER_COL_FDR]<=th,na.rm=T))
            }
         } else {
            eres = c(eres,"?")   
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
               #eres = c(eres,sum(crec$edgeR$block$fdr$table[,EDGER_COL_PVAL]<=th,na.rm=T))
               eres = c(eres,sum(topTags(crec$edgeR$block$LRT,
                                 nrow(crec$edgeR$block$counts))$table[,EDGER_COL_PVAL+1]<=th,na.rm=T))
            } else {
               #eres = c(eres,sum(crec$edgeR$block$fdr$table[,EDGER_COL_FDR]<=th,na.rm=T))
               eres = c(eres,sum(topTags(crec$edgeR$block$LRT,
                                 nrow(crec$edgeR$block$counts))$table[,EDGER_COL_FDR+1]<=th,na.rm=T))
            }
         } else {
            eres = c(eres,"?")   
         }
      }
      res = cbind(res,eres)
      cvec = c(cvec,'DB edgeR-block')
   }
   
   eres = NULL
   if(deseq) {
      for(crec in clist) {
         if(!is.null(names(crec$DESeq)) && (class(crec$DESeq) != "try-error") ){
            if(bUsePval) {
               eres = c(eres,sum(crec$DESeq$de$pval<=th,na.rm=T))
            } else {
               eres = c(eres,sum(crec$DESeq$de$padj<=th,na.rm=T))
            }
         } else {
            eres = c(eres,NA)   
         }
      }
      res = cbind(res,eres)
      cvec = c(cvec,'DB DESeq')
   }
      
   colnames(res) = cvec
   rownames(res) = 1:nrow(res)
   return(data.frame(res))
}
