#include <R.h>
#include <Rinternals.h>

#include "nodeGroup.h"
#include "croi_func.h"

extern "C" {

SEXP croi_count_reads(SEXP filename_r,SEXP insertLength_r,SEXP filetype_r,
                      SEXP bufferSize_r,SEXP chrom_r,SEXP left_r,SEXP right_r,
                      SEXP intervalCount_r,SEXP withoutDupes_r,SEXP counts_r) {
    Croi tree;
    bode::NodeGroup *ng;
    const char *filename;
    const char *chrom;
    int insertLength, filetype, bufferSize, intervalCount, withoutDupes;
    int *left,*right,*counts;
    int i,readCount,loadedReadCount;
    SEXP rv;
 
    filename = CHAR(STRING_ELT(filename_r,0));
    insertLength = INTEGER(insertLength_r)[0];
    filetype = INTEGER(filetype_r)[0];
    bufferSize = INTEGER(bufferSize_r)[0];
    intervalCount = INTEGER(intervalCount_r)[0];
    withoutDupes = LOGICAL(withoutDupes_r)[0];
    ng = new bode::NodeGroup(bufferSize);
    
    left = INTEGER(left_r);
    right = INTEGER(right_r);
    counts = INTEGER(counts_r);

    readCount = 0;
    loadedReadCount = 0;
    tree.open(filename,insertLength,filetype);
    loadedReadCount = tree.load(bufferSize,ng);
    readCount = loadedReadCount;
    for (i=0;i<intervalCount;i++) {
      chrom = CHAR(STRING_ELT(chrom_r,i));
      counts[i] = tree.count(chrom,left[i],right[i],withoutDupes);
    }
    ng->clear();

    while (loadedReadCount == bufferSize) {
      tree.clearCounts();
      loadedReadCount = tree.load(bufferSize,ng);
      readCount += loadedReadCount;
      for (i=0;i<intervalCount;i++) {
        chrom = CHAR(STRING_ELT(chrom_r,i));
        counts[i] += tree.count(chrom,left[i],right[i],withoutDupes);
      }
      ng->clear();
    }

    PROTECT(rv = allocVector(INTSXP,1));
    INTEGER(rv)[0] = readCount;
    UNPROTECT(1);
    return rv;
  }
                        
/*
  void croi_free_tree(SEXP tree_r) {
    Croi *tree;
    tree = (Croi *) R_ExternalPtrAddr(tree_r);
    delete tree;
    R_ClearExternalPtr(tree_r);
  }
*/

/*
  SEXP croi_load_reads(SEXP filename_r,SEXP insertLength_r,SEXP filetype_r) {
    Croi *tree;
    SEXP tree_r;
    const char *filename;
    int insertLength,filetype;
 
    insertLength = INTEGER(insertLength_r)[0];
    filetype = INTEGER(filetype_r)[0];
    filename = CHAR(STRING_ELT(filename_r,0));
    tree = new Croi();
    tree->load(filename,insertLength,filetype);
    tree_r = R_MakeExternalPtr((void *)tree,R_NilValue,R_NilValue);
    R_RegisterCFinalizerEx(tree_r,croi_free_tree,TRUE);
    return tree_r;
  }

  SEXP croi_tree_size(SEXP tree) {
    Croi *reads;
    SEXP rv;
    reads = (Croi *) R_ExternalPtrAddr(tree);
    PROTECT(rv = allocVector(INTSXP,1));
    INTEGER(rv)[0] = reads->size();
    UNPROTECT(1);
    return rv;
  }

  SEXP croi_count_reads(SEXP tree_r,SEXP chrom_r,SEXP left_r,SEXP right_r,SEXP ilen_r,SEXP withoutDupes_r) {
    Croi *tree;
    const char *chrom;
    int ilen,*left,*right,i,*count,withoutDupes;
    SEXP count_r;

    tree = (Croi *) R_ExternalPtrAddr(tree_r);
    ilen = INTEGER(ilen_r)[0];
    left = INTEGER(left_r);
    right = INTEGER(right_r);
    count_r = allocVector(INTSXP,ilen);
    count = INTEGER(count_r);
    withoutDupes = LOGICAL(withoutDupes_r)[0];

    for (i=0;i<ilen;i++) {
      chrom = CHAR(STRING_ELT(chrom_r,i));
      count[i] = tree->count(chrom,left[i],right[i],withoutDupes);
    }
    return count_r;
  }
*/

} // extern C
