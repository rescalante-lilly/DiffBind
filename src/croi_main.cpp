#include <R.h>
#include <Rinternals.h>

#include "croi_func.h"

extern "C" {

/*
  static Croi *reads = NULL;

  void croi_load_reads(char **filename,int *aligned_count) {
    if (reads != NULL) {
      delete reads;
    }
    reads = new Croi();
    *aligned_count = reads->load(*filename);
  }

  void croi_count_reads(char **chrom,int *left,int *right,int *ilen,int *counts) {
    int i;

    for (i=0;i<*ilen;i++) {
      counts[i] = reads->count(chrom[i],left[i],right[i]);
    }
  }
*/

  void croi_free_tree(SEXP tree_r) {
    Croi *tree;
    tree = (Croi *) R_ExternalPtrAddr(tree_r);
    delete tree;
    R_ClearExternalPtr(tree_r);
  }

  SEXP croi_load_reads(SEXP filename_r,SEXP insertLength_r) {
    Croi *tree;
    SEXP tree_r;
    const char *filename;
    int insertLength;
 
    insertLength = INTEGER(insertLength_r)[0];
    filename = CHAR(STRING_ELT(filename_r,0));
    tree = new Croi();
    tree->load(filename,insertLength);
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

  SEXP croi_count_reads(SEXP tree_r,SEXP chrom_r,SEXP left_r,SEXP right_r,SEXP ilen_r) {
    Croi *tree;
    const char *chrom;
    int ilen,*left,*right,i,*count;
    SEXP count_r;

    tree = (Croi *) R_ExternalPtrAddr(tree_r);
    ilen = INTEGER(ilen_r)[0];
    left = INTEGER(left_r);
    right = INTEGER(right_r);
    count_r = allocVector(INTSXP,ilen);
    count = INTEGER(count_r);

    for (i=0;i<ilen;i++) {
      chrom = CHAR(STRING_ELT(chrom_r,i));
      count[i] = tree->count(chrom,left[i],right[i]);
    }
    return count_r;
  }

} // extern C
