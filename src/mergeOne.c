#include <R.h>
#include <Rinternals.h>
#include <R_ext/RS.h>
#include <R_ext/Utils.h>

#define max(a,b) ((a)>(b)?(a):(b))
#define min(a,b) ((a)>(b)?(b):(a))

enum BedColumns { CHROM=0, LEFT=1, RIGHT=2 };

static int getIndex(SEXP list,char *name) {
  int i,llen;
  SEXP names = getAttrib(list,R_NamesSymbol);
  llen = length(list);
  for (i=0;i<llen;i++) {
    if (strcmp(CHAR(STRING_ELT(names,i)),name)==0) {
      return i;
    }
  }
  return -1;
}

SEXP mo_makeEmpty(int rows,int cols,SEXP colNames) {
  SEXP dest,chroms,lefts,rights,col;
  int i;

  // allocate the storage
  PROTECT(dest = allocVector(VECSXP,cols));
  PROTECT(chroms = allocVector(INTSXP,rows));
  PROTECT(lefts = allocVector(INTSXP,rows));
  PROTECT(rights = allocVector(INTSXP,rows));
  SET_VECTOR_ELT(dest,CHROM,chroms);
  SET_VECTOR_ELT(dest,LEFT,lefts);
  SET_VECTOR_ELT(dest,RIGHT,rights);
  for (i=RIGHT+1;i<cols;i++) {
    PROTECT(col = allocVector(REALSXP,rows));
    SET_VECTOR_ELT(dest,i,col);
  }
  UNPROTECT(cols);

  // make it a data frame
  SEXP class_attr;
  PROTECT(class_attr = allocVector(STRSXP,1));
  SET_STRING_ELT(class_attr,0,mkChar((char *)"data.frame"));
  setAttrib(dest,R_ClassSymbol,class_attr);
  UNPROTECT(1);

  // set up row names
  SEXP rnattr;
  PROTECT(rnattr = allocVector(INTSXP,rows));
  int *rnattrp = INTEGER(rnattr);
  for (i=0;i<rows;i++) {
    rnattrp[i] = i+1;
  }
  setAttrib(dest,R_RowNamesSymbol,rnattr);
  UNPROTECT(1);

  // install col names
  setAttrib(dest,R_NamesSymbol,colNames);

  return dest;
}

int mo_merge(SEXP dest,SEXP src,int keepAll,int minOverlap) {
  int rows,cols,si,di,i,wasMerged;
  int *dChrom,*sChrom;
  int *dLeft,*dRight,*sLeft,*sRight;
  double **srcScores,**destScores;
  cols = length(src);
  rows = length(VECTOR_ELT(src,0));
  di = 0;

  dChrom = INTEGER(VECTOR_ELT(dest,CHROM));
  dLeft = INTEGER(VECTOR_ELT(dest,LEFT));
  dRight = INTEGER(VECTOR_ELT(dest,RIGHT));
  sChrom = INTEGER(VECTOR_ELT(src,CHROM));
  sLeft = INTEGER(VECTOR_ELT(src,LEFT));
  sRight = INTEGER(VECTOR_ELT(src,RIGHT));

  srcScores = (double **) R_alloc(cols,sizeof(double *));
  destScores = (double **) R_alloc(cols,sizeof(double *));
  for (i=3;i<cols;i++) {
    srcScores[i] = REAL(VECTOR_ELT(src,i));
    destScores[i] = REAL(VECTOR_ELT(dest,i));
  }

  // set zero'th entry in dest
  /*SET_STRING_ELT(dChrom,0,STRING_ELT(sChrom,0)); */
  dChrom[0] = sChrom[0];
  dLeft[0] = sLeft[0];
  dRight[0] = sRight[0];
  for (i=RIGHT+1;i<cols;i++) {
    /* REAL(VECTOR_ELT(dest,i))[0] = REAL(VECTOR_ELT(src,i))[0]; */
    destScores[i][0] = srcScores[i][0];
  }

  wasMerged = 0;
  for (si=1;si<rows;si++) {
    // decide if we want to merge
    int wantMerge = (dChrom[di] == sChrom[si])
                     && (dRight[di] - sLeft[si] >= minOverlap);
    if (wantMerge) {
      dRight[di] = max(dRight[di],sRight[si]);
      for (i=RIGHT+1;i<cols;i++) {
        destScores[i][di] = max(srcScores[i][si],destScores[i][di]);
/*        REAL(VECTOR_ELT(dest,i))[di] = max(REAL(VECTOR_ELT(src,i))[si],
                                           REAL(VECTOR_ELT(dest,i))[di]);
*/
      }
      wasMerged = 1;
    } else {
      // copy src to dest
      if (wasMerged || keepAll) {
        di++;
      }
      /* SET_STRING_ELT(dChrom,di,STRING_ELT(sChrom,si)); */
      dChrom[di] = sChrom[si];
      dLeft[di] = sLeft[si];
      dRight[di] = sRight[si];
      for (i=RIGHT+1;i<cols;i++) {
        /* REAL(VECTOR_ELT(dest,i))[di] = REAL(VECTOR_ELT(src,i))[si]; */
        destScores[i][di] = srcScores[i][si];
      }
      wasMerged = 0;
    }
  }
  if (keepAll == 0 && wasMerged == 0) {
    di--;
  }
  return di+1;
}

SEXP mo_truncate(SEXP src,int len) {
  SEXP dest;
  int i,j,cols;
  double **srcScores,**destScores;

  cols = length(src);
  dest = mo_makeEmpty(len,cols,getAttrib(src,R_NamesSymbol));

  int *sChrom = INTEGER(VECTOR_ELT(src,CHROM));
  int *sLeft = INTEGER(VECTOR_ELT(src,LEFT));
  int *sRight = INTEGER(VECTOR_ELT(src,RIGHT));
  int *dChrom = INTEGER(VECTOR_ELT(dest,CHROM));
  int *dLeft = INTEGER(VECTOR_ELT(dest,LEFT));
  int *dRight = INTEGER(VECTOR_ELT(dest,RIGHT));

  srcScores = (double **) R_alloc(cols,sizeof(double *));
  destScores = (double **) R_alloc(cols,sizeof(double *));
  for (i=3;i<cols;i++) {
    srcScores[i] = REAL(VECTOR_ELT(src,i));
    destScores[i] = REAL(VECTOR_ELT(dest,i));
  }

  // copy the data
  for (i=0;i<len;i++) {
    /*SET_STRING_ELT(dChrom,i,STRING_ELT(sChrom,i)); */
    dChrom[i] = sChrom[i];
    dLeft[i] = sLeft[i];
    dRight[i] = sRight[i];
    for (j=RIGHT+1;j<cols;j++) {
      destScores[j][i] = srcScores[j][i];
/*      REAL(VECTOR_ELT(dest,j))[i] = REAL(VECTOR_ELT(src,j))[i];*/
    }
  }
  return dest;
}

int mo_validate(SEXP src) {
  int okay;
  int i,src_len;
  int chromInd,leftInd,rightInd;

  okay = 1;
  if (!isVectorList(src)) {
    error("Expecting a VectorList");
  }
  chromInd = getIndex(src,(char *)"CHR");
  leftInd = getIndex(src,(char *)"START");
  rightInd = getIndex(src,(char *)"END");
  if (chromInd != 0 || leftInd != 1 || rightInd != 2) {
    error("Expecting colnames 'chrom','left','right' in pos 1,2,3");
  }
  if (!isNumeric(VECTOR_ELT(src,0))) {
    error("Chrom column (1) should be numeric");
  }
  src_len = length(src);
  for (i=1;i<src_len;i++) {
    if (!isNumeric(VECTOR_ELT(src,i))) {
      error("Columns 2..n should be numeric");
    }
  }
  return okay;
}

SEXP mo_mergeOne(SEXP src,SEXP keepAll,SEXP minOverlap) {
  SEXP dest;
  int cols,rows,nnew;
  int cKeepAll,cMinOverlap;
  
  cKeepAll = INTEGER(keepAll)[0];
  cMinOverlap = INTEGER(minOverlap)[0];
  mo_validate(src); // throws an error if src fails to validate.
  cols = length(src);
  rows = length(VECTOR_ELT(src,0));
  R_CheckUserInterrupt();
  dest = mo_makeEmpty(rows,cols,getAttrib(src,R_NamesSymbol));
  R_CheckUserInterrupt();
  nnew = mo_merge(dest,src,cKeepAll,cMinOverlap);
  R_CheckUserInterrupt();

  // truncate the vectors
  dest = mo_truncate(dest,nnew);
  R_CheckUserInterrupt();

  UNPROTECT(2);
  return dest;
}

typedef struct ipset {
  int *chr;
  int *left;
  int *right;
  double **scores;
  int curr;
  int rows;
  int sWidth;
} *ipsetp;

int mo_isSorted(ipsetp s) {
  int i,ok;
  ok = 1;
  for (i=0;i<s->rows-1;i++) {
    if (s->chr[i] > s->chr[i+1]) {
      ok = 0;
    } else if (s->chr[i] == s->chr[i+1]) {
      if (s->left[i] > s->left[i+1]) {
        ok = 0;
      } else if (s->left[i] == s->left[i+1] && s->right[i] > s->right[i+1]) {
        ok = 0;
      }
    }
  }
  return ok;
}

int mo_cmp(ipsetp alpha,ipsetp bravo) {
  int ret;
  if (alpha->chr[alpha->curr] < bravo->chr[bravo->curr]) {
    ret = -1;
  } else if (alpha->chr[alpha->curr] > bravo->chr[bravo->curr]) {
    ret = 1;
  } else if (alpha->left[alpha->curr] < bravo->left[bravo->curr]) {
    ret = -1;
  } else if (alpha->left[alpha->curr] > bravo->left[bravo->curr]) {
    ret = 1;
  } else if (alpha->right[alpha->curr] < bravo->right[bravo->curr]) {
    ret = -1;
  } else if (alpha->right[alpha->curr] > bravo->right[bravo->curr]) {
    ret = 1;
  } else {
    ret = 0;
  }
  return ret;
}

int mo_overlap(ipsetp alpha,ipsetp bravo,int minOver) {
  int ret;
  if (alpha->chr[alpha->curr] != bravo->chr[bravo->curr]) {
    ret = 0;
  } else {
    int ocalc = min(alpha->right[alpha->curr],bravo->right[bravo->curr])
                - max(alpha->left[alpha->curr],bravo->left[bravo->curr]);
    ret = ocalc >= minOver;
  }
  return ret;
}

void mo_mergeInto(ipsetp dest,ipsetp src) {
  int i;

  dest->right[dest->curr] = max(dest->right[dest->curr],src->right[src->curr]);
  for (i=0;i<src->sWidth;i++) {
    dest->scores[i][dest->curr] = max(dest->scores[i][dest->curr],
                                      src->scores[i][src->curr]);
  }
  src->curr++;
}

void mo_initRow(ipsetp dest,ipsetp src,ipsetp zdest,double zero) {
  int i;

  dest->chr[dest->curr] = src->chr[src->curr];
  dest->left[dest->curr] = src->left[src->curr];
  dest->right[dest->curr] = src->right[src->curr];

  for (i=0;i<src->sWidth;i++) {
    dest->scores[i][dest->curr] = src->scores[i][src->curr];
  }
  for (i=0;i<zdest->sWidth;i++) {
    zdest->scores[i][zdest->curr] = zero;
  }
  src->curr++;
}

void mo_mergeOrIncrement(ipsetp dest,ipsetp src,int over,ipsetp z,double zero) {
  if (mo_overlap(dest,src,over)) {
    mo_mergeInto(dest,src);
  } else {
    dest->curr++;
    z->curr++;
    mo_initRow(dest,src,z,zero);
  }
}

ipsetp sexp2ipsetp(SEXP src) {
  ipsetp dest;
  int i;

  dest = (ipsetp) Calloc(1,struct ipset);
  dest->rows = length(VECTOR_ELT(src,0));
  dest->chr = INTEGER(VECTOR_ELT(src,0));
  dest->left = INTEGER(VECTOR_ELT(src,1));
  dest->right = INTEGER(VECTOR_ELT(src,2));
  dest->sWidth = length(src) - 3;
  dest->scores = (double **) Calloc(dest->sWidth,double *);
  for (i=0;i<dest->sWidth;i++) {
    dest->scores[i] = REAL(VECTOR_ELT(src,i+3));
  }
  dest->curr = 0;
  return dest;
}

void free_ipsetp(ipsetp *item) {
  Free((*item)->scores);
  Free((*item));
  *item = NULL;
}

void mungeTarget(ipsetp t,int l,int r) {
  for (int i=l;i<l+r;i++) {
    t->scores[i-l] = t->scores[i];
  }
  t->sWidth = r;
}

SEXP mo_mergeTwo(SEXP aexp,SEXP bexp,SEXP keep_s,SEXP overlap_s,SEXP zero_s) {
  int i,j,rowsT,colsT;
  SEXP target_s,t_names,src_names;
  ipsetp alpha,bravo,targetA,targetB;
  int keepAll,minOverlap;
  double zero;

  keepAll = INTEGER(keep_s)[0];
  minOverlap = INTEGER(overlap_s)[0];
  zero = REAL(zero_s)[0];

  /* the intervals for the left operand */
  alpha = sexp2ipsetp(aexp);
  bravo = sexp2ipsetp(bexp);

  /* sanity check: are the interval sets sorted? */
  if (!mo_isSorted(alpha) || !mo_isSorted(bravo)) {
    Rf_error("Attempt to merge unsorted interval sets.  Rejected.");
  }

  /* create target with appropriate size and column names*/
  rowsT = alpha->rows + bravo->rows;
  colsT = 3 + alpha->sWidth + bravo->sWidth;
  PROTECT(t_names = allocVector(STRSXP,colsT));
  src_names = getAttrib(aexp,R_NamesSymbol);
  for (i=0;i<length(src_names);i++) {
    SET_STRING_ELT(t_names,i,STRING_ELT(src_names,i));
  }
  j = length(src_names);
  src_names = getAttrib(bexp,R_NamesSymbol);
  for (i=3;i<length(src_names);i++) {
    SET_STRING_ELT(t_names,j,STRING_ELT(src_names,i));
    j++;
  }
  target_s = mo_makeEmpty(rowsT,colsT,t_names);
  UNPROTECT(1);
  
  /* set up pointers into the target */
  /* munge the results to create two parallel targets, except with different
     scoring regions. */
  targetA = sexp2ipsetp(target_s);
  targetB = sexp2ipsetp(target_s);
  mungeTarget(targetA,0,alpha->sWidth);
  mungeTarget(targetB,alpha->sWidth,bravo->sWidth);
  
  /* initialize target with smaller of two elements */
  if (mo_cmp(alpha,bravo) < 0) {
    /* initialize from left */
    mo_initRow(targetA,alpha,targetB,zero);
  } else {
    /* initialize from right */
    mo_initRow(targetB,bravo,targetA,zero);
  }
  
  /* iterate, merging or incrementing from left or right */
  while (alpha->curr < alpha->rows && bravo->curr < bravo->rows) {
    if (mo_cmp(alpha,bravo) < 0) {
      mo_mergeOrIncrement(targetA,alpha,minOverlap,targetB,zero);
    } else {
      mo_mergeOrIncrement(targetB,bravo,minOverlap,targetA,zero);
    }
  }
  
  /* finish off remaining of either left or right */
  while (alpha->curr < alpha->rows) {
    mo_mergeOrIncrement(targetA,alpha,minOverlap,targetB,zero);
  }
  while (bravo->curr < bravo->rows) {
    mo_mergeOrIncrement(targetB,bravo,minOverlap,targetA,zero);
  }

  /* truncate target to fit */
  target_s = mo_truncate(target_s,targetA->curr+1);

  free_ipsetp(&alpha);
  free_ipsetp(&bravo);
  free_ipsetp(&targetA);
  free_ipsetp(&targetB);

  UNPROTECT(2);
  return target_s;
}
