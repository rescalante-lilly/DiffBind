#include <unistd.h>
#include <R.h>
#include <Rinternals.h>

#include "nodeGroup.h"
#include "croi_func.h"
#include "iBucket.h"
#include "densitySet.h"

extern "C" {

SEXP croi_count_reads(SEXP filename_r,SEXP insertLength_r,SEXP filetype_r,
                      SEXP bufferSize_r,SEXP minMapQual_r,SEXP chrom_r,
                      SEXP left_r,SEXP right_r,
                      SEXP intervalCount_r,SEXP withoutDupes_r,
                      SEXP wantSummits_r,SEXP counts_r,SEXP summits_r,
                      SEXP heights_r) {
    Croi tree;
    bode::NodeGroup *ng;
    const char *filename;
    const char *chrom;
    int insertLength, filetype, bufferSize, intervalCount, withoutDupes;
    int withSummits,minMapQual;
    int *left,*right,*counts,*summits,*heights;
    int i,readCount,loadedReadCount;
    SEXP rv;
    IBucket *intervals;
    bode::DensitySet *densities;
    
 
    intervals = NULL;
    densities = NULL;
    filename = CHAR(STRING_ELT(filename_r,0));
    insertLength = INTEGER(insertLength_r)[0];
    filetype = INTEGER(filetype_r)[0];
    bufferSize = INTEGER(bufferSize_r)[0];
    minMapQual = INTEGER(minMapQual_r)[0];
    intervalCount = INTEGER(intervalCount_r)[0];
    withoutDupes = LOGICAL(withoutDupes_r)[0];
    withSummits = LOGICAL(wantSummits_r)[0];
    ng = new bode::NodeGroup(bufferSize);
    
    left = INTEGER(left_r);
    right = INTEGER(right_r);
    counts = INTEGER(counts_r);
    summits = INTEGER(summits_r);
    heights = INTEGER(heights_r);

    readCount = 0;
    loadedReadCount = 0;
    tree.open(filename,insertLength,filetype);
    if (withoutDupes) {
      intervals = new IBucket(intervalCount,tree.getIlength(),chrom_r,left,right);
    }
    if (withSummits) {
      std::string *cnames;
      cnames = new std::string[intervalCount];
      for (int i=0;i<intervalCount;i++) {
        cnames[i].assign(CHAR(STRING_ELT(chrom_r,i)));
      }
      densities = new bode::DensitySet(intervalCount,cnames,left,right);
      delete[] cnames;
    }
      
    loadedReadCount = tree.load(bufferSize,ng,intervals,densities,minMapQual);
    readCount = loadedReadCount;
    for (i=0;i<intervalCount;i++) {
      chrom = CHAR(STRING_ELT(chrom_r,i));
      counts[i] = tree.count(chrom,left[i],right[i],withoutDupes);
    }
    ng->clear();

    while (loadedReadCount == bufferSize) {
      tree.clearCounts();
      loadedReadCount = tree.load(bufferSize,ng,intervals,densities,minMapQual);
      readCount += loadedReadCount;
      for (i=0;i<intervalCount;i++) {
        chrom = CHAR(STRING_ELT(chrom_r,i));
        counts[i] += tree.count(chrom,left[i],right[i],withoutDupes);
      }
      ng->clear();
    }

    tree.close();

    /* get summits */
    if (withSummits) {
      for (i=0;i<intervalCount;i++) {
        int p;
        unsigned int h;
        densities->summit(i,p,h);
        summits[i] = p;
        heights[i] = h;
      }
    }

    /* clean up */
    if (intervals != NULL) {
      delete intervals;
    }
    if (densities != NULL) {
      delete densities;
    }
    delete ng;
    PROTECT(rv = allocVector(INTSXP,1));
    INTEGER(rv)[0] = readCount;
    UNPROTECT(1);
    return rv;
  }

} // extern C
