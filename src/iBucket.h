#ifndef __IBUCKET_H
#define __IBUCKET_H

#include <map>
#include <string>
#include <R.h>
#include <Rinternals.h>
#include "bitBucket.h"

struct pos {
  int chrom;
  int left;
  int right;
};

class IBucket {
  public:
    IBucket(int nInts,int readLength,SEXP chrom,int *left,int *right);
    ~IBucket(void);
    bool seen(std::string &chrom,int left,int right,int strand);

  private:
    struct pos *intervals;
    int icount;
    int readlen;
    BitBucket **tops;
    BitBucket **bots;
    std::map<std::string,int> chrom2num;
    int maxChrom;
    int c2n(std::string &c);
    int bsearch(int chrom,int left,int right);
    int cmp(int position,int chrom,int left,int right);
    std::string buffer;
};

#endif
