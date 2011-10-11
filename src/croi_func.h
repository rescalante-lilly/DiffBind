#ifndef __CROI_H
#define __CROI_H

#include "interval.h"
#include "intervalSetSplit.h"

class Croi {
  public:
    Croi(void);
    ~Croi(void);
    int load(const char *filename,int insertLength);
    int count(const char *chrom,int left,int right);
    int size();
  private:
    bode::IntervalSetSplit *isets;
    bode::Interval *iv;
};

#endif
