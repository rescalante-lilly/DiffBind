#ifndef __CROI_H
#define __CROI_H

#include "interval.h"
#include "intervalSet.h"

class Croi {
  public:
    Croi(void);
    ~Croi(void);
    int load(const char *filename,int insertLength);
    int count(const char *chrom,int left,int right,int withoutDupes);
    int size();
  private:
    bode::IntervalSet *isets;
    bode::Interval *iv;
};

#endif
