#ifndef __CROI_H
#define __CROI_H

#include "reader.h"
#include "nodeGroup.h"
#include "interval.h"
#include "intervalSet.h"
#include "iBucket.h"
#include "densitySet.h"

class Croi {
  public:
    Croi(void);
    ~Croi(void);
    void open(const char *filename,int insertLength,int filetype);
    int load(int maxReads,bode::NodeGroup *ng,IBucket *intervals,bode::DensitySet *densities);
    void close(void);
    int count(const char *chrom,int left,int right,int withoutDupes);
    int size(void);
    void clearCounts(void);
    int getIlength(void);
  private:
    bode::IntervalSet *isets;
    bode::Interval *iv;
    bode::Reader *rdr;
    int iLength;
    int getReadLength(const char *filename,int ftype);
};

#endif
