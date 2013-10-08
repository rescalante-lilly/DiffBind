#ifndef __CROI_H
#define __CROI_H

#include "reader.h"
#include "nodeGroup.h"
#include "interval.h"
#include "intervalSet.h"

class Croi {
  public:
    Croi(void);
    ~Croi(void);
    void open(const char *filename,int insertLength,int filetype);
    int load(int maxReads,bode::NodeGroup *ng);
    void close(void);
    int count(const char *chrom,int left,int right,int withoutDupes);
    int size(void);
    void clearCounts(void);
  private:
    bode::IntervalSet *isets;
    bode::Interval *iv;
    bode::Reader *rdr;
    int iLength;
};

#endif
