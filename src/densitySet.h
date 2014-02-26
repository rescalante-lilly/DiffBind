#ifndef __DENSITYSET_H
#define __DENSITYSET_H

#include <map>
#include <vector>
#include <string>

#include "intervalDensity.h"

namespace bode {

struct OffsetList {
  IntervalDensity *dmap;
  int position;
};

struct IntervalIndex {
  int left;
  int right;
  int index;
};

class DensitySet {
  public:
    DensitySet(int nInts,std::string *chrom,int *left,int *right);
    ~DensitySet(void);
    void add(std::string &chrom,int left,int right);
    bool summit(const int i,int &position,unsigned int &height);
    std::vector<struct IntervalIndex> getCVec(std::string &chrom);
    std::string density(const int i);

  private:
    std::vector<struct OffsetList> intervals;
    std::map<std::string,std::vector<struct IntervalIndex> > chrom2intervals;
    int nIntervals;
    int bsearch(std::vector<struct IntervalIndex> &v,int left,int right);
    int lsearch(std::vector<struct IntervalIndex> &v,int left,int right);
    bool olap(struct IntervalIndex &ii,int left,int right);
};

}

#endif
