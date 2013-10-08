#ifndef __INTERVALSET_H
#define __INTERVALSET_H

#include <map>
#include <string>
#include "interval.h"
#include "intervalTree.h"

namespace bode {

struct chromComp {
  bool operator()(const std::string &lhs,const std::string &rhs) const {
    int an = 0,bn = 0;
    const char *lc,*rc;
    lc = lhs.c_str();
    rc = rhs.c_str();
    an = atoi(lc+3);
    bn = atoi(rc+3);
    if (an > 0 && bn > 0) {
      return an < bn;
    } else {
      return lhs < rhs;
    }
  }
};

class IntervalSet {
  public:
    IntervalSet();
    ~IntervalSet();
    virtual void insert(Interval const *inter,bode::NodeGroup* ng);
    void clear(void);
    int coverage(std::string const &chrom,int point) const;
    int overlapping(Interval const *inter,int withoutDupes) const;
    int count() const;
    int realCount() const;
    int chromCount() const;
    std::map<std::string,IntervalTree *,chromComp>::iterator chromIter();

  protected:
    std::map<std::string,IntervalTree *,chromComp> *chroms;
};

}

#endif
