#ifndef __INTERVALTREESPLIT_H
#define __INTERVALTREESPLIT_H

#include "intervalNode.h"
#include "intervalTree.h"

namespace bode {

class IntervalTreeSplit: public IntervalTree {

  public:
    virtual ~IntervalTreeSplit(void) {};

  private:
    int i_countIntervals(IntervalNode *n,int left,int right);
};

}

#endif
