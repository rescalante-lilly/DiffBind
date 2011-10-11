#ifndef __INTERVALSETSPLIT_H
#define __INTERVALSETSPLIT_H

#include "interval.h"
#include "intervalSet.h"

namespace bode {

class IntervalSetSplit: public IntervalSet {
  public:
    void insert(Interval const *inter);
};

}

#endif
