#ifndef __NODEGROUP_H
#define __NODEGROUP_H

#include "intervalNode.h"

namespace bode {

class NodeGroup {
  public:
    NodeGroup(int size);
    ~NodeGroup(void);
    bode::IntervalNode *get(void);
    void pop(void);
    void clear(void);

  private:
    int nsize;
    int ncount;
    bode::IntervalNode *nodes;
};

}

#endif
