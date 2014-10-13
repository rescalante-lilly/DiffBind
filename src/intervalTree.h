#ifndef __INTERVALTREE_H
#define __INTERVALTREE_H

#include "intervalNode.h"
#include "nodeGroup.h"

namespace bode {

class IntervalTree {

  public:
    IntervalTree(void);
//    virtual ~IntervalTree(void);
    void insert(int left,int right,int strand,bode::NodeGroup *ng);
    int coverage(int point);
    int reads(int left,int right,int withoutDupes);
    int summit(int left,int right);
    int getCount(void);
    int realCount(void);
   
  private:
    IntervalNode *root;
    int count;

    IntervalNode *raw_insert(IntervalNode *node);
    void rebalance(IntervalNode *node);
    void leftRotate(IntervalNode *node);
    void rightRotate(IntervalNode *node);
    int i_coverage(IntervalNode *n,int point);
    virtual int i_countIntervals(IntervalNode *n,int left,int right,int withoutDupes);
    int i_realCount(IntervalNode *n);
};

}

#endif
