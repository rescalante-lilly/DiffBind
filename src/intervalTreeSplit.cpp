#include <algorithm>

#include "intervalTreeSplit.h"

int bode::IntervalTreeSplit::i_countIntervals(bode::IntervalNode *n,int left,int right) {
  int count = 0;

  if (n == NULL) {
    return 0;
  }
  if (left < n->r()) {
    count += i_countIntervals(n->getLeft(),left,right);
  }
  if (right >= n->l()) {
    count += i_countIntervals(n->getRight(),left,right);
  }
  int overlap = std::min(right,n->r()) - std::max(left,n->l());
/*  int nodeLen = std::min(n->r() - n->l(), right - left);
  if ((float)overlap >= (float)nodeLen / 2.0) {
    count += n->getCount();
  }
*/
  if (overlap == 0 && n->r() == n->l()) {
    if (left <= n->l() && right > n->l()) {
      overlap = 1;
    }
  }
  if (overlap > 0) {
    count += n->getCount();
  }
  return count;
}

