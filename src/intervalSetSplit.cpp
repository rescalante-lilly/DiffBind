#include <string>

#include "interval.h"
#include "intervalTreeSplit.h"
#include "intervalSetSplit.h"

void bode::IntervalSetSplit::insert(Interval const *inter) {
  const std::string chrom = inter->chrom();
  if (chroms->count(chrom) == 0) {
    (*chroms)[chrom] = new IntervalTreeSplit();
  }
  (*chroms)[chrom]->insert(inter->left(),inter->right());
}
