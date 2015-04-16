#include <map>
#include <string>
#include "interval.h"
#include "intervalSet.h"
using namespace std;

bode::IntervalSet::IntervalSet() {
  chroms = new map<string,IntervalTree *,chromComp>();
}

bode::IntervalSet::~IntervalSet() {
  map<string,IntervalTree *,chromComp>::iterator it;
  for (it=chroms->begin();it!=chroms->end();it++) {
    delete it->second;
  }

  delete chroms;
}

void bode::IntervalSet::insert(Interval const *inter,bode::NodeGroup* ng) {
  const string chrom = inter->chrom();
  if (chroms->count(chrom) == 0) {
    (*chroms)[chrom] = new IntervalTree();
  }
  (*chroms)[chrom]->insert(inter->left(),inter->right(),inter->strand(),ng);
}

int bode::IntervalSet::coverage(std::string const &chrom,int point) const {
  int rv;
  if (chroms->count(chrom) == 0) {
    rv = 0;
  } else {
    rv = (*chroms)[chrom]->coverage(point);
  }
  return rv;
}

int bode::IntervalSet::overlapping(Interval const *inter,int withoutDupes) const {
  int rv;
  const string chrom = inter->chrom();
  if (chroms->count(chrom) == 0) {
    rv = 0;
  } else {
    rv = (*chroms)[chrom]->reads(inter->left(),inter->right(),withoutDupes);
  }
  return rv;
}

int bode::IntervalSet::chromCount() const {
  return chroms->size();
}

int bode::IntervalSet::count() const {
  int reads = 0;
  map<string,IntervalTree *,chromComp>::iterator it;
  for (it=chroms->begin();it!=chroms->end();it++) {
    reads += it->second->getCount();
  }
  return reads;
}

int bode::IntervalSet::realCount() const {
  int reads = 0;
  map<string,IntervalTree *,chromComp>::iterator it;
  for (it=chroms->begin();it!=chroms->end();it++) {
    reads += it->second->realCount();
  }
  return reads;
}

std::map<std::string,bode::IntervalTree *,bode::chromComp>::iterator bode::IntervalSet::chromIter() {
  return chroms->begin();
}

void bode::IntervalSet::clear(void) {
  chroms->clear();
}
