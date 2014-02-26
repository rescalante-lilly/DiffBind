#include <map>
#include <string>
#include <iostream>

#include "intervalDensity.h"
#include "densitySet.h"

bool compareII(struct bode::IntervalIndex a,struct bode::IntervalIndex b) {
  return (a.left < b.left || (a.left == b.left && a.right < b.right));
}

bode::DensitySet::DensitySet(int nInts,std::string *chrom,int *left,int *right) {
  nIntervals = nInts;
  intervals.reserve(nIntervals);
  for (int i=0;i<nIntervals;i++) {
    struct OffsetList ol;
    ol.dmap = new bode::IntervalDensity(right[i] - left[i]);
    ol.position = left[i];
    intervals.push_back(ol);

    struct IntervalIndex ii;
    ii.left = left[i];
    ii.right = right[i];
    ii.index = i;
    chrom2intervals[chrom[i]].push_back(ii);
  }

  std::map<std::string,std::vector<struct IntervalIndex> >::iterator it;
  for (it=chrom2intervals.begin();it!=chrom2intervals.end();it++) {
    sort(it->second.begin(),it->second.end(),compareII);
  }
}

bode::DensitySet::~DensitySet(void) {
  for (int i=0;i<nIntervals;i++) {
    delete intervals[i].dmap;
  }
}

bool bode::DensitySet::olap(struct IntervalIndex &ii,int left,int right) {
  return (std::min(right,ii.right) - std::max(left,ii.left)) > 0;
}

int bode::DensitySet::bsearch(std::vector<struct IntervalIndex> &v,int left,int right) {
  int top,bot,mid;

  bot = 0;
  top = v.size() - 1;
  while (bot < top) {
    mid = (bot+top)/2;
    if (v[mid].right < left) {
      bot = mid+1;
    } else {
      top = mid;
    }
  }
  if (top == bot && olap(v[bot],left,right)) {
    return bot;
  } else {
    return -1;
  }
}

int bode::DensitySet::lsearch(std::vector<struct IntervalIndex> &v,int left,int right) {
  int pos,top;
  top = v.size();
  for (pos=0;pos<top;pos++) {
    if (olap(v[pos],left,right)) {
      return pos;
    }
  }
  return -1;
}

void bode::DensitySet::add(std::string &chrom,int left,int right) {
  int i;
  int top;

  if (chrom2intervals.count(chrom) == 0) {
    return;
  }

  std::vector<struct IntervalIndex>& vii = chrom2intervals[chrom];
  top = vii.size();

  i = bsearch(vii,left,right);
//  i = lsearch(vii,left,right);
  while (i != -1 && i < top && olap(vii[i],left,right)) {
    intervals[vii[i].index].dmap->set(left - vii[i].left,right - vii[i].left);
    i++;
  }
}

bool bode::DensitySet::summit(const int i,int &position,unsigned int &height) {
  int p;
  unsigned int h;
  intervals[i].dmap->summit(p,h);
  position = p + intervals[i].position;
  height = h;
  return true;
}

std::vector<struct bode::IntervalIndex> bode::DensitySet::getCVec(std::string &chrom) {
  return chrom2intervals[chrom];
}

std::string bode::DensitySet::density(const int i) {
  return intervals[i].dmap->str();
}
