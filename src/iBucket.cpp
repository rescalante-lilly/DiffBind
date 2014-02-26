#include <map>
#include <string>
#include <R.h>
#include <Rinternals.h>
#include "bitBucket.h"
#include "iBucket.h"
#include "string.h"

IBucket::IBucket(int nInts,int readLength,SEXP chrom,int *left,int *right) {
  maxChrom = 0;
  icount = nInts;
  intervals = new struct pos[icount];
  readlen = readLength;
  buffer.reserve(1024);
  tops = new BitBucket *[icount];
  bots = new BitBucket *[icount];
  for (int i=0;i<icount;i++) {
    buffer.assign(CHAR(STRING_ELT(chrom,i)));
    intervals[i].chrom = c2n(buffer);
    intervals[i].left = left[i] - readlen;
    intervals[i].right = right[i] + readlen;
    tops[i] = new BitBucket(right[i] - left[i] + readlen * 2);
    bots[i] = new BitBucket(right[i] - left[i] + readlen * 2);
  }
}

IBucket::~IBucket(void) {
  for (int i=0;i<icount;i++) {
    delete tops[i];
    delete bots[i];
  }
  delete intervals;
  delete tops;
  delete bots;
}

int IBucket::c2n(std::string &c) {
  int rv;
  if (chrom2num.count(c) == 0) {
    rv = maxChrom++;
    chrom2num[c] = rv;
  } else {
    rv = chrom2num[c];
  }
  return rv;
}

int IBucket::cmp(int position,int chrom,int left,int right) {
  int rv = 0;
  if (chrom < intervals[position].chrom) {
    rv = -1;
  } else if (chrom > intervals[position].chrom) {
    rv = 1;
  } else if (right <= intervals[position].left) {
    rv = -1;
  } else if (left >= intervals[position].right) {
    rv = 1;
  } // else they overlap, so cmp returns 0
  return rv;
}

int IBucket::bsearch(int chrom,int left,int right) {
  int top,bot,mid;
  bool found = false;
  int c;

  bot = 0;
  top = icount-1;
  mid = -1;
  while (found == false && top != bot && mid != bot) {
    mid = (top-bot)/2 + bot;
    c = cmp(mid,chrom,left,right);
    if (c < 0) {
      top = mid;
    } else if (c > 0) {
      bot = mid;
    } else {
      found = true;
    }
  }
  if (found == false) {
    mid = -1;
  }
  return mid;
}

bool IBucket::seen(std::string &chrom,int left,int right,int strand) {
  int cnum = c2n(chrom);
  int index = bsearch(cnum,left,right);
  if (index != -1) {
    if (strand > 0) {
      if (left >= intervals[index].left) {
        if (tops[index]->isSet(left - intervals[index].left)) {
          return true;
        } else {
          tops[index]->set(left - intervals[index].left);
        }
      }
    } else {
      if (right <= intervals[index].right) {
        if (bots[index]->isSet(intervals[index].right - right)) {
          return true;
        } else {
          bots[index]->set(intervals[index].right - right);
        }
      }
    }
  }
  return false;
}
