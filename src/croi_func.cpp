#include <R_ext/Utils.h>
#include "interval.h"
#include "nodeGroup.h"
#include "reader.h"
#include "croi_func.h"

Croi::Croi(void) {
  isets = new bode::IntervalSet();
  iv = new bode::Interval();
}

void Croi::open(const char *filename,int insertLength,int ftype) {
  rdr = bode::Reader::open(filename,ftype);
  iLength = insertLength;
}

int Croi::load(int maxReads,bode::NodeGroup *ng) {
  int read_count;
  bode::Interval *read_iv;

  read_count = 0;
  while (read_count < maxReads && (read_iv = rdr->next())) {
    if (read_iv->isMapped()) {
      if (iLength > 0) {
        read_iv->extend(iLength);
      }
      isets->insert(read_iv,ng);
      read_count++;
    }
    if (read_count % 10000 == 0) {
      R_CheckUserInterrupt();
    }
//    ng->clear();
  }
  return read_count;
}

void Croi::close(void) {
  rdr->close();
  delete rdr;
}

Croi::~Croi(void) {
  delete iv;
  delete isets;
}

void Croi::clearCounts(void) {
  isets->clear();
}

int Croi::count(const char *chrom,int left,int right,int withoutDupes) {
  iv->update(chrom,left,right);
  return isets->overlapping(iv,withoutDupes);
}

int Croi::size(void) {
  return isets->count();
}
