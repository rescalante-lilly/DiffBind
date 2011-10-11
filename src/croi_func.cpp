#include <R_ext/Utils.h>
#include "interval.h"
#include "reader.h"
#include "croi_func.h"

Croi::Croi(void) {
  isets = new bode::IntervalSetSplit();
  iv = new bode::Interval();
}

int Croi::load(const char *filename,int insertLength) {
  int read_count;
  std::string fn(filename);
  bode::Interval *read_iv;
  bode::Reader *rdr = bode::Reader::open(fn);

  read_count = 0;
  while ((read_iv = rdr->next())) {
    if (read_iv->isMapped()) {
      if (insertLength > 0) {
        read_iv->extend(insertLength);
      }
      isets->insert(read_iv);
      read_count++;
    }
    if (read_count % 10000 == 0) {
      R_CheckUserInterrupt();
    }
  }
  rdr->close();

  delete rdr;
  return read_count;
}

Croi::~Croi(void) {
  delete iv;
  delete isets;
}

int Croi::count(const char *chrom,int left,int right) {
  iv->update(chrom,left,right);
  return isets->overlapping(iv);
}

int Croi::size(void) {
  return isets->count();
}
