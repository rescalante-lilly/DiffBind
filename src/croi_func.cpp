#include <stdexcept>
#include <R_ext/Utils.h>
#include "interval.h"
#include "nodeGroup.h"
#include "reader.h"
#include "croi_func.h"
#include "iBucket.h"
#include "densitySet.h"

Croi::Croi(void) {
  isets = new bode::IntervalSet();
  iv = new bode::Interval();
}

int Croi::getReadLength(const char *filename,int ftype) {
  bode::Reader *fd;
  bode::Interval *iv = NULL;
  int rlen = -1;
  fd = bode::Reader::open(filename,ftype);
  iv = fd->next();
  while (iv && !(iv->isMapped())) {
    iv = fd->next();
  }
  if (iv) {
    rlen = iv->right() - iv->left();
  }
  fd->close();
  delete fd;
  return rlen;
}

void Croi::open(const char *filename,int insertLength,int ftype) {
  rdr = bode::Reader::open(filename,ftype);
  iLength = std::max(insertLength,getReadLength(filename,ftype));
}

int Croi::getIlength(void) {
  return iLength;
}

int Croi::load(int maxReads,bode::NodeGroup *ng,IBucket *intervals,bode::DensitySet *densities) {
  int read_count;
  bode::Interval *read_iv;
  std::string x(128,' ');

  read_count = 0;
  while (read_count < maxReads && (read_iv = rdr->next())) {
    if (read_iv->isMapped()) {
      if (iLength > 0) {
        read_iv->extend(iLength);
      }
      x.assign(read_iv->chrom());
      if (intervals == NULL || !intervals->seen(x,read_iv->left(),read_iv->right(),read_iv->strand())) {
        isets->insert(read_iv,ng);
        read_count++;
      }
      if (densities != NULL) {
        try {
          densities->add(x,read_iv->left(),read_iv->right());
        } catch (std::out_of_range oor) {
          warning("trapped exception from intervalDensity");
        }
      }
    }
    if (read_count % 10000 == 0) {
      R_CheckUserInterrupt();
    }
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
  iv->updatecstr(chrom,left,right);
  return isets->overlapping(iv,withoutDupes);
}

int Croi::size(void) {
  return isets->count();
}
