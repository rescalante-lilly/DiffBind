#include <string>
#include <samtools/sam.h>

#include "bam.h"
#include "interval.h"
#include "bamReader.h"

bode::BamReader::BamReader(std::string const &filename) {
  _fd = samopen(filename.c_str(),"rb",0);
  _seq = bam_init1();
  _bseq = new Bam();
  _bseq->setHeader(_fd->header);
  _eof = false;
}

bode::BamReader::~BamReader(void) {
  close();
  if (_seq != NULL) {
    if (_seq->data != NULL) {
      free(_seq->data);
      _seq->data = NULL;
    }
    free(_seq);
    _seq = NULL;
  }
  delete _bseq;
}

void bode::BamReader::close(void) {
  if (_fd != NULL) {
    samclose(_fd);
    _fd = NULL;
  }
}

bode::Interval *bode::BamReader::next(void) {
  int samrv;
  Interval *rv = NULL;

  samrv = samread(_fd,_seq);
  if (samrv > 0) {
    _bseq->update(_seq);
    rv = _bseq;
  } else {
    _bseq->setUnmapped();
    _eof = true;
  }
  return rv;
}

bode::BamReader *bode::BamReader::open(std::string const &filename) {
  BamReader *br = new BamReader(filename);
  return br;
}
