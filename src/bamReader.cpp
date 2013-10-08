#include <stdio.h>
#include <string>
#include <zlib.h>
#include <samtools/sam.h>
#include <samtools/bam.h>
#include <R.h>

#include "interval.h"
#include "bamReader.h"

bool bode::BamReader::isBam(std::string const &filename) {
  gzFile fd;
  char buf[5];
  int rd;

  fd = gzopen(filename.c_str(),"r");
  rd = gzread(fd,buf,4);
  gzclose(fd);
  return (buf[0] == 'B' && buf[1] == 'A' && buf[2] == 'M' && buf[3] == '\1');
}


bode::BamReader::BamReader(std::string const &filename) {
  if (!isBam(filename)) {
    error("file '%s' does not appear to be a BAM file (bad magic number)",filename.c_str());
  }
  _fd = samopen(filename.c_str(),"rb",0);
  _seq = bam_init1();
  _bseq = new Interval();
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
  int samrv,left,right,strand;
  std::string chrom;
  Interval *rv = NULL;

  samrv = samread(_fd,_seq);
  if (samrv > 0) {
    if (_seq->core.tid == -1) {
      _bseq->setUnmapped();
    } else {
      left = _seq->core.pos;
      right = bam_calend(&(_seq->core),bam1_cigar(_seq));
      chrom = _fd->header->target_name[_seq->core.tid];
      strand = bam1_strand(_seq) == 0 ? 1 : -1;
      _bseq->update(chrom,left,right,strand);
    }
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
