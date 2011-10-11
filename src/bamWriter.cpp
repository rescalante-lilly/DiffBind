#include "samtools/bam.h"
#include "bamWriter.h"
#include "bam.h"

bode::BamWriter::BamWriter(std::string const &filename,void *header) {
  _fd = samopen(filename.c_str(),"wb",header);
  _open = true;
}

bode::BamWriter::BamWriter(void) {
  _open = false;
  _fd = NULL;
}

bode::BamWriter::~BamWriter(void) {
  close();
}

void bode::BamWriter::close(void) {
  if (_fd != NULL) {
    samclose(_fd);
    _fd = NULL;
  }
}

void bode::BamWriter::write(Interval const &i) {
  std::string x;
  const Bam b = dynamic_cast<const Bam&>(i);
  samwrite(_fd,b._raw);
  b.seq(x);
}

bode::Writer *bode::BamWriter::open(std::string const &filename,void *header) {
  bam_header_t *bht;

  if (_fd != NULL) {
    samclose(_fd);
  }
  bht = (bam_header_t *) header;
  _fd = samopen(filename.c_str(),"wb",header);
  _open = true;
  return this;
}
