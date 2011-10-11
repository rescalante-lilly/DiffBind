#include <string.h>
#include "samtools/bam.h"

#include "sequence.h"
#include "interval.h"
#include "bam.h"

int bode::Bam::nucleotideMap[] = {'X','A','C','X','G','X','X','X','T',
                                  'X','X','X','X','X','X','X','X','N'};

bode::Bam::Bam(bam1_t *raw,bam_header_t *hdr) {
  _raw = raw;
  _hdr = hdr;
  setInterval();
  _name = std::string(bam1_qname(_raw));
  _seq = std::string();
  _null = false;
  seq(_seq);
}

bode::Bam::Bam(void) {
  _raw = NULL;
  _hdr = NULL;
  _mapped = false;
  _name = std::string();
  _seq = std::string();
  _null = true;
}

bode::Bam::Bam(bode::Bam const &b) {
  _raw = b._raw;
  _hdr = b._hdr;
  setInterval();
  _name = std::string(bam1_qname(_raw));
  _seq = std::string();
  _seq = b._seq;
  _null = b._null;
}

bode::Bam::~Bam(void) {
}

bode::Bam& bode::Bam::operator=(bode::Bam const &b) {
  _raw = b._raw;
  _hdr = b._hdr;
  setInterval();
  _seq = b._seq;
  _name = b._name;
  _null = b._null;
  return *this;
}

bool bode::operator==(bode::Bam const &l,bode::Bam const &r) {
  return l._raw == r._raw;
}

bool bode::operator<(bode::Bam const &l,bode::Bam const &r) {
  bool rv = false;

  if (l.isMapped() && r.isMapped()) {
    rv = ((Interval) l) < ((Interval) r);
  } else if (r.isMapped()) {
    rv = true; // mapped is greater than unmapped
  }
  return rv;
}

int bode::Bam::seq(std::string &dest) const {
  uint8_t *ch;
  int i;
  ch = bam1_seq(_raw);
  dest.clear();
  dest.reserve(_raw->core.l_qseq+1);
  for (i=0;i<_raw->core.l_qseq;i++) {
    char nuc = nucleotideMap[bam1_seqi(ch,i)];
    dest += nuc;
  }
  return 1;
}

void bode::Bam::update(bam1_t *raw) {
  _raw = raw;
  if (_name.length() == 0) {
    _name = std::string(bam1_qname(_raw));
  } else {
    _name.assign(bam1_qname(_raw));
  }
  _seq = std::string();
  _null = false;
  setInterval();
}

void bode::Bam::update(std::string const &chr,int l,int r) {
  int chrInd;

  _raw->core.pos = (int32_t) l;
  chrInd = chromIndex(chr);
  if (chrInd != -1) {
    _raw->core.tid = chrInd;
  }
  setInterval();
}

void bode::Bam::setUnmapped(void) {
  _mapped = false;
  if (_raw != NULL) {
    _raw->core.flag = _raw->core.flag|BAM_FUNMAP;
  }
}

int cigar2length(char *cigar) {
  char *p,*c,*end;
  int slen = 0;

  p = c = cigar;
  end = cigar + strlen(cigar);
  while (c != end) {
    while (isdigit(*c)) {
      c++;
    }
    if (c == end) {
      /* illegal CIGAR string... somebody should complain! */
      break;    }
    int n = atoi(p);
    switch (*c) {
      case 'M':
      case 'D': slen += n; break;
      case 'I': slen -= n; break;
    }
    c++;
    p = c;
  }
  return slen;
}

void bode::Bam::setInterval(void) {
  if (_raw != NULL) {
    _mapped = !(_raw->core.flag&BAM_FUNMAP);
    if (_mapped) {
      _left = _raw->core.pos;
/*      _right = _raw->core.pos + _raw->core.l_qseq; */
      _right = bam_calend(&(_raw->core),bam1_cigar(_raw));
      _chrom = _hdr->target_name[_raw->core.tid];
    }
  } else {
    _mapped = false;
  }
}

int32_t bode::Bam::chromIndex(std::string const &chr) {
  int32_t ind = -1;
  int i;

  for (i=0;i<_hdr->n_targets;i++) {
    if (chr.compare(_hdr->target_name[i]) == 0) {
      ind = i;
      break;
    }
  }
  return ind;
}

void bode::Bam::extend(int insertLen) {
  if (bam1_strand(_raw) == 0) { // forward
    _right = _left + insertLen;
  } else {
    _left = _right - insertLen;
  }
}
