#include <string>
#include <sstream>
#include "interval.h"

#define MAX_MAPQUAL 255

bode::Interval::Interval(std::string const &chr,int l,int r) {
  _chrom = chr;
  _left = l;
  _right = r;
  _mapped = true;
  _strand = 1;
  _mapqual = MAX_MAPQUAL;
}

bode::Interval::Interval(std::string const &chr,int l,int r,int s) {
  _chrom = chr;
  _left = l;
  _right = r;
  _mapped = true;
  _strand = s;
  _mapqual = MAX_MAPQUAL;
}

bode::Interval::Interval(std::string const &chr,int l,int r,int s,int m) {
  _chrom = chr;
  _left = l;
  _right = r;
  _mapped = true;
  _strand = s;
  _mapqual = m;
}

bode::Interval::Interval(Interval const &i) {
  _chrom = i._chrom;
  _left = i._left;
  _right = i._right;
  _mapped = i._mapped;
  _strand = i._strand;
  _mapqual = i._mapqual;
}

bode::Interval &bode::Interval::operator=(Interval const &i) {
  if (this != &i) {
    _chrom = i._chrom;
    _left = i._left;
    _right = i._right;
    _mapped = i._mapped;
    _strand = i._strand;
    _mapqual = i._mapqual;
  }
  return *this;
}

bool bode::operator==(bode::Interval const &l,bode::Interval const &r) {
  return (l._chrom == r._chrom &&
          l._left == r._left &&
          l._right == r._right &&
          l._strand == r._strand);
}

bool bode::operator<(Interval const &l,Interval const &r) {
  return (l._chrom < r._chrom ||
          (l._chrom == r._chrom && l._left < r._left) ||
          (l._chrom == r._chrom && l._left == r._left && l._right < r._right) ||
          (l._chrom == r._chrom && l._left == r._left && l._right == r._right && l._strand > r._strand));
}

void bode::Interval::updatecstr(const char *c_str,int l,int r) {
  _chrom.assign(c_str);
  _left = l;
  _right = r;
  _mapped = true;
  _mapqual = MAX_MAPQUAL;
}

void bode::Interval::update(std::string const &chr,int l,int r) {
  _chrom = chr;
  _left = l;
  _right = r;
  _mapped = true;
  _mapqual = MAX_MAPQUAL;
}

void bode::Interval::update(std::string const &chr,int l,int r,int s) {
  _chrom = chr;
  _left = l;
  _right = r;
  _mapped = true;
  _strand = s;
}

void bode::Interval::update(std::string const &chr,int l,int r,int s,int m) {
  _chrom = chr;
  _left = l;
  _right = r;
  _mapped = true;
  _strand = s;
  _mapqual = m;
}

std::string bode::Interval::format(void) const {
  std::ostringstream out;
  out << _chrom << ":" << _left << "-" << _right;
  return out.str();
}

void bode::Interval::extend(int insertLen) {
  if (_strand > 0) {
    _right = _left + insertLen;
  } else {
    _left = _right - insertLen;
  }
}
