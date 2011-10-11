#include <string>
#include <sstream>
#include "interval.h"

bode::Interval::Interval(std::string const &chr,int l,int r) {
  _chrom = chr;
  _left = l;
  _right = r;
  _mapped = true;
}

bode::Interval::Interval(Interval const &i) {
  _chrom = i._chrom;
  _left = i._left;
  _right = i._right;
  _mapped = i._mapped;
}

bode::Interval &bode::Interval::operator=(Interval const &i) {
  if (this != &i) {
    _chrom = i._chrom;
    _left = i._left;
    _right = i._right;
    _mapped = i._mapped;
  }
  return *this;
}

bool bode::operator==(bode::Interval const &l,bode::Interval const &r) {
  return (l._chrom == r._chrom &&
          l._left == r._left &&
          l._right == r._right);
}

bool bode::operator<(Interval const &l,Interval const &r) {
  return (l._chrom < r._chrom ||
          (l._chrom == r._chrom && l._left < r._left) ||
          (l._chrom == r._chrom && l._left == r._left && l._right < r._right));
}

void bode::Interval::update(std::string const &chr,int l,int r) {
  _chrom = chr;
  _left = l;
  _right = r;
  _mapped = true;
}

std::string bode::Interval::format(void) const {
  std::ostringstream out;
  out << _chrom << ":" << _left << "-" << _right;
  return out.str();
}

int bode::Interval::score(void) const {
  return 0;
}

void bode::Interval::extend(int insertLen) {
  _right = _left + insertLen;
}
