#include <string>
#include <sstream>

#include "interval.h"
#include "bed.h"

bode::Bed::Bed(std::string const &chr,int l,int r):Interval(chr,l,r) {
  _name = "";
  _strand = '.';
  _score = 0;
}

bode::Bed::Bed(std::string const &chr,int l,int r,std::string n,int sc,char s):Interval(chr,l,r) {
  _name = n;
  _strand = s;
  _score = sc;
}

bode::Bed::Bed(void) {
  _mapped = false;
}

bode::Bed::Bed(bode::Bed const &b) {
  _chrom = b._chrom;
  _left = b._left;
  _right = b._right;
  _mapped = b._mapped;
  _name = b._name;
  _strand = b._strand;
  _score = b._score;
}

void bode::Bed::update(std::string const &chr,int l,int r,std::string const &n,int score,char s) {
  _chrom = chr;
  _left = l;
  _right = r;
  _mapped = true;
  _name = n;
  if (s == '1') {
    s = '+';
  }
  _strand = s;
  _score = score;
}

void bode::Bed::update(std::string const &chr,int l,int r) {
  _chrom = chr;
  _left = l;
  _right = r;
  _mapped = true;
  _name = "";
  _strand = '.';
  _score = 0;
}

std::string bode::Bed::format(void) const {
  std::ostringstream out;
  out << _chrom << '\t' << _left << '\t' << _right;
  if (_name.length() > 0) {
    out << '\t' << _name << '\t' << _score << '\t' << _strand;
  }
  out << std::endl;
  return out.str();
}

void bode::Bed::extend(int insertLen) {
  if (_strand != '-') { // forward
    _right = _left + insertLen;
  } else {
    _left = _right - insertLen;
  }
}
