#include <string>
#include <ostream>

#include "sequence.h"

bode::Sequence::Sequence(std::string const &n,std::string const &s) {
  _name = n;
  _seq = s;
  _null = false;
}

void bode::Sequence::write(std::ostream &out) {
  out << _name << '\t' << _seq << std::endl;
}

std::string bode::Sequence::format(void) const {
  return _name + '\t' + _seq + '\n';
}
