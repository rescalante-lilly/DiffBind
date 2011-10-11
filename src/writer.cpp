#include <fstream>
#include <string>

#include "interval.h"
#include "writer.h"

bode::Writer::Writer(std::string const &filename,void *header) {
  _fd = new std::ofstream(filename.c_str());
  _open = true;
}

bode::Writer::Writer(void) {
  _open = false;
  _fd = NULL;
}

bode::Writer::~Writer(void) {
  close();
}

void bode::Writer::close(void) {
  if (_fd != NULL) {
    _fd->close();
    delete _fd;
    _fd = NULL;
  }
}

void bode::Writer::write(Interval const &i) {
  (*_fd) << i.format();
}

bode::Writer *bode::Writer::open(std::string const &filename,void *header) {
  _fd = new std::ofstream(filename.c_str());
  _open = true;
  return this;
}
