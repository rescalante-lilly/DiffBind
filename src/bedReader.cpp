#include <stdlib.h>
#include <string>
#include <fstream>
#include <zlib.h>
#include <iostream>

#include "interval.h"
#include "util.h"
#include "bedReader.h"

bode::BedReader::BedReader(std::string const &filename) {
  char *res;
  _fd = gzopen(filename.c_str(),"r");
  _buffer = new char[maxLine];

  res = gzgets(_fd,_buffer,maxLine);
  if (strncmp(_buffer,"track",5) == 0) {
    res = gzgets(_fd,_buffer,maxLine);
  }
  if (res == NULL) {
    _eof = true;
  }
  _bseq = new Interval();
}

bode::BedReader::~BedReader(void) {
  delete[] _buffer;
  delete _bseq;
  close();
}

bode::BedReader *bode::BedReader::open(std::string const &filename) {
  BedReader *br = new BedReader(filename);
  return br;
}

void bode::BedReader::close(void) {
/*  if (_fd != NULL && _fd->is_open()) { */
  if (_fd != NULL) {
    gzclose(_fd);
/*    delete _fd; */
    _fd = NULL;
  }
}

bode::Interval *bode::BedReader::next(void) {
  char *fields[12];
  int count,strand;
  bode::Interval *rv = NULL;

  if (_buffer[0] == '\0') {
    _eof = true;
    _bseq->setUnmapped();
    return rv;
  }
  bode::trimTrailing(_buffer);
  count = bode::splits(_buffer,fields,12);
  if (count < 6) {
    _bseq->update(fields[0],atoi(fields[1]),atoi(fields[2]));
    rv = _bseq;
  } else {
    if (fields[5][0] == '-') {
      strand = -1;
    } else {
      strand = 1;
    }
    _bseq->update(fields[0],atoi(fields[1]),atoi(fields[2]),strand);
    rv = _bseq;
  }
/*  if (_fd->eof()) { */
  if (gzeof(_fd)) {
    _buffer[0] = '\0';
  } else {
/*    _fd->getline(_buffer,maxLine); */
    char *rc;
    rc = gzgets(_fd,_buffer,maxLine);
    if (rc == NULL) {
      _eof = true;
    }
  }
  return rv;
}
