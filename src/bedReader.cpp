#include <stdlib.h>
#include <string>
#include <fstream>
#include <zlib.h>
#include <iostream>

#include "interval.h"
#include "bed.h"
#include "util.h"
#include "bedReader.h"

bode::BedReader::BedReader(std::string const &filename) {
  char *res;
/*  _fd = new std::ifstream(filename.c_str()); */
  _fd = gzopen(filename.c_str(),"r");
  _buffer = new char[maxLine];

/*  _fd->getline(_buffer,maxLine); */
  res = gzgets(_fd,_buffer,maxLine);
  if (strncmp(_buffer,"track",5) == 0) {
    res = gzgets(_fd,_buffer,maxLine);
/*    _fd->getline(_buffer,maxLine); */
  }
  if (res == NULL) {
    _eof = true;
  }
  _bseq = new Bed();
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
  int count;
  bode::Interval *rv = NULL;

  if (_buffer[0] == '\0') {
    _eof = true;
    _bseq->setUnmapped();
    return rv;
  }
  count = bode::splits(_buffer,fields,12);
  if (count == 3) {
    _bseq->update(fields[0],atoi(fields[1]),atoi(fields[2]));
    rv = _bseq;
  } else {
    if (fields[5][0] == '1') {
      fields[5][0] = '+';
    }
    _bseq->update(fields[0],atoi(fields[1]),atoi(fields[2]),fields[3],atoi(fields[4]),fields[5][0]);
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
