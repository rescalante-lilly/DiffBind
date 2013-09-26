#include <stdlib.h>
#include <string>
#include <fstream>
#include <zlib.h>
#include <iostream>
#include <ctype.h>
#include <R.h>

#include "interval.h"
#include "util.h"
#include "bedReader.h"

bool isDigits(char *s) {
  int slen = strlen(s);
  bool okay = true;
  for (int i=0;i<slen;i++) {
    if (!isdigit(s[i])) {
      okay = false;
      break;
    }
  }
  return okay;
}

bool bode::BedReader::isBed(std::string const &filename) {
  char *res;
  gzFile fd;
  char buffer[maxLine];
  bool okay = true;
  int count,lines;
  char *fields[12];
  bool digits;

  fd = gzopen(filename.c_str(),"r");
  for (lines = 0;lines < 10; lines++) { // check first 10 lines
    res = gzgets(fd,buffer,maxLine);
    if (res == NULL) {
      okay = false;
      break;
    }
    if (lines == 0 && strncmp(buffer,"track",5) == 0) {
      continue;
    }
    bode::trimTrailing(buffer);
    count = bode::splits(buffer,fields,12);
    if (count < 3 || (count > 6 && count != 12)) {
      okay = false;
      break;
    };
    if (!isDigits(fields[1]) || !isDigits(fields[2])) {
      okay = false;
      break;
    }
  }
  gzclose(fd);
  return okay;
}

bode::BedReader::BedReader(std::string const &filename) {
  char *res;

  if (!isBed(filename)) {
    error("file '%s' does not appear to be a BED file (coordinates are not integers)",filename.c_str());
  }
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
