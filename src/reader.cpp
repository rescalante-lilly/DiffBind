#include <string>
#include <R_ext/Error.h>

#include "interval.h"
#include "reader.h"
#include "bamReader.h"
#include "bedReader.h"

//class BamReader;
//class BedReader;

bode::Reader::~Reader(void) {
  /* this function exists solely to prevent a weird compiler error. */
}

bode::Reader *bode::Reader::open(std::string const &filename) {
  Reader *r = NULL;
  int flen;

  // Ugly.  Must fix.
  // See http://stackoverflow.com/questions/582331/c-is-there-a-way-to-instantiate-objects-from-a-string-holding-their-class-name
  // for a fix.

  flen = filename.length();
  if (filename.compare(flen-4,4,".bam") == 0) {
    r = new BamReader(filename);
  } else if (filename.compare(flen-4,4,".bed") == 0) {
    r = new BedReader(filename);
  } else if (filename.compare(flen-7,7,".bed.gz") == 0) {
    r = new BedReader(filename);
  } else {
    error("Unknown suffix in file '%s'.  Supported are: '.bam', '.bed', '.bed.gz'.",filename.c_str());
  }
  return r;
}
