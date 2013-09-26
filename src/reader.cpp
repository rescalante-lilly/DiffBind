#include <string>
#include <R_ext/Error.h>

#include "interval.h"
#include "reader.h"
#include "bamReader.h"
#include "bedReader.h"

//class BamReader;
//class BedReader;

#define NO_FILETYPE 0
#define BED_FILETYPE 1
#define BAM_FILETYPE 3

bode::Reader::~Reader(void) {
  /* this function exists solely to prevent a weird compiler error. */
}

bode::Reader *bode::Reader::open(std::string const &filename,int const &filetype) {
  Reader *r = NULL;
  int flen;

  // Ugly.  Must fix.
  // See http://stackoverflow.com/questions/582331/c-is-there-a-way-to-instantiate-objects-from-a-string-holding-their-class-name
  // for a fix.

  if (filetype == NO_FILETYPE) {
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
  } else {
    if (filetype == BED_FILETYPE) {
      r = new BedReader(filename);
    } else if (filetype == BAM_FILETYPE) {
      r = new BamReader(filename);
    } else {
      error("Unknown filetype %d in file '%s'.  Supported are 0 (use suffix), 1 (bed), 3 (bam).",filetype,filename.c_str());
    }
  }
  return r;
}
