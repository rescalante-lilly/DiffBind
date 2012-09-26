#ifndef __BEDREADER_H
#define __BEDREADER_H

#include <string.h>

#include "interval.h"
#include "reader.h"

namespace bode {

class BedReader: public Reader {
  public:
    BedReader(std::string const &filename);
    ~BedReader(void);

    Interval *next(void);
    void close(void);
    bool eof(void)                                             { return _eof; };
    static BedReader *open(std::string const &filename);

  private:
    static int const maxLine = 1024;
    gzFile _fd;
    Interval *_bseq;
    char *_buffer;
    bool _eof;
};

}

#endif
