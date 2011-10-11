#ifndef __BEDREADER_H
#define __BEDREADER_H

#include <string.h>

#include "interval.h"
#include "reader.h"
#include "bed.h"
#include "writer.h"

namespace bode {

class BedReader: public Reader {
  public:
    BedReader(std::string const &filename);
    ~BedReader(void);

    Interval *next(void);
    void close(void);
    bool eof(void)                                             { return _eof; };
    void *header(void)                                         { return NULL; };
    static BedReader *open(std::string const &filename);
    Writer *getWriter(void)                         { return new Writer(); };

  private:
    static int const maxLine = 1024;
    gzFile _fd;
    Bed *_bseq;
    char *_buffer;
    bool _eof;
};

}

#endif
