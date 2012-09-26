#ifndef __READER_H
#define __READER_H

#include <string>

#include "interval.h"

namespace bode {

class Reader {
  public:
    virtual ~Reader(void);
    virtual Interval *next(void) = 0;
    virtual void close(void) = 0;
    virtual bool eof(void) = 0;
    static Reader *open(std::string const &filename);
};

}

#endif
