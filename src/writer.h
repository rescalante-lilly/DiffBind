#ifndef __WRITER_H
#define __WRITER_H

#include <string>
#include <fstream>

#include "interval.h"

namespace bode {

class Writer {
  public:
    Writer(std::string const &filename,void *header);
    Writer(void);
    virtual ~Writer(void);

    virtual void close(void);
    virtual bool isOpen(void)                                 { return _open; };
    virtual void write(Interval const &i);
    virtual Writer *open(std::string const &filename,void *header);

  protected:
    std::ofstream *_fd;
    bool _open;
};

}

#endif
