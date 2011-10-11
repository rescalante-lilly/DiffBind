#ifndef __BAMWRITER_H
#define __BAMWRITER_H

#include <samtools/sam.h>

#include "writer.h"

namespace bode {

class BamWriter: public Writer {
  public:
    BamWriter(std::string const &filename,void *header);
    BamWriter(void);
    virtual ~BamWriter(void);

    void close(void);
    bool isOpen(void) { return _open; };
    void write(Interval const &i);
    Writer *open(std::string const &filename,void *header);

  protected:
    samfile_t *_fd;
    bool _open;
};

}

#endif
