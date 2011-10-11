#ifndef __BAMREADER_H
#define __BAMREADER_H

#include <samtools/sam.h>

#include "interval.h"
#include "reader.h"
#include "bam.h"
#include "writer.h"
#include "bamWriter.h"

namespace bode {

class BamReader: public Reader {
  public:
    BamReader(std::string const &filename);
    virtual ~BamReader(void);

    Interval *next(void);
    void close(void);
    bool eof(void)                                      { return _eof; };
    void *header(void)                                  { return _fd->header; };
    static BamReader *open(std::string const &filename);
    Writer *getWriter(void)                         { return new BamWriter(); };

  private:
    samfile_t *_fd;
    bam1_t *_seq;
    Bam *_bseq;
    bool _eof;
};

}

#endif
