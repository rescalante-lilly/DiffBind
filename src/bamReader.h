#ifndef __BAMREADER_H
#define __BAMREADER_H

#include <samtools/sam.h>

#include "interval.h"
#include "reader.h"

namespace bode {

class BamReader: public Reader {
  public:
    BamReader(std::string const &filename);
    virtual ~BamReader(void);

    Interval *next(void);
    void close(void);
    bool eof(void)                                      { return _eof; };
    static BamReader *open(std::string const &filename);

  private:
    samfile_t *_fd;
    bam1_t *_seq;
    Interval *_bseq;
    bool _eof;
};

}

#endif
