#ifndef __BAM_H
#define __BAM_H

/**
 * \file
 * \brief File containing the BAM class header.
 *
 * This is a test.
 */

#include "samtools/bam.h"

#include "interval.h"
#include "sequence.h"
#include "bamWriter.h"

namespace bode {

/**
 * \brief Bam -- holds a BAM sequence.
 *
 * This is more information about the BAM class. */
class Bam: public Interval, public Sequence {
  public:
    Bam(bam1_t *raw,bam_header_t *hdr);
    Bam(void);
    Bam(Bam const &b);
    Bam &operator=(Bam const &b);
    ~Bam(void);
    friend bool operator==(Bam const &l,Bam const &r);
    friend bool operator<(Bam const &l,Bam const &r);

    int seq(std::string &dest) const;
    void setHeader(bam_header_t *hdr)                           { _hdr = hdr; };
    void update(bam1_t *raw);
    void update(std::string const &chr,int l,int r);
    void setUnmapped(void);
    friend void bode::BamWriter::write(Interval const &i);
    virtual void extend(int insertLen);

  protected:
    bam1_t *_raw;
    bam_header_t *_hdr;

    void setInterval(void);
    int32_t chromIndex(std::string const &chr);

    static int nucleotideMap[];

};

}

#endif
