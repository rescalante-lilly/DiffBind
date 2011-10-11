#ifndef __BED_H
#define __BED_H

#include <string>

#include "interval.h"

namespace bode {

class Bed: public Interval {
  public:
    Bed(std::string const &chr,int l,int r);
    Bed(std::string const &chr,int l,int r,std::string n,int score,char s);
    Bed(void);
    Bed(Bed const &b);

    std::string const& name(void) const { return _name; };
    char strand(void) const             { return _strand; };
    int score(void) const               { return _score; };

    virtual void update(std::string const &chr,int l,int r);
    virtual void update(std::string const &chr,int l,int r,std::string const &n,int score,char s);
    virtual std::string format(void) const;
    virtual void extend(int insertLen);

  protected:
    std::string _name;
    char _strand;
    int _score;
};

}

#endif
