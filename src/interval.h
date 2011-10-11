#ifndef __INTERVAL_H
#define __INTERVAL_H

#include <string>

namespace bode {

class Interval {
  public:
    Interval(std::string const &chr,int l,int r);
    Interval(void)                                         { _mapped = false; };
    Interval(Interval const &i);
    virtual ~Interval(void) {};
    Interval &operator=(Interval const &i);
    friend bool operator==(Interval const &l,Interval const &r);
    friend bool operator<(Interval const &l,Interval const &r);

    std::string const& chrom(void) const                   { return _chrom; };
    int left(void) const                                   { return _left; }; 
    int right(void) const                                  { return _right; };
    bool isMapped(void) const                              { return _mapped; };

    virtual void update(std::string const &chr,int l,int r);
    virtual std::string format(void) const;
    virtual void setUnmapped(void)                         { _mapped = false; };
    virtual int score(void) const;
    virtual void extend(int insertLen);

  protected:
    int _left;
    int _right;
    std::string _chrom;
    bool _mapped;

};

bool operator==(Interval const &l,Interval const &r);
bool operator<(Interval const &l,Interval const &r);

}

#endif
