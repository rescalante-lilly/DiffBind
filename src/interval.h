#ifndef __INTERVAL_H
#define __INTERVAL_H

#include <string>

namespace bode {

class Interval {
  public:
    Interval(std::string const &chr,int l,int r);
    Interval(std::string const &chr,int l,int r,int strand);
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
    int strand(void) const                                 { return _strand; };

    virtual void updatecstr(const char *c_str,int l,int r);
    virtual void update(std::string const &chr,int l,int r);
    virtual void update(std::string const &chr,int l,int r,int s);
    virtual std::string format(void) const;
    virtual void setUnmapped(void)                         { _mapped = false; };
    virtual void extend(int insertLen);

  protected:
    int _left;
    int _right;
    std::string _chrom;
    bool _mapped;
    int _strand;

};

bool operator==(Interval const &l,Interval const &r);
bool operator<(Interval const &l,Interval const &r);

}

#endif
