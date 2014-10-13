#ifndef __SEQUENCE_H
#define __SEQUENCE_H

#include <ostream>
#include <string>

namespace bode {

class Sequence {
  public:
    Sequence(std::string const &n,std::string const &s);
    Sequence() { _null = true; _name = ""; _seq = ""; };
    Sequence(Sequence const &other);
    std::string seq(void) const            { return _seq; };
    std::string name(void) const           { return _name; };
    void setSeq(std::string const &s)      { _seq = s; };
    void setName(std::string const &n)     { _name = n; };
    virtual void write(std::ostream &out);
    virtual std::string format(void) const;
    bool isNull(void) const                { return _null; };
    void setNull(bool f)                   { _null = f; };
    friend bool operator==(Sequence const &a,Sequence const &b) {
      return a._seq == b._seq;
    };
    friend bool operator!=(Sequence const &a,Sequence const &b) {
      return a._seq != b._seq;
    };
    Sequence& operator=(Sequence const &other) {
      if (this != &other) {
        _name = other._name;
        _seq = other._seq;
        _null = other._null;
      }
      return *this;
    };
    bool operator!(void) const             { return !_null; };

  protected:
    std::string _name;
    std::string _seq;
    bool _null;
};

}

#endif
