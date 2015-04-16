#ifndef __INTERVALDENSITY_H
#define __INTERVALDENSITY_H

namespace bode {

class IntervalDensity {
  public:
    IntervalDensity(const int n);
    ~IntervalDensity(void);
    void set(const int left,const int right);
    unsigned int depth(const int i);
    void clear(void);
    unsigned int operator[](const int i);
    void summit(int &position, unsigned int &height);
    int nreads(void);
    std::string str(void);
  private:
    int size;
    unsigned int *counts;
    int reads;
};

}

#endif
