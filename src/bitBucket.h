#ifndef __BITBUCKET_H
#define __BITBUCKET_H

class BitBucket {
  public:
    BitBucket(const int n);
    ~BitBucket(void);
    void set(const int i);
    void reset(const int i);
    bool isSet(const int i);
    void clear(void);
    bool operator[](const int i);
  private:
    int size;
    unsigned char *bits;
};

#endif
