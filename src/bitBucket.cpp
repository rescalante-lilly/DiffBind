#include <R.h>

#include "bitBucket.h"

namespace bode {

BitBucket::BitBucket(const int n) {
  size = n;
  bits = new unsigned char[n];
}

BitBucket::~BitBucket(void) {
  delete bits;
}

void BitBucket::set(const int i) {
  if (i < 0 || i >= size) {
    error("Range error in BitBucket: %d (limit=%d)",i,size);
  }
  bits[i] = 1;
}

void BitBucket::reset(const int i) {
  if (i < 0 || i >= size) {
    error("Range error in BitBucket: %d (limit=%d)",i,size);
  }
  bits[i] = 0;
}

bool BitBucket::isSet(const int i) {
  if (i < 0 || i >= size) {
    error("Range error in BitBucket: %d (limit=%d)",i,size);
  }
  return bits[i] == 1;
}

void BitBucket::clear(void) {
  for (int i=0; i<size; i++) {
    bits[i] = 0;
  }
}

bool BitBucket::operator[](const int i) {
  if (i < 0 || i >= size) {
    error("Range error in BitBucket: %d (limit=%d)",i,size);
  }
  return bits[i] == 1;
}

}
