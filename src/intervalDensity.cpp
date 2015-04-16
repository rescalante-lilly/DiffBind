#include <exception>
#include <stdexcept>
#include <algorithm>
#include <sstream>

#include "intervalDensity.h"

bode::IntervalDensity::IntervalDensity(const int n) {
  size = n;
  counts = new unsigned int[n];
  for (int i=0;i<size;i++) {
    counts[i] = 0;
  }
  reads = 0;
}

bode::IntervalDensity::~IntervalDensity(void) {
  delete[] counts;
}

void bode::IntervalDensity::set(const int left,const int right) {
  if (left >= size || right <= 0) {
    throw std::out_of_range("invalid range on set");
  }
  int l = std::max(0,left);
  int r = std::min(size,right);
  for (int i=l;i<r;i++) {
    counts[i] += 1;
  }
  reads += 1;
}

unsigned int bode::IntervalDensity::depth(const int i) {
  if (i < 0 || i >= size) {
    throw std::out_of_range("invalid range on depth");
  }
  return counts[i];
}

void bode::IntervalDensity::clear(void) {
  for (int i=0; i<size; i++) {
    counts[i] = 0;
  }
}

unsigned int bode::IntervalDensity::operator[](const int i) {
  if (i < 0 || i >= size) {
    throw std::out_of_range("invalid range on array ref");
  }
  return counts[i];
}

void bode::IntervalDensity::summit(int &position,unsigned int &height) {
  unsigned int maxposf = 0,maxposr;
  unsigned int maxheightf = 0,maxheightr = 0;
  int posfwd,posrev;

  for (posfwd=0; posfwd<size; posfwd++) {
    if (counts[posfwd] > maxheightf) {
      maxheightf = counts[posfwd];
      maxposf = posfwd;
    }
  }
  maxposr = size-1;
  for (posrev=size-1;posrev>=0;posrev--) {
    if (counts[posrev] > maxheightr) {
      maxheightr = counts[posrev];
      maxposr = posrev;
    }
  }
  position = (maxposf + maxposr) / 2;
  height = maxheightf;
}

int bode::IntervalDensity::nreads(void) {
  return reads;
}

std::string bode::IntervalDensity::str(void) {
  std::ostringstream msg;
  for (int i=0;i<size;i++) {
    msg << counts[i] << " ";
  }
  return msg.str();
}
