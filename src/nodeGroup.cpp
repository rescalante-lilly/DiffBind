#include <R.h>
#include <Rinternals.h>
#include <stdio.h>
#include "intervalNode.h"
#include "nodeGroup.h"

bode::NodeGroup::NodeGroup(int size) {
  nsize = size;
  ncount = 0;
  nodes = new bode::IntervalNode[nsize];
}

bode::NodeGroup::~NodeGroup(void) {
  delete[] nodes;
  nodes = NULL;
  nsize = 0;
  ncount = 0;
}

void bode::NodeGroup::clear(void) {
  ncount = 0;
}

bode::IntervalNode *bode::NodeGroup::get(void) {
  bode::IntervalNode *n;
  n = &nodes[ncount];
  ncount += 1;
  return n;
}

void bode::NodeGroup::pop(void) {
  ncount -= 1;
}
