#include <cstddef>
#include "intervalNode.h"

/******************************************************************************/

/*
bode::IntervalNode::IntervalNode(int left,int right,int s) {
  leftEnd = left;
  rightEnd = right;
  if (s > 0) { // positive strand
    countFwd = 1;
    countRev = 0;
  } else {
    countFwd = 0;
    countRev = 1;
  }
  redFlag = 1;
  leftChild = NULL;
  rightChild = NULL;
  parent = NULL;
}
*/

void bode::IntervalNode::init(int left,int right,int s) {
  leftEnd = left;
  rightEnd = right;
  if (s > 0) { // positive strand
    countFwd = 1;
    countRev = 0;
  } else {
    countFwd = 0;
    countRev = 1;
  }
  redFlag = 1;
  leftChild = NULL;
  rightChild = NULL;
  parent = NULL;
}

/*
bode::IntervalNode::~IntervalNode() {
  if (leftChild != NULL) {
    delete leftChild;
  }
  if (rightChild != NULL) {
    delete rightChild;
  }
}
*/

/******************************************************************************/
