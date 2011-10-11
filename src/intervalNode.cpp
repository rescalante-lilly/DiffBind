#include <cstddef>
#include "intervalNode.h"

/******************************************************************************/

bode::IntervalNode::IntervalNode(int left,int right) {
  leftEnd = left;
  rightEnd = right;
  count = 1;
  redFlag = 1;
  leftChild = NULL;
  rightChild = NULL;
  parent = NULL;
}

bode::IntervalNode::~IntervalNode() {
  if (leftChild != NULL) {
    delete leftChild;
  }
  if (rightChild != NULL) {
    delete rightChild;
  }
}

/******************************************************************************/
