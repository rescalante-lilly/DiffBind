#include <cstddef>
#include "intervalNode.h"
#include "intervalTree.h"

// #define ISRED(x) ((x)->isRed==1)
// #define SETRED(x) ((x)->isRed = 1)
// #define SETBLACK(x) ((x)->isRed = 0)

/******************************************************************************/
/* interval set
 *
 */
bode::IntervalTree::IntervalTree(void) {
  root = NULL;
  count = 0;
}

bode::IntervalTree::~IntervalTree(void) {
  if (root != NULL) {
    delete root;
  }
}

#if 0
static int traverseRange(inodeP node,int left,int right) {
  int sum = 0;
  if (node->left != NULL && node->rightEnd >= left) {
    sum += traverseRange(node->left,left,right);
  }
  if (overlap(node,left,right)) {
    sum += node->count;
  }
  if (node->right != NULL && node->leftEnd <= right) {
    sum += traverseRange(node->right,left,right);
  }
  return sum;
}

int countIntervals(isetP tree,int left,int right) {
  return traverseRange(tree->root,left,right);
}
#endif

bode::IntervalNode *bode::IntervalTree::raw_insert(IntervalNode *node) {
  bode::IntervalNode *x,*y;
  x = root;
  y = NULL;
  while (x != NULL) {
    y = x;
    if (*node < *x) {
      x = x->getLeft();
    } else if (*node > *x) {
      x = x->getRight();
    } else { /* duplicate intervals */
      x->incrementCount();
      return x;
    }
  }
  node->setParent(y);
  if (y == NULL) {
    root = node;
  } else {
    if (*node < *y) {
      y->setLeft(node);
    } else {
      y->setRight(node);
    }
  }
  return node;
}

void bode::IntervalTree::leftRotate(IntervalNode *node) {
  bode::IntervalNode *y = node->getRight();
  node->setRight(y->getLeft());
  if (!y->nullLeft()) {
    y->getLeft()->setParent(node);
  }
  y->setParent(node->getParent());
  if (node->nullParent()) {
    root = y;
  } else if (node == node->getParent()->getLeft()) {
    node->getParent()->setLeft(y);
  } else {
    node->getParent()->setRight(y);
  }
  y->setLeft(node);
  node->setParent(y);
}

void bode::IntervalTree::rightRotate(IntervalNode *node) {
  bode::IntervalNode *y = node->getLeft();
  node->setLeft(y->getRight());
  if (!y->nullRight()) {
    y->getRight()->setParent(node);
  }
  y->setParent(node->getParent());
  if (node->nullParent()) {
    root = y;
  } else if (node == node->getParent()->getRight()) {
    node->getParent()->setRight(y);
  } else {
    node->getParent()->setLeft(y);
  }
  y->setRight(node);
  node->setParent(y);
}

void bode::IntervalTree::rebalance(IntervalNode *x) {
  while (x!=root && x->getParent()->isRed()) {
    if (x->getParent()->getParent()->getLeft() == x->getParent()) {
      bode::IntervalNode *y = x->getParent()->getParent()->getRight();
      if (y != NULL && y->isRed()) {
        x->getParent()->setBlack();
        y->setBlack();
        x->getParent()->getParent()->setRed();
        x = x->getParent()->getParent();
      } else {
        if (x == x->getParent()->getRight()) {
          x = x->getParent();
          leftRotate(x);
        }
        x->getParent()->setBlack();
        x->getParent()->getParent()->setRed();
        rightRotate(x->getParent()->getParent());
      }
    } else {
      bode::IntervalNode *y = x->getParent()->getParent()->getLeft();
      if (y != NULL && y->isRed()) {
        x->getParent()->setBlack();
        y->setBlack();
        x->getParent()->getParent()->setRed();
        x = x->getParent()->getParent();
      } else {
        if (x == x->getParent()->getLeft()) {
          x = x->getParent();
          rightRotate(x);
        }
        x->getParent()->setBlack();
        x->getParent()->getParent()->setRed();
        leftRotate(x->getParent()->getParent());
      }
    }
  }
  root->setBlack();
}

int bode::IntervalTree::i_coverage(bode::IntervalNode *n,int point) {
  if (n == NULL) {
    return 0;
  } else if (point >= n->r()) {
    return i_coverage(n->getRight(),point);
  } else if (point < n->l()) {
    return i_coverage(n->getLeft(),point);
  } else {
    return (  i_coverage(n->getLeft(),point)
            + n->getCount()
            + i_coverage(n->getRight(),point));
  }
}

int bode::IntervalTree::i_countIntervals(bode::IntervalNode *n,int left,int right) {
  if (n == NULL) {
    return 0;
  } else if (left >= n->r()) {
    return i_countIntervals(n->getRight(),left,right);
  } else if (right <= n->l()) {
    return i_countIntervals(n->getLeft(),left,right);
  } else {
    return (  i_countIntervals(n->getLeft(),left,right)
            + n->getCount()
            + i_countIntervals(n->getRight(),left,right));
  }
}

int bode::IntervalTree::i_realCount(bode::IntervalNode *n) {
  if (n == NULL) {
    return 0;
  } else {
    return i_realCount(n->getLeft())
           + n->getCount()
           + i_realCount(n->getRight());
  }
}

/******************************************************************************/

void bode::IntervalTree::insert(int left,int right) {
  bode::IntervalNode *x = new bode::IntervalNode(left,right);
  bode::IntervalNode *y;
  y = raw_insert(x);
  if (y != x) { // this means we found a duplicate interval and just incremented
                // its counter instead of inserting the new interval
    delete x;
    count++;
    return;
  }
  rebalance(x);
  count++;
}

int bode::IntervalTree::coverage(int point) {
  return i_coverage(root,point);
}

int bode::IntervalTree::reads(int left,int right) {
  return i_countIntervals(root,left,right);
}

int bode::IntervalTree::summit(int left,int right) {
  return 0;
}

int bode::IntervalTree::getCount(void) {
  return count;
}

int bode::IntervalTree::realCount(void) {
  return i_realCount(root);
}
