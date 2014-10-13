#ifndef __INTERVALNODE_H
#define __INTERVALNODE_H

#include <algorithm>
/* #define max(a,b) ((a)>(b)?(a):(b)) */
/* #define min(a,b) ((a)>(b)?(b):(a)) */

namespace bode {

class IntervalNode {
  public:
//    ~IntervalNode();
    void init(int left,int right,int s);
    void incrementCountF()          { countFwd++; };
    void incrementCountR()          { countRev++; };
    int getCount()                  { return countFwd+countRev; };
    int getCountFwd()               { return countFwd; };
    int getCountRev()               { return countRev; };
    void setParent(IntervalNode *p) { parent = p; };
    void setLeft(IntervalNode *l)   { leftChild = l; };
    void setRight(IntervalNode *r)  { rightChild = r; };
    bool nullParent(void)           { return parent == NULL; };
    bool nullLeft(void)             { return leftChild == NULL; };
    bool nullRight(void)            { return rightChild == NULL; };
    IntervalNode *getParent(void)   { return parent; };
    IntervalNode *getLeft(void)     { return leftChild; };
    IntervalNode *getRight(void)    { return rightChild; };
    int l(void)                     { return leftEnd; };
    int r(void)                     { return rightEnd; };
    void setBlack(void)             { redFlag = false; };
    void setRed(void)               { redFlag = true; };
    bool isRed(void)                { return redFlag; };
    bool overlap(int left,int right);
    bool operator==(const IntervalNode &other);
    bool operator<(const IntervalNode &other);
    bool operator>(const IntervalNode &other);

  private:
    int leftEnd;
    int rightEnd;
    int countFwd;
    int countRev;
    bool redFlag;
    IntervalNode *leftChild;
    IntervalNode *rightChild;
    IntervalNode *parent;
};

inline bool IntervalNode::overlap(int left,int right) {
  int ml = std::max(leftEnd,left);
  int mr = std::min(rightEnd,right);
  return (mr - ml) > 0;
}

inline bool IntervalNode::operator==(const IntervalNode &other) {
  return other.leftEnd == leftEnd
         && other.rightEnd == rightEnd;
}

inline bool IntervalNode::operator<(const IntervalNode &other) {
  return leftEnd < other.leftEnd ||
         (leftEnd == other.leftEnd && rightEnd < other.rightEnd);
}

inline bool IntervalNode::operator>(const IntervalNode &other) {
  return leftEnd > other.leftEnd ||
         (leftEnd == other.leftEnd && rightEnd > other.rightEnd);
}

}

#endif
