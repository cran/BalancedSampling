#ifndef KDTREECLASS_HEADER
#define KDTREECLASS_HEADER

#include <stddef.h>
#include <vector>

#include "KDNodeClass.h"
#include "KDStoreClass.h"

enum class KDTreeSplitMethod {
  variable = 0,
  maximalSpread = 1,
  midpointSlide = 2
};

KDTreeSplitMethod IntToKDTreeSplitMethod(const int);

class KDTree {
protected:
  double* data; // data array of length Np
  size_t N;
  size_t p;
  size_t bucketSize;
  KDTreeSplitMethod method = KDTreeSplitMethod::midpointSlide;
  size_t (KDTree::*SplitFindSplitUnit)(KDNode*, size_t*, const size_t) = nullptr;

public:
  KDNode* topNode = nullptr;

protected:
  KDTree();
public:
  KDTree(double*, const size_t, const size_t, const size_t, const KDTreeSplitMethod);
  ~KDTree();
  KDTree* Copy();
  void Prune();

protected:
  std::vector<double> liml = std::vector<double>(0);
  std::vector<double> limr = std::vector<double>(0);
  void SplitNode(KDNode*, size_t*, const size_t);
  size_t SplitByVariable(KDNode*, size_t*, const size_t);
  size_t SplitByMaximalSpread(KDNode*, size_t*, const size_t);
  size_t SplitByMidpointSlide(KDNode*, size_t*, const size_t);
  size_t SplitUnitsById(size_t*, const size_t, const size_t, const size_t);

public:
  KDNode* FindNode(const size_t);
  bool UnitExists(const size_t);
  void RemoveUnit(const size_t);
  double DistanceBetweenUnits(const size_t, const size_t);
private:
  double DistanceBetweenPointers(const double*, const double*);

public:
  void FindNeighbours(KDStore*, const size_t);
  void FindNeighbours(KDStore*, const double*);
private:
  void TraverseNodesForNeighbours(KDStore*, const size_t, const double*, KDNode*);
  void SearchNodeForNeighbour1(KDStore*, const size_t, const double*, KDNode*);
  void SearchNodeForNeighbours(KDStore*, const size_t, const double*, KDNode*);

public:
  void FindNeighboursCps(KDStore*, const std::vector<double>&, const size_t);
private:
  void TraverseNodesForNeighboursCps(KDStore*, const std::vector<double>&, const size_t, const double*, KDNode*, double*);
  void SearchNodeForNeighboursCps(KDStore*, const std::vector<double>&, const size_t, const double*, KDNode*, double*);
};

#endif
