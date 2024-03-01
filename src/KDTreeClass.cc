#include <algorithm>
#include <float.h>
#include <stddef.h>
#include <stdexcept>
#include <utility>
#include <vector>

#include "KDNodeClass.h"
#include "KDStoreClass.h"
#include "KDTreeClass.h"
#include "utils-matrix.h"

KDTreeSplitMethod IntToKDTreeSplitMethod(const int i) {
  if (0 <= i && i <= 2)
    return static_cast<KDTreeSplitMethod>(i);

  throw std::invalid_argument("split method does not exist");
  return KDTreeSplitMethod::midpointSlide;
}

// General constructor for KDTree
KDTree::KDTree(
  double* t_dt,
  const size_t t_N,
  const size_t t_p,
  const size_t t_bucketSize,
  const KDTreeSplitMethod t_method
) {
  Init(t_dt, t_N, t_p, t_bucketSize, t_method);

  size_t* splitUnits = new size_t[N];
  double* dt = data;
  for (size_t i = 0; i < N; i++) {
    splitUnits[i] = i;

    // Set limits
    for (size_t k = 0; k < p; k++, dt++) {
      if (*dt < liml[k])
        liml[k] = *dt;
      if (*dt > limr[k])
        limr[k] = *dt;
    }
  }

  if (N <= bucketSize) {
    // All units are in top
    topNode = new KDNode(nullptr, true);
    topNode->ReplaceUnits(splitUnits, N);
    delete[] splitUnits;
    return;
  }

  // We need to start splitting
  topNode = new KDNode(nullptr, false);
  SplitNode(topNode, splitUnits, N);

  delete[] splitUnits;
}

// Constructor for custom indices
KDTree::KDTree(
  double* t_dt,
  const size_t t_N,
  const size_t t_p,
  const size_t t_bucketSize,
  const KDTreeSplitMethod t_method,
  size_t* splitUnits,
  size_t splitUnitsN
) {
  Init(t_dt, t_N, t_p, t_bucketSize, t_method);

  // Set limits of current data
  for (size_t i = 0; i < splitUnitsN; i++) {
    double* dt = data + splitUnits[i] * p;

    for (size_t k = 0; k < p; k++, dt++) {
      if (*dt < liml[k])
        liml[k] = *dt;
      if (*dt > limr[k])
        limr[k] = *dt;
    }
  }

  // If all existing units fits in a bucket, we're done before we started ...
  if (N <= bucketSize) {
    topNode = new KDNode(nullptr, true);
    topNode->ReplaceUnits(splitUnits, splitUnitsN);
    return;
  }

  // ... otherwise, we need to start splitting
  topNode = new KDNode(nullptr, false);
  SplitNode(topNode, splitUnits, splitUnitsN);

  // delete[] splitUnits left for the caller
}

// Common initializers
void KDTree::Init(
  double* t_dt,
  const size_t t_N,
  const size_t t_p,
  const size_t t_bucketSize,
  const KDTreeSplitMethod t_method
) {
  data = t_dt;

  N = t_N;
  if (N < 1)
    throw std::invalid_argument("(init) N to small");

  p = t_p;
  if (p < 1)
    throw std::invalid_argument("(init) p to small");

  liml.resize(p, DBL_MAX);
  limr.resize(p, -DBL_MAX);

  bucketSize = t_bucketSize;
  if (bucketSize < 1)
    throw std::invalid_argument("(init) bucketSize to small");

  method = t_method;

  switch(method) {
  case KDTreeSplitMethod::variable:
    SplitFindSplitUnit = &KDTree::SplitByVariable;
    break;
  case KDTreeSplitMethod::maximalSpread:
    SplitFindSplitUnit = &KDTree::SplitByMaximalSpread;
    break;
  case KDTreeSplitMethod::midpointSlide:
    SplitFindSplitUnit = &KDTree::SplitByMidpointSlide;
    break;
  default:
    throw std::invalid_argument("split method does not exist");
    return;
  }
}

// Protected
KDTree::KDTree() {};

KDTree::~KDTree() {
  delete topNode;
}

KDTree* KDTree::Copy() {
  KDTree* newTree = new KDTree();
  newTree->data = data;
  newTree->N = N;
  newTree->p = p;
  newTree->bucketSize = bucketSize;
  newTree->method = method;
  newTree->SplitFindSplitUnit = SplitFindSplitUnit;
  newTree->liml.reserve(p);
  newTree->limr.reserve(p);

  for (size_t k = 0; k < p; k++) {
    newTree->liml.push_back(liml[k]);
    newTree->limr.push_back(limr[k]);
  }

  newTree->topNode = new KDNode(nullptr, topNode->IsTerminal());
  newTree->topNode->Copy(topNode);

  return newTree;
}
void KDTree::Prune() {
  if (topNode == nullptr)
    return;

  topNode->Prune(bucketSize);
  return;
}

void KDTree::SplitNode(KDNode* node, size_t* splitUnits, const size_t n) {
  // m should be the index of the last unit to not include in cleft
  // i.e. cleft = [0, m), cright = [m, n)
  size_t m = (this->*SplitFindSplitUnit)(node, splitUnits, n);

  // Sanity check
  if (m > n) {
    throw std::range_error("(SplitNodes) m > n");
    return;
  }

  // If m is 0 or n, we need to accept all units into the node
  if (m == 0 || m >= n) {
    node->SetTerminal(true);
    node->ReplaceUnits(splitUnits, n);
    return;
  }

  if (m <= bucketSize) {
    node->cleft = new KDNode(node, true);
    node->cleft->ReplaceUnits(splitUnits, m);
  } else {
    node->cleft = new KDNode(node, false);
    SplitNode(node->cleft, splitUnits, m);
  }

  size_t rSize = n - m;
  if (rSize <= bucketSize) {
    node->cright = new KDNode(node, true);
    node->cright->ReplaceUnits(splitUnits + m, rSize);
  } else {
    node->cright = new KDNode(node, false);
    SplitNode(node->cright, splitUnits + m, rSize);
  }

  return;
}

size_t KDTree::SplitByVariable(KDNode* node, size_t* splitUnits, const size_t n) {
  // We need to fix this function so it can handle when a split level is badly chosen
  // because it has no variability, but another split level does have variability.
  // Currently, it just stops when this happens, i.e. satisfied with creating a
  // possibly very large terminal node.
  size_t level = 0;
  for (KDNode* nd = node; nd->parent != nullptr; nd = nd->parent)
    level += 1;

  // Set split variable to the mod
  node->split = level % p;

  // n >> 1 tries to split the units by the median
  size_t m = SplitUnitsById(splitUnits, n, n >> 1, node->split);
  node->value = data[MatrixIdxRM(splitUnits[m - 1], node->split, p)];

  return m;
}

// splitSpread
size_t KDTree::SplitByMaximalSpread(KDNode* node, size_t* splitUnits, const size_t n) {
  double* mins = new double[p];
  double* maxs = new double[p];

  double* dt = data + splitUnits[0] * p;
  for (size_t k = 0; k < p; k++) {
    mins[k] = dt[k];
    maxs[k] = dt[k];
  }

  // Find the mins/maxs in curretn splitUnits
  for (size_t i = 1; i < n; i++) {
    dt = data + splitUnits[i] * p;
    for (size_t k = 0; k < p; k++) {
      if (dt[k] < mins[k])
        mins[k] = dt[k];
      else if (dt[k] > maxs[k])
        maxs[k] = dt[k];
    }
  }

  // Decide the splitting variable by finding the variable with the largest spread
  node->split = 0;
  double spread = maxs[0] - mins[0];
  for (size_t k = 1; k < p; k++) {
    double temp = maxs[k] - mins[k];
    if (temp > spread) {
      node->split = k;
      spread = temp;
    }
  }

  delete[] mins;
  delete[] maxs;

  // If there is no spread in any variable, we shouldn't split more
  if (spread == 0.0)
    return 0;

  // Find the median unit for splitting
  size_t m = SplitUnitsById(splitUnits, n, n >> 1, node->split);
  node->value = data[MatrixIdxRM(splitUnits[m - 1], node->split, p)];
  return m;
}

size_t KDTree::SplitByMidpointSlide(KDNode* node, size_t* splitUnits, const size_t n) {
  double* mins = new double[p];
  double* maxs = new double[p];

  for (size_t k = 0; k < p; k++) {
    mins[k] = liml[k];
    maxs[k] = limr[k];
  }

  // Get current window's limits by traversing up the tree
  for (KDNode* par = node; par->parent != nullptr; par = par->parent) {
    if (par->parent->cleft == par) {
      if (par->parent->value < maxs[par->parent->split])
        maxs[par->parent->split] = par->parent->value;
    } else {
      if (par->parent->value > mins[par->parent->split])
        mins[par->parent->split] = par->parent->value;
    }
  }

  // Decide the splitting variable by finding the variable with the largest window
  node->split = 0;
  double spread = maxs[0] - mins[0];
  for (size_t k = 1; k < p; k++) {
    double temp = maxs[k] - mins[k];
    if (temp > spread) {
      node->split = k;
      spread = temp;
    }
  }

  // Decide a candidate splitting value
  node->value = (maxs[node->split] + mins[node->split]) * 0.5;

  delete[] mins;
  delete[] maxs;

  // If there is no spread in any variable, we shouldn't split more
  if (spread == 0.0)
    return 0;

  double* dt = data + node->split;
  size_t l = 0;
  size_t r = n;
  double lbig = -DBL_MAX;
  double rsmall = DBL_MAX;

  // Sort splitUnits so that we have
  // x <= value is in range [0, l)
  // x > value is in range [r, n)
  // where value is the proposed split
  while (l < r) {
    double temp = *(dt + splitUnits[l] * p);
    if (temp <= node->value) {
      l += 1;

      if (temp > lbig) {
        lbig = temp;
        // If we know there are small units, we don't need to track the big
        rsmall = -DBL_MAX;
      }
    } else {
      r -= 1;
      std::swap(splitUnits[l], splitUnits[r]);

      if (temp < rsmall) {
        rsmall = temp;
        // If we know there are big units, we don't need to track the small
        lbig = DBL_MAX;
      }
    }
  }

  // If there exists units on both sides of the splitting value
  // we can be satisfied with the proposed split
  if (l > 0 && r < n)
    return l;

  // Now we have two cases: either
  // (1) all units are <= than the proposed splitting value, or
  // (2) all units are > than the proposed splitting value

  // We start with (2), we need to find all the smallest units, and move
  // the splitting value to these units
  if (l == 0) {
    for (size_t i = 0; i < n; i++) {
      double temp = *(dt + splitUnits[i] * p);
      if (temp == rsmall) {
        if (i != l)
          std::swap(splitUnits[i], splitUnits[l]);

        l += 1;
      }
    }

    // If it turns out that all units had the same value, maximum seperation has
    // been achieved
    if (l == n)
      return 0;

    // Otherwise, we can accept this candidate split as the splitting values
    node->value = rsmall;
    return l;
  }

  // We continue to the case (1), where we need to find all the biggest units
  // Here, we need to also find the next biggest value, as this will becom the
  // new splitting value.
  if (r == n) {
    rsmall = -DBL_MAX;

    for (size_t i = n; i-- > 0;) {
      double temp = *(dt + splitUnits[i] * p);
      if (temp == lbig) {
        r -= 1;

        if (i != r)
          std::swap(splitUnits[i], splitUnits[r]);
      } else {
        if (temp > rsmall)
          rsmall = temp;
      }
    }

    // If it turns out that all units had the same value, maximum seperation has
    // been achieved
    if (r == 0)
      return 0;

    // Otherwise, we can accept this candidate split as the splitting values
    node->value = rsmall;
    return r;
  }

  // If we managed to come here, something has went horribly worng
  throw std::runtime_error("(SplitByMidpointSlide) something went wrong in splitting");
  return 0;
}

size_t KDTree::SplitUnitsById(size_t* splitUnits, const size_t n, const size_t id, const size_t k) {
  size_t* tunits = new size_t[n];
  double* dt = data + k;
  size_t l = 0;
  size_t m = 0;
  size_t r = n;
  double value = *(dt + splitUnits[id] * p); // Proposed splitting value

  // Split units so that we get
  // x < value is in range [0, l)
  // x > value is in range [r, n)
  // x = value in tunits
  // where value is of the proposed index id
  for (size_t i = 0; i < r;) {
    double temp = *(dt + splitUnits[i] * p);
    if (temp < value) {
      if (i != l)
        splitUnits[l] = splitUnits[i];
      i += 1;
      l += 1;
    } else if (temp > value) {
      r -= 1;
      std::swap(splitUnits[i], splitUnits[r]);
    } else {
      tunits[m] = splitUnits[i];
      i += 1;
      m += 1;
    }
  }

  // Put back mid units into splitUnits
  for (size_t i = 0; i < m; i++)
    splitUnits[l + i] = tunits[i];

  delete[] tunits;

  // l + m = r
  // If the proposed id exists in [0, l), we need to sort this part further
  // If the proposed id exists in [r, n), we need to sort this part further
  // If the proposed id exists in [l, r), we are fine with choosing r
  // 0, ..., l, ..., r, ..., n
  if (id < l)
    return SplitUnitsById(splitUnits, l, id, k);
  if (id >= r)
    return r + SplitUnitsById(splitUnits + r, n - r, id - r, k);

  return r;
}


KDNode* KDTree::FindNode(const size_t id) {
  double* unit = data + id * p;
  KDNode* node = topNode;

  while (node != nullptr && !node->IsTerminal())
    node = unit[node->split] <= node->value ? node->cleft : node->cright;

  return node;
}

bool KDTree::UnitExists(const size_t id) {
  KDNode* node = FindNode(id);

  if (node == nullptr) {
    throw std::runtime_error("(UnitExists) node error");
    return false;
  }

  return node->UnitExists(id);
}

void KDTree::RemoveUnit(const size_t id) {
  KDNode* node = FindNode(id);

  if (node == nullptr) {
    throw std::runtime_error("(RemoveExists) node error");
    return;
  }

  node->RemoveUnit(id);
  return;
}

double KDTree::DistanceBetweenUnits(const size_t id1, const size_t id2) {
  return DistanceBetweenPointers(data + id1 * p, data + id2 * p);
}

double KDTree::DistanceBetweenPointers(const double* dt1, const double* dt2) {
  double distance = 0.0;

  for (size_t k = 0; k < p; k++) {
    double temp = dt1[k] - dt2[k];
    distance += temp * temp;
  }

  return distance;
}

void KDTree::FindNeighbours(KDStore* store, const size_t id) {
  store->Reset();

  if (topNode == nullptr) {
    throw std::runtime_error("(FindNeighbours) topNode is nullptr");
    return;
  }

  double* unit = data + id * p;

  TraverseNodesForNeighbours(store, id, unit, topNode);
  return;
}

void KDTree::FindNeighbours(KDStore* store, const double* unit) {
  store->Reset();

  if (topNode == nullptr) {
    throw std::runtime_error("(FindNeighbours) topNode is nullptr");
    return;
  }

  TraverseNodesForNeighbours(store, N + 1, unit, topNode);
  return;
}

void KDTree::TraverseNodesForNeighbours(
  KDStore* store,
  const size_t id,
  const double* unit,
  KDNode* node
) {
  if (node == nullptr) {
    throw std::runtime_error("(TraverseNodesForNeighbours) nullptr");
    return;
  }

  if (node->IsTerminal()) {
    if (store->maxSize == 1) {
      SearchNodeForNeighbour1(store, id, unit, node);
      return;
    }

    SearchNodeForNeighbours(store, id, unit, node);
    return;
  }

  double distance = unit[node->split] - node->value;
  KDNode* nextNode = distance <= 0.0 ? node->cleft : node->cright;

  TraverseNodesForNeighbours(store, id, unit, nextNode);

  // We only need to look at the wrong side of the tree if
  // (A) we have too few units
  // (B) the ball around the unit includes the other node
  if (!store->SizeFulfilled() || distance * distance <= store->MaximumDistance()) {
    TraverseNodesForNeighbours(store, id, unit, nextNode->GetSibling());
  }

  return;
}

void KDTree::SearchNodeForNeighbour1(
  KDStore* store,
  const size_t id,
  const double* unit,
  KDNode* node
) {
  size_t nodeSize = node->GetSize();
  double currentMinimum = store->MinimumDistance();

  for (size_t i = 0; i < nodeSize; i++) {
    size_t tid = node->units[i];
    // Skip if it is the same unit
    if (tid == id)
      continue;

    double distance = DistanceBetweenPointers(unit, data + tid * p);

    if (distance < currentMinimum) {
      store->AddUnitAndReset(tid);
      store->SetDistance(tid, distance);
      currentMinimum = distance;
    } else if (distance == currentMinimum) {
      store->AddUnit(tid);
      store->SetDistance(tid, distance);
    }
  }

  return;
}

void KDTree::SearchNodeForNeighbours(
  KDStore* store,
  const size_t id,
  const double* unit,
  KDNode* node
) {
  size_t nodeSize = node->GetSize();
  // Node is empty, we can skip
  if (nodeSize == 0)
    return;

  size_t originalSize = store->GetSize();
  bool originalFulfilled = store->SizeFulfilled();
  double currentMaximum = store->MaximumDistance();
  double nodeMinimum = DBL_MAX;
  // If we're full, set the nodeMax to currentMax, as we don't need to consider
  // units with larger distances. Otherwise, we set the nodeMax to 0.0
  double nodeMaximum = originalFulfilled ? currentMaximum : 0.0;

  // Search through all units in the node, and store the distances
  for (size_t i = 0; i < nodeSize; i++) {
    size_t tid = node->units[i];
    // Skip if it is the same unit
    if (tid == id)
      continue;

    double distance = DistanceBetweenPointers(unit, data + tid * p);

    // If we have a unit with distance larger than the nodeMax,
    // we continue if we're full,
    // if we're not full, we will add this unit and set the nodeMax to this new
    // distance.
    // If we were full before starting, we will just skip any units that are not
    // of interest. If we were not full before starting, we will add the first
    // units to fill up the size, after that we will only add units which are
    // smaller than the largest unit of these first added.
    if (distance > nodeMaximum) {
      if (store->SizeFulfilled())
        continue;
      else
        nodeMaximum = distance;
    }

    store->SetDistance(tid, distance);
    store->AddUnit(tid);

    if (distance < nodeMinimum)
      nodeMinimum = distance;
  }

  size_t storeSize = store->GetSize();

  // If we didn't add any units from this node, we have nothing to process
  if (storeSize == originalSize)
    return;

  // We now have two possibilities
  // - The node minimum is smaller than the current minimum
  // - The node minimum is larger than the current maximum
  // - The node minimum is somewhere inbetween
  size_t i;
  if (originalSize == 0 || nodeMinimum < store->MinimumDistance()) {
    i = 0;
  // } else if (nodeMinDistance >= currentMaxDistance) {
  //   i = originalSize;
  } else {
    for (i = originalSize; i-- > 0;) {
      if (nodeMinimum >= store->GetDistance(i))
        break;
    }

    i += 1;
  }

  // Sort the range [i, neighbours.size())
  store->SortNeighboursByDistance(i, storeSize);

  // When this loop breaks, we have an i that is at least as large as needed
  // 'i' will then be the smallest unit that is not to be included
  for(i += 1; i < storeSize; i++) {
    if (i < store->maxSize)
      continue;
    if (store->GetDistance(i - 1) < store->GetDistance(i))
      break;
  }

  store->neighbours.resize(i);
  return;
}

void KDTree::FindNeighboursCps(
  KDStore* store,
  const std::vector<double>& probabilities,
  const size_t id
) {
  store->Reset();

  if (topNode == nullptr) {
    throw std::runtime_error("(FindNeighbours) topNode is nullptr");
    return;
  }

  double* unit = data + id * p;
  double totalWeight = 0.0;

  TraverseNodesForNeighboursCps(store, probabilities, id, unit, topNode, &totalWeight);
  return;
}

void KDTree::TraverseNodesForNeighboursCps(
  KDStore* store,
  const std::vector<double>& probabilities,
  const size_t id,
  const double* unit,
  KDNode* node,
  double* totalWeight
) {
  if (node == nullptr) {
    throw std::runtime_error("(TraverseNodesForNeighbours) nullptr");
    return;
  }

  if (node->IsTerminal()) {
    SearchNodeForNeighboursCps(store, probabilities, id, unit, node, totalWeight);
    return;
  }

  double distance = unit[node->split] - node->value;
  KDNode* nextNode = distance <= 0.0 ? node->cleft : node->cright;

  TraverseNodesForNeighboursCps(store, probabilities, id, unit, nextNode, totalWeight);

  // We only need to look at the wrong side of the tree if
  // (A) we have too few units
  // (B) the ball around the unit includes the other node
  if (*totalWeight < 1.0 || distance * distance <= store->MaximumDistance()) {
    TraverseNodesForNeighboursCps(store, probabilities, id, unit, nextNode->GetSibling(), totalWeight);
  }

  return;
}

void KDTree::SearchNodeForNeighboursCps(
  KDStore* store,
  const std::vector<double>& probabilities,
  const size_t id,
  const double* unit,
  KDNode* node,
  double* totalWeight
) {
  size_t nodeSize = node->GetSize();
  // Node is empty, we can skip
  if (nodeSize == 0)
    return;

  size_t originalSize = store->GetSize();
  bool originalFulfilled = *totalWeight >= 1.0;
  double currentMaximum = store->MaximumDistance();
  double nodeMinimum = DBL_MAX;
  // If we're full, set the nodeMax to currentMax, as we don't need to consider
  // units with larger distances. Otherwise, we set the nodeMax to 0.0
  double nodeMaximum = originalFulfilled ? currentMaximum : 0.0;
  double nodeWeight = *totalWeight;

  // Search through all units in the node, and store the distances
  for (size_t i = 0; i < nodeSize; i++) {
    size_t tid = node->units[i];
    // Skip if it is the same unit
    if (tid == id)
      continue;

    double distance = DistanceBetweenPointers(unit, data + tid * p);

    // If we have a unit with distance larger than the nodeMax,
    // we continue if we're full,
    // if we're not full, we will add this unit and set the nodeMax to this new
    // distance.
    // If we were full before starting, we will just skip any units that are not
    // of interest. If we were not full before starting, we will add the first
    // units to fill up the size, after that we will only add units which are
    // smaller than the largest unit of these first added.
    if (distance > nodeMaximum) {
      if (nodeWeight >= 1.0)
        continue;
      else
        nodeMaximum = distance;
    }

    double weight = ((probabilities[id] + probabilities[tid]) <= 1.0)
      ? probabilities[tid] / (1.0 - probabilities[id])
      : (1.0 - probabilities[tid]) / probabilities[id];
    nodeWeight += weight;

    store->SetDistance(tid, distance);
    store->SetWeight(tid, weight);
    store->AddUnit(tid);


    if (distance < nodeMinimum)
      nodeMinimum = distance;
  }

  size_t storeSize = store->GetSize();

  // If we didn't add any units from this node, we have nothing to process
  if (storeSize == originalSize)
    return;

  // We now have two possibilities
  // - The node minimum is smaller than the current minimum
  // - The node minimum is larger than the current maximum
  // - The node minimum is somewhere inbetween
  size_t i;
  if (originalSize == 0 || nodeMinimum < store->GetDistance(0)) {
    i = 0;
    *totalWeight = 0.0;
  // } else if (nodeMinDistance >= currentMaxDistance) {
  //   i = originalSize;
  } else {
    for (i = originalSize; i-- > 0;) {
      if (nodeMinimum >= store->GetDistance(i))
        break;

      *totalWeight -= store->GetWeight(i);
    }

    i += 1;
    // "i" is now pointing to the biggest unit removed from totalWeight
    // Remember that the loop must have been break:ed, since otherwise
    // it would have chosen the other if-else path
  }

  // Sort the range [i, neighbours.size())
  store->SortNeighboursByDistance(i, storeSize);

  // When this loop breaks, we have an i that is at least as large as needed
  // 'i' will then be the smallest unit that is not to be included
  double prevDist = i == 0 ? -1.0 : store->GetDistance(i - 1);
  while(i < storeSize) {
    double thisDist = store->GetDistance(i);
    if (*totalWeight >= 1.0 && prevDist < thisDist)
      break;

    prevDist = thisDist;
    *totalWeight += store->GetWeight(i);
    i += 1;
  }

  store->neighbours.resize(i);
  return;
}
