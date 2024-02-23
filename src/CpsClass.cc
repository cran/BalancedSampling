#include <algorithm>
#include <float.h>
#include <stddef.h>
#include <stdexcept>
// #include <utility>
#include <vector>

#include "CpsClass.h"
#include "IndexListClass.h"
#include "KDStoreClass.h"
#include "KDTreeClass.h"
#include "uniform.h"
#include "utils.h"

CpsMethod IntToCpsMethod(const int i) {
  if (1 <= i && i <= 3)
    return static_cast<CpsMethod>(i);

  throw std::invalid_argument("cps-method does not exist");
  return CpsMethod::err;
}

Cps::Cps(
    const CpsMethod t_cpsMethod,
    const double* t_probabilities,
    double* xx,
    const size_t t_N,
    const size_t t_p,
    const double t_eps,
    const size_t t_treeBucketSize,
    const int t_treeMethod
) {
  set_direct = true;
  cpsMethod = t_cpsMethod;
  N = t_N;
  eps = t_eps;

  if (xx == nullptr) {
    throw std::invalid_argument("(Cps) x is nullptr");
    return;
  }

  idx = new IndexList(N);
  tree = new KDTree(xx, N, t_p, t_treeBucketSize, IntToKDTreeSplitMethod(t_treeMethod));
  store = new KDStore(N, 1);
  store->PrepareWeights();

  probabilities.resize(N);
  sample.reserve(N);

  switch(cpsMethod) {
  case CpsMethod::lcps:
    _Draw = &Cps::Draw_lcps;
    // candidates.reserve(32);
    candidates.reserve(N);
    break;
  case CpsMethod::scps:
    _Draw = &Cps::Draw_scps;
    break;
  case CpsMethod::scpscoord:
    _Draw = &Cps::Draw_scpscoord;
    break;
  default:
    throw std::invalid_argument("(Lpm::Init) no such LpmMethod");
    return;
  }

  set_draw = true;

  for (size_t i = N; i-- > 0; ) {
    probabilities[i] = t_probabilities[i];
    idx->Set(i);

    if (ProbabilityInt(probabilities[i], eps)) {
      EraseUnit(i);

      if (Probability1(probabilities[i], eps))
        AddUnitToSample(i);
    }
  }

  SetRandomStd();
}

Cps::~Cps() {
  if (set_direct) {
    delete idx;
    delete tree;
    delete store;
  }
}

void Cps::SetRandomStd() {
  _Random = &Cps::Random_std;
  set_random = true;
}

void Cps::SetRandomArr(double* t_rand) {
  randomValues = t_rand;
  _Random = &Cps::Random_arr;
  set_random = true;
}

void Cps::AddUnitToSample(const size_t id) {
  sample.push_back(id + 1);
  return;
}

void Cps::EraseUnit(const size_t id) {
  idx->Erase(id);

  // Needed like this as tree might be nullptr during landing
  if (tree != nullptr)
    tree->RemoveUnit(id);

  return;
}

void Cps::DecideUnit(const size_t id) {
  if (ProbabilityInt(probabilities[id], eps)) {
    EraseUnit(id);

    if (Probability1(probabilities[id], eps))
      AddUnitToSample(id);
  }

  return;
}

size_t Cps::Draw_lcps() {
  // Take care of edge cases
  if (idx->Length() <= 1) {
    if (idx->Length() == 1)
      return idx->Get(0);
    if (idx->Length() == 0)
      throw std::range_error("trying to find index in empty list");
  }

  double mindist = DBL_MAX;
  candidates.resize(0);

  // Loop through all remaining units.
  // Put the smallest distances in candidates
  for (size_t i = 0; i < idx->Length(); i++) {
    size_t id = idx->Get(i);
    tree->FindNeighboursCps(store, probabilities, id);
    double dist = store->MaximumDistance();

    if (dist < mindist) {
      candidates.resize(1);
      candidates[0] = id;
      mindist = dist;
    } else if (dist == mindist) {
      candidates.push_back(id);
    }
  }

  // Choose randomly from the units in candidates
  size_t k = sizeuniform(candidates.size());
  return candidates[k];
}

size_t Cps::Draw_scps() {
  return idx->Draw();
}

size_t Cps::Draw_scpscoord() {
  while(!idx->Exists(candidateIdx))
    candidateIdx += 1;

  size_t tunit = candidateIdx;
  candidateIdx += 1;
  return tunit;
}

double Cps::Random_std(const size_t id) {
  return stduniform();
}

double Cps::Random_arr(const size_t id) {
  return randomValues[id];
}

size_t Cps::Draw() {
  return (this->*_Draw)();
}

double Cps::Random(const size_t id) {
  return (this->*_Random)(id);
}

void Cps::Run() {
  if (!set_draw)
    throw std::runtime_error("_Draw is nullptr");
  if (!set_random)
    throw std::runtime_error("_Random is nullptr");

  while (idx->Length() > 1) {
    size_t id1 = Draw();

    // We need to remove the unit first, so that it is not searching itself
    // in the tree search
    EraseUnit(id1);

    // Find all neighbours
    tree->FindNeighboursCps(store, probabilities, id1);
    size_t len = store->GetSize();

    double slag = probabilities[id1];

    if (Random(id1) < probabilities[id1]) {
      slag -= 1.0;
      AddUnitToSample(id1);
      probabilities[id1] = 1.0;
    } else {
      probabilities[id1] = 0.0;
    }

    // The weight that remains to be put out to the neighbours
    double remweight = 1.0;

    // Loop through all found neighbours
    // The loop is conducted so that we take equal distance neighbours together
    for (size_t i = 0; i < len && remweight > eps;) {
      // First we need to find how many neighbours exists on the same distance
      // Initialize totweight to the first neighbour, then search through
      // until the distance differs from this first neighbour
      double totweight = store->GetWeight(i);

      size_t j = i + 1;
      for (; j < len; j++) {
        if (store->GetDistance(i) < store->GetDistance(j))
          break;

        totweight += store->GetWeight(j);
      }


      // If we only found one nearest neighbour, we resolve this and continue
      if (j - i == 1) {
        size_t id2 = store->neighbours[i];

        // Do not use more than the remaining weight
        double temp = remweight >= totweight ? totweight : remweight;

        probabilities[id2] += temp * slag;
        DecideUnit(id2);

        i += 1;
        remweight -= temp;
        continue;
      }

      // If we found multiple nearest neighbours
      if (remweight >= totweight) {
        // The remaining weight is larger than the total weight of the nearest neighbours
        // Loop through the nearest neighbours and update their probabilities
        for (; i < j; i++) {
          size_t id2 = store->neighbours[i];
          probabilities[id2] += store->weights[id2] * slag;
          DecideUnit(id2);
        }

        remweight -= totweight;
      } else {
        // The remaining weight is smaller than the total weight of the nearest neighbours
        // We need to sort this list, smallest weights first
        store->SortNeighboursByWeight(i, j);

        // Loop through all units, and update their weights
        // No unit can get more than a fair share
        for (; i < j; i++) {
          size_t id2 = store->neighbours[i];
          // Temp contains fair share
          double temp = remweight / (double)(j - i);
          // But we cannot update with more than the assigned weight
          if (store->weights[id2] < temp)
            temp = store->weights[id2];

          probabilities[id2] += temp * slag;
          DecideUnit(id2);
          remweight -= temp;
        }
      }
    }
  }

  // Sort out any remaining lone unit
  if (idx->Length() == 1) {
    size_t id1 = idx->Get(0);

    if (Random(id1) < probabilities[id1])
      AddUnitToSample(id1);

    EraseUnit(id1);
  }

  std::sort(sample.begin(), sample.end());
  return;
}
