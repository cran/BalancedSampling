#include <algorithm>
#include <float.h>
#include <stddef.h>
#include <stdexcept>
#include <utility>
#include <vector>

#include "CubeClass.h"
#include "CubeVectorInNullSpace.h"
#include "IndexListClass.h"
#include "KDStoreClass.h"
#include "KDTreeClass.h"
#include "ReducedRowEchelonForm.h"
#include "uniform.h"
#include "utils.h"
#include "utils-matrix.h"

//**********************************************
// Author: Wilmer Prentius
// Licence: GPL (>=2)
//**********************************************

CubeMethod IntToCubeMethod(const int i) {
  if (1 <= i && i <= 2)
    return static_cast<CubeMethod>(i);

  throw std::invalid_argument("cube-method does not exist");
  return CubeMethod::err;
}

// DIRECT CUBE
Cube::Cube(
  const double* t_probabilities,
  double* xxbalance,
  const size_t t_N,
  const size_t t_pbalance,
  const double t_eps
  ) {
  set_direct = true;
  cubeMethod = CubeMethod::cube;
  Init(t_probabilities, xxbalance, t_N, t_pbalance, t_eps);

  idx->Shuffle();
}

// DIRECT LCUBE
Cube::Cube(
  const double* t_probabilities,
  double* xxbalance,
  const size_t t_N,
  const size_t t_pbalance,
  const double t_eps,
  double* xxspread,
  const size_t t_pspread,
  const size_t t_treeBucketSize,
  const int t_treeMethod
) {
  set_direct = true;
  cubeMethod = CubeMethod::lcube;
  tree = new KDTree(xxspread, t_N, t_pspread, t_treeBucketSize, IntToKDTreeSplitMethod(t_treeMethod));

  Init(t_probabilities, xxbalance, t_N, t_pbalance, t_eps);
}

// INDIRECT CUBE
Cube::Cube(
  const CubeMethod t_cubeMethod,
  const size_t t_N,
  const size_t t_pbalance,
  const double t_eps
) {
  set_direct = false;
  cubeMethod = t_cubeMethod;
  InitIndirect(t_N, t_pbalance, t_eps);
}

void Cube::Init(
  const double* t_probabilities,
  double* xxbalance,
  const size_t t_N,
  const size_t t_pbalance,
  const double t_eps
) {
  InitIndirect(t_N, t_pbalance, t_eps);

  idx = new IndexList(N);

  for (size_t i = N; i-- > 0; ) {
    probabilities[i] = t_probabilities[i];
    idx->Set(i);

    if (ProbabilityInt(probabilities[i], eps)) {
      EraseUnit(i);

      if (Probability1(probabilities[i], eps))
        AddUnitToSample(i);

      continue;
    }

    for (size_t k = 0; k < pbalance; k++)
      amat[MatrixIdxCM(i, k, N)] =
        xxbalance[MatrixIdxCM(i, k, pbalance)] / probabilities[i];
  }

  return;
}

void Cube::InitIndirect(
  const size_t t_N,
  const size_t t_pbalance,
  const double t_eps
) {
  if (t_N == 0)
    return;

  N = t_N;
  pbalance = t_pbalance;
  eps = t_eps;

  probabilities.resize(N);
  sample.reserve(N);
  candidates.reserve(pbalance + 1);

  amat.resize(N * pbalance);
  uvec.resize(pbalance + 1);
  bmat.resize((pbalance + 1) * pbalance);

  switch(cubeMethod) {
  case CubeMethod::cube:
    _Draw = &Cube::Draw_cube;
    break;
  case CubeMethod::lcube:
    _Draw = &Cube::Draw_lcube;
    store = new KDStore(N, pbalance);
    break;
  default:
    throw std::invalid_argument("cubeMethod does not exist");
    break;
  }

  set_draw = true;
}

Cube::~Cube() {
  if (set_direct) {
    delete idx;
    delete tree;
  }

  delete store;
}

void Cube::AddUnitToSample(const size_t id) {
  sample.push_back(id + 1);
  return;
}

void Cube::EraseUnit(const size_t id) {
  idx->Erase(id);

  // Needed like this as tree might be nullptr during landing
  if (tree != nullptr)
    tree->RemoveUnit(id);

  return;
}

size_t Cube::MaxSize() {
  size_t il = idx->Length();
  return pbalance + 1 <= il ? pbalance + 1 : il;
}

void Cube::Draw_cube() {
  size_t maxSize = MaxSize();
  candidates.resize(0);

  for (size_t i = 0; i < maxSize; i++)
    candidates.push_back(idx->Get(i));

  return;
}

void Cube::Draw_lcube() {
  // maxSize - 1 since the first unit is drawn at random
  size_t maxSize = MaxSize() - 1;
  candidates.resize(1);

  // Set the first unit
  size_t id = idx->Draw();
  candidates[0] = id;

  // Prepare the store and run the algorithm
  store->maxSize = maxSize;
  tree->FindNeighbours(store, id);
  size_t size = store->GetSize();

  // If we have no equal dist units at the end, we can just return the candidates
  if (size == maxSize) {
    for (size_t i = 0; i < size; i++)
      candidates.push_back(store->neighbours[i]);
    return;
  }

  if (size < maxSize) {
    throw std::runtime_error("(Draw_lcube) size < maxSize - 1");
    return;
  }

  double maximumDistance = store->MaximumDistance();
  size_t i = 0;

  // Fill up all the units with distance less than maxdist
  // This algorithm will leave i to be one more than the size of candidates
  // when either maxSize has been filled (shouldn't happen),
  // or when the distance equals the maxdist
  for (; i < maxSize && store->GetDistance(i) < maximumDistance; i++)
    candidates.push_back(store->neighbours[i]);

  // Randomly select from the units at the end
  for (; i < maxSize; i++) {
    // Draw a random number from 0 to size - i.
    // If i = 2, size = 5, that means that we have already selected two units
    // and we have 3 more available in the store.
    size_t k = sizeuniform(size - i);
    candidates.push_back(store->neighbours[k + i]);
    if (k != 0)
      std::swap(store->neighbours[k + i], store->neighbours[i]);
  }

  return;
}

void Cube::Draw() {
  (this->*_Draw)();
  return;
}

void Cube::RunUpdate() {
  size_t maxSize = MaxSize();
  // bmat is transposed, and rref'ed
  ReducedRowEchelonForm(&bmat[0], maxSize - 1, maxSize);
  CubeVectorInNullSpace(&uvec[0], &bmat[0], maxSize);

  double lambda1 = DBL_MAX;
  double lambda2 = DBL_MAX;

  for (size_t i = 0; i < maxSize; i++) {
    double lval1 = std::abs(probabilities[candidates[i]] / uvec[i]);
    double lval2 = std::abs((1.0 - probabilities[candidates[i]]) / uvec[i]);

    if (uvec[i] >= 0.0) {
      if (lambda1 > lval2)
        lambda1 = lval2;
      if (lambda2 > lval1)
        lambda2 = lval1;
    } else {
      if (lambda1 > lval1)
        lambda1 = lval1;
      if (lambda2 > lval2)
        lambda2 = lval2;
    }
  }

  double lambda = stduniform(lambda1 + lambda2) < lambda2 ? lambda1 : -lambda2;

  for (size_t i = 0; i < maxSize; i++) {
    size_t id = candidates[i];
    probabilities[id] += lambda * uvec[i];

    if (ProbabilityInt(probabilities[id], eps)) {
      EraseUnit(id);

      if (Probability1(probabilities[id], eps))
        AddUnitToSample(id);
    }
  }

  return;
}

void Cube::RunFlight() {
  if (!set_draw)
    throw std::runtime_error("_Draw is nullptr");

  if (idx->Length() < pbalance + 1)
    return;

  // Cases:
  // - choose from tree and all
  // - choose from list and all

  size_t maxSize = MaxSize();

  while (idx->Length() >= maxSize) {
    Draw();

    // Prepare bmat
    // bmat is stored in Column Major order, but transposed and used in RM.
    for (size_t i = 0; i < maxSize; i++) {
      for (size_t k = 0; k < maxSize - 1; k++) {
        bmat[MatrixIdxCM(i, k, maxSize)] = amat[MatrixIdxCM(candidates[i], k, N)];
      }
    }

    RunUpdate();
  }

  return;
}

void Cube::RunLanding() {
  if (!set_draw)
    throw std::runtime_error("_Draw is nullptr");

  if (idx->Length() >= pbalance + 1)
    throw std::range_error("landingphase committed early");

  while (idx->Length() > 1) {
    size_t maxSize = idx->Length();
    candidates.resize(0);

    for (size_t i = 0; i < maxSize; i++) {
      size_t id = idx->Get(i);
      candidates.push_back(id);

      for (size_t k = 0; k < maxSize - 1; k++) {
        bmat[MatrixIdxCM(i, k, maxSize)] = amat[MatrixIdxCM(id, k, N)];
      }
    }

    RunUpdate();
  }

  if (idx->Length() == 1) {
    size_t id1 = idx->Get(0);

    if (stduniform() < probabilities[id1])
      AddUnitToSample(id1);

    EraseUnit(id1);
  }

  return;
}

void Cube::Run() {
  RunFlight();
  RunLanding();
  std::sort(sample.begin(), sample.end());

  return;
}
