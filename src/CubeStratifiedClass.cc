#include <algorithm>
#include <stddef.h>
#include <unordered_map>
#include <vector>

#include "CubeClass.h"
#include "CubeStratifiedClass.h"
#include "IndexListClass.h"
#include "KDTreeClass.h"
#include "utils.h"

// CUBE
CubeStratified::CubeStratified(
  int* t_strata,
  double* t_probabilities,
  double* t_xbalance,
  const size_t t_N,
  const size_t t_pbalance,
  const double t_eps
) {
  cubeMethod = CubeMethod::cube;

  Init(
    t_strata,
    t_probabilities,
    t_xbalance,
    t_N,
    t_pbalance,
    t_eps
  );
}

// LCUBE
CubeStratified::CubeStratified(
  int* t_strata,
  double* t_probabilities,
  double* t_xbalance,
  const size_t t_N,
  const size_t t_pbalance,
  const double t_eps,
  double* t_xspread,
  const size_t t_pspread,
  const size_t t_treeBucketSize,
  const int t_treeMethod
) {
  cubeMethod = CubeMethod::lcube;
  r_xspread = t_xspread;
  pspread = t_pspread;
  treeBucketSize = t_treeBucketSize;
  treeMethod = IntToKDTreeSplitMethod(t_treeMethod);

  Init(
    t_strata,
    t_probabilities,
    t_xbalance,
    t_N,
    t_pbalance,
    t_eps
  );
}

void CubeStratified::Init(
  int* t_strata,
  double* t_probabilities,
  double* t_xbalance,
  const size_t t_N,
  const size_t t_pbalance,
  const double t_eps
) {
  r_strata = t_strata;
  N = t_N;
  eps = t_eps;

  r_prob = t_probabilities;
  r_xbalance = t_xbalance;
  pbalance = t_pbalance;

  idx = new IndexList(N);

  probabilities.resize(N);
  sample.reserve(N);

  for (size_t i = N; i-- > 0; ) {
    idx->Set(i);

    if (ProbabilityInt(r_prob[i], eps)) {
      idx->Erase(i);

      if (Probability1(r_prob[i], eps))
        AddUnitToSample(i);

      continue;
    }

    if (stratumMap.count(r_strata[i]) == 0)
      stratumMap[r_strata[i]] = 1;
    else
      stratumMap[r_strata[i]] += 1;
  }

  stratumArr.resize(stratumMap.size());

  size_t stratumMax = 0;
  for (std::unordered_map<int, size_t>::iterator it = stratumMap.begin(); it != stratumMap.end(); ++it) {
    if (it->second > stratumMax)
      stratumMax = it->second;
  }

  size_t fullsizeMax = (pbalance + 1) * stratumMap.size();
  if (stratumMax < fullsizeMax)
    stratumMax = fullsizeMax;

  index.reserve(stratumMax);

  if (cubeMethod == CubeMethod::lcube)
    zspread.reserve(stratumMax * pspread);

  return;
}

CubeStratified::~CubeStratified() {
  delete idx;
}

void CubeStratified::AddUnitToSample(const size_t id) {
  sample.push_back(id + 1);
  return;
}

void CubeStratified::RunFlightPerStratum() {
  size_t maxSize = pbalance + 2;
  stratumArr.resize(0);

  for (std::unordered_map<int, size_t>::iterator it = stratumMap.begin(); it != stratumMap.end(); ++it) {
    // We can skip this completely if there are too few obs in the stratum
    if (it->second < maxSize)
      continue;

    size_t idxlen = idx->Length();
    index.resize(0);
    zspread.resize(0);
    Cube cube(cubeMethod, it->second, maxSize - 1, eps);
    IndexList* tidx = new IndexList(it->second);
    cube.idx = tidx;

    for (size_t i = 0; i < idxlen; i++) {
      size_t id = idx->Get(i);
      if (it->first != r_strata[id])
        continue;

      size_t indexSize = index.size();
      index.push_back(id);
      tidx->Set(indexSize);

      cube.probabilities[indexSize] = r_prob[id];

      cube.amat[MatrixIdxRow((size_t)0, indexSize, it->second)] = 1.0;

      for (size_t k = 0; k < pbalance; k++)
        cube.amat[MatrixIdxRow(k + 1, indexSize, it->second)] =
          r_xbalance[MatrixIdxCol(id, k, N)] / r_prob[id];

      if (cubeMethod == CubeMethod::lcube) {
        for (size_t k = 0; k < pspread; k++)
          zspread.push_back(r_xspread[MatrixIdxRow(id, k, pspread)]);
      }
    }

    if (cubeMethod == CubeMethod::lcube) {
      KDTree* tree = new KDTree(zspread.data(), it->second, pspread, treeBucketSize, treeMethod);
      cube.tree = tree;
      cube.RunFlight();
      cube.tree = nullptr;
      delete tree;
    } else {
      cube.RunFlight();
    }

    for (size_t i = 0; i < index.size(); i++) {
      size_t id = index[i];

      if (!tidx->Exists(i)) {
        idx->Erase(id);
        it->second -= 1;
      } else {
        probabilities[id] = cube.probabilities[i];
      }
    }

    for (size_t i = 0; i < cube.sample.size(); i++)
      AddUnitToSample(index[cube.sample[i] - 1]);

    if (it->second == 0)
      stratumArr.push_back(it->first);

    cube.idx = nullptr;
    delete tidx;
  }

  for (size_t k = 0; k < stratumArr.size(); k++)
    stratumMap.erase(stratumArr[k]);

  // Repopulate stratumArr with stratum numbers
  stratumArr.resize(0);
  for (std::unordered_map<int, size_t>::iterator it = stratumMap.begin(); it != stratumMap.end(); ++it)
    stratumArr.push_back(it->first);

  return;
}

void CubeStratified::RunFlightOnFull() {
  size_t strsize = stratumMap.size();
  size_t maxSize = pbalance + 1 + strsize;
  size_t idxlen = idx->Length();

  // If we can't run, we can't run
  if (idxlen < maxSize)
    return;

  index.resize(0);
  zspread.resize(0);
  Cube cube(cubeMethod, idxlen, maxSize - 1, eps);
  IndexList* tidx = new IndexList(idxlen);
  cube.idx = tidx;

  for (size_t i = 0; i < idxlen; i++) {
    size_t id = idx->Get(i);
    size_t indexSize = index.size();
    index.push_back(id);
    tidx->Set(indexSize);
    cube.probabilities[indexSize] = probabilities[id];

    for (size_t k = 0; k < strsize; k++)
      cube.amat[MatrixIdxRow(k, indexSize, idxlen)] = (r_strata[id] == stratumArr[k]) ? 1.0 : 0.0;

    for (size_t k = 0; k < pbalance; k++)
      cube.amat[MatrixIdxRow(strsize + k, indexSize, idxlen)]
        = r_xbalance[MatrixIdxCol(id, k, N)] / r_prob[id];

    if (cubeMethod == CubeMethod::lcube) {
      for (size_t k = 0; k < pspread; k++)
        zspread.push_back(r_xspread[MatrixIdxRow(id, k, pspread)]);
    }
  }

  if (cubeMethod == CubeMethod::lcube) {
    KDTree* tree = new KDTree(zspread.data(), idxlen, pspread, treeBucketSize, treeMethod);
    cube.tree = tree;
    cube.RunFlight();
    cube.tree = nullptr;
    delete tree;
  } else {
    cube.RunFlight();
  }

  for (size_t i = 0; i < index.size(); i++) {
    size_t id = index[i];

    if (!tidx->Exists(i)) {
      idx->Erase(id);
      stratumMap[r_strata[id]] -= 1;
    } else {
      probabilities[id] = cube.probabilities[i];
    }
  }

  for (size_t i = 0; i < cube.sample.size(); i++)
    AddUnitToSample(index[cube.sample[i] - 1]);

  cube.idx = nullptr;
  delete tidx;

  return;
}

void CubeStratified::RunLandingPerStratum() {
  size_t maxSize = pbalance + 2;

  for (std::unordered_map<int, size_t>::iterator it = stratumMap.begin(); it != stratumMap.end(); ++it) {
    // Skip if all decided
    if (it->second <= 0)
      continue;

    size_t idxlen = idx->Length();
    index.resize(0);
    Cube cube(CubeMethod::cube, it->second, maxSize - 1, eps);
    IndexList* tidx = new IndexList(it->second);
    cube.idx = tidx;

    for (size_t i = 0; i < idxlen; i++) {
      size_t id = idx->Get(i);
      if (it->first != r_strata[id])
        continue;

      size_t indexSize = index.size();
      index.push_back(id);
      tidx->Set(indexSize);
      cube.probabilities[indexSize] = probabilities[id];

      cube.amat[MatrixIdxRow((size_t)0, indexSize, it->second)] = 1.0;

      for (size_t k = 0; k < pbalance; k++)
        cube.amat[MatrixIdxRow(k + 1, indexSize, it->second)] =
          r_xbalance[MatrixIdxCol(id, k, N)] / probabilities[id];
    }

    cube.RunLanding();

    for (size_t i = 0; i < cube.sample.size(); i++)
      AddUnitToSample(index[cube.sample[i] - 1]);

    cube.idx = nullptr;
    delete tidx;
  }

  return;
}

void CubeStratified::Run() {
  RunFlightPerStratum();
  RunFlightOnFull();
  RunLandingPerStratum();

  std::sort(sample.begin(), sample.end());

  return;
}
