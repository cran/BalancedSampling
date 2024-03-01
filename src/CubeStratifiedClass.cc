#include <algorithm>
#include <stddef.h>
#include <unordered_map>

#include "CubeClass.h"
#include "CubeStratifiedClass.h"
#include "IndexListClass.h"
#include "KDTreeClass.h"
#include "utils.h"
#include "utils-matrix.h"

// CUBE
CubeStratified::CubeStratified(
  int* strata,
  double* probabilities,
  double* xbalance,
  const size_t N,
  const size_t pbalance,
  const double eps
) {
  cube_method_ = CubeMethod::cube;
  cube_ = new Cube(cube_method_, N, pbalance + 1, eps);

  Init(
    strata,
    probabilities,
    xbalance,
    N,
    pbalance,
    eps
  );
}

// LCUBE
CubeStratified::CubeStratified(
  int* strata,
  double* probabilities,
  double* xbalance,
  const size_t N,
  const size_t pbalance,
  const double eps,
  double* xspread,
  const size_t pspread,
  const size_t treeBucketSize,
  const int treeMethod
) {
  cube_method_ = CubeMethod::lcube;
  cube_ = new Cube(cube_method_, N, pbalance + 1, eps);

  Init(
    strata,
    probabilities,
    xbalance,
    N,
    pbalance,
    eps
  );

  rptr_xspread_ = xspread;
  p_spread_ = pspread;
  tree_bucket_size_ = treeBucketSize;
  tree_method_ = IntToKDTreeSplitMethod(treeMethod);
}

CubeStratified::~CubeStratified() {
  delete idx_;
  delete cube_;
}

void CubeStratified::Init(
  int* strata,
  double* probabilities,
  double* xbalance,
  const size_t N,
  const size_t pbalance,
  const double eps
) {
  rptr_strata_ = strata;
  rptr_probabilities_ = probabilities;
  rptr_xbalance_ = xbalance;

  eps_ = eps;
  N_ = N;
  p_balance_ = pbalance;

  idx_ = new IndexList(N);

  for (size_t i = N; i-- > 0; ) {
    idx_->Set(i);

    if (ProbabilityInt(rptr_probabilities_[i], eps_)) {
      idx_->Erase(i);

      if (Probability1(rptr_probabilities_[i], eps_))
        cube_->AddUnitToSample(i);

      continue;
    }

    if (stratum_map_.count(rptr_strata_[i]) == 0)
      stratum_map_[rptr_strata_[i]] = 1;
    else
      stratum_map_[rptr_strata_[i]] += 1;

    cube_->probabilities[i] = rptr_probabilities_[i];
  }

  stratum_arr_.resize(stratum_map_.size());

  size_t stratum_max = 0;
  for (StratumMap::iterator it = stratum_map_.begin(); it != stratum_map_.end(); ++it) {
    if (it->second > stratum_max)
      stratum_max = it->second;
  }

  size_t fullsize_max = (pbalance + 1) * stratum_map_.size();
  if (stratum_max < fullsize_max)
    stratum_max = fullsize_max;

  return;
}

void CubeStratified::PrepareAmat(const size_t id) {
  cube_->amat[MatrixIdxCM(id, (size_t)0, N_)] = 1.0;
  PrepareAmatAux(id, 1);
}

void CubeStratified::PrepareAmatAux(const size_t id, const size_t offset) {
  for (size_t k = 0; k < p_balance_; k++)
    cube_->amat[MatrixIdxCM(id, k + offset, N_)] =
      rptr_xbalance_[MatrixIdxCM(id, k, N_)] / rptr_probabilities_[id];
}

void CubeStratified::RunFlightPerStratum() {
  size_t max_competitors = p_balance_ + 2; // p + pi + 1
  stratum_arr_.resize(0);

  for (StratumMap::iterator it = stratum_map_.begin(); it != stratum_map_.end();) {
    // We can skip this completely if there are too few obs in the stratum
    if (it->second < max_competitors)
      continue;

    size_t remaining_units = idx_->Length();
    IndexList* tidx = idx_->CopyLen();
    cube_->idx = tidx;

    for (size_t i = remaining_units; i --> 0; ) {
      size_t id = idx_->Get(i);

      if (it->first == rptr_strata_[id]) {
        idx_->Erase(id);
      } else {
        tidx->Erase(id);
        continue;
      }

      PrepareAmat(id);
    }

    if (cube_method_ == CubeMethod::lcube) {
      size_t* list = tidx->CopyList();
      cube_->tree = new KDTree(rptr_xspread_, N_, p_spread_, tree_bucket_size_, tree_method_, list, tidx->Length());
      cube_->RunFlight();
      delete cube_->tree;
      cube_->tree = nullptr;
      delete[] list;
    } else {
      cube_->RunFlight();
    }

    if (tidx->Length() == 0) {
      // Remove stratum if it no longer has any units.
      it = stratum_map_.erase(it);
    } else {
      it->second = tidx->Length();
      stratum_arr_.push_back(it->first);

      // Add back any undecided units, i.e. units still in the tidx list.
      for (size_t i = 0, tn = tidx->Length(); i < tn; i++) {
        idx_->Add(tidx->Get(i));
      }

      ++it;
    }

    cube_->idx = nullptr;
    delete tidx;
  }
}

void CubeStratified::RunFlightOnFull() {
  size_t stratum_size = stratum_map_.size();
  size_t max_competitors = p_balance_ + 1 + stratum_size; // p + 1 + pi per stratum
  cube_->idx = idx_;

  // If we don't have enough units to run a flight step, skip it
  if (idx_->Length() < max_competitors)
    return;

  delete cube_->store;
  cube_->InitIndirect(N_, max_competitors - 1, eps_);

  for (size_t i = 0; i < idx_->Length(); i++) {
    size_t id = idx_->Get(i);

    for (size_t k = 0; k < stratum_size; k++)
      cube_->amat[MatrixIdxCM(id, k, N_)] = (rptr_strata_[id] == stratum_arr_[k]) ? 1.0 : 0.0;

    PrepareAmatAux(id, stratum_size);
  }

  if (cube_method_ == CubeMethod::lcube) {
    size_t* list = idx_->CopyList();
    cube_->tree = new KDTree(rptr_xspread_, N_, p_spread_, tree_bucket_size_, tree_method_, list, idx_->Length());
    cube_->RunFlight();
    delete cube_->tree;
    cube_->tree = nullptr;
    delete[] list;
  } else {
    cube_->RunFlight();
  }

  for (StratumMap::iterator it = stratum_map_.begin(); it != stratum_map_.end(); ++it)
    it->second = 0;

  for (size_t i = 0, n = idx_->Length(); i < n; i++)
    stratum_map_[rptr_strata_[idx_->Get(i)]] += 1;
}

void CubeStratified::RunLandingPerStratum() {
  size_t max_competitors = p_balance_ + 2; // p + 1 + pi
  delete cube_->store;
  cube_->InitIndirect(N_, max_competitors - 1, eps_);

  for (StratumMap::iterator it = stratum_map_.begin(); it != stratum_map_.end(); ++it) {
    // Skip if all decided
    if (it->second <= 0)
      continue;

    size_t remaining_units = idx_->Length();
    IndexList* tidx = idx_->CopyLen();
    cube_->idx = tidx;

    for (size_t i = remaining_units; i --> 0; ) {
      size_t id = idx_->Get(i);

      if (it->first == rptr_strata_[id]) {
        idx_->Erase(id);
      } else {
        tidx->Erase(id);
        continue;
      }

      PrepareAmat(id);
    }

    cube_->RunLanding();
    cube_->idx = nullptr;
    delete tidx;
  }
}

void CubeStratified::Run() {
  RunFlightPerStratum();
  RunFlightOnFull();
  RunLandingPerStratum();

  sample_ = cube_->sample;
  std::sort(sample_.begin(), sample_.end());

  return;
}

