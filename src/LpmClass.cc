#include <algorithm>
#include <stddef.h>
#include <stdexcept>
#include <utility>
#include <vector>

#include "IndexListClass.h"
#include "KDStoreClass.h"
#include "KDTreeClass.h"
#include "LpmClass.h"
#include "uniform.h"
#include "utils.h"

LpmMethod IntToLpmMethod(const int i) {
  if (1 <= i && i <= 5)
    return static_cast<LpmMethod>(i);

  throw std::invalid_argument("lpm-method does not exist");
  return LpmMethod::err;
}

// DIRECT
// -- DOUBLE
Lpm::Lpm(
  const LpmMethod t_lpMethod,
  const double* t_probabilities,
  double* xx,
  const size_t t_N,
  const size_t t_p,
  const double t_eps,
  const size_t t_treeBucketSize,
  const int t_treeMethod
) {
  set_direct = true;
  Init(xx, t_N, t_p, t_treeBucketSize, t_treeMethod, t_lpMethod);

  eps = t_eps;
  probabilities.resize(N);
  idx = new IndexList(N);

  for (size_t i = N; i-- > 0; ) {
    probabilities[i] = t_probabilities[i];
    idx->Set(i);

    if (ProbabilityInt(probabilities[i], eps)) {
      EraseUnit(i);

      if (Probability1(probabilities[i], eps))
        AddUnitToSample(i);
    }
  }

  _Run = &Lpm::Run_double;
  set_run = true;
}
// -- INT
Lpm::Lpm(
  const LpmMethod t_lpMethod,
  size_t t_probn,
  double* xx,
  const size_t t_N,
  const size_t t_p,
  const size_t t_treeBucketSize,
  const int t_treeMethod
) {
  set_direct = true;
  Init(xx, t_N, t_p, t_treeBucketSize, t_treeMethod, t_lpMethod);


  if (N == 0 || t_probn == 0) {
    idx = new IndexList(0);
  } else if (t_probn == N) {
    idx = new IndexList(0);
    for (size_t i = 0; i < N; i++)
      AddUnitToSample(i);
  } else {
    idx = new IndexList(N);
    idx->Fill();
    // iprobabilities.resize(0);
    iprobabilities.resize(N, t_probn);
  }

  _Run = &Lpm::Run_int;
  set_run = true;
}

void Lpm::Init(
  double* xx,
  const size_t t_N,
  const size_t t_p,
  const size_t t_bucketSize,
  const int t_method,
  const LpmMethod t_lpMethod
) {
  N = t_N;
  lpMethod = t_lpMethod;
  sample.reserve(N);

  switch(lpMethod) {
  case LpmMethod::lpm1:
    _Draw = &Lpm::Draw_lpm1;
    candidates.reserve(16);
    break;
  case LpmMethod::lpm2:
    _Draw = &Lpm::Draw_lpm2;
    break;
  case LpmMethod::lpm1search:
    _Draw = &Lpm::Draw_lpm1search;
    candidates.reserve(16);
    history.reserve(N);
    break;
  case LpmMethod::rpm:
    _Draw = &Lpm::Draw_rpm;
    break;
  case LpmMethod::spm:
    _Draw = &Lpm::Draw_spm;
    break;
  default:
    throw std::invalid_argument("(Lpm::Init) no such LpmMethod");
    return;
  }

  set_draw = true;

  if (xx != nullptr) {
    tree = new KDTree(xx, t_N, t_p, t_bucketSize, IntToKDTreeSplitMethod(t_method));
    store = new KDStore(N, 1);
  }

  return;
}

Lpm::~Lpm() {
  if (set_direct) {
    delete idx;
    delete tree;
    delete store;
  }
}

void Lpm::AddUnitToSample(const size_t id) {
  sample.push_back(id + 1);
  return;
}

void Lpm::EraseUnit(const size_t id) {
  idx->Erase(id);

  // Needed like this as tree might be nullptr during landing
  if (tree != nullptr)
    tree->RemoveUnit(id);

  return;
}

void Lpm::Draw_lpm1() {
  while (true) {
    pair[0] = idx->Draw();
    tree->FindNeighbours(store, pair[0]);
    size_t len = store->GetSize();

    candidates.reserve(len);
    candidates.resize(0);

    // We need to reuse the store to check for mutuality, so we put the results
    // in candidates
    for (size_t i = 0; i < len; i++)
      candidates.push_back(store->neighbours[i]);

    for (size_t i = 0; i < len;) {
      tree->FindNeighbours(store, candidates[i]);
      size_t tlen = store->GetSize();

      bool found = false;

      // Loop through, and set if pair[0] was found
      for (size_t j = 0; j < tlen; j++) {
        if (store->neighbours[j] == pair[0]) {
          found = true;
          break;
        }
      }

      // If we found something, we look through the next candidate.
      // Otherwise, we discard this candidate.
      if (found) {
        i += 1;
      } else {
        len -= 1;
        candidates[i] = candidates[len];
      }
    }

    // If there are still units in candidates, we select randomly amongst these
    // candidates, and return. Otherwise, we continue looking.
    if (len > 0) {
      pair[1] = candidates[sizeuniform(len)];
      return;
    }
  }
}

void Lpm::Draw_lpm2() {
  pair[0] = idx->Draw();
  tree->FindNeighbours(store, pair[0]);
  size_t k = sizeuniform(store->GetSize());
  pair[1] = store->neighbours[k];
  return;
}

void Lpm::Draw_lpm1search() {
  // Go back in the history and remove units that does not exist
  while (history.size() > 0) {
    if (idx->Exists(history.back()))
      break;

    history.pop_back();
  }

  // If there is no history, we draw a unit at random
  if (history.size() == 0)
    history.push_back(idx->Draw());

  while (true) {
    // Set the first unit to the last in history
    pair[0] = history.back();
    // Find this units nearest neighbours
    tree->FindNeighbours(store, pair[0]);
    size_t len = store->GetSize();

    candidates.reserve(len);
    candidates.resize(0);

    // We need to reuse the store to check for mutuality, so we put the results
    // in history
    for (size_t i = 0; i < len; i++)
      candidates.push_back(store->neighbours[i]);

    // Go through all nearest neighbours
    for (size_t i = 0; i < len;) {
      // Find the neighbours nearest neighbours
      tree->FindNeighbours(store, candidates[i]);
      size_t tlen = store->GetSize();

      bool found = false;

      // Check if any of these are the history-unit
      for (size_t j = 0; j < tlen; j++) {
        if (store->neighbours[j] == pair[0]) {
          found = true;
          break;
        }
      }

      // If the history-unit exists among the nearest neighbours, we continue
      // to see if any other of the history-units neighbours also are mutual.
      // Otherwise, the history-unit is not among the nearest neighbours,
      // we swap places and continue the search.
      if (found) {
        i += 1;
      } else {
        len -= 1;
        if (i != len)
          std::swap(candidates[i], candidates[len]);
      }
    }

    // If we found one or more mutual neighbours, we select one at random
    if (len > 0) {
      pair[1] = candidates[sizeuniform(len)];
      return;
    }

    // If we come here, no mutual neighbours exist

    // We might need to clear the history if the search has been going on for
    // too long. This can probably? happen if there is a long history, and
    // updates has affected previous units.
    if (history.size() == N) {
      history.resize(0);
      history.push_back(pair[0]);
    }

    // We select a unit at random to become the next history unit, and traverse
    // one step further.
    size_t k = sizeuniform(candidates.size());
    history.push_back(candidates[k]);
  }
}

void Lpm::Draw_rpm() {
  // Draw a unit i at random from 0 to N
  // Draw a unit j at random from 0 to N-1
  // If unit i == j, then set j = N
  pair[0] = idx->Draw();
  size_t len = idx->Length() - 1;
  pair[1] = idx->Get(sizeuniform(len));

  if (pair[0] == pair[1])
    pair[1] = idx->Get(len);

  return;
}

void Lpm::Draw_spm() {
  // If the first unit does no longer exist,
  // set the first unit to be the second unit.
  if (!idx->Exists(pair[0])) {
    pair[0] = pair[1];

    // If the second unit also doesn't exist, increment by 1
    while (!idx->Exists(pair[0])) {
      pair[0] += 1;

      if (pair[0] >= N)
        throw std::range_error("invalid value of pair 0");
    }

    pair[1] = pair[0] + 1;
  }

  while (!idx->Exists(pair[1])) {
    pair[1] += 1;

    if (pair[1] >= N)
      throw std::range_error("invalid value of pair 1");
  }

  return;
}

void Lpm::Run_double() {
  while (idx->Length() > 1) {
    Draw();
    size_t id1 = pair[0];
    size_t id2 = pair[1];

    double* p1 = &probabilities[id1];
    double* p2 = &probabilities[id2];
    double psum = *p1 + *p2;

    if (psum > 1.0) {
      if (1.0 - *p2 > stduniform(2.0 - psum)) {
        *p1 = 1.0;
        *p2 = psum - 1.0;
      } else {
        *p1 = psum - 1.0;
        *p2 = 1.0;
      }
    } else {
      if (*p2 > stduniform(psum)) {
        *p1 = 0.0;
        *p2 = psum;
      } else {
        *p1 = psum;
        *p2 = 0.0;
      }
    }

    if (ProbabilityInt(*p1, eps)) {
      EraseUnit(id1);

      if (Probability1(*p1, eps))
        AddUnitToSample(id1);
    }

    if (ProbabilityInt(*p2, eps)) {
      EraseUnit(id2);

      if (Probability1(*p2, eps))
        AddUnitToSample(id2);
    }
  }

  if (idx->Length() == 1) {
    size_t id1 = idx->Get(0);

    if (stduniform() < probabilities[id1])
      AddUnitToSample(id1);

    EraseUnit(id1);
  }

  return;
}

void Lpm::Run_int() {
  while (idx->Length() > 1) {
    Draw();
    size_t id1 = pair[0];
    size_t id2 = pair[1];

    size_t* p1 = &iprobabilities[id1];
    size_t* p2 = &iprobabilities[id2];
    size_t psum = *p1 + *p2;

    if (psum > N) {
      if (N - *p2 > sizeuniform(2 * N - psum)) {
        *p1 = N;
        *p2 = psum - N;
      } else {
        *p1 = psum - N;
        *p2 = N;
      }
    } else {
      if (*p2 > sizeuniform(psum)) {
        *p1 = 0;
        *p2 = psum;
      } else {
        *p1 = psum;
        *p2 = 0;
      }
    }

    if (ProbabilityInt(*p1, N)) {
      EraseUnit(id1);

      if (Probability1(*p1, N))
        AddUnitToSample(id1);
    }

    if (ProbabilityInt(*p2, N)) {
      EraseUnit(id2);

      if (Probability1(*p2, N))
        AddUnitToSample(id2);
    }
  }

  if (idx->Length() == 1) {
    size_t id1 = idx->Get(0);

    if (sizeuniform(N) < iprobabilities[id1])
      AddUnitToSample(id1);

    EraseUnit(id1);
  }

  return;
}

void Lpm::Draw() {
  (this->*_Draw)();
  return;
}

void Lpm::Run() {
  if (!set_run)
    throw std::runtime_error("_run is nullptr");
  if (!set_draw)
    throw std::runtime_error("_draw is nullptr");

  (this->*_Run)();

  std::sort(sample.begin(), sample.end());
  return;
}
