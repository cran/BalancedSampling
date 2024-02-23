#include <algorithm>
#include <float.h>
#include <stdexcept>
#include <vector>
#include <stddef.h>

#include "KDStoreClass.h"

KDStore::KDStore() {}
KDStore::KDStore(const size_t t_N, const size_t t_maxSize) {
  Set(t_N, t_maxSize);
}
KDStore::~KDStore() {}

void KDStore::Set(const size_t t_N, const size_t t_maxSize) {
  if (t_maxSize == 0) {
    throw std::range_error("(Set) size must be > 0");
    return;
  }

  if (t_N == 0) {
    throw std::range_error("(Set) N must be > 0");
    return;
  }

  N = t_N;
  maxSize = t_maxSize;

  neighbours.reserve(N);
  distances.resize(N);

  Reset();
  return;
}

void KDStore::PrepareWeights() {
  weights.resize(N);
  return;
}

void KDStore::Reset() {
  neighbours.resize(0);
  return;
}

size_t KDStore::GetSize() {
  return neighbours.size();
}

bool KDStore::SizeFulfilled() {
  return neighbours.size() >= maxSize;
}

double KDStore::MinimumDistance() {
  if (neighbours.size() == 0)
    return DBL_MAX;

  return distances[neighbours[0]];
}

double KDStore::MaximumDistance() {
  if (neighbours.size() == 0)
    return DBL_MAX;

  return distances[neighbours.back()];
}

double KDStore::GetDistance(const size_t i) {
  return distances[neighbours.at(i)];
}

void KDStore::SetDistance(const size_t id, const double distance) {
  distances[id] = distance;
  return;
}

double KDStore::GetWeight(const size_t i) {
  return weights[neighbours.at(i)];
}

void KDStore::SetWeight(const size_t id, const double weight) {
  weights[id] = weight;
  return;
}

void KDStore::AddUnit(const size_t id) {
  neighbours.push_back(id);
}

void KDStore::AddUnitAndReset(const size_t id) {
  neighbours.resize(1);
  neighbours[0] = id;
}

void KDStore::SortNeighboursByDistance(const size_t from, const size_t to) {
  size_t size = GetSize();
  if (to <= from || size < to) {
    throw std::range_error("(SortNeighboursByDistance) bad input");
    return;
  }

  size_t* tn = neighbours.data();
  std::sort(
    tn + from,
    tn + to,
    [this](size_t a, size_t b) { return distances[a] < distances[b]; }
  );

  return;
}

void KDStore::SortNeighboursByWeight(const size_t from, const size_t to) {
  size_t size = GetSize();
  if (to <= from || size < to) {
    throw std::range_error("(SortNeighboursByDistance) bad input");
    return;
  }

  size_t* tn = neighbours.data();
  std::sort(
    tn + from,
    tn + to,
    [this](size_t a, size_t b) { return weights[a] < weights[b]; }
  );

  return;
}

