#ifndef KDSTORECLASS_HEADER
#define KDSTORECLASS_HEADER

#include <vector>
#include <stddef.h>

class KDStore {
public:
  size_t N;
  size_t maxSize;
  std::vector<size_t> neighbours = std::vector<size_t>(0);
  std::vector<double> distances = std::vector<double>(0);
  std::vector<double> weights = std::vector<double>(0);

  KDStore();
  KDStore(const size_t, const size_t);
  ~KDStore();
  void Set(const size_t, const size_t);
  void PrepareWeights();
  void Reset();

  size_t GetSize();
  bool SizeFulfilled();
  double MinimumDistance();
  double MaximumDistance();

  double GetDistance(const size_t);
  void SetDistance(const size_t, const double);
  double GetWeight(const size_t);
  void SetWeight(const size_t, const double);
  void AddUnit(const size_t);
  void AddUnitAndReset(const size_t);

  void SortNeighboursByDistance(const size_t, const size_t);
  void SortNeighboursByWeight(const size_t, const size_t);
};

#endif
