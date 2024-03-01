#ifndef CUBESTRATIFIEDCLASS_HEADER
#define CUBESTRATIFIEDCLASS_HEADER

#include <stddef.h>
#include <unordered_map>
#include <vector>

#include "CubeClass.h"
#include "IndexListClass.h"
#include "KDTreeClass.h"

using StratumMap = std::unordered_map<int, size_t>;

class CubeStratified {
public:
  CubeMethod cube_method_;
  size_t N_;
  size_t p_balance_;
  size_t p_spread_;
  double eps_ = 1e-12;

  IndexList* idx_ = nullptr;
  Cube* cube_;

protected:
  double* rptr_probabilities_ = nullptr; // NO DELETE
  double* rptr_xbalance_ = nullptr; // NO DELETE
  double* rptr_xspread_ = nullptr; // NO DELETE
  int* rptr_strata_ = nullptr; // NO DELETE

  size_t tree_bucket_size_ = 40;
  KDTreeSplitMethod tree_method_ = KDTreeSplitMethod::midpointSlide;

  StratumMap stratum_map_;
  std::vector<int> stratum_arr_;

public:
  std::vector<size_t> sample_;

public:
  // CUBE
  CubeStratified(
    int* strata,
    double* probabilities,
    double* xbalance,
    const size_t N,
    const size_t pbalance,
    const double eps
  );

  // LCUBE
  CubeStratified(
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
  );

  ~CubeStratified();

  void Init(
    int* strata,
    double* probabilities,
    double* xbalance,
    const size_t N,
    const size_t pbalance,
    const double eps
  );

private:
  void PrepareAmat(const size_t);
  void PrepareAmatAux(const size_t, const size_t);
  void RunFlightPerStratum();
  void RunFlightOnFull();
  void RunLandingPerStratum();

public:
  void Run();
};

#endif

