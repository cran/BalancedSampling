#ifndef CUBESTRATIFIEDCLASS_HEADER
#define CUBESTRATIFIEDCLASS_HEADER

#include <stddef.h>
#include <unordered_map>
#include <vector>

#include "CubeClass.h"
#include "IndexListClass.h"
#include "KDTreeClass.h"

class CubeStratified {
public:
  CubeMethod cubeMethod;
  size_t N;
  size_t pbalance;
  size_t pspread;
  double eps = 1e-12;

  IndexList* idx = nullptr;

  std::vector<double> probabilities = std::vector<double>(0);

protected:
  double* r_prob = nullptr; // NO DELETE
  double* r_xbalance = nullptr; // NO DELETE
  double* r_xspread = nullptr; // NO DELETE
  int* r_strata = nullptr; // NO DELETE

  std::vector<double> zspread = std::vector<double>(0);

  size_t treeBucketSize = 40;
  KDTreeSplitMethod treeMethod = KDTreeSplitMethod::midpointSlide;

  std::unordered_map<int, size_t> stratumMap;
  std::vector<int> stratumArr = std::vector<int>(0);
  std::vector<size_t> index = std::vector<size_t>(0);

public:
  std::vector<int> sample;

public:
  CubeStratified(
    int* t_strata,
    double* t_probabilities,
    double* t_xbalance,
    const size_t t_N,
    const size_t t_pbalance,
    const double t_eps
  );

  CubeStratified(
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
  );

  void Init(
    int* t_strata,
    double* t_probabilities,
    double* t_xbalance,
    const size_t t_N,
    const size_t t_pbalance,
    const double t_eps
  );

  ~CubeStratified();

public:
  void AddUnitToSample(const size_t);

private:
  void RunFlightPerStratum();
  void RunFlightOnFull();
  void RunLandingPerStratum();

public:
  void Run();
};

#endif
