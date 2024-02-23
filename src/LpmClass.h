#ifndef LPMCLASS_HEADER
#define LPMCLASS_HEADER

#include <stddef.h>
#include <vector>

#include "IndexListClass.h"
#include "KDStoreClass.h"
#include "KDTreeClass.h"

enum class LpmMethod {
  err = 0,
  lpm1 = 1,
  lpm2 = 2,
  lpm1search = 3,
  rpm = 4,
  spm = 5
};

LpmMethod IntToLpmMethod(const int i);

class Lpm {
protected:
  bool set_direct = false;
  bool set_draw = false;
  bool set_run = false;

  void (Lpm::*_Draw)() = nullptr;
  void (Lpm::*_Run)() = nullptr;

public:
  LpmMethod lpMethod;
  size_t N;
  double eps = 1e-12;

  IndexList* idx = nullptr;
  KDTree* tree = nullptr;
  KDStore* store = nullptr;

  std::vector<double> probabilities = std::vector<double>(0);
  std::vector<size_t> iprobabilities = std::vector<size_t>(0);

protected:
  size_t pair [2] = {0, 1};
  std::vector<size_t> candidates = std::vector<size_t>(0);
  std::vector<size_t> history = std::vector<size_t>(0);

public:
  std::vector<size_t> sample = std::vector<size_t>(0);

  // DIRECT
  // -- DOUBLE
  Lpm(
    const LpmMethod t_lpMethod,
    const double* t_probabilities,
    double* xx,
    const size_t t_N,
    const size_t t_p,
    const double t_eps,
    const size_t t_treeBucketSize,
    const int t_treeMethod
  );
  // -- INT
  Lpm(
    const LpmMethod t_lpMethod,
    size_t t_probn,
    double* xx,
    const size_t t_N,
    const size_t t_p,
    const size_t t_treeBucketSize,
    const int t_treeMethod
    );

  void Init(
    double* xx,
    const size_t t_N,
    const size_t t_p,
    const size_t t_bucketSize,
    const int t_method,
    const LpmMethod t_lpMethod
  );

  ~Lpm();

public:
  void AddUnitToSample(const size_t id);
  void EraseUnit(const size_t id);

protected:
  void Draw_lpm1();
  void Draw_lpm2();
  void Draw_lpm1search();
  void Draw_rpm();
  void Draw_spm();

  void Run_double();
  void Run_int();

  void Draw();

public:
  void Run();
};

#endif
