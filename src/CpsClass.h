#ifndef CPSCLASS_HEADER
#define CPSCLASS_HEADER

#include <stddef.h>
#include <vector>

#include "IndexListClass.h"
#include "KDStoreClass.h"
#include "KDTreeClass.h"

enum class CpsMethod {
  err = 0,
  lcps = 1,
  scps = 2,
  scpscoord = 3
};

CpsMethod IntToCpsMethod(const int i);

class Cps {
protected:
  bool set_direct = false;
  bool set_draw = false;
  bool set_random = false;

  size_t (Cps::*_Draw)() = nullptr;
  double (Cps::*_Random)(const size_t) = nullptr;

public:
  CpsMethod cpsMethod;
  size_t N;
  double eps = 1e-12;

  IndexList* idx = nullptr;
  KDTree* tree = nullptr;
  KDStore* store = nullptr;

  std::vector<double> probabilities = std::vector<double>(0);
  double* randomValues = nullptr;

protected:
  std::vector<size_t> candidates = std::vector<size_t>(0);
  size_t candidateIdx = 0;

public:
  std::vector<size_t> sample = std::vector<size_t>(0);

  // DIRECT
  Cps(
    const CpsMethod t_cpsMethod,
    const double* t_probabilities,
    double* xx,
    const size_t t_N,
    const size_t t_p,
    const double t_eps,
    const size_t t_treeBucketSize,
    const int t_treeMethod
  );

  ~Cps();

  void SetRandomStd();
  void SetRandomArr(double*);

public:
  void AddUnitToSample(const size_t id);
  void EraseUnit(const size_t id);
  void DecideUnit(const size_t id);

protected:
  size_t Draw_lcps();
  size_t Draw_scps();
  size_t Draw_scpscoord();

  double Random_std(const size_t);
  double Random_arr(const size_t);

  size_t Draw();
  double Random(const size_t);

public:
  void Run();
};

#endif
