#ifndef CUBECLASS_HEADER
#define CUBECLASS_HEADER

#include <stddef.h>
#include <vector>

#include "IndexListClass.h"
#include "KDStoreClass.h"
#include "KDTreeClass.h"

enum class CubeMethod {
  err = 0,
  cube = 1,
  lcube = 2
};

CubeMethod IntToCubeMethod(const int i);

class Cube {
protected:
  bool set_direct = false;
  bool set_draw = false;

  void (Cube::*_Draw)() = nullptr;

public:
  CubeMethod cubeMethod;
  size_t N;
  size_t pbalance;
  double eps = 1e-12;

  IndexList* idx = nullptr;
  KDTree* tree = nullptr;
  KDStore* store = nullptr;

  std::vector<double> probabilities = std::vector<double>(0);
  std::vector<double> amat = std::vector<double>(0);

protected:
  std::vector<size_t> candidates = std::vector<size_t>(0);
  std::vector<double> bmat = std::vector<double>(0);
  std::vector<double> uvec = std::vector<double>(0);

public:
  std::vector<size_t> sample = std::vector<size_t>(0);

  // DIRECT CUBE
  Cube(
    const double* t_probabilities,
    double* xxbalance,
    const size_t t_N,
    const size_t t_pbalance,
    const double t_eps
  );
  // DIRECT LCUBE
  Cube(
    const double* t_probabilities,
    double* xxbalance,
    const size_t t_N,
    const size_t t_pbalance,
    const double t_eps,
    double* xxspread,
    const size_t t_pspread,
    const size_t t_treeBucketSize,
    const int t_treeMethod
  );
  // INDIRECT CUBE
  Cube(const CubeMethod, const size_t t_N, const size_t t_pbalance, const double t_eps);

  void Init(
    const double* t_probabilities,
    double* xxbalance,
    const size_t t_N,
    const size_t t_pbalance,
    const double t_eps
  );

  void InitIndirect(const size_t t_N, const size_t t_pbalance, const double t_eps);

  ~Cube();

public:
  void AddUnitToSample(const size_t id);
  void EraseUnit(const size_t id);

protected:
  size_t MaxSize();

  void Draw_cube();
  void Draw_lcube();

  void Draw();
  void RunUpdate();

public:
  void RunFlight();
  void RunLanding();
  void Run();
};

#endif
