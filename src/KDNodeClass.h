#ifndef KDNODECLASS_HEADER
#define KDNODECLASS_HEADER

#include <stddef.h>
#include <vector>

class KDNode {
  // REGULAR NODE
public:
  KDNode* parent = nullptr;
  KDNode* cleft = nullptr; // <=
  KDNode* cright = nullptr; // >

  size_t split;
  double value = 0.0;

private:
  bool terminal = false;
public:
  std::vector<size_t> units = std::vector<size_t>(0);

public:
  KDNode(KDNode*, const size_t);
  KDNode(KDNode*, const bool);
  ~KDNode();

  void Copy(const KDNode*);
  void Prune(const size_t);

  void SetTerminal(const size_t);
  void SetTerminal(const bool);
  bool IsTerminal();

  KDNode* GetSibling();

  void AddUnit(const size_t);
  void ReplaceUnits(const size_t*, const size_t);
  void ReplaceUnits(const std::vector<size_t>&);
  void RemoveUnit(const size_t);
  bool UnitExists(const size_t);
  size_t GetSize();
};

#endif
