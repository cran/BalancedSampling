#include <algorithm>
#include <stddef.h>
#include <utility>
#include <vector>

#include "KDNodeClass.h"

KDNode::KDNode(KDNode* t_parent, const size_t t_size) {
  parent = t_parent;
  SetTerminal(t_size);
}

KDNode::KDNode(KDNode* t_parent, const bool t_bool) {
  parent = t_parent;
  SetTerminal(t_bool);
}

KDNode::~KDNode() {
  delete cleft;
  delete cright;
}

void KDNode::Copy(const KDNode* original) {
  if (terminal) {
    ReplaceUnits(original->units);
    return;
  }

  split = original->split;
  value = original->value;

  cleft = new KDNode(this, original->cleft->IsTerminal());
  cleft->Copy(original->cleft);

  cright = new KDNode(this, original->cright->IsTerminal());
  cright->Copy(original->cright);

  return;
}

void KDNode::Prune(const size_t bucketSize) {
  // We Can't prune a terminal node
  if (terminal)
    return;

// If a node is not terminal, but is missing a child, something went wrong
  if (cleft == nullptr || cright == nullptr)
    return;

  if (!cleft->terminal)
    cleft->Prune(bucketSize);

  if (!cright->terminal)
    cleft->Prune(bucketSize);

  if (!cleft->terminal || !cright->terminal)
    return;

  // Both nodes are terminal
  size_t nunits = cleft->GetSize() + cright->GetSize();
  if (nunits > bucketSize)
    return;

  // The nodes are too small, prune
  units.reserve(nunits);
  units.resize(0);

  size_t* tarr = cleft->units.data();
  for (size_t i = 0; i < cleft->GetSize(); i++)
    units.push_back(tarr[i]);

  tarr = cright->units.data();
  for (size_t i = 0; i < cright->GetSize(); i++)
    units.push_back(tarr[i]);

  terminal = true;

  delete cleft;
  delete cright;
  cleft = nullptr;
  cright = nullptr;

  return;
}

void KDNode::SetTerminal(const size_t t_size) {
  terminal = t_size > 0;
  return;
}

void KDNode::SetTerminal(const bool t_bool) {
  terminal = t_bool;
  return;
}

bool KDNode::IsTerminal() {
  return terminal;
}

KDNode* KDNode::GetSibling() {
  if (parent == nullptr)
    return nullptr;

  return this == parent->cleft ? parent->cright : parent->cleft;
}

void KDNode::AddUnit(const size_t id) {
  if (!terminal)
    return;

  units.push_back(id);
  return;
}

void KDNode::ReplaceUnits(const size_t* t_units, const size_t nunits) {
  units.reserve(nunits);
  units.resize(0);

  for (size_t i = 0; i < nunits; i++)
    units.push_back(t_units[i]);

  return;
}

void KDNode::ReplaceUnits(const std::vector<size_t>& t_units) {
  ReplaceUnits(t_units.data(), t_units.size());
  return;
}


void KDNode::RemoveUnit(const size_t id) {
  size_t nunits = units.size();

  // If there are no units here, nothing to remove
  if (nunits == 0)
    return;

  size_t* dt = units.data();

  // Find the pointer to id in units
  size_t* it = std::find(dt, dt + nunits, id);

  // If nothing found, return
  if (it == dt + nunits)
    return;

  // Something found, thus we have one less unit
  nunits -= 1;

  // If the iterator is pointing to the position at the last pos,
  // we don't need to switch
  if (it != dt + nunits) {
    *it = dt[nunits];
  }

  units.pop_back();
  return;
}

bool KDNode::UnitExists(const size_t id) {
  if (units.size() == 0)
    return false;

  std::vector<size_t>::iterator it = std::find(units.begin(), units.end(), id);

  return it != units.end();
}

size_t KDNode::GetSize() {
  if (!terminal)
    return 0;

  return units.size();
}
