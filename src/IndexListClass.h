#ifndef INDEXLISTCLASS_HEADER
#define INDEXLISTCLASS_HEADER

#include <stddef.h>

class IndexList {
private:
  size_t* list = nullptr;
  size_t* reverse = nullptr;
  size_t len = 0;
  size_t capacity = 0;
public:
  IndexList();
  IndexList(const size_t);
  ~IndexList();
  IndexList* Copy();
  IndexList* CopyLen();
  size_t* CopyList();

  size_t Length();
  void Fill();
  void Reset();
  void Resize(const size_t);
  void Shuffle();

  void Set(const size_t);
  void Add(const size_t);
  size_t Get(const size_t);
  size_t GetK(const size_t);
  size_t GetLast();
  bool Exists(const size_t);
  size_t Draw();
  void Erase(const size_t);
};

#endif
