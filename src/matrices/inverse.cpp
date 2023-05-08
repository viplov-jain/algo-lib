#include <stdio.h>
#include <vector>

#define AA(i, j) A[n * (i) + j]
#define BB(i, j) B[n * (i) + j]

template <typename T>
void swap(int i, int j, T *A, T *B, const int n) {
  if (i == j) return;
  for (int k = 0; k < n; k++) {
    std::swap(AA(i, k), AA(j, k));
    std::swap(BB(i, k), BB(j, k));
  }
}

template <typename T>
void scale(const int i, const T a, T *A, T *B, const int n) {
  if (a == 1) return;
  for (int k = 0; k < n; k++) {
    AA(i, k) *=  a;
    BB(i, k) *=  a;
  }
}

template <typename T>
void add(const int i, const int j, const T a, T *A, T *B, const int n) {
  if (a == 0) return;
  for (int k = 0; k < n; k++) {
    AA(j, k) = (AA(j, k) + a * AA(i, k));
    BB(j, k) = (BB(j, k) + a * BB(i, k));
  }
}

template <typename T>
bool invert(T *A, T *B, const int n) {
  for (int i = 0; i < n; i++) {
    for (int j = i; j < n; j++) {
      if (abs(AA(j, i)) > 1e-6) {
        swap(i, j, A, B, n);
        break;
      }
    }
    if (abs(AA(i, i)) < 1e-6) {
      return false;
    }
    scale(i, 1.0 / AA(i, i), A, B, n);
    for (int j = i + 1; j < n; j++) {
      add(i, j, -AA(j, i), A, B, n);
    }
  }
  for (int i = 0; i < n; i++) {
    for (int j = i + 1; j < n; j++) {
      add(j, i, -AA(i, j), A, B, n);
    }
  }
  return true;
}

int main() {
  return 0;
}
