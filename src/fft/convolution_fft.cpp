#include <stdio.h>
#include <vector>
#include <complex>
#include <iostream>
#include <bitset>
#define int64_t long long
#define uint64_t unsigned long long
#define N_LIM 20
#define Complex std::complex<double>
#define MAX(a, b) (a) > (b) ? (a) : (b)

double PI = std::acos(-1);
Complex W[(1 << N_LIM)];
int REV[(1 << N_LIM)];

void init_fft(int L) {
  for (int r = 0; !(r >> L); r++) {
    double th = 2 * PI * r / pow(2.0, L);
    W[r] = Complex(std::cos(th), std::sin(th));
    REV[r] = (REV[r >> 1] >> 1) | ((r & 1) << (L - 1));
  }
}

void fft(Complex* A, int L) {
  for (int i = 0; !(i >> L); i++) {
    if (i > REV[i]) {
      std::swap(A[i], A[REV[i]]);
    }
  }
  for (int s = 0; s < L; s++) {
    int m = 1 << s, m2 = m << 1, rm2 = 1 << (L - s - 1);
    for (int k = 0; !(k >> L); k += m2) {
      for (int j = 0; !(j >> s); j++) {
        Complex &As = A[k + j], &At = A[k + j + m];
        Complex t = W[rm2 * j] * At;
        At = As - t;
        As = As + t;
      }
    }
  }
}

Complex AA[1 << N_LIM], BB[1 << N_LIM];
std::vector<int64_t> convolution(const std::vector<int>& A, const std::vector<int>& B) {
  std::vector<int64_t> res(A.size() + B.size() - 1);
  if (A.size() < 50 || B.size() < 50) {
    for (int i = 0; i < A.size(); i++) {
      for (int j = 0; j < B.size(); j++) {
        res[i + j] += (int64_t)A[i] * (int64_t)B[j];
      }
    }
    return res;
  }
  int L = 1;
  while ((A.size() + B.size() - 1) >> L) L++;
  init_fft(L);
  for (int i = 0; !(i >> L); i++) {
    AA[i] = i < A.size() ? A[i] : 0;
  }
  for (int i = 0; !(i >> L); i++) {
    BB[i] = i < B.size() ? B[i] : 0;
  }
  fft(AA, L);
  fft(BB, L);
  for (int i = 0; !(i >> L); i++) {
    AA[i] *= BB[i];
  }
  for (int i = 0; !(i >> L); i++) {
    W[i] = std::conj(W[i]);
  }
  fft(AA, L);
  for (int i = 0; i < res.size(); i++) {
    res[i] = (int64_t)((std::real(AA[i]) / (double)(1 << L)) + 0.5);
  }
  return res;
}

std::vector<int64_t> self_convolution(const std::vector<int>& A) {
  std::vector<int64_t> res(2 * A.size() - 1);
  if (A.size() < 50) {
    for (int i = 0; i < A.size(); i++) {
      for (int j = 0; j < A.size(); j++) {
        res[i + j] += (int64_t)A[i] * (int64_t)A[j];
      }
    }
    return res;
  }
  int L = 1;
  while ((2 * A.size() - 1) >> L) L++;
  init_fft(L);
  for (int i = 0; !(i >> L); i++) {
    AA[i] = i < A.size() ? A[i] : 0;
  }
  fft(AA, L);
  for (int i = 0; !(i >> L); i++) {
    AA[i] *= AA[i];
  }
  for (int i = 0; !(i >> L); i++) {
    W[i] = std::conj(W[i]);
  }
  fft(AA, L);
  for (int i = 0; i < res.size(); i++) {
    res[i] = (int64_t)((std::real(AA[i]) / (double)(1 << L)) + 0.5);
  }
  return res;
}

#define BASE_MASK 65535
#define BASE_L 16

std::vector<int> mul(const std::vector<int>& a, const std::vector<int>& b) {
  if (a.size() < b.size()) return mul(b, a);
  std::vector<int64_t> R = &a == &b ? self_convolution(a) : convolution(a, b);
  std::vector<int> res;
  res.reserve(a.size() + b.size() - 1);
  int64_t carry = 0, val;
  for (int i = 0; i < R.size(); i++) {
    val = R[i] + carry;
    res.push_back(val & BASE_MASK);
    carry = val >> BASE_L;
  }
  while (carry) {
    res.push_back(carry & BASE_MASK);
    carry >>= BASE_L;
  }
  return res;
}

std::vector<int> sq(const std::vector<int>& a) {
  std::vector<int64_t> R = self_convolution(a);
  std::vector<int> res;
  res.reserve(2 * a.size() - 1);
  int64_t carry = 0, val;
  for (int i = 0; i < R.size(); i++) {
    val = R[i] + carry;
    res.push_back(val & BASE_MASK);
    carry = val >> BASE_L;
  }
  while (carry) {
    res.push_back(carry & BASE_MASK);
    carry >>= BASE_L;
  }
  return res;
}

void mul(std::vector<int>& inp, int v) {
  uint64_t carry = 0, val;
  for (int i = 0; i < inp.size(); i++) {
    val = (uint64_t)inp[i] * (uint64_t)v + carry;
    inp[i] = val & BASE_MASK;
    carry = val >> BASE_L;
  }
  while (carry) {
    inp.push_back(carry & BASE_MASK);
    carry >>= BASE_L;
  }
}
