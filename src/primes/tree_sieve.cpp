#include <stdio.h>
#include <algorithm>
#include <vector>
#define PRE_LIM 100

typedef unsigned long i64;
typedef int i32;

std::vector<i32> primes;
std::vector<i32> lst;

void sieve(
  std::vector<bool>& is_prime,
  const std::vector<i32>& vals,
  i32 l,
  i32 h,
  i32 m) {
  if (h == l + 1) {
    is_prime[l] = false;
    return;
  }
  i32 mid = (l + h) / 2;
  i32 v_mid = vals[mid] / m;

  if (vals[l] / m != v_mid) {
    sieve(is_prime, vals, l, mid, m);
  }
  if (vals[h - 1] / m != v_mid) {
    sieve(is_prime, vals, mid, h, m);
  }
}

void init() {
  primes.reserve(6542);
  std::vector<bool> is_prime(65536, true);
  for (int i = 2; i < 256; i++) {
    if (is_prime[i]) {
      for (int j = i * i; j < 65536; j += i) {
        is_prime[j] = false;
      }
    }
  }
  for (int i = 2; i < 65536; i++) {
    if (is_prime[i]) {
      primes.push_back(i);
    }
  }
  int pl = 0;
  while (primes[pl] < PRE_LIM) {
    pl++;
  }
  // printf("%lu\n", primes.size());
  // for (int i = 0; i < 20; i++) {
  //   printf("(%i %i)", i, primes[i]);
  // }
  // double val = 0;
  // for (int i = 9; i < primes.size(); i++) {
  //   val += 1.0 / primes[i];
  // }
  // printf("%lf\n", val);
  // return;
  // lst.reserve(33333333);
  lst.reserve(1 << 23);
  i32 a = 1234567891, mask = (1U << 31) - 1;
  for (int i = 0; i < 33333333; i++) {
    a += 1234567890;
    a &= mask;

    bool br = false;
    for (int i = 1; i < pl; i++) {
      if (a % primes[i] == 0) {
        br = true;
        break;
      }
    }
    if (!br) {
      lst.push_back(a);
    }
  }
  printf("%li\n", lst.size());
  std::sort(lst.begin(), lst.end());
  for (int i = 0; i < 100; i++) {
    printf("%i ", lst[i]);
  }
  fflush(stdout);
  is_prime = std::vector<bool>(1 << 23, true);
  for (const auto& p: primes) {
    if (p >= 10000) {
      sieve(is_prime, lst, 0, lst.size(), p);
      // printf("%i ", p);
      // fflush(stdout);
    }
  }

  return;
}

int main() {

}
