#pragma GCC optimize ("O3")
#include <stdio.h>
#include <vector>
#include <math.h>
#include <algorithm>
using uint8_v = u_int8_t __attribute__((vector_size(8)));
using uint32_v = u_int32_t __attribute__((vector_size(32)));

#define CASE(M, I) \
  case M:\
    sv[I] = (s * (s + a) - low) / 30; \
    break;

uint8_v MASK = {254, 253, 251, 247, 239, 223, 191, 127};
uint8_v BIT = {1, 2, 4, 8, 16, 32, 64, 128};
uint32_v MOD_V = {1, 7, 11, 13, 17, 19, 23, 29};
uint32_v ZERO = {0, 0, 0, 0, 0, 0, 0, 0};

inline void print(uint32_v val) {
    printf("{");
    for (int i = 0; i < 7; i++) printf("%u ", val[i]);
    printf("%u}\n", val[7]);
}

inline __attribute__((always_inline)) bool less(uint32_v v, u_int32_t u) {
    for (int i = 0; i < 8; i++) if (v[i] >= u) return false;
    return true;
}

int main() {
  long sum = 129; // 2 + 3 + 5 + 7 + 11 + 13 + 17 + 19 + 23 + 29
  const u_int32_t lim = 2147483640; // 30 | (lim + 1)
  const u_int32_t segment_size = 262080 / 8; // 30 | segment_size
  u_int32_t lim_sqrt = sqrt(lim);

  u_int32_t i = 3;
  u_int32_t n = 1, n30;
  u_int32_t s = 3;

  std::vector<u_int8_t> sieve(segment_size);
  std::vector<bool> is_prime(lim_sqrt + 1, true);
  std::vector<bool> wheel_pos(lim_sqrt + 1);
  std::vector<u_int32_t> small_primes;
  // 1, 7, 11, 13, 17, 19, 23, 29
  std::vector<uint32_v> mod_grp;

  u_int32_t count = 10; // 2, 3, 5, 7, 11, 13, 17, 19, 23, 29

  u_int32_t low30;
  for (u_int32_t low = 0; low <= lim; low += segment_size * 30) {
    std::fill(sieve.begin(), sieve.end(), 255);
    u_int32_t high = low + segment_size * 30 - 1;
    high = std::min(high, lim);
    low30 = low / 30;

    for (; i * i <= high; i += 2) {
      if (is_prime[i]) {
        for (u_int32_t j = i * i; j <= lim_sqrt; j += i) {
          is_prime[j] = false;
        }
      }
    }

    for (; s * s <= high; s += 2) {
      if (is_prime[s] && s > 5) {
        uint32_v sv;
        for (u_int32_t a = 0; a < 30; a += 2) {
          switch (s * (s + a) % 30) {
            CASE(01, 0)
            CASE(07, 1)
            CASE(11, 2)
            CASE(13, 3)
            CASE(17, 4)
            CASE(19, 5)
            CASE(23, 6)
            CASE(29, 7)
          }
        }
        small_primes.push_back(s);
        mod_grp.push_back(sv);
      }
    }
    for (u_int32_t i = 0; i < small_primes.size(); i++) {
      const u_int32_t k = small_primes[i];
      uint32_v j = mod_grp[i];
    //   print(j);
      for (; less(j, segment_size); j += k) {
        #pragma GCC ivdep
        for (int jj = 0; jj < 8; jj ++) sieve[j[jj]] &= MASK[jj];
      }
      for (int jj = 0; jj < 8; jj ++) {
        for (; j[jj] < segment_size; j[jj] += k) sieve[j[jj]] &= MASK[jj];
      }
      mod_grp[i] = j - segment_size;
    }
    for (; (n30 = 30 * n) <= high; n++) {
        uint32_v is_pr;
        uint8_v is_pr_s = sieve[n - low30] & BIT;
        #pragma GCC ivdep
        for (int jj = 0; jj < 8; jj ++) is_pr[jj] = is_pr_s[jj];

        uint32_v pr_vec = is_pr ? (n30 + MOD_V) : ZERO;
        // #pragma GCC ivdep
        // for (int jj = 0; jj < 8; jj ++) {
        //     sum += pr_vec[jj];
        //     if (pr_vec[jj]) count++;
        // }
    }
  }
  printf("%i %li\n", count, sum);
  return 0;
}
