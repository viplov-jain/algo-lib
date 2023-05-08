#include <stdio.h>
#include <vector>
#include <math.h>
#include <algorithm>

#define SEGMENTED_LOOP(X, Y) \
  j = X[i];\
  for (; j < segment_size; j += k) {\
    sieve[j] &= Y;\
  }\
  X[i] = j - segment_size;

#define CASE(X, M) \
  case M:\
    X.push_back((s * (s + a) - low) / 30); \
    break;

#define IF_BIT(B, M) \
  if (sieve[n - low30] & B) { \
    count++;\
    sum += n30 + M;\
  }

u_int8_t MASK0 = 254, MASK1 = 253, MASK2 = 251, MASK3 = 247, MASK4 = 239, MASK5 = 223, MASK6 = 191, MASK7 = 127;
u_int8_t BIT0 = 1, BIT1 = 2, BIT2 = 4, BIT3 = 8, BIT4 = 16, BIT5 = 32, BIT6 = 64, BIT7 = 128;

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
  std::vector<u_int32_t> mod_grp01, mod_grp07, mod_grp11, mod_grp13, mod_grp17, mod_grp19, mod_grp23, mod_grp29;

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
        small_primes.push_back(s);
        for (u_int32_t a = 0; a < 30; a += 2) {
          switch (s * (s + a) % 30) {
            CASE(mod_grp01, 1)
            CASE(mod_grp07, 7)
            CASE(mod_grp11, 11)
            CASE(mod_grp13, 13)
            CASE(mod_grp17, 17)
            CASE(mod_grp19, 19)
            CASE(mod_grp23, 23)
            CASE(mod_grp29, 29)
          }
        }
      }
    }

    for (u_int32_t i = 0; i < small_primes.size(); i++) {
      u_int32_t j;
      const u_int32_t k = small_primes[i];
      SEGMENTED_LOOP(mod_grp01, MASK0)
      SEGMENTED_LOOP(mod_grp07, MASK1)
      SEGMENTED_LOOP(mod_grp11, MASK2)
      SEGMENTED_LOOP(mod_grp13, MASK3)
      SEGMENTED_LOOP(mod_grp17, MASK4)
      SEGMENTED_LOOP(mod_grp19, MASK5)
      SEGMENTED_LOOP(mod_grp23, MASK6)
      SEGMENTED_LOOP(mod_grp29, MASK7)
    }
    for (; (n30 = 30 * n) <= high; n++) {
      IF_BIT(BIT0, 1)
      IF_BIT(BIT1, 7)
      IF_BIT(BIT2, 11)
      IF_BIT(BIT3, 13)
      IF_BIT(BIT4, 17)
      IF_BIT(BIT5, 19)
      IF_BIT(BIT6, 23)
      IF_BIT(BIT7, 29)
    }
  }
  printf("%i %li\n", count, sum);
  return 0;
}
