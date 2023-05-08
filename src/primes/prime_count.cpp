#include <iostream>
#include <math.h>
#include <vector>
#include <unordered_map>
#include <map>
#include <cstring>
#include <unordered_set>
#define PR_LIM0 40000000
#define PR_LIM 2500000
std::vector<int32_t> pr_arr, pi_arr;
std::unordered_map<int64_t, int64_t> pi_vals;
std::vector<int64_t> phi_vals(510510 * 7, 0);
std::vector<int64_t> phi_vals2(169000, 0);
std::unordered_map<int64_t, int64_t> phi_vals3;
int64_t phi(int64_t x, int32_t a);
int64_t pi(int64_t x);

void init();

int64_t phi(int64_t x, int32_t a) {
  if (!a || !x) {
    return x;
  }
  if (x < pr_arr[a]) {
    return 1;
  }
  if (a < 8) {
    int64_t v1 = phi_vals[(a - 1) * 510510 + (x % 510510)];
    return 92160 * (x / 510510) + v1;
  }
  if (x < 1000 && a < 169) {
    return phi_vals2[a * 1000 + x];
  }
  std::unordered_map<int64_t, int64_t>::iterator find = phi_vals3.find(x * 169 + a);
  if (find != phi_vals3.end()) {
    return find->second;
  }

  int64_t val = phi(x, a - 1) - phi(x / pr_arr[a], a - 1);
  phi_vals3[x * 169 + a] = val;
  return val;
}

int64_t pi(int64_t x) {
  if (x < PR_LIM0) {
    return pi_arr[x];
  }
  std::unordered_map<int64_t, int64_t>::iterator find = pi_vals.find(x);
  if (find != pi_vals.end()) {
    return find->second;
  }
  int64_t a = sqrt(sqrt(x)), b = sqrt(x), c = cbrt(x);
  a = pi_arr[a];
  b = pi_arr[b];
  c = pi_arr[c];
  int64_t val = phi(x, a) + ((a + b - 2) * (b - a + 1)) / 2;
  for (int i = a + 1; i <= b; i++) {
    val -= pi(x / pr_arr[i]);
  }
  for (int i = a + 1; i <= c; i++) {
    int32_t bi = sqrt(x / pr_arr[i]);
    bi = pi_arr[bi];
    for (int j = i; j <= bi; j++) {
      val -= pi(x / ((int64_t)pr_arr[i] * (int64_t)pr_arr[j])) - j + 1;
    }
  }
  pi_vals[x] = val;
  return val;
}

void sieve(int n, bool is_prime[], std::vector<int> &l) {
  // std::memset(is_prime, true, sizeof(is_prime));
  for (int p = 2; p * p < n; p++) {
    if (is_prime[p]) {
      for (int i = p * 2; i < n; i += p) {
        is_prime[i] = false;
      }
    }
  }
  for (int p = 2; p < n; p++) {
    if (is_prime[p]) {
      l.push_back(p);
    }
  }
}

void segmented_sieve(int n, std::vector<int> &pi, std::vector<int> &pr) {
  pr.push_back(0);
  int segment_s = floor(sqrt(PR_LIM0)) + 1;
  std::vector<int32_t> primes;
  bool is_prime[segment_s + 1];
  std::memset(is_prime, true, sizeof(is_prime));
  sieve(segment_s, is_prime, primes);
  for (int i = 2; i < segment_s; i++) {
    if (is_prime[i]) {
      pr.push_back(i);
    }
    pi[i] = pr.size() - 1;
  }
  int low = segment_s;
  int high = 2 * segment_s;

  while (low < n) {
    if (high >= n) {
      high = n;
    }
    std::memset(is_prime, true, sizeof(is_prime));
    for (int i = 0; i < primes.size(); i++) {
      int j = floor(low / primes[i]) * primes[i];
      if (j < low) {
        j += primes[i];
      }
      for (; j < high; j += primes[i]) {
        is_prime[j - low] = false;
      }
    }
    for (int i = low; i < high; i++) {
      if (is_prime[i - low] == true) {
        pr.push_back(i);
      }
      pi[i] = pr.size() - 1;
    }
    low = low + segment_s;
    high = high + segment_s;
  }
}

void init() {
  pi_arr.resize(PR_LIM0);
  pr_arr.reserve(PR_LIM);
  segmented_sieve(PR_LIM0, pi_arr, pr_arr);

  for (int n = 0; n < 510510; n++) {
    phi_vals[n] = n - (n / 2);
  }
  for (int a = 2; a < 8; a++) {
    for (int n = 0; n < 510510; n++) {
      phi_vals[(a - 1) * 510510 + n] = phi_vals[(a - 2) * 510510 + n] - phi_vals[(a - 2) * 510510 + n / pr_arr[a]];
    }
  }
  for (int n = 0; n < 1000; n++) {
    phi_vals2[n] = n;
  }
  for (int a = 1; a < 169; a++) {
    for (int n = 0; n < 1000; n++) {
      phi_vals2[a * 1000 + n] = phi_vals2[(a - 1) * 1000 + n] - phi_vals2[(a - 1) * 1000 + n / pr_arr[a]];
    }
  }
}

int main() {
  init();
  for (int i = 1; i < sqrt(235711131719L); i++) {
    pi(235711131719L / i);
  }
  return 0;
}
