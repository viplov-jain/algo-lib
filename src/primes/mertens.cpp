#include <stdio.h>
#include <unordered_map>
#include <vector>
#include <string>
#include <math.h>
#define PRE_LIM 8000000
#define int128_t __int128

std::vector<int32_t> mertens_s(PRE_LIM, 0);
std::unordered_map<int64_t, int32_t> mertens_map;

void init();
int32_t mertens_compute(int64_t n);

int main() {
  init();
  int64_t val = 235711131719L;
  for (int i = 1; i < 10; i++) {
    printf("%i ", mertens_compute(val - i));
  }
  printf("\n");
}

void init() {
  int pre_sqrt = sqrt(PRE_LIM);
  std::vector<bool> is_pr(PRE_LIM, true), mu_p(PRE_LIM, false), is_sqrfr(PRE_LIM, true);
  is_pr[0] = is_pr[1] = false;
  mu_p[1] = true;
  for (int i = 2; i < PRE_LIM; i++) {
    if (is_pr[i]) {
      for (int j = i * 2; j < PRE_LIM; j += i) {
        is_pr[j] = false;
        mu_p[j] = !mu_p[j / i];
      }
      if (i <= pre_sqrt) {
        int i2 = i * i;
        for (int j = i2; j < PRE_LIM; j += i2) {
          is_sqrfr[j] = false;
        }
      }
    }
  }
  for (int i = 1; i < PRE_LIM; i++) {
    mertens_s[i] = mertens_s[i - 1] + (is_sqrfr[i] ? (mu_p[i] ? 1 : -1) : 0);
  }
}

int32_t mertens_compute(int64_t n) {
  if (n < PRE_LIM) {
    return mertens_s[n];
  }
  auto val = mertens_map.find(n);
  if (val != mertens_map.end()) {
    return val->second;
  }

  int64_t ret = 1, v = sqrt(n) + 1;
  for (int i = 1; i <= v; i++) {
    ret -= mertens_compute(i) * ((n / i) - (n / (i + 1)));
  }

  for (int i = 2; v < n / i; i++) {
    ret -= mertens_compute(n / i);
  }

  mertens_map[n] = ret;
  return ret;
}
