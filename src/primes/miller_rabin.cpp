#include <stdio.h>
#include <vector>
#include <random>
#include <unordered_set>

typedef __int128 int128_t;
#define MOD 1000000007
#define LIM 1048576

std::vector<bool> is_pr(LIM, true);

int64_t gcd(int64_t a, int64_t b) {
  while (b) std::swap(a = a % b, b);
  return a;
}

int64_t exp_(int64_t a, int64_t e, int64_t m) {
	if (!e) return 1;
	int64_t t = exp_(a, e >> 1, m);
	t = ((int128_t)t * (int128_t)t) % m;
	if (e & 1) return ((int128_t)t * (int128_t)a) % m;
	return t;
}

bool mr_pass(int64_t a, int s, int64_t d, int64_t n) {
	int64_t apow = exp_(a, d, n);
	if (apow == 1) return true;
	for (int i = 1 ; i < s; i++) {
		if (apow == n - 1) return true;
		apow = ((int128_t)apow * (int128_t)apow) % n;
  }
	return apow == n - 1;
}

bool is_prime_mr(int64_t n) {
  if (n < LIM) {
    return is_pr[n];
  }
	int64_t d = n - 1;
	int s = 0;
	while (!(d & 1)) {d >>= 1; s++;}
  std::vector<int> test_list;
  if (n < 2047) {
    test_list = {2};
  } else if (n < 9080191) {
    test_list = {31, 73};
  } else if (n < 4759123141L) {
    test_list = {2, 7, 61};
  } else if (n <  1122004669633L) {
    test_list = {2, 13, 23, 1662803};
  } else if (n < 2152302898747L) {
    test_list = {2, 3, 5, 7, 11};
  } else if (n < 3474749660383L) {
    test_list = {2, 3, 5, 7, 11, 13};
  } else if (n < 341550071728321L) {
    test_list = {2, 3, 5, 7, 11, 13, 17};
  } else if (n < 3825123056546413051L) {
    test_list = {2, 3, 5, 7, 11, 13, 17, 19, 23};
  } else {
    test_list = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};
  }
	for (int64_t tst: test_list) {
  	if (!mr_pass(tst, s, d, n)) {
  		return false;
  	}
	}
	return true;
}

void init() {
  is_pr[0] = is_pr[1] = false;
  for (int i = 2; i < 1024; i++) {
    if (!is_pr[i]) {
      continue;
    }
    for (int j = i * i; j < LIM; j += i) {
      is_pr[j] = false;
    }
  }
}

int main() {
  init();
  std::vector<int64_t> tsts{2, 3, 5, 4321, 4321343, 4324321};
  for (auto tst: tsts) {
    if (is_prime_mr(tst)) {
      printf("%lli ", tst);
    }
  }
}
