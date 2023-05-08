#include <stdio.h>
#include <vector>
#include <random>
#include <unordered_set>

typedef __int128 int128_t;
#define MOD 1000000007
#define LIM 1048576

std::vector<bool> is_pr(LIM, true);
std::vector<int> fctr(LIM);

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

void factorize(int64_t N, std::unordered_set<int64_t>& factors) {
  if (N == 1) {
    return;
  }
  if (!(N & 1)) {
    factors.insert(2);
    factorize(N >> 1, factors);
    return;
  }
  if (is_prime_mr(N)) {
    factors.insert(N);
    return;
  }
  if (N < LIM) {
    factors.insert(fctr[N]);
    factorize(N / fctr[N], factors);
    return;
  }
  int128_t x, y, c, q, ys, k;
  int64_t m, g, r;

  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<int64_t> distr(2, N - 1);

  y = distr(gen); c = distr(gen); m = distr(gen);
  r = q = 1;
  g = gcd(N, sqrt(N));

  while (g == 1) {
    x = y;
    for (int i = 0; i < r; i++) {
      y = (y * y + c) % N;
    }
    k = 0;
    while (k < r && g == 1) {
      ys = y;
      for (int i = 0; i < m && i < r - k; i++) {
        y = (y * y + c) % N;
        q = (q * (x > y ? x - y : y - x)) % N;
      }
      g = gcd(q, N);
      k += m;
    }
    r <<= 1;
  }

  if (g == N) {
    do {
      ys = (ys * ys + c) % N;
      g = gcd((x > ys ? x - ys : ys - x), N);
    } while (g < 2);
  }
  factorize(g, factors);
  factorize(N / g, factors);
}

void init() {
  for (int i = 0 ; i < LIM; i++) {
    fctr[i] = i;
  }
  is_pr[0] = is_pr[1] = false;
  for (int i = 2; i < 1024; i++) {
    if (!is_pr[i]) {
      continue;
    }
    for (int j = i * i; j < LIM; j += i) {
      is_pr[j] = false;
      fctr[j] = i;
    }
  }
}

int main() {
  init();
  std::unordered_set<int64_t> factors;
  factorize(4732892, factors);
  for (auto fc: factors) {
    printf("%lli ", fc);
  }
}
