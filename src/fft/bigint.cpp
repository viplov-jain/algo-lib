#include <stdio.h>
#include <assert.h>
#include <vector>
#include <map>
#include <algorithm>
#include <cmath>
typedef __int128 i128;
typedef unsigned int i32;
typedef unsigned long i64;

#define POLY_MAX_SIZE 10000000 // 10^7
std::vector<i32> bit_rev_table(POLY_MAX_SIZE), factors(1 << 22);


struct bigint {
  std::vector<uint64_t> parts;
  bool negative = false;

  bigint(std::vector<uint64_t> vals, bool neg = false) {
    parts = vals;
    negative = neg;
  }
};


void add_h(std::vector<uint32_t> x, const std::vector<uint32_t>& y) {
  int xlen = x.size(), ylen = y.size(), d = 0;
  if (ylen > xlen) {
    x.reserve(ylen);
  }
  while (xlen < ylen) {
    x.push_back(0);
    xlen++;
  }
  for (int i = 0; i < xlen; i++) {
    x[i] += y[i] + d;
    if ((x[i] == y[i] && d == 1) || x[i] < y[i]) {
      d = 1;
    }
  }
  if (d == 1) {
    x.push_back(1);
  }
}

bool diff_h(std::vector<uint32_t> x, const std::vector<uint32_t>& y) {
  int xlen = x.size(), ylen = y.size(), d = 0;
  if (ylen > xlen) {
    x.reserve(ylen);
  }
  while (xlen < ylen) {
    x.push_back(0);
    xlen++;
  }
  for (int i = 0; i < xlen; i++) {
    x[i] -= y[i] + d;
    if ((x[i] == y[i] && d == 1) || x[i] > y[i]) {
      d = 1;
    }
  }
  if (d == 1) {
    for (int i = 0; i < xlen; i++) {
      x[i] = -x[i];
    }
    return true;
  }
  return false;
}



void init() {
  // Assuming 8 | BIT_REV_TABLE_K
  std::vector<uint16_t> lookup(1 << 8);
  for (uint16_t i = 0; !(i >> 8); i++) {
    uint16_t rev = 0;
    for (uint16_t j = 0; j < 8; j++) {
      if (i & (1 << j)) {
        rev |= 1 << (7 - j);
      }
    }
    lookup[i] = rev;
  }

  for (i32 i = 0; i < POLY_MAX_SIZE; i++) {
    i32 rev = 0;
    rev |= ((i32) lookup[i & 0xff]) << 24;
    rev |= ((i32) lookup[(i >> 8) & 0xff]) << 16;
    rev |= ((i32) lookup[(i >> 16) & 0xff]) << 8;
    rev |= ((i32) lookup[(i >> 24)]);
    bit_rev_table[i] = rev;
  }

  for (int i = 0; !(i >> 22); i++) {
    factors[i] = i;
  }
  for (int i = 2; !(i >> 11); i++) {
    if (factors[i] == i) {
      for (int j = i * i; !(j >> 22); j += i) {
        factors[j] = i;
      }
    }
  }
}

i32 mod_pow(i32 a, i32 r, i32 p) {
  if (r == 0) {
    return 1;
  }
  i32 val = mod_pow(a, r >> 1, p);
  val = ((i64) val * val) % p;
  if (r & 1) {
    return ((i64) val * a) % p;
  } else {
    return val;
  }
}

bool mr_pass(i32 a, i32 s, i32 d, i32 n) {
  i32 apow = mod_pow(a, d, n);
  if (apow == 1)
    return true;
  for (int i = 1; i < s; i++) {
    if (apow == n - 1)
      return true;
    apow = ((i64) apow * apow) % n;
  }
  return apow == n - 1;
}

bool is_prime(i32 n) {
  if (n == 2 || n == 7 || n == 61) {
    return true;
  }
  i32 d = n - 1;
  int s = 0;
  while (!(d & 1)) {
    d >>= 1;
    s++;
  }
  return mr_pass(2, s, d, n) && mr_pass(7, s, d, n) && mr_pass(61, s, d, n);
}

i32 primitive_root_check(i32 a, i32 p) {
  i32 t = p - 1;
  while (!(t & 1)) {
    t >>= 1;
  }
  if (mod_pow(a, (p - 1) / 2, p) == 1) {
    return false;
  }
  while (t != 1) {
    i32 pi = factors[t];
    while (!(t % pi)) {
      t /= pi;
    }
    if (mod_pow(a, (p - 1) / pi, p) == 1) {
      return false;
    }
  }
  return true;
}

i32 primitive_root(i32 p) {
  for (int a = 2;; a++) {
    if (primitive_root_check(a, p)) {
      return a;
    }
  }
}

i32 gcd(i32 a, i32 b) {
  if (b == 0) {
    return a;
  }
  return gcd(b, a % b);
}

std::pair<int, int> extended_gcd(int a, int b) {
  int x = 0, y = 1;
  int u = 1, v = 0;
  int m, n, q, r;
  while (a != 0) {
    q = b / a; r = b % a;
    m = (x - u * q), n = (y - v * q);
    b = a; a = r;
    x = u; y = v;
    u = m; v = n;
  }
  // assert(b == 1);
  return std::make_pair(x, y);
}

i32 mod_inv(i32 a, i32 mod) {
  int ret = extended_gcd(a, mod).first;
  if (ret < 0) {
    return ret + mod;
  } else {
    return ret;
  }
}

i32 cr_extend(std::vector<std::pair<i32, i32>> a, i32 mod) {
  i128 val = 0, prod_all = 1;
  for (int i = 0; i < a.size(); i++) {
    i32 prod = 1;
    i128 prod_m = 1;
    for (int j = 0; j < a.size(); j++) {
      if (i != j) {
        prod = ((i64) prod * a[j].second) % a[i].second;
        prod_m = (i128) prod_m * a[j].second;
      }
    }
    prod_m = (i128) prod_m * mod_inv(prod, a[i].second);
    val += (i128) prod_m * a[i].first;
    prod_all = (i128) prod_all * a[i].second;
  }
  return (val % prod_all) % mod;
}

std::vector<i32> dft(
  std::vector<i32>&& a,
  i32 p,  // Prime < 2^30
  i32 pr, // Radix of 2^k in Z_p
  int16_t n) {
    // assert((p - 1) & ((1 << (n + 1)) - 1));
    // assert(a.size() == (1 << n));

    i32 x, y, x1, y1, r, w;
    std::vector<i32>::iterator s, t;
    std::vector<i32> wvals(1 << n, 1);
    for (int i = 1; !(i >> n); i++) {
      wvals[i] = ((i64) wvals[i - 1] * pr) % p;
    }
    // assert(wvals[1 << (n - 1)] == p - 1);
    for (int i = 0; i < n; i++) {
      int i2 = n - i - 1, i2_p = 1 << i2;
      s = a.begin();
      t = a.begin() + i2_p;
      for (int j = 0; !(j >> i); j++) {
        for (int k = 0; !(k >> i2); k++) {
          w = wvals[k << i];
          x = *s;
          y = *t;
          x1 = x + y;
          if (x1 > p) {
            x1 -= p;
          }
          r = x + p - y;
          y1 = ((i64) w * r) % p;
          *s = x1;
          *t = y1;
          s++;
          t++;
        }
        s += i2_p;
        t += i2_p;
      }
    }
    std::vector<i32> ret(1 << n);
    for (int i = 0; !(i >> n); i++) {
      ret[bit_rev_table[i] >> (32 - n)] = a[i];
    }
    return ret;
}

std::vector<i32> convolute_h(
  std::vector<i32>&& a,
  std::vector<i32>&& b,
  i32 p,  // Prime
  i32 pr, // Radix of 2^k in Z_p
  int16_t n) {
    // assert((p - 1) & ((1 << (n + 1)) - 1));
    // assert(!((a.size() + b.size() - 1) >> n));
    a.resize(1 << n, 0);
    b.resize(1 << n, 0);
    std::vector<i32> x = dft(std::move(a), p, pr, n);
    std::vector<i32> y = dft(std::move(b), p, pr, n);

    for (int i = 0; !(i >> n); i++) {
      x[i] = ((i64) x[i] * y[i]) % p;
    }
    y.clear();

    i32 v2 = mod_inv(1 << n, p);

    std::vector<i32> z = dft(std::move(x), p, mod_inv(pr, p), n);
    for (int i = 0; !(i >> n); i++) {
      z[i] = ((i64) z[i] * v2) % p;
    }
    return z;
}

std::vector<i32> convolute(
  std::vector<i32>&& a,
  std::vector<i32>&& b,
  i32 mod) {
    if (a.size() + b.size() - 1 < 128) {
      std::vector<i32> c(a.size() + b.size() - 1, 0);
      int i = 0, j = 0;
      for (const auto& av: a) {
        for (const auto& bv: b) {
          c[i + j] += ((i64) av * bv) % mod;
          c[i + j] = c[i + j] > mod ? c[i + j] - mod : c[i + j];
          j++;
        }
        j = 0;
        i++;
      }
      return c;
    }
    int16_t n = log2(a.size() + b.size() - 1) + 1;
    double t = 2 * log2(mod) + log2(a.size() + b.size()) - 1;

    std::vector<i32> primes, prs;
    i32 p = (1 << 30) + (1 << n) + 1;
    while (t > 0) {
      do {
        p -= (1 << (n + 1));
      } while (!is_prime(p));
      t -= log2(p);
      primes.push_back(p);
      prs.push_back(mod_pow(primitive_root(p), (p - 1) >> n, p));
    }
    std::vector<std::vector<i32>> convs;
    for (int i = 0; i < primes.size(); i++) {
      convs.push_back(convolute_h(
        std::vector<i32>(a),
        std::vector<i32>(b),
        primes[i],
        prs[i],
        n));
    }

    std::vector<i32> ret(1 << n);
    for (int i = 0; !(i >> n); i++) {
      std::vector<std::pair<i32, i32>> coeff(primes.size());
      for (int j = 0; j < primes.size(); j++) {
        coeff[j] = std::make_pair(convs[j][i], primes[j]);
      }
      ret[i] = cr_extend(coeff, mod);
    }
    return ret;
}



i32 mod_p = 1000000007;
std::vector<i32> mod_p_invs(3000000);

void test() {
  std::vector<i32> p1, p2, p3, p4;
  for (int i = 0; i < 1000000; i++) {
    p1.push_back((749784329 + i * 329743892 + i * i * 297483297) % mod_p);
    p2.push_back((758943279 + i * 578943784 + i * i * 438217433) % mod_p);
  }
  // p3 = std::vector<i32>(p1.size() + p2.size() - 1, 0);
  // for (int i = 0; i < 1000000; i++) {
  //   for (int j = 0; i + j < 100; j++) {
  //     p3[i + j] += ((i64) p1[i] * (i64) p2[j]) % mod_p;
  //     p3[i + j] = p3[i + j] > mod_p ? p3[i + j] - mod_p : p3[i + j];
  //   }
  // }
  p4 = convolute(std::move(p1), std::move(p2), mod_p);
  // for (int i = 0; i < 100; i++)  {
  //   if (p3[i] != p4[i]) {
  //     printf(".");
  //   }
  // }
}

void init2() {
  mod_p_invs[1] = 1;
  for (int i = 2; i < 3000000; i++) {
    i32 r = mod_p_invs[mod_p % i], q = (mod_p / i);
    mod_p_invs[i] = mod_p - ((i64) q * r) % mod_p;
  }
}

std::vector<i32> stirling_2(i32 n) {
  std::vector<i32> fac_inv(n + 1), pw(n + 1);
  fac_inv[0] = 1;
  pw[0] = 0;
  for (int i = 1; i <= n; i++) {
    fac_inv[i] = ((i64) fac_inv[i - 1] * mod_p_invs[i]) % mod_p;
    pw[i] = ((i64) fac_inv[i] * mod_pow(i, n, mod_p)) % mod_p;
  }
  for (int i = 1; i <= n; i += 2) {
    fac_inv[i] = mod_p - fac_inv[i];
  }
  std::vector<i32> ret = convolute(std::move(fac_inv), std::move(pw), mod_p);
  ret.resize(n + 1);

  return ret;
}

std::vector<i32> bin_coeffs(i32 n, i32 a, i32 l) {
  i32 ainv = mod_inv(a, mod_p), val = mod_pow(a, n, mod_p);
  std::vector<i32> ret(l);
  for (int i = 0; i < l; i++) {
    ret[i] = val;
    val = ((i64) val * ainv) % mod_p;
    val = ((i64) val * (n - i)) % mod_p;
    val = ((i64) val * mod_p_invs[i + 1]) % mod_p;
  }
  return ret;
}

i32 summatory(i64 N, i64 a, i64 r) {
  if (r == 0) {
    if (a == 1) {
      return N;
    } else {
      return (
        (i64) (mod_pow(a, N + 1, mod_p) - 1) * mod_inv(a - 1, mod_p) - 1
      ) % mod_p;
    }
  }
  std::vector<i32> coeff;
  if (a == 1) {
    coeff = bin_coeffs(N + 1, 1, r + 2);
    coeff.erase(coeff.begin());
  } else {
    i32 a1inv = mod_inv(a - 1, mod_p);
    std::vector<i32> v1 = bin_coeffs(N + 1, a, r + 1), v2(r + 1);
    v1[0] -= 1;
    v2[0] = a1inv;
    for (int i = 1; i <= r; i++) {
      v2[i] = ((i64) v2[i - 1] * a1inv) % mod_p;
      v2[i] = mod_p - v2[i];
    }
    coeff = convolute(std::move(v1), std::move(v2), mod_p);
    coeff.resize(r + 1);
  }
  return 0;
  std::vector<i32> c2 = stirling_2(r);
  i32 fact_v = 1, pow_v = 1, ret = 0, val = 1;
  for (int i = 0; i <= r; i++) {

    val = ((i64) c2[i] * coeff[i]) % mod_p;
    val = ((i64) val * pow_v) % mod_p;
    val = ((i64) val * fact_v) % mod_p;
    ret += val;
    ret = ret > mod_p ? ret - mod_p : ret;

    pow_v = ((i64) a * pow_v) % mod_p;
    fact_v = ((i64) (i + 1) * fact_v) % mod_p;
  }
  return ret;
}

int main() {
  init();
  // init2();
  test();
  return 0;
  int y = 1;
  while (y--) {
    summatory(123456789, 214365879, 1100000 - y);
  }
  return 0;

  i32 t, N, a, r;
  scanf("%u", &t);
  while (t--) {
    scanf("%u %u %u", &N, &a, &r);
    printf("%u\n", summatory(N, a, r));
  }
  return 0;
}
