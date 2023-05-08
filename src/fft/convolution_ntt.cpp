#include <iostream>
#include <vector>
#include <cmath>
typedef __int128 int128_t;
typedef long long int64_t;
typedef unsigned int uint32_t;
typedef unsigned long long uint64_t;
#define N_LIM 20
uint32_t REV[1 << N_LIM];


#define P_BITS 32
static inline uint32_t Redc(uint64_t n, uint32_t mod, uint32_t minv) {
    uint32_t x = n * minv;
    x = (n + (uint64_t)mod * x) >> P_BITS;
    return x < mod ? x : x - mod;
}

constexpr uint32_t PrepRedc(uint32_t P) {
    uint64_t minv = -P;
    for (int b = 3; b < P_BITS; b *= 2) {
      minv = 2*minv + P * minv * minv;
    }

    return minv;
}

template<uint32_t P>
class ZP {
  uint32_t x;

  static constexpr uint32_t MOD = P;
  static constexpr uint32_t MINV = PrepRedc(P);

  ZP(uint32_t x_) {
    x = x_;
  }
  public:
  ZP() {
    x = 0;
  }

  static ZP unsafe(uint32_t x_) {
    return ZP(x_);
  }
  static ZP from_int(uint32_t x_) {
    return ZP(((uint64_t) x_ << 32) % P);
  }

  uint32_t to_int() {
    return Redc(x, P, MINV);
  }

  inline ZP operator*(const ZP& b) const {
    return ZP(Redc((uint64_t) x * b.x, P, MINV));
  }

  inline void operator*=(const ZP& b) {
    x = Redc((uint64_t) x * b.x, P, MINV);
  }

  inline ZP operator+(const ZP& b) const {
    uint32_t v = x + b.x;
    return ZP(v > P ? v - P : v);
  }

  inline void operator+=(const ZP& b) {
    x += b.x;
    x = x > P ? x - P : x;
  }

  inline ZP operator-(const ZP& b) const {
    uint32_t v = x + (P - b.x);
    return ZP(v > P ? v - P : v);
  }

  inline void operator-=(const ZP& b) {
    x += P - b.x;
    x = x > P ? x - P : x;
  }
};

template<typename T>
T pow_(T a, uint32_t r) {
  if (r == 0) {
    return T::from_int(1);
  }
  T val = pow_(a, r >> 1);
  val = val * val;
  if (r & 1) {
    return val * a;
  } else {
    return val;
  }
}

template<uint32_t P>
std::vector<ZP<P>> init_w(ZP<P> root, int L) {
  for (int r = 0; !(r >> L); r++) {
    REV[r] = (REV[r >> 1] >> 1) | ((r & 1) << (L - 1));
  }
  auto w = pow_(root, (P - 1) >> L);
  std::vector<ZP<P>> W(1 << L);
  W[0] = ZP<P>::from_int(1);
  for (int i = 1; !(i >> L); i++) {
    W[i] = W[i - 1] * w;
  }

  return W;
}

template<typename T>
void inv_w(std::vector<T>& W, int L) {
  for (int i = 1; !(i >> (L - 1)); i++) {
    std::swap(W[i], W[(1 << L) - i]);
  }
}


uint32_t mod_pow(uint32_t a, uint32_t r, uint32_t p) {
  if (r == 0) return 1;
  if (r & 1) return ((uint64_t) a * mod_pow(a, r - 1, p)) % p;
  return mod_pow((uint64_t) a * a % p, r, p);
}


std::pair<int64_t, int64_t> extended_gcd(int64_t a, int64_t b) {
  int64_t x = 0, y = 1;
  int64_t u = 1, v = 0;
  int64_t m, n, q, r;
  while (a) {
    q = b / a; r = b % a;
    m = (x - u * q), n = (y - v * q);
    b = a; a = r;
    x = u; y = v;
    u = m; v = n;
  }
  // b == gcd(a_init, b_init)
  return std::make_pair(x, y);
}

uint32_t mod_inv(uint32_t a, uint32_t mod) {
  int ret = extended_gcd(a, mod).first;
  if (ret < 0) {
    return ret + mod;
  } else {
    return ret;
  }
}

template<typename T>
void ntt(std::vector<T>& A, const std::vector<T> W,int L) {
  for (int i = 0; !(i >> L); i++) {
    if (i > REV[i]) {
      std::swap(A[i], A[REV[i]]);
    }
  }

  for (int s = 0; s < L; s++) {
    int m = 1 << s, m2 = m << 1, rm2 = 1 << (L - s - 1);
    for (int k = 0; !(k >> L); k += m2) {
      for (int j = 0; !(j >> s); j++) {
        auto &As = A[k + j], &At = A[k + j + m];
        auto t = W[rm2 * j] * At;
        At = As - t;
        As = As + t;
      }
    }
  }
}

template<typename T>
std::vector<T> convolute_h(
  std::vector<T>&& a,
  std::vector<T>&& b,
  std::vector<T>&& W,
  const int16_t n) {
    // assert((p - 1) & ((1 << (n + 1)) - 1));
    // assert(!((a.size() + b.size() - 1) >> n));

    a.resize(1 << n, T());
    b.resize(1 << n, T());

    ntt(a, W, n);
    ntt(b, W, n);

    for (int i = 0; !(i >> n); i++) {
      a[i] = a[i] * b[i];
    }
    inv_w(W, n);
    ntt(a, W, n);
    auto v2 = T::unsafe(1 << (32 - n));
    for (int i = 0; !(i >> n); i++) {
      a[i] = a[i] * v2;
    }
    return std::move(a);
}

template<typename T>
std::vector<T> convolute_bf(
  std::vector<T>&& a,
  std::vector<T>&& b) {

  std::vector<T> c(a.size() + b.size() - 1);
  int i, j;

  i = 0;
  for (const auto& av: a) {
    j = 0;
    for (const auto& bv: b) {
      c[i + j] += av * bv;
      j++;
    }
    i++;
  }
  return c;
}

const uint32_t MOD = 1000000007;
std::vector<ZP<MOD>> mod_p_invs(2000000);

template<typename T>
std::vector<T> convolute(
  std::vector<T>&& a,
  std::vector<T>&& b) {
    // if (a.size() < 64 || b.size() < 64) {
    //   return convolute_bf(std::move(a), std::move(b), MOD);
    // }
    int16_t n = log2(a.size() + b.size() - 1) + 1;

    // pr < 2^31, 2^25 | pr - 1
    std::vector<uint32_t> primes = {{1107296257, 1711276033, 1811939329}};
    std::vector<uint32_t> primitive_roots = {{10, 29, 13}};

    std::vector<ZP<1107296257>> A1(a.size()), B1(b.size());
    std::vector<ZP<1711276033>> A2(a.size()), B2(b.size());
    std::vector<ZP<1811939329>> A3(a.size()), B3(b.size());

    for (int i = 0; i < a.size(); i++) {
      uint32_t v = a[i].to_int();
      A1[i] = ZP<1107296257>::from_int(v);
      A2[i] = ZP<1711276033>::from_int(v);
      A3[i] = ZP<1811939329>::from_int(v);
    }
    // return std::vector<uint32_t>(0);

    for (int i = 0; i < b.size(); i++) {
      uint32_t v = b[i].to_int();
      B1[i] = ZP<1107296257>::from_int(v);
      B2[i] = ZP<1711276033>::from_int(v);
      B3[i] = ZP<1811939329>::from_int(v);
    }
    // auto O1 = convolute_bf(std::vector<ZP<1107296257>>(A1), std::vector<ZP<1107296257>>(B1));
    // auto O2 = convolute_bf(std::vector<ZP<1711276033>>(A2), std::vector<ZP<1711276033>>(B2));
    // auto O3 = convolute_bf(std::vector<ZP<1811939329>>(A3), std::vector<ZP<1811939329>>(B3));

    auto R1 = convolute_h(std::move(A1), std::move(B1), init_w(ZP<1107296257>::from_int(10), n), n);
    auto R2 = convolute_h(std::move(A2), std::move(B2), init_w(ZP<1711276033>::from_int(29), n), n);
    auto R3 = convolute_h(std::move(A3), std::move(B3), init_w(ZP<1811939329>::from_int(13), n), n);

    std::vector<T> ret(1 << n);
    std::vector<int128_t> C = {{
      ((int128_t)1711276033 * 1811939329) * 448191345,
      ((int128_t)1107296257 * 1811939329) * 285212624,
      ((int128_t)1711276033 * 1107296257) * 776545473
    }};
    int128_t C0 = ((int128_t)1107296257 * 1711276033) * 1811939329;
    for (int i = 0; !(i >> n); i++) {
      // if (i < 10) std::cout << "(" <<
      //   R1[i].to_int() - O1[i].to_int() << " " <<
      //   R2[i].to_int() - O2[i].to_int() << " " <<
      //   R3[i].to_int() - O3[i].to_int() << ")";
      ret[i] = T::from_int(((C[0] * R1[i].to_int() + C[1] * R2[i].to_int() + C[2] * R3[i].to_int()) % C0) % T::MOD);
    }
    ret.resize(a.size() + b.size() - 1);
    return ret;
}

void init() {
  std::vector<uint32_t> temp(mod_p_invs.size());
  temp[1] = 1;
  mod_p_invs[1] = ZP<MOD>::from_int(1);
  for (int i = 2; i < mod_p_invs.size(); i++) {
    temp[i] = MOD - ((uint64_t) (MOD / i) * temp[MOD % i]) % MOD;
    mod_p_invs[i] = ZP<MOD>::from_int(temp[i]);
  }
}

std::vector<ZP<MOD>> stirling_2(uint32_t n) {
  std::vector<ZP<MOD>> fac_inv(n + 1), pw(n + 1);
  fac_inv[0] = ZP<MOD>::from_int(1);
  pw[0] = ZP<MOD>();
  for (int i = 1; i <= n; i++) {
    fac_inv[i] = fac_inv[i - 1] * mod_p_invs[i];
    pw[i] = fac_inv[i] * pow_(ZP<MOD>::from_int(i), n);
  }
  for (int i = 1; i <= n; i += 2) {
    fac_inv[i] = ZP<MOD>() - fac_inv[i];
  }
  auto ret = convolute(std::move(fac_inv), std::move(pw));
  ret.resize(n + 1);
  return ret;
}

// std::vector<uint32_t> bin_coeffs(uint32_t n, uint32_t a, uint32_t l) {
//   uint32_t ainv = mod_inv(a, MOD), val = mod_pow(a, n, MOD);
//   std::vector<uint32_t> ret(l);
//   for (int i = 0; i < l; i++) {
//     ret[i] = val;
//     val = ((uint64_t) val * ainv) % MOD;
//     val = ((uint64_t) val * (n - i)) % MOD;
//     val = ((uint64_t) val * mod_p_invs[i + 1]) % MOD;
//   }
//   return ret;
// }

int main() {
  init();
  //
  while (1) {
  auto res = stirling_2(200000);
  for (int i = 0; i < 100; i++) {
    std::cout << res[i].to_int() << " ";
  }
  std::cout << std::endl;
  return 0;
  }
  return 0;
}
