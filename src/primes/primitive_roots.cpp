#include <iostream>
#include <vector>
using namespace std;


int is_prime(int p) {
    for (int i = 2; i * i <= p; i++) {
        if (p % i == 0) return false;
    }
    return true;
}

int64_t powmod(const int64_t a, const int64_t k, const int64_t n) {
    if (!k) return 1;
    int64_t ret = powmod(a, k / 2, n);
    ret = (ret * ret) % n;
    return k & 1 ? (ret * a) % n : ret;
}

void get_factors(int n, vector<int>& F) {
    if (n == 1) return;
    int p = 2;
    while (n % p) {
        if (p * p > n) {F.push_back(n); return;}
        p++;
    }
    F.push_back(p);
    for (n /= p; n % p == 0; n /= p);
    get_factors(n, F);
}

int order(int x, int p, const vector<int>& F) {
   int res = p - 1;

   // try to remove factors from m until we can't remove any more
   for (const auto& q: F) {
      while (res % q == 0) {
         int n = res / q;
         if (powmod(x, n, p) != 1) break;
         res = n;
      }
   }
   return res;
}

int primitive_root(int p, const vector<int>& F) {
    if (p == 2) return 1;
    for (int g = 2;; g++) if (order(g, p, F) == p - 1) return g;
}

int primitive_root(int p) {
    vector<int> F;
    get_factors(p - 1, F);
    return primitive_root(p, F);
}


int main() {
    int p = (1 << 30) + 1;
    while (!(p >> 31)) {
        while (!is_prime(p)) {
            p += (1 << 25);
        }
        cout << p << " " << primitive_root(p) << endl;
        p += (1 << 25);
    }
}
