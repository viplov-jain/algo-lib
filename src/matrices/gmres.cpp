#include <stdio.h>
#include <math.h>
#include <vector>
#include <unordered_map>

#define EPS 1e-14
#define AA(i, j) A[N * (i) + (j)]
#define BB(i, j) B[N * (i) + (j)]
#define CC(i, j) C[N * (i) + (j)]

template <typename T>
void matrmul(T *A, T *B, T *C, const int N) {
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      CC(i, j) = 0;
    }
  }
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      for (int k = 0; k < N; k++) {
        CC(i, k) += AA(i, j) * BB(j, k);
      }
    }
  }
}

void printmat(double *A, const int N) {
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      printf("%16.12f ", AA(i, j));
    }
    printf("\n");
  }
  printf("\n");
}

template <typename T>
void swap(int i, int j, T *A, T *B, const int N) {
  if (i == j) return;
  for (int k = 0; k < N; k++) {
    std::swap(AA(i, k), AA(j, k));
  }
  std::swap(B[i], B[j]);
}

template <typename T>
void scale(const int i, const T a, T *A, T *B, const int N) {
  if (abs(a - 1) < EPS) return;
  for (int k = 0; k < N; k++) {
    AA(i, k) *= a;
  }
  B[i] *= a;
}

template <typename T>
void add(const int i, const int j, const T a, T *A, T *B, const int N) {
  if (abs(a) < EPS) return;
  for (int k = 0; k < N; k++) {
    AA(j, k) += a * AA(i, k);
  }
  B[j] += a * B[i];
}

// A X = (1, 0, 0...)'

#define HH(i, j) H[(i) * (N + 1) + (j)]
#define QQ(i, j) Q[(i) * (N + 1) + (j)]

template <typename T>
void solve_hess(T *H, T* X, int k, int N) {
  // Something wrong here
  printf("\n[[[[--------]]]]\n");
  for (int i = 0; i < k; i++) {
    printf("{{%f}}", X[i]);
  }
  printf("\n[[[[--------]]]]\n");
  std::vector<double> Q((N + 1) * (N + 1)), test(N + 1, 0.0);
  for (int i = 0; i < k; i++) {
    for (int j = 0; j < k; j++) {
      QQ(i, j) = HH(i, j);
    }
  }
  for (int i = 0; i < k; i++) {
    printf("%16.12f ", X[i]);
  }
  printf("\n");
  for (int i = 0; i < k; i++) {
    const double h = 1.0 / hypot(HH(i, i), HH(i + 1, i));
    const double c = HH(i, i) * h;
    const double s = HH(i + 1, i) * h;
    double a, b;
    for (int j = i; j < k; j++) {
      a = HH(i, j), b = HH(i + 1, j);
      HH(i    , j) = c * a + s * b;
      HH(i + 1, j) = s * a - c * b;
    }
    a = X[i]; b = X[i + 1];
    X[i    ] = c * a + s * b;
    X[i + 1] = s * a - c * b;
    const double r = HH(i, i);
    for (int j = i; j < k; j++) {
      HH(i, j) /= r;
    }
    printmat(H, N + 1);
    X[i] /= r;
    for (int i = 0; i < k; i++) {
      printf("%16.12f ", X[i]);
    }
    printf("\n");
  }

  for (int i = k - 1; i >= 0; i--) {
    for (int j = i - 1; j >= 0; j--) {
      X[j] -= HH(j, i) * X[i];
      HH(j, i) -= HH(j, i) * HH(i, i);
    }
    printmat(H, N + 1);
    for (int i = 0; i < k; i++) {
      printf("%16.12f ", X[i]);
    }
    printf("\n");
  }

  printf("\n[[[[--------]]]]\n");
  for (int i = 0; i < k; i++) {
    for (int j = 0; j < k; j++) {
      test[i] += QQ(i, j) * X[j];
      // printf("%f ", HH(i, j));
    }
    printf("((%f))", test[i]);
  }
  printf("\n[[[[--------]]]]\n");


}
// Numerical Implementations of the Generalized Minimal Residual Method (GMRES)
template <typename T>
void solve_eq(const T *A, T *X, const int N) {
  // printmat(A, N);
  for (int i = 0; i < N; i++) {
    X[i] = 10;
  }
  std::vector<std::vector<int>> places(N);
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      if (abs(AA(i, j)) > EPS) {
        places[i].push_back(j);
      }
    }
  }
  printf("\n\n---------\n\n");
  for (int iter = 0; iter < 200; iter++) {
    for (int i = 0; i < N; i++) {
      double val = 0.0;
      for (auto j: places[i]) {
        if (j == i) continue;
        val += AA(i, j) * X[j];
      }
      X[i] = (i == 0 ? 1.0 : 0.0) - val;
    }
  }
  printf("\n\n---------\n\n");

  for (int i = 0; i < N; i++) {
    printf("%f ", X[i]);
  }

  double test = 0.0;
  for (int i = 0; i < N; i++) {
    double val = i == 0 ? 1 : 0;
    for (auto j: places[i]) {
      val -= AA(i, j) * X[j];
    }
    test += val * val;
  }
  printf("Epsilon: %.18f\n", test);
  printf("\n\n---------\n\n");


  std::vector<double> B(N * N, 0), H((N + 1) * (N + 1), 0), Q((N + 1) * (N + 1), 0);
  std::vector<double> R0(N);
  double rnorm = 0.0;
  for (int i = 0; i < N; i++) {
    double val = (i == 0 ? 1 : 0);
    for (auto j: places[i]) {
      val -= AA(i, j) * X[j];
    }
    R0[i] = val;
    rnorm += val * val;
  }
  rnorm = sqrt(rnorm);
  if (rnorm < EPS) {
    return;
  }
  for (int i = 0; i < N; i++) {
    BB(i, 0) = R0[i] / rnorm;
  }

  int k_max = std::min(100, N - 1);
  std::vector<double> Vk(N), Hk;
  for (int k = 0; k < k_max; k++) {
    Vk = std::vector<double>(N, 0);
    for (int i = 0; i < N; i++) {
      for (auto j: places[k]) {
        Vk[i] += AA(i, j) * BB(j, k);
      }
    }
    Hk = std::vector<double>(k, 0);
    for (int j = 0; j < k; j++) {
      double val = 0;
      for (int i = 0; i < N; i++) {
        val += BB(i, j) * Vk[i];
      }
      Hk[j] = val;
    }
    for (int i = 0; i < N; i++) {
      double val = 0;
      for (int j = 0; j < k; j++) {
        val += BB(i, j) * Hk[j];
      }
      Vk[i] -= val;
    }
    double Vkn = 0;
    for (int i = 0; i < N; i++) {
      Vkn += Vk[i] * Vk[i];
    }
    Vkn = sqrt(Vkn);
    if (Vkn < EPS) {
      k_max = k;
      printf("FDSAFDSAFDSA\n\n");
      break;
    }
    for (int i = 0; i < N; i++) {
      Vk[i] /= Vkn;
    }
    for (int i = 0; i < N; i++) {
      BB(i, k + 1) = Vk[i];
    }
    for (int i = 0; i < k; i++) {
      HH(i, k) = Hk[i];
    }
    HH(k + 1, k) = Vkn;
  }
  for (int i = 0; i < k_max; i++) {
    R0[i] = i == 0 ? rnorm : 0;
  }
  solve_hess<double>(H.data(), R0.data(), k_max, N);

  for (int i = 0; i < N; i++) {
    for (int j = 0; j < k_max; j++) {
      X[i] = X[i] + BB(i, j) * R0[j];
    }
  }
  for (int i = 0; i < k_max; i++) {
    // printf("{{%f}} ", R0[i]);
  }

  double test2 = 0.0;
  for (int i = 0; i < N; i++) {
    double val = i == 0 ? 1 : 0;
    for (auto j: places[i]) {
      val -= AA(i, j) * X[j];
    }
    test2 += val * val;
  }
  printf("[[(%f)]]\n", test2);
  printf("\n\n---------\n\n");
}
