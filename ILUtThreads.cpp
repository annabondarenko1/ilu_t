#include <iostream>
#include <chrono>
#include <thread>

#include "CSR_threads.cpp"

const int N = 1000;
double t{0.01};
const double p{5};

void FillIn(Matrix& A, double p){
  std::srand(std::time(0));
  int size = A.m_rows;
  int nonZeroElements = size* size * p / 100;

  for (int i = 0; i < A.m_rows; ++i) {
    for (int j = 0; j < A.m_cols; ++j) {
      A[i][j] = 0;
    }
  }

    int count = 0;
    while (count < nonZeroElements) {
        int index1 = std::rand() % (size);
        int index2 = std::rand() % (size);

        if (A[index1][index2] == 0) {
            A[index1][index2] = std::rand() % 100 + 1; 
            count++;
        }
    }
}

int main() {
  Matrix A{N, N}, L{N, N}, U{N, N};
  auto start = std::chrono::steady_clock::now();
  FillIn(A, p);
  CSRMatrix l{L}, a{A}, u{U};

  
  for (int i = 1; i < N; ++i) {
    auto start = std::chrono::steady_clock::now();
    SparseVector w = a.GetSparseRow(i); // w = a_{i*}



    for (int k = 0; k < i; ++k) {
      if (std::abs(w[k]) < kEps) {
        continue;
      }

      w.values[k] = w[k]/a.Get(k, k); // w_k = w_k / a_{kk}
      w.Drop(t, k);
      if (std::abs(w[k]) >= kEps) {
        w -= Multiply(w[k], u.GetSparseRow(k));
      }
    }

    for (int k = 0; k < N; ++k) {
      w.Drop(t, k);
    }
    l.AppendRow(w);
    auto end = std::chrono::steady_clock::now();

    std::cout << "iteration " << i << " time (ms) = " << 
    std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count() << std::endl;
  }

  auto end = std::chrono::steady_clock::now();
  std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << std::endl;
  // L*U ~ A , LU != A
  // Распараллелить внешний и внутренний циклы через #pragma omp parallel
  // Распараллелить внешний и внутренний циклы через threads

  return 0;
}