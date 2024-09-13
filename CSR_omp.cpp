#include <cmath>
#include <limits>
#include <vector>
#include <unordered_map>
#include <omp.h>

const double kEps{1e-6};

struct Matrix {
  Matrix(int rows, int cols) {
    m_mtrx = new double[rows*cols]{0};
    m_rows = rows;
    m_cols = cols;
  };
  Matrix(double* mtrx, int rows, int cols) {
    m_mtrx = mtrx;
    m_rows = rows;
    m_cols = cols
  };
  double* operator[](int i) {return m_mtrx + (m_rows * i);}

  void Print() {
    for (int i = 0; i < m_rows; ++i) {
      for (int j = 0; j < m_cols; ++j) {
        std::cout << m_mtrx[i*m_rows + j] << " ";
      }
      std::cout << std::endl;
    }
  }

  int m_rows;
  int m_cols;
  double* m_mtrx;
};

struct SparseVector {

  SparseVector() = default;

  // Check in-place
  void operator/=(double k) {
    
    #pragma omp parallel for
    for (int i = 0; i < len; ++i) {
      if (values.find(i) == values.end()) {
        continue;
      }

      values[i] /= k;
    }
  }

  // Check
  void operator-=(const SparseVector& vec) {
    #pragma omp parallel for
    for (int i = 0; i < len; ++i) {
      if (vec.values.find(i) == vec.values.end())
        continue;

      values[i] -= vec.values.at(i);
    }

  }

  void Add(int i, double val) {
    values[i] = val;
  }

  double operator[] (int k) {
    if (values.find(k) == values.end())
      return 0;

    return values[k];
  }

  void Drop(double t, int k) {
    auto it = values.find(k);
    if (it == values.end()) {
      return;
    }
      // auto start{std::chrono::steady_clock::now()};
      auto n = Get2Norm();
      // auto end{std::chrono::steady_clock::now()};
    // std::cout << "Norm time (ns) = "
              // << std::chrono::duration_cast<std::chrono::nanoseconds>(end-start).count() << std::endl;

    if (std::abs(it->second - n*t) < kEps) {
      values.erase(it);
    }


  }

  // double Get2Norm() {
  //   auto threads_num = omp_get_num_threads();
  //   std::vector<double> norms(threads_num, 0.0);

  //   #pragma omp parallel for
  //   for (int i = 0; i < len; ++i) {
  //     if (values.find(i) == values.end())
  //       continue;
    
  //     int thread_idx = omp_get_thread_num();
  //     norms[thread_idx] += values[i]*values[i];
  //   }

  //   double norm = 0;
  //   for (int i = 0; i < threads_num; ++i) {
  //     norm += norms[i];
  //   }

  //   return std::sqrt(norm);
  // }
  // double Get2Norm() {
  //   double norm = 0;
  //   std::vector<double> values_vec{};

  //   for (auto it = values.begin(); it != values.end(); ++it) {
  //     values_vec.emplace_back(it->second);
  //   }
  //   #pragma omp parallel for shared(values_vec) reduction(+: norm)
  //   for (int i = 0; i < values_vec.size(); ++i) {
  //     norm += values_vec[i]*values_vec[i];
  //   }

  //   return std::sqrt(norm);
  // }

  double Get2Norm() {
    double norm = 0;

    #pragma omp parallel reduction(+: norm)
    for (int i = 0; i < len; ++i) {
      if (values.find(i) == values.end()) {
        continue;
      }

      norm += values[i]*values[i];
    }

    return std::sqrt(norm);
  }

  SparseVector(const std::vector<int>& idxs_, const std::vector<double>& values_, int size) {

    for (int i = 0; i < idxs_.size(); ++i) {
      values[i] = values_[i];
    }

    len = size;
  }
  int len{};
  std::unordered_map<int, double> values;
};

SparseVector Multiply(double k, SparseVector vec) {
  for (auto it = vec.values.begin(); it != vec.values.end(); ++it) {
    it->second *= k;
  }

  return vec;
}

struct CSRMatrix {

  double Get(int i, int j) {
    int start = m_row_idx[i];
    int end = m_row_idx[i+1];

    for (int k = start; k < end; ++k) {
      if (m_col_idx[k] == j) {
        return m_values[k];
      }
    }

    return 0;
  }
  // можно распараллелить с std::threads
  CSRMatrix (Matrix mtrx) {
    int rows = mtrx.m_rows;
    int cols = mtrx.m_cols;

    std::vector<SparseVector> sparse_vectors(rows);
    // Acceleration - OK from 300ms to 150ms
    #pragma omp parallel for
    for (int i = 0; i < rows; ++i) {
      for (int j = 0; j < cols; ++j) {
        auto elem = mtrx[i][j];
        if (std::abs(elem) > kEps) {
          sparse_vectors[i].values[j] = elem;
        }
      }
    }

    for (int i = 0; i < rows; ++i) {
      AppendRow(sparse_vectors[i]);
    }

  }

  Matrix GetMatrix() {
    auto rows = GetRowCount();
    auto cols = GetColCount();
    Matrix res{rows, cols};

    // No acceleration (from 6(no omp) to 3 ms (omp) )
    for (int i = 0; i < rows; ++i) {
      auto start_idx = m_row_idx[i];
      auto stop_idx = m_row_idx[i+1];
      for (int j = start_idx; j < stop_idx; ++j) {
       // m_rows_idx[start : stop] contains idxs for v

        int row = i;
        int col = m_col_idx[j];
        auto value = m_values[j];

        res[row][col] = value;
      }
    }

    return res;
  }

  int GetRowCount() {
    return m_row_idx.size() - 1;
  }

  int GetColCount() {
    // Count of rows equals to count of cols
    return GetRowCount();
  }

  void Print() {
    std::cout << "m_col_idx = [ ";
    for (const auto& x: m_col_idx) {
      std::cout << x << " ";
    } 
    std::cout << "]\n";

    std::cout << "m_values = [ ";
    for (const auto& x: m_values) {
      std::cout << x << " ";
    }
    std::cout << "]\n";

    std::cout << "m_row_idx = [ ";
    for (const auto& x : m_row_idx) {
      std::cout << x << " ";
    }
    std::cout << "]\n";
  }

  SparseVector GetSparseRow(int i) {
    int start_idx = m_row_idx[i];
    int end_idx = m_row_idx[i+1];

    std::vector<double> values(end_idx-start_idx);
    std::vector<int> idxs(end_idx - start_idx);
    // NO acceleration
    for (int i = start_idx; i < end_idx; ++i) {
      values[i-start_idx] = m_values[i];
      idxs[i-start_idx] = m_col_idx[i];
    }
    
    return SparseVector{idxs, values, GetColCount()};
  }

  void AppendRow(const SparseVector& vec) {
    for (auto it = vec.values.begin(); it != vec.values.end(); ++it) {
      m_values.push_back(it->second);
      m_col_idx.push_back(it->first);
    }

    m_row_idx.push_back(m_row_idx.back() + vec.values.size());
   
  }
  std::vector<double> m_values;
  std::vector<int> m_col_idx;
  std::vector<int> m_row_idx = {0};
};