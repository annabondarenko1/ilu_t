#include <cmath>
#include <limits>
#include <vector>
#include <unordered_map>
#include <thread>
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
    
    
    for (int i = 0; i < len; ++i) {
      if (values.find(i) == values.end()) {
        continue;
      }

      values[i] /= k;
    }
  }

  // Check
  void operator-=(const SparseVector& vec) {
    
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
      auto start{std::chrono::steady_clock::now()};
      auto n = Get2Norm();
      auto end{std::chrono::steady_clock::now()};
    // std::cout << "Norm time (ns) = "
              // << std::chrono::duration_cast<std::chrono::nanoseconds>(end-start).count() << std::endl;

    if (std::abs(it->second - n*t) < kEps) {
      values.erase(it);
    }


  }


  // double Get2Norm() {
  //   auto threads_num = std::thread::hardware_concurrency();
  //   std::vector<double> norms(threads_num, 0.0);
  //   std::vector<std::thread> threads{};
  //   double norm = 0;
  //   for (int i = 0; i < threads_num; ++i) {
  //     threads.emplace_back(&SparseVector::GetSubNorm, *this, i, threads_num, std::ref(norms[i]));
  //   }

  //   for (int i = 0; i < threads_num; ++i) {
  //     norm += norms[i];
  //   }

  //   std::cout << 1 << std::endl;
  //   return std::sqrt(norm);
  // }

  void GetSubNorm(int idx, int threads_num, std::vector<double>& norms) {
    auto start{std::chrono::steady_clock::now()};
    int l = (len/threads_num)*idx;
    int r = idx == (threads_num-1) ? len: (len/threads_num)*(idx+1);

    for (int i = l; i < r; i += threads_num) {
      norms[idx] += values[i]*values[i];
    }
    auto end{std::chrono::steady_clock::now()};


  }

  // double Get2Norm() {
  //   double norm = 0;

  //   int threads_num = 4;
  //   auto start{std::chrono::steady_clock::now()};
  //   std::vector<std::thread> threads{};
  //   std::vector<double> norms(threads_num);

  //   for (int i = 0; i < threads_num; ++i) {
  //     threads.emplace_back(
  //       &SparseVector::GetSubNorm, *this, i, threads_num, std::ref(norms));
  //   }
  //   auto end{std::chrono::steady_clock::now()};
  //   std::cout << "GetSubNorm (nsec) = " << std::chrono::duration_cast<std::chrono::nanoseconds>(end-start).count() << std::endl;


  //   for (int i = 0; i < threads_num; ++i) {
  //     norm += norms[i];
  //   }

  //   for (int i = 0; i < threads_num; ++i) {
  //     threads[i].join();
  //   }
  //   auto res = std::sqrt(norm);

  //   return res;
  // }

  double Get2Norm() {
    double norm = 0;

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

    SparseVector(const double *vec, int len) {
    for (int i = 0; i < len; ++i) {
      if (std::abs(vec[i]) >= kEps) {
        values[i] = vec[i];
      }
    }
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

  void SparseRows(Matrix &mtrx, int thread_idx, int thread_num,
                  std::vector<SparseVector> &thread_results) {

    for (int row_idx = thread_idx; row_idx < mtrx.m_rows;
          row_idx += thread_num) {
      thread_results[row_idx] =
          SparseVector(mtrx[row_idx], mtrx.m_cols);
    }
  }

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
  CSRMatrix(Matrix& mtrx) {
    auto start = std::chrono::steady_clock::now();
    std::vector<SparseVector> thread_results(mtrx.m_rows);
    std::vector<std::thread> threads{};
    unsigned int thread_num{std::thread::hardware_concurrency()};

    for (int i = 0; i < thread_num; ++i) {
      threads.emplace_back(
          &CSRMatrix::SparseRows, *this, std::ref(mtrx), i, thread_num, std::ref(thread_results));
    }

    for (int i = 0; i < thread_num; ++i) {
      threads[i].join();
    }
    for (int i = 0; i < mtrx.m_rows; ++i) {
      AppendRow(thread_results[i]);
    }
    auto end = std::chrono::steady_clock::now();
    std::cout << "CSR format time (ms) = " << std::chrono::duration_cast<std::chrono::milliseconds>(end -start).count() << std::endl;
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