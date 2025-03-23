#include <chrono>
#include <complex>
#include <iostream>
#include <vector>
using namespace std;

typedef complex<double> Complex;
typedef double Real;

/**
 * Operator overloading for printing vectors.
 * @tparam T
 * @param os
 * @param v
 * @return
 */
template <typename T> ostream &operator<<(ostream &os, const vector<T> &v) {
  for (size_t i = 0; i < v.size(); i++) {
    os << v[i] << " ";
  }
  return os;
}

vector<Real> random_real_vector(int size) {
  auto result = vector<Real>();
  for (int i = 0; i < size; i++) {
    result.push_back((double)rand() / (double)RAND_MAX);
  }
  return result;
}

vector<Complex> vector_as_complex(vector<Real> &v) {
  auto result = vector<Complex>(v.size());
  for (size_t i = 0; i < v.size(); i++) {
    result[i] = Complex(v[i], 0);
  }
  return result;
}

template <typename R> vector<R> poly_add(vector<R> &a, vector<R> &b) {
  auto res = vector<R>(max(a.size(), b.size()));
  for (size_t i = 0; i < a.size(); i++) {
    res[i] = a[i];
  }
  for (size_t i = 0; i < b.size(); i++) {
    res[i] += b[i];
  }
  return res;
}

template <typename R> vector<R> poly_sub(vector<R> &a, vector<R> &b) {
  auto res = vector<R>(max(a.size(), b.size()));
  for (size_t i = 0; i < a.size(); i++) {
    res[i] = a[i];
  }
  for (size_t i = 0; i < b.size(); i++) {
    res[i] -= b[i];
  }
  return res;
}

// Shift up degrees by n
template <typename R> vector<R> poly_shift_up(int n, vector<R> &p) {
  auto result = p;
  for (int i = 0; i < n; i++) {
    result.insert(result.begin(), 0);
  }
  return result;
}

// Basic polynomial multiplication
template <typename R> vector<R> poly_mult_basic(vector<R> &a, vector<R> &b) {
  if (a.size() == 0 && b.size() == 0)
    return vector<R>(0);
  auto res = vector<R>(a.size() + b.size() - 1, 0);
  for (size_t i = 0; i < a.size(); i++) {
    // Start with i 0s
    auto tmp = vector<R>(i, 0);
    for (R bj : b) {
      tmp.push_back(a[i] * bj);
    }
    res = poly_add(res, tmp);
  }
  return res;
}

// TODO Reduce allocations

/**
 * A step of the Karatsuba function.
 * @param deg_bnd power-of-2 degree bound
 */
template <typename R>
vector<R> poly_mult_Karatsuba_step(const size_t deg_bnd, vector<R> &a,
                                   vector<R> &b) {
  if (deg_bnd <= 1)
    return poly_mult_basic(a, b);

  const auto next_bnd = deg_bnd >> 1;
  const auto next_bnd_in_a = min(next_bnd, a.size());
  const auto next_bnd_in_b = min(next_bnd, b.size());
  auto a0 = vector(a.begin(), a.begin() + next_bnd_in_a);
  auto a1 = vector(a.begin() + next_bnd_in_a, a.end());
  auto b0 = vector(b.begin(), b.begin() + next_bnd_in_b);
  auto b1 = vector(b.begin() + next_bnd_in_b, b.end());

  auto prod0 = poly_mult_Karatsuba_step(next_bnd, a0, b0);
  auto prod1 = poly_mult_Karatsuba_step(next_bnd, a1, b1);
  auto a01 = poly_add(a0, a1);
  auto b01 = poly_add(b0, b1);
  auto prod_add = poly_mult_Karatsuba_step(next_bnd, a01, b01);

  auto tmp1 = poly_sub(prod_add, prod0);
  auto tmp2 = poly_sub(tmp1, prod1);
  auto mid_term = poly_shift_up(next_bnd, tmp2);
  auto high_term = poly_shift_up(deg_bnd, prod1);
  auto higher = poly_add(mid_term, high_term);
  return poly_add(prod0, higher);
}

template <typename R>
vector<R> poly_mult_Karatsuba(vector<R> &a, vector<R> &b) {
  size_t deg_bound = 1;
  while (deg_bound < max(a.size(), b.size()))
    deg_bound = deg_bound << 1;

  return poly_mult_Karatsuba_step(deg_bound, a, b);
}

void basic_vs_Karatsuba(size_t size) {
  auto p = random_real_vector(size);
  auto q = random_real_vector(size);

  cout << "Degree " << size - 1 << endl;

  auto begin = chrono::high_resolution_clock::now();
  poly_mult_basic(p, q);
  auto end = chrono::high_resolution_clock::now();
  auto spent = chrono::duration<double>(end - begin);
  cout << "Basic took " << spent.count() << "s" << endl;

  begin = chrono::high_resolution_clock::now();
  poly_mult_Karatsuba(p, q);
  end = chrono::high_resolution_clock::now();
  spent = chrono::duration<double>(end - begin);
  cout << "Karatsuba took " << spent.count() << "s" << endl;
}

int main() {
  {
    auto p = vector<Real>{1.0, 2.0};
    auto q = vector<Real>{3.0, 4.0, 5.0};
    cout << "P: " << p << endl;
    cout << "Q: " << q << endl;
    cout << "P + Q: " << poly_add(p, q) << endl;
    cout << "basic P * Q: " << poly_mult_basic(p, q) << endl;
    cout << "Karatsuba P * Q: " << poly_mult_Karatsuba(p, q) << endl;
    cout << endl;
  }

  {
    auto p = random_real_vector(6);
    auto q = random_real_vector(8);
    cout << "P: " << p << endl;
    cout << "Q: " << q << endl;
    cout << "basic P * Q: " << poly_mult_basic(p, q) << endl;
    cout << "Karatsuba P * Q: " << poly_mult_Karatsuba(p, q) << endl;
  }

  basic_vs_Karatsuba(125);
  basic_vs_Karatsuba(250);
  basic_vs_Karatsuba(500);
  basic_vs_Karatsuba(1000);
  basic_vs_Karatsuba(2000);
  // Approaches at this point
  basic_vs_Karatsuba(4000);

  // {
  //   auto p = random_real_vector(4000);
  //   auto q = random_real_vector(4000);
  //   auto begin = chrono::high_resolution_clock::now();
  //   poly_mult_Karatsuba(p, q);
  //   auto end = chrono::high_resolution_clock::now();
  //   auto spent = chrono::duration<double>(end - begin);
  //   cout << "Karatsuba took " << spent.count() << "s" << endl;
  // }

  return 0;
}
