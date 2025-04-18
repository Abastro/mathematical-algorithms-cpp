#include <algorithm>
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

template <typename T> ostream &operator<<(ostream &os, const span<T> &v) {
  for (auto x : v) {
    os << x << " ";
  }
  return os;
}

vector<int> random_int_vector(size_t size) {
  auto result = vector<int>();
  for (size_t i = 0; i < size; i++) {
    result.push_back(rand() % 100);
  }
  return result;
}

vector<Real> random_real_vector(size_t size) {
  auto result = vector<Real>();
  for (size_t i = 0; i < size; i++) {
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

// Add polynomials in-place assuming the size is allocated
template <typename R>
void add_inplace(span<R> &a, span<R> &b, span<R> &result) {
  for (size_t i = 0; i < a.size(); i++)
    result[i] = a[i] + b[i];
}

// Subtract polynomials in-place assuming the size is allocated
template <typename R>
void sub_inplace(span<R> &a, span<R> &b, span<R> &result) {
  for (size_t i = 0; i < a.size(); i++)
    result[i] = a[i] - b[i];
}

template <typename R> vector<R> poly_normalize(vector<R> &a) {
  int i;
  for (i = a.size() - 1; i >= 0; i--) {
    if (a[i] != 0)
      break;
  }
  return vector(a.begin(), a.begin() + i + 1);
}

// Computes basic multiplication, assuming result is initialized as 0 and has enough space
template <typename R> void poly_mult_basic_span(span<R> &a, span<R> &b, span<R> &result) {
  for (size_t i = 0; i < a.size(); i++) {
    // Add at i-th position
    auto at_ith = result.subspan(i, b.size());
    for (size_t j = 0; j < b.size(); j++) {
      at_ith[j] = at_ith[j] + a[i] * b[j];
    }
  }
}

// Basic polynomial multiplication.
template <typename R> vector<R> poly_mult_basic(vector<R> &a, vector<R> &b) {
  if (a.empty() && b.empty())
    return vector<R>(0);

  vector<R> result = vector<R> (a.size() + b.size() - 1);
  auto span_a = span(a);
  auto span_b = span(b);
  auto span_result = span(result);
  poly_mult_basic_span(span_a, span_b, span_result);

  return result;
}

#define THRESHOLD 16

/**
 * A step of the Karatsuba function.
 *
 * NOTE: interestingly, the basic case is quite a performance bottleneck.
 * Hence, the basic case needs to be implemented well.
 *
 * @param deg_bnd power-of-2 degree bound
 * @param buffer the buffer which is used only throughout the invocation
 */
template <typename R>
void poly_mult_Karatsuba_step(const size_t deg_bnd, span<R> &a, span<R> &b,
                              span<R> &result, span<R> &buffer) {
  // Result may be reused, so this needs clearing
  for (auto &entry : result)
    entry = 0;

  if (deg_bnd <= THRESHOLD) {
    poly_mult_basic_span(a, b, result);
    return;
  }

  const auto next_bnd = deg_bnd >> 1;
  auto a0 = a.subspan(0, next_bnd);
  auto a1 = a.subspan(next_bnd, next_bnd);
  auto b0 = b.subspan(0, next_bnd);
  auto b1 = b.subspan(next_bnd, next_bnd);

  auto a01 = buffer.subspan(0, next_bnd);
  auto b01 = buffer.subspan(next_bnd, next_bnd);

  auto prod0 = result.subspan(0, deg_bnd);
  auto prod1 = result.subspan(deg_bnd, deg_bnd);
  auto prod_add = buffer.subspan(next_bnd * 2, next_bnd * 2);

  auto next_buffer = buffer.subspan(next_bnd * 4, next_bnd * 4);

  // correctly put into prod0 and prod1 position
  poly_mult_Karatsuba_step(next_bnd, a0, b0, prod0, next_buffer);
  poly_mult_Karatsuba_step(next_bnd, a1, b1, prod1, next_buffer);

  add_inplace(a0, a1, a01);
  add_inplace(b0, b1, b01);
  poly_mult_Karatsuba_step(next_bnd, a01, b01, prod_add, next_buffer);

  // adjust prod_add
  sub_inplace(prod_add, prod0, prod_add);
  sub_inplace(prod_add, prod1, prod_add);

  // Add middle term at X^next_bnd position
  auto result_mid = result.subspan(next_bnd, deg_bnd);
  add_inplace(prod_add, result_mid, result_mid);
}

template <typename R>
vector<R> poly_mult_Karatsuba(vector<R> &a, vector<R> &b) {
  size_t deg_bound = 1;
  while (deg_bound < max(a.size(), b.size()))
    deg_bound = deg_bound << 1;

  a.resize(deg_bound);
  b.resize(deg_bound);
  auto result = vector<R>(deg_bound << 1);
  auto buffer = vector<R>(deg_bound * 4);

  auto span_a = span(a);
  auto span_b = span(b);
  auto span_result = span(result);
  auto span_buffer = span(buffer);
  poly_mult_Karatsuba_step(deg_bound, span_a, span_b, span_result, span_buffer);
  return result;
}

void basic_vs_Karatsuba(size_t size) {
  auto p = random_int_vector(size);
  auto q = random_int_vector(size);

  cout << "Degree " << size - 1 << endl;

  auto begin = chrono::high_resolution_clock::now();
  auto basic = poly_mult_basic(p, q);
  auto end = chrono::high_resolution_clock::now();
  auto spent = chrono::duration<double>(end - begin);
  cout << "Basic took " << spent.count() << "s" << endl;

  begin = chrono::high_resolution_clock::now();
  auto karat = poly_mult_Karatsuba(p, q);
  end = chrono::high_resolution_clock::now();
  spent = chrono::duration<double>(end - begin);
  cout << "Karatsuba took " << spent.count() << "s" << endl;

  if (poly_normalize(basic) == poly_normalize(karat))
    cout << "+ Match" << endl;
  else
    cout << "- Mismatch" << endl;
  cout << endl;
}

void only_Karatsuba(size_t size) {
  auto p = random_int_vector(size);
  auto q = random_int_vector(size);

  cout << "Degree " << size - 1 << endl;

  auto begin = chrono::high_resolution_clock::now();
  poly_mult_Karatsuba(p, q);
  auto end = chrono::high_resolution_clock::now();
  auto spent = chrono::duration<double>(end - begin);
  cout << "Karatsuba took " << spent.count() << "s" << endl;
}

int main() {
  {
    auto p = vector<int>{1, 2};
    auto q = vector<int>{3, 4, 5};
    cout << "P: " << p << endl;
    cout << "Q: " << q << endl;
    cout << "P + Q: " << poly_add(p, q) << endl;
    cout << "basic P * Q: " << poly_mult_basic(p, q) << endl;
    auto karat = poly_mult_Karatsuba(p, q);
    cout << "Karatsuba P * Q: " << karat << endl;
    cout << "normalized: " << poly_normalize(karat) << endl;
    cout << endl;
  }

  {
    auto p = random_int_vector(6);
    auto q = random_int_vector(8);
    cout << "P: " << p << endl;
    cout << "Q: " << q << endl;
    auto basic = poly_mult_basic(p, q);
    cout << "basic P * Q: " << basic << endl;
    auto karat = poly_mult_Karatsuba(p, q);
    cout << "Karatsuba P * Q: " << karat << endl;
    cout << endl;
  }

  basic_vs_Karatsuba(128);
  basic_vs_Karatsuba(256);
  basic_vs_Karatsuba(512);
  basic_vs_Karatsuba(1024);
  basic_vs_Karatsuba(2048);
  basic_vs_Karatsuba(4096);
  basic_vs_Karatsuba(8192);
  basic_vs_Karatsuba(16384);

  only_Karatsuba(32768);
  only_Karatsuba(65536);
  only_Karatsuba(131072);

  only_Karatsuba(1 << 20);

  return 0;
}
