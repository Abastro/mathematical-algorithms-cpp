#include <list>
#include <vector>
#include <version>

#include "types.hpp"

/**
 * Constraint: degrees should have the same size across the board
 */
template<typename R>
struct monomial {
  R coefficient;
  vector<size_t> degrees;
};

template <typename R>
std::ostream& operator<<(std::ostream& os, const monomial<R>& r) {
  os << r.coefficient;
  for (size_t i = 0; i < r.degrees.size(); i++) {
    os << " " << "X_" << i << "^" << r.degrees[i];
  }

  return os;
}

template<typename R>
int lexi_compare(const monomial<R> &lhs, const monomial<R> &rhs) {
  for (int i = 0; i < lhs.degrees.size(); i++) {
    if (lhs.degrees[i] > rhs.degrees[i]) {
      return 1;
    }
    if (lhs.degrees[i] < rhs.degrees[i]) {
      return -1;
    }
  }
  // All degree terms are the same
  return 0;
}

/**
 * Ordered according to lexicographic order
 * @tparam R ring
 */
template<typename R>
struct multipoly {
  list<monomial<R>> monomials;
};

template <typename R>
std::ostream& operator<<(std::ostream& os, const multipoly<R>& poly) {
  if (poly.monomials.empty()) {
    os << "0";
  } else {
    auto ite = poly.monomials.begin();
    os << *ite;
    ++ite;
    for (; ite != poly.monomials.end(); ++ite) {
      os << " + " << *ite;
    }
  }

  return os;
}

template <typename R>
multipoly<R> &add_monomial(multipoly<R> &poly, const monomial<R> &monomial) {
  auto ite = poly.monomials.begin();
  for (; ite != poly.monomials.end(); ++ite) {
    int cmp = lexi_compare(monomial, *ite);
    if (cmp < 0) {
      poly.monomials.insert(ite, monomial);
      return poly;
    }
    if (cmp == 0) {
      ite->coefficient = ite->coefficient + monomial.coefficient;
      if (ite->coefficient == 0) {
        poly.monomials.erase(ite);
      }
      return poly;
    }
  }

  // Fail to insert before the very back
  poly.monomials.insert(ite, monomial);
  return poly;
}

template <typename R>
multipoly<R> operator+(const multipoly<R> &lhs, const multipoly<R> &rhs) {
  auto result = lhs;
  for (auto &rhs_mono : rhs.monomials) {
    add_monomial(result, rhs_mono);
  }
  return result;
}

template <typename R>
multipoly<R> operator-(const multipoly<R> &lhs, const multipoly<R> &rhs) {
  auto result = lhs;
  for (auto &rhs_mono : rhs.monomials) {
    auto negate_rhs = rhs_mono;
    negate_rhs.coefficient = -negate_rhs.coefficient;
    add_monomial(result, negate_rhs);
  }
  return result;
}

template <typename R>
multipoly<R> operator*(const multipoly<R> &lhs, const multipoly<R> &rhs) {
  auto result = lhs;
  // TODO Implement multiplication
  return result;
}

template <typename R>
bool is_zero(multipoly<R> &poly) {
  return poly.monomials.empty();
}

int main() {
  auto m1 = monomial{2, vector<size_t>{3, 4}};
  auto m2 = monomial{2, vector<size_t>{4, 2}};
  auto m3 = monomial{3, vector<size_t>{3, 2}};
  cout << "m1, m2, m3: " << m1 << ", " << m2 << ", " << m3 << endl;
  cout << "compare m1, m2: " << lexi_compare(m1, m2) << endl;
  cout << "compare m1, m3: " << lexi_compare(m1, m3) << endl;
  cout << "compare m1, m1: " << lexi_compare(m1, m1) << endl;

  cout << endl;
  auto poly = multipoly<int>();
  cout << "poly - step 0: " << poly << endl;
  add_monomial(poly, m1);
  cout << "poly - step 1: " << poly << endl;
  add_monomial(poly, m2);
  cout << "poly - step 2: " << poly << endl;
  add_monomial(poly, m3);
  cout << "poly - step 3: " << poly << endl;
  add_monomial(poly, m1);
  cout << "poly - step 4: " << poly << endl;

  cout << "poly + poly: " << poly + poly << endl;
  auto sub = poly - poly;
  cout << "poly - poly: " << sub << endl;
  cout << "poly - poly is zero: " << is_zero(sub) << endl;

  return 0;
}