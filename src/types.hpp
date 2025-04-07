#pragma once
#include <cmath>
#include <iostream>
using namespace std;

/**
 * A generic function for computing the GCD
 * @tparam R : this template variable should be a type for which all the necessary
 * operations for Euclidean domains have been implemented
 * @param a
 * @param b
 * @return
 */
template <typename R>
R gcd(R a, R b) {
  if (b == 0)
    return a;
  R r = a % b;
  while (r != 0) {
    a = b;
    b = r;
    r = a % b;
  }
  return b;
}

/**
 * Rational on R which is always normalized.
 * @tparam R : a type for Euclidean domain
 */
template<typename R>
struct ratio {
  R num;
  R denom;
};

template<typename R>
ratio<R> &rational_of(R n) {
  return ratio<R>(n, 1);
}

template<typename R>
ratio<R> &normalize(ratio<R> &n) {
  const R gcd = gcd(n.num, n.denom);
  n.num = n.num / gcd;
  n.denom = n.denom / gcd;
  return n;
}

template<typename R>
ratio<R> &operator+(const ratio<R> &l, const ratio<R> &r) {
  return normalize(rational(l.num * r.denom + r.num * l.denom, l.denom * r.denom));
}

template<typename R>
ratio<R> &operator-(const ratio<R> &l, const ratio<R> &r) {
  return normalize(rational(l.num * r.denom - r.num * l.denom, l.denom * r.denom));
}

template<typename R>
ratio<R> &operator*(const ratio<R> &l, const ratio<R> &r) {
  return normalize(rational(l.num * r.num, l.denom * r.denom));
}

template<typename R>
ratio<R> &operator/(const ratio<R> &l, const ratio<R> &r) {
  return normalize(rational(l.num * r.denom, r.denom * l.num));
}

template<typename R>
bool operator==(const ratio<R> &l, const ratio<R> &r) {
  return l.num == r.num && l.denom == r.denom;
}

template<typename R>
bool operator!=(const ratio<R> &l, const ratio<R> &r) {
  return l.num != r.num || l.denom != r.denom;
}

template <typename R>
std::ostream& operator<<(std::ostream& os, const ratio<R>& r) {
  os << r.a << "/" << r.b;

  return os;
}

/**
 * Prime field where number is kept to be in range [0, p).
 */
template<unsigned int>
struct prime_field {
  size_t number;
};

template<unsigned int p>
prime_field<p> &operator+(const prime_field<p> &l, const prime_field<p> &r) {
  return prime_field((l.number + r.number) % p);
}

template<unsigned int p>
prime_field<p> &operator-(const prime_field<p> &l, const prime_field<p> &r) {
  return prime_field((l.number - r.number) % p);
}

template<unsigned int p>
prime_field<p> &operator*(const prime_field<p> &l, const prime_field<p> &r) {
  return prime_field((l.number * r.number) % p);
}

// TODO prime_field /

template<unsigned int p>
bool operator==(const prime_field<p> &l, const prime_field<p> &r) {
  return l.number == r.number;
}

template<unsigned int p>
bool operator!=(const prime_field<p> &l, const prime_field<p> &r) {
  return l.number != r.number;
}