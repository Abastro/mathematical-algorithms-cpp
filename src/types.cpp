/*
 * Until next time, let's boilerplate later
 */
#include <cmath>
#include <iostream>
using namespace std;

template<typename R>
struct rational {
  R num;
  R denom;
};

template<typename R>
rational<R> operator+(rational<R> &l, rational<R> &r) {
  return rational(l.num * r.denom + r.num * l.denom, l.denom * r.denom);
}

template<typename R>
rational<R> operator-(rational<R> &l, rational<R> &r) {
  return rational(l.num * r.denom - r.num * l.denom, l.denom * r.denom);
}

template<typename R>
rational<R> operator*(rational<R> &l, rational<R> &r) {
  return rational(l.num * r.num, l.denom * r.denom);
}

template<typename R>
rational<R> operator/(rational<R> &l, rational<R> &r) {
  return rational(l.num * r.denom, r.denom * l.num);
}

// TODO Take care of normalization
template<typename R>
bool operator==(rational<R> &l, rational<R> &r) {
  return l.num == r.num && l.denom == r.denom;
}

template<typename R>
bool operator!=(rational<R> &l, rational<R> &r) {
  return l.num != r.num || l.denom != r.denom;
}

/*
template<typename R>
ostream operator<< (ostream &os, rational<R> &r) {

}
*/

template<unsigned int>
struct prime_field {
  size_t number;
};

template<unsigned int p>
prime_field<p> operator+(const prime_field<p> &l, const prime_field<p> &r) {
  return prime_field(l.number + r.number % p);
}

template<unsigned int p>
prime_field<p> operator-(const prime_field<p> &l, const prime_field<p> &r) {
  return prime_field(l.number - r.number % p);
}

template<unsigned int p>
prime_field<p> operator*(const prime_field<p> &l, const prime_field<p> &r) {
  return prime_field(l.number * r.number % p);
}

// TODO /, ==, !=

int main() {
  cout << "Hello World!" << endl;
}
