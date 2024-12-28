#pragma once

#include <iostream>

namespace math {
class Rational {
 public:
  Rational(int, int);
  Rational(int);
  Rational();

  Rational operator-() const;

  Rational operator+(const Rational&) const;
  Rational operator-(const Rational&) const;
  Rational operator*(const Rational&) const;
  Rational operator/(const Rational&) const;

  Rational abs() const;
  Rational pow(uint32_t) const;

  Rational& operator+=(const Rational&);
  Rational& operator-=(const Rational&);
  Rational& operator*=(const Rational&);
  Rational& operator/=(const Rational&);

  bool operator==(const Rational&) const;
  bool operator!=(const Rational&) const;
  bool operator<(const Rational&) const;
  bool operator<=(const Rational&) const;
  bool operator>(const Rational&) const;
  bool operator>=(const Rational&) const;

  friend std::ostream& operator<<(std::ostream&, const Rational&);
  friend std::istream& operator>>(std::istream&, Rational&);

 private:
  int numerator_;
  int denominator_;

  void Reduce();
};
}  // namespace math