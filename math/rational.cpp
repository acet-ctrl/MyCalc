#include "rational.h"

#include <numeric>
#include <stdexcept>

namespace math {
Rational::Rational(int numerator, int denominator)
    : numerator_(numerator), denominator_(denominator) {
  if (!denominator_) throw std::invalid_argument("Denominator cannot be zero");
  Reduce();
}

Rational::Rational(int numerator) : numerator_(numerator), denominator_(1) {}

Rational Rational::operator-() const {
  return Rational(-numerator_, denominator_);
}

Rational Rational::operator+(const Rational& rational) const {
  return Rational(
      numerator_ * rational.denominator_ + rational.numerator_ * denominator_,
      denominator_ * rational.denominator_);
}

Rational Rational::operator-(const Rational& rational) const {
  return Rational(
      numerator_ * rational.denominator_ - rational.numerator_ * denominator_,
      denominator_ * rational.denominator_);
}

Rational Rational::operator*(const Rational& rational) const {
  return Rational(numerator_ * rational.numerator_,
                  denominator_ * rational.denominator_);
}

Rational Rational::operator/(const Rational& rational) const {
  if (!rational.numerator_) throw std::invalid_argument("Division by zero");
  return Rational(numerator_ * rational.denominator_,
                  denominator_ * rational.numerator_);
}

Rational& Rational::operator+=(const Rational& rational) {
  numerator_ =
      numerator_ * rational.denominator_ + rational.numerator_ * denominator_;
  denominator_ *= rational.denominator_;
  Reduce();
  return *this;
}

Rational& Rational::operator-=(const Rational& rational) {
  numerator_ =
      numerator_ * rational.denominator_ - rational.numerator_ * denominator_;
  denominator_ *= rational.denominator_;
  Reduce();
  return *this;
}

Rational& Rational::operator*=(const Rational& rational) {
  numerator_ *= rational.numerator_;
  denominator_ *= rational.denominator_;
  Reduce();
  return *this;
}

Rational& Rational::operator/=(const Rational& rational) {
  if (!rational.numerator_) throw std::invalid_argument("Division by zero");
  numerator_ *= rational.denominator_;
  denominator_ *= rational.numerator_;
  Reduce();
  return *this;
}

bool Rational::operator==(const Rational& rational) const {
  return numerator_ == rational.numerator_ &&
         denominator_ == rational.denominator_;
}

bool Rational::operator!=(const Rational& rational) const {
  return !(*this == rational);
}

bool Rational::operator<(const Rational& rational) const {
  return numerator_ * rational.denominator_ <
         rational.numerator_ * denominator_;
}

bool Rational::operator<=(const Rational& rational) const {
  return numerator_ * rational.denominator_ <=
         rational.numerator_ * denominator_;
}

bool Rational::operator>(const Rational& rational) const {
  return numerator_ * rational.denominator_ >
         rational.numerator_ * denominator_;
}

bool Rational::operator>=(const Rational& rational) const {
  return numerator_ * rational.denominator_ >=
         rational.numerator_ * denominator_;
}

void Rational::Reduce() {
  if (!numerator_) {
    denominator_ = 1;
    return;
  }
  int gcd = std::gcd(numerator_, denominator_);
  numerator_ /= gcd;
  denominator_ /= gcd;
  if (denominator_ < 0) {
    numerator_ = -numerator_;
    denominator_ = -denominator_;
  }
}

std::ostream& operator<<(std::ostream& os, const Rational& rational) {
  os << rational.numerator_;
  if (rational.denominator_ != 1) os << '/' << rational.denominator_;
  return os;
}

std::istream& operator>>(std::istream& is, Rational& rational) {
  is >> rational.numerator_;
  if (is.peek() == '/') {
    is.ignore();
    is >> rational.denominator_;
    rational.Reduce();
  } else
    rational.denominator_ = 1;
  return is;
}
}  // namespace math