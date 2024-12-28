#pragma once

#include <sys/stat.h>

#include <memory>

#include "rational.h"

namespace math {
class Matrix {
 public:
  Matrix(int, int);
  Matrix(int);
  Matrix() = delete;
  Matrix(const Matrix&) = delete;
  Matrix(Matrix&&) = default;
  ~Matrix();

  static Matrix Identity(int);

  Rational& operator()(int, int);
  Rational operator()(int, int) const;

  friend std::ostream& operator<<(std::ostream&, const Matrix&);
  friend std::istream& operator>>(std::istream&, Matrix&);

  friend class ElementaryOperation;
  friend class SwapOperation;
  friend class MultiplyOperation;
  friend class AddOperation;

  int rows() const { return rows_; }
  int cols() const { return cols_; }

 private:
  int rows_;
  int cols_;
  Rational* data_;
};

class ElementaryOperation {
 public:
  virtual ~ElementaryOperation() = default;

  virtual ElementaryOperation& rowOperate(Matrix& matrix) = 0;
  virtual ElementaryOperation& colOperate(Matrix& matrix) = 0;

  virtual ElementaryOperation& inverse() = 0;

 protected:
  ElementaryOperation(const Rational& value) : value_(value) {};

  Rational value_;
};

class SwapOperation : public ElementaryOperation {
 public:
  SwapOperation(int, int);

  SwapOperation& rowOperate(Matrix& matrix) override;
  SwapOperation& colOperate(Matrix& matrix) override;

  SwapOperation& inverse() override;

 private:
  int i_;
  int j_;
};

class MultiplyOperation : public ElementaryOperation {
 public:
  MultiplyOperation(int, const Rational&);

  MultiplyOperation& rowOperate(Matrix& matrix) override;
  MultiplyOperation& colOperate(Matrix& matrix) override;

  MultiplyOperation& inverse() override;

 private:
  int i_;
  Rational scalar_;
};

class AddOperation : public ElementaryOperation {
 public:
  AddOperation(int, int, const Rational&);

  AddOperation& rowOperate(Matrix& matrix) override;
  AddOperation& colOperate(Matrix& matrix) override;

  AddOperation& inverse() override;

 private:
  int i_;
  int j_;
  Rational scalar_;
};
}  // namespace math