#pragma once

#include <memory>

#include "rational.h"

namespace math {
class Matrix {
 public:
  Matrix(int, int);
  Matrix(int);
  Matrix() = delete;

  static Matrix Identity(int);

  Rational& operator()(int, int);
  Rational operator()(int, int) const;

  Rational dotRow(int, int) const;
  Rational dotCol(int, int) const;

  friend std::ostream& operator<<(std::ostream&, const Matrix&);
  friend std::istream& operator>>(std::istream&, Matrix&);

  friend class SwapOperation;
  friend class MultiplyOperation;
  friend class AddOperation;

  int rows() const { return rows_; }
  int cols() const { return cols_; }

 private:
  int rows_;
  int cols_;
  std::shared_ptr<Rational[]> data_;
};

class ElementaryOperation {
 public:
  virtual ~ElementaryOperation() = default;

  ElementaryOperation& operate(Matrix& matrix);
  virtual ElementaryOperation& inverse() = 0;

 protected:
  ElementaryOperation(bool isRowOp) : isRowOp_(isRowOp) {};

  virtual ElementaryOperation& rowOperate(Matrix& matrix) = 0;
  virtual ElementaryOperation& colOperate(Matrix& matrix) = 0;

  bool isRowOp_;
};

class SwapOperation : public ElementaryOperation {
 public:
  SwapOperation(bool, int, int);
  SwapOperation(int, int);

  SwapOperation& inverse() override;

 private:
  SwapOperation& rowOperate(Matrix& matrix) override;
  SwapOperation& colOperate(Matrix& matrix) override;

  int i_;
  int j_;
};

class MultiplyOperation : public ElementaryOperation {
 public:
  MultiplyOperation(bool, int, const Rational&);
  MultiplyOperation(int, const Rational&);

  MultiplyOperation& inverse() override;

 private:
  int i_;
  Rational scalar_;

  MultiplyOperation& rowOperate(Matrix& matrix) override;
  MultiplyOperation& colOperate(Matrix& matrix) override;
};

class AddOperation : public ElementaryOperation {
 public:
  AddOperation(bool, int, int, const Rational&);
  AddOperation(int, int, const Rational&);

  AddOperation& inverse() override;

 private:
  int i_;
  int j_;
  Rational scalar_;

  AddOperation& rowOperate(Matrix& matrix) override;
  AddOperation& colOperate(Matrix& matrix) override;
};
}  // namespace math