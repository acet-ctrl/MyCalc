#include "matrix.h"

#include <stdexcept>

#include "rational.h"

namespace math {
Matrix::Matrix(int rows, int cols) : rows_(rows), cols_(cols) {
  if (rows_ <= 0 || cols_ <= 0)
    throw std::invalid_argument("Invalid matrix size");
  data_ = new Rational[rows * cols];
}

Matrix::~Matrix() { delete[] data_; }

Rational& Matrix::operator()(int i, int j) {
  if (i < 0 || i >= rows_ || j < 0 || j >= cols_)
    throw std::out_of_range("Index out of range");
  return data_[i * cols_ + j];
}

Rational Matrix::operator()(int i, int j) const {
  if (i < 0 || i >= rows_ || j < 0 || j >= cols_)
    throw std::out_of_range("Index out of range");
  return data_[i * cols_ + j];
}

std::ostream& operator<<(std::ostream& os, const Matrix& matrix) {
  int c = 0;
  for (int i = 0; i < matrix.rows_ * matrix.cols_; ++i) {
    os << matrix.data_[i] << ' ';
    c = (c + 1) % matrix.cols_;
    os << (c ? " " : "\n");
  }
  return os;
}

std::istream& operator>>(std::istream& is, Matrix& matrix) {
  for (int i = 0; i < matrix.rows_ * matrix.cols_; ++i) is >> matrix.data_[i];
  return is;
}

SwapOperation::SwapOperation(int i, int j)
    : ElementaryOperation(-1), i_(i), j_(j) {
  if (i_ == j_) throw std::invalid_argument("Indices must be different");
}

SwapOperation& SwapOperation::rowOperate(Matrix& matrix) {
  if (i_ < 0 || i_ >= matrix.rows_ || j_ < 0 || j_ >= matrix.rows_)
    throw std::out_of_range("Index out of range");
  int cols = matrix.cols_;
  Rational* i = matrix.data_ + i_ * cols;
  Rational* j = matrix.data_ + j_ * cols;
  while (cols--) std::swap(*i++, *j++);
  return *this;
}

SwapOperation& SwapOperation::colOperate(Matrix& matrix) {
  if (i_ < 0 || i_ >= matrix.cols_ || j_ < 0 || j_ >= matrix.cols_)
    throw std::out_of_range("Index out of range");
  int rows = matrix.rows_;
  Rational* i = matrix.data_ + i_;
  Rational* j = matrix.data_ + j_;
  while (rows--) {
    std::swap(*i, *j);
    i += matrix.cols_;
    j += matrix.cols_;
  }
  return *this;
}

SwapOperation& SwapOperation::inverse() { return *this; }

MultiplyOperation::MultiplyOperation(int i, const Rational& scalar)
    : ElementaryOperation(scalar), i_(i), scalar_(scalar) {}

MultiplyOperation& MultiplyOperation::rowOperate(Matrix& matrix) {
  if (i_ < 0 || i_ >= matrix.rows_)
    throw std::out_of_range("Index out of range");
  int cols = matrix.cols_;
  Rational* row = matrix.data_ + i_ * cols;
  while (cols--) *row++ *= value_;
  return *this;
}

MultiplyOperation& MultiplyOperation::colOperate(Matrix& matrix) {
  if (i_ < 0 || i_ >= matrix.cols_)
    throw std::out_of_range("Index out of range");
  int rows = matrix.rows_;
  Rational* col = matrix.data_ + i_;
  while (rows--) {
    *col *= value_;
    col += matrix.cols_;
  }
  return *this;
}

MultiplyOperation& MultiplyOperation::inverse() {
  scalar_ = Rational(1) / scalar_;
  value_ = Rational(1) / value_;
  return *this;
}

AddOperation::AddOperation(int i, int j, const Rational& scalar)
    : ElementaryOperation(1), i_(i), j_(j), scalar_(scalar) {
  if (i_ == j_) throw std::invalid_argument("Indices must be different");
}

AddOperation& AddOperation::rowOperate(Matrix& matrix) {
  if (i_ < 0 || i_ >= matrix.rows_ || j_ < 0 || j_ >= matrix.rows_)
    throw std::out_of_range("Index out of range");
  int cols = matrix.cols_;
  Rational* i = matrix.data_ + i_ * cols;
  Rational* j = matrix.data_ + j_ * cols;
  while (cols--) *i++ += value_ * *j++;
  return *this;
}

AddOperation& AddOperation::colOperate(Matrix& matrix) {
  if (i_ < 0 || i_ >= matrix.cols_ || j_ < 0 || j_ >= matrix.cols_)
    throw std::out_of_range("Index out of range");
  int rows = matrix.rows_;
  Rational* i = matrix.data_ + i_;
  Rational* j = matrix.data_ + j_;
  while (rows--) {
    *i += value_ * *j;
    i += matrix.cols_;
    j += matrix.cols_;
  }
  return *this;
}

AddOperation& AddOperation::inverse() {
  scalar_ = -scalar_;
  return *this;
}
}  // namespace math