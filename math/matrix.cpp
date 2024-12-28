#include "matrix.h"

#include <memory>
#include <stdexcept>
#include <utility>

#include "rational.h"

namespace math {
Matrix::Matrix(int rows, int cols) : rows_(rows), cols_(cols) {
  if (rows_ <= 0 || cols_ <= 0)
    throw std::invalid_argument("Invalid matrix size");
  data_ = std::shared_ptr<Rational[]>(new Rational[rows * cols]);
}

Matrix::Matrix(int n) : Matrix(n, n) {}

Matrix Matrix::Identity(int n) {
  Matrix matrix(n);
  for (int i = 0; i < n; ++i) matrix(i, i) = 1;
  return matrix;
}

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

Rational Matrix::dotRow(int i, int j) const {
  if (i < 0 || i >= rows_ || j < 0 || j >= rows_)
    throw std::out_of_range("Index out of range");
  Rational result;
  Rational* row1 = data_.get() + i * cols_;
  Rational* row2 = data_.get() + j * cols_;
  for (int k = 0; k < cols_; ++k) result += *row1++ * *row2++;
  return result;
}

Rational Matrix::dotCol(int i, int j) const {
  if (i < 0 || i >= cols_ || j < 0 || j >= cols_)
    throw std::out_of_range("Index out of range");
  Rational result;
  Rational* col1 = data_.get() + i;
  Rational* col2 = data_.get() + j;
  for (int k = 0; k < rows_; ++k) {
    result += *col1 * *col2;
    col1 += cols_;
    col2 += cols_;
  }
  return result;
}

std::ostream& operator<<(std::ostream& os, const Matrix& matrix) {
  int c = 0;
  for (int i = 0; i < matrix.rows_ * matrix.cols_; ++i) {
    os << matrix.data_[i];
    c = (c + 1) % matrix.cols_;
    os << (c ? " " : "\n");
  }
  return os;
}

std::istream& operator>>(std::istream& is, Matrix& matrix) {
  for (int i = 0; i < matrix.rows_ * matrix.cols_; ++i) is >> matrix.data_[i];
  return is;
}

ElementaryOperation& ElementaryOperation::operate(Matrix& matrix) {
  return isRowOp_ ? rowOperate(matrix) : colOperate(matrix);
}

SwapOperation::SwapOperation(bool isRowOp, int i, int j)
    : ElementaryOperation(isRowOp), i_(i), j_(j) {
  if (i_ == j_) throw std::invalid_argument("Indices must be different");
}

SwapOperation::SwapOperation(int i, int j) : SwapOperation(true, i, j) {}

SwapOperation& SwapOperation::rowOperate(Matrix& matrix) {
  if (i_ < 0 || i_ >= matrix.rows_ || j_ < 0 || j_ >= matrix.rows_)
    throw std::out_of_range("Index out of range");
  int cols = matrix.cols_;
  Rational* i = matrix.data_.get() + i_ * cols;
  Rational* j = matrix.data_.get() + j_ * cols;
  while (cols--) std::swap(*i++, *j++);
  return *this;
}

SwapOperation& SwapOperation::colOperate(Matrix& matrix) {
  if (i_ < 0 || i_ >= matrix.cols_ || j_ < 0 || j_ >= matrix.cols_)
    throw std::out_of_range("Index out of range");
  int rows = matrix.rows_;
  Rational* i = matrix.data_.get() + i_;
  Rational* j = matrix.data_.get() + j_;
  while (rows--) {
    std::swap(*i, *j);
    i += matrix.cols_;
    j += matrix.cols_;
  }
  return *this;
}

SwapOperation& SwapOperation::inverse() {
  isRowOp_ = !isRowOp_;
  return *this;
}

MultiplyOperation::MultiplyOperation(bool isRowOp, int i,
                                     const Rational& scalar)
    : ElementaryOperation(isRowOp), i_(i), scalar_(scalar) {
  if (scalar_ == 0) throw std::invalid_argument("Scalar cannot be zero");
}

MultiplyOperation::MultiplyOperation(int i, const Rational& scalar)
    : MultiplyOperation(true, i, scalar) {}

MultiplyOperation& MultiplyOperation::rowOperate(Matrix& matrix) {
  if (i_ < 0 || i_ >= matrix.rows_)
    throw std::out_of_range("Index out of range");
  int cols = matrix.cols_;
  Rational* row = matrix.data_.get() + i_ * cols;
  while (cols--) *row++ *= scalar_;
  return *this;
}

MultiplyOperation& MultiplyOperation::colOperate(Matrix& matrix) {
  if (i_ < 0 || i_ >= matrix.cols_)
    throw std::out_of_range("Index out of range");
  int rows = matrix.rows_;
  Rational* col = matrix.data_.get() + i_;
  while (rows--) {
    *col *= scalar_;
    col += matrix.cols_;
  }
  return *this;
}

MultiplyOperation& MultiplyOperation::inverse() {
  isRowOp_ = !isRowOp_;
  scalar_ = Rational(1) / scalar_;
  return *this;
}

AddOperation::AddOperation(bool isRowOp, int i, int j, const Rational& scalar)
    : ElementaryOperation(isRowOp), i_(i), j_(j), scalar_(scalar) {
  if (i_ == j_) throw std::invalid_argument("Indices must be different");
}

AddOperation::AddOperation(int i, int j, const Rational& scalar)
    : AddOperation(true, i, j, scalar) {}

AddOperation& AddOperation::rowOperate(Matrix& matrix) {
  if (scalar_ == 0) return *this;
  if (i_ < 0 || i_ >= matrix.rows_ || j_ < 0 || j_ >= matrix.rows_)
    throw std::out_of_range("Index out of range");
  int cols = matrix.cols_;
  Rational* i = matrix.data_.get() + i_ * cols;
  Rational* j = matrix.data_.get() + j_ * cols;
  while (cols--) *i++ += scalar_ * *j++;
  return *this;
}

AddOperation& AddOperation::colOperate(Matrix& matrix) {
  if (scalar_ == 0) return *this;
  if (i_ < 0 || i_ >= matrix.cols_ || j_ < 0 || j_ >= matrix.cols_)
    throw std::out_of_range("Index out of range");
  int rows = matrix.rows_;
  Rational* i = matrix.data_.get() + i_;
  Rational* j = matrix.data_.get() + j_;
  while (rows--) {
    *i += scalar_ * *j;
    i += matrix.cols_;
    j += matrix.cols_;
  }
  return *this;
}

AddOperation& AddOperation::inverse() {
  isRowOp_ = !isRowOp_;
  scalar_ = -scalar_;
  std::swap(i_, j_);
  return *this;
}
}  // namespace math