#include "mycalc.h"

#include <iostream>
#include <ostream>

#include "matrix.h"
#include "rational.h"

namespace math {
void plu_decomposition(Matrix& A) {
  int rows = A.rows(), cols = A.cols();
  Matrix P = Matrix::Identity(rows), L = Matrix::Identity(rows), U = A;
  int pivCol = 0;
  for (int curRow = 0; curRow < rows && pivCol < cols; ++curRow, ++pivCol) {
    int pivRow;
    for (pivRow = curRow; pivCol < cols; ++pivCol) {
      for (int tryRow = curRow + 1; tryRow < rows; ++tryRow)
        if (U(tryRow, pivCol).abs() > U(pivRow, pivCol).abs()) pivRow = tryRow;
      if (U(pivRow, pivCol) != 0) break;
    }
    if (pivCol == cols) break;
    if (pivRow > curRow)
      SwapOperation(pivRow, curRow)
          .operate(U)
          .operate(L)
          .inverse()
          .operate(L)
          .operate(P);
    Rational pivEle = U(curRow, pivCol);
    for (int modRow = curRow + 1; modRow < rows; ++modRow) {
      AddOperation(modRow, curRow, -U(modRow, pivCol) / pivEle)
          .operate(U)
          .inverse()
          .operate(L);
    }
  }
  std::cout << "A = PLU, where\n"
            << "P =\n"
            << P << "L =\n"
            << L << "U =\n"
            << U;
}

void determinant(Matrix& A) {
  int rows = A.rows();
  if (rows != A.cols()) {
    std::cerr << "Matrix is not square\n";
    return;
  }
  Rational det = 1;
  for (int curRow = 0; curRow < rows; ++curRow) {
    int tryRow;
    for (tryRow = curRow; tryRow < rows; ++tryRow)
      if (A(tryRow, curRow) != 0) break;
    if (tryRow == rows) {
      det = 0;
      break;
    }
    if (tryRow > curRow) {
      det = -det;
      SwapOperation(tryRow, curRow).operate(A);
    }
    Rational pivEle = A(curRow, curRow);
    det *= pivEle;
    MultiplyOperation(curRow, Rational(1) / pivEle).operate(A);
    for (int modRow = curRow + 1; modRow < rows; ++modRow)
      AddOperation(modRow, curRow, -A(modRow, curRow)).operate(A);
  }
  std::cout << "det(A) = " << det << '\n';
}

void inverse(Matrix& A) {
  int rows = A.rows(), cols = A.cols();
  if (rows != cols) {
    std::cerr << "Matrix is not square\n";
    return;
  }
  Matrix B = Matrix::Identity(rows);
  for (int curRow = 0; curRow < rows; ++curRow) {
    int tryRow;
    for (tryRow = curRow; tryRow < rows; ++tryRow)
      if (A(tryRow, curRow) != 0) break;
    if (tryRow == rows) {
      std::cerr << "Matrix is not invertible\n";
      return;
    }
    if (tryRow > curRow) SwapOperation(tryRow, curRow).operate(A).operate(B);
    MultiplyOperation(curRow, Rational(1) / A(curRow, curRow))
        .operate(A)
        .operate(B);
    for (int modRow = 0; modRow < rows; ++modRow)
      if (modRow != curRow)
        AddOperation(modRow, curRow, -A(modRow, curRow)).operate(A).operate(B);
  }
  std::cout << "A^(-1) =\n" << B;
}

void qr_decomposition(Matrix& A) {
  int cols = A.cols();
  Matrix Q = A, R = Matrix::Identity(cols);
  Rational* selfDot = new Rational[cols];
  for (int curCol = 0; curCol < cols; ++curCol) {
    for (int refCol = 0; refCol < curCol; ++refCol) {
      if (selfDot[refCol] == 0) continue;
      AddOperation(false, curCol, refCol,
                   -Q.dotCol(curCol, refCol) / selfDot[refCol])
          .operate(Q)
          .inverse()
          .operate(R);
    }
    selfDot[curCol] = Q.dotCol(curCol, curCol);
  }
  std::cout << "A = QR, where\n"
            << "Q =\n"
            << Q << "R =\n"
            << R;
}

void svd_decomposition(Matrix& A) {}

void jordan_canonical(Matrix& A) {}
}  // namespace math