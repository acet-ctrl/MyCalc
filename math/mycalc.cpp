#include "mycalc.h"

#include "matrix.h"

namespace math {
void plu_decomposition(Matrix& U) {
  int rows = U.rows();
  Matrix P = Matrix::Identity(rows), L(rows);
  for (int i = 0; i < rows; ++i) {
    int pivot = i;
    for (int j = i + 1; j < rows; ++j)
      if (U(j, i).abs() > U(pivot, i).abs()) pivot = j;
    if (pivot != i)
      SwapOperation(pivot, i).rowOperate(U).inverse().rowOperate(P);
    for (int j = i + 1; j < rows; ++j) {
      Rational factor = U(j, i) / U(i, i);
      AddOperation(j, i, -factor).rowOperate(U);
      U(j, i) = factor;
    }
  }
}
}  // namespace math