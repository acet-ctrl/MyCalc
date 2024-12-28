#pragma once

#include "matrix.h"

namespace math {
void plu_decomposition(Matrix&);
void determinant(Matrix&);
void inverse(Matrix&);
void qr_decomposition(Matrix&);
void svd_decomposition(Matrix&);
void jordan_canonical(Matrix&);
}  // namespace math