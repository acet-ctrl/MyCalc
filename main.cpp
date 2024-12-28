#include <iostream>

#include "math/mycalc.h"

int main() {
  while (true) {
    std::cout << "Welcome to MyCalc!\n"
              << "    1: PLU decomposition\n"
              << "    2: Determinant\n"
              << "    3: Inverse\n"
              << "    4: QR decomposition\n"
              << "    5: SVD decomposition\n"
              << "    6: Jordan canonical form\n"
              << "Other: Exit\n"
              << "Enter your choice: ";
    char choice;
    std::cin >> choice;
    if (choice < '1' || choice > '6') {
      std::cout << "Goodbye!\n";
      break;
    }
    int rows, cols;
    std::cout << "Enter the number of rows and columns: ";
    std::cin >> rows >> cols;
    math::Matrix A(rows, cols);
    std::cout << "Enter the matrix:\n";
    std::cin >> A;
    switch (choice) {
      case '1':
        math::plu_decomposition(A);
        break;
      case '2':
        math::determinant(A);
        break;
      case '3':
        math::inverse(A);
        break;
      case '4':
        math::qr_decomposition(A);
        break;
      case '5':
        math::svd_decomposition(A);
        break;
      default:
        math::jordan_canonical(A);
        break;
    }
  }
}
