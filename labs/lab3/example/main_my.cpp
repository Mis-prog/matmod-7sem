#include "symplectic.h"
#include "verle.h"

int main() {
    Verle verleMethod;
    double alpha, beta;
    std::cout << "Введите альфа и бетта:\n";
    std::cin >> alpha, beta;
    verleMethod.SetAlphaBetta(alpha, beta);
    verleMethod.Calculate();

    return 0;
}