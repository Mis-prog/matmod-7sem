//#include "symplectic.h"
//#include "verle.h"
//
//int main() {
//    Verle verleMethod;
//    double alpha, my_beta;
//    std::cout << "Введите альфа и бетта:\n";
//    std::cin >> alpha, my_beta;
//    verleMethod.SetAlphaBetta(alpha, my_beta);
//    verleMethod.Calculate();
//
//    return 0;
//}


#include "verle.h"
#include "symplectic.h"
#include <iostream>

int main(){
    Verle methodVerle;
    methodVerle.start();
    return 0;
}
