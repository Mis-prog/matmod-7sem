#include "include/speed.h"
#include "include/symplectic.h"

int main() {
    Symplectic speedMethods;
    speedMethods.setAlphaBeta(0, 300);
    speedMethods.setupdateStep(10);
    speedMethods.setnumberOfSteps(1e6);
    speedMethods.start();
    return 0;
}