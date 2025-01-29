#include "include/speed.h"
#include "include/symplectic.h"

int main() {
    Symplectic speedMethods;
    speedMethods.setAlphaBeta(0, 600);
    speedMethods.setupdateStep(50);
    speedMethods.setnumberOfSteps(5e4);
    speedMethods.start();
    return 0;
}

// 2 cолитона 500 сипмлекс
// 1 солитона 600 сиплекс 655, 550
// 3 солитона 700 симплекс 695