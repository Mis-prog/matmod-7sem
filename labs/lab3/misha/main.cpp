#include "include/speed.h"
#include "include/symplectic.h"

int main() {
    Speed speedMethods;
    speedMethods.setAlphaBeta(0, 500);
    speedMethods.setupdateStep(500);
    speedMethods.setnumberOfSteps(1e5);
    speedMethods.start();
    return 0;
}

// 2 cолитона 500 сипмлекс
// 1 солитона 600 сиплекс 655, 550
// 3 солитона 700 симплекс 695