#include "../include/symplectic.h"

Symplectic::Symplectic() : Chart("symplectic") {}

void Symplectic::start() {
    isStarted = true;
    calculate();
}

void Symplectic::calculate() {
    for (int i = currentStep; i < numberOfSteps; ++i) {
        for (int j = 0; j < numberOfParticles; ++j) {
            offsets[j] = offsets[j] + speeds[j] * Dtheta * tau;
        }

        accelerations = calcCommonAccelerations();

        for (int j = 0; j < numberOfParticles; ++j) {
            speeds[j] = speeds[j] + accelerations[j] / 2.0 * tau;
            offsets[j] = offsets[j] + speeds[j] * (1 - 2 * Dtheta) * tau;
        }

        accelerations = calcCommonAccelerations();

        for (int j = 0; j < numberOfParticles; ++j) {
            speeds[j] = speeds[j] + accelerations[j] / 2.0 * tau;
            offsets[j] = offsets[j] + speeds[j] * Dtheta * tau;
        }

        finiteHamiltonian = calcHamiltonian();
        currentStep++;



//        if (!isStarted) break;

        if (currentStep % updateStep == 0) {
            callDataChanged();
        }
    }
    std::cout << "Конечный гамильтон\n";
    std::cout << finiteHamiltonian << std::endl;
    isStarted = false;
    callDataChanged();
}
