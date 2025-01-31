#include "../include/speed.h"

Speed::Speed() : Chart("speed") {}

void Speed::start() {
    isStarted = true;
    initialHamiltonian = calcHamiltonian();
    std::cout << "Начальный гамильтон\n";
    std::cout << initialHamiltonian << std::endl;
    this->calculate();
}

void Speed::calculate() {
    for (int i = currentStep; i < numberOfSteps; ++i) {
        for (int j = 0; j < numberOfParticles; ++j) {
            offsets[j] = offsets[j] + speeds[j] * tau +
                         accelerations[j] / 2.0 * tau * tau;
            speeds[j] = speeds[j] + accelerations[j] / 2.0 * tau;
        }

        accelerations = calcCommonAccelerations();

        for (int j = 0; j < numberOfParticles; ++j) {
            speeds[j] = speeds[j] + accelerations[j] / 2.0 * tau;
        }

        finiteHamiltonian = calcHamiltonian();
        currentStep++;


        if (currentStep % updateStep == 0) {
            callDataChanged();
        }
    }
    std::cout << "Конечный гамильтон\n";
    std::cout << finiteHamiltonian << std::endl;

    isStarted = false;
    callDataChanged();
}