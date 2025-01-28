//#include "chart.h"
//
//class Symplectic : public Chart {
//    const double Dtheta = 0.1931833275037836;
//
//public:
//    Symplectic() : Chart("symplectic") {}
//
//    void Calculate() {
//        for (int i = CurrentStep; i < NumberOfSteps; ++i) {
//            for (size_t j = 0; j < Offsets.size(); ++j) {
//                Offsets[j] += Speeds[j] * Dtheta * Tau;
//            }
//
//            Accelerations = CalcCommonAccelerations();
//
//            for (size_t j = 0; j < Speeds.size(); ++j) {
//                Speeds[j] += Accelerations[j] / 2.0 * Tau;
//                Offsets[j] += Speeds[j] * (1 - 2 * Dtheta) * Tau;
//            }
//
//            Accelerations = CalcCommonAccelerations();
//
//            for (size_t j = 0; j < Speeds.size(); ++j) {
//                Speeds[j] += Accelerations[j] / 2.0 * Tau;
//                Offsets[j] += Speeds[j] * Dtheta * Tau;
//            }
//
//            FiniteHamiltonian = CalcHamiltonian();
//            ++CurrentStep;
//
//            if (i % UpdateStep == 0){
//                this->SaveResult();
//            }
//        }
//    }
//};


// Symplectic.hpp
#pragma once
#include "chart.h"
#include <thread>

class Symplectic : public Chart {
private:
    static constexpr long  double Dtheta = 0.1931833275037836;
    void calculate();

public:
    Symplectic();
    void start() override;
};

// Symplectic.cpp

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
