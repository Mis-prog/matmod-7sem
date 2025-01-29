//#include "chart.h"
//
//class Verle : public Chart {
//public:
//    Verle() : Chart("verle") {}
//
//    void Calculate() {
//        for (int i = CurrentStep; i < NumberOfSteps; ++i) {
//            for (size_t j = 0; j < Offsets.size(); ++j) {
//                Offsets[j] += Speeds[j] * Tau + Accelerations[j] / 2.0 * Tau * Tau;
//                Speeds[j] += Accelerations[j] / 2.0 * Tau;
//            }
//
//            Accelerations = CalcCommonAccelerations();
//
//            for (size_t j = 0; j < Speeds.size(); ++j) {
//                Speeds[j] += Accelerations[j] / 2.0 * Tau;
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

// Verle.hpp
#pragma once

#include "chart.h"
#include <thread>

class Verle : public Chart {
private:
    void calculate();

public:
    Verle();

    void start() override;
};

// Verle.cpp

Verle::Verle() : Chart("verle") {}

void Verle::start() {
    isStarted = true;
    this->calculate();
//    std::thread calcThread(&Verle::calculate, this);
//    calcThread.detach();
}

void Verle::calculate() {
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