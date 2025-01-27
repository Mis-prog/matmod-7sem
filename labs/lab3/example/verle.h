#include "chart.h"

class Verle : public Chart {
public:
    Verle() : Chart("verle") {}

    void Calculate() {
        for (int i = CurrentStep; i < NumberOfSteps; ++i) {
            for (size_t j = 0; j < Offsets.size(); ++j) {
                Offsets[j] += Speeds[j] * Tau + Accelerations[j] / 2.0 * Tau * Tau;
                Speeds[j] += Accelerations[j] / 2.0 * Tau;
            }

            Accelerations = CalcCommonAccelerations();

            for (size_t j = 0; j < Speeds.size(); ++j) {
                Speeds[j] += Accelerations[j] / 2.0 * Tau;
            }

            FiniteHamiltonian = CalcHamiltonian();
            ++CurrentStep;

            if (i % UpdateStep == 0){
                this->SaveResult();
            }
        }
    }
};