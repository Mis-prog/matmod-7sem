#include "chart.h"

class Symplectic : public Chart {
    const double Dtheta = 0.1931833275037836;

public:
    Symplectic() : Chart("symplectic") {}

    void Calculate() {
        for (int i = CurrentStep; i < NumberOfSteps; ++i) {
            for (size_t j = 0; j < Offsets.size(); ++j) {
                Offsets[j] += Speeds[j] * Dtheta * Tau;
            }

            Accelerations = CalcCommonAccelerations();

            for (size_t j = 0; j < Speeds.size(); ++j) {
                Speeds[j] += Accelerations[j] / 2.0 * Tau;
                Offsets[j] += Speeds[j] * (1 - 2 * Dtheta) * Tau;
            }

            Accelerations = CalcCommonAccelerations();

            for (size_t j = 0; j < Speeds.size(); ++j) {
                Speeds[j] += Accelerations[j] / 2.0 * Tau;
                Offsets[j] += Speeds[j] * Dtheta * Tau;
            }

            FiniteHamiltonian = CalcHamiltonian();
            ++CurrentStep;

            if (i % UpdateStep == 0){
                this->SaveResult();
            }
        }
    }
};