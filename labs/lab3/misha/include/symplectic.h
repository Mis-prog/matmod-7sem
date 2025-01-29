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