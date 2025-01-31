#pragma once

#include "chart.h"
#include <thread>

class Speed : public Chart {
private:
    void calculate();

public:
    Speed();

    void start() override;
};