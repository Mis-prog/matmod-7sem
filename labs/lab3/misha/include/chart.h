#pragma once

#include <vector>
#include <functional>
#include <cmath>
#include <string>
#include <fstream>
#include <iostream>

class Chart {
protected:
    std::string header;
    bool isStarted = false;
    long double alpha = 0.6;
    long double beta = 0;
    long double tau = 0.005;
    long double mass = 0.5;
    int numberOfSteps = 100000;
    int currentStep = 0;
    int updateStep = 500;
    int numberOfParticles = 1000;
    long double initialDeviation = 0.5;
    long double initialHamiltonian = 0;
    long double finiteHamiltonian = 0;
    std::ofstream out;

    std::vector<long double> offsets;
    std::vector<long double> speeds;
    std::vector<long double> accelerations;

    // Callback for data changes
    std::function<void()> dataChangedCallback;

    std::vector<long double> calcCommonAccelerations();

    long double calcHamiltonian();

    std::vector<long double> calcGradV();

    long double calcGradV(long double qM1,long double q,long double qP1);

    std::vector<long double> calcV();

    long double calcV(long double q,long double qP1);

    void callDataChanged();

public:
    Chart(const std::string &header);

    virtual void start() = 0;

    virtual void stop();

    virtual void clear();

    void setAlphaBeta(long double alpha, long double beta);
    void setnumberOfSteps(int numberOfSteps);
    void setupdateStep(int updateStep);

    ~Chart();
};