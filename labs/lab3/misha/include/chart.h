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
    double alpha = 0.6;
    double beta = 0;
    double tau = 0.01;
    double mass = 0.5;
    int numberOfSteps = 100000;
    int currentStep = 0;
    int updateStep = 500;
    int numberOfParticles = 500;
    double initialDeviation = 0.5;
    double initialHamiltonian = 0;
    double finiteHamiltonian = 0;
    std::ofstream out;

    std::vector<double> offsets;
    std::vector<double> speeds;
    std::vector<double> accelerations;

    // Callback for data changes
    std::function<void()> dataChangedCallback;

    std::vector<double> calcCommonAccelerations();

    double calcHamiltonian();

    std::vector<double> calcGradV();

    double calcGradV(double qM1,double q,double qP1);

    std::vector<double> calcV();

    double calcV(double q,double qP1);

    void callDataChanged();

public:
    Chart(const std::string &header);

    virtual void start() = 0;

    virtual void stop();

    virtual void clear();

    void setAlphaBeta(double alpha, double beta);
    void setnumberOfSteps(int numberOfSteps);
    void setupdateStep(int updateStep);

    ~Chart();
};