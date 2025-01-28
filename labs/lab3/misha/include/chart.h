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

    virtual ~Chart() = default;

    virtual void start() = 0;

    virtual void stop();

    virtual void clear();

    // Getters and setters
    bool getIsStarted() const { return isStarted; }

    void setIsStarted(bool value) { isStarted = value; }

    double getAlpha() const { return alpha; }

    void setAlpha(double value) { alpha = value; }

    double getBeta() const { return beta; }

    void setBeta(double value) { beta = value; }

    double getTau() const { return tau; }

    void setTau(double value) { tau = value; }

    double getMass() const { return mass; }

    void setMass(double value) { mass = value; }

    int getNumberOfSteps() const { return numberOfSteps; }

    void setNumberOfSteps(int value) { numberOfSteps = value; }

    int getCurrentStep() const { return currentStep; }

    void setCurrentStep(int value) { currentStep = value; }

    int getUpdateStep() const { return updateStep; }

    void setUpdateStep(int value) { updateStep = value; }

    int getNumberOfParticles() const { return numberOfParticles; }

    void setNumberOfParticles(int value) {
        numberOfParticles = value;
        clear();
    }

    double getInitialDeviation() const { return initialDeviation; }

    void setInitialDeviation(double value) { initialDeviation = value; }

    void setDataChangedCallback(std::function<void()> callback) {
        dataChangedCallback = callback;
    }
};