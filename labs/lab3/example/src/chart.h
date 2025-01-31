//#pragma once
//#include <vector>
//#include <string>
//#include <functional>
//#include <thread>
//#include <cmath>
//#include <fstream>
//#include <iomanip>
//#include <iostream>
//
//class Chart {
//private:
//    void Clear() {
//            CurrentStep = 0;
//
//            Offsets.resize(NumberOfParticles, 0.0);
//            Offsets[NumberOfParticles / 2 - 1] = - InitialDeviation;
//            Offsets[NumberOfParticles / 2] =  InitialDeviation;
//
//            Speeds.resize(NumberOfParticles, 0.0);
//            Accelerations = CalcCommonAccelerations();
//        }
//protected:
//    std::string Header;
//    bool IsStarted = false;
//    double Alpha = 0.0;
//    double Beta = 0.7;
//    double Tau = 0.01;
//    double Mass = 0.5;
//    int NumberOfSteps = 1'000'000;
//    int CurrentStep = 0;
//    int UpdateStep = 10;
//    int NumberOfParticles = 500;
//    double InitialDeviation = 0.5;
//    double InitalHamiltonian = 0.0;
//    double FiniteHamiltonian = 0.0;
//
//    std::vector<double> Offsets;
//    std::vector<double> Speeds;
//    std::vector<double> Accelerations;
//
//    std::ofstream fout;
//
//public:
//    Chart(const std::string& header) : Header(header) {
//        Clear();
//        fout.open("../labs/lab3/example/result/" + header + ".txt");
//    }
//
//    void SetAlphaBetta(double Alpha,double Beta){
//        this -> Alpha = Alpha;
//        this -> Beta = Beta;
//        InitalHamiltonian = CalcHamiltonian();
//        std :: cout << "Начальный гамильтолиан для "<< Header << ": "<< std::fixed << std::setprecision(15) << InitalHamiltonian  << std::endl;
//    }
//
//    ~Chart() {
//        fout.close();
//        std :: cout << "Конечный гамильтолиан для " << Header << ": " << std::fixed << std::setprecision(15) << FiniteHamiltonian << std::endl;
//    }
//
//protected:
//
//    void SaveResult(){
//        for (auto value : Offsets ){
//            fout << value << " ";
//        }
//        fout << std::endl;
//    }
//
//    std::vector<double> CalcCommonAccelerations() {
//        auto gradV = CalcGradV();
//        for (auto& acc : gradV) {
//            acc *= -1.0 / Mass;
//        }
//        return gradV;
//    }
//
//    double CalcHamiltonian() {
//        double result = 0.0;
//        auto V = CalcV();
//        for (size_t i = 0; i < V.size(); ++i) {
//            result += Speeds[i] * Speeds[i] * Mass / 2 + V[i];
//        }
//        return result;
//    }
//
//    std::vector<double> CalcGradV() {
//        std::vector<double> result(NumberOfParticles);
//
//        result[0] = CalcGradV(Offsets[NumberOfParticles - 1], Offsets[0], Offsets[1]);
//        result[NumberOfParticles - 1] = CalcGradV(Offsets[NumberOfParticles - 2],
//                                                  Offsets[NumberOfParticles - 1],
//                                                  Offsets[0]);
//
//        for (int i = 1; i < NumberOfParticles - 1; ++i) {
//            result[i] = CalcGradV(Offsets[i - 1], Offsets[i], Offsets[i + 1]);
//        }
//
//        return result;
//    }
//
//    double CalcGradV(double qM1, double q, double qP1) {
//        return -(qP1 - q) + (q - qM1) +
//               Alpha * (-std::pow(qP1 - q, 2) + std::pow(q - qM1, 2)) +
//               Beta * (-std::pow(qP1 - q, 3) + std::pow(q - qM1, 3));
//    }
//
//    std::vector<double> CalcV() {
//        std::vector<double> result(Offsets.size());
//        for (size_t i = 0; i < result.size() - 1; ++i) {
//            result[i] = CalcV(Offsets[i], Offsets[i + 1]);
//        }
//        result[result.size() - 1] = CalcV(Offsets[result.size() - 1], Offsets[0]);
//        return result;
//    }
//
//    double CalcV(double q, double qP1) {
//        double r = qP1 - q;
//        return r * r / 2.0 + Alpha * r * r * r / 3.0 + Beta * r * r * r * r / 4.0;
//    }
//};

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

Chart::Chart(const std::string &header) : header(header) {
    out.open("../labs/lab3/example/result/" + header + ".txt");
    clear();
}

void Chart::stop() {
    isStarted = false;
    callDataChanged();
}

void Chart::clear() {
    isStarted = false;
    currentStep = 0;

    offsets.resize(numberOfParticles, 0.0);
    offsets[numberOfParticles / 2 - 1] = -initialDeviation;
    offsets[numberOfParticles / 2] = initialDeviation;

    speeds.resize(numberOfParticles, 0.0);
    accelerations = calcCommonAccelerations();

    initialHamiltonian = calcHamiltonian();
    std::cout << "Начальный гамильтон\n";
    std::cout << initialHamiltonian << std::endl;
    finiteHamiltonian = initialHamiltonian;

    callDataChanged();
}

std::vector<double> Chart::calcCommonAccelerations() {
    std::vector<double> result = calcGradV();
    for (double &val: result) {
        val *= -1.0 / mass;
    }
    return result;
}

double Chart::calcHamiltonian() {
    double result = 0;
    std::vector<double> V = calcV();
    for (size_t i = 0; i < V.size(); ++i) {
        result += speeds[i] * speeds[i] * mass / 2 + V[i];
    }
    return result;
}

std::vector<double> Chart::calcGradV() {
    std::vector<double> result(numberOfParticles);

    result[0] = calcGradV(offsets[numberOfParticles - 1], offsets[0], offsets[1]);
    result[numberOfParticles - 1] = calcGradV(offsets[numberOfParticles - 2],
                                              offsets[numberOfParticles - 1],
                                              offsets[0]);

    for (int i = 1; i < numberOfParticles - 1; ++i) {
        result[i] = calcGradV(offsets[i - 1], offsets[i], offsets[i + 1]);
    }

    return result;
}

double Chart::calcGradV(double qM1,double q,double qP1) {
    return -(qP1 - q) + (q - qM1) +
           alpha * (-std::pow(qP1 - q, 2) + std::pow(q - qM1, 2)) +
           beta * (-std::pow(qP1 - q, 3) + std::pow(q - qM1, 3));
}

std::vector<double> Chart::calcV() {
    std::vector<double> result(offsets.size());
    for (size_t i = 0; i < result.size() - 1; ++i) {
        result[i] = calcV(offsets[i], offsets[i + 1]);
    }
    result[result.size() - 1] = calcV(offsets[result.size() - 1], offsets[0]);
    return result;
}

double Chart::calcV(double q,double qP1) {
    double r = qP1 - q;
    return r * r / 2.0 + alpha * r * r * r / 3.0 + beta * r * r * r * r / 4.0;
}

void Chart::callDataChanged() {
//    if (dataChangedCallback) {
//        dataChangedCallback();
//    }
    for (auto value: speeds) {
        out << value << " ";
    }
    out << std::endl;
}