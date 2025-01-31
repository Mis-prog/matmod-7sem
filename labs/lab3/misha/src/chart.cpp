#include "../include/chart.h"

Chart::Chart(const std::string &header) : header(header) {
    out.open("../labs/lab3/misha/result/" + header + ".txt");
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

    finiteHamiltonian = initialHamiltonian;

    callDataChanged();
}

std::vector<long double> Chart::calcCommonAccelerations() {
    std::vector<long double> result = calcGradV();
    for (long double &val: result) {
        val *= -1.0 / mass;
    }
    return result;
}

long double Chart::calcHamiltonian() {
    long double result = 0;
    std::vector<long double> V = calcV();
    for (size_t i = 0; i < V.size(); ++i) {
        result += speeds[i] * speeds[i] * mass / 2 + V[i];
    }
    return result;
}

std::vector<long double> Chart::calcGradV() {
    std::vector<long double> result(numberOfParticles);

    result[0] = calcGradV(offsets[numberOfParticles - 1], offsets[0], offsets[1]);
    result[numberOfParticles - 1] = calcGradV(offsets[numberOfParticles - 2],
                                              offsets[numberOfParticles - 1],
                                              offsets[0]);

    for (int i = 1; i < numberOfParticles - 1; ++i) {
        result[i] = calcGradV(offsets[i - 1], offsets[i], offsets[i + 1]);
    }

    return result;
}

long double Chart::calcGradV(long double qM1,long  double q,long  double qP1) {
    return -(qP1 - q) + (q - qM1) +
           alpha * (-std::pow(qP1 - q, 2) + std::pow(q - qM1, 2)) +
           beta * (-std::pow(qP1 - q, 3) + std::pow(q - qM1, 3));
}

std::vector<long double> Chart::calcV() {
    std::vector<long double> result(offsets.size());
    for (size_t i = 0; i < result.size() - 1; ++i) {
        result[i] = calcV(offsets[i], offsets[i + 1]);
    }
    result[result.size() - 1] = calcV(offsets[result.size() - 1], offsets[0]);
    return result;
}

long double Chart::calcV(long double q, long double qP1) {
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


void Chart::setAlphaBeta(long double alpha,long  double beta) {
    this->alpha = alpha;
    this->beta = beta;
}

Chart::~Chart() {
    out.close();
}

void Chart::setnumberOfSteps(int numberOfSteps) {
    this->numberOfSteps = numberOfSteps;
}

void Chart::setupdateStep(int updateStep) {
    this->updateStep = updateStep;
}