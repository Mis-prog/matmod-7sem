#include <iostream>
#include <vector>
#include <cmath>
#include <fstream> 
#include <string>

double a;
double b;

std::vector<double> system(const std::vector<double>& state) {
    std::vector<double> derivatives(3);
    derivatives[0] = state[1];
    derivatives[1] = a*state[1] + (1- state[2])* state[0];
    derivatives[2] = b*state[2] + state[0] * state[0];
    return derivatives;
}

std::vector<double> rungeKutta4(const std::vector<double>& state, double dt) {
    std::vector<double> k1 = system(state);
    std::vector<double> stateK2(3), stateK3(3), stateK4(3);

    for (int i = 0; i < 3; ++i) stateK2[i] = state[i] + 0.5 * dt * k1[i];
    std::vector<double> k2 = system(stateK2);

    for (int i = 0; i < 3; ++i) stateK3[i] = state[i] + 0.5 * dt * k2[i];
    std::vector<double> k3 = system(stateK3);

    for (int i = 0; i < 3; ++i) stateK4[i] = state[i] + dt * k3[i];
    std::vector<double> k4 = system(stateK4);

    std::vector<double> nextState(3);
    for (int i = 0; i < 3; ++i) {
        nextState[i] = state[i] + (dt / 6.0) * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
    }

    return nextState;
}

int main() {
    cout << "Введите значения a и b\n"; 
    cin >> a >> b;
    std::vector<double> state = { 0, 0, 0 };
    double dt = 0.01;
    int numSteps = 4000;
    std::ofstream file("../labs/lab4/misha/result/trajectory"+ std::to_string(a) +  "_" +  std::to_string(b) ".csv");
    std::ofstream << "x y z";
    for (int step = 0; step < numSteps; ++step) {
        state = rungeKutta4(state, dt);
       
        file << state[0] << " " << state[1] << " " << state[2] << "\n";
    }
    file.close();

    return 0;
}
