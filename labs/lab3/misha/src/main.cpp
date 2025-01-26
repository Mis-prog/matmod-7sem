#include <fstream>
#include <iostream>
#include <omp.h>
#include <time.h>
#include <vector>
#include <iomanip>
#include <sstream>

using namespace std;

// 
double F(double q_last, double q, double q_next, double m, double alpha, double beta) {
    return q_next - 2 * q + q_last +
           alpha * ((q_next - q) * (q_next - q) - (q - q_last) * (q - q_last)) +
           beta * ((q_next - q) * (q_next - q) * (q_next - q) -
                   (q - q_last) * (q - q_last) * (q - q_last));
}

// Потенциал
double V(double alpha, double beta, double r) {
    return r * r / 2 + alpha * r * r * r / 3 + beta * r * r * r * r / 4;
}

// Полная энергия
double H(vector<double> q, vector<double> v, double m, double alpha,
         double beta) {
    double summ = 0;
    for (int i = 0; i < q.size() - 1; i++) {
        summ += (m * v[i] * v[i] / 2 + V(alpha, beta, q[i + 1] - q[i]));
    }
    summ += (m * v[q.size() - 1] * v[q.size() - 1] / 2 +
             V(alpha, beta, (q[0] - q[q.size() - 1])));
    return summ;
}

// Скорострой алгоритм Верле
vector<vector<double>> SpeedVerle(vector<double> q, vector<double> v,
                                  double alpha, double beta, double tau, int N,
                                  double m) {
    double l = 1;
    string title = "../labs/lab3/misha/result/speed_" + std::to_string(beta);+".txt";
    ofstream f(title);
    vector<double> a(q.size(), 0);
    int s = q.size() - 1;
    for (int t = 0; t < N; t++) {

        for (int i = 0; i < s + 1; i++) {
            q[i] = q[i] + v[i] * tau + 0.5 * a[i] * tau * tau;
            v[i] = v[i] + 0.5 * a[i] * tau;
        }
        a[0] = F(l * q[s], l * q[0], l * q[1], m, alpha, beta) / l;
        v[0] = v[0] + 0.5 * a[0] * tau;
        for (int i = 1; i < s; i++) {
            a[i] = F(q[i - 1], q[i], q[i + 1], m, alpha, beta);
            v[i] = v[i] + 0.5 * a[i] * tau;
        }
        a[s] = F(q[s - 1], q[s], q[0], m, alpha, beta) / l;
        v[s] = v[s] + 0.5 * a[s] * tau;
        if (t % 1000 == 0) {
            for (int i = 0; i < v.size(); i++)
                f << v[i] << " ";
            f << endl;
        }
    }

    f.close();
    return {q, v};
}

// Симплектический алгоритм Верле
double ksi = 0.1931833275037836;
vector<vector<double>> SimplexVerle(vector<double> q, vector<double> v,
                                    double alpha, double beta, double tau,
                                    int N, double m) {
    vector<double> a(q.size(), 0);
    int s = q.size() - 1;
    ostringstream stream;
    stream << fixed << setprecision(2) << beta;
    string title = "../labs/lab3/misha/result/symplectic_" + stream.str() + ".txt";
    ofstream f(title);
    for (int t = 0; t < N; t++) {

        for (int i = 0; i < s + 1; i++) {
            q[i] = q[i] + v[i] * tau * ksi;
        }
        a[0] = F(q[s], q[0], q[1], m, alpha, beta);
        for (int i = 1; i < s; i++) {
            a[i] = F(q[i - 1], q[i], q[i + 1], m, alpha, beta);
        }
        a[s] = F(q[s - 1], q[s], q[0], m, alpha, beta);
        for (int i = 0; i < s + 1; i++) {
            v[i] = v[i] + 0.5 * a[i] * tau;
            q[i] = q[i] + v[i] * tau * (1 - 2 * ksi);
        }
        a[0] = F(q[s], q[0], q[1], m, alpha, beta);
        for (int i = 1; i < s; i++) {
            a[i] = F(q[i - 1], q[i], q[i + 1], m, alpha, beta);
        }
        a[s] = F(q[s - 1], q[s], q[0], m, alpha, beta);

        for (int i = 0; i < s + 1; i++) {
            v[i] = v[i] + 0.5 * a[i] * tau;
            q[i] = q[i] + v[i] * tau * ksi;
        }
        if (t % 1000 == 0) {
            for (int i = 0; i < v.size(); i++)
                f << v[i] << " ";
            f << endl;
        }
    }

    f.close();
    return {q, v};
}

int main() {
    int N = 1000;
    vector<double> q(N, 0);
    vector<double> v(N, 0);
    double m = 1;
    q[N / 2 - 1] = 1;
    q[N / 2] = -1;
    double a = 0.7;
    double b = 0;
    double H0 = H(q, v, m, a, b);
    cout << "H0: " << H0 << endl;
    double start = clock();
    // auto res = SpeedVerle(q, v, a, b, 0.02, 2000000, m);
    auto res = SimplexVerle(q, v, a, b, 0.01,2000000,m);
    double finish = clock();
    cout << "Time: " << (finish - start) / CLOCKS_PER_SEC << endl;
    cout << (H(res[0], res[1], m, a, b) - H0) / H0 << " ";
}