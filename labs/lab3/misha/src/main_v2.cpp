// Matmod3.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//
#include <omp.h>
#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

double F(double q_last, double q, double q_next, double m, double alpha, double beta) {

    return 1.0 / (m) * ((q_next - 2 * q + q_last) + alpha * (q_next - q) * (q_next - q) +
                        beta * (q_next - q) * (q_next - q) * (q_next - q)
                        - alpha * (q - q_last) * (q - q_last) - beta * (q - q_last) * (q - q_last) * (q - q_last));
}

double Vf(double alpha, double beta, double r) {
    return r * r / 2 + alpha * r * r * r / 3 + beta * r * r * r * r / 4;
}

double H(vector<double> q, vector<double> v, double m, double alpha, double beta) {
    double summ = 0;
    for (int i = 0; i < q.size() - 1; i++) {
        summ += (m * v[i] * v[i] / 2 + Vf(alpha, beta, q[i + 1] - q[i]));
    }
    summ += (m * v[q.size() - 1] * v[q.size() - 1] / 2 + Vf(alpha, beta, (q[0] - q[q.size() - 1])));
    return summ;
}

vector<vector<double>>
SpeedVerle(vector<double> q, vector<double> v, double alpha, double beta, double tau, int N, double m) {
    // double m = 2.0;
    double l = 1;
    ofstream f("../labs/lab3/misha/result/speed.txt");
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
            for (int i = 0; i < v.size(); i++) f << v[i] << " ";
            f << endl;
        }
    }

    f.close();
    return {q, v};
}

double ksi = 0.1931833275037836;

vector<vector<double>>
SimplexVerle(vector<double> q, vector<double> v, double alpha, double beta, double tau, int N, double m) {
    vector<double> a(q.size(), 0);
    int s = q.size() - 1;
    ofstream f("../labs/lab3/misha/result/simplex.txt");
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

//        if (t > N / 1000 and t < N / 800)
        for (int i = s / 2; i >= 0; i--) {
            v[i] = -v[s - i];
            q[i] = -q[s - i];
//                v[s - i] = -v[i];
//                q[s - i] = -q[i];
        }

        if (t % 1000 == 0) {
            //  cout << t + 1 << endl;
            //    cout << v[s / 2] << " " << v[s / 2 + 1] << endl;
            for (int i = 0; i < v.size(); i++) f << v[i] << " ";
            f << endl;
        }
    }
    f.close();
    return {q, v};
}

int main() {
    // main1();
    int N = 1000;
    vector<double> q(N, 0);
    vector<double> v(N, 0);
    double m = 1;
    double l = 1;
    q[N / 2 - 1] = 0.5;
    q[N / 2] = -0.5;
    //  v[N / 2 - 1] = -0.7;
    //  v[N / 2] = 0.7;
    double a = 0.0;
    double b = 100.0;
    std::cout << "Введите a и b:\n";
    std::cin >> a >> b;
    double H0 = H(q, v, m, a, b);
    cout << H0 << endl;
    double start = clock();
    auto res = SimplexVerle(q, v, a, b, 0.01, 1000000, m);
    double finish = clock();
    cout << (finish - start) / CLOCKS_PER_SEC << endl;
    //  for (int i = 0; i < res[0].size(); i++) cout << res[0][i] << " " << res2[0][i] << endl;
    cout << (H(res[0], res[1], m, a, b) - H0) / H0 << " ";// << H(res2[0], res2[1], 1, a, b) - H0 << endl;
}


// 1
// 2 127(+-) 129.1
// 3