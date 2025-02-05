#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

using namespace std;

// Проверка выполнения закона сохранения энергии
const int n = 500; // Число частиц
const int N_time = 1e6; // Число слоев по времени
const double alpha = 1, my_beta = 1.;
const double tau = 0.01, m = 1.0;
const double xi = 0.1931833275037836;
double q0 = 0.06;

//1 солитон alpha
//const int n = 1000; // число частиц
//const int N_time = 1e5; // число слоев по времени
//const double alpha = 0.7, my_beta = 0.;
//const double tau = 0.01, m = 1.0;
//const double xi = 0.1931833275037836;
//double q0 = 0.5;

// 1 солитон my_beta
//const int n = 1000; // Число частиц
//const int N_time = 1e5; // Число слоев по времени
//const double alpha = 0., my_beta = 50;
//const double tau = 0.01, m = 1.0;
//const double xi = 0.1931833275037836;
//double q0 = 0.9;
const double v0 = 1.;

// 2 солитона my_beta
//const int n = 1000; // Число частиц
//const int N_time = 3e6; // Число слоев по времени
//const double alpha = 0., my_beta = 106.;
//const double tau = 0.01, m = 1.0;
//const double xi = 0.1931833275037836;
//const double q0 = 0.9;

// 3 солитона my_beta
//const int n = 1000; // Число частиц
//const int N_time = 1e5; // Число слоев по времени
//const double alpha = 0., my_beta = 200;
//const double tau = 0.01, m = 1.0;
//const double xi = 0.1931833275037836;
//double q0 = 0.9;

std::vector<double> GradV(std::vector<double>& q)
{
    std::vector<double> GradV_vector(q.size());
    for (int i = 1; i < n - 1; i++)
    {
        GradV_vector[i] = (q[i + 1] - 2 * q[i] + q[i - 1]) +
                          alpha * (pow(q[i + 1] - q[i], 2) - pow(q[i] - q[i - 1], 2)) +
                          my_beta * (pow(q[i + 1] - q[i], 3) - pow(q[i] - q[i - 1], 3));
    }
    return GradV_vector;
}

void Verle(std::vector<double>& q, std::vector<double>& v, std::vector<double>& a)
{
    std::ofstream f_out;
    f_out.open("data//Verle.dat");
    clock_t t_solv_end = clock();
    for (int i = 0; i < N_time; i++)
    {
        clock_t t = clock();
        if (i % 100 == 0)
        {
            for (int j = 0; j < n; j++)
            {
                f_out << v[j] << " ";
            }
            f_out << "\n";
        }

        for (int j = 1; j < n - 1; j++)
        {
            q[j] = q[j] + v[j] * tau + 0.5 * a[j] * tau * tau;

            v[j] = v[j] + 0.5 * a[j] * tau;
        }

        q[0] = q[n - 2];
        v[0] = v[n - 2];

        q[n - 1] = q[1];
        v[n - 1] = v[1];

        std::vector<double> GradV_vector = GradV(q);

        for (int j = 1; j < n - 1; j++)
        {
            a[j] = 1. / m * GradV_vector[j];

            v[j] = v[j] + 0.5 * a[j] * tau;
        }
        t = clock() - t;
#ifndef NDEBUG
        std::cout << "Time for one iteration" << t << " clicks (" << ((double)t) / CLOCKS_PER_SEC << " seconds).\n";
        std::cout << "End calculation for t = " << (i * 100) / (N_time) << "%\n";
#endif //  NDEBUG
        if (i % 100 == 0) std::cout << "End calculation Verle for t = " << (i * 100) / (N_time) << "%\n";
    }
    t_solv_end = clock() - t_solv_end;
    std::cout << "Time Solve Verle: " << t_solv_end << " clicks (" << ((double)t_solv_end) / CLOCKS_PER_SEC << " seconds).\n";
    for (int j = 0; j < n; j++)
    {
        f_out << v[j] << " ";
    }
    f_out << "\n";
    f_out.close();
}

void simplexVerle(std::vector<double>& q, std::vector<double>& v, std::vector<double>& a)
{
    std::ofstream f_out;
    f_out.open("data//simplexVerle.dat");
    clock_t t_solv_end = clock();

    for (int i = 0; i < N_time; i++)
    {
        for (int j = 0; j < n; j++)
        {
            f_out << v[j] << " ";
        }
        f_out << "\n";
        clock_t t = clock();
        for (int j = 0; j < n; j++)
        {
            q[j] = q[j] + v[j] * xi * tau;
        }

        std::vector<double> GradV_vector = GradV(q);

        for (int j = 1; j < n - 1; j++)
        {
            a[j] = 1. / m * GradV_vector[j];

            v[j] = v[j] + 0.5 * a[j] * tau;
        }

        for (int j = 0; j < n; j++)
        {
            q[j] = q[j] + v[j] * (1. - 2. * xi) * tau;
        }

        GradV_vector = GradV(q);
        for (int j = 1; j < n - 1; j++)
        {
            a[j] = 1. / m * GradV_vector[j];

            v[j] = v[j] + 0.5 * a[j] * tau;

            q[j] = q[j] + v[j] * xi * tau;
        }
        int s = q.size() / 2;
        for (int i = s + 1; i < q.size(); i++) {
            v[i] = -v[q.size() - 1 - i];
            a[i] = -a[q.size() - 1 - i];
            q[i] = -q[q.size() - 1 - i];
        }


        t = clock() - t;
#ifndef NDEBUG
        //std::cout << "Time for one iteration" << t << " clicks (" << ((double)t) / CLOCKS_PER_SEC << " seconds).\n";
        //std::cout << "End calculation for t = " << (i * 100) / (N_time) << "%\n";
#endif //  NDEBUG
        //if (i % 100 == 0)  std::cout << "End calculation simplexVerle for t = " << (i * 100) / (N_time) << "%\n";
    }
    t_solv_end = clock() - t_solv_end;
    std::cout << "Time Solve simplexVerle: " << t_solv_end << " clicks (" << ((double)t_solv_end) / CLOCKS_PER_SEC << " seconds).\n";


    for (int j = 0; j < n; j++)
    {
        f_out << v[j] << " ";
    }
    f_out << "\n";
    f_out.close();
}

double V(double r)
{
    return (r * r / 2. + alpha * r * r * r / 3. + my_beta * r * r * r * r / 4.);
}

double Gamiltonian(std::vector<double>& v, std::vector<double>& q)
{
    double P = 0.0;
    // Левая точка выколота
    for (int i = 1; i < n - 1; i++)
    {
        P += m * v[i] * v[i] / 2. + V(q[i + 1] - q[i]);
    }
    return P;
}

double F(double q_last, double q, double q_next, double m, double alpha, double beta) {

    return 1.0 / (m) * ((q_next - 2 * q + q_last) + alpha * (q_next - q) * (q_next - q) + beta * (q_next - q) * (q_next - q) * (q_next - q)
                        - alpha * (q - q_last) * (q - q_last) - beta * (q - q_last) * (q - q_last) * (q - q_last));
}
vector<vector<double>> SimplexVerle2(vector<double> q, vector<double>v, double alpha, double beta, double tau, int N, double m) {
    vector<double> a(q.size(), 0);
    int s = q.size() / 2 - 1;
    ofstream f("data//speed.txt");
    for (int t = 0; t < N; t++) {

        for (int i = 0; i < s + 1; i++) {
            q[i] = q[i] + v[i] * tau * xi;
        }
        a[0] = F(-q[0], q[0], q[1], m, alpha, beta);
        for (int i = 1; i < s; i++) {
            a[i] = F(q[i - 1], q[i], q[i + 1], m, alpha, beta);
        }
        a[s] = F(q[s - 1], q[s], -q[s], m, alpha, beta);
        for (int i = 0; i < s + 1; i++) {
            v[i] = v[i] + 0.5 * a[i] * tau;
            q[i] = q[i] + v[i] * tau * (1 - 2 * xi);
        }
        a[0] = F(-q[0], q[0], q[1], m, alpha, beta);
        for (int i = 1; i < s; i++) {
            a[i] = F(q[i - 1], q[i], q[i + 1], m, alpha, beta);
        }
        a[s] = F(q[s - 1], q[s], -q[s], m, alpha, beta);

        for (int i = 0; i < s + 1; i++) {
            v[i] = v[i] + 0.5 * a[i] * tau;
            q[i] = q[i] + v[i] * tau * xi;
        }
        for (int i = s + 1; i < q.size(); i++) {
            v[i] = -v[q.size() - 1 - i];
            a[i] = -a[q.size() - 1 - i];
            q[i] = -q[q.size() - 1 - i];
        }
        if (t % 300 == 0 && t > 000000) {
            //  cout << t + 1 << endl;
            //    cout << v[s / 2] << " " << v[s / 2 + 1] << endl;
            for (int i = 0; i < v.size(); i++) f << v[i] << " ";
            f << endl;
        }
    }
    for (int i = 0; i < v.size(); i++) f << v[i] << " ";
    // f << endl;
    f.close();
    return { q,v };
}


int main()
{
    std::vector<double> q(n, 0), v(n, 0), a(n, 0);
    std::vector<double> GradV_vector;
    for (int i = 0; i < n; i++)
    {
        q[i] = v[i] = a[i] = 0.0;
    }
    v[n / 2] = 1.;
    v[n / 2 + 1] = -v0;
    //q[n / 2 - 1] = q0;
    //q[n / 2] = -q0;
    GradV_vector = GradV(q);
    for (int j = 1; j < n - 1; j++)
    {
        a[j] = 1. / m * GradV_vector[j];
    }
    double start = Gamiltonian(v, q);
    std::cout << "first H = " << start << "\n";
    //simplexVerle(q, v, a);
    simplexVerle(q, v, a);
    double end = Gamiltonian(v, q);
    std::cout << "after simplexVerle H = " << end << "\n";
    std::cout << "simplexVerle уrror: " << abs(start - end) << "\n";

    return 0;
}

// Начальные условия
//q[n / 2 + 1] = q0;
//q[n / 2] = -q0;
//GradV_vector = GradV(q);
//for (int j = 1; j < n - 1; j++)
//{
//	a[j] = 1. / m * GradV_vector[j];
//}

//double start = Gamiltonian(v, q);
//std::cout << "first H = " << start << "\n";

//Verle(q, v, a);
//double end = Gamiltonian(v, q);
//std::cout << "after Verle H = " << end << "\n";
//std::cout << "Verle уrror: " << abs(start - end) << "\n";

//Начальные условия
