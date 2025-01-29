#include <iostream>
#include <stdio.h>
#include <vector>
#include <iomanip>
#include <fstream>
#include <time.h>
#include <cmath>

using namespace std;

const double xi = 0.1931833275037836;
double alpha = 0.;
double my_beta = 760.;
const double m = 1;
const int I = 12; // номер варианта

template <typename T>
vector<T> add_vector2D(vector<T> a, vector<T> b) {
    vector<T> result = vector<T>(a.size());

    if (a.size() != b.size())
        throw "the sizes dont match";

    for (auto i = 0; i < a.size(); i++) {
        result[i] = a[i] + b[i];
    }

    return result;
}


long double FPU(long double r)
{
    return r * r * (0.5 + r * (alpha / 3 + my_beta * r * 0.25));
}

long double FPU_der(long double r)
{
    return r * (1 + r * (alpha + my_beta * r));
}

void write_in_txt(vector<long double> a, ofstream& out) {

    for (auto o : a) {
        out << o << " ";
    }

    out << endl;
}

void verle(vector<long double>& q_cur,
           vector<long double>& v_cur,
           long long k,
           long long n,
           long double dt,
           long long write_count,
           ofstream& out_q,
           ofstream& out_v)
{
    cout << "COMMON VERLE METHOD" << endl;

    // как часто будем записывать в файл
    long long write_step = k / write_count;

    // считаем начальный потенциал
    long double h0 = 0;
    for (int j = 1; j < 2 * n + 1; j++) {
        h0 += FPU(q_cur[j] - q_cur[j - 1]);
        h0 += 0.5 * v_cur[j] * v_cur[j];
    }

    cout << "t = " << setw(10) << setprecision(10) << 0 << ": H0 = " << h0 << endl;

    // определяем a
    vector<long double> a_cur = std::vector<long double>(v_cur);

    for (int i = 1; i < 2 * n + 1; i++) {
        a_cur[i] = -(1 / m) * (-FPU_der(q_cur[i + 1] - q_cur[i]) + FPU_der(q_cur[i] - q_cur[i - 1]));
    }
    a_cur[0] = a_cur[2 * n];
    a_cur[2 * n + 1] = a_cur[1];

    int order_i = 0;
    long double h;

    // Собственно сам алгоритм

    for (long long i = 0; i < k; i++) {

        for (int j = 0; j < 2 * n + 2; j++) {
            q_cur[j] = q_cur[j] + dt * (v_cur[j] + 0.5 * a_cur[j] * dt);
        }

        for (int j = 0; j < 2 * n + 2; j++) {
            v_cur[j] = v_cur[j] + 0.5 * a_cur[j] * dt;
        }

        for (int j = 1; j < 2 * n + 1; j++) {
            a_cur[j] = - (1 / m) * (-FPU_der(q_cur[j + 1] - q_cur[j]) + FPU_der(q_cur[j] - q_cur[j - 1]));
        }

        a_cur[0] = a_cur[2 * n];
        a_cur[2 * n + 1] = a_cur[1];

        for (int j = 0; j < 2 * n + 2; j++) {
            v_cur[j] = v_cur[j] + 0.5 * a_cur[j] * dt;
        }

        // запись в файл
        if ((i) % 10 == 0 and i < 100000) {
            write_in_txt(q_cur, out_q);
            write_in_txt(v_cur, out_v);
        }

        // на степенях десятки считаем и выводим погрешность
        if (int((i + 1) * dt) == pow(10, order_i)) {

            h = 0;

            for (int j = 1; j < 2 * n + 1; j++) {
                h += FPU(q_cur[j] - q_cur[j - 1]);
                h += 0.5 * v_cur[j] * v_cur[j] * m;
            }

            cout << "t = " << setw(10) << setprecision(10) << dt * (i + 1) << ": H" << (order_i + 1) << " = " << h << "\t";
            cout << "dH = " << setw(10) << setprecision(10) << abs(h0 - h) << endl;

            order_i += 1;
        }
    }
}

void simplex(vector<long double>& q_cur,
             vector<long double>& v_cur,
             long long k,
             long long n,
             long double dt,
             long long write_count,
             ofstream& out_q,
             ofstream& out_v)
{
    cout << "SIMPLEX VERLE METHOD" << endl;
    // считаем начальный потенциал
    long double h0 = 0;

    for (int j = 1; j < 2 * n + 1; j++) {
        h0 += FPU(q_cur[j] - q_cur[j - 1]);
        h0 += 0.5 * v_cur[j] * v_cur[j];
    }

    cout << "t = " << setw(10) << setprecision(10) << 0 << ": H0 = " << h0 << endl;


    // определяем a
    vector<long double> a_cur = std::vector<long double>(v_cur);

    for (int i = 1; i < 2 * n + 1; i++) {
        a_cur[i] = -(1 / m) * (-FPU_der(q_cur[i + 1] - q_cur[i]) + FPU_der(q_cur[i] - q_cur[i - 1]));
    }
    a_cur[0] = a_cur[2 * n];
    a_cur[2 * n + 1] = a_cur[1];

    int order_i = 0;
    long double h;

    // сам алгоритм
    for (int i = 0; i < k; i++)
    {
        // шаг 1
        for (int j = 1; j < q_cur.size() - 1; j++)
        {
            q_cur[j] = q_cur[j] + v_cur[j] * dt * xi;
        }
        q_cur[0] = q_cur[2 * n]; q_cur[2 * n + 1] = q_cur[1];

        // шаг 2
        for (int j = 1; j < 2 * n + 1; j++) {
            a_cur[j] = -(1 / m) * (-FPU_der(q_cur[j + 1] - q_cur[j]) + FPU_der(q_cur[j] - q_cur[j - 1]));
        }
        a_cur[0] = a_cur[2 * n];
        a_cur[2 * n + 1] = a_cur[1];

        // шаг 3
        for (int j = 1; j < q_cur.size() - 1; j++)
        {
            v_cur[j] = v_cur[j] + 0.5 * a_cur[j] * dt;
        }
        v_cur[0] = v_cur[2 * n]; v_cur[2 * n + 1] = v_cur[1];

        // шаг 4
        for (int j = 1; j < q_cur.size() - 1; j++)
        {
            q_cur[j] = q_cur[j] + v_cur[j] * dt * (1 - 2 * xi);
        }
        q_cur[0] = q_cur[2 * n]; q_cur[2 * n + 1] = q_cur[1];

        // шаг 5
        for (int j = 1; j < 2 * n + 1; j++) {
            a_cur[j] = -(1 / m) * (-FPU_der(q_cur[j + 1] - q_cur[j]) + FPU_der(q_cur[j] - q_cur[j - 1]));
        }
        a_cur[0] = a_cur[2 * n];
        a_cur[2 * n + 1] = a_cur[1];

        // шаг 6
        for (int j = 1; j < q_cur.size() - 1; j++)
        {
            v_cur[j] = v_cur[j] + 0.5 * a_cur[j] * dt;
        }
        v_cur[0] = v_cur[2 * n]; v_cur[2 * n + 1] = v_cur[1];

        // шаг 7
        for (int j = 1; j < q_cur.size() - 1; j++)
        {
            q_cur[j] = q_cur[j] + v_cur[j] * xi * dt;
        }
        q_cur[0] = q_cur[2 * n]; q_cur[2 * n + 1] = q_cur[1];

        // запись в файл
        if ((i) % 5 == 0 and i < 100000) {
            write_in_txt(q_cur, out_q);
            write_in_txt(v_cur, out_v);
        }

        // на степенях десятки считаем и выводим погрешность
        if (int((i + 1) * dt) == pow(10, order_i)) {

            h = 0;

            for (int j = 1; j < 2 * n + 1; j++) {
                h += FPU(q_cur[j] - q_cur[j - 1]);
                h += 0.5 * v_cur[j] * v_cur[j] * m;
            }

            cout << "t = " << setw(10) << setprecision(10) << dt * (i + 1) << ": H" << (order_i + 1) << " = " << h << "\t";
            cout << "dH = " << setw(10) << setprecision(10) << abs(h0 - h) << endl;

            order_i += 1;
        }
    }
}


int main() {
    int n = 100; // количество частиц
    //long double q_start = 0.01 * I; // init начальных значений
    long double q_start = 0.02;
    //long double q_start = 0.5; // init начальных значений


    long long N = 1e6; // количество временных шагов
    long double dt = 0.01; // шаг по времени
    long double t_end = dt * N; // конечное время
    ofstream out_q_verle("../labs/lab3/example/result/out_verle_q.txt");
    ofstream out_v_verle("../labs/lab3/example/result/out_verle_v.txt");
//    ofstream out_q_simplex("../labs/lab3/example/result/out_simplex_q.txt");
//    ofstream out_v_simplex("../labs/lab3/example/result/out_simplex_v.txt");


    vector<long double> q_cur = vector<long double>(2 * n + 2);
    vector<long double> v_cur = vector<long double>(2 * n + 2);

    cout << "N = " << N << endl;
    cout << "n = " << n << endl;
    cout << "dt = " << dt << endl;
    // начальные значения для q (остальные 0)

    q_cur[n] = -q_start;
    q_cur[n + 1] = q_start;
//    out_q_verle << alpha << " " << my_beta << endl;
//    out_v_verle << alpha << " " << my_beta << endl;
    verle(q_cur, v_cur, N, n, dt, 1000, out_q_verle, out_v_verle);

//    q_cur = vector<long double>(2 * n + 2);
//    v_cur = vector<long double>(2 * n + 2);
//    //начальные значения для q (остальные 0)
//    q_cur[n] = -q_start;
//    q_cur[n + 1] = q_start;
//    out_q_simplex << alpha << " " << my_beta << endl;
//    out_v_simplex << alpha << " " << my_beta << endl;
//    simplex(q_cur, v_cur, N, n, dt, 1000, out_q_simplex, out_v_simplex);

    return 0;
}