#include <iostream>
#include <cmath>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>
#include <fstream>
#include <gsl/gsl_matrix.h>

using namespace std;

const double G = 6.67e-11;
const double M1 = 2.0e30; // масса тела 1 (звезда), кг
const double M2 = 6.4e23; // масса тела 2 (планета), кг
const double M3 = 1.1e16; // масса астероида
const double R1 = 696340e3; // радиус солнца
const double R2 = 3390e3; //  радиус тела 2 (планета), км
const double R12 = 228e9; // расстояние между телом 1 и телом 2, км
const double U2 = 24e3; // начальная скорость тела 2, км/с
const double R3 = 11.1e3; // радиус тела 3 (астероид), км
const double R23 = 9.4e6; // расстояние между телом 2 и телом 3, км
const double U3 = 2.14e3; //  начальная скорость тела 3, км/с

double r12(double x1, double y1, double x2, double y2) {
    double diff_x = x2 - x1;
    double diff_y = y2 - y1;
    return sqrt(diff_x * diff_x + diff_y * diff_y);
}

double r13(double x1, double y1, double x3, double y3) {
    double diff_x = x3 - x1;
    double diff_y = y3 - y1;
    return sqrt(diff_x * diff_x + diff_y * diff_y);
}

double r23(double x2, double y2, double x3, double y3) {
    double diff_x = x3 - x2;
    double diff_y = y3 - y2;
    return sqrt(diff_x * diff_x + diff_y * diff_y);
}

int system(double t, const double y[], double f[], void *params) {
    double x1 = y[0], vx1 = y[1];
    double y1 = y[2], vy1 = y[3];
    double x2 = y[4], vx2 = y[5];
    double y2 = y[6], vy2 = y[7];
    double x3 = y[8], vx3 = y[9];
    double y3 = y[10], vy3 = y[11];


    double r_12 = r12(x1, y1, x2, y2);
    double r_13 = r13(x1, y1, x3, y3);
    double r_23 = r23(x2, y2, x3, y3);

    if (r_12 == 0 || r_13 == 0 || r_23 == 0) {
        cerr << "Ошибка: деление на ноль!" << endl;
        return GSL_FAILURE;
    }

    f[0] = vx1;
    f[1] = G * M2 * (x2 - x1) / pow(r_12, 3) + G * M3 * (x3 - x1) / pow(r_13, 3); // ax1

    f[2] = vy1;
    f[3] = G * M2 * (y2 - y1) / pow(r_12, 3) + G * M3 * (y3 - y1) / pow(r_13, 3); //ay1

    f[4] = vx2;
    f[5] = G * M1 * (x1 - x2) / pow(r_12, 3) + G * M3 * (x3 - x2) / pow(r_23, 3); // ax2

    f[6] = vy2;
    f[7] = G * M1 * (y1 - y2) / pow(r_12, 3) + G * M3 * (y3 - y2) / pow(r_23, 3); //ay2

    f[8] = vx3;
    f[9] = G * M1 * (x1 - x3) / pow(r_13, 3) + G * M2 * (x2 - x3) / pow(r_23, 3); //ax3

    f[10] = vy3;
    f[11] = G * M1 * (y1 - y3) / pow(r_13, 3) + G * M2 * (y2 - y3) / pow(r_23, 3); //ay3

    return GSL_SUCCESS;
}


int main() {
    double y[12]; // Вектор состояния
    y[0] = 0;
    y[1] = 0;
    y[2] = 0;
    y[3] = 0;
    y[4] = (R1 + R12 + R2); // начальное положение x тела 2
    y[5] = 0; // начальная скорость по x тела 2
    y[6] = 0; // начальное положение y тела 2
    y[7] = U2; // начальная скорость по y тела 2 в м/с
    y[8] = (R1 + R12 + 2 * R2 + R23 + R3); // начальное положение x тела 3
    y[9] = 0; // начальная скорость по x тела 3
    y[10] = 0; // начальное положение y тела 3
    y[11] = (U3 + U2); // начальная скорость по y тела 3 в м/с


    double t = 0.0; // начальное время
    unsigned int t_end = 60 * 60 * 24 * 687 * 1000; // конечное время
    int n = 1000000;
    double h = t_end / n;
    cout << "h: " << h << endl;
//    cout << "turnover steps: " << (60 * 60 * 24 * 687)/h;
    gsl_odeiv2_system sys = {system, nullptr, 12, nullptr};
    gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk8pd, h, 1e-2, 1e-2);

    ofstream fout("../../../../../labs/lab1/misha/plot/output.txt");
    fout << "x1 y1 x2 y2 x3 y3\n";
    double t_curr = t;
    while (t_curr < t_end) {

        int status = gsl_odeiv2_driver_apply(d, &t_curr, (t_curr + h), y);

        if (status != GSL_SUCCESS) {
            cerr << "Ошибка интегрирования: " << gsl_strerror(status) << endl;
            break;
        }

        if ((int) (t_curr / 10000) % 10 == 0) {
            fout << y[0] << " " << y[2] << " " << y[4] << " " << y[6] << " " << y[8] << " " << y[10] << endl;
        }
        t_curr += h;
    }

    gsl_odeiv2_driver_free(d);

    return 0;
}
