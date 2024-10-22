#include <iostream>
#include <cmath>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>
#include <vector>
#include <iomanip>
#include <fstream>


const double G = 6.67e-11;
const double m1 = 2e30, m2 = 1.9e27, m3 = 1.5e23;

// Правые части системы уравнений
int func(double t, const double y[], double f[], void *params) {
    (void)(t); // Время не используется напрямую
    double r12x = y[0], r12y = y[1], r13x = y[2], r13y = y[3];
    double v2x = y[4], v2y = y[5], v3x = y[6], v3y = y[7];

    double r2 = sqrt(r12x * r12x + r12y * r12y);
    double r3 = sqrt(r13x * r13x + r13y * r13y);

    double r23x = r13x - r12x;
    double r23y = r13y - r12y;
    double r = sqrt(r23x * r23x + r23y * r23y);

    // Уравнения движения
    f[0] = v2x;
    f[1] = v2y;
    f[2] = v3x;
    f[3] = v3y;
    f[4] = G * (-(m1 / (r2 * r2 * r2) * r12x) + (m3 / (r * r * r) * r23x));
    f[5] = G * (-(m1 / (r2 * r2 * r2) * r12y) + (m3 / (r * r * r) * r23y));
    f[6] = G * (-(m1 / (r3 * r3 * r3) * r13x) - (m2 / (r * r * r) * r23x));
    f[7] = G * (-(m1 / (r3 * r3 * r3) * r13y) - (m2 / (r * r * r) * r23y));

    return GSL_SUCCESS;
}

// Событие для остановки при заданном условии
int event_stop(double t, const double y[], double* f, void* params) {
    // Условие остановки — это пересечение траекторий, похожее на то, что было в Python.
    double threshold = 0.001; // Можем добавить условие для остановки симуляции
    double r1 = sqrt(y[0] * y[0] + y[1] * y[1]);
    double r2 = sqrt(y[2] * y[2] + y[3] * y[3]);

    return (r1 - r2) < threshold ? GSL_SUCCESS : GSL_CONTINUE;
}

int main() {
    // Инициализация условий и системы
    gsl_odeiv2_system sys = {func, nullptr, 8, nullptr};

    gsl_odeiv2_driver *driver = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk8pd,
                                                              1e-6, 1e-6, 0.0);

    // Начальные условия
    double t = 0.0;
    double t1 = 3000000000.0;
    double y[8] = {780e9, 0, 780e9 + 1070e6, 0, 0, 13000, 0, 23900};

    std::vector<std::vector<double>> orbitsx(1200), orbitsy(1200);
    std::ofstream fout("../../../../../labs/lab1/misha/plot/output_2.txt");
    fout << "x2 y2\n";

    for (int i = 0; i < 2000; ++i) {
        double ti = i * 1000.0;
        int status = gsl_odeiv2_driver_apply(driver, &t, ti, y);

        if (status != GSL_SUCCESS) {
            std::cerr << "Error: " << gsl_strerror(status) << std::endl;
            break;
        }

        // Записываем орбиты
        orbitsx[i].push_back(y[2] - y[0]);
        orbitsy[i].push_back(y[3] - y[1]);

        // Печатаем результат
        fout << y[0] << " " << y[1] << std::endl;
    }

    // Освобождаем память
    gsl_odeiv2_driver_free(driver);

    return 0;
}

