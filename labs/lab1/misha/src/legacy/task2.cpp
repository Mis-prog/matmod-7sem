#include <cmath>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>
#include <iostream>

struct Constants {
    static constexpr double G = 6.67e-11; // гравитационная постоянная
    static constexpr double M1 = 2.0e30; // масса звезды (кг)
    static constexpr double M2 = 6.4e23; // масса планеты (кг)
    static constexpr double M3 = 1.1e16; // масса астероида (кг)
    static constexpr double R1 = 696340e3; // радиус звезды (м)
    static constexpr double R2 = 3390e3; // радиус планеты (м)
    static constexpr double R3 = 11.1e3; // радиус астероида (м)
    static constexpr double R12 = 228e9; // начальное расстояние звезда-планета (м)
    static constexpr double R23 = 9.4e6; // начальное расстояние планета-астероид (м)
    static constexpr double U2 = 24e3; // начальная скорость планеты (м/с)
    static constexpr double U3 = 2.14e3; // начальная скорость астероида (м/с)

    static constexpr double U_oil = 3040; //скорость истечения
    static constexpr double H = 200e3; // высота орбиты
    static constexpr double T = 1200; // время работы ракеты
    static constexpr double M0 = 10; // полезная нагрузка
};

struct Params {
    double G; // гравитационная постоянная
    double H; // высота орбиты
    double U4; // начальная скорость ракеты (первая космическая)
    double U_oil; // скорость истечения
    double T; // время полного израсходования топлива
    double M0; // масса полезной нагрузки
    double Mt; // масса топлива
    double Mk; // масса конструкции ракеты
    double M1, M2, M3; // массы объектов
    double r12x, r12y; // координаты второго объекта
    double r13x, r13y; // координаты третьего объекта
    double alpha; // начальный угол
};

double m(double t, const Params *p) {
    if (t < p->T) {
        return (p->M0 + p->Mt) / (1 - p->alpha) - (p->Mt * t) / p->T;
    } else {
        return p->M0;
    }
}

double dm(double t, const Params *p) {
    if (t < p->T) {
        return -p->Mt / p->T;
    } else {
        return 0.0;
    }
}

int system_of_odes(double t, const double y[], double f[], void *params) {
    Params *p = static_cast<Params *>(params);

    // Переменные состояния
    double rx = y[0], ry = y[1], r13x = y[2], r13y = y[3];
    double vx = y[4], vy = y[5], v3x = y[6], v3y = y[7];

    // Вычисление расстояний
    double r = std::sqrt(rx * rx + ry * ry);
    double r2 = std::sqrt((rx - p->r12x) * (rx - p->r12x) + (ry - p->r12y) * (ry - p->r12y));
    double r3 = std::sqrt((rx - r13x) * (rx - r13x) + (ry - r13y) * (ry - r13y));
    double r13 = std::sqrt(r13x * r13x + r13y * r13y);
    double r23x = r13x - p->r12x;
    double r23y = r13y - p->r12y;
    double r23 = std::sqrt(r23x * r23x + r23y * r23y);

    // Вычисление массы и ее изменения
    double mass = m(t, p);
    double dmass = dm(t, p);
    double velocity = std::sqrt(vx * vx + vy * vy);

    // Уравнения движения
    f[0] = vx;
    f[1] = vy;
    f[2] = v3x;
    f[3] = v3y;

    f[4] = (-p->U_oil * dmass / mass * vx / velocity +
            p->G * (-p->M1 * rx / std::pow(r, 3) - p->M2 * (rx - p->r12x) / std::pow(r2, 3) - p->M3 * (rx - r13x) /
                    std::pow(r3, 3)));

    f[5] = (-p->U_oil * dmass / mass * vy / velocity +
            p->G * (-p->M1 * ry / std::pow(r, 3) - p->M2 * (ry - p->r12y) / std::pow(r2, 3) - p->M3 * (ry - r13y) /
                    std::pow(r3, 3)));

    f[6] = p->G * (-(p->M1 / std::pow(r13, 3)) * r13x - (p->M2 / std::pow(r23, 3)) * r23x);
    f[7] = p->G * (-(p->M1 / std::pow(r13, 3)) * r13y - (p->M2 / std::pow(r23, 3)) * r23y);

    return GSL_SUCCESS;
}

int main() {
    double y[8] = {
        Constants::R1 + Constants::R12 + Constants::R2, 0, 0,
        Constants::U2, // астероид (x,vx,y,vy)
        Constants::R1 + Constants::R12 + 2 * Constants::R2 + Constants::R23 + Constants::R3, 0, 0,
        Constants::U3 + Constants::U2 // ракета (x,vx,y,vy)
    };

    Params params = {
        .G = Constants::G,
        .M1 = Constants::M1,
        .M2 = Constants::M2,
        .M3 = Constants::M3,

        .M0 = Constants::M0,
        .U_oil = Constants::U_oil,
        .T = Constants::T,
        .H = Constants::H
    };


    gsl_odeiv2_system sys = {system_of_odes, nullptr, 8, &params};

    gsl_odeiv2_driver *driver = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk8pd, 1e-6, 1e-6, 0.0);

    double t = 0.0;
    double t1 = 1e5;

    while (t < t1) {
        int status = gsl_odeiv2_driver_apply(driver, &t, t1, y);

        if (status != GSL_SUCCESS) {
            std::cerr << "Error, return value= " << status << std::endl;
            break;
        }

        std::cout << "t = " << t << " x = " << y[0] << " y = " << y[1] << " vx = " << y[4] << " vy = " << y[5] <<
                std::endl;
    }

    gsl_odeiv2_driver_free(driver);
    return 0;
}
