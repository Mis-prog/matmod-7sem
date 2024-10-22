#include <iostream>
#include <cmath>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>
#include <fstream>
#include <string>
#include <stdexcept>

// Физические константы
struct Constants {
    static constexpr double G = 6.67e-11;    // гравитационная постоянная
    static constexpr double M1 = 2.0e30;     // масса звезды (кг)
    static constexpr double M2 = 6.4e23;     // масса планеты (кг)
    static constexpr double M3 = 1.1e16;     // масса астероида (кг)
    static constexpr double R1 = 696340e3;   // радиус звезды (м)
    static constexpr double R2 = 3390e3;     // радиус планеты (м)
    static constexpr double R3 = 11.1e3;     // радиус астероида (м)
    static constexpr double R12 = 228e9;     // начальное расстояние звезда-планета (м)
    static constexpr double R23 = 9.4e6;     // начальное расстояние планета-астероид (м)
    static constexpr double U2 = 24e3;       // начальная скорость планеты (м/с)
    static constexpr double U3 = 2.14e3;     // начальная скорость астероида (м/с)
};

// Вспомогательные функции
class Physics {
public:
    static double distance(double x1, double y1, double x2, double y2) {
        double dx = x2 - x1;
        double dy = y2 - y1;
        return std::sqrt(dx * dx + dy * dy);
    }

    static int calculateForces(double t, const double y[], double f[], void *params) {
        try {
            // Координаты и скорости планеты и астероида
            // Солнце находится в начале координат (0,0)
            double x2 = y[0], vx2 = y[1], y2 = y[2], vy2 = y[3];  // планета
            double x3 = y[4], vx3 = y[5], y3 = y[6], vy3 = y[7];  // астероид

            // Расчёт расстояний
            double r_12 = std::sqrt(x2 * x2 + y2 * y2);  // расстояние от солнца до планеты
            double r_13 = std::sqrt(x3 * x3 + y3 * y3);  // расстояние от солнца до астероида
            double r_23 = distance(x2, y2, x3, y3);  // расстояние между планетой и астероидом

            // Скорости
            f[0] = vx2;
            f[2] = vy2;  // планета
            f[4] = vx3;
            f[6] = vy3;  // астероид

            // Ускорения для планеты (тело 2)
            f[1] = -Constants::G * Constants::M1 * x2 / std::pow(r_12, 3) +
                   Constants::G * Constants::M3 * (x3 - x2) / std::pow(r_23, 3);

            f[3] = -Constants::G * Constants::M1 * y2 / std::pow(r_12, 3) +
                   Constants::G * Constants::M3 * (y3 - y2) / std::pow(r_23, 3);

            // Ускорения для астероида (тело 3)
            f[5] = -Constants::G * Constants::M1 * x3 / std::pow(r_13, 3) +
                   Constants::G * Constants::M2 * (x2 - x3) / std::pow(r_23, 3);

            f[7] = -Constants::G * Constants::M1 * y3 / std::pow(r_13, 3) +
                   Constants::G * Constants::M2 * (y2 - y3) / std::pow(r_23, 3);

            return GSL_SUCCESS;
        }
        catch (const std::exception &e) {
            std::cerr << "Ошибка: " << e.what() << std::endl;
            return GSL_FAILURE;
        }
    }
};

bool checkFullOrbit(double x_curr, double y_curr) {
    double dx = std::abs(x_curr - (Constants::R1 + Constants::R12 + Constants::R2));
    double dy = std::abs(y_curr);
    return (dx < 100 && dy < 100);
}

int main() {
    double y[8] = {
            Constants::R1 + Constants::R12 + Constants::R2, 0, 0,
            Constants::U2,                    // планета (x,vx,y,vy)
            Constants::R1 + Constants::R12 + 2 * Constants::R2 + Constants::R23 + Constants::R3, 0, 0,
            Constants::U3 + Constants::U2  // астероид (x,vx,y,vy)
    };

    // Параметры интегрирования
    double t = 0.0;                           // начальное время
    double t_end = 60. * 60 * 24 * 365 ;     // конечное время (1000 лет)
    int n_steps = 1000;                   // количество шагов
    double h = t_end / n_steps;// шаг интегрирования

    std::cout << "Шаг: " << h << std::endl;

    gsl_odeiv2_system sys = {Physics::calculateForces, nullptr, 8, nullptr};
    gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new(
            &sys, gsl_odeiv2_step_rkf45, h, 1e-6, 1e-6
    );

    // Открытие файла для записи результатов
    std::ofstream fout("../../../../../labs/lab1/misha/plot/output.csv");
    if (!fout.is_open()) {
        throw std::runtime_error("Не удалось открыть файл для записи");
    }

    fout << "x2 y2 x3 y3\n";  // заголовок файла

    // Основной цикл интегрирования
    double t_curr = t;
    int i = 0;
    while (t_curr < t_end) {
        int status = gsl_odeiv2_driver_apply(d, &t_curr, t_curr + h, y);

        if (status != GSL_SUCCESS) {
            throw std::runtime_error("Ошибка интегрирования: " +
                                     std::string(gsl_strerror(status)));
        }

        if (checkFullOrbit(y[0], y[2])) {
            std::cout << "Полная орбита астероида достигнута на  " << i << " шаге\n";
            break;
        }


        fout <<
             y[0] << " " << y[2] << " "  // координаты планеты
             << y[4] << " " << y[6] << "\n"; // координаты астероида
        t_curr += h;
        i++;
    }

    // Освобождение ресурсов
    gsl_odeiv2_driver_free(d);
    fout.close();

    std::cout << "Расчет успешно завершен\n";
    return 0;
}