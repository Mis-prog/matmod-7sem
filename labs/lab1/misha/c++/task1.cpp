#include <iostream>
#include <cmath>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>
#include <fstream>
#include <string>
#include <stdexcept>
#include <vector>
#include "progressbar.hpp"
#include "omp.h"

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


bool check_interseption_start_point(double x, double y) {
    double dx = Constants::R1 + Constants::R12 + 2 * Constants::R2 + Constants::R23 + Constants::R3 - x;
    double dy = y;
    return (sqrt(dx * dx + dy * dy) < 2e8);
}

bool check_interseption_start_point_udp2(double x, double y) {
    double initial_x = Constants::R1 + Constants::R12 + 2 * Constants::R2 + Constants::R23 + Constants::R3;
    double initial_y = 0.0;

    double angle_current = std::atan2(y, x);
    double angle_initial = std::atan2(initial_y, initial_x);

    double angle_difference = std::abs(angle_current - angle_initial);
    if (angle_difference > M_PI) angle_difference = 2 * M_PI - angle_difference;


    double threshold_angle = 0.000085 * M_PI;

    return angle_difference < threshold_angle;
}


int main() {
    double y[8] = {
            Constants::R1 + Constants::R12 + Constants::R2, 0, 0,
            Constants::U2,                    // планета (x,vx,y,vy)
            Constants::R1 + Constants::R12 + 2 * Constants::R2 + Constants::R23 + Constants::R3, 0, 0,
            Constants::U3 + Constants::U2  // астероид (x,vx,y,vy)
    };

    double t = 0.0;// начальное время
    double t_circle_end = 60. * 60 * 24 * 365;
    double t_end = t_circle_end * 2000;     // конечное время
    double h = 5000;// шаг интегрирования
    int n_steps = round(t_end / h);// количество шагов
    progressbar bar(n_steps);
    bar.set_todo_char(" ");
    bar.set_done_char("█");
    bar.set_opening_bracket_char("{");
    bar.set_closing_bracket_char("}");


    std::cout << "Шаг: " << h << ",кол-во шагов: " << n_steps << std::endl;

    gsl_odeiv2_system sys = {Physics::calculateForces, nullptr, 8, nullptr};
    gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new(
            &sys, gsl_odeiv2_step_rk8pd, h, 1e-6, 1e-6
    );

    std::ofstream fout_main("../../../../../labs/lab1/misha/plot/path_full.csv");
    std::ofstream fout_log("../../../../../labs/lab1/misha/plot/log.txt");

    fout_main << "x2 y2 x3 y3\n";

    std::vector<std::vector<double>> orbitsX_Spytnik(4);
    std::vector<std::vector<double>> orbitsY_Spytnik(4);

    std::vector<std::vector<double>> orbitsX_Planeta(4);
    std::vector<std::vector<double>> orbitsY_Planeta(4);


    double t_curr = t;
    int curr_i = 0;
    int prev_i = 0;
    int count = 0;
    while (t_curr < t_end) {
        bar.update();

        int status = gsl_odeiv2_driver_apply(d, &t_curr, t_curr + h, y);

        if (status != GSL_SUCCESS) {
            throw std::runtime_error("Ошибка интегрирования: " +
                                     std::string(gsl_strerror(status)));
        }

        if (check_interseption_start_point_udp2(y[4], y[6])) {
            fout_log << "Пересечение на шаге " << curr_i << std::endl;
            count++;
        }

        switch (count) {
            case 1:
                orbitsX_Spytnik[0].push_back(y[4]);
                orbitsY_Spytnik[0].push_back(y[6]);

                orbitsX_Planeta[0].push_back(y[0]);
                orbitsY_Planeta[0].push_back(y[2]);
                break;
            case 100:
                orbitsX_Spytnik[1].push_back(y[4]);
                orbitsY_Spytnik[1].push_back(y[6]);

                orbitsX_Planeta[1].push_back(y[0]);
                orbitsY_Planeta[1].push_back(y[2]);
                break;
            case 500:
                orbitsX_Spytnik[2].push_back(y[4]);
                orbitsY_Spytnik[2].push_back(y[6]);

                orbitsX_Planeta[2].push_back(y[0]);
                orbitsY_Planeta[2].push_back(y[2]);
                break;
            case 1000:
                orbitsX_Spytnik[3].push_back(y[4]);
                orbitsY_Spytnik[3].push_back(y[6]);

                orbitsX_Planeta[3].push_back(y[0]);
                orbitsY_Planeta[3].push_back(y[2]);
                break;
        }

        fout_main <<
                  y[0] << " " << y[2] << " "  // координаты планеты
                  << y[4] << " " << y[6] << "\n"; // координаты астероида

        curr_i++;
    }

    gsl_odeiv2_driver_free(d);
    fout_main.close();
    fout_log.close();


    for (int i = 0; i < orbitsX_Spytnik.size(); i++) {
        std::ofstream fout_spytnik("../../../../../labs/lab1/misha/plot/path_spytnik_" + std::to_string(i) + ".csv");
        fout_spytnik << "x y\n";
        for (int j = 0; j < orbitsX_Spytnik[i].size(); j++) {
            fout_spytnik << orbitsX_Spytnik[i][j] << " " << orbitsY_Spytnik[i][j] << std::endl;
        }
        fout_spytnik.close();
    }

    for (int i = 0; i < orbitsX_Spytnik.size(); i++) {
        std::ofstream fout_planeta("../../../../../labs/lab1/misha/plot/path_planeta_" + std::to_string(i) + ".csv");
        fout_planeta << "x y\n";
        for (int j = 0; j < orbitsX_Planeta[i].size(); j++) {
            fout_planeta << orbitsX_Planeta[i][j] << " " << orbitsY_Planeta[i][j] << std::endl;
        }
        fout_planeta.close();
    }

    orbitsX_Spytnik.clear();
    orbitsY_Spytnik.clear();

    std::cout << "\nРасчет успешно завершен\n" << "Количество пересечений " << count;
//    <<"\nКол-во шагов: " << curr_i << "\nВремя: " << t_end - t_curr;
    return 0;
}