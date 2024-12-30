#include <iostream>
#include <cmath>
#include <boost/numeric/odeint.hpp>
#include <fstream>
#include <chrono>
#include <iomanip>

using namespace boost::numeric::odeint;
using state_type = std::array<double, 8>;
using namespace std;

struct Constants {
    static constexpr double G = 6.67e-11;    // гравитационная постоянная
    static constexpr double M1 = 2.0e30;     // масса звезды (кг)
    static constexpr double M2 = 6.0e24;     // масса планеты (кг)
    static constexpr double M3 = 7.3e22;     // масса астероида (кг)
    static constexpr double R1 = 696340e3;   // радиус звезды (м)
    static constexpr double R2 = 6378e3;     // радиус планеты (м)
    static constexpr double R3 = 1737e3;     // радиус астероида (м)
    static constexpr double R12 = 150e9;     // начальное расстояние звезда-планета (м)
    static constexpr double R23 = 384e6;     // начальное расстояние планета-астероид (м)
    static constexpr double U2 = 30e3;       // начальная скорость планеты (м/с)
    static constexpr double U3 = 1e3;     // начальная скорость астероида (м/с)
};

// Вспомогательные функции
class Physics {
public:
    static bool intersection;
    static bool fast_intersection;

    static double distance(double x1, double y1, double x2, double y2) {
        double dx = x2 - x1;
        double dy = y2 - y1;
        return std::sqrt(dx * dx + dy * dy);
    }

    static void calculateForces(const state_type &y, state_type &f, double /* t */) {
        // Координаты и скорости планеты и астероида
        // Солнце находится в начале координат (0,0)
        double x2 = y[0], vx2 = y[1], y2 = y[2], vy2 = y[3]; // планета
        double x3 = y[4], vx3 = y[5], y3 = y[6], vy3 = y[7]; // астероид

        // Расчёт расстояний
        double r_12 = std::sqrt(x2 * x2 + y2 * y2); // расстояние от солнца до планеты
        double r_13 = std::sqrt(x3 * x3 + y3 * y3); // расстояние от солнца до астероида
        double r_23 = distance(x2, y2, x3, y3); // расстояние между планетой и астероидом

        // Скорости
        f[0] = vx2;
        f[2] = vy2; // планета
        f[4] = vx3;
        f[6] = vy3; // астероид

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
    }

    // нормализация
    static pair<double, double> normalize_vector(double x, double y) {
        double length = sqrt(x * x + y * y);
        return {x / length, y / length};
    }

    // проверка на коллианеарность
    static vector<double> coordinates_on_line(double x2, double y2, double x3, double y3) {
        auto planeta_norm = normalize_vector(x2, y2);
        auto spytnik_norm = normalize_vector(x3, y3);

        double vector_product = planeta_norm.first * spytnik_norm.second - planeta_norm.second * spytnik_norm.first;
        if (vector_product < 1e-7 and abs(x3) > abs(x2)) return vector<double>{x2, y2, x3, y3};
        else return vector<double>{};
    }

    static double coordinates_on_line_coef(double x2, double y2, double x3, double y3) {
        auto planeta_norm = normalize_vector(x2, y2);
        auto spytnik_norm = normalize_vector(x3, y3);

        double planeta_coef = planeta_norm.second / spytnik_norm.first;
        double spytnik_coef = spytnik_norm.second / spytnik_norm.first;
        if (abs(planeta_coef - spytnik_coef) < 1e-5 and abs(x3) > abs(x2)) return abs(planeta_coef - spytnik_coef);
    }

    // проверка на пересечение
    static bool intersection_angle(double x, double y) {
        auto norm_current = normalize_vector(x, y);
        if (norm_current.first >= 0.99 and !intersection) {
            intersection = true;
            return true;
        } else if (norm_current.first < 0 and intersection) {
            intersection = false;
        }
        return false;
    }

    // вычисление нач координат для второй задачи
    static bool fast_crossing_check(double x, double y) {
        auto norm_current = normalize_vector(x, y);
        if (norm_current.first >= 0.97)
            return true;
        return false;
    }
};

bool Physics::intersection = false;


int main() {
    state_type y = {
        Constants::R1 + Constants::R12 + Constants::R2, 0, 0,
        Constants::U2, // планета (x,vx,y,vy)
        Constants::R1 + Constants::R12 + 2 * Constants::R2 + Constants::R23 + Constants::R3, 0, 0,
        Constants::U3 + Constants::U2 // астероид (x,vx,y,vy)
    };

    double t = 0.0;
    double t_circle_end = 60. * 60 * 24 * 365;
    double t_end = t_circle_end * 1;
    double h = 5000;

    runge_kutta_cash_karp54<state_type> stepper;

    std::ofstream fout_main("../labs/lab1/alia/result/task1/path_full.csv");
    fout_main << "x2 y2 x3 y3\n";

    std::vector<std::vector<double> > orbitsX_Spytnik(4);
    std::vector<std::vector<double> > orbitsY_Spytnik(4);

    std::vector<std::vector<double> > orbitsX_Planeta(4);
    std::vector<std::vector<double> > orbitsY_Planeta(4);

    int curr_i = 0;
    int count = 0;
    double t_curr = t;
    double point_collinear;

    double y_min = 100;
    while (t_curr < t_end) {
        stepper.do_step(Physics::calculateForces, y, t, h);

        if (Physics::intersection_angle(y[0], y[2])) {
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
                // h = 1000; // меняем шаг для уточнения нач данных для второй задачиw
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

        // if (count > 1000) {
        //     point_collinear = Physics::coordinates_on_line_coef(y[0], y[2], y[4], y[6]);
        //     if (point_collinear && point_collinear < y_min) {
        //         y_min = point_collinear;
        //         cout << endl << "Возможное решение на: " << count - 1 << " шаге" << " ,коэффицент: " << fixed <<
        //                 setprecision(10) << point_collinear << endl;
        //         for (auto value: y) {
        //             cout << fixed << setprecision(2) << value << " ";
        //         }
        //     }
        // }


        fout_main <<
                y[0] << " " << y[2] << " " // координаты планеты
                << y[4] << " " << y[6] << "\n"; // координаты астероида
        t_curr += h;
        curr_i++;
    }
    fout_main.close();

    cout << "Кол-во пересечений: " << count - 1 << endl;

    for (int i = 0; i < orbitsX_Spytnik.size(); i++) {
        std::ofstream fout_spytnik(
            "../labs/lab1/alia/result/task1/path_spytnik_" + std::to_string(i) + ".csv");
        fout_spytnik << "x y\n";
        for (int j = 0; j < orbitsX_Spytnik[i].size(); j++) {
            fout_spytnik << orbitsX_Spytnik[i][j] << " " << orbitsY_Spytnik[i][j] << std::endl;
        }
        fout_spytnik.close();
    }

    for (int i = 0; i < orbitsX_Spytnik.size(); i++) {
        std::ofstream fout_planeta(
            "../labs/lab1/alia/result/task1/path_planeta_" + std::to_string(i) + ".csv");
        fout_planeta << "x y\n";
        for (int j = 0; j < orbitsX_Planeta[i].size(); j++) {
            fout_planeta << orbitsX_Planeta[i][j] << " " << orbitsY_Planeta[i][j] << std::endl;
        }
        fout_planeta.close();
    }

    orbitsX_Spytnik.clear();
    orbitsY_Spytnik.clear();
    return 0;
}
