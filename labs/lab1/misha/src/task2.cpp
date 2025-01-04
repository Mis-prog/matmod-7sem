#include <iostream>
#include <cmath>
#include <boost/numeric/odeint.hpp>
#include <fstream>
#include <chrono>

using namespace boost::numeric::odeint;
using state_type = std::array<double, 8>;
using namespace std;

// Физические константы
struct Constants {
    static constexpr double G = 6.67e-11; // гравитационная постоянная
    static constexpr double M1 = 2.0e30; // масса звезды (кг)
    static constexpr double M2 = 6.4e23; // масса планеты (кг)
    static constexpr double M3 = 1.1e16; // масса спутника (кг)
    static constexpr double R1 = 696340e3; // радиус звезды (м)
    static constexpr double R2 = 3390e3; // радиус планеты (м)
    static constexpr double R3 = 11.1e3; // радиус спутника (м)
    static constexpr double R12 = 228e9; // начальное расстояние звезда-планета (м)
    static constexpr double R23 = 9.4e6; // начальное расстояние планета-спутника (м)
    static constexpr double U2 = 24e3; // начальная скорость планеты (м/с)
    static constexpr double U3 = 2.14e3; // начальная скорость спутника (м/с)

    static constexpr double T = 1200.0; // время работы двигателя (с)
    static constexpr double H = 200e3; // высота орбиты (м)
    static constexpr double M0 = 10.0; // масса полезной назрузки (кг)
    static constexpr double U = 3040.0; // скорость истечения (м/c)
    static constexpr double koef = 0.025;
};

// Вспомогательные функции
class Physics {
public:
    static double mt;
    static double r12x; // координаты планеты
    static double r12y;


    static double distance(double x1, double y1, double x2, double y2) {
        double dx = x2 - x1;
        double dy = y2 - y1;
        return std::sqrt(dx * dx + dy * dy);
    }

    static double m(double t) {
        if (t >= Constants::T) {
            return Constants::M0;
        } else {
            return (Constants::M0 + mt) / (1 - Constants::koef) - (mt * t) / Constants::T;
        }
    }

    static double dm(double t) {
        if (t > Constants::T) {
            return 0.0;
        } else {
            return -mt / Constants::T;
        }
    }

    static void calculateForces(const state_type &y, state_type &f, double t) {
        // Координаты и скорости спутника и ракеты
    double rx = y[0], ry = y[1], vx = y[4], vy = y[5];  // ракета
    double r3x = y[2], r3y = y[3], v3x = y[6], v3y = y[7];  // спутник

    // Проверка столкновения ракеты со спутником
    if (sqrt((rx - r3x) * (rx - r3x) + (ry - r3y) * (ry - r3y)) <= Constants::R3) {
        throw std::runtime_error("Ракета попала в спутник.");
    }

    // Расчёт расстояний
    double r1 = std::sqrt(rx * rx + ry * ry);  // расстояние от ракеты до планеты
    double r2 = distance(rx, ry, r3x, r3y);    // расстояние между ракетой и спутником
    double r3 = std::sqrt(r3x * r3x + r3y * r3y);  // расстояние от спутника до планеты

    // Расчет скорости ракеты
    double v = std::sqrt(vx * vx + vy * vy);

    // Скорости
    f[0] = vx;
    f[1] = vy;
    f[2] = v3x;
    f[3] = v3y;

    // Ускорения для ракеты
    std::pair<double, double> f_rocket = {
    -(Constants::U * Physics::dm(t) * vx) / (v * Physics::m(t)),
    -(Constants::U * Physics::dm(t) * vy) / (v * Physics::m(t))
    };

    std::pair<double, double> f_gravity = {
        Constants::G * (-Constants::M2 * rx / std::pow(r1, 3) - Constants::M3 * (rx - r3x) / std::pow(r2, 3)),
        Constants::G * (-Constants::M2 * ry / std::pow(r1, 3) - Constants::M3 * (ry - r3y) / std::pow(r2, 3))
    };

    // Ускорения для ракеты
    f[4] = f_rocket.first + f_gravity.first;
    f[5] = f_rocket.second + f_gravity.second;


    // Ускорения для спутника
    f[6] = -Constants::G * Constants::M2 * r3x / std::pow(r3, 3);  // сила от планеты
    f[7] = -Constants::G * Constants::M2 * r3y / std::pow(r3, 3);  // сила от планеты
    }
};

double Physics::mt;
double Physics::r12x;
double Physics::r12y;

int main() {
    double r12x0 = -50368219856.43, r12y0 = -219503615669.65,
            v2x0=23688.42, v2y0=-5739.71,
            r13x0 = -50370393498.29, r13y0 = -219513089385.80,
            v3x0 = 25513.15, v3y0 = -6666.22; // нач координаты планеты и спутника

    double angle_input;
    cout << "Введите общую массу и угол: \n";
    cin >> Physics::mt >> angle_input;
    double angle = angle_input * M_PI / 180;

    // Переводим спутник в относительные координаты
    double r3x = r13x0 - r12x0;
    double r3y = r13y0 - r12y0;
    double v3x = v3x0 - v2x0;
    double v3y = v3y0 - v2y0;

    // Расстояние до спутника
    double r3 = std::sqrt(r3x * r3x + r3y * r3y);

    // Рассчитываем начальное положение ракеты
    double v0 = std::sqrt(Constants::G * Constants::M2 / (Constants::R2 + Constants::H));
    double rx0 = (Constants::R2 + Constants::H) * (r3x * cos(angle) - r3y * sin(angle)) / r3;
    double ry0 = (Constants::R2 + Constants::H) * (r3x * sin(angle) + r3y * cos(angle)) / r3;
    
    // Расстояние от ракеты до центра
    double r0 = sqrt(rx0 * rx0 + ry0 * ry0);
    
    // Начальная скорость ракеты (перпендикулярная радиус-вектору)
    double vx0 = -v0 * ry0 / r0;
    double vy0 = v0 * rx0 / r0;

    // Планета в центре координат
    Physics::r12x = 0;
    Physics::r12y = 0;



    std::ofstream fout_main("../labs/lab1/misha/result/task2/full_trajectory.csv");
    fout_main << "x y x3 y3\n";

     state_type y = {
        rx0, ry0,  // координаты ракеты
        r3x, r3y,  // координаты спутника
        vx0, vy0,  // скорость ракеты
        v3x, v3y   // скорость спутника
    };

    double t = 0.0;
    double t_circle_end = 60. * 60 * 10;
    double t_end = t_circle_end;
    double h = 0.1;

    double t_curr = t;
    runge_kutta_cash_karp54<state_type> stepper;

    integrate_adaptive(stepper, Physics::calculateForces, y, t, t_end, h,
                       [&](const state_type &state, double t) {
                           fout_main << y[0] << " " << y[1] << " " << y[2] << " " << y[3] << endl;
                       });
    fout_main.close();

    return 0;
}
