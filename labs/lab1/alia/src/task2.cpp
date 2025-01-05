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

    static constexpr double T = 2400.0; // время работы двигателя (с)
    static constexpr double H = 300e3; // высота орбиты (м)
    static constexpr double M0 = 120.0; // масса полезной назрузки (кг)
    static constexpr double U = 3060.0; // скорость истечения (м/c)
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
    double r13x = y[2], r13y = y[3], v3x = y[6], v3y = y[7];  // спутник

    // Проверка столкновения ракеты со спутником
    if (sqrt((rx - r13x) * (rx - r13x) + (ry - r13y) * (ry - r13y)) <= Constants::R3) {
        throw std::runtime_error("Ракета попала в спутник.");
    }

    // Расчёт расстояний
    double r = std::sqrt(rx * rx + ry * ry); // расстояние от солнца до ракеты
    double r2 = distance(rx, ry, r12x, r12y); // расстояние от ракеты до планеты
    double r3 = distance(rx, ry, r13x, r13y); // расстояние от ракеты до спутника
    double r13 = std::sqrt(r13x * r13x + r13y * r13y); // расстояние от спутника до солнца
    double r23 = distance(r13x, r13y, r12x, r12y); // расстояние от спутника до планеты

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
         Constants::G * (
            // -Constants::M1 * rx / std::pow(r, 3)
         - Constants::M2 * (rx - r12x) / std::pow(r2, 3)- Constants::M3 * (rx - r13x) / std::pow(r3, 3)),
        Constants::G * (
            // -Constants::M1 * ry / std::pow(r, 3)
        - Constants::M2 * (ry - r12y) / std::pow(r2, 3)- Constants::M3 * (ry - r13y) / std::pow(r3, 3))
    };

    // Ускорения для ракеты
    f[4] = f_gravity.first+f_rocket.first;
    f[5] = f_gravity.second+f_rocket.second;


    // Ускорения для спутника
    f[6] =
        //  -Constants::G * Constants::M1 * r13x / std::pow(r13, 3) 
         -  Constants::G * Constants::M2 * (r13x - r12x) / std::pow(r23, 3);  // сила от планеты и солнца
    f[7] =
            //  -Constants::G * Constants::M1 * r13y / std::pow(r13, 3) 
             -Constants::G * Constants::M2 * (r13y - r12y) / std::pow(r23, 3);  // сила от планеты и солнца
    }
};

double Physics::mt;
double Physics::r12x;
double Physics::r12y;

int main() {
    double r12x0 = -63051822463.82, r12y0 =141005424242.41,
            v2x0= -26912.62, v2y0= -11515.37,
            r13x0 = -63198811221.89, r13y0 = 141334141912.81,
            v3x0 =  -27897.52, v3y0 = -11974.71; // нач координаты планеты и спутника

    double angle_input;
    cout << "Введите общую массу и угол: \n";
    cin >> Physics::mt >> angle_input;
    double angle = angle_input * M_PI / 180;

    double v3x = v3x0 - v2x0;
    double v3y = v3y0 - v2y0;

    // Планета в центре координат
    Physics::r12x = r12x0;
    Physics::r12y = r12y0;

    // Расстояние до спутника
    double r3x = r13x0 - Physics::r12x, r3y = r13y0 - Physics::r12y;
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

    rx0+=r12x0;
    ry0+=r12y0;


    std::ofstream fout_main("../labs/lab1/alia/result/task2/full_trajectory.csv");
    fout_main << "x y x3 y3\n";

    state_type y = {
        rx0, ry0,  // координаты ракеты
        r13x0, r13y0,  // координаты спутника
        vx0, vy0,  // скорость ракеты
        v3x, v3y   // скорость спутника
    };

    double t = 0.0;
    double t_circle_end = 60. * 60 * 24 * 28;
    double t_end = t_circle_end;
    double h = 5;

    double t_curr = t;
    runge_kutta_cash_karp54<state_type> stepper;

    integrate_adaptive(stepper, Physics::calculateForces, y, t, t_end, h,
                       [&](const state_type &state, double t) {
                           fout_main << y[0] << " " << y[1] << " " << y[2] << " " << y[3] << endl;
                       });
    fout_main.close();

    return 0;
}
