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
        double rx = y[0], ry = y[1], vx = y[4], vy = y[5]; // ракета
        double r13x = y[2], r13y = y[3], v3x = y[6], v3y = y[7]; // спутник

        // Проверка столкновения ракеты со спутником
        // if (sqrt((rx - r3x) * (rx - r3x) + (ry - r3y) * (ry - r3y)) <= Constants::R3) {
        //     throw std::runtime_error("Ракета попала в спутник.");
        // }

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
            Constants::G * (-Constants::M1 * rx / std::pow(r, 3) - Constants::M2 * (rx - r12x) / std::pow(r2, 3) -
                            Constants::M3 * (rx - r13x) / std::pow(r3, 3)),
            Constants::G * (-Constants::M1 * ry / std::pow(r, 3) - Constants::M2 * (ry - r12y) / std::pow(r2, 3) -
                            Constants::M3 * (ry - r13y) / std::pow(r3, 3))
        };

        // Ускорения для ракеты
        f[4] = f_rocket.first + f_gravity.first;
        f[5] = f_rocket.second + f_gravity.second;


        // Ускорения для спутника
        f[6] = -Constants::G * Constants::M1 * r13x / std::pow(r13, 3) -
               Constants::G * Constants::M2 * (r13x - r12x) / std::pow(r23, 3); // сила от планеты и солнца
        f[7] = -Constants::G * Constants::M1 * r13y / std::pow(r13, 3) -
               Constants::G * Constants::M2 * (r13y - r12y) / std::pow(r23, 3); // сила от планеты и солнца
    }
};

double Physics::mt;
double Physics::r12x;
double Physics::r12y;


int init(double mt, double angle, bool output = false) {
    ofstream fout_main;

    if (output) {
        std::ostringstream oss;
        oss << std::fixed << std::setprecision(2) << mt; // Ограничиваем 2 знаками после запятой
        std::string mt_str = oss.str();

        std::ostringstream oss_angle;
        oss_angle << std::fixed << std::setprecision(2) << angle;
        std::string angle_str = oss_angle.str();

        fout_main.open("../labs/lab1/misha/result/task2/full_trajectory_mt-" + mt_str + "_angle-" + angle_str + ".csv");
        fout_main << "x y x3 y3\n";
    }

    double r12x0 = -50368219856.43, r12y0 = -219503615669.65,
            v2x0 = 23688.42, v2y0 = -5739.71,
            r13x0 = -50370393498.29, r13y0 = -219513089385.80,
            v3x0 = 25513.15, v3y0 = -6666.22; // нач координаты планеты и спутника

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

    rx0 += r12x0;
    ry0 += r12y0;


    std::ofstream fout_main("../labs/lab1/misha/result/task2/full_trajectory.csv");
    fout_main << "x y x3 y3\n";

    state_type y = {
        rx0, ry0, // координаты ракеты
        r13x0, r13y0, // координаты спутника
        vx0, vy0, // скорость ракеты
        v3x, v3y // скорость спутника
    };

    double t = 0.0;
    double t_circle_end = 60. * 60 * 10;
    double t_end = t_circle_end;
    double h = 0.1;

    double t_curr = t;
    runge_kutta_cash_karp54<state_type> stepper;

    int status = 0; //0 - нет столкновения в пределах орбиты спутника,1 - cтолкновение со спутником

    while (t < t_end) {
        stepper.do_step(Physics::calculateForces, y, t, h);

        if (output) {
            fout_main << y[0] << " " << y[1] << " " << y[2] << " " << y[3] << std::endl;
        }

        // Проверка на столкновение ракеты и спутника
        double dx_rocket_sat = y[0] - y[2];
        double dy_rocket_sat = y[1] - y[3];
        double r_rocket_sat = std::sqrt(dx_rocket_sat * dx_rocket_sat + dy_rocket_sat * dy_rocket_sat);
        if (r_rocket_sat <= Constants::R3) {
            status = 1;
            break;
        }
        t += h;
    }

    if (output) {
        fout_main.close();
    }
}