#include <iostream>
#include <cmath>
#include <boost/numeric/odeint.hpp>
#include <fstream>
#include <chrono>

using namespace boost::numeric::odeint;
using state_type = std::array<double, 12>;
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
    // static double r12x; // координаты планеты
    // static double r12y;


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
        double rx = y[0], ry = y[1], vx = y[6], vy = y[7]; // ракета
        double r13x = y[2], r13y = y[3], v3x = y[8], v3y = y[9]; // cпутник
        double r12x = y[4], r12y = y[5], v2x = y[10], v2y = y[11];

        if (sqrt((rx - r13x) * (rx - r13x) + (ry - r13y) * (ry - r13y)) <= Constants::R3) {
            throw std::runtime_error("Ракета попала в планету.");
        }
        // Расчёт расстояний
        double r = std::sqrt(rx * rx + ry * ry); // расстояние от солнца до ракеты
        double r2 = distance(rx, ry, r12x, r12y); // расстояние от ракеты до планеты
        double r3 = distance(rx, ry, r13x, r13y); // расстояние от ракеты до спутника
        double r12= std::sqrt(r12x*r12x+r12y*r12y); // расстояние от планеты до солнца
        double r13 = std::sqrt(r13x * r13x + r13y * r13y); // расстояние от спутника до солнца
        double r23 = distance(r13x, r13y, r12x, r12y); // расстояние от спутника до планеты

        // Расчет скорости
        double v = std::sqrt(vx * vx + vy * vy); // скорость ракеты


        // Скорости
        f[0] = vx;
        f[1] = vy; // ракета
        f[2] = v3x;
        f[3] = v3y; // спутник
        f[4] = v2x;
        f[5] = v2y; // планета

        // Ускорения для ракеты
        f[6] = -(Constants::U * dm(t)* vx )/(v*m(t)) +
                    Constants::G * (
                   -Constants::M1 * rx / std::pow(r, 3)
                   - Constants::M2 * (rx - r12x) / std::pow(r2, 3)
                   - Constants::M3 * (rx - r13x) / std::pow(r3, 3)
               );

        f[7] = - (Constants::U * dm(t)* vy) / (v*m(t)) +
                Constants::G * (
                   -Constants::M1 * ry / std::pow(r, 3)
                   - Constants::M2 * (ry - r12y) / std::pow(r2, 3)
                   - Constants::M3 * (ry - r13y) / std::pow(r3, 3)
               );

        // Ускорения для спутника
        f[8] = -Constants::G * Constants::M1 * r13x / std::pow(r13, 3) -
               Constants::G * Constants::M2 * (r13x - r12x) / std::pow(r23, 3);

        f[9] = -Constants::G * Constants::M1 * r13y / std::pow(r13, 3) -
               Constants::G * Constants::M2 * (r13y - r12y) / std::pow(r23, 3);

        // Ускорение для планеты
        f[10] = -Constants::G * Constants::M1 *  r12x/ std::pow(r12, 3) +
                Constants::G * Constants::M3 * (r12x - r13x) / std::pow(r23, 3);

        f[11]=-Constants::G * Constants::M1 *  r12y/ std::pow(r12, 3) +
                Constants::G * Constants::M3 * (r12y - r13y) / std::pow(r23, 3);
    }
};

double Physics::mt;
// double Physics::r12x;
// double Physics::r12y;

int main() {
    double r12x0 = -50368219856.43, r12y0 = -219503615669.65,
            v2x0=23688.42, v2y0=-5739.71,
            r13x0 = -50370393498.29, r13y0 = -219513089385.80,
            v3x0 = 25513.15, v3y0 = -6666.22; // нач координаты планеты и спутника
    Physics::mt = 5000;
    double angle = 180. * M_PI / 180;

    double rx0, ry0, vx0, vy0;

    // Расчет нач значений
    // Physics::r12x = r12x0;
    // Physics::r12y = r12y0;
    double r3x = r13x0 - r12x0, r3y = r13y0 - r12y0; // координаты ракеты относительно планеты
    double r3 = std::sqrt(r3x * r3x + r3y * r3y); // расстояние между планетой и центром

    double v0 = 1.0 * std::sqrt(Constants::G * Constants::M2 / (Constants::R2 + Constants::H));
    rx0 = (Constants::R2 + Constants::H) * (r3x * cos(angle) - r3y * sin(angle)) / r3;
    ry0 = (Constants::R2 + Constants::H) * (r3x * sin(angle) + r3y * cos(angle)) / r3;
    double r0 = sqrt(rx0 * rx0 + ry0 * ry0);
    vx0 = -v0 * ry0 / r0;
    vy0 = v0 * rx0 / r0;

    vx0+=v2x0;
    vy0+=v2y0;

    rx0 += r12x0;
    ry0 += r12y0;


    std::ofstream fout_main("../labs/lab1/misha/res_task2/full_trajectory.csv");
    fout_main << "x y x3 y3 x2 y2\n";

    state_type y = {
        rx0, ry0, r13x0, r13y0,r12x0, r12y0,
        vx0, vy0, v3x0, v3y0,v2x0,v2y0 
    };

    double t = 0.0;
    double t_circle_end = 60. * 60 * 24 * 800;
    double t_end = t_circle_end;
    double h = 1000;

    double t_curr = t;
    runge_kutta_dopri5<state_type> stepper;

    integrate_adaptive(stepper, Physics::calculateForces, y, t, t_end, h,
                       [&](const state_type &state, double t) {
                           fout_main << y[0] << " " << y[1] << " " << y[2] << " " << y[3] << " " << y[4] << " " << y[5] <<  endl;
                       });
    fout_main.close();
    // try {
    //     while (t_curr < t_end) {
    //         stepper.do_step(Physics::calculateForces, y, t, h);
    //
    //         t_curr += h;
    //     }
    //     fout_main.close();
    // } catch (exception &e) {
    //     std::cerr << e.what() << std::endl;
    // }

    return 0;
}
