#include <iostream>
#include <cmath>
#include <boost/numeric/odeint.hpp>
#include <fstream>
#include <chrono>
#include <exception>
#include <string>
#include "task2_const.h"

using namespace boost::numeric::odeint;
using state_type = std::array<double, 8>;
using namespace std;

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

    // // Проверка столкновения ракеты со спутником
    // if (sqrt((rx - r13x) * (rx - r13x) + (ry - r13y) * (ry - r13y)) <= Constants::R3) {
    //     cout << "Ракета попала в спутник\n";
    //     return ;
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


int init(double mt,double angle,bool output=false){
    ofstream fout_main;

    if (output){
        std::ostringstream oss;
        oss << std::fixed << std::setprecision(2) << mt; // Ограничиваем 2 знаками после запятой
        std::string mt_str = oss.str();

        std::ostringstream oss_angle;
        oss_angle << std::fixed << std::setprecision(2) << angle;
        std::string angle_str = oss_angle.str();

        fout_main.open("../labs/lab1/alia/result/task2/full_trajectory_mt-" + mt_str + "_angle-" + angle_str + ".csv");
        fout_main << "x y x3 y3\n";
    }

    double r12x0 = -63051822463.82, r12y0 =141005424242.41,
            v2x0= -26912.62, v2y0= -11515.37,
            r13x0 = -63198811221.89, r13y0 = 141334141912.81,
            v3x0 =  -27897.52, v3y0 = -11974.71; // нач координаты планеты и спутника

    double _angle = angle * M_PI / 180;

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
    double rx0 = (Constants::R2 + Constants::H) * (r3x * cos(_angle) - r3y * sin(_angle)) / r3;
    double ry0 = (Constants::R2 + Constants::H) * (r3x * sin(_angle) + r3y * cos(_angle)) / r3;
    
    // Расстояние от ракеты до центра
    double r0 = sqrt(rx0 * rx0 + ry0 * ry0);
    
    // Начальная скорость ракеты (перпендикулярная радиус-вектору)
    double vx0 = -v0 * ry0 / r0;
    double vy0 = v0 * rx0 / r0;

    rx0+=r12x0;
    ry0+=r12y0;

    double t = 0.0;
    double t_circle_end = 60. * 60 * 24 * 28 * 10;
    double t_end = t_circle_end;
    double h = 7;

    state_type y = {
        rx0, ry0,  // координаты ракеты
        r13x0, r13y0,  // координаты спутника
        vx0, vy0,  // скорость ракеты
        v3x, v3y   // скорость спутника
    };
    
    double t_curr = t;
    runge_kutta_cash_karp54<state_type> stepper;


    int status=0; //0 - нет столкновения,1 - cтолкновение со спутником, 2 - столкновение с планетой

    while (t < t_end) {
        stepper.do_step(Physics::calculateForces, y, t, h);

        if (output){
            fout_main << y[0] << " " << y[1] << " " << y[2] << " " << y[3] << std::endl;
        }

        // Проверка на столкновение ракеты и планеты
        double r_rocket_planet = std::sqrt(y[0] * y[0] + y[1] * y[1]);
        if (r_rocket_planet <= Constants::R2) {
            // std::cout << "Столкновение ракеты с планетой на времени t = " << t << " секунд." << std::endl;
            status = 2;
            break;
        }

        // Проверка на столкновение ракеты и спутника
        double dx_rocket_sat = y[0] - y[2];
        double dy_rocket_sat = y[1] - y[3];
        double r_rocket_sat = std::sqrt(dx_rocket_sat * dx_rocket_sat + dy_rocket_sat * dy_rocket_sat);
        if (r_rocket_sat <= Constants::R3) {
            // std::cout << "Столкновение ракеты со спутником на времени t = " << t << " секунд." << std::endl;
            status = 1;
            break;
        }
        t+=h;
    }
    if (output) {
        fout_main.close();
    }

    return status;
}