#include "task2.h"

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
    double t_circle_end = 60. * 60 * 24 * 28 * 5;
    double t_end = t_circle_end;
    double h = 7;

    double t_curr = t;
    runge_kutta_cash_karp54<state_type> stepper;

    integrate_adaptive(stepper, Physics::calculateForces, y, t, t_end, h,
                       [&](const state_type &state, double t) {
                           fout_main << y[0] << " " << y[1] << " " << y[2] << " " << y[3] << endl;
                            
                            // Проверка на столкновение ракеты и планеты
                            double r_rocket_planet = std::sqrt(y[0] * y[0] + y[1] * y[1]);
                            if (r_rocket_planet <= Constants::R2) {
                                std::cout << "Столкновение ракеты с планетой на времени t = " << t << " секунд." << std::endl;
                                fout_main.close();
                                exit(0); 
                            }

                            // Проверка на столкновение ракеты и спутника
                            double dx_rocket_sat = y[0] - y[2];
                            double dy_rocket_sat = y[1] - y[3];
                            double r_rocket_sat = std::sqrt(dx_rocket_sat * dx_rocket_sat + dy_rocket_sat * dy_rocket_sat);
                            if (r_rocket_sat <= Constants::R3) {
                                std::cout << "Столкновение ракеты со спутником на времени t = " << t << " секунд." << std::endl;
                                fout_main.close();
                                exit(0);
                            }
                       });
    fout_main.close();

    return 0;
}
