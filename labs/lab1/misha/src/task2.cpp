#include "task2.h"


int main() {
    double r12x0 = -45779049843.29, r12y0 = -220565040943.48,
            v2x0 = 23796.90, v2y0 = -5243.20,
            r13x0 = -45779234115.16, r13y0 = -220575669506.13,
            v3x0 = 26035.80, v3y0 = -5315.35; // нач координаты планеты и спутника
    Physics::mt = 40;
    double angle = 270. * M_PI / 180;

    double rx0, ry0, vx0, vy0;

    // Расчет нач значений
    Physics::r12x = r12x0;
    Physics::r12y = r12y0;
    double r3x = r13x0 - Physics::r12x, r3y = r13y0 - Physics::r12y; // координаты ракеты относительно планеты
    double r3 = std::sqrt(r3x * r3x + r3y * r3y); // расстояние между планетой и центром

    double v0 = 1.0 * std::sqrt(Constants::G * Constants::M2 / (Constants::R2 + Constants::H));
    rx0 = (Constants::R2 + Constants::H) * (r3x * cos(angle) - r3y * sin(angle)) / r3;
    ry0 = (Constants::R2 + Constants::H) * (r3x * sin(angle) + r3y * cos(angle)) / r3;
    double r0 = sqrt(rx0 * rx0 + ry0 * ry0);
    vx0 = -v0 * ry0 / r0;
    vy0 = v0 * rx0 / r0;

    // vx0+=v2x0;
    // vy0+=v2y0;

    rx0 += Physics::r12x;
    ry0 += Physics::r12y;


    std::ofstream fout_main("../../../../../labs/lab1/misha/res_task2/full_trajectory.csv");
    fout_main << "x y x3 y3\n";

    state_type y = {
            rx0, ry0, r13x0, r13y0,
            vx0, vy0, v3x0 - v2x0, v3y0 - v2y0
    };

    double t = 0.0;
    double t_circle_end = 60. * 60 * 24;
    double t_end = t_circle_end;
    double h = 0.1;

    double t_curr = t;
    runge_kutta_dopri5<state_type> stepper;

    integrate_adaptive(stepper, Physics::calculateForces, y, t, t_end, h,
                       [&](const state_type &state, double t) {
                           fout_main << y[0] << " " << y[1] << " " << y[2] << " " << y[3] << endl;
                       });
    fout_main.close();

//    try {
//        while (t_curr < t_end) {
//            stepper.do_step(Physics::calculateForces, y, t, h);
//            fout_main << y[0] << " " << y[1] << " " << y[2] << " " << y[3] << endl;
//            t_curr += h;
//        }
//        fout_main.close();
//    } catch (exception &e) {
//        std::cerr << e.what() << std::endl;
//    }

    return 0;
}


//const int max_iterations = 100;
//const double tolerance = 1e-6;
//
//// Начальные точки для оптимизации
//std::vector<std::vector<double>> initial_points = {
//        {M_PI/2, 200.0},     // Первая точка
//        {M_PI/3, 250.0},     // Вторая точка
//        {2*M_PI/3, 150.0}    // Третья точка
//};
//
//// Создаем оптимизатор Нелдера-Мида
//OptimizationProblem problem;
//
//// Запускаем оптимизацию
//std::vector<double> result = {M_PI/2, 200.0};  // Начальное приближение
//double minimum = std::numeric_limits<double>::max();
//
//for (int i = 0; i < max_iterations; ++i) {
//std::vector<double> new_point = result;
//
//// Случайное возмущение текущей точки
//new_point[0] += (rand() / double(RAND_MAX) - 0.5) * M_PI / 10;
//new_point[1] += (rand() / double(RAND_MAX) - 0.5) * 50;
//
//double new_value = problem(new_point);
//if (new_value < minimum) {
//minimum = new_value;
//result = new_point;
//
//cout << "Iteration " << i << ": angle = " << result[0] * 180/M_PI
//<< " degrees, fuel = " << result[1]
//<< " kg, distance = " << minimum << " m" << endl;
//}
//}
//
//cout << "\nOptimal parameters found:" << endl;
//cout << "Angle: " << result[0] * 180/M_PI << " degrees" << endl;
//cout << "Fuel mass: " << result[1] << " kg" << endl;
//cout << "Minimum distance: " << minimum << " m" << endl;
//
//// Запуск финальной симуляции с оптимальными параметрами
//std::ofstream fout_main("optimal_trajectory.csv");
//fout_main << "x y x3 y3\n";
//
//// ... остальной код main() для сохранения траектории ...
//
//return 0;