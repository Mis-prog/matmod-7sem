#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <cvode/cvode.h>
#include <nvector/nvector_serial.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <omp.h>
#include "progressbar.hpp"


struct Constants {
    static constexpr double G = 6.67e-11; // гравитационная постоянная
    static constexpr double M1 = 2.0e30; // масса звезды (кг)
    static constexpr double M2 = 6.4e23; // масса планеты (кг)
    static constexpr double M3 = 1.1e16; // масса астероида (кг)
    static constexpr double R1 = 696340e3; // радиус звезды (м)
    static constexpr double R2 = 3390e3; // радиус планеты (м)
    static constexpr double R3 = 11.1e3; // радиус астероида (м)
    static constexpr double R12 = 228e9; // начальное расстояние звезда-планета (м)
    static constexpr double R23 = 9.4e6; // начальное расстояние планета-астероид (м)
    static constexpr double U2 = 24e3; // начальная скорость планеты (м/с)
    static constexpr double U3 = 2.14e3; // начальная скорость астероида (м/с)
};


static int orbitDerivatives(double t, N_Vector y, N_Vector ydot, void *user_data) {
    double x2 = NV_Ith_S(y, 0), vx2 = NV_Ith_S(y, 1);
    double y2 = NV_Ith_S(y, 2), vy2 = NV_Ith_S(y, 3);
    double x3 = NV_Ith_S(y, 4), vx3 = NV_Ith_S(y, 5);
    double y3 = NV_Ith_S(y, 6), vy3 = NV_Ith_S(y, 7);

    double r_12 = std::sqrt(x2 * x2 + y2 * y2);
    double r_13 = std::sqrt(x3 * x3 + y3 * y3);
    double r_23 = std::sqrt(std::pow(x3 - x2, 2) + std::pow(y3 - y2, 2));


    NV_Ith_S(ydot, 0) = vx2;
    NV_Ith_S(ydot, 2) = vy2;

    NV_Ith_S(ydot, 4) = vx3;
    NV_Ith_S(ydot, 6) = vy3;

    NV_Ith_S(ydot, 1) = -Constants::G * Constants::M1 * x2 / std::pow(r_12, 3) +
                        Constants::G * Constants::M3 * (x3 - x2) / std::pow(r_23, 3);


    NV_Ith_S(ydot, 3) = -Constants::G * Constants::M1 * y2 / std::pow(r_12, 3) +
                        Constants::G * Constants::M3 * (y3 - y2) / std::pow(r_23, 3);


    NV_Ith_S(ydot, 5) = -Constants::G * Constants::M1 * x3 / std::pow(r_13, 3) +
                        Constants::G * Constants::M2 * (x2 - x3) / std::pow(r_23, 3);


    NV_Ith_S(ydot, 7) = -Constants::G * Constants::M1 * y3 / std::pow(r_13, 3) +
                        Constants::G * Constants::M2 * (y2 - y3) / std::pow(r_23, 3);
    return 0;
}

int main() {
    const int N = 8; // Число переменных
    N_Vector y = N_VNew_Serial(N);

    // Начальные условия
    NV_Ith_S(y, 0) = Constants::R1 + Constants::R12 + Constants::R2;
    NV_Ith_S(y, 1) = 0;
    NV_Ith_S(y, 2) = 0;
    NV_Ith_S(y, 3) = Constants::U2; // планета (x,vx,y,vy)
    NV_Ith_S(y, 4) = Constants::R1 + Constants::R12 + 2 * Constants::R2 + Constants::R23 + Constants::R3;
    NV_Ith_S(y, 5) = 0;
    NV_Ith_S(y, 6) = 0;
    NV_Ith_S(y, 7) = Constants::U3 + Constants::U2; // астероид (x,vx,y,vy)


    void *cvode_mem = CVodeCreate(CV_BDF);

    CVodeInit(cvode_mem, orbitDerivatives, 0, y);
    CVodeSStolerances(cvode_mem, 1e-4, 1e-4);

    SUNMatrix A = SUNDenseMatrix(N, N);
    SUNLinearSolver LS = SUNLinSol_Dense(y, A);
    CVodeSetLinearSolver(cvode_mem, LS, A);

    double t = 0.0;
    double t_end = 60.0 * 60 * 24 * 365 * 500;
    double h = 4000;

    std::ofstream fout_main("../../../../../labs/lab1/misha/res_task1/path_full.csv");
    fout_main << "x2 y2 x3 y3\n";

    while (t < t_end) {
        int flag = CVode(cvode_mem, t + h, y, &t, CV_NORMAL);

        if (flag < 0) {
            std::cerr << "Ошибка интегрирования" << std::endl;
            break;
        }

        fout_main << NV_Ith_S(y, 0) << " " << NV_Ith_S(y, 2) << " "
                << NV_Ith_S(y, 4) << " " << NV_Ith_S(y, 6) << "\n";
    }


    N_VDestroy(y);
    CVodeFree(&cvode_mem);
    SUNLinSolFree(LS);
    SUNMatDestroy(A);

    std::cout << "Расчет успешно завершен" << std::endl;
    return 0;
}
