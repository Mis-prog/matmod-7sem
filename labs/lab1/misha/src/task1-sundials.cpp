#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <omp.h>
#include <cvode/cvode.h>
#include <nvector/nvector_serial.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <sunmatrix/sunmatrix_dense.h>

struct Constants {
    static constexpr double G = 6.67e-11;
    static constexpr double M1 = 2.0e30;
    static constexpr double M2 = 6.4e23;
    static constexpr double M3 = 1.1e16;
    static constexpr double R1 = 696340e3;
    static constexpr double R2 = 3390e3;
    static constexpr double R3 = 11.1e3;
    static constexpr double R12 = 228e9;
    static constexpr double R23 = 9.4e6;
    static constexpr double U2 = 24e3;
    static constexpr double U3 = 2.14e3;
};

struct UserData {
    double G, M1, M2, M3;
};

static int orbitDerivatives(double t, N_Vector y, N_Vector ydot, void* user_data) {
    UserData* data = static_cast<UserData*>(user_data);

    #pragma omp parallel sections
    {
        #pragma omp section
        {
            // Расчет для планеты
            double x2 = NV_Ith_S(y, 0), vx2 = NV_Ith_S(y, 1);
            double y2 = NV_Ith_S(y, 2), vy2 = NV_Ith_S(y, 3);
            double x3 = NV_Ith_S(y, 4), r_12 = std::sqrt(x2 * x2 + y2 * y2);
            double r_23 = std::sqrt(std::pow(x3 - x2, 2) + std::pow(NV_Ith_S(y, 6) - y2, 2));

            NV_Ith_S(ydot, 0) = vx2;
            NV_Ith_S(ydot, 2) = vy2;
            NV_Ith_S(ydot, 1) = -data->G * data->M1 * x2 / std::pow(r_12, 3) +
                                data->G * data->M3 * (x3 - x2) / std::pow(r_23, 3);
            NV_Ith_S(ydot, 3) = -data->G * data->M1 * y2 / std::pow(r_12, 3) +
                                data->G * data->M3 * (NV_Ith_S(y, 6) - y2) / std::pow(r_23, 3);
        }

        #pragma omp section
        {
            // Расчет для астероида
            double x3 = NV_Ith_S(y, 4), vx3 = NV_Ith_S(y, 5);
            double y3 = NV_Ith_S(y, 6), vy3 = NV_Ith_S(y, 7);
            double x2 = NV_Ith_S(y, 0), r_13 = std::sqrt(x3 * x3 + y3 * y3);
            double r_23 = std::sqrt(std::pow(x3 - x2, 2) + std::pow(y3 - NV_Ith_S(y, 2), 2));

            NV_Ith_S(ydot, 4) = vx3;
            NV_Ith_S(ydot, 6) = vy3;
            NV_Ith_S(ydot, 5) = -data->G * data->M1 * x3 / std::pow(r_13, 3) +
                                data->G * data->M2 * (x2 - x3) / std::pow(r_23, 3);
            NV_Ith_S(ydot, 7) = -data->G * data->M1 * y3 / std::pow(r_13, 3) +
                                data->G * data->M2 * (NV_Ith_S(y, 2) - y3) / std::pow(r_23, 3);
        }
    }

    return 0;
}

int main() {
    const int N = 8;
    N_Vector y = N_VNew_Serial(N);

    // Точные начальные условия
    NV_Ith_S(y, 0) = Constants::R1 + Constants::R12;  // x2
    NV_Ith_S(y, 1) = Constants::U2;                   // vx2
    NV_Ith_S(y, 2) = 0.0;                             // y2
    NV_Ith_S(y, 3) = 0.0;                             // vy2

    NV_Ith_S(y, 4) = Constants::R1 + Constants::R12 + Constants::R23;  // x3
    NV_Ith_S(y, 5) = Constants::U3;                   // vx3
    NV_Ith_S(y, 6) = 0.0;                             // y3
    NV_Ith_S(y, 7) = 0.0;                             // vy3

    // Остальной код без изменений
    UserData data{Constants::G, Constants::M1, Constants::M2, Constants::M3};



    void* cvode_mem = CVodeCreate(CV_BDF);
    CVodeSetMaxErrTestFails(cvode_mem, 20);  // Увеличьте допустимое число провалов теста ошибок
    CVodeSetMinStep(cvode_mem, 1e-10);       // Очень маленький минимальный шаг
    CVodeSetMaxStep(cvode_mem, 1e5);
    CVodeInit(cvode_mem, orbitDerivatives, 0, y);
    CVodeSStolerances(cvode_mem, 1e-6, 1e-6);
    CVodeSetUserData(cvode_mem, &data);
    CVodeSetMaxNumSteps(cvode_mem, 100000);

    // Дополнительные настройки solver
    CVodeSetInitStep(cvode_mem, 100.0);
    CVodeSetMinStep(cvode_mem, 1.0);
    CVodeSetMaxStep(cvode_mem, 10000.0);

    // Остальной код без изменений


    SUNMatrix A = SUNDenseMatrix(N, N);
    SUNLinearSolver LS = SUNLinSol_Dense(y, A);
    CVodeSetLinearSolver(cvode_mem, LS, A);

    double t = 0.0, t_end = 60.0 * 60 * 24 * 365 * 2000;
    double h = 5000.0;

    std::ofstream fout_main("path_full.csv");
    fout_main << "x2 y2 x3 y3\n";

    while (t < t_end) {
        int flag = CVode(cvode_mem, t + h, y, &t, CV_NORMAL);

        if (flag < 0) break;

        fout_main << NV_Ith_S(y, 0) << " " << NV_Ith_S(y, 2) << " "
                  << NV_Ith_S(y, 4) << " " << NV_Ith_S(y, 6) << "\n";
    }

    // Освобождение ресурсов
    N_VDestroy(y);
    CVodeFree(&cvode_mem);
    SUNLinSolFree(LS);
    SUNMatDestroy(A);

    return 0;
}