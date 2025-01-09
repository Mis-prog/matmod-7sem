#pragma once
#include "task2.h"
#include "omp.h"


struct Result {
    double angle;
    double min_fuel;
};


double findMinFuel(double angle, double epsilon = 0.05) {
    double left = 275;
    double right = 290;

    double success = 0;

    while (right - left > epsilon) {
        double mid = (left + right) / 2.0;
        int result = init(mid, angle, true);

        if (result == 1) {
            success = mid;
            right = mid;
        } else if (result == 3) {
            right = mid;
        } else if (result == 0) {
            left = mid;
        }
    }
    cout << "Оптимальное количество топлива: " << right << " для угла: " << angle << endl;
    return success;
}

void analyzeAllAngles(double start_angle = 0.0, double end_angle = 360.0, double angle_step = 1.0) {
    int num_angles = static_cast<int>((end_angle - start_angle) / angle_step) + 1;
    std::vector<Result> results(num_angles);

    int num_threads = omp_get_max_threads();
    std::cout << "Запуск вычислений на " << num_threads << " потоках" << std::endl;

#pragma omp parallel for schedule(dynamic) num_threads(num_threads)
    for (int i = 0; i < num_angles; i++) {
        double angle = start_angle + i * angle_step;
        double min_fuel = findMinFuel(angle);

        results[i] = {angle, min_fuel};

#pragma omp critical
        {
            std::cout << "Поток " << omp_get_thread_num()
                    << " | Угол: " << std::fixed << std::setprecision(2) << angle
                    << "° | Минимальное топливо: " << min_fuel << std::endl;
        }
    }

    // Сохраняем результаты в файл
    std::ofstream output("../labs/lab1/alia/result/min_fuel_results.csv");
    output << "angle,min_fuel\n";

    for (const auto &result: results) {
        output << std::fixed << std::setprecision(2)
                << result.angle << "," << result.min_fuel << "\n";
    }

    output.close();
}
