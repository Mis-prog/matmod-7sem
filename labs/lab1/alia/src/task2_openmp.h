#pragma once
#include "task2.h"
#include "omp.h"


struct Result {
    double angle;
    double min_fuel;
};


double findMinFuel(double angle, double epsilon = 0.05) {
    double left = 272.6;
    double right = 290;

    double success = 0;

    double mt = left;
    while (true) {
        success = init(mt, angle, false);
        if (success == 1 && mt < right) {
            // cout << "Оптимальное количество топлива: " << mt << " для угла: " << angle << endl;
            break;
        }
        mt += 0.1;
    }
    if (mt > right) {
        return -1;
    }
    return mt;
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

        if (min_fuel == -1) {
            continue; // Пропускаем неудачные расчеты
        }
        results[i] = {angle, min_fuel};

#pragma omp critical
        {
            std::cout << "Поток " << omp_get_thread_num()
                    << " | Угол: " << std::fixed << std::setprecision(2) << angle
                    << "° | Минимальное топливо: " << min_fuel << std::endl;
        }
    }

    // Сохраняем результаты в файл
    std::ofstream output("../labs/lab1/alia/result/t/min_fuel_results.csv");
    output << "angle,min_fuel\n";

    for (const auto &result: results) {
        output << std::fixed << std::setprecision(2)
                << result.angle << "," << result.min_fuel << "\n";
    }

    output.close();
}
