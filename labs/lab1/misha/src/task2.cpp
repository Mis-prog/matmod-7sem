#include "task2_openmp.h"

int main() {
    // init(2,270,true);

    // findMinFuel(0);

    double start_time = omp_get_wtime();

    analyzeAllAngles(0,360 , 2.0);

    double end_time = omp_get_wtime();
    std::cout << "Общее время выполнения: " << (end_time - start_time) << " секунд" << std::endl;
    return 0;
}
