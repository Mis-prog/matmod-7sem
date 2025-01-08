#include "task2_openmp.h"

int main() {
    double start_time = omp_get_wtime();


    analyzeAllAngles(0.0, 1.0 , 1.0);
    
    double end_time = omp_get_wtime();
    std::cout << "Общее время выполнения: " << (end_time - start_time) << " секунд" << std::endl;
    return 0;
}
