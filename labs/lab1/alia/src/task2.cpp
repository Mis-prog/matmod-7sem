#include "task2_openmp.h"
// #include "task2.h"

int main() {
    // double start_time = omp_get_wtime();


    // analyzeAllAngles(0.0, 1.0 , 1.0);
    
    // double end_time = omp_get_wtime();
    // std::cout << "Общее время выполнения: " << (end_time - start_time) << " секунд" << std::endl;
    // double angle,mt;
    // cout << "Введите угол и массу топлива: \n";  
    // cin >> angle >> mt;
    // init(284,0,true);

    double mt=findMinFuel(0);
    cout << "Итоговое кол-во топливо: " <<  mt;
    return 0;
}
