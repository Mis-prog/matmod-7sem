#include <iostream>
#include <stdio.h>
#include <vector>
#include <iomanip>
#include <fstream>
#include <time.h>
#include <math.h>

double f(double t, double x, double y, double z, double a, double b) {
    return a * (y - x);
}

double g(double t, double x, double y, double z, double a, double b) {
    return b * y  - x * z;
}

double h(double t, double x, double y, double z, double a, double b) {
    return -3 * z + x * y;
}

std::vector<double> new_point(double f(double, double, double, double, double, double),
                            double g(double, double, double, double, double, double),
                            double h(double, double, double, double, double, double),
                            double t, double x, double y, double z, 
                            double step, double a, double b)
{
    double kx0, ky0, kz0, kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3;

    kx0 = f(t, x, y ,z, a, b);
    ky0 = g(t, x, y, z, a, b);
    kz0 = h(t, x, y, z, a, b);

    kx1 = f(t + (step / 2), x + (step / 2) * kx0, y + (step / 2) * ky0, z + (step / 2) * kz0, a, b);
    ky1 = g(t + (step / 2), x + (step / 2) * kx0, y + (step / 2) * ky0, z + (step / 2) * kz0, a, b);
    kz1 = h(t + (step / 2), x + (step / 2) * kx0, y + (step / 2) * ky0, z + (step / 2) * kz0, a, b);

    kx2 = f(t + (step / 2), x + (step / 2) * kx1, y + (step / 2) * ky1, z + (step / 2) * kz1, a, b);
    ky2 = g(t + (step / 2), x + (step / 2) * kx1, y + (step / 2) * ky1, z + (step / 2) * kz1, a, b);
    kz2 = h(t + (step / 2), x + (step / 2) * kx1, y + (step / 2) * ky1, z + (step / 2) * kz1, a, b);

    kx3 = f(t + step, x + step * kx2, y + step * ky2, z + step * kz2, a, b);
    ky3 = g(t + step, x + step * kx2, y + step * ky2, z + step * kz2, a, b);
    kz3 = h(t + step, x + step * kx2, y + step * ky2, z + step * kz2, a, b);

    x = x + (step / 6) * (kx0 + 2 * kx1 + 2 * kx2 + kx3);
    y = y + (step / 6) * (ky0 + 2 * ky1 + 2 * ky2 + ky3);
    z = z + (step / 6) * (kz0 + 2 * kz1 + 2 * kz2 + kz3);

    return std::vector<double>{x, y, z};
}


void write_file(double a, double b,
                std::vector<double> t,
                std::vector<double> x, 
                std::vector<double> y, 
                std::vector<double> z)
{
    std::ofstream fout;
    std::ofstream fout_stat;
    fout.open("../labs/lab4/misha/result/result.txt");
    fout_stat.open("../labs/lab4/misha/result/stat_points.txt");

    if (t.size() != x.size() || t.size() != y.size() || t.size() != z.size()) {
        std::cout << "Dimension t is not equal with other dimensions" << std::endl;
        fout.close();
        throw "dimensions error";
    }

    fout_stat << a << " " << b << std::endl;
    fout_stat << sqrt(3*b) << " " << sqrt(3*b) << " " << b << std::endl;
    fout_stat << -sqrt(3*b) << " " << -sqrt(3*b) << " " << b << std::endl;
    fout_stat << 0 << " " << 0 << " " << 0 << std::endl;
    fout_stat.close();

    for (int i = 0; i < x.size(); i++) {
        fout << t[i] << " " << x[i] << " " << y[i] << " " << z[i] << "\n";
    }

    fout.close();
}


int main()
{
    double step = 0.01;
    int n = 100000;
    double t0 = 0, t1 = t0 + n * step;
    std::cout << "Введите  a и b: " << std::endl;
    double a = - 0.01 , b = 2.9; 
    std::cin >> a >> b; 
    double x0 = 1, y0 = 1, z0 = 1;

    std::vector<double> t(n), x(n), y(n), z(n);
    std::vector<double> point(3);
    double t_;

    x[0] = x0; y[0] = y0; z[0] = z0;
    for (int i = 0; i < n - 1; i++)
    {
        t_ = i * step;
        point = new_point(f, g, h, t_, x[i], y[i], z[i], step, a, b);
        t[i + 1] = t_;
        x[i + 1] = point[0];
        y[i + 1] = point[1];
        z[i + 1] = point[2];
    }

    write_file(a, b, t, x, y, z);

    // system("python ../labs/lab4/misha/src/draw.py");
    
    return 0;
}
