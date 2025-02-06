#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <windows.h>

double F(double q_last, double q, double q_next, double m, double alpha, double beta) {

    return 1.0 / (m) * ((q_next - 2 * q + q_last) + alpha * (q_next - q) * (q_next - q) +
                        beta * (q_next - q) * (q_next - q) * (q_next - q)
                        - alpha * (q - q_last) * (q - q_last) - beta * (q - q_last) * (q - q_last) * (q - q_last));
}

double Vf(double alpha, double beta, double r) {
    return r * r / 2 + alpha * r * r * r / 3 + beta * r * r * r * r / 4;
}

double H(double *q, double *v, int size, double m, double alpha, double beta) {
    double summ = 0;
    for (int i = 0; i < size - 1; i++) {
        summ += (m * v[i] * v[i] / 2 + Vf(alpha, beta, q[i + 1] - q[i]));
    }
    summ += (m * v[size - 1] * v[size - 1] / 2 + Vf(alpha, beta, (q[0] - q[size - 1])));
    return summ;
}

void SpeedVerle(double *q, double *v, int size, double alpha, double beta, double tau, int N, double m) {
    double l = 1;
    double *a = (double *) malloc(size * sizeof(double));
    int s = size - 1;
    FILE *f = fopen("../labs/lab3/misha/result/verle.txt", "w");

    for (int t = 0; t < N; t++) {

        for (int i = 0; i < s + 1; i++) {
            q[i] = q[i] + v[i] * tau + 0.5 * a[i] * tau * tau;
            v[i] = v[i] + 0.5 * a[i] * tau;
        }
        a[0] = F(l * q[s], l * q[0], l * q[1], m, alpha, beta) / l;
        v[0] = v[0] + 0.5 * a[0] * tau;
        for (int i = 1; i < s; i++) {
            a[i] = F(q[i - 1], q[i], q[i + 1], m, alpha, beta);
            v[i] = v[i] + 0.5 * a[i] * tau;

        }
        a[s] = F(q[s - 1], q[s], q[0], m, alpha, beta) / l;
        v[s] = v[s] + 0.5 * a[s] * tau;
        if (t % 1000 == 0) {
            for (int i = 0; i < size; i++) {
                fprintf(f, "%f ", v[i]);
            }
            fprintf(f, "\n");
        }
    }
    fclose(f);
    free(a);
}

double ksi = 0.1931833275037836;

void SimplexVerle(double *q, double *v, int size, double alpha, double beta, double tau, int N, double m) {
    double *a = (double *) malloc(size * sizeof(double));
    int s = size - 1;
    FILE *f = fopen("../labs/lab3/misha/result/simplex.txt", "w");

    for (int t = 0; t < N; t++) {

        for (int i = 0; i < s + 1; i++) {
            q[i] = q[i] + v[i] * tau * ksi;
        }
        a[0] = F(q[s], q[0], q[1], m, alpha, beta);
        for (int i = 1; i < s; i++) {
            a[i] = F(q[i - 1], q[i], q[i + 1], m, alpha, beta);
        }
        a[s] = F(q[s - 1], q[s], q[0], m, alpha, beta);
        for (int i = 0; i < s + 1; i++) {
            v[i] = v[i] + 0.5 * a[i] * tau;
            q[i] = q[i] + v[i] * tau * (1 - 2 * ksi);
        }
        a[0] = F(q[s], q[0], q[1], m, alpha, beta);
        for (int i = 1; i < s; i++) {
            a[i] = F(q[i - 1], q[i], q[i + 1], m, alpha, beta);
        }
        a[s] = F(q[s - 1], q[s], q[0], m, alpha, beta);

        for (int i = 0; i < s + 1; i++) {
            v[i] = v[i] + 0.5 * a[i] * tau;
            q[i] = q[i] + v[i] * tau * ksi;
        }

//        for (int i = s / 2; i >= 0; i--) {
//            v[i] = -v[s - i];
//            q[i] = -q[s - i];
//        }

        if (t % 1000 == 0) {
            for (int i = 0; i < size; i++) {
                fprintf(f, "%f ", v[i]);
            }
            fprintf(f, "\n");
        }
    }
    fclose(f);
    free(a);
}

int main() {
    SetConsoleOutputCP(CP_UTF8);
    int N = 1000;
    int size = N;
    double *q = (double *) malloc(size * sizeof(double));
    double *v = (double *) malloc(size * sizeof(double));
    double m = 1;
    double alpha = 0.0;
    double beta = 100.0;
//    printf("Введите значение для alpha: ");
//    scanf("%lf", &alpha);

    printf("Введите значение для beta: ");
    scanf("%lf", &beta);
    for (int i = 0; i < N; i++) {
        q[i] = 0.0;
        v[i] = 0.0;
    }

    q[N / 2 - 1] = -0.5;
    q[N / 2] = 0.5;

    double H0 = H(q, v, size, m, alpha, beta);
    printf("H0: %f\n", H0);

    SimplexVerle(q, v, size, alpha, beta, 0.01, 5e5, m);

    printf("Final H: %f\n", H(q, v, size, m, alpha, beta));

    free(q);
    free(v);
    return 0;
}