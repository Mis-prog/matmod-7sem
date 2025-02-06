

double F(double q_last, double q, double q_next, double m, double alpha, double beta);

double Vf(double alpha, double beta, double r);

double H(double *q, double *v, int size, double m, double alpha, double beta);

void SpeedVerle(double *q, double *v, int size, double alpha, double beta, double tau, int N, double m);

void SimplexVerle(double *q, double *v, int size, double alpha, double beta, double tau, int N, double m);

void SimplexVerleNew(double *q, double *v, double *a, int size, double alpha, double beta, double tau, double m);