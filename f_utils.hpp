#include <cstdlib>
#include <cuda_runtime.h>
#include <cmath>

// Physical constants and some thrash like a const iterations
const int N_X = 20;
const int N = 10000;
const int POINT_NUMBER = 54;
const double PI = 3.141592653589793;
const double K = 36.75;
const double H = 1.0/N;
const double A = 107.0; // tmp constant, must be as parameter
const double Z = 47.0; // tmp constant, must be as parameter

double*
create_physic_array(int size);

void
get_init_assumption
(
        double **x,
        double **fi,
        double z,
        double tet,
        double r0,
        double h,
        double nu_0,
        int n
);

void
calculate_entrope(
        double* x2int32,
        double* se_array,
        double* fi,
        double h,
        double n
);


__host__
__device__
double
rect(
        const double* x,
        const double* f,
        int n
);


__host__
__device__
double
fint_neg12(double x);


__host__
__device__
double
fint_neg12_der(double x);


__host__
__device__
double
fint_12(double x);


__host__
__device__
double
fint_32(double x);



__host__
__device__
double Y(double x);

