#include <cstdlib>

// Physical constants and some thrash like a const iterations
const int N_X = 20;
const int N = 100000;
const int POINT_NUMBER = 54;
const double PI = 3.141592653589793;
const double K = 36.75;
const double H = 1.0/N;
const double A = 107.0; // tmp constant, must be as parameter
const double Z = 47.0; // tmp constant, must be as parameter

double* create_physic_array(int size);

void
get_init_assumption
(
        double **dfi,
        double **x,
        double **fi,
        double z,
        double tet,
        double r0,
        double h,
        double nu_0,
        int n
);


