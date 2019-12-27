#include <cuda_runtime.h>

// Physical constants and some thrash like a const iterations
const int N_X = 20;
const int N = 1000;
const int POINT_NUMBER = 54;
const double PI = 3.141592653589793;
const double K = 36.75;
const double H = 1.0/N;
const double A = 107.0;
const double Z = 47.0;

double* create_physic_array(int size);


