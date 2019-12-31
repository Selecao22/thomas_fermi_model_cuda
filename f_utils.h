#ifndef _F_UTILS_H
#define _F_UTILS_H

#include <stdlib.h>
#include <math.h>

// Physical constants and some thrash like a const iterations
//const int N_X = 20;
#define POINT_NUMBER 54
#define PI 3.141592653589793
#define K 36.75
#define A 107.0 // tmp constant, must be as parameter
#define Z 47.0 // tmp constant, must be as parameter

double*
create_physic_array(int size);

void
get_init_assumption
(
        double *x,
        double *fi,
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


double
rect(
        const double* x,
        const double* f,
        int n
);


double
fint_neg12(double x);


double
fint_neg12_der(double x);


double
fint_12(double x);


double
fint_32(double x);


double Y(double x);


#endif