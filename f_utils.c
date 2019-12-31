//
// Created by nikky on 26.12.2019.
//
#include "f_utils.h"

double* create_physic_array(int size)
{
    double *array = (double*)calloc(size, sizeof(double));
    double point_step = 0.1;

    array[0] = 0.1;

    for (int i = 1; i < POINT_NUMBER; ++i) {
        array[i] = array[i - 1] + point_step;
        if (((i % 10 == 9) && i < 11 ) || ((i % 9) == 0 && i > 11))
            point_step*=10;
    }

    return array;
}


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
        )
{
    for (int block = 0; block < n; ++block) {
        double hi = h * block;
        fi[block] = (z / tet / r0) * (1.0 - 1.5 * pow(hi, 2.0) + 0.5 * pow(hi, 6.0)) - nu_0 * pow(hi, 2.0);
        x[block] = hi;
    }
}


void calculate_entrope(
        double* x2int32,
        double* se_array,
        double* fi,
        double h,
        double n
        )
{
    for (int block = 0; block < n; ++block) {
        x2int32[block] = 2.0 * pow(block * h, 5.0) * fint_32(fi[block] / pow(block * h, 2.0));
        se_array[block] = 2.0 * pow(block * h, 5.0) * ((5.0 / 3.0) * fint_32(fi[block] / pow(block * h, 2.0)) -
                                                       (fi[block] / pow(block * h, 2.0)) * fint_12(fi[block] / pow(block * h, 2.0)));
    }
}


// IMPROVE - REDUCTION OPERATION

double
rect(
        const double* x,
        const double* f,
        int n
        )
{
    double s = 0.0;

    for (int i = 1; i < n; ++i)
        s = s + f[i] * (x[i] - x[i - 1]);

    return s;
}


double
fint_neg12(double x)
{
    if ( x >= 100.0)
    {
        return 2.0 * pow(x, 0.5);
    }

    if (x <= - 50.0)
    {
        return sqrt(PI) * exp(x);
    }

    double pi6_pow_1_3 = pow(PI / 6.0, 1.0 / 3.0);
    double exp_of_x = exp(2.0 * x / 3.0);

    // quesioning about my existance
    double res = 2.0 * 1.5 * pow(1.5, 0.5) *
            pow(log(1.0 + pi6_pow_1_3) * exp_of_x, 0.5) /
            (1.0 + pi6_pow_1_3 * exp_of_x) * (-pi6_pow_1_3) * exp_of_x * 1.5;

    return res;
}


double
fint_12(double x) {
    if (x >= 100.0) {
        return 1.5 * pow(x, 1.5);
    }

    if (x <= -50.0)
    {
        return 0.5 * sqrt(PI) * exp(x);
    }

    double pi6_pow_1_3 = pow(PI / 6.0, 1.0 / 3.0);
    double exp_of_x = exp(2.0 * x / 3.0);

    return pow(1.5, 0.5) * pow(log(pi6_pow_1_3 * exp_of_x), 1.5);
}


double
fint_32(double x)
{
    if (x >= 100.0)
    {
        return (2.0 / 5.0) * pow(x, 2.5);
    }

    if (x <= -50.0)
    {
        return 0.75 * sqrt(PI) * exp(x);
    }

    double pi6_pow_1_3 = pow(PI / 6.0, 1.0 / 3.0);
    double exp_of_x = exp(2.0 * x / 3.0);

    double f12 = pow(1.5, 0.5) * pow(log(1.0 + pi6_pow_1_3 * exp_of_x), 1.5);
    return 0.3 * f12 * pow(185.0 * f12 + 18.0 * pow(f12, 2.0), 1.0 / 3.0);
}



double
fint_neg12_der(double x)
{
    if (x >= 120.0)
        return 1 / sqrt(1.0*x) + 12.4298 / pow(x, 4.5+1.2337) / pow(x, 2.5);

    if (x <= -50.0)
        return sqrt(PI) * exp(x);

    return 2.0 * sqrt(3.0/2.0) * (-pow(2.0, (1.0 / 3.0)) * pow(PI, 2.0 / 3.0) * exp(4.0 * x / 3.0) *
        sqrt(log(pow(PI / 6.0, 1.0 / 3.0) * exp(2.0 * x /3.0) + 1.0))) / (3.0 * pow(3.0, 2.0 / 3.0) *
        pow(pow(PI/6.0, 1.0/3.0) * exp(2.0 * x / 3.0) + 1.0, 2.0) +
        (pow(PI, 2.0 / 3.0) * exp( 4.0 * x / 3.0) / ( 3.0 * pow(6.0, 2.0 / 3.0) * (pow(PI / 6.0, 1.0 / 3.0) *
        pow(exp(2.0 * x / 3.0) + 1.0, 2.0) * sqrt(log(pow(PI / 6.0, 1.0 / 3.0) * exp(2.0 * x / 3.0) + 1.0)))) +
        (pow(2.0,2.0 / 3.0) * pow(PI / 3.0, 1.0 / 3.0) * exp(2.0 * x / 3.0)*sqrt(log(pow(PI / 6.0,1.0 / 3.0) *
        exp((2.0 * x) / 3.0) + 1.0))) / (3.0 * (pow(PI / 6.0,1.0 / 3.0) * exp(2.0 * x / 3.0) + 1.0))));
}



double Y(double x)
{
    double low_bound = -10.0;
    const int N_grid = 1000;
    double h;
    double x_array[N_grid];
    double y_array[N_grid];

    if (x > 100.0)
        return (11.0 / 3.0) * pow(x, 2.0);

    if (x < -50.0)
        return  (6.0 / 4.0) * (PI / 2.0) * exp(2.0 * x) + fint_neg12(x) * fint_12(x) / 2.0;

    h = (x - low_bound) / N_grid;
    for (int i = 0; i < N_grid; ++i) {
        x_array[i] = low_bound + i * h;
        y_array[i] = pow(fint_neg12(x_array[i]), 2.0);
    }

    return  (6.0 / 4.0) * rect(x_array, y_array, N_grid) + fint_neg12(x) * fint_12(x) / 2.0;

}



