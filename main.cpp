#include <iostream>
#include <cmath>

#include "f_utils.hpp"

int main()
{
    double *rho_array = create_physic_array(N);
    double *t_array = create_physic_array(N);

    //linear matrix
    double *delta_array = new double[POINT_NUMBER * POINT_NUMBER];

    for (int i = 0; i < POINT_NUMBER * POINT_NUMBER; ++i) {
        delta_array[i] = 0.0;
    }

    for (int i = 0; i < POINT_NUMBER; ++i) {
        for (int j = 0; j < POINT_NUMBER; ++j) {
            auto te = t_array[i];
            auto t = te * 1e-3;
            auto tet = K * t;

            auto r0 = 1.388 * pow((A / rho_array[j]), 1/3);
            auto q = 2.795e-3 * Z * rho_array[j] / A / pow(t, 1.5); // what will be in power?
            auto v = (4.0 / 3.0) * PI * pow(r0, 3);
            auto aa = 4.0 / PI * sqrt(2.0 * tet) * pow(r0, 2.0);



        }
    }

    delete [] rho_array;
    delete [] t_array;
    delete [] delta_array;
}