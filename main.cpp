#include <iostream>
#include <cmath>

#include "f_utils.hpp"

int main()
{
    double *rho_array = create_physic_array(N);
    double *t_array = create_physic_array(N);

    //linear matrix
    auto *delta_array = (double*)calloc(POINT_NUMBER * POINT_NUMBER, sizeof(double));

    for (int i = 0; i < POINT_NUMBER * POINT_NUMBER; ++i) {
        delta_array[i] = 0.0;
    }

    for (int i = 0; i < POINT_NUMBER; ++i) {
        auto t = t_array[i] * 1e-3;

        for (int j = 0; j < POINT_NUMBER; ++j) {
            auto tet = K * t;

            auto r0 = 1.388 * pow((A / rho_array[j]), 1.0/3.0);
            auto q = 2.795e-3 * Z * rho_array[j] / A / pow(t, 1.5);
            auto v = (4.0 / 3.0) * PI * pow(r0, 3);
            auto aa = 4.0 / PI * sqrt(2.0 * tet) * pow(r0, 2.0);
            double nu_0 = 0.0;

            if (q > 100)
            {
                nu_0 = -pow(q* 3.0 / 2.0, 2.0 / 3.0);
            } else if (q < -100)
            {
                nu_0 = log(sqrt(PI) / 2.0 / q);
            } else {
                nu_0 = 0.5 * log(PI / 6.0) - 1.5 * log(exp(pow(2.0 * pow(q, 2.0 / 3.0), 1.0 / 3.0)) - 1.0);
            }

            double *fi;
            double *x;
            double *dfi;


            get_init_assumption(
                    &dfi,
                    &x,
                    &fi,
                    Z,
                    tet,
                    r0,
                    H,
                    nu_0,
                    N
                    );



            free(fi);
            free(x);
            free(dfi);
        }
    }

    free(rho_array);
    free(t_array);
    free(delta_array);
}