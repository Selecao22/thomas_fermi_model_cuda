#include <iostream>

#include "f_utils.hpp"

int main()
{
    double *rho_array = create_physic_array(N);
    double *t_array = create_physic_array(N);

    //linear matrix
    auto *delta_array = (double*)calloc(POINT_NUMBER * POINT_NUMBER, sizeof(double));

    for (int i = 0; i < POINT_NUMBER; ++i) {
        auto t = t_array[i] * 1e-3;

        for (int j = 0; j < POINT_NUMBER; ++j) {
            auto tet = K * t;

            auto r0 = 1.388 * pow((A / rho_array[j]), 1.0/3.0);
            auto q = 2.795e-3 * Z * rho_array[j] / A / pow(t, 1.5);
            auto v = (4.0 / 3.0) * PI * pow(r0, 3);
            auto aa = 4.0 / PI * sqrt(2.0 * tet) * pow(r0, 2.0);
            double nu_0;

            if (q > 100)
            {
                nu_0 = -pow(q* 3.0 / 2.0, 2.0 / 3.0);
            } else if (q < -100)
            {
                nu_0 = log(sqrt(PI) / 2.0 / q);
            } else {
                nu_0 = 0.5 * log(PI / 6.0) - 1.5 *
                        log(exp(pow(2.0 * pow(q, 2.0 / 3.0), 1.0 / 3.0)) - 1.0);
            }

            double *fi;
            double *x;

            get_init_assumption(
                    &x,
                    &fi,
                    Z,
                    tet,
                    r0,
                    H,
                    nu_0,
                    N
                    );

            auto *alpha = (double*)calloc(N - 1, sizeof(double));
            auto *beta = (double*)calloc(N - 1, sizeof(double));

            for(int lp = 0; lp < 10; lp++)
            {
                double h_pow2 = pow(H, 2.0);
                alpha[N - 2] = 1.0 / (1.0 - 2.0 * H + h_pow2 *
                        (1.0 + aa* fint_neg12((fi[N - 1]))));
                beta[N - 2] = -aa * h_pow2 *
                        (fint_12(fi[N - 1]) - fi[N - 1] * fint_neg12(fi[N - 1])) /
                        (1.0 - 2.0 * H + h_pow2 * (1.0 + aa * fint_neg12(fi[N - 1])));

                for (int k = N - 2; k > 0; k--) {
                    double a_c = H * (2.0 * k + 1.0);
                    double c_c = H * (2.0 * k - 1.0);
                    double b_c = -4.0 * H * k * (1.0 + aa * h_pow2 *
                            pow(H * k, 2.0) * fint_neg12(fi[k] / pow(H * k, 2.0)));
                    double d_c = 4.0 * aa * h_pow2 * pow(H * k, 3.0) *
                            (2.0 * h_pow2 * fint_12(fi[k] / pow(H * k, 2.0)) -
                            fi[k] * fint_neg12(fi[k] / pow(H * k, 2.0)));

                    alpha[N - 2] = -a_c / (b_c + c_c * alpha[k]);
                    beta[N - 2] = (d_c - c_c * beta[k]) / (b_c + c_c * alpha[k]);
                }

                for (int k = 0; k < N - 1; ++k) {
                    fi[k + 1] = alpha[k] * fi[k] + beta[k];
                }

            }



            auto *dfi = (double*)calloc(N, sizeof(double));

            free(fi);
            free(x);
            free(dfi);
            free(alpha);
            free(beta);

//            std::cout << (int)(((double)(i * POINT_NUMBER + j)) / (double)(POINT_NUMBER * POINT_NUMBER) * 100.0) << "%\r";
//            std::flush(std::cout);
            std::cout << i * POINT_NUMBER + j << std::endl;
        }
    }

    free(rho_array);
    free(t_array);
    free(delta_array);
}