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
                nu_0 = -pow(q * 3.0 / 2.0, 2.0 / 3.0);
            else if (q < -100)
                nu_0 = log(sqrt(PI) / 2.0 / q);

            else
                nu_0 = 0.5 * log(PI / 6.0) - 1.5 *
                        log(exp(pow(2.0 * pow(q, 2.0 / 3.0), 1.0 / 3.0)) - 1.0);

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


            // classic model
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

                    alpha[k - 1] = -a_c / (b_c + c_c * alpha[k]);
                    beta[k - 1] = (d_c - c_c * beta[k]) / (b_c + c_c * alpha[k]);
                }

                for (int k = 0; k < N - 1; k++)
                    fi[k + 1] = alpha[k] * fi[k] + beta[k];

            }

            nu_0 = -fi[N - 1];

            double pe = pow(2.0 * tet, 2.5) / (6.0 * pow(PI, 2)) * fint_32(-nu_0);
            double psum = 29420.0 * (pe + tet / v);

            auto x2int32 = (double*)calloc(N, sizeof(double));
            auto se_array = (double*)calloc(N, sizeof(double));

            calculate_entrope(x2int32, se_array, fi, H, N);

            auto dif0 = (fi[1] - fi[0]) / pow(H, 2.0);
            auto int_01_calc = rect(x, x2int32, N);

            auto se = (4.0 * sqrt(2.0)* pow(tet, 1.5) * pow(r0, 3.0) / PI) * rect(x, se_array, N);
            auto s = 96.48 * (se + 1.5 * log(1836 * A * tet * pow(v, 2.0 / 3.0) / (2.0 * PI)) + 2.5) / A;

            auto epe = (2.0 * sqrt(2.0) * v * pow(tet, 2.5) / pow(PI, 2.0)) *
                    (fint_32(-nu_0) - 3.0 * int_01_calc);

            auto eke = (3.0 * sqrt(2.0) * v * pow(tet, 2.5) / pow(PI, 2.0)) * int_01_calc;
            auto ee = epe + eke;
            auto e0 = -0.76874512422 * pow(Z, 7.0 / 3.0);
            auto e = 2626.0 * (ee - e0 + 1.5 * tet) / A;

            auto q_quant = (2.0 * sqrt(2.0 * tet) / PI) * pow(r0, 2.0) * fint_neg12(fi[N - 1]);
            auto f_quant = (4.0 * sqrt(2.0 * tet) / PI) * pow(r0, 2.0) *
                    ((7.0 / 4.0) * pow(fint_neg12(fi[N - 1]), 2.5) * fint_neg12_der(fi[N - 1]));
            alpha[N - 2] = 1.0 / (1.0 - 2.0 * H + pow(H, 2.0) * (1.0 + 2.0 * q_quant));
            beta[N - 2] = -2.0 * pow(H, 2.0) * f_quant / (1.0 - 2.0 * H + pow(H, 2.0) * (1.0 + 2.0 * q_quant));


            // quantum correlation
            for (int k = N - 1; k > 0; k--) {
                q_quant = (2.0 * sqrt(2.0 * tet) / PI) * pow(r0, 2.0) *
                        fint_neg12(fi[k] / pow(H * k, 2.0));

                f_quant = (4.0 * sqrt(2.0 * tet) / PI) * pow(r0, 2.0) *
                          ((7.0 / 4.0) * pow(fint_neg12(fi[k] / pow(H * k, 2.0)), 2.5) *
                          fint_12(fi[k] / pow(H * k, 2.0)) * fint_neg12_der(fi[k] / pow(H * k, 2.0)));
                auto a_c = H * (2.0 * k + 1.0);
                auto c_c = H * (2.0 * k - 1.0);
                auto b_c = -4.0 * k * H * (1.0 + 2.0 * pow(H * k, 2.0) * pow(H, 2.0) * q_quant);
                auto d_c = f_quant * 8.0 * pow(H * k, 3.0) * pow(H, 3.0);
                alpha[k - 1] = -a_c / (b_c + c_c * alpha[k]);
                beta[k - 1] = (d_c - c_c * beta[k]) / (b_c + c_c * alpha[k]);
            }

            auto *dfi = (double*)calloc(N, sizeof(double));

            for (int k = 0; k < N - 1; ++k)
                dfi[k + 1] = alpha[k] * dfi[k] + beta[k];


            auto dpe = pow(tet, 2.0) / (3.0 * pow(PI, 3.0)) * (dfi[N - 1] * fint_12(fi[N - 1]) + Y(fi[N]));

            auto dee_array = (double*)calloc(N, sizeof(double));
            auto dse_array = (double*)calloc(N, sizeof(double));

            for (int k = 0; k < N; ++k) {

                dse_array[k] = 2.0 * (k * H) * (pow(k * H, 2.0) * dfi[k] *
                        fint_12(fi[k] / pow(k * H, 2.0)) + 2.0 *
                        pow(k * H, 4.0) * Y(fi[k] / pow(k * H, 2.0)));

                dee_array[k] = 2.0 * (k * H) * (pow(k * H, 2.0) * dfi[k] *
                        fint_12(fi[k] / pow(k * H, 2.0)) + 2.0 *
                        pow(k * H, 4.0) * Y(fi[k] / pow(k * H, 2.0)));
            }

            auto dee = (2.0 * pow(tet, 2.0) * pow(r0, 2.0)) / (3.0 * pow(PI, 2.0)) * rect(x, dee_array, N) +
                    (sqrt(2.0 * tet) * Z / (6.0 * PI)) * (dfi[1] - dfi[0]) /
                    pow(H, 2.0) + 0.269900170 * pow(Z, 5.0 / 3.0);

            auto dmu = sqrt(2.0) * sqrt(tet) / (6.0 * PI) * (0.5 * fint_neg12(fi[N - 1] + dfi[N - 1]));

            delta_array[i * POINT_NUMBER + j] = dee / (ee - e0);

            free(fi);
            free(x);
            free(dfi);
            free(alpha);
            free(beta);
            free(x2int32);
            free(se_array);
            free(dee_array);
            free(dse_array);

            std::cout << (int)(((double)(i * POINT_NUMBER + j)) / (double)(POINT_NUMBER * POINT_NUMBER) * 100.0) << "%\r";
            std::flush(std::cout);
//            std::cout << i * POINT_NUMBER + j << std::endl;
        }
    }

    free(rho_array);
    free(t_array);
    free(delta_array);
}