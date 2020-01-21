#include <omp.h>
#include <iostream>

#include "f_utils.h"

int main(int argc, char** argv)
{
    int N = 10000;
    double H;

    if (argc == 1){
        std::cout << "usage: " << argv[0] << " N\n";
        return 0;
    }

    if (argc == 2) {
        N = atoi(argv[1]);
        if (N <= 0){
            std::cout << "error: N must be from 1 to 100 000 000\n";
            return 1;
        }
    }
    H = 1.0 / N;
    N++;
    double *rho_array = create_physic_array(N);
    double *t_array = create_physic_array(N);

    //linear matrix
    double *delta_array = (double*)calloc(POINT_NUMBER * POINT_NUMBER, sizeof(double));
    int progress = 0;
    double start = omp_get_wtime();
// single tmp
#pragma omp parallel shared(rho_array, t_array, delta_array, progress, H, N, std::cout) default(none)
    {
        int thr = omp_get_thread_num();
        int thr_count = omp_get_num_threads();

        for (int i = thr; i < POINT_NUMBER; i+=thr_count) {
            double t = t_array[i] * 0.001;

            for (int j = 0; j < POINT_NUMBER; ++j) {


                double *alpha = (double*)calloc(N - 1, sizeof(double));
                double *beta = (double*)calloc(N - 1, sizeof(double));

                double *fi = (double*)calloc(N, sizeof(double));
                double *x = (double*)calloc(N, sizeof(double));

                double *x2int32 = (double*)calloc(N, sizeof(double));
                double *se_array = (double*)calloc(N, sizeof(double));

                double *dfi = (double*)calloc(N, sizeof(double));

                double *dee_array = (double*)calloc(N, sizeof(double));
                double *dse_array = (double*)calloc(N, sizeof(double));

                double tet = K * t;

                double r0 = 1.388 * pow(A / rho_array[j], 1.0 / 3.0);
                double q = 0.002795 * Z * rho_array[j] / A / pow(t, 1.5);
                double v = (4.0 / 3.0) * PI * pow(r0, 3.0);
                double aa = 4.0 / PI * sqrt(2.0 * tet) * pow(r0, 2.0);
                double nu_0;

                if (q > 100.0)
                    nu_0 = -pow(q * 3.0 / 2.0, 2.0 / 3.0);
                else if (q < -100.0)
                    nu_0 = log(sqrt(PI) / 2.0 / q);

                else
                    nu_0 = 0.5 * log(PI / 6.0) - 1.5 *
                            log(exp(pow(2.0 * pow(q, 2.0 / 3.0), 1.0 / 3.0)) - 1.0);

                get_init_assumption(
                        x,
                        fi,
                        Z,
                        tet,
                        r0,
                        H,
                        nu_0,
                        N
                );


                // classic model
                for(int lp = 0; lp < 10; lp++)
                {
                    double h_pow2 = pow(H, 2.0);
                    alpha[N - 2] = 1.0 / (1.0 - 2.0 * H + h_pow2 *
                                                          (1.0 + aa * fint_neg12(fi[N - 1])));
                    beta[N - 2] = -aa * h_pow2 *
                                  (2.0 * fint_12(fi[N - 1]) - fi[N - 1] * fint_neg12(fi[N - 1])) /
                                  (1.0 - 2.0 * H + h_pow2 * (1.0 + aa * fint_neg12(fi[N - 1])));

                    for (int k = N - 2; k > 0; k--) {
                        double a_c = H * (2.0 * k + 1.0);
                        double c_c = H * (2.0 * k - 1.0);
                        double b_c = -4.0 * H * k * (1.0 + aa * h_pow2 *
                                                           pow(H * k, 2.0) * fint_neg12(fi[k] / pow(H * k, 2.0)));
                        double d_c = 4.0 * aa * h_pow2 * pow(H * k, 3.0) *
                                     (2.0 * pow(H * k, 2.0) * fint_12(fi[k] / pow(H * k, 2.0)) -
                                      fi[k] * fint_neg12(fi[k] / pow(H * k, 2.0)));

                        alpha[k - 1] = -a_c / (b_c + c_c * alpha[k]);
                        beta[k - 1] = (d_c - c_c * beta[k]) / (b_c + c_c * alpha[k]);
                    }

                    for (int k = 0; k < N - 1; k++)
                        fi[k + 1] = alpha[k] * fi[k] + beta[k];

                }

                nu_0 = -fi[N - 1];

//                double pe = pow(2.0 * tet, 2.5) / (6.0 * pow(PI, 2)) * fint_32(-nu_0);
//                double psum = 29420.0 * (pe + tet / v);

                calculate_entrope(x2int32, se_array, fi, H, N);

//                double dif0 = (fi[1] - fi[0]) / pow(H, 2.0);
                double int_01_calc = rect(x, x2int32, N);

//                double se = (4.0 * sqrt(2.0)* pow(tet, 1.5) * pow(r0, 3.0) / PI) * rect(x, se_array, N);
//                double s = 96.48 * (se + 1.5 * log(1836 * A * tet * pow(v, 2.0 / 3.0) / (2.0 * PI)) + 2.5) / A;

                double epe = (2.0 * sqrt(2.0) * v * pow(tet, 2.5) / pow(PI, 2.0)) *
                           (fint_32(-nu_0) - 3.0 * int_01_calc);

                double eke = (3.0 * sqrt(2.0) * v * pow(tet, 2.5) / pow(PI, 2.0)) * int_01_calc;
                double ee = epe + eke;
                double e0 = -0.76874512422 * pow(Z, 7.0 / 3.0);
//                double e = 2626.0 * (ee - e0 + 1.5 * tet) / A;

                // quantum correlation
                double q_quant = (2.0 * sqrt(2.0 * tet) / PI) * pow(r0, 2.0) * fint_neg12(fi[N - 1]);
                double f_quant = (4.0 * sqrt(2.0 * tet) / PI) * pow(r0, 2.0) *
                               ((7.0 / 4.0) * pow(fint_neg12(fi[N - 1]), 2.0) + 0.5 * fint_12(fi[N - 1]) *
                               fint_neg12_der(fi[N - 1]));

                alpha[N - 2] = 1.0 / (1.0 - 2.0 * H + pow(H, 2.0) * (1.0 + 2.0 * q_quant));
                beta[N - 2] = -2.0 * pow(H, 2.0) * f_quant / (1.0 - 2.0 * H + pow(H, 2.0) * (1.0 + 2.0 * q_quant));


                for (int k = N - 1; k > 0; k--) {
                    q_quant = (2.0 * sqrt(2.0 * tet) / PI) * pow(r0, 2.0) *
                              fint_neg12(fi[k] / pow(H * k, 2.0));

                    f_quant = (4.0 * sqrt(2.0 * tet) / PI) * pow(r0, 2.0) * pow(H * k, 2.0) *
                               ((7.0 / 4.0) * pow(fint_neg12(fi[k] / pow(H * k, 2.0)), 2.0) + 0.5 *
                               fint_12(fi[k] / pow(H * k, 2.0)) * fint_neg12_der(fi[k] / pow(H * k, 2.0)));

                    double a_c = H * (2.0 * k + 1.0);
                    double c_c = H * (2.0 * k - 1.0);
                    double b_c = -4.0 * k * H * (1.0 + 2.0 * pow(H * k, 2.0) * pow(H, 2.0) * q_quant);
                    double d_c = f_quant * 8.0 * pow(H * k, 3.0) * pow(H, 2.0);
                    alpha[k - 1] = -a_c / (b_c + c_c * alpha[k]);
                    beta[k - 1] = (d_c - c_c * beta[k]) / (b_c + c_c * alpha[k]);
                }

                for (int k = 0; k < N - 1; ++k)
                    dfi[k + 1] = alpha[k] * dfi[k] + beta[k];


//                double dpe = pow(tet, 2.0) / (3.0 * pow(PI, 3.0)) * (dfi[N - 1] * fint_12(fi[N - 1]) + Y(fi[N]));

                calculate_dee_and_dse(dse_array, dee_array, fi, dfi, H, N);

                double dee = (2.0 * pow(tet, 2.0) * pow(r0, 3.0)) / (3.0 * pow(PI, 2.0)) * rect(x, dee_array, N) +
                           (sqrt(2.0 * tet) * Z / (6.0 * PI)) * (dfi[1] - dfi[0]) /
                           pow(H, 2.0) + 0.269900170 * pow(Z, 5.0 / 3.0);

//                double dmu = sqrt(2.0) * sqrt(tet) / (6.0 * PI) * (0.5 * fint_neg12(fi[N - 1] + dfi[N - 1]));

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

#pragma omp critical
                {
                    progress++;
                    std::cout << (int)((double)progress / (double)(POINT_NUMBER * POINT_NUMBER) * 100.0) << "%\r";
                    std::cout.flush();
                }


//            std::cout << i * POINT_NUMBER + j << std::endl;
            }
        }
    }
    printf("time: %f\n", omp_get_wtime() - start);


    FILE* rha;
    FILE* ta;
    FILE* da;

    rha = fopen("data3ddy.csv", "wt");
    ta = fopen("data3ddx.csv", "wt");
    da = fopen("test_data_vol.csv", "wt");

    fseek(rha, 0, SEEK_SET);
    fseek(ta, 0, SEEK_SET);
    fseek(da, 0, SEEK_SET);

    for (int i = 0; i < POINT_NUMBER; ++i) {
        fprintf(rha, "%f\t", rho_array[i]);
        fprintf(ta, "%f\t", t_array[i]);
    }

    for (int i = 0; i < POINT_NUMBER; ++i) {
        for (int j = 0; j < POINT_NUMBER; ++j) {
            fprintf(da, "%f\t", log10(fabs(delta_array[i * POINT_NUMBER + j])));
        }
        fprintf(da, "\n");
    }

    free(rho_array);
    free(t_array);
    free(delta_array);

    fclose(rha);
    fclose(ta);
    fclose(da);

    return 0;
}