//
// Created by nikky on 26.12.2019.
//
#include "f_utils.hpp"

double* create_physic_array(int size)
{
    auto *array = (double*)calloc(size, sizeof(double));
    double point_step = 0.1;

    array[0] = 0.1;

    for (int i = 1; i < POINT_NUMBER; ++i) {
        array[i] = array[i - 1] + point_step;
        if (((i % 10 == 9) && i < 11 ) || ((i % 9) == 0 && i > 11))
            point_step*=10;
    }

    return array;
}


__global__
void
get_init_assumpiton_cuda
(
        double* x,
        double* fi,
        double z,
        double tet,
        double r0,
        double h,
        double nu_0,
        int n
        )
{
    unsigned int block = blockDim.x * blockIdx.x + threadIdx.x;
    double hi = block * h;

    if (block < n)
    {
        fi[block] = (z / tet / r0) * (1.0 - 1.5 * pow(hi, 2.0) + 0.5 * pow(hi, 6.0)) - nu_0 * pow(hi, 2.0);
        x[block] = hi;
    }

}

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
        )
{
    *fi = (double*)calloc(n, sizeof(double));
    *x = (double*)calloc(n, sizeof(double));

    double* d_fi;
    double* d_x;

    cudaMalloc(&d_fi, n * sizeof(double));
    cudaMalloc(&d_x, n * sizeof(double));

    get_init_assumpiton_cuda<<<(n / 1024) + 1, 1024>>>(
            d_x,
            d_fi,
            z,
            tet,
            r0,
            h,
            nu_0,
            n
    );

    cudaMemcpy(*fi, d_fi, sizeof(double) * n, cudaMemcpyDeviceToHost);
    cudaMemcpy(*x, d_x, sizeof(double) * n, cudaMemcpyDeviceToHost);

    cudaFree(d_fi);
    cudaFree(d_x);
}

__global__
void
calculate_approximation_cuda(
        double *alpha,
        double *beta,
        double *fi,
        int n
        )
{
    int block = blockIdx.x * blockDim.x + threadIdx.x;

    if (block < n - 1){
        fi[block + 1] = alpha[block] * fi[block] + beta[block];
    }
}

void
calculate_approximation(
        double *alpha,
        double *beta,
        double *fi,
        int n
        )
{
    double *d_alpha;
    double *d_beta;
    double *d_fi;

    auto size_of_array_ab = sizeof(double) * (n - 1);
    auto size_of_array_fi = sizeof(double) * n;

    cudaMalloc(&d_alpha, size_of_array_ab);
    cudaMalloc(&d_beta, size_of_array_ab);
    cudaMalloc(&d_fi, size_of_array_fi);

    cudaMemcpy(d_alpha, alpha, size_of_array_ab, cudaMemcpyHostToDevice);
    cudaMemcpy(d_beta, beta, size_of_array_ab, cudaMemcpyHostToDevice);
    cudaMemcpy(d_fi, fi, size_of_array_fi, cudaMemcpyHostToDevice);

    calculate_approximation_cuda<<<(n / 1024) + 1, 1024>>>(d_alpha, d_beta, d_fi, n);

    cudaMemcpy(fi, d_fi, size_of_array_fi, cudaMemcpyDeviceToHost);

    cudaFree(d_alpha);
    cudaFree(d_beta);
    cudaFree(d_fi);
}

__global__
void
calculate_entrope_cuda(
        double* x2int32,
        double* se_array,
        double* fi,
        double h,
        double n
)
{
    int block = blockDim.x * blockIdx.x + threadIdx.x;

    if (block < n)
    {
        x2int32[block] = 2.0 * pow(block * h, 5.0) * fint_32(fi[block] / pow(block * h, 2.0));
        se_array[block] = 2.0 * pow(block * h, 5.0) * ((5.0 / 3.0) * fint_32(fi[block] / pow(block * h, 2.0)) -
                (fi[block] / pow(block * h, 2.0)) * fint_12(fi[block] / pow(block * h, 2.0)));
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
    double *d_x2int32;
    double *d_se_array;
    double *d_fi;
    auto size_of_array = sizeof(double) * n;

    cudaMalloc(&d_x2int32, size_of_array);
    cudaMalloc(&d_se_array, size_of_array);
    cudaMalloc(&d_fi, size_of_array);

    cudaMemcpy(d_fi, fi, size_of_array, cudaMemcpyHostToDevice);

    calculate_entrope_cuda <<< (n / 1024) + 1, 1024 >>> (
            d_x2int32,
            d_se_array,
            d_fi,
            h,
            n
            );

    cudaMemcpy(x2int32, d_x2int32, size_of_array, cudaMemcpyDeviceToHost);
    cudaMemcpy(se_array, d_se_array, size_of_array, cudaMemcpyDeviceToHost);

    cudaFree(d_x2int32);
    cudaFree(d_se_array);
    cudaFree(d_fi);
}

__host__
__device__
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

__host__
__device__
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

__host__
__device__
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