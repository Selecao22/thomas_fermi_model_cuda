//
// Created by nikky on 26.12.2019.
//

#include "f_utils.hpp"
#include <cuda_runtime.h>

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
        double* dfi,
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
        dfi[block] = 0;
    }

}

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
        )
{
//    *fi = (double*)calloc(n, sizeof(double));
    *fi = (double*)malloc(n * sizeof(double));
    *x = (double*)calloc(n, sizeof(double));
    *dfi = (double*)calloc(n, sizeof(double));
    double* d_fi;
    double* d_x;
    double* d_dfi;

    cudaMalloc(&d_fi, n * sizeof(double));
    cudaMalloc(&d_x, n * sizeof(double));
    cudaMalloc(&d_dfi, n * sizeof(double));

    get_init_assumpiton_cuda<<<n / 1024, 1024>>>(
            d_dfi,
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
    cudaMemcpy(*dfi, d_dfi, sizeof(double) * n, cudaMemcpyDeviceToHost);

    cudaFree(d_fi);
    cudaFree(d_x);
    cudaFree(d_dfi);
}

