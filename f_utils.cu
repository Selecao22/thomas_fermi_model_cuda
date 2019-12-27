//
// Created by nikky on 26.12.2019.
//

#include <new>
#include "f_utils.hpp"

double* create_physic_array(int size)
{
    double *array;

    try {
         array = new double [size];
    }
    catch (std::bad_alloc) {
        return nullptr;
    }
    double point_step = 0.1;

    array[0] = 0.1;

    for (int i = 1; i < POINT_NUMBER; ++i) {
        array[i] = array[i - 1] + point_step;
        if (((i % 10 == 9) && i < 11 ) || ((i % 9) == 0 && i > 11))
            point_step*=10;
    }

    return array;
}

