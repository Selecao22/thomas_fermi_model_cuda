cmake_minimum_required(VERSION 3.11)
project(thomas_fermi_model)
set(CMAKE_C_STANDARD 11)

set(THOMAS_FERMI_SOURCES
        main.c
        )

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fopenmp -lm")

add_executable(real_model ${THOMAS_FERMI_SOURCES})

add_library(tf_cuda_improve STATIC f_utils.c)
target_link_libraries(real_model tf_cuda_improve)