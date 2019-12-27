cmake_minimum_required(VERSION 3.11)
project(thomas_fermi_model CXX CUDA)
set(CMAKE_CXX_STANDARD 11)

enable_language(CUDA)

set(THOMAS_FERMI_SOURCES
        main.cpp
        )

add_library(tf_cuda_improve STATIC f_utils.cu)

set_target_properties(tf_cuda_improve
        PROPERTIES CUDA_SEPARABLE_COMPILATION ON
        POSITION_INDEPENDENT_CODE ON)

add_executable(real_model ${THOMAS_FERMI_SOURCES})
target_include_directories(real_model PUBLIC ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})