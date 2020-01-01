cmake_minimum_required(VERSION 3.11)
project(thomas_fermi_model CXX CUDA)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CUDA_STANDARD 11)

enable_language(CUDA)


set(THOMAS_FERMI_SOURCES
        main.cpp
        )

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")

add_library(tf_cuda_improve STATIC f_utils.cu)
set_target_properties(tf_cuda_improve PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON
        POSITION_INDEPENDENT_CODE ON)
target_include_directories(tf_cuda_improve PUBLIC ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})

add_executable(real_model ${THOMAS_FERMI_SOURCES})

target_link_libraries(real_model tf_cuda_improve)