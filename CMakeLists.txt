cmake_minimum_required(VERSION 3.29)

project(QuestQMC LANGUAGES C CXX Fortran)

find_package(OpenBLAS REQUIRED CONFIG)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED On)

add_compile_options(
    -fopenmp -O3 -Wall
    $<$<COMPILE_LANGUAGE:Fortran>:-std=legacy>
)
add_link_options(-fopenmp)

add_subdirectory(EXAMPLE)
add_subdirectory(SRC)
