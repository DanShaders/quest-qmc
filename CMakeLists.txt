cmake_minimum_required(VERSION 3.29)

project(QuestQMC LANGUAGES C CXX Fortran)

find_package(OpenBLAS REQUIRED CONFIG)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED On)
set(CMAKE_EXPORT_COMPILE_COMMANDS On)

add_compile_options(-fPIE)
add_link_options(-fPIE)

if (DQMC_USE_LIBCXX)
    add_compile_options($<$<COMPILE_LANGUAGE:CXX>:-stdlib=libc++>)
    add_link_options(-stdlib=libc++)
endif()
if (NOT DQMC_LINK_RPATH STREQUAL "")
    add_link_options(-rpath "${DQMC_LINK_RPATH}")
endif()

add_compile_options(
    -fopenmp -O3 -Wall
    $<$<COMPILE_LANGUAGE:Fortran>:-std=legacy>
)
add_link_options(-fopenmp)

add_subdirectory(EXAMPLE)
add_subdirectory(SRC)
