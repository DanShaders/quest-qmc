set(SOURCES
    dqmc_2dperl.F90
    dqmc_bonds.F90
        detour_dqmc_cfg.F90
        detour_dqmc_cfg_defs.F90
        dqmc_cfg.cpp
        dqmc_cfg.F90
    dqmc_geom_param.F90
    dqmc_geom_wrap.F90
    dqmc_gfun.F90
    dqmc_gtau.F90
    dqmc_hamilt.F90
    dqmc_hubbard.F90
    dqmc_kbonds.F90
    dqmc_latt.F90
    dqmc_matb.F90
    dqmc_mpi.F90
    dqmc_omp.F90
    dqmc_phy0.F90
    dqmc_phy2.F90
    dqmc_reclatt.F90
    dqmc_seqb.F90
    dqmc_struct.F90
    dqmc_symm.F90
    dqmc_tdm1.F90
    dqmc_tdm2.F90
        detour_dqmc_util.F90
        detour_dqmc_util_defs.F90
        dqmc_util.cpp
        dqmc_util.F90
    dqmc_wspace.F90

    binding.F90
    config.cpp

    blaslapack.c++
    computeg.c++
    kernel_cpu.c++
    profile.c++
)

add_library(dqmc STATIC ${SOURCES})
target_compile_options(dqmc PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:-Wno-c-binding-type>)
target_link_libraries(dqmc PUBLIC OpenBLAS::OpenBLAS FastFloat::fast_float)
target_include_directories(dqmc INTERFACE ${CMAKE_CURRENT_BINARY_DIR})
