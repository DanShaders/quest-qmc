#pragma once

#include "SRC/common/Types.h"

namespace dqmc::simulation {

struct HubbardModelParameters {
    // Number of time slices to split time evolution operator to for Trotter's approximation.
    u32 time_slices;

    // Imaginary time step per time slice (beta = l * dtau).
    f64 dtau;

    // FIXME: More control for auxiliary field type.
    // FIXME: Allow (un)fixing randomness seed.

    // FIXME: Deduce parameters for Monte-Carlo based on some target accuracy.
    // Number of sweeps to perform to equilibrate the Hubbard-Strotonovich field.
    u64 warm_up_steps;

    // Number of sweeps to perform to measure observables.
    u64 measurement_steps;
};

} // namespace dqmc::simulation
