#pragma once

#include "SRC/common/Types.h"
#include "SRC/parser/DiagnosticEngine.h"
#include "SRC/parser/Forward.h"
#include "SRC/sema/Geometry.h"

namespace dqmc::sema {

inline constexpr f64 epsilon = 1e-12;

struct Context {
    parser::DiagnosticEngine& diag;
    parser::ParsedFreeformGeometry const& geometry;
    parser::ParsedFreeformGeometryParameters const& parameters;

    Lattice& lattice;
    Hamiltonian& hamiltonian;
    EquivalenceClasses& equivalence_classes;
};

parser::DiagnosticOr<void> build_lattice(Context& ctx);
parser::DiagnosticOr<void> build_hamiltonian(Context& ctx);
parser::DiagnosticOr<void> find_equivalence_classes(Context& ctx);

inline parser::DiagnosticOr<Geometry> build_geometry(
    parser::DiagnosticEngine& diag,
    parser::ParsedFreeformGeometry const& geometry,
    parser::ParsedFreeformGeometryParameters const& parameters)
{
    Geometry result;
    Context ctx {
        .diag = diag,
        .geometry = geometry,
        .parameters = parameters,
        .lattice = result.lattice,
        .hamiltonian = result.hamiltonian,
        .equivalence_classes = result.equivalence_classes,
    };
    TRY(build_lattice(ctx));
    TRY(build_hamiltonian(ctx));
    TRY(find_equivalence_classes(ctx));
    return result;
}

} // namespace dqmc::sema
