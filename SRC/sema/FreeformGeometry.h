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
    std::optional<LatticeBipartition>& bipartition;
};

parser::DiagnosticOr<void> build_lattice(Context& ctx);
parser::DiagnosticOr<void> build_hamiltonian(Context& ctx);
parser::DiagnosticOr<void> find_equivalence_classes(Context& ctx);
void find_bipartition_if_exists(Context& ctx);

inline parser::DiagnosticOr<std::unique_ptr<Geometry>> build_geometry(
    parser::DiagnosticEngine& diag,
    parser::ParsedFreeformGeometry const& geometry,
    parser::ParsedFreeformGeometryParameters const& parameters)
{
    auto result = std::make_unique<Geometry>();
    Context ctx {
        .diag = diag,
        .geometry = geometry,
        .parameters = parameters,
        .lattice = result->lattice,
        .hamiltonian = result->hamiltonian,
        .equivalence_classes = result->equivalence_classes,
        .bipartition = result->bipartition,
    };
    TRY(build_lattice(ctx));
    TRY(build_hamiltonian(ctx));
    TRY(find_equivalence_classes(ctx));
    find_bipartition_if_exists(ctx);
    return result;
}

} // namespace dqmc::sema
