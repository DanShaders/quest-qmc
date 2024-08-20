#pragma once

#include "SRC/common/Eigen.h"
#include "SRC/common/Types.h"
#include "SRC/parser/DiagnosticEngine.h"
#include "SRC/parser/Forward.h"

namespace dqmc::sema {

inline constexpr f64 epsilon = 1e-12;

struct Site {
    std::string label;
    Vector3d cartesian_position;
    Vector3d fractionary_position;
};

struct Lattice {
    // Formats structure the same way as old Fortran code did. Useful for debugging purposes.
    void legacy_compatible_format_into(std::ostream& out) const;

    int dimensions;
    // Cartesian components of primitive cell (*)
    Matrix3d basis;

    // Cartesian components of supercell (*)
    Matrix3d supercell_cartesian_basis;
    // Fractionary components of supercell (*)
    Matrix3i supercell_fractionary_basis;

    std::vector<Vector3i> primitive_cells_in_supercell;

    std::vector<Site> primitive_cell_sites;

    std::vector<Site> sites;

    // (*) -- vectors are columns of the matrix
};

struct Hamiltonian {
    struct OnSiteInteration {
        f64 mu_up;
        f64 mu_down;
        f64 u;
    };

    MatrixXd hoppings[2];
    std::vector<OnSiteInteration> interactions;
};

struct Context {
    parser::DiagnosticEngine& diag;
    parser::ParsedFreeformGeometry const& geometry;
    parser::ParsedFreeformGeometryParameters const& parameters;

    Lattice& lattice;
    Hamiltonian& hamiltonian;
};

parser::DiagnosticOr<void> build_lattice(Context& ctx);
parser::DiagnosticOr<void> build_hamiltonian(Context& ctx);

struct Geometry {
    Lattice lattice;
    Hamiltonian hamiltonian;
};

inline parser::DiagnosticOr<Geometry> build_geometry(
    parser::DiagnosticEngine& diag,
    parser::ParsedFreeformGeometry const& geometry,
    parser::ParsedFreeformGeometryParameters const& parameters)
{
    Geometry result;
    Context ctx { diag, geometry, parameters, result.lattice, result.hamiltonian };
    TRY(build_lattice(ctx));
    TRY(build_hamiltonian(ctx));
    return result;
}

} // namespace dqmc::sema
