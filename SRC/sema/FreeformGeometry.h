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
    Vector3d fractional_position;
};

struct Lattice {
    struct CellLookupResult {
        int primitive_cell;
        Vector3i supercell;
    };

    // Returns struct `result` such that `primitive_cells_in_supercell[result.primitive_cell] +
    // supercell_fractional_basis * result.supercell == coords`.
    CellLookupResult fractional_coords_to_cell(Vector3i const& coords) const;

    // Formats structure the same way as old Fortran code did. Useful for debugging purposes.
    void legacy_compatible_format_into(std::ostream& out) const;

    int dimensions;
    // Cartesian components of primitive cell (*)
    Matrix3d basis;

    // Cartesian components of supercell (*)
    Matrix3d supercell_cartesian_basis;
    // Fractional components of supercell (*)
    Matrix3i supercell_fractional_basis;

    std::vector<Vector3i> primitive_cells_in_supercell;

    std::vector<Site> primitive_cell_sites;

    std::vector<Site> sites;

    // abs(det(supercell_fractional_basis)), same as primitive_cells_in_supercell.size()
    int supercell_size;
    // supercell_fractional_basis^(-1) * supercell_size
    Matrix3i supercell_basis_inverse;
    // Map from primitive_cells_in_supercell elements into their index
    std::map<Vector3i, int, Vector3iComparator> cell_by_fractional_coords;

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
