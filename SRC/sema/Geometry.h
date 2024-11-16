#pragma once

#include <map>

#include "SRC/common/Eigen.h"
#include "SRC/common/Types.h"

namespace dqmc::sema {

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
    Matrix3i supercell_basis_adjugate;
    // Map from primitive_cells_in_supercell elements into their index
    std::map<Vector3i, int, Vector3iComparator> cell_by_fractional_coords;

    // (*) -- vectors are columns of the matrix
};

struct Hamiltonian {
    struct OnSiteInteration {
        f64 mu_up;
        f64 mu_down;
        f64 u;

        bool operator==(OnSiteInteration const&) const = default;
    };

    MatrixXd hoppings[2];
    std::vector<OnSiteInteration> interactions;
};

struct EquivalenceClasses {
    // Equivalence classes of Green's function. For every equivalence class k, there exists a real
    // number v such that for all sites i and j such that `pair_class[i, j] == k`,
    // `G_ij == greens_fn_phase[i, j] * v`.
    MatrixXi pair_class;

    MatrixXi greens_fn_phase; // 1, -1, or 0
};

struct Geometry {
    Lattice lattice;
    Hamiltonian hamiltonian;
    EquivalenceClasses equivalence_classes;
};

} // namespace dqmc::sema
