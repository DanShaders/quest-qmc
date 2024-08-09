#pragma once

#include <Eigen/Core>

#include "SRC/common/Types.h"
#include "SRC/parser/DiagnosticEngine.h"
#include "SRC/parser/Forward.h"

namespace dqmc::sema {

using Eigen::Matrix3d, Eigen::Matrix3i, Eigen::MatrixXd,
    Eigen::Vector3d, Eigen::Vector3i;

inline constexpr f64 epsilon = 1e-12;

struct Site {
    std::string label;
    Vector3d cartesian_position;
    Vector3d fractionary_position;
};

struct Lattice {
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

class FreeformGeometry {
public:
    static parser::DiagnosticOr<FreeformGeometry> create(
        parser::DiagnosticEngine& diag,
        parser::ParsedFreeformGeometry const& geometry,
        parser::ParsedFreeformGeometryParameters const& parameters);

    void legacy_compatible_format_into(std::ostream& stream) const;

private:
    FreeformGeometry() = default;

    parser::DiagnosticOr<void> initialize(
        parser::DiagnosticEngine& diag,
        parser::ParsedFreeformGeometry const& geometry,
        parser::ParsedFreeformGeometryParameters const& parameters);

    Lattice m_lattice;
};

} // namespace dqmc::sema
