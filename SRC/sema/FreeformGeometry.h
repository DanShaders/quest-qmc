#pragma once

#include <Eigen/Core>

#include "SRC/common/Error.h"
#include "SRC/common/Types.h"
#include "SRC/parser/Forward.h"

namespace dqmc::sema {

using Eigen::Matrix3d, Eigen::Matrix3i, Eigen::Vector3d, Eigen::Vector3i;

inline constexpr f64 epsilon = 1e-12;

class FreeformGeometry {
public:
    class Site {
    public:
        Site(std::string label, Vector3d cartesian_position, Vector3d fractionary_position)
            : m_label(std::move(label))
            , m_cartesian_position(std::move(cartesian_position))
            , m_fractionary_position(std::move(fractionary_position))
        {
        }

        std::string const& label() const& { return m_label; }
        Vector3d const& cartesian_position() const& { return m_cartesian_position; }
        Vector3d const& fractionary_position() const& { return m_fractionary_position; }

    private:
        std::string m_label;
        Vector3d m_cartesian_position;
        Vector3d m_fractionary_position;
    };

    static std::expected<FreeformGeometry, Empty> create(
        parser::DiagnosticEngine& diag,
        parser::ParsedFreeformGeometry const& geometry,
        parser::ParsedFreeformGeometryParameters const& parameters);

    void legacy_compatible_format_into(std::ostream& stream) const;

    int dimensions() const { return m_dimensions; }
    Matrix3d const& lattice_basis() const& { return m_lattice_basis; }
    Matrix3d const& supercell_cartesian_basis() const& { return m_supercell_cartesian_basis; }
    Matrix3i const& supercell_fractionary_basis() const& { return m_supercell_fractionary_basis; }

    int primitive_cells_in_supercell_count() const { return m_primitive_cells_in_supercell.size(); }
    std::vector<Vector3i> const& primitive_cells_in_supercell() const& { return m_primitive_cells_in_supercell; }

    int primitive_cell_sites_count() const { return m_primitive_cell_sites.size(); }
    std::vector<Site> const& primitive_cell_sites() const { return m_primitive_cell_sites; }

    int sites_count() const { return m_sites.size(); }
    std::vector<Site> const& sites() const { return m_sites; }

private:
    FreeformGeometry() = default;

    std::expected<void, Empty> initialize(
        parser::DiagnosticEngine& diag,
        parser::ParsedFreeformGeometry const& geometry,
        parser::ParsedFreeformGeometryParameters const& parameters);

    int m_dimensions = 0;
    // Cartesian components of primitive cell (*)
    Matrix3d m_lattice_basis;

    // Cartesian components of supercell (*)
    Matrix3d m_supercell_cartesian_basis;
    // Fractionary components of supercell (*)
    Matrix3i m_supercell_fractionary_basis;

    std::vector<Vector3i> m_primitive_cells_in_supercell;

    std::vector<Site> m_primitive_cell_sites;

    std::vector<Site> m_sites;

    // (*) -- vectors are columns of the matrix
};

} // namespace dqmc::sema
