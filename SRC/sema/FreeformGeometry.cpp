#include <Eigen/LU>
#include <set>

#include "SRC/parser/FreeformGeometryParser.h"
#include "SRC/sema/FreeformGeometry.h"

namespace dqmc::sema {

namespace {

std::vector<Vector3i> compute_primitive_cells_in_supercell(FreeformGeometry const& geometry)
{
    // Compute the number of points with integer fractionary coordinates in a supercell. These will
    // be precisely points we are interested to find explicitly.
    int det = geometry.supercell_fractionary_basis().determinant();

    Matrix3i supercell_basis_inverse; // = (supercell_fractionary_basis)^-1 * det(supercell_fractionary_basis)
    {
        auto const& m = geometry.supercell_fractionary_basis();
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                supercell_basis_inverse(j, i)
                    = m((i + 1) % 3, (j + 1) % 3) * m((i + 2) % 3, (j + 2) % 3)
                    - m((i + 2) % 3, (j + 1) % 3) * m((i + 1) % 3, (j + 2) % 3);
            }
        }
    }
    if (det < 0) {
        supercell_basis_inverse = -supercell_basis_inverse;
        det = -det;
    }

    std::vector<Vector3i> result = { { 0, 0, 0 } };
    std::set<
        Vector3i,
        decltype([](Vector3i const& a, Vector3i const& b) {
            return std::tie(a[0], a[1], a[2]) < std::tie(b[0], b[1], b[2]);
        })>
        seen_vertices;
    seen_vertices.emplace(0, 0, 0);

    // To find points themselves, we want to do BFS on the graph where vertices are points with
    // integer fractionary coordinates and edges connect neighboring (distance = 1) points. This
    // simple algorithm would have worked if integer points were always connected in every possible
    // supercell. However, this is not the case (consider basis {(3, 1, 0), (5, 1, 0), (0, 0, 1)}).
    // Fortunately, there is an easy fix: if one uses a quotient of the full (infinite) primary
    // lattice graph under equivalence classes induced by superlattice, then every equivalence class
    // (and, therefore, every integer fractionary point) will be reached.
    for (size_t i = 0; i < result.size(); ++i) {
        auto site = result[i];

        for (int direction = 0; direction < 3; ++direction) {
            for (int sign : { -1, 1 }) {
                Vector3i neighbor = site + supercell_basis_inverse.col(direction) * sign;

                for (int component = 0; component < 3; ++component) {
                    int& coordinate = neighbor[component];
                    coordinate -= coordinate / det * det;
                    if (coordinate < 0) {
                        coordinate += det;
                    }
                }

                if (!seen_vertices.contains(neighbor)) {
                    seen_vertices.insert(neighbor);
                    result.push_back(neighbor);
                }
            }
        }
    }

    VERIFY(result.size() == det);

    for (auto& primitive_cell : result) {
        primitive_cell = geometry.supercell_fractionary_basis() * primitive_cell / det;
    }

    // Order primitive cells in the same way as in the legacy code.
    std::ranges::sort(result, [](Vector3i const& a, Vector3i const& b) {
        return std::tie(a[2], a[1], a[0]) < std::tie(b[2], b[1], b[0]);
    });

    return result;
}

std::vector<FreeformGeometry::Site> compute_sites(FreeformGeometry const& geometry)
{
    using Site = FreeformGeometry::Site;

    auto translate = [&](Site const& site, Vector3i const& translation) {
        Vector3d fractionary_coordinates = site.fractionary_position() + translation.cast<f64>();
        Vector3d cartesian_coordinates = geometry.lattice_basis() * fractionary_coordinates;
        return Site { site.label(), cartesian_coordinates, fractionary_coordinates };
    };

    std::vector<Site> result;
    for (Vector3i const& cell : geometry.primitive_cells_in_supercell()) {
        for (Site const& site : geometry.primitive_cell_sites()) {
            result.push_back(translate(site, cell));
        }
    }
    return result;
}

} // namespace

std::expected<FreeformGeometry, Empty> FreeformGeometry::create(parser::DiagnosticEngine& diag, parser::ParsedFreeformGeometry const& geometry)
{
    FreeformGeometry result;
    TRY(result.initialize(diag, geometry));
    return result;
}

void FreeformGeometry::legacy_compatible_format_into(std::ostream& stream) const
{
    std::print(stream, " ================================================================\n"
                       " Basic real space geometry info\n"
                       "\n"
                       " Crystal atomic basis\n");
    for (int i = 0; i < primitive_cell_sites_count(); ++i) {
        auto const& coords = primitive_cell_sites()[i].fractionary_position();
        std::println(stream, "{:3}{:14.7f}{:14.7f}{:14.7f}", i, coords[0], coords[1], coords[2]);
    }

    std::print(stream, "\n"
                       " Basis cell vectors\n");
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            std::print(stream, "{:14.7f}", lattice_basis()(j, i));
        }
        std::print(stream, "\n");
    }

    std::print(stream, "\n"
                       "\n"
                       " Supercell vectors (fractionary unit)\n");
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            std::print(stream, "{:5}", supercell_fractionary_basis()(j, i));
        }
        std::print(stream, "\n");
    }

    std::print(stream, "\n"
                       " Super-Lattice vectors (cartesian)\n");
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            std::print(stream, "{:14.7f}", supercell_cartesian_basis()(j, i));
        }
        std::print(stream, "\n");
    }

    std::print(stream, "\n"
                       " ================================================================\n"
                       " Real space lattice\n"
                       "\n");
    std::print(stream, " Number of orbitals in primitive cell: {:12}\n", primitive_cell_sites_count());
    std::print(stream, " Total number of orbitals:             {:12}\n", sites_count());
    std::print(stream, " index  label   type       X           Y         Z   \n");
    for (int i = 0; i < sites_count(); ++i) {
        auto const& site = sites()[i];
        auto const& coords = site.cartesian_position();
        std::print(stream, "{:3} {:3} {:3}{:14.5f}{:14.5f}{:14.5f}\n",
            i, site.label(), i % primitive_cell_sites_count(), coords[0], coords[1], coords[2]);
    }
}

std::expected<void, Empty> FreeformGeometry::initialize(parser::DiagnosticEngine& diag, parser::ParsedFreeformGeometry const& geometry)
{
    // dimensions
    VERIFY(geometry.dimensions >= 1 && geometry.dimensions <= 3);
    m_dimensions = geometry.dimensions;

    // lattice_basis
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            f64 coefficient = geometry.lattice_basis[i][j];
            m_lattice_basis(j, i) = coefficient;
            if (i >= dimensions() || j >= dimensions()) {
                VERIFY((i == j && coefficient == 1e3) || (i != j && coefficient == 0));
            }
        }
    }

    if (std::abs(m_lattice_basis.determinant()) < epsilon) {
        diag.error(geometry.lattice_basis_location,
            "lattice basis is not linearly independent");
        return Empty::error();
    }

    // supercell_fractionary_basis
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            int coefficient = geometry.supercell_basis[i][j];
            m_supercell_fractionary_basis(j, i) = coefficient;
            if (i >= dimensions() || j >= dimensions()) {
                VERIFY((i == j && coefficient == 1) || (i != j && coefficient == 0));
            }
        }
    }

    if (supercell_fractionary_basis().determinant() == 0) {
        diag.error(geometry.supercell_basis_location,
            "supercell does not contain any primary cells");
        return Empty::error();
    }

    // supercell_cartesian_basis
    m_supercell_cartesian_basis = lattice_basis() * supercell_fractionary_basis().cast<f64>();

    // primitive_cells_in_supercell
    m_primitive_cells_in_supercell = compute_primitive_cells_in_supercell(*this);

    // primitive_cell_sites
    Matrix3d lattice_basis_inverse = lattice_basis().inverse();
    for (auto const& site : geometry.primitive_cell_sites) {
        Vector3d cartesian_position { site.displacement[0], site.displacement[1], site.displacement[2] };
        m_primitive_cell_sites.emplace_back(
            site.label,
            std::move(cartesian_position),
            lattice_basis_inverse * cartesian_position);
    }

    // sites
    m_sites = compute_sites(*this);

    return {};
}

} // namespace dqmc::sema
