#include <Eigen/LU>
#include <set>

#include "SRC/common/Enumerate.h"
#include "SRC/parser/FreeformGeometryParser.h"
#include "SRC/sema/FreeformGeometry.h"

namespace dqmc::sema {

namespace {

struct LatticeBuildingContext : Context {
    LatticeBuildingContext(Context& ctx)
        : Context(ctx)
    {
    }

    void compute_primitive_cells_in_supercell();
    void compute_sites();
    parser::DiagnosticOr<void> build_lattice();
};

// Computes primitive cells (points with integer fractional coordinates) in a parallelogram
// specified by supercell basis.
void LatticeBuildingContext::compute_primitive_cells_in_supercell()
{
    // Compute the number of points with integer fractional coordinates in a supercell. These will
    // be precisely points we are interested to find explicitly.
    std::tie(lattice.supercell_size, lattice.supercell_basis_adjugate)
        = compute_adjugate(lattice.supercell_fractional_basis);

    std::vector<Vector3i> result = { { 0, 0, 0 } };
    std::set<Vector3i, Vector3iComparator> seen_vertices;
    seen_vertices.emplace(0, 0, 0);

    // To find points themselves, we want to do BFS on the graph where vertices are points with
    // integer fractional coordinates and edges connect neighboring (distance = 1) points. This
    // simple algorithm would have worked if integer points were always connected in every possible
    // supercell. However, this is not the case (consider basis {(3, 1, 0), (5, 1, 0), (0, 0, 1)}).
    // Fortunately, there is an easy fix: if one uses a quotient of the full (infinite) primary
    // lattice graph under equivalence classes induced by superlattice, then every equivalence class
    // (and, therefore, every integer fractional point) will be reached. In practice, this means
    // that if step out of the supercell when considering an edge, we apply a linear combination of
    // supercell basis vectors to the point to get back into the supercell.
    for (size_t i = 0; i < result.size(); ++i) {
        auto site = result[i];

        for (int direction = 0; direction < 3; ++direction) {
            Vector3i neighbor = site + lattice.supercell_basis_adjugate.col(direction);

            for (int component = 0; component < 3; ++component) {
                int& coordinate = neighbor[component];
                coordinate %= lattice.supercell_size;
                if (coordinate < 0) {
                    coordinate += lattice.supercell_size;
                }
            }

            if (!seen_vertices.contains(neighbor)) {
                seen_vertices.insert(neighbor);
                result.push_back(neighbor);
            }
        }
    }

    VERIFY(result.size() == lattice.supercell_size);

    for (auto& primitive_cell : result) {
        primitive_cell = lattice.supercell_fractional_basis * primitive_cell / lattice.supercell_size;
    }

    // Order primitive cells in the same way as in the legacy code.
    std::ranges::sort(result, [](Vector3i const& a, Vector3i const& b) {
        return std::tie(a[2], a[1], a[0]) < std::tie(b[2], b[1], b[0]);
    });

    lattice.primitive_cells_in_supercell = result;

    for (auto [i, site] : enumerate(lattice.primitive_cells_in_supercell)) {
        lattice.cell_by_fractional_coords[site] = i;
    }
}

// Computes lattice sites by copying and translating sites specified in #ORB section to all
// primitive cells found in compute_primitive_cells_in_supercell.
void LatticeBuildingContext::compute_sites()
{
    auto translate = [&](Site const& site, Vector3i const& translation) {
        Vector3d fractional_coordinates = site.fractional_position + translation.cast<f64>();
        Vector3d cartesian_coordinates = lattice.basis * fractional_coordinates;
        return Site { site.label, cartesian_coordinates, fractional_coordinates };
    };

    for (Vector3i const& cell : lattice.primitive_cells_in_supercell) {
        for (Site const& site : lattice.primitive_cell_sites) {
            lattice.sites.push_back(translate(site, cell));
        }
    }
}

// Populates fields of `Context::lattice` field.
parser::DiagnosticOr<void> LatticeBuildingContext::build_lattice()
{
    // dimensions
    VERIFY(geometry.dimensions >= 1 && geometry.dimensions <= 3);
    lattice.dimensions = geometry.dimensions;
    int dims = lattice.dimensions;

    // basis
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            f64 coefficient = geometry.lattice_basis[i][j];
            lattice.basis(j, i) = coefficient;
            if (i >= dims || j >= dims) {
                VERIFY((i == j && coefficient == 1e3) || (i != j && coefficient == 0));
            }
        }
    }

    if (std::abs(lattice.basis.determinant()) < epsilon) {
        return diag.error(geometry.lattice_basis_location,
            "lattice basis is not linearly independent");
    }

    // supercell_fractional_basis
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            int coefficient = geometry.supercell_basis[i][j];
            lattice.supercell_fractional_basis(j, i) = coefficient;
            if (i >= dims || j >= dims) {
                VERIFY((i == j && coefficient == 1) || (i != j && coefficient == 0));
            }
        }
    }

    if (lattice.supercell_fractional_basis.determinant() == 0) {
        return diag.error(geometry.supercell_basis_location,
            "supercell does not contain any primary cells");
    }

    // supercell_cartesian_basis
    lattice.supercell_cartesian_basis = lattice.basis * lattice.supercell_fractional_basis.cast<f64>();

    // primitive_cells_in_supercell
    compute_primitive_cells_in_supercell();

    // primitive_cell_sites
    Matrix3d lattice_basis_inverse = lattice.basis.inverse();
    for (auto const& site : geometry.primitive_cell_sites) {
        Vector3d cartesian_position { site.displacement[0], site.displacement[1], site.displacement[2] };
        lattice.primitive_cell_sites.emplace_back(
            site.label,
            std::move(cartesian_position),
            lattice_basis_inverse * cartesian_position);
    }

    // sites
    compute_sites();

    return {};
}

} // namespace

parser::DiagnosticOr<void> build_lattice(Context& ctx)
{
    return LatticeBuildingContext { ctx }.build_lattice();
}

auto Lattice::fractional_coords_to_cell(Vector3i const& coords) const -> CellLookupResult
{
    CellLookupResult result;

    Vector3i coords_to_lookup = supercell_basis_adjugate * coords;
    for (auto [i, component] : enumerate(coords_to_lookup)) {
        int coordinate_inside_supercell = component % supercell_size;
        if (coordinate_inside_supercell < 0) {
            coordinate_inside_supercell += supercell_size;
        }
        result.supercell[i] = (coords_to_lookup[i] - coordinate_inside_supercell) / supercell_size;
        component = coordinate_inside_supercell;
    }

    coords_to_lookup = supercell_fractional_basis * coords_to_lookup / supercell_size;

    result.primitive_cell = cell_by_fractional_coords.at(coords_to_lookup);

    return result;
}

void Lattice::legacy_compatible_format_into(std::ostream& stream) const
{
    std::print(stream, " ================================================================\n"
                       " Basic real space geometry info\n"
                       "\n"
                       " Crystal atomic basis\n");
    for (auto [i, site] : enumerate(primitive_cell_sites)) {
        auto const& coords = site.fractional_position;
        std::println(stream, "{:3}{:14.7f}{:14.7f}{:14.7f}", i, coords[0], coords[1], coords[2]);
    }

    std::print(stream, "\n"
                       " Basis cell vectors\n");
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            std::print(stream, "{:14.7f}", basis(j, i));
        }
        std::print(stream, "\n");
    }

    std::print(stream, "\n"
                       "\n"
                       " Supercell vectors (fractional unit)\n");
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            std::print(stream, "{:5}", supercell_fractional_basis(j, i));
        }
        std::print(stream, "\n");
    }

    std::print(stream, "\n"
                       " Super-Lattice vectors (cartesian)\n");
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            std::print(stream, "{:14.7f}", supercell_cartesian_basis(j, i));
        }
        std::print(stream, "\n");
    }

    std::print(stream, "\n"
                       " ================================================================\n"
                       " Real space lattice\n"
                       "\n");
    std::print(stream, " Number of orbitals in primitive cell: {:12}\n", primitive_cell_sites.size());
    std::print(stream, " Total number of orbitals:             {:12}\n", sites.size());
    std::print(stream, " index  label   type       X           Y         Z   \n");
    for (auto [i, site] : enumerate(sites)) {
        auto const& coords = site.cartesian_position;
        std::print(stream, "{:3} {:3} {:3}{:14.5f}{:14.5f}{:14.5f}\n",
            i, site.label, i % primitive_cell_sites.size(), coords[0], coords[1], coords[2]);
    }
}

} // namespace dqmc::sema
