#include <Eigen/LU>
#include <set>

#include "SRC/parser/FreeformGeometryParametersParser.h"
#include "SRC/parser/FreeformGeometryParser.h"
#include "SRC/sema/FreeformGeometry.h"

namespace dqmc::sema {

namespace {

using Vector3iComparator = decltype([](Vector3i const& a, Vector3i const& b) {
    return std::tie(a[0], a[1], a[2]) < std::tie(b[0], b[1], b[2]);
});

// Returns { abs(det(basis)), (basis)^-1 * abs(det) }.
std::pair<int, Matrix3i> compute_premultiplied_inverse(Matrix3i const& basis)
{
    int det = basis.determinant();

    Matrix3i inverse;
    {
        auto const& m = basis;
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                inverse(j, i)
                    = m((i + 1) % 3, (j + 1) % 3) * m((i + 2) % 3, (j + 2) % 3)
                    - m((i + 2) % 3, (j + 1) % 3) * m((i + 1) % 3, (j + 2) % 3);
            }
        }
    }
    if (det < 0) {
        inverse = -inverse;
        det = -det;
    }

    return { det, inverse };
}

std::vector<Vector3i> compute_primitive_cells_in_supercell(Lattice const& lattice)
{
    // Compute the number of points with integer fractionary coordinates in a supercell. These will
    // be precisely points we are interested to find explicitly.
    auto [det, supercell_basis_inverse] = compute_premultiplied_inverse(lattice.supercell_fractionary_basis);

    std::vector<Vector3i> result = { { 0, 0, 0 } };
    std::set<Vector3i, Vector3iComparator> seen_vertices;
    seen_vertices.emplace(0, 0, 0);

    // To find points themselves, we want to do BFS on the graph where vertices are points with
    // integer fractionary coordinates and edges connect neighboring (distance = 1) points. This
    // simple algorithm would have worked if integer points were always connected in every possible
    // supercell. However, this is not the case (consider basis {(3, 1, 0), (5, 1, 0), (0, 0, 1)}).
    // Fortunately, there is an easy fix: if one uses a quotient of the full (infinite) primary
    // lattice graph under equivalence classes induced by superlattice, then every equivalence class
    // (and, therefore, every integer fractionary point) will be reached. In practice, this means
    // that if step out of the supercell when considering an edge, we apply a linear combination of
    /// supercell basis vectors to the point to get back into the supercell.
    for (size_t i = 0; i < result.size(); ++i) {
        auto site = result[i];

        for (int direction = 0; direction < 3; ++direction) {
            for (int sign : { -1, 1 }) {
                Vector3i neighbor = site + supercell_basis_inverse.col(direction) * sign;

                for (int component = 0; component < 3; ++component) {
                    int& coordinate = neighbor[component];
                    coordinate %= det;
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
        primitive_cell = lattice.supercell_fractionary_basis * primitive_cell / det;
    }

    // Order primitive cells in the same way as in the legacy code.
    std::ranges::sort(result, [](Vector3i const& a, Vector3i const& b) {
        return std::tie(a[2], a[1], a[0]) < std::tie(b[2], b[1], b[0]);
    });

    return result;
}

std::vector<Site> compute_sites(Lattice const& lattice)
{
    auto translate = [&](Site const& site, Vector3i const& translation) {
        Vector3d fractionary_coordinates = site.fractionary_position + translation.cast<f64>();
        Vector3d cartesian_coordinates = lattice.basis * fractionary_coordinates;
        return Site { site.label, cartesian_coordinates, fractionary_coordinates };
    };

    std::vector<Site> result;
    for (Vector3i const& cell : lattice.primitive_cells_in_supercell) {
        for (Site const& site : lattice.primitive_cell_sites) {
            result.push_back(translate(site, cell));
        }
    }
    return result;
}

parser::DiagnosticOr<void> initialize_lattice(
    parser::DiagnosticEngine& diag,
    Lattice& lattice,
    parser::ParsedFreeformGeometry const& geometry)
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

    // supercell_fractionary_basis
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            int coefficient = geometry.supercell_basis[i][j];
            lattice.supercell_fractionary_basis(j, i) = coefficient;
            if (i >= dims || j >= dims) {
                VERIFY((i == j && coefficient == 1) || (i != j && coefficient == 0));
            }
        }
    }

    if (lattice.supercell_fractionary_basis.determinant() == 0) {
        return diag.error(geometry.supercell_basis_location,
            "supercell does not contain any primary cells");
    }

    // supercell_cartesian_basis
    lattice.supercell_cartesian_basis = lattice.basis * lattice.supercell_fractionary_basis.cast<f64>();

    // primitive_cells_in_supercell
    lattice.primitive_cells_in_supercell = compute_primitive_cells_in_supercell(lattice);

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
    lattice.sites = compute_sites(lattice);

    return {};
}

parser::DiagnosticOr<void> initialize_hamiltonian(
    parser::DiagnosticEngine& diag,
    Hamiltonian& hamiltonian,
    Lattice const& lattice,
    parser::ParsedFreeformGeometry const& geometry,
    parser::ParsedFreeformGeometryParameters const& parameters)
{
    Matrix3d lattice_basis_inverse = lattice.basis.inverse();
    int cell_site_count = lattice.primitive_cell_sites.size();
    int total_site_count = lattice.sites.size();
    int cells_count = lattice.primitive_cells_in_supercell.size();

    hamiltonian.interactions.resize(total_site_count);
    for (auto& site : hamiltonian.interactions) {
        site.mu_up = parameters.mu_up;
        site.mu_down = parameters.mu_down;
    }

    hamiltonian.hoppings[0] = MatrixXd::Zero(total_site_count, total_site_count);
    hamiltonian.hoppings[1] = MatrixXd::Zero(total_site_count, total_site_count);

    auto [supercell_size, supercell_basis_inverse] = compute_premultiplied_inverse(lattice.supercell_fractionary_basis);

    std::map<Vector3i, int, Vector3iComparator> index_of_primitive_cell;
    for (int i = 0; i < cells_count; ++i) {
        index_of_primitive_cell[lattice.primitive_cells_in_supercell[i]] = i;
    }

    auto process_on_site_interaction = [&](parser::ParsedFreeformGeometry::OnSiteInteraction const& interaction) -> parser::DiagnosticOr<void> {
        for (int cell = 0; cell < cells_count; ++cell) {
            int site_index = cell * cell_site_count + interaction.site;
            auto& site = hamiltonian.interactions[site_index];
            site.u += interaction.u;
            site.mu_up -= interaction.mu_up_offset;
            site.mu_down -= interaction.mu_down_offset;
        }
        return {};
    };

    auto process_hopping = [&](parser::ParsedFreeformGeometry::Hopping const& hopping) -> parser::DiagnosticOr<void> {
        Vector3d fractionary_shift_fp = lattice.sites[hopping.to].cartesian_position
            - (lattice.sites[hopping.from].cartesian_position + Vector3d { hopping.coordinate_delta });
        fractionary_shift_fp = lattice_basis_inverse * fractionary_shift_fp;

        // FIXME: We should not allow shifts that have more dimensions than the primary cell lattice
        //        itself. We should also split dimension number of the lattice and of the points
        //        inside it.
        Vector3i fractionary_shift = fractionary_shift_fp.cast<int>();

        double badness = (fractionary_shift.cast<f64>() - fractionary_shift_fp).squaredNorm();

        if (badness > epsilon) {
            return diag.error(hopping.location,
                "distance {:.1e} (in supercell lattice basis) to the best candidate for the "
                "mentioned site is larger than allowed {:.1e} rounding error",
                sqrt(badness), sqrt(epsilon));
        }

        for (int from_cell_index = 0; from_cell_index < cells_count; ++from_cell_index) {
            auto from_cell = lattice.primitive_cells_in_supercell[from_cell_index];
            auto to_cell = from_cell + fractionary_shift;

            Vector3i to_coordinates = supercell_basis_inverse * to_cell;
            bool should_negate_phase = false;
            for (int component = 0; component < 3; ++component) {
                int coordinate = to_coordinates[component];
                int coordinate_inside_supercell = coordinate % supercell_size;
                if (coordinate_inside_supercell < 0) {
                    coordinate_inside_supercell += supercell_size;
                }
                to_coordinates[component] = coordinate_inside_supercell;

                int phase_power = (coordinate - coordinate_inside_supercell) / supercell_size;
                should_negate_phase ^= (phase_power * parameters.should_negate_phase[component]) & 1;
            }

            to_coordinates = lattice.supercell_fractionary_basis * to_coordinates / supercell_size;
            int to_cell_index = index_of_primitive_cell.at(to_coordinates);

            int from_site = from_cell_index * cell_site_count + hopping.from;
            int to_site = to_cell_index * cell_site_count + hopping.to;

            if (from_site == to_site) {
                return diag.error(hopping.location,
                    "hopping terms with equal endpoints in quotient of the lattice by "
                    "supercell are not implemented, shift chemical potential instead");
            }

            int phase_shift = should_negate_phase ? -1 : 1;

            hamiltonian.hoppings[0](from_site, to_site) += hopping.t_up * phase_shift;
            hamiltonian.hoppings[0](to_site, from_site) += hopping.t_up * phase_shift;
            hamiltonian.hoppings[1](to_site, from_site) += hopping.t_down * phase_shift;
            hamiltonian.hoppings[1](from_site, to_site) += hopping.t_down * phase_shift;
        }

        return {};
    };

    for (auto const& term : geometry.hamiltonian) {
        TRY(visit(term, process_hopping, process_on_site_interaction));
    }
    return {};
}

} // namespace

parser::DiagnosticOr<FreeformGeometry> FreeformGeometry::create(
    parser::DiagnosticEngine& diag,
    parser::ParsedFreeformGeometry const& geometry,
    parser::ParsedFreeformGeometryParameters const& parameters)
{
    FreeformGeometry result;
    TRY(result.initialize(diag, geometry, parameters));
    return result;
}

void FreeformGeometry::legacy_compatible_format_into(std::ostream& stream) const
{
    std::print(stream, " ================================================================\n"
                       " Basic real space geometry info\n"
                       "\n"
                       " Crystal atomic basis\n");
    for (int i = 0; i < m_lattice.primitive_cell_sites.size(); ++i) {
        auto const& coords = m_lattice.primitive_cell_sites[i].fractionary_position;
        std::println(stream, "{:3}{:14.7f}{:14.7f}{:14.7f}", i, coords[0], coords[1], coords[2]);
    }

    std::print(stream, "\n"
                       " Basis cell vectors\n");
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            std::print(stream, "{:14.7f}", m_lattice.basis(j, i));
        }
        std::print(stream, "\n");
    }

    std::print(stream, "\n"
                       "\n"
                       " Supercell vectors (fractionary unit)\n");
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            std::print(stream, "{:5}", m_lattice.supercell_fractionary_basis(j, i));
        }
        std::print(stream, "\n");
    }

    std::print(stream, "\n"
                       " Super-Lattice vectors (cartesian)\n");
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            std::print(stream, "{:14.7f}", m_lattice.supercell_cartesian_basis(j, i));
        }
        std::print(stream, "\n");
    }

    std::print(stream, "\n"
                       " ================================================================\n"
                       " Real space lattice\n"
                       "\n");
    std::print(stream, " Number of orbitals in primitive cell: {:12}\n", m_lattice.primitive_cell_sites.size());
    std::print(stream, " Total number of orbitals:             {:12}\n", m_lattice.sites.size());
    std::print(stream, " index  label   type       X           Y         Z   \n");
    for (int i = 0; i < m_lattice.sites.size(); ++i) {
        auto const& site = m_lattice.sites[i];
        auto const& coords = site.cartesian_position;
        std::print(stream, "{:3} {:3} {:3}{:14.5f}{:14.5f}{:14.5f}\n",
            i, site.label, i % m_lattice.primitive_cell_sites.size(), coords[0], coords[1], coords[2]);
    }
}

parser::DiagnosticOr<void> FreeformGeometry::initialize(
    parser::DiagnosticEngine& diag,
    parser::ParsedFreeformGeometry const& geometry,
    parser::ParsedFreeformGeometryParameters const& parameters)
{
    TRY(initialize_lattice(diag, m_lattice, geometry));
    TRY(initialize_hamiltonian(diag, m_hamiltonian, m_lattice, geometry, parameters));
    return {};
}

} // namespace dqmc::sema
