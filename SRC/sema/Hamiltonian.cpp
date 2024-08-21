#include <Eigen/LU>

#include "SRC/parser/FreeformGeometryParametersParser.h"
#include "SRC/parser/FreeformGeometryParser.h"
#include "SRC/sema/FreeformGeometry.h"

namespace dqmc::sema {

namespace {

using RawGeometry = parser::ParsedFreeformGeometry;

struct HamiltonianBuildingContext : Context {
    HamiltonianBuildingContext(Context& ctx)
        : Context(ctx)
    {
    }

    void add_on_site_interaction(RawGeometry::OnSiteInteraction const& interaction);
    parser::DiagnosticOr<void> add_hopping(RawGeometry::Hopping const& hopping);
    parser::DiagnosticOr<void> build_hamiltonian();

    Matrix3d lattice_basis_inverse = lattice.basis.inverse();
    int cell_site_count = lattice.primitive_cell_sites.size();
    int total_site_count = lattice.sites.size();
    int cells_count = lattice.primitive_cells_in_supercell.size();
};

// Adds an interaction term (on-site U and chemical potential shifts) to the hamiltonian.
void HamiltonianBuildingContext::add_on_site_interaction(RawGeometry::OnSiteInteraction const& interaction)
{
    for (int cell = 0; cell < cells_count; ++cell) {
        int site_index = cell * cell_site_count + interaction.site;
        auto& site = hamiltonian.interactions[site_index];
        site.u += interaction.u;
        site.mu_up -= interaction.mu_up_offset;
        site.mu_down -= interaction.mu_down_offset;
    }
}

// Adds a hopping term to the hamiltonian.
parser::DiagnosticOr<void> HamiltonianBuildingContext::add_hopping(RawGeometry::Hopping const& hopping)
{
    // First, we figure out the sum of lattice basis vectors (modulo superlattice) that gives the
    // delta specified.
    Vector3d fractional_shift_fp = lattice.sites[hopping.to].cartesian_position
        - (lattice.sites[hopping.from].cartesian_position + Vector3d { hopping.coordinate_delta });
    fractional_shift_fp = lattice_basis_inverse * fractional_shift_fp;

    // FIXME: We should not allow shifts that have more dimensions than the primary cell lattice
    //        itself. We should also split dimension number of the lattice and of the points
    //        inside it.
    Vector3i fractional_shift = fractional_shift_fp.cast<int>();

    double badness = (fractional_shift.cast<f64>() - fractional_shift_fp).squaredNorm();

    if (badness > epsilon) {
        return diag.error(hopping.location,
            "distance {:.1e} (in supercell lattice basis) to the best candidate for the "
            "mentioned site is larger than allowed {:.1e} rounding error",
            sqrt(badness), sqrt(epsilon));
    }

    for (int from_cell_index = 0; from_cell_index < cells_count; ++from_cell_index) {
        // Next, for every primitive cell in a supercell, we figure out the primitive cell
        // to which fractional_shift would lead us.
        auto from_cell = lattice.primitive_cells_in_supercell[from_cell_index];
        auto to_cell = from_cell + fractional_shift;

        auto [to_cell_index, supercell] = lattice.fractional_coords_to_cell(to_cell);
        bool should_negate_phase = false;
        for (int i = 0; i < 3; ++i) {
            int phase_power = supercell[i];
            should_negate_phase ^= (phase_power * parameters.should_negate_phase[i]) & 1;
        }

        // Lastly, we compute site indices and add hoppings to the hamiltonian.
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
}

// Iterates over parsed entries specified in "#HAMILT" and calculates hopping integrals, on-site
// interaction terms, and chemical potentials.
parser::DiagnosticOr<void> HamiltonianBuildingContext::build_hamiltonian()
{
    hamiltonian.interactions.resize(total_site_count);
    for (auto& site : hamiltonian.interactions) {
        site.mu_up = parameters.mu_up;
        site.mu_down = parameters.mu_down;
    }

    hamiltonian.hoppings[0] = MatrixXd::Zero(total_site_count, total_site_count);
    hamiltonian.hoppings[1] = MatrixXd::Zero(total_site_count, total_site_count);

    for (auto const& term : geometry.hamiltonian) {
        TRY(visit(
            term,
            [&](RawGeometry::OnSiteInteraction const& interaction) {
                add_on_site_interaction(interaction);
                return parser::DiagnosticOr<void> {};
            },
            [&](RawGeometry::Hopping const& hopping) {
                return add_hopping(hopping);
            }));
    }
    return {};
}

} // namespace

parser::DiagnosticOr<void> build_hamiltonian(Context& ctx)
{
    return HamiltonianBuildingContext { ctx }.build_hamiltonian();
}

} // namespace dqmc::sema
