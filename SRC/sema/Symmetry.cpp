#include <Eigen/Geometry>
#include <numbers>
#include <ranges>
#include <set>

#include "SRC/common/Enumerate.h"
#include "SRC/parser/FreeformGeometryParametersParser.h"
#include "SRC/parser/FreeformGeometryParser.h"
#include "SRC/sema/FreeformGeometry.h"

namespace {

// Used for pretty-printing sites in diagnostic messages.
struct FormattableCoordinates {
    Eigen::Vector3d coords;
};

struct FormattableSite {
    dqmc::sema::Context ctx;
    int site_index;
};

} // namespace

template<>
struct std::formatter<FormattableCoordinates> : std::formatter<std::string> {
    auto format(FormattableCoordinates const& wrapper, auto& ctx) const
    {
        return formatter<string>::format(std::format("({}, {}, {})", wrapper.coords.x(), wrapper.coords.y(), wrapper.coords.z()), ctx);
    }
};

template<>
struct std::formatter<FormattableSite> : std::formatter<std::string> {
    auto format(FormattableSite const& site, auto& ctx) const
    {
        int index = site.site_index % site.ctx.lattice.primitive_cell_sites.size();
        auto const& coords = site.ctx.lattice.sites[site.site_index].cartesian_position;

        return formatter<string>::format(std::format("site {} at {}", index, FormattableCoordinates { coords }), ctx);
    }
};

namespace dqmc::sema {

namespace {

using RawGeometry = parser::ParsedFreeformGeometry;
using Transform = Eigen::Transform<f64, 3, Eigen::TransformTraits::Affine>;

struct SymmetryProcessingContext : Context {
    SymmetryProcessingContext(Context& ctx)
        : Context(ctx)
    {
    }

    FormattableSite as_site(int index)
    {
        return { *this, index };
    }

    parser::DiagnosticOr<int> find_mapped_site(int site);
    parser::DiagnosticOr<std::vector<int>> compute_site_permutation();
    void find_equivalence_classes();

    // Current symmetry source range
    parser::SourceRange location;
    // Symmetry operation (in Cartesian coordinates) and switch to the fractional coordinates
    Transform action;
};

// Constructs a Transform that applies a given symmetry to a point. Transform does not change basis.
Transform construct_transform(RawGeometry::SymmetryAction const& action)
{
    using Eigen::Matrix4d;

    return Transform { visit(
        action,
        [&](RawGeometry::Rotation const& rotation) -> Matrix4d {
            Vector3d o { rotation.center[0], rotation.center[1], rotation.center[2] };
            // (-1)^(rotation.is_rotoreflection) * \hat{R} (P - O) + O
            //                                       (2)     (1)  (3)

            Matrix4d operator1 {
                { 1, 0, 0, -o.x() },
                { 0, 1, 0, -o.y() },
                { 0, 0, 1, -o.z() },
                { 0, 0, 0, 1 },
            };

            // See https://en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle
            Vector3d u { rotation.direction[0], rotation.direction[1], rotation.direction[2] };
            u.normalize();
            Matrix3d outer_product = u * u.transpose();
            Matrix3d cross_product_matrix {
                { 0, -u.z(), u.y() },
                { u.z(), 0, -u.x() },
                { -u.y(), u.x(), 0 },
            };
            double angle = 2 * std::numbers::pi_v<f64> / rotation.order;
            Matrix3d rhat
                = std::cos(angle) * Matrix3d::Identity()
                + (1 - std::cos(angle)) * outer_product
                + std::sin(angle) * cross_product_matrix;

            Matrix4d operator2 = Matrix4d::Identity();
            operator2.block<3, 3>(0, 0) = rhat;

            double coefficient = rotation.is_rotoreflection ? -1 : 1;
            Matrix4d operator3 {
                { coefficient, 0, 0, o.x() },
                { 0, coefficient, 0, o.y() },
                { 0, 0, coefficient, o.z() },
                { 0, 0, 0, 1 },
            };

            return operator3 * operator2 * operator1;
        },
        [&](RawGeometry::Reflection const& reflection) -> Matrix4d {
            Vector3d a { reflection.point[0], reflection.point[1], reflection.point[2] };
            Vector3d n { reflection.normal[0], reflection.normal[1], reflection.normal[2] };
            // P - 2 * (n * ((P - A) * n)) / abs(n)^2
            //           (3)   (1)  (2)

            Matrix4d operator1 {
                { 1, 0, 0, -a.x() },
                { 0, 1, 0, -a.y() },
                { 0, 0, 1, -a.z() },
                { 0, 0, 0, 1 },
            };

            Eigen::Matrix<double, 2, 4> operator2 {
                { n.x(), n.y(), n.z(), 0 },
                { 0, 0, 0, 1 },
            };

            double coefficient = 2 / n.squaredNorm();
            Eigen::Matrix<double, 4, 2> operator3 {
                { coefficient * n.x(), 0 },
                { coefficient * n.y(), 0 },
                { coefficient * n.z(), 0 },
                { 0, 1 },
            };

            return Matrix4d::Identity() - operator3 * operator2 * operator1;
        },
        [&](RawGeometry::Inversion const& inversion) {
            // 2 * inversion.center - P
            return Matrix4d {
                { -1, 0, 0, 2 * inversion.center[0] },
                { 0, -1, 0, 2 * inversion.center[1] },
                { 0, 0, -1, 2 * inversion.center[2] },
                { 0, 0, 0, 1 },
            };
        }) };
}

constexpr std::format_string<> warn_lattice_not_symmetric = "lattice is not symmetric under specified space transformation";

// Applies the current symmetry to the site with index `site_index` and finds a site that the
// symmetry maps the given point to. Only sites with the same label as the site `site_index` are
// considered. For function to succeed, best candidate for the mapped site must be close enough to
// the transformed point and be far enough from the second best candidate.
parser::DiagnosticOr<int> SymmetryProcessingContext::find_mapped_site(int site_index)
{
    auto const& site = lattice.sites[site_index];
    auto const& site_coords = site.cartesian_position;

    Vector3d map_fractionary = action * site_coords;

    std::vector<std::pair<f64, size_t>> candidates;

    for (auto [i, candidate] : enumerate(lattice.primitive_cell_sites)) {
        if (candidate.label != site.label) {
            continue;
        }

        Vector3d shift = candidate.fractional_position - map_fractionary;

        candidates.push_back({ (shift - shift.cast<int>().cast<f64>()).squaredNorm(), i });
    }
    if (candidates.size() > 1) {
        std::ranges::nth_element(candidates, candidates.begin() + 1);
    }

    int primitive_index = site_index % lattice.primitive_cell_sites.size();
    auto site_location = geometry.primitive_cell_sites[primitive_index].location;
    Vector3d map_cartesian = lattice.basis * map_fractionary;

    if (candidates[0].first > epsilon) {
        diag.warning(location, warn_lattice_not_symmetric);
        return diag.note(site_location,
            "because {} is mapped to point {} which is too far away "
            "(badness={:.1e} > tolerance={:.1e}) from the nearest site with the same label",
            as_site(site_index), FormattableCoordinates { map_cartesian },
            sqrt(candidates[0].first), sqrt(epsilon));
    }

    if (candidates.size() > 1 && candidates[1].first <= candidates[0].first * 100) {
        diag.warning(location, "symmetry cannot be transformed to site permutation");
        return diag.note(site_location,
            "because {} is mapped to point {} which can be snapped "
            "to either site {} or {}, disambiguate sites using label field",
            as_site(site_index), FormattableCoordinates { map_cartesian },
            candidates[0].second, candidates[1].second);
    }

    int mapped_site_remainder = candidates[0].second;
    auto const& best_candidate = lattice.primitive_cell_sites[mapped_site_remainder];
    Vector3i shift = (map_fractionary - best_candidate.fractional_position).cast<int>();
    int mapped_site_quotient = lattice.fractional_coords_to_cell(shift).primitive_cell;

    return mapped_site_remainder + mapped_site_quotient * lattice.primitive_cell_sites.size();
}

// Verifies that the given integer vector is a permutation.
void verify_is_permutation(std::vector<int> const& permutation)
{
    auto indices = permutation | std::ranges::to<std::set<int>>();
    VERIFY(indices.size() == permutation.size()
        && *indices.begin() == 0 && *indices.rbegin() == permutation.size() - 1);
}

// Computes the mapping between original sites and sites mapped by the current symmetry. Checks that
// Hamiltonian stays the same after symmetry application.
parser::DiagnosticOr<std::vector<int>> SymmetryProcessingContext::compute_site_permutation()
{
    std::vector<int> permutation;
    for (int i = 0; i < lattice.sites.size(); ++i) {
        permutation.push_back(TRY(find_mapped_site(i)));
    }
    verify_is_permutation(permutation);

    for (auto [from, to] : enumerate(permutation)) {
        if (hamiltonian.interactions[from] != hamiltonian.interactions[to]) {
            diag.warning(location, warn_lattice_not_symmetric);
            return diag.note(location,
                "because on-site interaction term of {} is not equivalent to one of {}",
                as_site(from), as_site(to));
        }
    }

    constexpr std::array<std::string_view, 2> spin_name = { "up", "down" };

    for (auto [orig_from, mapped_from] : enumerate(permutation)) {
        for (auto [orig_to, mapped_to] : enumerate(permutation)) {
            for (int spin = 0; spin < 2; ++spin) {
                if (hamiltonian.hoppings[spin](orig_from, orig_to) != hamiltonian.hoppings[spin](mapped_from, mapped_to)) {
                    diag.warning(location, warn_lattice_not_symmetric);
                    return diag.note(location,
                        "because hopping integral from {} to {} for spin-{} electron is not equal "
                        "to one from {} to {}",
                        as_site(orig_from), as_site(orig_to), spin_name[spin],
                        as_site(mapped_from), as_site(mapped_to));
                }
            }
        }
    }

    return permutation;
}

// During construction of equivalence classes, we store them in a disjoint set union which allows us
// to easily and efficiently merge classes found equivalent. To support arbitrary intra-class phase
// shifts, every edge in the structure stores a relative phase difference between nodes it connects.
// Additionally, each tree can be marked inconsistent. This happens if for some pair of sites
// (i, j), we deduce (either directly or not) conflicting phase shifts. The latter means that value
// of a function equivalence class represents (in practice, Green's function value) must be 0.
class DisjointSetUnion {
public:
    struct RootLookupResult {
        int root;
        int phase = 0;
        bool is_inconsistent = false;

        bool operator==(RootLookupResult const&) const = default;

        RootLookupResult contract(RootLookupResult other)
        {
            root = other.root;
            phase ^= other.phase;
            is_inconsistent |= other.is_inconsistent;
            return *this;
        }
    };

    DisjointSetUnion(int n)
        : parent(n)
        , rank(n)
    {
        for (int i = 0; i < n; ++i) {
            parent[i] = { i };
        }
    }

    RootLookupResult get(int i) const
    {
        if (parent[i].root == i) {
            return parent[i];
        }
        return parent[i].contract(get(parent[i].root));
    }

    void merge(int a, int b, int relative_phase)
    {
        auto [a_root, a_phase, a_is_inconsistent] = get(a);
        auto [b_root, b_phase, b_is_inconsistent] = get(b);
        a = a_root;
        b = b_root;
        relative_phase ^= a_phase ^ b_phase;

        if (a == b) {
            if (relative_phase != 0) {
                parent[a].is_inconsistent = true;
            }
            return;
        }
        if (rank[a] == rank[b]) {
            ++rank[a];
        }
        if (rank[a] < rank[b]) {
            std::swap(a, b);
            std::swap(a_is_inconsistent, b_is_inconsistent);
        }
        VERIFY((parent[b] == RootLookupResult { b, 0, b_is_inconsistent }));
        parent[b] = { a, relative_phase, a_is_inconsistent || b_is_inconsistent };
        parent[a].is_inconsistent = parent[b].is_inconsistent;
    }

private:
    mutable std::vector<RootLookupResult> parent;
    std::vector<int> rank;
};

// Finds equivalence classes of pairs of sites based on the symmetries provided in the input and
// some (pseudo)symmetries that are always present.
void SymmetryProcessingContext::find_equivalence_classes()
{
    int sites_count = lattice.sites.size();
    int primitive_sites_count = lattice.primitive_cell_sites.size();

    DisjointSetUnion dsu(sites_count * sites_count);

    auto as_dsu_index = [&](int from, int to) {
        return from * lattice.sites.size() + to;
    };

    // First, add symmetries provided in '#SYMM' section.
    for (auto const& symmetry : geometry.symmetries) {
        location = symmetry.location;
        action = lattice.basis.inverse() * construct_transform(symmetry.action);
        auto permutation_or_error = compute_site_permutation();
        if (!permutation_or_error.has_value()) {
            continue;
        }

        auto permutation = permutation_or_error.value();
        for (auto [orig_from, mapped_from] : enumerate(permutation)) {
            for (auto [orig_to, mapped_to] : enumerate(permutation)) {
                dsu.merge(as_dsu_index(orig_from, orig_to), as_dsu_index(mapped_from, mapped_to), 0);
            }
        }
    }

    // Then, acknowledge that we would have G_ij == G_ji because hopping integrals are real. (In
    // general, G_ij^* == G_ji.)
    for (int i = 0; i < sites_count; ++i) {
        for (int j = i + 1; j < sites_count; ++j) {
            dsu.merge(as_dsu_index(i, j), as_dsu_index(j, i), 0);
        }
    }

    // Apply translations by lattice basis vectors. Strictly speaking, these might not be
    // symmetries because of twisted boundary conditions. However, in that case, it is possible to
    // distribute boundary phase shift without disrupting Hamiltonian much. See section II.A of
    // https://arxiv.org/pdf/2211.07494 for details. H_TPG is trivially symmetric under lattice
    // basis vector translation and it turns out that
    // G_jl_TPG = exp(1j * \vec{k} * (\vec{j} - \vec{l})) * G_kl_TBC with suitably defined
    // coordinates of sites and wave-vector \vec{k}.
    Vector3i k {
        parameters.should_negate_phase[0],
        parameters.should_negate_phase[1],
        parameters.should_negate_phase[2]
    };
    k = lattice.supercell_basis_adjugate.transpose() * k;

    for (auto translation : std::initializer_list<Vector3i> { { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } }) {
        std::vector<int> permutation;
        for (auto const& cell : lattice.primitive_cells_in_supercell) {
            permutation.push_back(lattice.fractional_coords_to_cell(cell + translation).primitive_cell);
        }
        verify_is_permutation(permutation);

        auto compute_phase_shift = [&](int from_cell, int to_cell) {
            auto from = lattice.primitive_cells_in_supercell[from_cell];
            auto to = lattice.primitive_cells_in_supercell[to_cell];
            return k.dot(from - to);
        };

        for (auto [orig_from, mapped_from] : enumerate(permutation)) {
            for (auto [orig_to, mapped_to] : enumerate(permutation)) {
                int scaled_shift = compute_phase_shift(orig_from, orig_to) - compute_phase_shift(mapped_from, mapped_to);
                VERIFY(scaled_shift % lattice.supercell_size == 0);

                int shift = scaled_shift / lattice.supercell_size % 2 == 0 ? 0 : 1;

                for (int from_site = 0; from_site < primitive_sites_count; ++from_site) {
                    for (int to_site = 0; to_site < primitive_sites_count; ++to_site) {
                        int orig_from_site = orig_from * primitive_sites_count + from_site;
                        int orig_to_site = orig_to * primitive_sites_count + to_site;
                        int mapped_from_site = mapped_from * primitive_sites_count + from_site;
                        int mapped_to_site = mapped_to * primitive_sites_count + to_site;
                        dsu.merge(
                            as_dsu_index(orig_from_site, orig_to_site),
                            as_dsu_index(mapped_from_site, mapped_to_site),
                            shift);
                    }
                }
            }
        }
    }

    // Lastly, iterate over pairs of sites, find their final equivalence classes and remap indices
    // to be a bit more sensible in debug output.
    std::map<int, int> pair_classes_remap;

    equivalence_classes.pair_class = MatrixXi::Zero(sites_count, sites_count);
    equivalence_classes.greens_fn_phase = MatrixXi::Zero(sites_count, sites_count);

    for (int i = 0; i < sites_count; ++i) {
        for (int j = 0; j < sites_count; ++j) {
            auto [pair_class, phase, is_inconsistent] = dsu.get(i * sites_count + j);
            int current_class;
            if (auto it = pair_classes_remap.find(pair_class); it != pair_classes_remap.end()) {
                current_class = it->second;
            } else {
                current_class = pair_classes_remap.size();
                pair_classes_remap[pair_class] = current_class;
            }

            if (!is_inconsistent) {
                phase = phase == 0 ? 1 : -1;
            } else {
                phase = 0;
            }

            equivalence_classes.pair_class(i, j) = current_class;
            equivalence_classes.greens_fn_phase(i, j) = phase;
        }
    }
}

} // namespace

parser::DiagnosticOr<void> find_equivalence_classes(Context& ctx)
{
    SymmetryProcessingContext { ctx }.find_equivalence_classes();
    return {};
}

} // namespace dqmc::sema
