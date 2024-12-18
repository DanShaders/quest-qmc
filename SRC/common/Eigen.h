#include <Eigen/Core>

namespace dqmc {

using Eigen::Matrix3d, Eigen::Matrix3i, Eigen::MatrixXd, Eigen::MatrixXi,
    Eigen::Vector3d, Eigen::Vector3i;

using Vector3iComparator = decltype([](Vector3i const& a, Vector3i const& b) {
    return std::tie(a[0], a[1], a[2]) < std::tie(b[0], b[1], b[2]);
});

// Returns { abs(det(basis)), (basis)^-1 * abs(det) }.
std::pair<int, Matrix3i> compute_adjugate(Matrix3i const& basis);

inline Vector3i round_to_nearest_integer(Vector3d const& v)
{
    return v.array().round().cast<int>();
}

} // namespace dqmc
