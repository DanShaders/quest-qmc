#include <Eigen/Core>

namespace dqmc {

using Eigen::Matrix3d, Eigen::Matrix3i, Eigen::MatrixXd,
    Eigen::Vector3d, Eigen::Vector3i;

using Vector3iComparator = decltype([](Vector3i const& a, Vector3i const& b) {
    return std::tie(a[0], a[1], a[2]) < std::tie(b[0], b[1], b[2]);
});

// Returns { abs(det(basis)), (basis)^-1 * abs(det) }.
std::pair<int, Matrix3i> compute_premultiplied_inverse(Matrix3i const& basis);

} // namespace dqmc
