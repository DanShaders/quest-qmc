#include <Eigen/LU>

#include "SRC/common/Eigen.h"

namespace dqmc {

std::pair<int, Matrix3i> compute_adjugate(Matrix3i const& basis)
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

} // namespace dqmc
