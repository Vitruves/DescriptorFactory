#pragma once

#include <Eigen/Dense>
#include <cmath>

namespace desfact {
namespace descriptors {
namespace detail {

// Helper function implementations
inline bool isMatrixSizeValid(const Eigen::MatrixXd& matrix, unsigned long maxDim = 2000) {
    return matrix.rows() > 0 && matrix.cols() > 0 && 
           static_cast<unsigned long>(matrix.rows()) < maxDim && 
           static_cast<unsigned long>(matrix.cols()) < maxDim;
}

inline bool hasInvalidValues(const Eigen::MatrixXd& matrix) {
    for (int i = 0; i < matrix.rows(); ++i) {
        for (int j = 0; j < matrix.cols(); ++j) {
            if (!std::isfinite(matrix(i, j))) {
                return true;
            }
        }
    }
    return false;
}

inline Eigen::MatrixXd safeMatrixMultiply(const Eigen::MatrixXd& a, const Eigen::MatrixXd& b) {
    if (!isMatrixSizeValid(a) || !isMatrixSizeValid(b) || a.cols() != b.rows()) {
        return Eigen::MatrixXd(0, 0);
    }
    
    Eigen::MatrixXd result = a * b;
    
    if (hasInvalidValues(result)) {
        return Eigen::MatrixXd(0, 0);
    }
    
    return result;
}

} // namespace detail
} // namespace descriptors
} // namespace desfact 