#include "descriptors/eigen.hpp"
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <GraphMol/GraphMol.h>
#include <GraphMol/BondIterators.h>
#include <GraphMol/AtomIterators.h>
#include <cmath>
#include "utils.hpp"

namespace desfact {
namespace descriptors {

// Initialize thread_local cache variables
thread_local std::unordered_map<const RDKit::ROMol*, EigenDescriptor::MatrixCacheEntry> EigenDescriptor::adjacencyCache;
thread_local std::unordered_map<const RDKit::ROMol*, EigenDescriptor::MatrixCacheEntry> EigenDescriptor::degreeCache;
thread_local std::unordered_map<const RDKit::ROMol*, EigenDescriptor::MatrixCacheEntry> EigenDescriptor::laplacianCache;
thread_local std::unordered_map<const RDKit::ROMol*, EigenDescriptor::MatrixCacheEntry> EigenDescriptor::normalizedLaplacianCache;
thread_local std::unordered_map<const RDKit::ROMol*, EigenDescriptor::MatrixCacheEntry> EigenDescriptor::signlessLaplacianCache;
thread_local std::unordered_map<const RDKit::ROMol*, EigenDescriptor::MatrixCacheEntry> EigenDescriptor::weightedAdjacencyCache;
thread_local std::unordered_map<const RDKit::ROMol*, EigenDescriptor::MatrixCacheEntry> EigenDescriptor::distanceMatrixCache;
thread_local std::unordered_map<Eigen::MatrixXd, EigenDescriptor::EigenvalueCacheEntry, EigenDescriptor::MatrixHash, EigenDescriptor::MatrixEqual> EigenDescriptor::eigenvalueCache;

// EigenDescriptor base class implementation
EigenDescriptor::EigenDescriptor(const std::string& name, const std::string& description)
    : Descriptor(name, description) {}
    
// Static method to clear all caches
void EigenDescriptor::clearAllCaches() {
    adjacencyCache.clear();
    degreeCache.clear();
    laplacianCache.clear();
    normalizedLaplacianCache.clear();
    signlessLaplacianCache.clear();
    weightedAdjacencyCache.clear();
    distanceMatrixCache.clear();
    eigenvalueCache.clear();
}

Eigen::MatrixXd EigenDescriptor::buildAdjacencyMatrix(const RDKit::ROMol* mol) const {
    int numAtoms = mol->getNumAtoms();
    Eigen::MatrixXd adjacency = Eigen::MatrixXd::Zero(numAtoms, numAtoms);
    
    for (const RDKit::Bond* bond : mol->bonds()) {
        int startIdx = bond->getBeginAtomIdx();
        int endIdx = bond->getEndAtomIdx();
        adjacency(startIdx, endIdx) = 1.0;
        adjacency(endIdx, startIdx) = 1.0;  // Undirected graph
    }
    
    return adjacency;
}

Eigen::MatrixXd EigenDescriptor::buildDegreeMatrix(const Eigen::MatrixXd& adjacency) const {
    int size = adjacency.rows();
    Eigen::MatrixXd degree = Eigen::MatrixXd::Zero(size, size);
    
    // Sum each row to get the degree
    for (int i = 0; i < size; ++i) {
        degree(i, i) = adjacency.row(i).sum();
    }
    
    return degree;
}

Eigen::MatrixXd EigenDescriptor::buildLaplacianMatrix(const Eigen::MatrixXd& degree, const Eigen::MatrixXd& adjacency) const {
    return degree - adjacency;
}

Eigen::MatrixXd EigenDescriptor::buildNormalizedLaplacian(const Eigen::MatrixXd& degree, const Eigen::MatrixXd& adjacency) const {
    int size = adjacency.rows();
    Eigen::MatrixXd normalizedLaplacian = Eigen::MatrixXd::Identity(size, size);
    
    // Compute D^(-1/2)
    Eigen::MatrixXd dInvSqrt = Eigen::MatrixXd::Zero(size, size);
    for (int i = 0; i < size; ++i) {
        if (degree(i, i) > 0) {
            dInvSqrt(i, i) = 1.0 / std::sqrt(degree(i, i));
        }
    }
    
    // L_norm = I - D^(-1/2) * A * D^(-1/2)
    normalizedLaplacian -= dInvSqrt * adjacency * dInvSqrt;
    
    return normalizedLaplacian;
}

Eigen::VectorXd EigenDescriptor::computeEigenvalues(const Eigen::MatrixXd& matrix) const {
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(matrix);
    return solver.eigenvalues();
}

double EigenDescriptor::computeVariance(const Eigen::VectorXd& values) const {
    double mean = values.mean();
    double variance = 0.0;
    
    for (int i = 0; i < values.size(); ++i) {
        double diff = values(i) - mean;
        variance += diff * diff;
    }
    
    return variance / values.size();
}

double EigenDescriptor::computeSkewness(const Eigen::VectorXd& values) const {
    double mean = values.mean();
    double variance = computeVariance(values);
    double stddev = std::sqrt(variance);
    double skewness = 0.0;
    
    if (stddev < 1e-10) return 0.0;
    
    for (int i = 0; i < values.size(); ++i) {
        double diff = values(i) - mean;
        skewness += diff * diff * diff;
    }
    
    return skewness / (values.size() * std::pow(stddev, 3));
}

double EigenDescriptor::computeKurtosis(const Eigen::VectorXd& values) const {
    double mean = values.mean();
    double variance = computeVariance(values);
    double kurtosis = 0.0;
    
    if (variance < 1e-10) return 0.0;
    
    for (int i = 0; i < values.size(); ++i) {
        double diff = values(i) - mean;
        kurtosis += std::pow(diff, 4);
    }
    
    return kurtosis / (values.size() * variance * variance) - 3.0;  // Excess kurtosis
}

// Helper to compute singular values of a matrix
Eigen::VectorXd EigenDescriptor::computeSingularValues(const Eigen::MatrixXd& matrix) const {
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(matrix, Eigen::ComputeThinU | Eigen::ComputeThinV);
    return svd.singularValues();
}

// Helper for building normalized adjacency matrix
Eigen::MatrixXd EigenDescriptor::buildNormalizedAdjacency(const Eigen::MatrixXd& degree, const Eigen::MatrixXd& adjacency) const {
    int size = adjacency.rows();
    
    // D^(-1/2)
    Eigen::MatrixXd dInvSqrt = Eigen::MatrixXd::Zero(size, size);
    for (int i = 0; i < size; ++i) {
        if (degree(i, i) > 0) {
            dInvSqrt(i, i) = 1.0 / std::sqrt(degree(i, i));
        }
    }
    
    // D^(-1/2) * A * D^(-1/2)
    return dInvSqrt * adjacency * dInvSqrt;
}

// Helper for building signless Laplacian matrix: Q = D + A
Eigen::MatrixXd EigenDescriptor::buildSignlessLaplacian(const Eigen::MatrixXd& degree, const Eigen::MatrixXd& adjacency) const {
    return degree + adjacency;
}

// Helper for building weighted adjacency matrix
Eigen::MatrixXd EigenDescriptor::buildWeightedAdjacency(const RDKit::ROMol* mol) const {
    int numAtoms = mol->getNumAtoms();
    Eigen::MatrixXd weightedAdjacency = Eigen::MatrixXd::Zero(numAtoms, numAtoms);
    
    for (const RDKit::Bond* bond : mol->bonds()) {
        int startIdx = bond->getBeginAtomIdx();
        int endIdx = bond->getEndAtomIdx();
        
        // Weight based on bond type: single=1, double=2, triple=3, aromatic=1.5
        double weight = 1.0;
        switch (bond->getBondType()) {
            case RDKit::Bond::SINGLE:
                weight = 1.0;
                break;
            case RDKit::Bond::DOUBLE:
                weight = 2.0;
                break;
            case RDKit::Bond::TRIPLE:
                weight = 3.0;
                break;
            case RDKit::Bond::AROMATIC:
                weight = 1.5;
                break;
            default:
                weight = 1.0;
        }
        
        weightedAdjacency(startIdx, endIdx) = weight;
        weightedAdjacency(endIdx, startIdx) = weight; // Undirected graph
    }
    
    return weightedAdjacency;
}



std::variant<double, int, std::string> SignlessLaplacianZeroEigenvalues::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return 0.0;
    
    const RDKit::ROMol* rdkitMol = mol.getMolecule().get();
    if (!rdkitMol || rdkitMol->getNumAtoms() == 0) return 0.0;
    
    // Get signless Laplacian matrix from cache or build it
    Eigen::MatrixXd adjacency = getAdjacencyMatrix(rdkitMol);
    Eigen::MatrixXd degree = getDegreeMatrix(rdkitMol);
    Eigen::MatrixXd signlessLaplacian = getSignlessLaplacianMatrix(rdkitMol);
    
    // Get eigenvalues from cache or compute them
    Eigen::VectorXd eigenvalues = getEigenvalues(signlessLaplacian);
    
    // Count zero eigenvalues (with tolerance)
    int zeroCount = 0;
    for (int i = 0; i < eigenvalues.size(); ++i) {
        if (std::abs(eigenvalues(i)) < 1e-10) zeroCount++;
    }
    
    return zeroCount;
}

std::variant<double, int, std::string> SignlessLaplacianPositiveEigenvalues::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return 0;
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol || rdkMol->getNumAtoms() == 0) return 0;
    
    // Get signless Laplacian matrix from cache or build it
    Eigen::MatrixXd signlessLaplacian = getSignlessLaplacianMatrix(rdkMol);
    
    // Get eigenvalues from cache or compute them
    Eigen::VectorXd eigenvalues = getEigenvalues(signlessLaplacian);
    
    int count = 0;
    for (int i = 0; i < eigenvalues.size(); ++i) {
        if (eigenvalues(i) > 1e-10) {
            count++;
        }
    }
    
    return count;
}

std::variant<double, int, std::string> SignlessLaplacianNegativeEigenvalues::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return 0;
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol || rdkMol->getNumAtoms() == 0) return 0;
    
    // Get signless Laplacian matrix from cache or build it
    Eigen::MatrixXd signlessLaplacian = getSignlessLaplacianMatrix(rdkMol);
    
    // Get eigenvalues from cache or compute them
    Eigen::VectorXd eigenvalues = getEigenvalues(signlessLaplacian);
    
    int count = 0;
    for (int i = 0; i < eigenvalues.size(); ++i) {
        if (eigenvalues(i) < -1e-10) {
            count++;
        }
    }
    
    return count;
}

std::variant<double, int, std::string> SignlessLaplacianEigenvalueVariance::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return 0.0;
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol || rdkMol->getNumAtoms() == 0) return 0.0;
    
    // Get signless Laplacian matrix from cache or build it
    Eigen::MatrixXd signlessLaplacian = getSignlessLaplacianMatrix(rdkMol);
    
    // Get eigenvalues from cache or compute them
    Eigen::VectorXd eigenvalues = getEigenvalues(signlessLaplacian);
    
    return computeVariance(eigenvalues);
}

std::variant<double, int, std::string> SignlessLaplacianEigenvalueSkewness::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return 0.0;
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol || rdkMol->getNumAtoms() == 0) return 0.0;
    
    // Get signless Laplacian matrix from cache or build it
    Eigen::MatrixXd signlessLaplacian = getSignlessLaplacianMatrix(rdkMol);
    
    // Get eigenvalues from cache or compute them
    Eigen::VectorXd eigenvalues = getEigenvalues(signlessLaplacian);
    
    return computeSkewness(eigenvalues);
}

std::variant<double, int, std::string> WeightedAdjacencySpectralRadius::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return 0.0;
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (rdkMol->getNumAtoms() == 0) return 0.0;
    
    Eigen::MatrixXd weightedAdj = buildWeightedAdjacency(rdkMol);
    Eigen::VectorXd eigenvalues = computeEigenvalues(weightedAdj);
    
    return eigenvalues.maxCoeff();
}

std::variant<double, int, std::string> MeanFirstPassageTime::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return 0.0;
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (rdkMol->getNumAtoms() <= 1) return 0.0;
    
    Eigen::MatrixXd adjacency = buildAdjacencyMatrix(rdkMol);
    Eigen::MatrixXd degree = buildDegreeMatrix(adjacency);
    Eigen::MatrixXd laplacian = buildLaplacianMatrix(degree, adjacency);
    
    // Compute pseudo-inverse (Moore-Penrose)
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(laplacian, Eigen::ComputeFullU | Eigen::ComputeFullV);
    Eigen::VectorXd singularValues = svd.singularValues();
    Eigen::MatrixXd U = svd.matrixU();
    Eigen::MatrixXd V = svd.matrixV();
    
    // Build pseudo-inverse
    Eigen::MatrixXd pinv = Eigen::MatrixXd::Zero(laplacian.rows(), laplacian.cols());
    for (int i = 0; i < singularValues.size(); ++i) {
        if (singularValues(i) > 1e-10) {
            pinv += (1.0 / singularValues(i)) * V.col(i) * U.col(i).transpose();
        }
    }
    
    // Mean first passage time = trace(pinv) / (n-1)
    int n = laplacian.rows();
    return pinv.trace() / (n - 1);
}

std::variant<double, int, std::string> CommuteTimeDistance::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return 0.0;
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (rdkMol->getNumAtoms() <= 1) return 0.0;
    
    Eigen::MatrixXd adjacency = buildAdjacencyMatrix(rdkMol);
    Eigen::MatrixXd degree = buildDegreeMatrix(adjacency);
    Eigen::MatrixXd laplacian = buildLaplacianMatrix(degree, adjacency);
    
    // Compute pseudo-inverse (Moore-Penrose)
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(laplacian, Eigen::ComputeFullU | Eigen::ComputeFullV);
    Eigen::VectorXd singularValues = svd.singularValues();
    Eigen::MatrixXd U = svd.matrixU();
    Eigen::MatrixXd V = svd.matrixV();
    
    // Build pseudo-inverse
    Eigen::MatrixXd pinv = Eigen::MatrixXd::Zero(laplacian.rows(), laplacian.cols());
    for (int i = 0; i < singularValues.size(); ++i) {
        if (singularValues(i) > 1e-10) {
            pinv += (1.0 / singularValues(i)) * V.col(i) * U.col(i).transpose();
        }
    }
    
    // Average commute time distance
    double sum = 0.0;
    int count = 0;
    int n = laplacian.rows();
    
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            sum += (pinv(i, i) + pinv(j, j) - 2 * pinv(i, j));
            count++;
        }
    }
    
    return count > 0 ? sum / count : 0.0;
}

std::variant<double, int, std::string> KirchhoffIndexVariance::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return 0.0;
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (rdkMol->getNumAtoms() <= 1) return 0.0;
    
    Eigen::MatrixXd adjacency = buildAdjacencyMatrix(rdkMol);
    Eigen::MatrixXd degree = buildDegreeMatrix(adjacency);
    Eigen::MatrixXd laplacian = buildLaplacianMatrix(degree, adjacency);
    
    // Calculate Kirchhoff indices for each vertex-deleted subgraph
    std::vector<double> kirchhoffIndices;
    int n = laplacian.rows();
    
    for (int k = 0; k < n; ++k) {
        // Create the subgraph Laplacian by removing row/col k
        Eigen::MatrixXd subLaplacian = Eigen::MatrixXd::Zero(n-1, n-1);
        int row_idx = 0;
        for (int i = 0; i < n; ++i) {
            if (i == k) continue;
            int col_idx = 0;
            for (int j = 0; j < n; ++j) {
                if (j == k) continue;
                subLaplacian(row_idx, col_idx) = laplacian(i, j);
                col_idx++;
            }
            row_idx++;
        }
        
        // Compute eigenvalues of subgraph Laplacian
        Eigen::VectorXd subEigenvalues = computeEigenvalues(subLaplacian);
        
        // Calculate Kirchhoff index as the sum of reciprocals of non-zero eigenvalues
        double kirchhoffIndex = 0.0;
        for (int i = 0; i < subEigenvalues.size(); ++i) {
            if (subEigenvalues(i) > 1e-10) {
                kirchhoffIndex += 1.0 / subEigenvalues(i);
            }
        }
        kirchhoffIndex *= (n-1);
        kirchhoffIndices.push_back(kirchhoffIndex);
    }
    
    // Calculate variance of Kirchhoff indices
    double mean = 0.0;
    for (double k : kirchhoffIndices) {
        mean += k;
    }
    mean /= kirchhoffIndices.size();
    
    double variance = 0.0;
    for (double k : kirchhoffIndices) {
        double diff = k - mean;
        variance += diff * diff;
    }
    variance /= kirchhoffIndices.size();
    
    return variance;
}

std::variant<double, int, std::string> EffectiveGraphResistanceDistribution::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return 0.0;
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (rdkMol->getNumAtoms() <= 1) return 0.0;
    
    Eigen::MatrixXd adjacency = buildAdjacencyMatrix(rdkMol);
    Eigen::MatrixXd degree = buildDegreeMatrix(adjacency);
    Eigen::MatrixXd laplacian = buildLaplacianMatrix(degree, adjacency);
    
    // Compute pseudo-inverse (Moore-Penrose)
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(laplacian, Eigen::ComputeFullU | Eigen::ComputeFullV);
    Eigen::VectorXd singularValues = svd.singularValues();
    Eigen::MatrixXd U = svd.matrixU();
    Eigen::MatrixXd V = svd.matrixV();
    
    // Build pseudo-inverse
    Eigen::MatrixXd pinv = Eigen::MatrixXd::Zero(laplacian.rows(), laplacian.cols());
    for (int i = 0; i < singularValues.size(); ++i) {
        if (singularValues(i) > 1e-10) {
            pinv += (1.0 / singularValues(i)) * V.col(i) * U.col(i).transpose();
        }
    }
    
    // Calculate effective resistances
    std::vector<double> resistances;
    int n = laplacian.rows();
    
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            double resistance = pinv(i, i) + pinv(j, j) - 2 * pinv(i, j);
            resistances.push_back(resistance);
        }
    }
    
    // Calculate standard deviation of resistance distribution as a measure
    double mean = 0.0;
    for (double r : resistances) {
        mean += r;
    }
    mean /= resistances.size();
    
    double variance = 0.0;
    for (double r : resistances) {
        double diff = r - mean;
        variance += diff * diff;
    }
    variance /= resistances.size();
    
    return std::sqrt(variance); // Return standard deviation
}

std::variant<double, int, std::string> LocalClusteringCoefficientDistribution::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return 0.0;
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (rdkMol->getNumAtoms() < 3) return 0.0;
    
    Eigen::MatrixXd adjacency = buildAdjacencyMatrix(rdkMol);
    Eigen::MatrixXd squared = adjacency * adjacency;
    Eigen::MatrixXd cubed = squared * adjacency;
    
    // Calculate local clustering coefficients
    std::vector<double> clusteringCoefficients;
    int n = adjacency.rows();
    
    for (int i = 0; i < n; ++i) {
        double degree = adjacency.row(i).sum();
        if (degree < 2) {
            // For isolated or pendant vertices, clustering coefficient is 0
            clusteringCoefficients.push_back(0.0);
        } else {
            // Number of triangles that include vertex i is (A³)_ii / 2
            double triangles = cubed(i, i) / 2.0;
            // Maximum possible triangles is d_i*(d_i-1)/2
            double maxTriangles = degree * (degree - 1) / 2.0;
            double coef = triangles / maxTriangles;
            clusteringCoefficients.push_back(coef);
        }
    }
    
    // Calculate standard deviation of clustering coefficients
    double mean = 0.0;
    for (double c : clusteringCoefficients) {
        mean += c;
    }
    mean /= clusteringCoefficients.size();
    
    double variance = 0.0;
    for (double c : clusteringCoefficients) {
        double diff = c - mean;
        variance += diff * diff;
    }
    variance /= clusteringCoefficients.size();
    
    return std::sqrt(variance); // Return standard deviation
}

std::variant<double, int, std::string> GraphRobustnessIndex::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return 0.0;
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (rdkMol->getNumAtoms() <= 1) return 0.0;
    
    Eigen::MatrixXd adjacency = buildAdjacencyMatrix(rdkMol);
    Eigen::MatrixXd degree = buildDegreeMatrix(adjacency);
    Eigen::MatrixXd laplacian = buildLaplacianMatrix(degree, adjacency);
    Eigen::VectorXd eigenvalues = computeEigenvalues(laplacian);
    
    // Sort eigenvalues in ascending order
    std::vector<double> sortedEigenvalues(eigenvalues.data(), eigenvalues.data() + eigenvalues.size());
    std::sort(sortedEigenvalues.begin(), sortedEigenvalues.end());
    
    // Second smallest eigenvalue (algebraic connectivity)
    double algebraicConnectivity = sortedEigenvalues.size() > 1 ? sortedEigenvalues[1] : 0.0;
    double maxDegree = degree.diagonal().maxCoeff();
    
    // Robustness index definition based on algebraic connectivity and max degree
    return algebraicConnectivity / maxDegree;
}

std::variant<double, int, std::string> NormalizedEstradaIndex::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return 0.0;
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (rdkMol->getNumAtoms() == 0) return 0.0;
    
    Eigen::MatrixXd adjacency = buildAdjacencyMatrix(rdkMol);
    Eigen::MatrixXd degree = buildDegreeMatrix(adjacency);
    Eigen::MatrixXd normAdjacency = buildNormalizedAdjacency(degree, adjacency);
    Eigen::VectorXd eigenvalues = computeEigenvalues(normAdjacency);
    
    double estradaIndex = 0.0;
    for (int i = 0; i < eigenvalues.size(); ++i) {
        estradaIndex += std::exp(eigenvalues(i));
    }
    
    // Normalize by number of vertices
    return estradaIndex / eigenvalues.size();
}

std::variant<double, int, std::string> GraphBipartivityIndex::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return 0.0;
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (rdkMol->getNumAtoms() == 0) return 0.0;
    
    Eigen::MatrixXd adjacency = buildAdjacencyMatrix(rdkMol);
    Eigen::VectorXd eigenvalues = computeEigenvalues(adjacency);
    
    double sumSquaredEigs = eigenvalues.squaredNorm();
    if (sumSquaredEigs < 1e-10) return 0.0;
    
    // Calculate the sum of squares of only negative eigenvalues
    double sumSquaredNegative = 0.0;
    for (int i = 0; i < eigenvalues.size(); ++i) {
        if (eigenvalues(i) < 0) {
            sumSquaredNegative += eigenvalues(i) * eigenvalues(i);
        }
    }
    
    // Bipartivity index = sum(λ_i^2 where λ_i < 0) / sum(λ_i^2 for all i)
    return sumSquaredNegative / sumSquaredEigs;
}

std::variant<double, int, std::string> SpanningTreeEntropy::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return 0.0;
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (rdkMol->getNumAtoms() <= 1) return 0.0;
    
    Eigen::MatrixXd adjacency = buildAdjacencyMatrix(rdkMol);
    Eigen::MatrixXd degree = buildDegreeMatrix(adjacency);
    Eigen::MatrixXd laplacian = buildLaplacianMatrix(degree, adjacency);
    
    // Get eigenvalues of the Laplacian matrix
    Eigen::VectorXd eigenvalues = computeEigenvalues(laplacian);
    
    // Sort them to find non-zero eigenvalues
    std::vector<double> sortedEigenvalues(eigenvalues.data(), eigenvalues.data() + eigenvalues.size());
    std::sort(sortedEigenvalues.begin(), sortedEigenvalues.end());
    
    // Calculate the geometric mean of non-zero eigenvalues (μ)
    double logProduct = 0.0;
    int count = 0;
    for (double val : sortedEigenvalues) {
        if (val > 1e-10) {
            logProduct += std::log(val);
            count++;
        }
    }
    
    if (count == 0) return 0.0;
    
    double geometricMean = std::exp(logProduct / count);
    
    // Calculate entropy using formula: -Σ(λ_i/(n*μ) * log(λ_i/(n*μ)))
    double entropy = 0.0;
    int n = eigenvalues.size();
    
    for (double val : sortedEigenvalues) {
        if (val > 1e-10) {
            double p = val / (n * geometricMean);
            entropy -= p * std::log2(p);
        }
    }
    
    return entropy;
}


// Adjacency matrix descriptors implementation
AdjacencyNonZeroEntries::AdjacencyNonZeroEntries()
    : EigenDescriptor("AdjNonZero", "Number of non-zero entries in adjacency matrix (edges)") {}

std::variant<double, int, std::string> AdjacencyNonZeroEntries::calculate(const Molecule& mol) const {
    if (!mol.isValid()) {
        return 0;
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    
    return static_cast<int>(rdkMol->getNumBonds());
}

AdjacencyMatrixTrace::AdjacencyMatrixTrace()
    : EigenDescriptor("AdjTrace", "Trace of adjacency matrix") {}

std::variant<double, int, std::string> AdjacencyMatrixTrace::calculate(const Molecule& mol) const {
    if (!mol.isValid()) {
        return 0.0;
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    Eigen::MatrixXd adjacency = buildAdjacencyMatrix(rdkMol);
    
    return adjacency.trace();
}

AdjacencyFrobeniusNorm::AdjacencyFrobeniusNorm()
    : EigenDescriptor("AdjFrobNorm", "Frobenius norm of adjacency matrix") {}

std::variant<double, int, std::string> AdjacencyFrobeniusNorm::calculate(const Molecule& mol) const {
    if (!mol.isValid()) {
        return 0.0;
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    Eigen::MatrixXd adjacency = buildAdjacencyMatrix(rdkMol);
    
    return adjacency.norm();
}

AdjacencySpectralRadius::AdjacencySpectralRadius()
    : EigenDescriptor("AdjSpecRad", "Spectral radius (largest eigenvalue) of adjacency matrix") {}

std::variant<double, int, std::string> AdjacencySpectralRadius::calculate(const Molecule& mol) const {
    if (!mol.isValid()) {
        return 0.0;
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (rdkMol->getNumAtoms() == 0) {
        return 0.0;
    }
    
    Eigen::MatrixXd adjacency = buildAdjacencyMatrix(rdkMol);
    Eigen::VectorXd eigenvalues = computeEigenvalues(adjacency);
    
    return eigenvalues.maxCoeff();
}

AdjacencySmallestEigenvalue::AdjacencySmallestEigenvalue()
    : EigenDescriptor("AdjMinEig", "Smallest eigenvalue of adjacency matrix") {}

std::variant<double, int, std::string> AdjacencySmallestEigenvalue::calculate(const Molecule& mol) const {
    if (!mol.isValid()) {
        return 0.0;
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (rdkMol->getNumAtoms() == 0) {
        return 0.0;
    }
    
    Eigen::MatrixXd adjacency = buildAdjacencyMatrix(rdkMol);
    Eigen::VectorXd eigenvalues = computeEigenvalues(adjacency);
    
    return eigenvalues.minCoeff();
}

AdjacencySumEigenvalues::AdjacencySumEigenvalues()
    : EigenDescriptor("AdjEigSum", "Sum of eigenvalues of adjacency matrix") {}

std::variant<double, int, std::string> AdjacencySumEigenvalues::calculate(const Molecule& mol) const {
    if (!mol.isValid()) {
        return 0.0;
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (rdkMol->getNumAtoms() == 0) {
        return 0.0;
    }
    
    Eigen::MatrixXd adjacency = buildAdjacencyMatrix(rdkMol);
    Eigen::VectorXd eigenvalues = computeEigenvalues(adjacency);
    
    return eigenvalues.sum();
}

AdjacencySumSquaresEigenvalues::AdjacencySumSquaresEigenvalues()
    : EigenDescriptor("AdjEigSumSq", "Sum of squares of eigenvalues of adjacency matrix") {}

std::variant<double, int, std::string> AdjacencySumSquaresEigenvalues::calculate(const Molecule& mol) const {
    if (!mol.isValid()) {
        return 0.0;
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (rdkMol->getNumAtoms() == 0) {
        return 0.0;
    }
    
    Eigen::MatrixXd adjacency = buildAdjacencyMatrix(rdkMol);
    Eigen::VectorXd eigenvalues = computeEigenvalues(adjacency);
    
    return eigenvalues.squaredNorm();
}

AdjacencyEigenvalueVariance::AdjacencyEigenvalueVariance()
    : EigenDescriptor("AdjEigVar", "Variance of eigenvalues of adjacency matrix") {}

std::variant<double, int, std::string> AdjacencyEigenvalueVariance::calculate(const Molecule& mol) const {
    if (!mol.isValid()) {
        return 0.0;
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (rdkMol->getNumAtoms() == 0) {
        return 0.0;
    }
    
    Eigen::MatrixXd adjacency = buildAdjacencyMatrix(rdkMol);
    Eigen::VectorXd eigenvalues = computeEigenvalues(adjacency);
    
    return computeVariance(eigenvalues);
}

AdjacencyEigenvalueSkewness::AdjacencyEigenvalueSkewness()
    : EigenDescriptor("AdjEigSkew", "Skewness of eigenvalues of adjacency matrix") {}

std::variant<double, int, std::string> AdjacencyEigenvalueSkewness::calculate(const Molecule& mol) const {
    if (!mol.isValid()) {
        return 0.0;
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (rdkMol->getNumAtoms() < 3) {
        return 0.0;
    }
    
    Eigen::MatrixXd adjacency = buildAdjacencyMatrix(rdkMol);
    Eigen::VectorXd eigenvalues = computeEigenvalues(adjacency);
    
    return computeSkewness(eigenvalues);
}

AdjacencyEigenvalueKurtosis::AdjacencyEigenvalueKurtosis()
    : EigenDescriptor("AdjEigKurt", "Kurtosis of eigenvalues of adjacency matrix") {}

std::variant<double, int, std::string> AdjacencyEigenvalueKurtosis::calculate(const Molecule& mol) const {
    if (!mol.isValid()) {
        return 0.0;
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (rdkMol->getNumAtoms() < 4) {
        return 0.0;
    }
    
    Eigen::MatrixXd adjacency = buildAdjacencyMatrix(rdkMol);
    Eigen::VectorXd eigenvalues = computeEigenvalues(adjacency);
    
    return computeKurtosis(eigenvalues);
}

AdjacencyPositiveEigenvalues::AdjacencyPositiveEigenvalues()
    : EigenDescriptor("AdjPosEig", "Number of positive eigenvalues of adjacency matrix") {}

std::variant<double, int, std::string> AdjacencyPositiveEigenvalues::calculate(const Molecule& mol) const {
    if (!mol.isValid()) {
        return 0;
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (rdkMol->getNumAtoms() == 0) {
        return 0;
    }
    
    Eigen::MatrixXd adjacency = buildAdjacencyMatrix(rdkMol);
    Eigen::VectorXd eigenvalues = computeEigenvalues(adjacency);
    
    int count = 0;
    for (int i = 0; i < eigenvalues.size(); ++i) {
        if (eigenvalues(i) > 1e-10) {
            count++;
        }
    }
    
    return count;
}

AdjacencyNegativeEigenvalues::AdjacencyNegativeEigenvalues()
    : EigenDescriptor("AdjNegEig", "Number of negative eigenvalues of adjacency matrix") {}

std::variant<double, int, std::string> AdjacencyNegativeEigenvalues::calculate(const Molecule& mol) const {
    if (!mol.isValid()) {
        return 0;
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (rdkMol->getNumAtoms() == 0) {
        return 0;
    }
    
    Eigen::MatrixXd adjacency = buildAdjacencyMatrix(rdkMol);
    Eigen::VectorXd eigenvalues = computeEigenvalues(adjacency);
    
    int count = 0;
    for (int i = 0; i < eigenvalues.size(); ++i) {
        if (eigenvalues(i) < -1e-10) {
            count++;
        }
    }
    
    return count;
}

AdjacencyZeroEigenvalues::AdjacencyZeroEigenvalues()
    : EigenDescriptor("AdjZeroEig", "Number of zero eigenvalues of adjacency matrix") {}

std::variant<double, int, std::string> AdjacencyZeroEigenvalues::calculate(const Molecule& mol) const {
    if (!mol.isValid()) {
        return 0;
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (rdkMol->getNumAtoms() == 0) {
        return 0;
    }
    
    Eigen::MatrixXd adjacency = buildAdjacencyMatrix(rdkMol);
    Eigen::VectorXd eigenvalues = computeEigenvalues(adjacency);
    
    int count = 0;
    for (int i = 0; i < eigenvalues.size(); ++i) {
        if (std::abs(eigenvalues(i)) < 1e-10) {
            count++;
        }
    }
    
    return count;
}

AdjacencyMaxDegree::AdjacencyMaxDegree()
    : EigenDescriptor("AdjMaxDeg", "Maximum degree (max row sum of adjacency matrix)") {}

std::variant<double, int, std::string> AdjacencyMaxDegree::calculate(const Molecule& mol) const {
    if (!mol.isValid()) {
        return 0;
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (rdkMol->getNumAtoms() == 0) {
        return 0;
    }
    
    Eigen::MatrixXd adjacency = buildAdjacencyMatrix(rdkMol);
    
    double maxDegree = 0;
    for (int i = 0; i < adjacency.rows(); ++i) {
        double rowSum = adjacency.row(i).sum();
        maxDegree = std::max(maxDegree, rowSum);
    }
    
    return static_cast<int>(maxDegree);
}

AdjacencyMinDegree::AdjacencyMinDegree()
    : EigenDescriptor("AdjMinDeg", "Minimum degree (min row sum of adjacency matrix)") {}

std::variant<double, int, std::string> AdjacencyMinDegree::calculate(const Molecule& mol) const {
    if (!mol.isValid()) {
        return 0;
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (rdkMol->getNumAtoms() == 0) {
        return 0;
    }
    
    Eigen::MatrixXd adjacency = buildAdjacencyMatrix(rdkMol);
    
    double minDegree = std::numeric_limits<double>::max();
    for (int i = 0; i < adjacency.rows(); ++i) {
        double rowSum = adjacency.row(i).sum();
        if (rowSum > 0) { // Ignore isolated vertices (degree 0)
            minDegree = std::min(minDegree, rowSum);
        }
    }
    
    return static_cast<int>(minDegree == std::numeric_limits<double>::max() ? 0 : minDegree);
}

AdjacencyMeanDegree::AdjacencyMeanDegree()
    : EigenDescriptor("AdjMeanDeg", "Mean degree (average row sum of adjacency matrix)") {}

std::variant<double, int, std::string> AdjacencyMeanDegree::calculate(const Molecule& mol) const {
    if (!mol.isValid()) {
        return 0.0;
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (rdkMol->getNumAtoms() == 0) {
        return 0.0;
    }
    
    Eigen::MatrixXd adjacency = buildAdjacencyMatrix(rdkMol);
    double totalDegree = 0.0;
    
    for (int i = 0; i < adjacency.rows(); ++i) {
        totalDegree += adjacency.row(i).sum();
    }
    
    return totalDegree / adjacency.rows();
}

AdjacencyDegreeVariance::AdjacencyDegreeVariance()
    : EigenDescriptor("AdjDegVar", "Variance of degrees in adjacency matrix") {}

std::variant<double, int, std::string> AdjacencyDegreeVariance::calculate(const Molecule& mol) const {
    if (!mol.isValid()) {
        return 0.0;
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (rdkMol->getNumAtoms() == 0) {
        return 0.0;
    }
    
    Eigen::MatrixXd adjacency = buildAdjacencyMatrix(rdkMol);
    Eigen::VectorXd degrees(adjacency.rows());
    
    for (int i = 0; i < adjacency.rows(); ++i) {
        degrees(i) = adjacency.row(i).sum();
    }
    
    return computeVariance(degrees);
}

AdjacencyDegreeStdDev::AdjacencyDegreeStdDev()
    : EigenDescriptor("AdjDegStd", "Standard deviation of degrees in adjacency matrix") {}

std::variant<double, int, std::string> AdjacencyDegreeStdDev::calculate(const Molecule& mol) const {
    if (!mol.isValid()) {
        return 0.0;
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (rdkMol->getNumAtoms() == 0) {
        return 0.0;
    }
    
    Eigen::MatrixXd adjacency = buildAdjacencyMatrix(rdkMol);
    Eigen::VectorXd degrees(adjacency.rows());
    
    for (int i = 0; i < adjacency.rows(); ++i) {
        degrees(i) = adjacency.row(i).sum();
    }
    
    return std::sqrt(computeVariance(degrees));
}

AdjacencyMatrixRank::AdjacencyMatrixRank()
    : EigenDescriptor("AdjRank", "Rank of adjacency matrix") {}

std::variant<double, int, std::string> AdjacencyMatrixRank::calculate(const Molecule& mol) const {
    if (!mol.isValid()) {
        return 0;
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (rdkMol->getNumAtoms() == 0) {
        return 0;
    }
    
    Eigen::MatrixXd adjacency = buildAdjacencyMatrix(rdkMol);
    Eigen::FullPivLU<Eigen::MatrixXd> lu(adjacency);
    
    return static_cast<int>(lu.rank());
}

GraphEnergy::GraphEnergy()
    : EigenDescriptor("GraphEnergy", "Energy of the graph (sum of absolute eigenvalues)") {}

std::variant<double, int, std::string> GraphEnergy::calculate(const Molecule& mol) const {
    if (!mol.isValid()) {
        return 0.0;
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (rdkMol->getNumAtoms() == 0) {
        return 0.0;
    }
    
    Eigen::MatrixXd adjacency = buildAdjacencyMatrix(rdkMol);
    Eigen::VectorXd eigenvalues = computeEigenvalues(adjacency);
    
    return eigenvalues.cwiseAbs().sum();
}

// Laplacian matrix descriptors implementation
LaplacianSpectralRadius::LaplacianSpectralRadius()
    : EigenDescriptor("LapSpecRad", "Spectral radius of Laplacian matrix") {}

std::variant<double, int, std::string> LaplacianSpectralRadius::calculate(const Molecule& mol) const {
    if (!mol.isValid()) {
        return 0.0;
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (rdkMol->getNumAtoms() == 0) {
        return 0.0;
    }
    
    Eigen::MatrixXd adjacency = buildAdjacencyMatrix(rdkMol);
    Eigen::MatrixXd degree = buildDegreeMatrix(adjacency);
    Eigen::MatrixXd laplacian = buildLaplacianMatrix(degree, adjacency);
    Eigen::VectorXd eigenvalues = computeEigenvalues(laplacian);
    
    return eigenvalues.maxCoeff();
}

LaplacianAlgebraicConnectivity::LaplacianAlgebraicConnectivity()
    : EigenDescriptor("LapAlgConn", "Algebraic connectivity (second-smallest eigenvalue of Laplacian)") {}

std::variant<double, int, std::string> LaplacianAlgebraicConnectivity::calculate(const Molecule& mol) const {
    if (!mol.isValid()) {
        return 0.0;
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (rdkMol->getNumAtoms() <= 1) {
        return 0.0;
    }
    
    Eigen::MatrixXd adjacency = buildAdjacencyMatrix(rdkMol);
    Eigen::MatrixXd degree = buildDegreeMatrix(adjacency);
    Eigen::MatrixXd laplacian = buildLaplacianMatrix(degree, adjacency);
    Eigen::VectorXd eigenvalues = computeEigenvalues(laplacian);
    
    // Sort eigenvalues
    std::vector<double> sortedEigenvalues(eigenvalues.data(), eigenvalues.data() + eigenvalues.size());
    std::sort(sortedEigenvalues.begin(), sortedEigenvalues.end());
    
    // Second smallest eigenvalue (first one is usually close to zero for connected graphs)
    double secondSmallest = (sortedEigenvalues.size() > 1) ? sortedEigenvalues[1] : 0.0;
    
    return secondSmallest;
}

LaplacianZeroEigenvalues::LaplacianZeroEigenvalues()
    : EigenDescriptor("LapZeroEig", "Number of zero eigenvalues in Laplacian (number of components)") {}

std::variant<double, int, std::string> LaplacianZeroEigenvalues::calculate(const Molecule& mol) const {
    if (!mol.isValid()) {
        return 0;
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (rdkMol->getNumAtoms() == 0) {
        return 0;
    }
    
    Eigen::MatrixXd adjacency = buildAdjacencyMatrix(rdkMol);
    Eigen::MatrixXd degree = buildDegreeMatrix(adjacency);
    Eigen::MatrixXd laplacian = buildLaplacianMatrix(degree, adjacency);
    Eigen::VectorXd eigenvalues = computeEigenvalues(laplacian);
    
    int count = 0;
    for (int i = 0; i < eigenvalues.size(); ++i) {
        if (std::abs(eigenvalues(i)) < 1e-10) {
            count++;
        }
    }
    
    return count;
}

LaplacianEnergy::LaplacianEnergy()
    : EigenDescriptor("LapEnergy", "Laplacian energy (sum of squares of eigenvalues)") {}

std::variant<double, int, std::string> LaplacianEnergy::calculate(const Molecule& mol) const {
    if (!mol.isValid()) {
        return 0.0;
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (rdkMol->getNumAtoms() == 0) {
        return 0.0;
    }
    
    Eigen::MatrixXd adjacency = buildAdjacencyMatrix(rdkMol);
    Eigen::MatrixXd degree = buildDegreeMatrix(adjacency);
    Eigen::MatrixXd laplacian = buildLaplacianMatrix(degree, adjacency);
    Eigen::VectorXd eigenvalues = computeEigenvalues(laplacian);
    
    return eigenvalues.squaredNorm();
}

LaplacianMatrixTrace::LaplacianMatrixTrace()
    : EigenDescriptor("LapTrace", "Trace of Laplacian matrix") {}

std::variant<double, int, std::string> LaplacianMatrixTrace::calculate(const Molecule& mol) const {
    if (!mol.isValid()) {
        return 0.0;
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (rdkMol->getNumAtoms() == 0) {
        return 0.0;
    }
    
    Eigen::MatrixXd adjacency = buildAdjacencyMatrix(rdkMol);
    Eigen::MatrixXd degree = buildDegreeMatrix(adjacency);
    Eigen::MatrixXd laplacian = buildLaplacianMatrix(degree, adjacency);
    
    return laplacian.trace();
}

LaplacianMatrixDeterminant::LaplacianMatrixDeterminant()
    : EigenDescriptor("LapDet", "Determinant of Laplacian matrix") {}

std::variant<double, int, std::string> LaplacianMatrixDeterminant::calculate(const Molecule& mol) const {
    if (!mol.isValid()) {
        return 0.0;
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (rdkMol->getNumAtoms() == 0) {
        return 0.0;
    }
    
    Eigen::MatrixXd adjacency = buildAdjacencyMatrix(rdkMol);
    Eigen::MatrixXd degree = buildDegreeMatrix(adjacency);
    Eigen::MatrixXd laplacian = buildLaplacianMatrix(degree, adjacency);
    
    // For connected graphs, Laplacian is singular, so use the product of non-zero eigenvalues
    Eigen::VectorXd eigenvalues = computeEigenvalues(laplacian);
    double product = 1.0;
    int zeros = 0;
    
    for (int i = 0; i < eigenvalues.size(); ++i) {
        if (std::abs(eigenvalues(i)) < 1e-10) {
            zeros++;
        } else {
            product *= eigenvalues(i);
        }
    }
    
    // If all eigenvalues are zero, return 0
    if (zeros == eigenvalues.size()) {
        return 0.0;
    }
    
    return product;
}

LaplacianTotalEffectiveResistance::LaplacianTotalEffectiveResistance()
    : EigenDescriptor("LapTotEffRes", "Total effective resistance (sum of inverse nonzero Laplacian eigenvalues)") {}

std::variant<double, int, std::string> LaplacianTotalEffectiveResistance::calculate(const Molecule& mol) const {
    if (!mol.isValid()) {
        return 0.0;
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (rdkMol->getNumAtoms() <= 1) {
        return 0.0;
    }
    
    Eigen::MatrixXd adjacency = buildAdjacencyMatrix(rdkMol);
    Eigen::MatrixXd degree = buildDegreeMatrix(adjacency);
    Eigen::MatrixXd laplacian = buildLaplacianMatrix(degree, adjacency);
    Eigen::VectorXd eigenvalues = computeEigenvalues(laplacian);
    
    double sum = 0.0;
    for (int i = 0; i < eigenvalues.size(); ++i) {
        if (eigenvalues(i) > 1e-10) {
            sum += 1.0 / eigenvalues(i);
        }
    }
    
    return sum;
}

LaplacianKirchhoffIndex::LaplacianKirchhoffIndex()
    : EigenDescriptor("LapKirchhoff", "Kirchhoff index (sum of reciprocal Laplacian eigenvalues)") {}

std::variant<double, int, std::string> LaplacianKirchhoffIndex::calculate(const Molecule& mol) const {
    if (!mol.isValid()) {
        return 0.0;
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (rdkMol->getNumAtoms() <= 1) {
        return 0.0;
    }
    
    Eigen::MatrixXd adjacency = buildAdjacencyMatrix(rdkMol);
    Eigen::MatrixXd degree = buildDegreeMatrix(adjacency);
    Eigen::MatrixXd laplacian = buildLaplacianMatrix(degree, adjacency);
    Eigen::VectorXd eigenvalues = computeEigenvalues(laplacian);
    
    double sum = 0.0;
    int n = eigenvalues.size();
    
    for (int i = 0; i < n; ++i) {
        if (eigenvalues(i) > 1e-10) {
            sum += 1.0 / eigenvalues(i);
        }
    }
    
    return n * sum;
}

LaplacianEigenvalueVariance::LaplacianEigenvalueVariance()
    : EigenDescriptor("LapEigVar", "Variance of Laplacian eigenvalues") {}

std::variant<double, int, std::string> LaplacianEigenvalueVariance::calculate(const Molecule& mol) const {
    if (!mol.isValid()) {
        return 0.0;
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (rdkMol->getNumAtoms() == 0) {
        return 0.0;
    }
    
    Eigen::MatrixXd adjacency = buildAdjacencyMatrix(rdkMol);
    Eigen::MatrixXd degree = buildDegreeMatrix(adjacency);
    Eigen::MatrixXd laplacian = buildLaplacianMatrix(degree, adjacency);
    Eigen::VectorXd eigenvalues = computeEigenvalues(laplacian);
    
    return computeVariance(eigenvalues);
}

LaplacianEigenvalueSkewness::LaplacianEigenvalueSkewness()
    : EigenDescriptor("LapEigSkew", "Skewness of Laplacian eigenvalues") {}

std::variant<double, int, std::string> LaplacianEigenvalueSkewness::calculate(const Molecule& mol) const {
    if (!mol.isValid()) {
        return 0.0;
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (rdkMol->getNumAtoms() < 3) {
        return 0.0;
    }
    
    Eigen::MatrixXd adjacency = buildAdjacencyMatrix(rdkMol);
    Eigen::MatrixXd degree = buildDegreeMatrix(adjacency);
    Eigen::MatrixXd laplacian = buildLaplacianMatrix(degree, adjacency);
    Eigen::VectorXd eigenvalues = computeEigenvalues(laplacian);
    
    return computeSkewness(eigenvalues);
}

// Normalized Laplacian descriptors
NormalizedLaplacianSpectralRadius::NormalizedLaplacianSpectralRadius()
    : EigenDescriptor("NormLapSpecRad", "Spectral radius of normalized Laplacian") {}

std::variant<double, int, std::string> NormalizedLaplacianSpectralRadius::calculate(const Molecule& mol) const {
    if (!mol.isValid()) {
        return 0.0;
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (rdkMol->getNumAtoms() == 0) {
        return 0.0;
    }
    
    Eigen::MatrixXd adjacency = buildAdjacencyMatrix(rdkMol);
    Eigen::MatrixXd degree = buildDegreeMatrix(adjacency);
    Eigen::MatrixXd normLaplacian = buildNormalizedLaplacian(degree, adjacency);
    Eigen::VectorXd eigenvalues = computeEigenvalues(normLaplacian);
    
    return eigenvalues.maxCoeff();
}

NormalizedLaplacianSmallestNonzero::NormalizedLaplacianSmallestNonzero()
    : EigenDescriptor("NormLapMinNonZero", "Smallest nonzero normalized Laplacian eigenvalue") {}

std::variant<double, int, std::string> NormalizedLaplacianSmallestNonzero::calculate(const Molecule& mol) const {
    if (!mol.isValid()) {
        return 0.0;
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (rdkMol->getNumAtoms() <= 1) {
        return 0.0;
    }
    
    Eigen::MatrixXd adjacency = buildAdjacencyMatrix(rdkMol);
    Eigen::MatrixXd degree = buildDegreeMatrix(adjacency);
    Eigen::MatrixXd normLaplacian = buildNormalizedLaplacian(degree, adjacency);
    Eigen::VectorXd eigenvalues = computeEigenvalues(normLaplacian);
    
    // Sort eigenvalues
    std::vector<double> sortedEigenvalues(eigenvalues.data(), eigenvalues.data() + eigenvalues.size());
    std::sort(sortedEigenvalues.begin(), sortedEigenvalues.end());
    
    // Find first non-zero eigenvalue
    double smallestNonZero = 0.0;
    for (double val : sortedEigenvalues) {
        if (std::abs(val) > 1e-10) {
            smallestNonZero = val;
            break;
        }
    }
    
    return smallestNonZero;
}

NormalizedLaplacianLargestEigenvalue::NormalizedLaplacianLargestEigenvalue()
    : EigenDescriptor("NormLapMaxEig", "Largest normalized Laplacian eigenvalue") {}

std::variant<double, int, std::string> NormalizedLaplacianLargestEigenvalue::calculate(const Molecule& mol) const {
    if (!mol.isValid()) {
        return 0.0;
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (rdkMol->getNumAtoms() == 0) {
        return 0.0;
    }
    
    Eigen::MatrixXd adjacency = buildAdjacencyMatrix(rdkMol);
    Eigen::MatrixXd degree = buildDegreeMatrix(adjacency);
    Eigen::MatrixXd normLaplacian = buildNormalizedLaplacian(degree, adjacency);
    Eigen::VectorXd eigenvalues = computeEigenvalues(normLaplacian);
    
    return eigenvalues.maxCoeff();
}

NormalizedLaplacianEnergy::NormalizedLaplacianEnergy()
    : EigenDescriptor("NormLapEnergy", "Normalized Laplacian energy") {}

std::variant<double, int, std::string> NormalizedLaplacianEnergy::calculate(const Molecule& mol) const {
    if (!mol.isValid()) {
        return 0.0;
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (rdkMol->getNumAtoms() == 0) {
        return 0.0;
    }
    
    Eigen::MatrixXd adjacency = buildAdjacencyMatrix(rdkMol);
    Eigen::MatrixXd degree = buildDegreeMatrix(adjacency);
    Eigen::MatrixXd normLaplacian = buildNormalizedLaplacian(degree, adjacency);
    Eigen::VectorXd eigenvalues = computeEigenvalues(normLaplacian);
    
    int n = eigenvalues.size();
    double sum = 0.0;
    double avg = n > 0 ? n / static_cast<double>(n) : 0;
    
    for (int i = 0; i < n; ++i) {
        sum += std::abs(eigenvalues(i) - avg);
    }
    
    return sum;
}

NormalizedLaplacianTrace::NormalizedLaplacianTrace()
    : EigenDescriptor("NormLapTrace", "Normalized Laplacian trace") {}

std::variant<double, int, std::string> NormalizedLaplacianTrace::calculate(const Molecule& mol) const {
    if (!mol.isValid()) {
        return 0.0;
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (rdkMol->getNumAtoms() == 0) {
        return 0.0;
    }
    
    Eigen::MatrixXd adjacency = buildAdjacencyMatrix(rdkMol);
    Eigen::MatrixXd degree = buildDegreeMatrix(adjacency);
    Eigen::MatrixXd normLaplacian = buildNormalizedLaplacian(degree, adjacency);
    
    return normLaplacian.trace();
}

// Degree matrix descriptors
DegreeMatrixMaxDegree::DegreeMatrixMaxDegree()
    : EigenDescriptor("DegMatMaxDeg", "Maximum degree from diagonal of degree matrix") {}

std::variant<double, int, std::string> DegreeMatrixMaxDegree::calculate(const Molecule& mol) const {
    if (!mol.isValid()) {
        return 0;
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (rdkMol->getNumAtoms() == 0) {
        return 0;
    }
    
    Eigen::MatrixXd adjacency = buildAdjacencyMatrix(rdkMol);
    Eigen::MatrixXd degree = buildDegreeMatrix(adjacency);
    
    double maxDegree = 0.0;
    for (int i = 0; i < degree.rows(); ++i) {
        maxDegree = std::max(maxDegree, degree(i, i));
    }
    
    return static_cast<int>(maxDegree);
}

DegreeMatrixMinDegree::DegreeMatrixMinDegree()
    : EigenDescriptor("DegMatMinDeg", "Minimum degree from diagonal of degree matrix") {}

std::variant<double, int, std::string> DegreeMatrixMinDegree::calculate(const Molecule& mol) const {
    if (!mol.isValid()) {
        return 0;
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (rdkMol->getNumAtoms() == 0) {
        return 0;
    }
    
    Eigen::MatrixXd adjacency = buildAdjacencyMatrix(rdkMol);
    Eigen::MatrixXd degree = buildDegreeMatrix(adjacency);
    
    double minDegree = std::numeric_limits<double>::max();
    for (int i = 0; i < degree.rows(); ++i) {
        if (degree(i, i) > 0.0) {
            minDegree = std::min(minDegree, degree(i, i));
        }
    }
    
    return static_cast<int>(minDegree == std::numeric_limits<double>::max() ? 0 : minDegree);
}

DegreeMatrixAvgDegree::DegreeMatrixAvgDegree()
    : EigenDescriptor("DegMatAvgDeg", "Average degree from diagonal of degree matrix") {}

std::variant<double, int, std::string> DegreeMatrixAvgDegree::calculate(const Molecule& mol) const {
    if (!mol.isValid()) {
        return 0.0;
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (rdkMol->getNumAtoms() == 0) {
        return 0.0;
    }
    
    Eigen::MatrixXd adjacency = buildAdjacencyMatrix(rdkMol);
    Eigen::MatrixXd degree = buildDegreeMatrix(adjacency);
    
    double sumDegrees = degree.trace();
    return sumDegrees / degree.rows();
}

DegreeMatrixVariance::DegreeMatrixVariance()
    : EigenDescriptor("DegMatVar", "Variance of degrees from diagonal of degree matrix") {}

std::variant<double, int, std::string> DegreeMatrixVariance::calculate(const Molecule& mol) const {
    if (!mol.isValid()) {
        return 0.0;
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (rdkMol->getNumAtoms() == 0) {
        return 0.0;
    }
    
    Eigen::MatrixXd adjacency = buildAdjacencyMatrix(rdkMol);
    Eigen::MatrixXd degree = buildDegreeMatrix(adjacency);
    
    Eigen::VectorXd degrees(degree.rows());
    for (int i = 0; i < degree.rows(); ++i) {
        degrees(i) = degree(i, i);
    }
    
    return computeVariance(degrees);
}

DegreeMatrixSkewness::DegreeMatrixSkewness()
    : EigenDescriptor("DegMatSkew", "Skewness of degrees from diagonal of degree matrix") {}

std::variant<double, int, std::string> DegreeMatrixSkewness::calculate(const Molecule& mol) const {
    if (!mol.isValid()) {
        return 0.0;
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (rdkMol->getNumAtoms() < 3) {
        return 0.0;
    }
    
    Eigen::MatrixXd adjacency = buildAdjacencyMatrix(rdkMol);
    Eigen::MatrixXd degree = buildDegreeMatrix(adjacency);
    
    Eigen::VectorXd degrees(degree.rows());
    for (int i = 0; i < degree.rows(); ++i) {
        degrees(i) = degree(i, i);
    }
    
    return computeSkewness(degrees);
}

DegreeMatrixEntropy::DegreeMatrixEntropy()
    : EigenDescriptor("DegMatEnt", "Entropy of degrees from diagonal of degree matrix") {}

std::variant<double, int, std::string> DegreeMatrixEntropy::calculate(const Molecule& mol) const {
    if (!mol.isValid()) {
        return 0.0;
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (rdkMol->getNumAtoms() == 0) {
        return 0.0;
    }
    
    Eigen::MatrixXd adjacency = buildAdjacencyMatrix(rdkMol);
    Eigen::MatrixXd degree = buildDegreeMatrix(adjacency);
    
    // Count occurrences of each degree value
    std::map<int, int> degreeCounts;
    for (int i = 0; i < degree.rows(); ++i) {
        int deg = static_cast<int>(degree(i, i));
        degreeCounts[deg]++;
    }
    
    // Calculate entropy
    double entropy = 0.0;
    int numAtoms = degree.rows();
    for (const auto& [deg, count] : degreeCounts) {
        double prob = static_cast<double>(count) / numAtoms;
        if (prob > 0.0) {
            entropy -= prob * std::log2(prob);
        }
    }
    
    return entropy;
}


// Higher-order adjacency powers descriptors
std::variant<double, int, std::string> NumberOf2Walks::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return 0;
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (rdkMol->getNumAtoms() == 0) return 0;
    
    Eigen::MatrixXd adjacency = buildAdjacencyMatrix(rdkMol);
    Eigen::MatrixXd squared = adjacency * adjacency;
    
    return static_cast<int>(squared.trace());
}

std::variant<double, int, std::string> NumberOf3Walks::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return 0;
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (rdkMol->getNumAtoms() == 0) return 0;
    
    Eigen::MatrixXd adjacency = buildAdjacencyMatrix(rdkMol);
    Eigen::MatrixXd squared = adjacency * adjacency;
    Eigen::MatrixXd cubed = squared * adjacency;
    
    return static_cast<int>(cubed.trace());
}

std::variant<double, int, std::string> NumberOf4Walks::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return 0;
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (rdkMol->getNumAtoms() == 0) return 0;
    
    Eigen::MatrixXd adjacency = buildAdjacencyMatrix(rdkMol);
    Eigen::MatrixXd squared = adjacency * adjacency;
    Eigen::MatrixXd fourth = squared * squared;
    
    return static_cast<int>(fourth.trace());
}

std::variant<double, int, std::string> MeanClosed3WalksPerNode::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return 0.0;
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    int numAtoms = rdkMol->getNumAtoms();
    if (numAtoms == 0) return 0.0;
    
    Eigen::MatrixXd adjacency = buildAdjacencyMatrix(rdkMol);
    Eigen::MatrixXd squared = adjacency * adjacency;
    Eigen::MatrixXd cubed = squared * adjacency;
    
    return cubed.trace() / numAtoms;
}

std::variant<double, int, std::string> MeanClosed4WalksPerNode::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return 0.0;
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    int numAtoms = rdkMol->getNumAtoms();
    if (numAtoms == 0) return 0.0;
    
    Eigen::MatrixXd adjacency = buildAdjacencyMatrix(rdkMol);
    Eigen::MatrixXd squared = adjacency * adjacency;
    Eigen::MatrixXd fourth = squared * squared;
    
    return fourth.trace() / numAtoms;
}

std::variant<double, int, std::string> Walk2Energy::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return 0.0;
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (rdkMol->getNumAtoms() == 0) return 0.0;
    
    Eigen::MatrixXd adjacency = buildAdjacencyMatrix(rdkMol);
    Eigen::MatrixXd squared = adjacency * adjacency;
    Eigen::VectorXd singularValues = computeSingularValues(squared);
    
    return singularValues.sum();
}

std::variant<double, int, std::string> Walk3Energy::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return 0.0;
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (rdkMol->getNumAtoms() == 0) return 0.0;
    
    Eigen::MatrixXd adjacency = buildAdjacencyMatrix(rdkMol);
    Eigen::MatrixXd squared = adjacency * adjacency;
    Eigen::MatrixXd cubed = squared * adjacency;
    Eigen::VectorXd singularValues = computeSingularValues(cubed);
    
    return singularValues.sum();
}

std::variant<double, int, std::string> Walk4Energy::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return 0.0;
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (rdkMol->getNumAtoms() == 0) return 0.0;
    
    Eigen::MatrixXd adjacency = buildAdjacencyMatrix(rdkMol);
    Eigen::MatrixXd squared = adjacency * adjacency;
    Eigen::MatrixXd fourth = squared * squared;
    Eigen::VectorXd singularValues = computeSingularValues(fourth);
    
    return singularValues.sum();
}

std::variant<double, int, std::string> GraphIrregularityWalkCount::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return 0.0;
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (rdkMol->getNumAtoms() <= 1) return 0.0;
    
    Eigen::MatrixXd adjacency = buildAdjacencyMatrix(rdkMol);
    Eigen::MatrixXd squared = adjacency * adjacency;
    
    std::vector<double> walkCounts;
    for (int i = 0; i < squared.rows(); ++i) {
        walkCounts.push_back(squared(i, i)); // Diagonal elements are closed walk counts
    }
    
    double irregularity = 0.0;
    for (size_t i = 0; i < walkCounts.size(); ++i) {
        for (size_t j = i + 1; j < walkCounts.size(); ++j) {
            irregularity += std::abs(walkCounts[i] - walkCounts[j]);
        }
    }
    
    return irregularity;
}

std::variant<double, int, std::string> TrianglesToPathsRatio::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return 0.0;
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (rdkMol->getNumAtoms() < 3) return 0.0;
    
    Eigen::MatrixXd adjacency = buildAdjacencyMatrix(rdkMol);
    Eigen::MatrixXd squared = adjacency * adjacency;
    Eigen::MatrixXd cubed = squared * adjacency;
    
    double paths = squared.trace();
    double triangles = cubed.trace() / 6.0;
    
    return paths > 0 ? triangles / paths : 0.0;
}

// SVD-based descriptors
std::variant<double, int, std::string> MaxSingularValue::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return 0.0;
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (rdkMol->getNumAtoms() == 0) return 0.0;
    
    Eigen::MatrixXd adjacency = buildAdjacencyMatrix(rdkMol);
    Eigen::VectorXd singularValues = computeSingularValues(adjacency);
    
    return singularValues.size() > 0 ? singularValues(0) : 0.0;
}

std::variant<double, int, std::string> MinNonZeroSingularValue::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return 0.0;
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (rdkMol->getNumAtoms() == 0) return 0.0;
    
    Eigen::MatrixXd adjacency = buildAdjacencyMatrix(rdkMol);
    Eigen::VectorXd singularValues = computeSingularValues(adjacency);
    
    double minNonZero = 0.0;
    for (int i = 0; i < singularValues.size(); ++i) {
        if (singularValues(i) > 1e-10 && (minNonZero == 0.0 || singularValues(i) < minNonZero)) {
            minNonZero = singularValues(i);
        }
    }
    
    return minNonZero;
}

std::variant<double, int, std::string> ConditionNumber::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return 0.0;
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (rdkMol->getNumAtoms() == 0) return 0.0;
    
    Eigen::MatrixXd adjacency = buildAdjacencyMatrix(rdkMol);
    Eigen::VectorXd singularValues = computeSingularValues(adjacency);
    
    if (singularValues.size() == 0) return 0.0;
    
    double maxSingular = singularValues(0);
    double minNonZero = 0.0;
    for (int i = 0; i < singularValues.size(); ++i) {
        if (singularValues(i) > 1e-10 && (minNonZero == 0.0 || singularValues(i) < minNonZero)) {
            minNonZero = singularValues(i);
        }
    }
    
    return minNonZero > 0 ? maxSingular / minNonZero : 0.0;
}

std::variant<double, int, std::string> SumSingularValues::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return 0.0;
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (rdkMol->getNumAtoms() == 0) return 0.0;
    
    Eigen::MatrixXd adjacency = buildAdjacencyMatrix(rdkMol);
    Eigen::VectorXd singularValues = computeSingularValues(adjacency);
    
    return singularValues.sum();
}

std::variant<double, int, std::string> FrobeniusNormSingularValues::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return 0.0;
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (rdkMol->getNumAtoms() == 0) return 0.0;
    
    Eigen::MatrixXd adjacency = buildAdjacencyMatrix(rdkMol);
    Eigen::VectorXd singularValues = computeSingularValues(adjacency);
    
    return std::sqrt(singularValues.squaredNorm());
}

std::variant<double, int, std::string> SingularValueEntropy::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return 0.0;
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (rdkMol->getNumAtoms() == 0) return 0.0;
    
    Eigen::MatrixXd adjacency = buildAdjacencyMatrix(rdkMol);
    Eigen::VectorXd singularValues = computeSingularValues(adjacency);
    
    double sum = singularValues.sum();
    if (sum < 1e-10) return 0.0;
    
    double entropy = 0.0;
    for (int i = 0; i < singularValues.size(); ++i) {
        if (singularValues(i) > 1e-10) {
            double p = singularValues(i) / sum;
            entropy -= p * std::log2(p);
        }
    }
    
    return entropy;
}

std::variant<double, int, std::string> SingularValueVariance::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return 0.0;
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (rdkMol->getNumAtoms() == 0) return 0.0;
    
    Eigen::MatrixXd adjacency = buildAdjacencyMatrix(rdkMol);
    Eigen::VectorXd singularValues = computeSingularValues(adjacency);
    
    return computeVariance(singularValues);
}

std::variant<double, int, std::string> SingularValueSkewness::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return 0.0;
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (rdkMol->getNumAtoms() < 3) return 0.0;
    
    Eigen::MatrixXd adjacency = buildAdjacencyMatrix(rdkMol);
    Eigen::VectorXd singularValues = computeSingularValues(adjacency);
    
    return computeSkewness(singularValues);
}

std::variant<double, int, std::string> SpectralEffectiveRank::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return 0.0;
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (rdkMol->getNumAtoms() == 0) return 0.0;
    
    Eigen::MatrixXd adjacency = buildAdjacencyMatrix(rdkMol);
    Eigen::VectorXd singularValues = computeSingularValues(adjacency);
    
    double sum = singularValues.sum();
    if (sum < 1e-10) return 0.0;
    
    double entropy = 0.0;
    for (int i = 0; i < singularValues.size(); ++i) {
        if (singularValues(i) > 1e-10) {
            double p = singularValues(i) / sum;
            entropy -= p * std::log2(p);
        }
    }
    
    return std::pow(2, entropy);
}

std::variant<double, int, std::string> NuclearNorm::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return 0.0;
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (rdkMol->getNumAtoms() == 0) return 0.0;
    
    Eigen::MatrixXd adjacency = buildAdjacencyMatrix(rdkMol);
    Eigen::VectorXd singularValues = computeSingularValues(adjacency);
    
    return singularValues.sum();
}

// Normalized adjacency matrix descriptors
std::variant<double, int, std::string> NormalizedAdjacencySpectralRadius::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return 0.0;
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (rdkMol->getNumAtoms() == 0) return 0.0;
    
    Eigen::MatrixXd adjacency = buildAdjacencyMatrix(rdkMol);
    Eigen::MatrixXd degree = buildDegreeMatrix(adjacency);
    Eigen::MatrixXd normAdjacency = buildNormalizedAdjacency(degree, adjacency);
    Eigen::VectorXd eigenvalues = computeEigenvalues(normAdjacency);
    
    return eigenvalues.cwiseAbs().maxCoeff();
}

std::variant<double, int, std::string> NormalizedEigenvalueGap::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return 0.0;
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (rdkMol->getNumAtoms() <= 1) return 0.0;
    
    Eigen::MatrixXd adjacency = buildAdjacencyMatrix(rdkMol);
    Eigen::MatrixXd degree = buildDegreeMatrix(adjacency);
    Eigen::MatrixXd normAdjacency = buildNormalizedAdjacency(degree, adjacency);
    Eigen::VectorXd eigenvalues = computeEigenvalues(normAdjacency);
    
    std::vector<double> absEigenvalues;
    for (int i = 0; i < eigenvalues.size(); ++i) {
        absEigenvalues.push_back(std::abs(eigenvalues(i)));
    }
    
    std::sort(absEigenvalues.begin(), absEigenvalues.end(), std::greater<double>());
    
    if (absEigenvalues.size() >= 2) {
        return absEigenvalues[0] - absEigenvalues[1];
    }
    
    return 0.0;
}

std::variant<double, int, std::string> SumNormalizedEigenvalues::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return 0.0;
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (rdkMol->getNumAtoms() == 0) return 0.0;
    
    Eigen::MatrixXd adjacency = buildAdjacencyMatrix(rdkMol);
    Eigen::MatrixXd degree = buildDegreeMatrix(adjacency);
    Eigen::MatrixXd normAdjacency = buildNormalizedAdjacency(degree, adjacency);
    Eigen::VectorXd eigenvalues = computeEigenvalues(normAdjacency);
    
    return eigenvalues.sum();
}

std::variant<double, int, std::string> VarianceNormalizedEigenvalues::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return 0.0;
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (rdkMol->getNumAtoms() == 0) return 0.0;
    
    Eigen::MatrixXd adjacency = buildAdjacencyMatrix(rdkMol);
    Eigen::MatrixXd degree = buildDegreeMatrix(adjacency);
    Eigen::MatrixXd normAdjacency = buildNormalizedAdjacency(degree, adjacency);
    Eigen::VectorXd eigenvalues = computeEigenvalues(normAdjacency);
    
    return computeVariance(eigenvalues);
}

std::variant<double, int, std::string> CountNormalizedEigenvaluesAboveHalf::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return 0;
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (rdkMol->getNumAtoms() == 0) return 0;
    
    Eigen::MatrixXd adjacency = buildAdjacencyMatrix(rdkMol);
    Eigen::MatrixXd degree = buildDegreeMatrix(adjacency);
    Eigen::MatrixXd normAdjacency = buildNormalizedAdjacency(degree, adjacency);
    Eigen::VectorXd eigenvalues = computeEigenvalues(normAdjacency);
    
    int count = 0;
    for (int i = 0; i < eigenvalues.size(); ++i) {
        if (eigenvalues(i) > 0.5) {
            count++;
        }
    }
    
    return count;
}

std::variant<double, int, std::string> NormalizedEnergy::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return 0.0;
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (rdkMol->getNumAtoms() == 0) return 0.0;
    
    Eigen::MatrixXd adjacency = buildAdjacencyMatrix(rdkMol);
    Eigen::MatrixXd degree = buildDegreeMatrix(adjacency);
    Eigen::MatrixXd normAdjacency = buildNormalizedAdjacency(degree, adjacency);
    Eigen::VectorXd eigenvalues = computeEigenvalues(normAdjacency);
    
    return eigenvalues.cwiseAbs().sum();
}

std::variant<double, int, std::string> LargestNormalizedEigenvectorCentrality::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return 0.0;
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (rdkMol->getNumAtoms() == 0) return 0.0;
    
    Eigen::MatrixXd adjacency = buildAdjacencyMatrix(rdkMol);
    Eigen::MatrixXd degree = buildDegreeMatrix(adjacency);
    Eigen::MatrixXd normAdjacency = buildNormalizedAdjacency(degree, adjacency);
    
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(normAdjacency);
    Eigen::VectorXd eigenvalues = solver.eigenvalues();
    Eigen::MatrixXd eigenvectors = solver.eigenvectors();
    
    int maxIndex = 0;
    for (int i = 1; i < eigenvalues.size(); ++i) {
        if (std::abs(eigenvalues(i)) > std::abs(eigenvalues(maxIndex))) {
            maxIndex = i;
        }
    }
    
    Eigen::VectorXd dominantEigenvector = eigenvectors.col(maxIndex);
    return dominantEigenvector.cwiseAbs().maxCoeff();
}

std::variant<double, int, std::string> AverageNormalizedEigenvectorCentrality::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return 0.0;
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (rdkMol->getNumAtoms() == 0) return 0.0;
    
    Eigen::MatrixXd adjacency = buildAdjacencyMatrix(rdkMol);
    Eigen::MatrixXd degree = buildDegreeMatrix(adjacency);
    Eigen::MatrixXd normAdjacency = buildNormalizedAdjacency(degree, adjacency);
    
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(normAdjacency);
    Eigen::VectorXd eigenvalues = solver.eigenvalues();
    Eigen::MatrixXd eigenvectors = solver.eigenvectors();
    
    int maxIndex = 0;
    for (int i = 1; i < eigenvalues.size(); ++i) {
        if (std::abs(eigenvalues(i)) > std::abs(eigenvalues(maxIndex))) {
            maxIndex = i;
        }
    }
    
    Eigen::VectorXd dominantEigenvector = eigenvectors.col(maxIndex);
    return dominantEigenvector.cwiseAbs().mean();
}

std::variant<double, int, std::string> NormalizedAdjacencyMatrixRank::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return 0;
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (rdkMol->getNumAtoms() == 0) return 0;
    
    Eigen::MatrixXd adjacency = buildAdjacencyMatrix(rdkMol);
    Eigen::MatrixXd degree = buildDegreeMatrix(adjacency);
    Eigen::MatrixXd normAdjacency = buildNormalizedAdjacency(degree, adjacency);
    
    Eigen::FullPivLU<Eigen::MatrixXd> lu(normAdjacency);
    return static_cast<int>(lu.rank());
}

std::variant<double, int, std::string> NormalizedAdjacencyEntropy::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return 0.0;
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (rdkMol->getNumAtoms() == 0) return 0.0;
    
    Eigen::MatrixXd adjacency = buildAdjacencyMatrix(rdkMol);
    Eigen::MatrixXd degree = buildDegreeMatrix(adjacency);
    Eigen::MatrixXd normAdjacency = buildNormalizedAdjacency(degree, adjacency);
    Eigen::VectorXd eigenvalues = computeEigenvalues(normAdjacency);
    
    double entropy = 0.0;
    for (int i = 0; i < eigenvalues.size(); ++i) {
        double absVal = std::abs(eigenvalues(i));
        if (absVal > 1e-10) {
            entropy -= absVal * std::log2(absVal);
        }
    }
    
    return entropy;
}

// Signless Laplacian matrix descriptors
std::variant<double, int, std::string> SignlessLaplacianSpectralRadius::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return 0.0;
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (rdkMol->getNumAtoms() == 0) return 0.0;
    
    Eigen::MatrixXd adjacency = buildAdjacencyMatrix(rdkMol);
    Eigen::MatrixXd degree = buildDegreeMatrix(adjacency);
    Eigen::MatrixXd signlessLaplacian = buildSignlessLaplacian(degree, adjacency);
    Eigen::VectorXd eigenvalues = computeEigenvalues(signlessLaplacian);
    
    return eigenvalues.maxCoeff();
}

std::variant<double, int, std::string> SignlessLaplacianSmallestEigenvalue::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return 0.0;
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (rdkMol->getNumAtoms() == 0) return 0.0;
    
    Eigen::MatrixXd adjacency = buildAdjacencyMatrix(rdkMol);
    Eigen::MatrixXd degree = buildDegreeMatrix(adjacency);
    Eigen::MatrixXd signlessLaplacian = buildSignlessLaplacian(degree, adjacency);
    Eigen::VectorXd eigenvalues = computeEigenvalues(signlessLaplacian);
    
    return eigenvalues.minCoeff();
}

std::variant<double, int, std::string> SignlessLaplacianEnergy::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return 0.0;
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (rdkMol->getNumAtoms() == 0) return 0.0;
    
    Eigen::MatrixXd adjacency = buildAdjacencyMatrix(rdkMol);
    Eigen::MatrixXd degree = buildDegreeMatrix(adjacency);
    Eigen::MatrixXd signlessLaplacian = buildSignlessLaplacian(degree, adjacency);
    Eigen::VectorXd eigenvalues = computeEigenvalues(signlessLaplacian);
    
    return eigenvalues.cwiseAbs().sum();
}

std::variant<double, int, std::string> SignlessLaplacianTrace::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return 0.0;
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (rdkMol->getNumAtoms() == 0) return 0.0;
    
    Eigen::MatrixXd adjacency = buildAdjacencyMatrix(rdkMol);
    Eigen::MatrixXd degree = buildDegreeMatrix(adjacency);
    Eigen::MatrixXd signlessLaplacian = buildSignlessLaplacian(degree, adjacency);
    
    return signlessLaplacian.trace();
}

std::variant<double, int, std::string> SignlessLaplacianDeterminant::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return 0.0;
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (rdkMol->getNumAtoms() == 0) return 0.0;
    
    Eigen::MatrixXd adjacency = buildAdjacencyMatrix(rdkMol);
    Eigen::MatrixXd degree = buildDegreeMatrix(adjacency);
    Eigen::MatrixXd signlessLaplacian = buildSignlessLaplacian(degree, adjacency);
    
    return signlessLaplacian.determinant();
}


// Other linear algebra properties
GraphIrregularity::GraphIrregularity()
    : EigenDescriptor("GraphIrreg", "Graph irregularity (sum of absolute differences of degrees)") {}

std::variant<double, int, std::string> GraphIrregularity::calculate(const Molecule& mol) const {
    if (!mol.isValid()) {
        return 0.0;
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (rdkMol->getNumAtoms() <= 1) {
        return 0.0;
    }
    
    Eigen::MatrixXd adjacency = buildAdjacencyMatrix(rdkMol);
    
    std::vector<double> degrees;
    for (int i = 0; i < adjacency.rows(); ++i) {
        degrees.push_back(adjacency.row(i).sum());
    }
    
    double irregularity = 0.0;
    for (size_t i = 0; i < degrees.size(); ++i) {
        for (size_t j = i + 1; j < degrees.size(); ++j) {
            irregularity += std::abs(degrees[i] - degrees[j]);
        }
    }
    
    return irregularity;
}

WienerIndex::WienerIndex()
    : EigenDescriptor("WienerIdx", "Wiener index (sum of all shortest-path distances)") {}

std::variant<double, int, std::string> WienerIndex::calculate(const Molecule& mol) const {
    if (!mol.isValid()) {
        return 0.0;
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (rdkMol->getNumAtoms() <= 1) {
        return 0.0;
    }
    
    int n = rdkMol->getNumAtoms();
    Eigen::MatrixXd adjacency = buildAdjacencyMatrix(rdkMol);
    
    // Initialize distance matrix with adjacency
    Eigen::MatrixXd dist = Eigen::MatrixXd::Constant(n, n, std::numeric_limits<double>::infinity());
    for (int i = 0; i < n; ++i) {
        dist(i, i) = 0.0;  // Distance to self is 0
        for (int j = 0; j < n; ++j) {
            if (adjacency(i, j) > 0.5) {  // If there's an edge
                dist(i, j) = 1.0;
            }
        }
    }
    
    // Floyd-Warshall algorithm for all-pairs shortest paths
    for (int k = 0; k < n; ++k) {
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                if (dist(i, k) + dist(k, j) < dist(i, j)) {
                    dist(i, j) = dist(i, k) + dist(k, j);
                }
            }
        }
    }
    
    // Sum up all distances
    double wienerIndex = 0.0;
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            if (dist(i, j) != std::numeric_limits<double>::infinity()) {
                wienerIndex += dist(i, j);
            }
        }
    }
    
    return wienerIndex;
}

EstradaIndex::EstradaIndex()
    : EigenDescriptor("EstradaIdx", "Estrada index (sum of exp(eigenvalues of adjacency matrix))") {}

std::variant<double, int, std::string> EstradaIndex::calculate(const Molecule& mol) const {
    if (!mol.isValid()) {
        return 0.0;
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (rdkMol->getNumAtoms() == 0) {
        return 0.0;
    }
    
    Eigen::MatrixXd adjacency = buildAdjacencyMatrix(rdkMol);
    Eigen::VectorXd eigenvalues = computeEigenvalues(adjacency);
    
    double estradaIndex = 0.0;
    for (int i = 0; i < eigenvalues.size(); ++i) {
        estradaIndex += std::exp(eigenvalues(i));
    }
    
    return estradaIndex;
}

NumberSpanningTrees::NumberSpanningTrees()
    : EigenDescriptor("NumSpanTrees", "Number of spanning trees (via Laplacian minor)") {}

std::variant<double, int, std::string> NumberSpanningTrees::calculate(const Molecule& mol) const {
    if (!mol.isValid()) {
        return 0;
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (rdkMol->getNumAtoms() <= 1) {
        return 0;
    }
    
    Eigen::MatrixXd adjacency = buildAdjacencyMatrix(rdkMol);
    Eigen::MatrixXd degree = buildDegreeMatrix(adjacency);
    Eigen::MatrixXd laplacian = buildLaplacianMatrix(degree, adjacency);
    
    // Get the first minor (remove last row and column)
    int n = laplacian.rows();
    if (n <= 1) return 0;
    
    Eigen::MatrixXd minor = laplacian.topLeftCorner(n-1, n-1);
    double det = minor.determinant();
    
    // Return rounded value as integer
    return static_cast<int>(std::round(std::abs(det)));
}

GraphEccentricity::GraphEccentricity()
    : EigenDescriptor("GraphEcc", "Graph eccentricity (maximum eccentricity of any vertex)") {}

std::variant<double, int, std::string> GraphEccentricity::calculate(const Molecule& mol) const {
    if (!mol.isValid()) {
        return 0;
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (rdkMol->getNumAtoms() <= 1) {
        return 0;
    }
    
    int n = rdkMol->getNumAtoms();
    Eigen::MatrixXd adjacency = buildAdjacencyMatrix(rdkMol);
    
    // Initialize distance matrix with adjacency
    Eigen::MatrixXd dist = Eigen::MatrixXd::Constant(n, n, std::numeric_limits<double>::infinity());
    for (int i = 0; i < n; ++i) {
        dist(i, i) = 0.0;  // Distance to self is 0
        for (int j = 0; j < n; ++j) {
            if (adjacency(i, j) > 0.5) {  // If there's an edge
                dist(i, j) = 1.0;
            }
        }
    }
    
    // Floyd-Warshall algorithm for all-pairs shortest paths
    for (int k = 0; k < n; ++k) {
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                if (dist(i, k) + dist(k, j) < dist(i, j)) {
                    dist(i, j) = dist(i, k) + dist(k, j);
                }
            }
        }
    }
    
    // Find eccentricity of each vertex and return maximum
    double maxEccentricity = 0.0;
    for (int i = 0; i < n; ++i) {
        double eccentricity = 0.0;
        for (int j = 0; j < n; ++j) {
            if (dist(i, j) != std::numeric_limits<double>::infinity()) {
                eccentricity = std::max(eccentricity, dist(i, j));
            }
        }
        maxEccentricity = std::max(maxEccentricity, eccentricity);
    }
    
    return static_cast<int>(maxEccentricity);
}

SpectralGap::SpectralGap()
    : EigenDescriptor("SpecGap", "Spectral gap (difference between largest and second-largest eigenvalues)") {}

std::variant<double, int, std::string> SpectralGap::calculate(const Molecule& mol) const {
    if (!mol.isValid()) {
        return 0.0;
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (rdkMol->getNumAtoms() <= 1) {
        return 0.0;
    }
    
    Eigen::MatrixXd adjacency = buildAdjacencyMatrix(rdkMol);
    Eigen::VectorXd eigenvalues = computeEigenvalues(adjacency);
    
    // Sort eigenvalues in descending order
    std::vector<double> sortedEigenvalues(eigenvalues.data(), eigenvalues.data() + eigenvalues.size());
    std::sort(sortedEigenvalues.begin(), sortedEigenvalues.end(), std::greater<double>());
    
    // Gap between the two largest eigenvalues
    double gap = 0.0;
    if (sortedEigenvalues.size() >= 2) {
        gap = sortedEigenvalues[0] - sortedEigenvalues[1];
    }
    
    return gap;
}

TraceMatrixPower2::TraceMatrixPower2()
    : EigenDescriptor("TraceMatPow2", "Trace of adjacency matrix squared (related to number of 2-walks)") {}

std::variant<double, int, std::string> TraceMatrixPower2::calculate(const Molecule& mol) const {
    if (!mol.isValid()) {
        return 0;
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (rdkMol->getNumAtoms() == 0) {
        return 0;
    }
    
    Eigen::MatrixXd adjacency = buildAdjacencyMatrix(rdkMol);
    Eigen::MatrixXd squared = adjacency * adjacency;
    
    return static_cast<int>(squared.trace());
}

NumberTriangles::NumberTriangles()
    : EigenDescriptor("NumTriangles", "Number of triangles in the graph (from trace of A³/6)") {}

std::variant<double, int, std::string> NumberTriangles::calculate(const Molecule& mol) const {
    if (!mol.isValid()) {
        return 0;
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (rdkMol->getNumAtoms() < 3) {
        return 0;
    }
    
    Eigen::MatrixXd adjacency = buildAdjacencyMatrix(rdkMol);
    Eigen::MatrixXd squared = adjacency * adjacency;
    Eigen::MatrixXd cubed = squared * adjacency;
    
    // Number of triangles is trace(A³)/6
    double triangles = cubed.trace() / 6.0;
    
    return static_cast<int>(std::round(triangles));
}

GraphDiameter::GraphDiameter()
    : EigenDescriptor("GraphDiam", "Graph diameter (maximum distance between any two vertices)") {}

std::variant<double, int, std::string> GraphDiameter::calculate(const Molecule& mol) const {
    if (!mol.isValid()) {
        return 0;
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (rdkMol->getNumAtoms() <= 1) {
        return 0;
    }
    
    int n = rdkMol->getNumAtoms();
    Eigen::MatrixXd adjacency = buildAdjacencyMatrix(rdkMol);
    
    // Initialize distance matrix with adjacency
    Eigen::MatrixXd dist = Eigen::MatrixXd::Constant(n, n, std::numeric_limits<double>::infinity());
    for (int i = 0; i < n; ++i) {
        dist(i, i) = 0.0;  // Distance to self is 0
        for (int j = 0; j < n; ++j) {
            if (adjacency(i, j) > 0.5) {  // If there's an edge
                dist(i, j) = 1.0;
            }
        }
    }
    
    // Floyd-Warshall algorithm for all-pairs shortest paths
    for (int k = 0; k < n; ++k) {
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                if (dist(i, k) + dist(k, j) < dist(i, j)) {
                    dist(i, j) = dist(i, k) + dist(k, j);
                }
            }
        }
    }
    
    // Find the maximum finite distance (diameter)
    double diameter = 0.0;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (dist(i, j) != std::numeric_limits<double>::infinity()) {
                diameter = std::max(diameter, dist(i, j));
            }
        }
    }
    
    return static_cast<int>(diameter);
}

// Higher-order adjacency powers descriptors
NumberOf2Walks::NumberOf2Walks()
    : EigenDescriptor("Num2Walks", "Number of 2-walks (Tr(A²))") {}

NumberOf3Walks::NumberOf3Walks()
    : EigenDescriptor("Num3Walks", "Number of 3-walks (Tr(A³))") {}

NumberOf4Walks::NumberOf4Walks()
    : EigenDescriptor("Num4Walks", "Number of 4-walks (Tr(A⁴))") {}

MeanClosed3WalksPerNode::MeanClosed3WalksPerNode()
    : EigenDescriptor("MeanClosed3Walks", "Mean number of closed 3-walks per node") {}

MeanClosed4WalksPerNode::MeanClosed4WalksPerNode()
    : EigenDescriptor("MeanClosed4Walks", "Mean number of closed 4-walks per node") {}

Walk2Energy::Walk2Energy()
    : EigenDescriptor("Walk2Energy", "2-walk energy (sum of singular values of A²)") {}

Walk3Energy::Walk3Energy()
    : EigenDescriptor("Walk3Energy", "3-walk energy (sum of singular values of A³)") {}

Walk4Energy::Walk4Energy()
    : EigenDescriptor("Walk4Energy", "4-walk energy (sum of singular values of A⁴)") {}

GraphIrregularityWalkCount::GraphIrregularityWalkCount()
    : EigenDescriptor("GraphIrregWalk", "Graph irregularity based on walk counts") {}

TrianglesToPathsRatio::TrianglesToPathsRatio()
    : EigenDescriptor("TriPathRatio", "Ratio Tr(A³)/Tr(A²) (triangles to paths)") {}

// SVD-based descriptors
MaxSingularValue::MaxSingularValue()
    : EigenDescriptor("MaxSingVal", "Maximum singular value of adjacency matrix") {}

MinNonZeroSingularValue::MinNonZeroSingularValue()
    : EigenDescriptor("MinNonZeroSingVal", "Minimum non-zero singular value of adjacency matrix") {}

ConditionNumber::ConditionNumber()
    : EigenDescriptor("CondNumber", "Condition number (max/min singular value)") {}

SumSingularValues::SumSingularValues()
    : EigenDescriptor("SumSingVal", "Sum of singular values of adjacency matrix") {}

FrobeniusNormSingularValues::FrobeniusNormSingularValues()
    : EigenDescriptor("FrobNormSingVal", "Frobenius norm from singular values") {}

SingularValueEntropy::SingularValueEntropy()
    : EigenDescriptor("SingValEntropy", "Entropy of singular values after normalization") {}

SingularValueVariance::SingularValueVariance()
    : EigenDescriptor("SingValVar", "Variance of singular values") {}

SingularValueSkewness::SingularValueSkewness()
    : EigenDescriptor("SingValSkew", "Skewness of singular values") {}

SpectralEffectiveRank::SpectralEffectiveRank()
    : EigenDescriptor("SpecEffRank", "Spectral effective rank (Shannon entropy of normalized singular values)") {}

NuclearNorm::NuclearNorm()
    : EigenDescriptor("NuclearNorm", "Nuclear norm (sum of singular values)") {}

NormalizedEigenvalueGap::NormalizedEigenvalueGap()
    : EigenDescriptor("NormEigGap", "Gap between first two normalized eigenvalues") {}

SumNormalizedEigenvalues::SumNormalizedEigenvalues()
    : EigenDescriptor("SumNormEig", "Sum of normalized eigenvalues") {}

VarianceNormalizedEigenvalues::VarianceNormalizedEigenvalues()
    : EigenDescriptor("VarNormEig", "Variance of normalized eigenvalues") {}

CountNormalizedEigenvaluesAboveHalf::CountNormalizedEigenvaluesAboveHalf()
    : EigenDescriptor("CountNormEigAboveHalf", "Number of normalized eigenvalues > 0.5") {}

NormalizedEnergy::NormalizedEnergy()
    : EigenDescriptor("NormEnergy", "Normalized energy (sum of absolute normalized eigenvalues)") {}

LargestNormalizedEigenvectorCentrality::LargestNormalizedEigenvectorCentrality()
    : EigenDescriptor("MaxNormEigCent", "Largest normalized eigenvector centrality") {}

AverageNormalizedEigenvectorCentrality::AverageNormalizedEigenvectorCentrality()
    : EigenDescriptor("AvgNormEigCent", "Average normalized eigenvector centrality") {}

NormalizedAdjacencyMatrixRank::NormalizedAdjacencyMatrixRank()
    : EigenDescriptor("NormAdjRank", "Normalized adjacency matrix rank") {}

NormalizedAdjacencyEntropy::NormalizedAdjacencyEntropy()
    : EigenDescriptor("NormAdjEntropy", "Normalized adjacency entropy") {}

// Signless Laplacian matrix descriptors
SignlessLaplacianSpectralRadius::SignlessLaplacianSpectralRadius()
    : EigenDescriptor("SignLapSpecRad", "Signless Laplacian spectral radius") {}

SignlessLaplacianSmallestEigenvalue::SignlessLaplacianSmallestEigenvalue()
    : EigenDescriptor("SignLapMinEig", "Smallest eigenvalue of signless Laplacian") {}

SignlessLaplacianEnergy::SignlessLaplacianEnergy()
    : EigenDescriptor("SignLapEnergy", "Energy of signless Laplacian (sum of absolute eigenvalues)") {}

SignlessLaplacianTrace::SignlessLaplacianTrace()
    : EigenDescriptor("SignLapTrace", "Trace of signless Laplacian") {}

SignlessLaplacianDeterminant::SignlessLaplacianDeterminant()
    : EigenDescriptor("SignLapDet", "Determinant of signless Laplacian") {}

SignlessLaplacianZeroEigenvalues::SignlessLaplacianZeroEigenvalues()
    : EigenDescriptor("SignLapZeroEig", "Number of zero eigenvalues of signless Laplacian") {}

SignlessLaplacianPositiveEigenvalues::SignlessLaplacianPositiveEigenvalues()
    : EigenDescriptor("SignLapPosEig", "Number of positive eigenvalues of signless Laplacian") {}

SignlessLaplacianNegativeEigenvalues::SignlessLaplacianNegativeEigenvalues()
    : EigenDescriptor("SignLapNegEig", "Number of negative eigenvalues of signless Laplacian") {}

SignlessLaplacianEigenvalueVariance::SignlessLaplacianEigenvalueVariance()
    : EigenDescriptor("SignLapEigVar", "Signless Laplacian eigenvalue variance") {}

SignlessLaplacianEigenvalueSkewness::SignlessLaplacianEigenvalueSkewness()
    : EigenDescriptor("SignLapEigSkew", "Signless Laplacian eigenvalue skewness") {}

// Miscellaneous graph matrices descriptors
WeightedAdjacencySpectralRadius::WeightedAdjacencySpectralRadius()
    : EigenDescriptor("WtAdjSpecRad", "Weighted adjacency spectral radius") {}

MeanFirstPassageTime::MeanFirstPassageTime()
    : EigenDescriptor("MeanFirstPassage", "Mean first passage time approximated by pseudo-inverse of Laplacian") {}

CommuteTimeDistance::CommuteTimeDistance()
    : EigenDescriptor("CommuteTime", "Commute time distance (using pseudo-inverse of Laplacian)") {}

KirchhoffIndexVariance::KirchhoffIndexVariance()
    : EigenDescriptor("KirchhoffVar", "Kirchhoff index variance among subgraphs") {}

EffectiveGraphResistanceDistribution::EffectiveGraphResistanceDistribution()
    : EigenDescriptor("EffGraphResDist", "Effective graph resistance distribution (subgraphs)") {}

LocalClusteringCoefficientDistribution::LocalClusteringCoefficientDistribution()
    : EigenDescriptor("LocalClustCoefDist", "Local clustering coefficient distribution (approximated with A³)") {}

GraphRobustnessIndex::GraphRobustnessIndex()
    : EigenDescriptor("GraphRobustIdx", "Graph robustness index (based on spectral gap)") {}

NormalizedEstradaIndex::NormalizedEstradaIndex()
    : EigenDescriptor("NormEstradaIdx", "Estrada index from normalized adjacency matrix") {}

GraphBipartivityIndex::GraphBipartivityIndex()
    : EigenDescriptor("GraphBipart", "Graph bipartivity index (based on eigenvalues)") {}

SpanningTreeEntropy::SpanningTreeEntropy()
    : EigenDescriptor("SpanTreeEntropy", "Spanning tree entropy (related to Laplacian eigenvalues)") {}

} // namespace descriptors
} // namespace desfact