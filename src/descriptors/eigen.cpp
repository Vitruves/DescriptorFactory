#include "descriptors/eigen.hpp"
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <GraphMol/GraphMol.h>
#include <GraphMol/BondIterators.h>
#include <GraphMol/AtomIterators.h>
#include <cmath>
#include <limits>
#include <vector>
#include <algorithm>
#include <map>
#include "utils.hpp"
#include <GraphMol/MolOps.h>

namespace desfact {
namespace descriptors {

// Initialize thread_local cache variables
thread_local std::unordered_map<const RDKit::ROMol*, EigenDescriptor::MatrixCacheEntry> EigenDescriptor::adjacencyCache;
thread_local std::unordered_map<const RDKit::ROMol*, EigenDescriptor::MatrixCacheEntry> EigenDescriptor::degreeCache;
thread_local std::unordered_map<const RDKit::ROMol*, EigenDescriptor::MatrixCacheEntry> EigenDescriptor::laplacianCache;
thread_local std::unordered_map<const RDKit::ROMol*, EigenDescriptor::MatrixCacheEntry> EigenDescriptor::normalizedLaplacianCache;
thread_local std::unordered_map<const RDKit::ROMol*, EigenDescriptor::MatrixCacheEntry> EigenDescriptor::signlessLaplacianCache;
thread_local std::unordered_map<const RDKit::ROMol*, EigenDescriptor::MatrixCacheEntry> EigenDescriptor::adjPow2Cache;
thread_local std::unordered_map<const RDKit::ROMol*, EigenDescriptor::MatrixCacheEntry> EigenDescriptor::adjPow3Cache;
thread_local std::unordered_map<const RDKit::ROMol*, EigenDescriptor::MatrixCacheEntry> EigenDescriptor::adjPow4Cache;
thread_local std::unordered_map<const RDKit::ROMol*, EigenDescriptor::MatrixCacheEntry> EigenDescriptor::weightedAdjacencyCache;
thread_local std::unordered_map<const RDKit::ROMol*, EigenDescriptor::MatrixCacheEntry> EigenDescriptor::distanceMatrixCache;
thread_local std::unordered_map<const RDKit::ROMol*, EigenDescriptor::MatrixCacheEntry> EigenDescriptor::atomicNumWeightedAdjCache;
thread_local std::unordered_map<const RDKit::ROMol*, EigenDescriptor::MatrixCacheEntry> EigenDescriptor::atomicNumWeightedLapCache;
thread_local std::unordered_map<const RDKit::ROMol*, EigenDescriptor::MatrixCacheEntry> EigenDescriptor::normalizedAdjCache;
thread_local std::unordered_map<const RDKit::ROMol*, EigenDescriptor::MatrixCacheEntry> EigenDescriptor::lapPow2Cache;

thread_local std::unordered_map<Eigen::MatrixXd, EigenDescriptor::EigenvalueCacheEntry, EigenDescriptor::MatrixHash, EigenDescriptor::MatrixEqual> EigenDescriptor::eigenvalueCache;
thread_local std::unordered_map<Eigen::MatrixXd, EigenDescriptor::EigenvalueCacheEntry, EigenDescriptor::MatrixHash, EigenDescriptor::MatrixEqual> EigenDescriptor::singularValueCache;
thread_local std::unordered_map<Eigen::MatrixXd, EigenDescriptor::EigenvalueCacheEntry, EigenDescriptor::MatrixHash, EigenDescriptor::MatrixEqual> EigenDescriptor::dominantEigenvectorCache;

// EigenDescriptor base class implementation
EigenDescriptor::EigenDescriptor(const std::string& name, const std::string& description)
    : Descriptor(name, description) {}

// Static method to clear all caches - Implementation belongs here
void EigenDescriptor::clearThreadLocalCaches() {
    adjacencyCache.clear();
    degreeCache.clear();
    laplacianCache.clear();
    normalizedLaplacianCache.clear();
    signlessLaplacianCache.clear();
    weightedAdjacencyCache.clear();
    distanceMatrixCache.clear();
    atomicNumWeightedAdjCache.clear();
    atomicNumWeightedLapCache.clear();
    normalizedAdjCache.clear();
    adjPow2Cache.clear();
    adjPow3Cache.clear();
    adjPow4Cache.clear();
    lapPow2Cache.clear();
    eigenvalueCache.clear();
    singularValueCache.clear();
    dominantEigenvectorCache.clear();
}

// --- Matrix Building Methods - Implementations belong here ---
Eigen::MatrixXd EigenDescriptor::buildAdjacencyMatrix(const RDKit::ROMol* mol) const {
    int numAtoms = mol->getNumAtoms();
    Eigen::MatrixXd adjacency = Eigen::MatrixXd::Zero(numAtoms, numAtoms);
    for (const RDKit::Bond* bond : mol->bonds()) {
        int startIdx = bond->getBeginAtomIdx();
        int endIdx = bond->getEndAtomIdx();
        adjacency(startIdx, endIdx) = 1.0;
        adjacency(endIdx, startIdx) = 1.0; // Undirected graph
    }
    return adjacency;
}

Eigen::MatrixXd EigenDescriptor::buildDegreeMatrix(const Eigen::MatrixXd& adjacency) const {
    int size = adjacency.rows();
    Eigen::MatrixXd degree = Eigen::MatrixXd::Zero(size, size);
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
    Eigen::MatrixXd dInvSqrt = Eigen::MatrixXd::Zero(size, size);
    for (int i = 0; i < size; ++i) {
        if (degree(i, i) > 1e-10) { // Avoid division by zero for isolated atoms
            dInvSqrt(i, i) = 1.0 / std::sqrt(degree(i, i));
        }
    }
    normalizedLaplacian -= dInvSqrt * adjacency * dInvSqrt;
    return normalizedLaplacian;
}

Eigen::MatrixXd EigenDescriptor::buildSignlessLaplacian(const Eigen::MatrixXd& degree, const Eigen::MatrixXd& adjacency) const {
    return degree + adjacency;
}

Eigen::MatrixXd EigenDescriptor::buildWeightedAdjacency(const RDKit::ROMol* mol) const {
    int numAtoms = mol->getNumAtoms();
    Eigen::MatrixXd weightedAdjacency = Eigen::MatrixXd::Zero(numAtoms, numAtoms);
    for (const RDKit::Bond* bond : mol->bonds()) {
        int startIdx = bond->getBeginAtomIdx();
        int endIdx = bond->getEndAtomIdx();
        double weight = 1.0;
        switch (bond->getBondType()) {
            case RDKit::Bond::SINGLE:   weight = 1.0; break;
            case RDKit::Bond::DOUBLE:   weight = 2.0; break;
            case RDKit::Bond::TRIPLE:   weight = 3.0; break;
            case RDKit::Bond::AROMATIC: weight = 1.5; break;
            default: weight = 1.0; // Fallback for other types
        }
        weightedAdjacency(startIdx, endIdx) = weight;
        weightedAdjacency(endIdx, startIdx) = weight;
    }
    return weightedAdjacency;
}

Eigen::MatrixXd EigenDescriptor::buildDistanceMatrix(const RDKit::ROMol* mol) const {
     int n = mol->getNumAtoms();
     if (n == 0) return Eigen::MatrixXd(0, 0);

     Eigen::MatrixXd dist = Eigen::MatrixXd::Constant(n, n, std::numeric_limits<double>::infinity());
     const double* distMatRaw = RDKit::MolOps::getDistanceMat(*mol, false, false, true); // useBO=false, useAtomWts=false, force=true

     for (int i = 0; i < n; ++i) {
         dist(i, i) = 0.0;
         for (int j = i + 1; j < n; ++j) {
             // Check bounds for safety, although RDKit should provide n*n
             if (i * n + j >= n * n) continue;
             double d = distMatRaw[i * n + j];
             if (d >= 0) { // RDKit returns -1 for unconnected pairs
                 dist(i, j) = d;
                 dist(j, i) = d;
             }
         }
     }
     return dist;
}

Eigen::MatrixXd EigenDescriptor::buildAtomicNumberWeightedAdjacency(const RDKit::ROMol* mol) const {
    int numAtoms = mol->getNumAtoms();
    Eigen::MatrixXd weightedAdj = Eigen::MatrixXd::Zero(numAtoms, numAtoms);
    for (const RDKit::Bond* bond : mol->bonds()) {
        int startIdx = bond->getBeginAtomIdx();
        int endIdx = bond->getEndAtomIdx();
        double weight = static_cast<double>(bond->getBeginAtom()->getAtomicNum()) *
                        static_cast<double>(bond->getEndAtom()->getAtomicNum());
        weightedAdj(startIdx, endIdx) = weight;
        weightedAdj(endIdx, startIdx) = weight;
    }
    return weightedAdj;
}

Eigen::MatrixXd EigenDescriptor::buildAtomicNumberWeightedLaplacian(const RDKit::ROMol* mol) const {
    Eigen::MatrixXd weightedAdj = buildAtomicNumberWeightedAdjacency(mol);
    int size = weightedAdj.rows();
    Eigen::MatrixXd degree = Eigen::MatrixXd::Zero(size, size);
    for (int i = 0; i < size; ++i) {
        degree(i, i) = weightedAdj.row(i).sum(); // Use weighted degree
    }
    return degree - weightedAdj;
}

Eigen::MatrixXd EigenDescriptor::buildNormalizedAdjacency(const Eigen::MatrixXd& degree, const Eigen::MatrixXd& adjacency) const {
    int size = adjacency.rows();
    Eigen::MatrixXd dInvSqrt = Eigen::MatrixXd::Zero(size, size);
    for (int i = 0; i < size; ++i) {
         if (degree(i, i) > 1e-10) {
            dInvSqrt(i, i) = 1.0 / std::sqrt(degree(i, i));
         }
    }
    return dInvSqrt * adjacency * dInvSqrt;
}


// --- Compute Eigenvalues/vectors/SVD - Implementations belong here ---
Eigen::VectorXd EigenDescriptor::computeEigenvalues(const Eigen::MatrixXd& matrix) const {
    if (matrix.rows() == 0) return Eigen::VectorXd(0);
    if (matrix.isApprox(matrix.transpose())) {
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(matrix);
        return solver.eigenvalues();
    } else {
         Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(matrix); // Assume symmetric for descriptor context
         return solver.eigenvalues();
    }
}

Eigen::VectorXd EigenDescriptor::computeSingularValues(const Eigen::MatrixXd& matrix) const {
    if (matrix.rows() == 0 || matrix.cols() == 0) return Eigen::VectorXd(0);
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(matrix);
    return svd.singularValues();
}

Eigen::VectorXd EigenDescriptor::computeDominantEigenvector(const Eigen::MatrixXd& matrix) const {
     if (matrix.rows() == 0) return Eigen::VectorXd(0);
     if (matrix.isApprox(matrix.transpose())) {
         Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(matrix);
         if (solver.info() != Eigen::Success) return Eigen::VectorXd(0); // Check solver success
         if (solver.eigenvalues().size() == 0) return Eigen::VectorXd(0); // Check eigenvalues available
         Eigen::Index maxIndex;
         solver.eigenvalues().cwiseAbs().maxCoeff(&maxIndex);
         if (maxIndex >= solver.eigenvectors().cols()) return Eigen::VectorXd(0); // Check index validity
         return solver.eigenvectors().col(maxIndex);
     } else {
         return Eigen::VectorXd::Zero(matrix.rows()); // Assume symmetric context
     }
}

// --- Statistical Helpers - Implementations belong here ---
double EigenDescriptor::computeVariance(const Eigen::VectorXd& values) const {
    if (values.size() == 0) return 0.0;
    double mean = values.mean();
    return (values.array() - mean).square().mean();
}

double EigenDescriptor::computeSkewness(const Eigen::VectorXd& values) const {
    if (values.size() < 3) return 0.0; // Skewness undefined for < 3 points
    double mean = values.mean();
    double variance = computeVariance(values);
    if (variance < 1e-10) return 0.0;
    double stddev = std::sqrt(variance);
    if (stddev < 1e-10) return 0.0; // Avoid division by zero
    double m3 = (values.array() - mean).pow(3).mean();
    return m3 / (stddev * stddev * stddev);
}

double EigenDescriptor::computeKurtosis(const Eigen::VectorXd& values) const {
    if (values.size() < 4) return 0.0; // Kurtosis undefined for < 4 points
    double mean = values.mean();
    double variance = computeVariance(values);
    if (variance < 1e-10) return 0.0;
    double var_sq = variance * variance;
    if (var_sq < 1e-10) return 0.0; // Avoid division by zero
    double m4 = (values.array() - mean).pow(4).mean();
    // Return excess kurtosis
    return m4 / var_sq - 3.0;
}

// --- Implementations for all calculate() methods ---

// Adjacency based
AdjacencyNonZeroEntries::AdjacencyNonZeroEntries() : EigenDescriptor("AdjNonZero", "Number of non-zero entries (edges)") {}
std::variant<double, int, std::string> AdjacencyNonZeroEntries::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    return static_cast<int>(rdkMol->getNumBonds() * 2);
}

AdjacencyMatrixTrace::AdjacencyMatrixTrace() : EigenDescriptor("AdjTrace", "Trace of adjacency matrix") {}
std::variant<double, int, std::string> AdjacencyMatrixTrace::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    return 0.0; // Trace of simple graph adjacency is 0
}

AdjacencyFrobeniusNorm::AdjacencyFrobeniusNorm() : EigenDescriptor("AdjFrobNorm", "Frobenius norm of adjacency matrix") {}
std::variant<double, int, std::string> AdjacencyFrobeniusNorm::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    if (rdkMol->getNumAtoms() == 0) return 0.0;
    const Eigen::MatrixXd& adjacency = getAdjacencyMatrix(rdkMol);
    return adjacency.norm();
}

AdjacencySpectralRadius::AdjacencySpectralRadius() : EigenDescriptor("AdjSpecRad", "Spectral radius (largest abs eigenvalue) of adjacency matrix") {}
std::variant<double, int, std::string> AdjacencySpectralRadius::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    if (rdkMol->getNumAtoms() == 0) return 0.0;
    const Eigen::MatrixXd& adjacency = getAdjacencyMatrix(rdkMol);
    const Eigen::VectorXd& eigenvalues = getEigenvalues(adjacency);
    if (eigenvalues.size() == 0) return 0.0; // Handle case where eigenvalues couldn't be computed
    return eigenvalues.cwiseAbs().maxCoeff();
}

AdjacencySmallestEigenvalue::AdjacencySmallestEigenvalue() : EigenDescriptor("AdjMinEig", "Smallest eigenvalue of adjacency matrix") {}
std::variant<double, int, std::string> AdjacencySmallestEigenvalue::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    if (rdkMol->getNumAtoms() == 0) return 0.0;
    const Eigen::MatrixXd& adjacency = getAdjacencyMatrix(rdkMol);
    const Eigen::VectorXd& eigenvalues = getEigenvalues(adjacency);
    if (eigenvalues.size() == 0) return 0.0;
    return eigenvalues.minCoeff();
}

AdjacencySumEigenvalues::AdjacencySumEigenvalues() : EigenDescriptor("AdjEigSum", "Sum of eigenvalues of adjacency matrix (Trace)") {}
std::variant<double, int, std::string> AdjacencySumEigenvalues::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    return 0.0; // Trace is 0
}

AdjacencySumSquaresEigenvalues::AdjacencySumSquaresEigenvalues() : EigenDescriptor("AdjEigSumSq", "Sum of squares of eigenvalues of adjacency matrix (Frobenius norm squared)") {}
std::variant<double, int, std::string> AdjacencySumSquaresEigenvalues::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    if (rdkMol->getNumAtoms() == 0) return 0.0;
    return static_cast<double>(rdkMol->getNumBonds() * 2);
}

AdjacencyEigenvalueVariance::AdjacencyEigenvalueVariance() : EigenDescriptor("AdjEigVar", "Variance of eigenvalues of adjacency matrix") {}
std::variant<double, int, std::string> AdjacencyEigenvalueVariance::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    if (rdkMol->getNumAtoms() == 0) return 0.0;
    const Eigen::MatrixXd& adjacency = getAdjacencyMatrix(rdkMol);
    const Eigen::VectorXd& eigenvalues = getEigenvalues(adjacency);
    if (eigenvalues.size() == 0) return 0.0;
    return computeVariance(eigenvalues);
}

AdjacencyEigenvalueSkewness::AdjacencyEigenvalueSkewness() : EigenDescriptor("AdjEigSkew", "Skewness of eigenvalues of adjacency matrix") {}
std::variant<double, int, std::string> AdjacencyEigenvalueSkewness::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    if (rdkMol->getNumAtoms() < 3) return 0.0;
    const Eigen::MatrixXd& adjacency = getAdjacencyMatrix(rdkMol);
    const Eigen::VectorXd& eigenvalues = getEigenvalues(adjacency);
    if (eigenvalues.size() < 3) return 0.0;
    return computeSkewness(eigenvalues);
}

AdjacencyEigenvalueKurtosis::AdjacencyEigenvalueKurtosis() : EigenDescriptor("AdjEigKurt", "Kurtosis of eigenvalues of adjacency matrix") {}
std::variant<double, int, std::string> AdjacencyEigenvalueKurtosis::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    if (rdkMol->getNumAtoms() < 4) return 0.0;
    const Eigen::MatrixXd& adjacency = getAdjacencyMatrix(rdkMol);
    const Eigen::VectorXd& eigenvalues = getEigenvalues(adjacency);
    if (eigenvalues.size() < 4) return 0.0;
    return computeKurtosis(eigenvalues);
}

AdjacencyPositiveEigenvalues::AdjacencyPositiveEigenvalues() : EigenDescriptor("AdjPosEig", "Number of positive eigenvalues of adjacency matrix") {}
std::variant<double, int, std::string> AdjacencyPositiveEigenvalues::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    if (rdkMol->getNumAtoms() == 0) return 0;
    const Eigen::MatrixXd& adjacency = getAdjacencyMatrix(rdkMol);
    const Eigen::VectorXd& eigenvalues = getEigenvalues(adjacency);
    if (eigenvalues.size() == 0) return 0;
    return static_cast<int>((eigenvalues.array() > 1e-10).count());
}

AdjacencyNegativeEigenvalues::AdjacencyNegativeEigenvalues() : EigenDescriptor("AdjNegEig", "Number of negative eigenvalues of adjacency matrix") {}
std::variant<double, int, std::string> AdjacencyNegativeEigenvalues::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    if (rdkMol->getNumAtoms() == 0) return 0;
    const Eigen::MatrixXd& adjacency = getAdjacencyMatrix(rdkMol);
    const Eigen::VectorXd& eigenvalues = getEigenvalues(adjacency);
    if (eigenvalues.size() == 0) return 0;
    return static_cast<int>((eigenvalues.array() < -1e-10).count());
}

AdjacencyZeroEigenvalues::AdjacencyZeroEigenvalues() : EigenDescriptor("AdjZeroEig", "Number of zero eigenvalues of adjacency matrix") {}
std::variant<double, int, std::string> AdjacencyZeroEigenvalues::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    if (rdkMol->getNumAtoms() == 0) return 0;
    const Eigen::MatrixXd& adjacency = getAdjacencyMatrix(rdkMol);
    const Eigen::VectorXd& eigenvalues = getEigenvalues(adjacency);
    if (eigenvalues.size() == 0) return 0;
    return static_cast<int>((eigenvalues.array().abs() < 1e-10).count());
}

AdjacencyMaxDegree::AdjacencyMaxDegree() : EigenDescriptor("AdjMaxDeg", "Maximum degree") {}
std::variant<double, int, std::string> AdjacencyMaxDegree::calculate(const Molecule& mol) const {
     if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    if (rdkMol->getNumAtoms() == 0) return 0;
    const Eigen::MatrixXd& degreeMat = getDegreeMatrix(rdkMol);
    if (degreeMat.rows() == 0) return 0;
    return static_cast<int>(degreeMat.diagonal().maxCoeff());
}

AdjacencyMinDegree::AdjacencyMinDegree() : EigenDescriptor("AdjMinDeg", "Minimum degree") {}
std::variant<double, int, std::string> AdjacencyMinDegree::calculate(const Molecule& mol) const {
     if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    if (rdkMol->getNumAtoms() == 0) return 0;
    const Eigen::MatrixXd& degreeMat = getDegreeMatrix(rdkMol);
    if (degreeMat.rows() == 0) return 0;
    Eigen::VectorXd degrees = degreeMat.diagonal();
    if (degrees.size() == 0) return 0;
    // Find min degree > 0, handle case where all degrees are 0
    double minDeg = std::numeric_limits<double>::max();
    bool nonZeroFound = false;
    for(int i=0; i<degrees.size(); ++i) {
        if(degrees[i] > 1e-10) {
            minDeg = std::min(minDeg, degrees[i]);
            nonZeroFound = true;
        }
    }
    return nonZeroFound ? static_cast<int>(minDeg) : 0;
}

AdjacencyMeanDegree::AdjacencyMeanDegree() : EigenDescriptor("AdjMeanDeg", "Mean degree") {}
std::variant<double, int, std::string> AdjacencyMeanDegree::calculate(const Molecule& mol) const {
     if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    if (rdkMol->getNumAtoms() == 0) return 0.0;
    const Eigen::MatrixXd& degreeMat = getDegreeMatrix(rdkMol);
    if (degreeMat.rows() == 0) return 0.0;
    return degreeMat.diagonal().mean();
}

AdjacencyDegreeVariance::AdjacencyDegreeVariance() : EigenDescriptor("AdjDegVar", "Variance of degrees") {}
std::variant<double, int, std::string> AdjacencyDegreeVariance::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    if (rdkMol->getNumAtoms() == 0) return 0.0;
    const Eigen::MatrixXd& degreeMat = getDegreeMatrix(rdkMol);
    if (degreeMat.rows() == 0) return 0.0;
    return computeVariance(degreeMat.diagonal());
}

AdjacencyDegreeStdDev::AdjacencyDegreeStdDev() : EigenDescriptor("AdjDegStd", "Standard deviation of degrees") {}
std::variant<double, int, std::string> AdjacencyDegreeStdDev::calculate(const Molecule& mol) const {
    auto variance_val = AdjacencyDegreeVariance().calculate(mol);
    if (std::holds_alternative<double>(variance_val)) {
        double variance = std::get<double>(variance_val);
        return (variance >= 0.0) ? std::sqrt(variance) : 0.0;
    } else {
        return variance_val; // Propagate error string
    }
}

AdjacencyMatrixRank::AdjacencyMatrixRank() : EigenDescriptor("AdjRank", "Rank of adjacency matrix") {}
std::variant<double, int, std::string> AdjacencyMatrixRank::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    if (rdkMol->getNumAtoms() == 0) return 0;
    const Eigen::MatrixXd& adjacency = getAdjacencyMatrix(rdkMol);
    Eigen::FullPivLU<Eigen::MatrixXd> lu(adjacency);
    return static_cast<int>(lu.rank());
}

GraphEnergy::GraphEnergy() : EigenDescriptor("GraphEnergy", "Energy of the graph (sum of absolute adjacency eigenvalues)") {}
std::variant<double, int, std::string> GraphEnergy::calculate(const Molecule& mol) const {
     if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    if (rdkMol->getNumAtoms() == 0) return 0.0;
    const Eigen::MatrixXd& adjacency = getAdjacencyMatrix(rdkMol);
    const Eigen::VectorXd& eigenvalues = getEigenvalues(adjacency);
    if (eigenvalues.size() == 0) return 0.0;
    return eigenvalues.cwiseAbs().sum();
}


// Laplacian based
LaplacianSpectralRadius::LaplacianSpectralRadius() : EigenDescriptor("LapSpecRad", "Spectral radius of Laplacian matrix") {}
std::variant<double, int, std::string> LaplacianSpectralRadius::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    if (rdkMol->getNumAtoms() == 0) return 0.0;
    const Eigen::MatrixXd& laplacian = getLaplacianMatrix(rdkMol);
    const Eigen::VectorXd& eigenvalues = getEigenvalues(laplacian);
    if (eigenvalues.size() == 0) return 0.0;
    return eigenvalues.maxCoeff();
}

LaplacianAlgebraicConnectivity::LaplacianAlgebraicConnectivity() : EigenDescriptor("LapAlgConn", "Algebraic connectivity (second-smallest eigenvalue of Laplacian)") {}
std::variant<double, int, std::string> LaplacianAlgebraicConnectivity::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    if (rdkMol->getNumAtoms() <= 1) return 0.0;
    const Eigen::MatrixXd& laplacian = getLaplacianMatrix(rdkMol);
    const Eigen::VectorXd& eigenvalues = getEigenvalues(laplacian);
    if (eigenvalues.size() < 2) return 0.0;
    return eigenvalues(1);
}

LaplacianZeroEigenvalues::LaplacianZeroEigenvalues() : EigenDescriptor("LapZeroEig", "Number of zero eigenvalues (connected components)") {}
std::variant<double, int, std::string> LaplacianZeroEigenvalues::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    if (rdkMol->getNumAtoms() == 0) return 0;
    const Eigen::MatrixXd& laplacian = getLaplacianMatrix(rdkMol);
    const Eigen::VectorXd& eigenvalues = getEigenvalues(laplacian);
    if (eigenvalues.size() == 0) return 0;
    return static_cast<int>((eigenvalues.array().abs() < 1e-10).count());
}

LaplacianEnergy::LaplacianEnergy() : EigenDescriptor("LapEnergy", "Laplacian energy (sum |lambda_i - mean_degree|)") {}
std::variant<double, int, std::string> LaplacianEnergy::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    int n = rdkMol->getNumAtoms();
    if (n == 0) return 0.0;
    const Eigen::MatrixXd& laplacian = getLaplacianMatrix(rdkMol);
    const Eigen::VectorXd& eigenvalues = getEigenvalues(laplacian);
    if (eigenvalues.size() == 0) return 0.0;
    const Eigen::MatrixXd& degreeMat = getDegreeMatrix(rdkMol);
    if (degreeMat.rows() == 0) return 0.0;
    double mean_degree = degreeMat.diagonal().sum() / n;
    return (eigenvalues.array() - mean_degree).abs().sum();
}

LaplacianMatrixTrace::LaplacianMatrixTrace() : EigenDescriptor("LapTrace", "Trace of Laplacian matrix (sum of degrees)") {}
std::variant<double, int, std::string> LaplacianMatrixTrace::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    if (rdkMol->getNumAtoms() == 0) return 0.0;
    const Eigen::MatrixXd& laplacian = getLaplacianMatrix(rdkMol);
    return laplacian.trace();
}

LaplacianMatrixDeterminant::LaplacianMatrixDeterminant() : EigenDescriptor("LapDet", "Determinant of Laplacian matrix (related to spanning trees)") {}
std::variant<double, int, std::string> LaplacianMatrixDeterminant::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    if (rdkMol->getNumAtoms() == 0) return 1.0;
    return 0.0; // Determinant is always 0 for N > 0
}

LaplacianTotalEffectiveResistance::LaplacianTotalEffectiveResistance() : EigenDescriptor("LapTotEffRes", "Total effective resistance (sum 1/lambda_i for non-zero)") {}
std::variant<double, int, std::string> LaplacianTotalEffectiveResistance::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    if (rdkMol->getNumAtoms() <= 1) return 0.0;
    const Eigen::MatrixXd& laplacian = getLaplacianMatrix(rdkMol);
    const Eigen::VectorXd& eigenvalues = getEigenvalues(laplacian);
    if (eigenvalues.size() < 2) return 0.0;
    double sum_inv_eig = 0.0;
    for (int i = 1; i < eigenvalues.size(); ++i) {
        if (eigenvalues(i) > 1e-10) {
            sum_inv_eig += 1.0 / eigenvalues(i);
        } else {
            return "Error: Non-positive Laplacian eigenvalue"; // Should not happen for connected graph lambda_1
        }
    }
    return sum_inv_eig;
}

LaplacianKirchhoffIndex::LaplacianKirchhoffIndex() : EigenDescriptor("LapKirchhoff", "Kirchhoff index (N * sum 1/lambda_i for non-zero)") {}
std::variant<double, int, std::string> LaplacianKirchhoffIndex::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    int n = rdkMol->getNumAtoms();
    if (n <= 1) return 0.0;

    auto totalResistanceResult = LaplacianTotalEffectiveResistance().calculate(mol);
    if (std::holds_alternative<double>(totalResistanceResult)) {
        return n * std::get<double>(totalResistanceResult);
    } else {
        return totalResistanceResult; // Propagate error
    }
}

LaplacianEigenvalueVariance::LaplacianEigenvalueVariance() : EigenDescriptor("LapEigVar", "Variance of Laplacian eigenvalues") {}
std::variant<double, int, std::string> LaplacianEigenvalueVariance::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    if (rdkMol->getNumAtoms() == 0) return 0.0;
    const Eigen::MatrixXd& laplacian = getLaplacianMatrix(rdkMol);
    const Eigen::VectorXd& eigenvalues = getEigenvalues(laplacian);
    if (eigenvalues.size() == 0) return 0.0;
    return computeVariance(eigenvalues);
}

LaplacianEigenvalueSkewness::LaplacianEigenvalueSkewness() : EigenDescriptor("LapEigSkew", "Skewness of Laplacian eigenvalues") {}
std::variant<double, int, std::string> LaplacianEigenvalueSkewness::calculate(const Molecule& mol) const {
     if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    if (rdkMol->getNumAtoms() < 3) return 0.0;
    const Eigen::MatrixXd& laplacian = getLaplacianMatrix(rdkMol);
    const Eigen::VectorXd& eigenvalues = getEigenvalues(laplacian);
    if (eigenvalues.size() < 3) return 0.0;
    return computeSkewness(eigenvalues);
}


// Normalized Laplacian based
NormalizedLaplacianSpectralRadius::NormalizedLaplacianSpectralRadius() : EigenDescriptor("NormLapSpecRad", "Spectral radius of normalized Laplacian") {}
std::variant<double, int, std::string> NormalizedLaplacianSpectralRadius::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    if (rdkMol->getNumAtoms() == 0) return 0.0;
    const Eigen::MatrixXd& normLap = getNormalizedLaplacianMatrix(rdkMol);
    const Eigen::VectorXd& eigenvalues = getEigenvalues(normLap);
    if (eigenvalues.size() == 0) return 0.0;
    return eigenvalues.maxCoeff();
}

NormalizedLaplacianSmallestNonzero::NormalizedLaplacianSmallestNonzero() : EigenDescriptor("NormLapMinNonZero", "Smallest non-zero eigenvalue of normalized Laplacian") {}
std::variant<double, int, std::string> NormalizedLaplacianSmallestNonzero::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    if (rdkMol->getNumAtoms() <= 1) return 0.0;
    const Eigen::MatrixXd& normLap = getNormalizedLaplacianMatrix(rdkMol);
    const Eigen::VectorXd& eigenvalues = getEigenvalues(normLap);
    if (eigenvalues.size() < 2) return 0.0;
    for (int i = 1; i < eigenvalues.size(); ++i) {
        if (eigenvalues(i) > 1e-10) {
            return eigenvalues(i);
        }
    }
    return 0.0; // All non-zero eigenvalues were effectively zero (or only one eigenvalue)
}

NormalizedLaplacianLargestEigenvalue::NormalizedLaplacianLargestEigenvalue() : EigenDescriptor("NormLapMaxEig", "Largest eigenvalue of normalized Laplacian") {}
std::variant<double, int, std::string> NormalizedLaplacianLargestEigenvalue::calculate(const Molecule& mol) const {
    return NormalizedLaplacianSpectralRadius().calculate(mol);
}

NormalizedLaplacianEnergy::NormalizedLaplacianEnergy() : EigenDescriptor("NormLapEnergy", "Normalized Laplacian energy sum(|lambda_i - 1|)") {}
std::variant<double, int, std::string> NormalizedLaplacianEnergy::calculate(const Molecule& mol) const {
     if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    if (rdkMol->getNumAtoms() == 0) return 0.0;
    const Eigen::MatrixXd& normLap = getNormalizedLaplacianMatrix(rdkMol);
    const Eigen::VectorXd& eigenvalues = getEigenvalues(normLap);
    if (eigenvalues.size() == 0) return 0.0;
    return (eigenvalues.array() - 1.0).abs().sum();
}

NormalizedLaplacianTrace::NormalizedLaplacianTrace() : EigenDescriptor("NormLapTrace", "Trace of normalized Laplacian (number of atoms)") {}
std::variant<double, int, std::string> NormalizedLaplacianTrace::calculate(const Molecule& mol) const {
     if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    return static_cast<double>(rdkMol->getNumAtoms());
}


// Degree matrix based (delegates ensure consistency)
DegreeMatrixMaxDegree::DegreeMatrixMaxDegree() : EigenDescriptor("DegMatMaxDeg", "Max degree from Degree Matrix") {}
std::variant<double, int, std::string> DegreeMatrixMaxDegree::calculate(const Molecule& mol) const {
    return AdjacencyMaxDegree().calculate(mol);
}

DegreeMatrixMinDegree::DegreeMatrixMinDegree() : EigenDescriptor("DegMatMinDeg", "Min degree from Degree Matrix") {}
std::variant<double, int, std::string> DegreeMatrixMinDegree::calculate(const Molecule& mol) const {
    return AdjacencyMinDegree().calculate(mol);
}

DegreeMatrixAvgDegree::DegreeMatrixAvgDegree() : EigenDescriptor("DegMatAvgDeg", "Avg degree from Degree Matrix") {}
std::variant<double, int, std::string> DegreeMatrixAvgDegree::calculate(const Molecule& mol) const {
    return AdjacencyMeanDegree().calculate(mol);
}

DegreeMatrixVariance::DegreeMatrixVariance() : EigenDescriptor("DegMatVar", "Variance of degrees from Degree Matrix") {}
std::variant<double, int, std::string> DegreeMatrixVariance::calculate(const Molecule& mol) const {
    return AdjacencyDegreeVariance().calculate(mol);
}

DegreeMatrixSkewness::DegreeMatrixSkewness() : EigenDescriptor("DegMatSkew", "Skewness of degrees from Degree Matrix") {}
std::variant<double, int, std::string> DegreeMatrixSkewness::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    if (rdkMol->getNumAtoms() < 3) return 0.0;
    const Eigen::MatrixXd& degreeMat = getDegreeMatrix(rdkMol);
    if (degreeMat.rows() < 3) return 0.0;
    return computeSkewness(degreeMat.diagonal());
}

DegreeMatrixEntropy::DegreeMatrixEntropy() : EigenDescriptor("DegMatEnt", "Entropy of degree distribution") {}
std::variant<double, int, std::string> DegreeMatrixEntropy::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    int n = rdkMol->getNumAtoms();
    if (n == 0) return 0.0;
    const Eigen::MatrixXd& degreeMat = getDegreeMatrix(rdkMol);
    if (degreeMat.rows() == 0) return 0.0;
    Eigen::VectorXd degrees = degreeMat.diagonal();
    if (degrees.size() == 0) return 0.0;

    std::map<int, int> counts;
    for (int i = 0; i < n; ++i) {
        counts[static_cast<int>(round(degrees(i)))]++;
    }

    double entropy = 0.0;
    for (auto const& [deg, count] : counts) {
        if (count <= 0) continue;
        double p = static_cast<double>(count) / n;
        if (p > 1e-10) {
            entropy -= p * std::log2(p);
        }
    }
    return entropy;
}


// Other graph properties
GraphIrregularity::GraphIrregularity() : EigenDescriptor("GraphIrreg", "Graph irregularity (sum |deg_i - deg_j|)") {}
std::variant<double, int, std::string> GraphIrregularity::calculate(const Molecule& mol) const {
     if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    int n = rdkMol->getNumAtoms();
    if (n <= 1) return 0.0;
    const Eigen::MatrixXd& degreeMat = getDegreeMatrix(rdkMol);
    if (degreeMat.rows() == 0) return 0.0;
    Eigen::VectorXd degrees = degreeMat.diagonal();
    if (degrees.size() != n) return "Error: Degree vector size mismatch"; // Safety check
    double irregularity = 0.0;
    for(int i=0; i<n; ++i) {
        for (int j=i+1; j<n; ++j) {
            irregularity += std::abs(degrees(i) - degrees(j));
        }
    }
    return irregularity;
}

WienerIndex::WienerIndex() : EigenDescriptor("WienerIdx", "Wiener index (sum of shortest path distances)") {}
std::variant<double, int, std::string> WienerIndex::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    int n = rdkMol->getNumAtoms();
    if (n <= 1) return 0.0;

    const double* distMatRaw = RDKit::MolOps::getDistanceMat(*rdkMol, false, false, true);
    double wienerSum = 0.0;
    bool connected = true;

    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            if (i * n + j >= n * n) return "Error: Distance matrix index out of bounds"; // Safety
            double d = distMatRaw[i * n + j];
            if (d < 0) {
                 connected = false;
                 break;
            }
             wienerSum += d;
        }
         if(!connected) break;
    }

    if (!connected) {
         return "Error: Molecule has disconnected components (Wiener)";
    }
    return wienerSum;
}

EstradaIndex::EstradaIndex() : EigenDescriptor("EstradaIdx", "Estrada index (sum exp(adj_eigenvalues))") {}
std::variant<double, int, std::string> EstradaIndex::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    if (rdkMol->getNumAtoms() == 0) return 0.0;
    const Eigen::MatrixXd& adjacency = getAdjacencyMatrix(rdkMol);
    const Eigen::VectorXd& eigenvalues = getEigenvalues(adjacency);
    if (eigenvalues.size() == 0) return 0.0;
    return eigenvalues.array().exp().sum();
}

NumberSpanningTrees::NumberSpanningTrees() : EigenDescriptor("NumSpanTrees", "Number of spanning trees (Laplacian minor determinant)") {}
std::variant<double, int, std::string> NumberSpanningTrees::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    int n = rdkMol->getNumAtoms();
    if (n == 0) return 0;
    if (n == 1) return 1;

    auto zeroEigResult = LaplacianZeroEigenvalues().calculate(mol);
    if (std::holds_alternative<int>(zeroEigResult)) {
        if (std::get<int>(zeroEigResult) != 1) {
             return 0; // Disconnected graph has 0 spanning trees
        }
    } else {
         return zeroEigResult; // Propagate error
    }

    const Eigen::MatrixXd& laplacian = getLaplacianMatrix(rdkMol);
    if (n <= 1 || laplacian.rows() < n-1 || laplacian.cols() < n-1) {
         // This should not happen if n > 1 and connected, but check anyway
         return "Error: Laplacian matrix size invalid for minor calculation";
    }
    Eigen::MatrixXd minor = laplacian.topLeftCorner(n - 1, n - 1);
    double det = minor.determinant();
    return static_cast<int>(std::round(std::abs(det)));
}

GraphEccentricity::GraphEccentricity() : EigenDescriptor("GraphEcc", "Maximum eccentricity (max shortest path from a vertex)") {}
std::variant<double, int, std::string> GraphEccentricity::calculate(const Molecule& mol) const {
     if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    int n = rdkMol->getNumAtoms();
    if (n == 0) return 0;
    if (n == 1) return 0;

    const double* distMatRaw = RDKit::MolOps::getDistanceMat(*rdkMol, false, false, true);
    double maxEccentricity = 0.0;

    for (int i = 0; i < n; ++i) {
        double currentEccentricity = 0.0;
        for (int j = 0; j < n; ++j) {
             if (i * n + j >= n * n) return "Error: Distance matrix index out of bounds (Ecc)"; // Safety
             double d = distMatRaw[i * n + j];
             if (d < 0) {
                 return "Error: Molecule has disconnected components (Ecc)";
             }
             currentEccentricity = std::max(currentEccentricity, d);
        }
        maxEccentricity = std::max(maxEccentricity, currentEccentricity);
    }
    return static_cast<int>(maxEccentricity);
}

SpectralGap::SpectralGap() : EigenDescriptor("SpecGap", "Adjacency spectral gap (lambda_1 - lambda_2)") {}
std::variant<double, int, std::string> SpectralGap::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    if (rdkMol->getNumAtoms() <= 1) return 0.0;
    const Eigen::MatrixXd& adjacency = getAdjacencyMatrix(rdkMol);
    const Eigen::VectorXd& eigenvalues = getEigenvalues(adjacency);
    int n = eigenvalues.size();
    if (n < 2) return 0.0;
    return eigenvalues(n - 1) - eigenvalues(n - 2); // Assumes sorted ascending by Eigen solver
}

TraceMatrixPower2::TraceMatrixPower2() : EigenDescriptor("TraceMatPow2", "Trace(A^2) (twice number of edges)") {}
std::variant<double, int, std::string> TraceMatrixPower2::calculate(const Molecule& mol) const {
     if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    return static_cast<double>(rdkMol->getNumBonds() * 2);
}

NumberTriangles::NumberTriangles() : EigenDescriptor("NumTriangles", "Number of triangles (Trace(A^3) / 6)") {}
std::variant<double, int, std::string> NumberTriangles::calculate(const Molecule& mol) const {
     if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    if (rdkMol->getNumAtoms() < 3) return 0.0;
    const Eigen::MatrixXd& adj3 = getAdjPow3(rdkMol);
    return adj3.trace() / 6.0;
}

GraphDiameter::GraphDiameter() : EigenDescriptor("GraphDiam", "Graph diameter (max shortest path)") {}
std::variant<double, int, std::string> GraphDiameter::calculate(const Molecule& mol) const {
    // Diameter is the maximum eccentricity
    auto eccentricity_val = GraphEccentricity().calculate(mol);
    if (std::holds_alternative<int>(eccentricity_val)) {
        return std::get<int>(eccentricity_val);
    } else {
        return eccentricity_val; // Propagate error string
    }
}


// Higher-order adjacency powers
NumberOf2Walks::NumberOf2Walks() : EigenDescriptor("Num2Walks", "Number of 2-walks (Tr(A^2))") {}
std::variant<double, int, std::string> NumberOf2Walks::calculate(const Molecule& mol) const {
    return TraceMatrixPower2().calculate(mol); // Delegate
}

NumberOf3Walks::NumberOf3Walks() : EigenDescriptor("Num3Walks", "Number of 3-walks (Tr(A^3))") {}
std::variant<double, int, std::string> NumberOf3Walks::calculate(const Molecule& mol) const {
     if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    if (rdkMol->getNumAtoms() == 0) return 0.0;
    const Eigen::MatrixXd& adj3 = getAdjPow3(rdkMol);
    return adj3.trace();
}

NumberOf4Walks::NumberOf4Walks() : EigenDescriptor("Num4Walks", "Number of 4-walks (Tr(A^4))") {}
std::variant<double, int, std::string> NumberOf4Walks::calculate(const Molecule& mol) const {
     if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    if (rdkMol->getNumAtoms() == 0) return 0.0;
    const Eigen::MatrixXd& adj4 = getAdjPow4(rdkMol);
    return adj4.trace();
}

MeanClosed3WalksPerNode::MeanClosed3WalksPerNode() : EigenDescriptor("MeanClosed3Walks", "Mean number of closed 3-walks per node (Tr(A^3)/N)") {}
std::variant<double, int, std::string> MeanClosed3WalksPerNode::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    int n = rdkMol->getNumAtoms();
    if (n == 0) return 0.0;
    auto trace3_val = NumberOf3Walks().calculate(mol);
    if (std::holds_alternative<double>(trace3_val)) {
         return std::get<double>(trace3_val) / n;
    } else {
         return trace3_val; // Propagate error or int
    }
}

MeanClosed4WalksPerNode::MeanClosed4WalksPerNode() : EigenDescriptor("MeanClosed4Walks", "Mean number of closed 4-walks per node (Tr(A^4)/N)") {}
std::variant<double, int, std::string> MeanClosed4WalksPerNode::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    int n = rdkMol->getNumAtoms();
    if (n == 0) return 0.0;
     auto trace4_val = NumberOf4Walks().calculate(mol);
     if (std::holds_alternative<double>(trace4_val)) {
          return std::get<double>(trace4_val) / n;
     } else {
          return trace4_val; // Propagate error or int
     }
}

Walk2Energy::Walk2Energy() : EigenDescriptor("Walk2Energy", "Sum of absolute eigenvalues of A^2") {}
std::variant<double, int, std::string> Walk2Energy::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    if (rdkMol->getNumAtoms() == 0) return 0.0;
    const Eigen::MatrixXd& adj2 = getAdjPow2(rdkMol);
    const Eigen::VectorXd& eigenvalues = getEigenvalues(adj2); // A^2 is symmetric PSD, eigs >= 0
    if (eigenvalues.size() == 0) return 0.0;
    return eigenvalues.sum(); // sum(|lambda_i|) = sum(lambda_i) = Trace(A^2)
}

Walk3Energy::Walk3Energy() : EigenDescriptor("Walk3Energy", "Sum of absolute eigenvalues of A^3") {}
std::variant<double, int, std::string> Walk3Energy::calculate(const Molecule& mol) const {
     if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    if (rdkMol->getNumAtoms() == 0) return 0.0;
    const Eigen::MatrixXd& adj3 = getAdjPow3(rdkMol);
    const Eigen::VectorXd& eigenvalues = getEigenvalues(adj3); // A^3 is symmetric
    if (eigenvalues.size() == 0) return 0.0;
    return eigenvalues.cwiseAbs().sum();
}

Walk4Energy::Walk4Energy() : EigenDescriptor("Walk4Energy", "Sum of absolute eigenvalues of A^4") {}
std::variant<double, int, std::string> Walk4Energy::calculate(const Molecule& mol) const {
     if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    if (rdkMol->getNumAtoms() == 0) return 0.0;
    const Eigen::MatrixXd& adj4 = getAdjPow4(rdkMol);
    const Eigen::VectorXd& eigenvalues = getEigenvalues(adj4); // A^4 is symmetric PSD, eigs >= 0
    if (eigenvalues.size() == 0) return 0.0;
    return eigenvalues.sum(); // sum(|lambda_i|) = sum(lambda_i) = Trace(A^4)
}

GraphIrregularityWalkCount::GraphIrregularityWalkCount() : EigenDescriptor("GraphIrregWalk", "Irregularity based on diagonal of A^2") {}
std::variant<double, int, std::string> GraphIrregularityWalkCount::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    int n = rdkMol->getNumAtoms();
    if (n <= 1) return 0.0;
    const Eigen::MatrixXd& adj2 = getAdjPow2(rdkMol);
    if (adj2.rows() != n) return "Error: A^2 matrix size mismatch"; // Safety
    Eigen::VectorXd diag = adj2.diagonal(); // diag(A^2)_i = degree_i
    if (diag.size() != n) return "Error: A^2 diagonal size mismatch"; // Safety
    double irregularity = 0.0;
    for(int i=0; i<n; ++i) {
        for (int j=i+1; j<n; ++j) {
            irregularity += std::abs(diag(i) - diag(j));
        }
    }
    return irregularity;
}

TrianglesToPathsRatio::TrianglesToPathsRatio() : EigenDescriptor("TriPathRatio", "Ratio Tr(A^3)/Tr(A^2)") {}
std::variant<double, int, std::string> TrianglesToPathsRatio::calculate(const Molecule& mol) const {
     if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    if (rdkMol->getNumAtoms() < 3) return 0.0;

    auto trace3_val = NumberOf3Walks().calculate(mol);
    auto trace2_val = NumberOf2Walks().calculate(mol);

    if (!std::holds_alternative<double>(trace3_val)) return trace3_val; // Propagate error
    if (!std::holds_alternative<double>(trace2_val)) return trace2_val; // Propagate error

    double trace3 = std::get<double>(trace3_val);
    double trace2 = std::get<double>(trace2_val);

    return (trace2 > 1e-10) ? (trace3 / trace2) : 0.0;
}


// SVD based
MaxSingularValue::MaxSingularValue() : EigenDescriptor("MaxSingVal", "Maximum singular value of adjacency matrix") {}
std::variant<double, int, std::string> MaxSingularValue::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    if (rdkMol->getNumAtoms() == 0) return 0.0;
    const Eigen::MatrixXd& adjacency = getAdjacencyMatrix(rdkMol);
    const Eigen::VectorXd& singularValues = getSingularValues(adjacency); // Assumes sorted descending
    return singularValues.size() > 0 ? singularValues(0) : 0.0;
}

MinNonZeroSingularValue::MinNonZeroSingularValue() : EigenDescriptor("MinNonZeroSingVal", "Minimum non-zero singular value of adjacency matrix") {}
std::variant<double, int, std::string> MinNonZeroSingularValue::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    if (rdkMol->getNumAtoms() == 0) return 0.0;
    const Eigen::MatrixXd& adjacency = getAdjacencyMatrix(rdkMol);
    const Eigen::VectorXd& singularValues = getSingularValues(adjacency);
    if (singularValues.size() == 0) return 0.0;
    double min_sv = 0.0;
    for (int i = singularValues.size() - 1; i >= 0; --i) {
        if (singularValues(i) > 1e-10) {
            min_sv = singularValues(i);
            break;
        }
    }
    return min_sv;
}

ConditionNumber::ConditionNumber() : EigenDescriptor("CondNumber", "Condition number (max_sv / min_non_zero_sv)") {}
std::variant<double, int, std::string> ConditionNumber::calculate(const Molecule& mol) const {
    auto max_sv_val = MaxSingularValue().calculate(mol);
    auto min_sv_val = MinNonZeroSingularValue().calculate(mol);

    if (!std::holds_alternative<double>(max_sv_val)) return max_sv_val; // Propagate error
    if (!std::holds_alternative<double>(min_sv_val)) return min_sv_val; // Propagate error

    double max_sv = std::get<double>(max_sv_val);
    double min_sv = std::get<double>(min_sv_val);

    return (min_sv > 1e-10) ? (max_sv / min_sv) : std::numeric_limits<double>::infinity();
}

SumSingularValues::SumSingularValues() : EigenDescriptor("SumSingVal", "Sum of singular values (Nuclear Norm)") {}
std::variant<double, int, std::string> SumSingularValues::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    if (rdkMol->getNumAtoms() == 0) return 0.0;
    const Eigen::MatrixXd& adjacency = getAdjacencyMatrix(rdkMol);
    const Eigen::VectorXd& singularValues = getSingularValues(adjacency);
    return singularValues.sum();
}

FrobeniusNormSingularValues::FrobeniusNormSingularValues() : EigenDescriptor("FrobNormSingVal", "Frobenius norm from singular values (sqrt(sum(sv^2)))") {}
std::variant<double, int, std::string> FrobeniusNormSingularValues::calculate(const Molecule& mol) const {
    return AdjacencyFrobeniusNorm().calculate(mol);
}

SingularValueEntropy::SingularValueEntropy() : EigenDescriptor("SingValEntropy", "Entropy of normalized singular values") {}
std::variant<double, int, std::string> SingularValueEntropy::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    if (rdkMol->getNumAtoms() == 0) return 0.0;
    const Eigen::MatrixXd& adjacency = getAdjacencyMatrix(rdkMol);
    const Eigen::VectorXd& singularValues = getSingularValues(adjacency);
    if (singularValues.size() == 0) return 0.0;
    double sum_sv = singularValues.sum();
    if (sum_sv < 1e-10) return 0.0;

    double entropy = 0.0;
    for(int i=0; i<singularValues.size(); ++i) {
        if (singularValues(i) > 1e-10) {
            double p = singularValues(i) / sum_sv;
            entropy -= p * std::log2(p);
        }
    }
    return entropy;
}

SingularValueVariance::SingularValueVariance() : EigenDescriptor("SingValVar", "Variance of singular values") {}
std::variant<double, int, std::string> SingularValueVariance::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    if (rdkMol->getNumAtoms() == 0) return 0.0;
    const Eigen::MatrixXd& adjacency = getAdjacencyMatrix(rdkMol);
    const Eigen::VectorXd& singularValues = getSingularValues(adjacency);
    if (singularValues.size() == 0) return 0.0;
    return computeVariance(singularValues);
}

SingularValueSkewness::SingularValueSkewness() : EigenDescriptor("SingValSkew", "Skewness of singular values") {}
std::variant<double, int, std::string> SingularValueSkewness::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    if (rdkMol->getNumAtoms() < 3) return 0.0;
    const Eigen::MatrixXd& adjacency = getAdjacencyMatrix(rdkMol);
    const Eigen::VectorXd& singularValues = getSingularValues(adjacency);
    if (singularValues.size() < 3) return 0.0;
    return computeSkewness(singularValues);
}

SpectralEffectiveRank::SpectralEffectiveRank() : EigenDescriptor("SpecEffRank", "Spectral effective rank (exp(entropy of normalized sv))") {}
std::variant<double, int, std::string> SpectralEffectiveRank::calculate(const Molecule& mol) const {
    auto entropy_val = SingularValueEntropy().calculate(mol);
    if (std::holds_alternative<double>(entropy_val)) {
        // Using log2 entropy S = -sum(p*log2(p)), effective rank is 2^S
        return std::pow(2.0, std::get<double>(entropy_val));
    } else {
        return entropy_val; // Propagate error
    }
}

NuclearNorm::NuclearNorm() : EigenDescriptor("NuclearNorm", "Nuclear norm (Sum of singular values)") {}
std::variant<double, int, std::string> NuclearNorm::calculate(const Molecule& mol) const {
    return SumSingularValues().calculate(mol); // Delegate
}


// Normalized Adjacency based
NormalizedAdjacencySpectralRadius::NormalizedAdjacencySpectralRadius() : EigenDescriptor("NormAdjSpecRad", "Spectral radius of normalized adjacency matrix") {}
std::variant<double, int, std::string> NormalizedAdjacencySpectralRadius::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    if (rdkMol->getNumAtoms() == 0) return 0.0;
    const Eigen::MatrixXd& normAdj = getNormalizedAdjacency(rdkMol);
    const Eigen::VectorXd& eigenvalues = getEigenvalues(normAdj);
    if (eigenvalues.size() == 0) return 0.0;
    return eigenvalues.cwiseAbs().maxCoeff();
}

NormalizedEigenvalueGap::NormalizedEigenvalueGap() : EigenDescriptor("NormEigGap", "Gap between largest two abs normalized adjacency eigenvalues") {}
std::variant<double, int, std::string> NormalizedEigenvalueGap::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    if (rdkMol->getNumAtoms() <= 1) return 0.0;
    const Eigen::MatrixXd& normAdj = getNormalizedAdjacency(rdkMol);
    const Eigen::VectorXd& eigenvalues = getEigenvalues(normAdj);
    if (eigenvalues.size() < 2) return 0.0;
    std::vector<double> abs_eigs;
    for (int i = 0; i < eigenvalues.size(); ++i) abs_eigs.push_back(std::abs(eigenvalues[i]));
    std::sort(abs_eigs.begin(), abs_eigs.end(), std::greater<double>());
    return abs_eigs[0] - abs_eigs[1];
}

SumNormalizedEigenvalues::SumNormalizedEigenvalues() : EigenDescriptor("SumNormEig", "Sum of normalized adjacency eigenvalues (Trace)") {}
std::variant<double, int, std::string> SumNormalizedEigenvalues::calculate(const Molecule& mol) const {
     if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    if (rdkMol->getNumAtoms() == 0) return 0.0;
    const Eigen::MatrixXd& normAdj = getNormalizedAdjacency(rdkMol);
    const Eigen::VectorXd& eigenvalues = getEigenvalues(normAdj);
    if (eigenvalues.size() == 0) return 0.0;
    return eigenvalues.sum();
}

VarianceNormalizedEigenvalues::VarianceNormalizedEigenvalues() : EigenDescriptor("VarNormEig", "Variance of normalized adjacency eigenvalues") {}
std::variant<double, int, std::string> VarianceNormalizedEigenvalues::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    if (rdkMol->getNumAtoms() == 0) return 0.0;
    const Eigen::MatrixXd& normAdj = getNormalizedAdjacency(rdkMol);
    const Eigen::VectorXd& eigenvalues = getEigenvalues(normAdj);
    if (eigenvalues.size() == 0) return 0.0;
    return computeVariance(eigenvalues);
}

CountNormalizedEigenvaluesAboveHalf::CountNormalizedEigenvaluesAboveHalf() : EigenDescriptor("CountNormEigAboveHalf", "Number of normalized adjacency eigenvalues > 0.5") {}
std::variant<double, int, std::string> CountNormalizedEigenvaluesAboveHalf::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    if (rdkMol->getNumAtoms() == 0) return 0;
    const Eigen::MatrixXd& normAdj = getNormalizedAdjacency(rdkMol);
    const Eigen::VectorXd& eigenvalues = getEigenvalues(normAdj);
    if (eigenvalues.size() == 0) return 0;
    return static_cast<int>((eigenvalues.array() > 0.5).count());
}

NormalizedEnergy::NormalizedEnergy() : EigenDescriptor("NormEnergy", "Normalized adjacency energy (sum |lambda_i|)") {}
std::variant<double, int, std::string> NormalizedEnergy::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    if (rdkMol->getNumAtoms() == 0) return 0.0;
    const Eigen::MatrixXd& normAdj = getNormalizedAdjacency(rdkMol);
    const Eigen::VectorXd& eigenvalues = getEigenvalues(normAdj);
    if (eigenvalues.size() == 0) return 0.0;
    return eigenvalues.cwiseAbs().sum();
}

LargestNormalizedEigenvectorCentrality::LargestNormalizedEigenvectorCentrality() : EigenDescriptor("MaxNormEigCent", "Max component of eigenvector for largest abs normAdj eigenvalue") {}
std::variant<double, int, std::string> LargestNormalizedEigenvectorCentrality::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    if (rdkMol->getNumAtoms() == 0) return 0.0;
    const Eigen::MatrixXd& normAdj = getNormalizedAdjacency(rdkMol);
    const Eigen::VectorXd& eigenvector = getDominantEigenvector(normAdj);
    if (eigenvector.size() == 0) return 0.0;
    return eigenvector.cwiseAbs().maxCoeff();
}

AverageNormalizedEigenvectorCentrality::AverageNormalizedEigenvectorCentrality() : EigenDescriptor("AvgNormEigCent", "Mean abs component of eigenvector for largest abs normAdj eigenvalue") {}
std::variant<double, int, std::string> AverageNormalizedEigenvectorCentrality::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    if (rdkMol->getNumAtoms() == 0) return 0.0;
    const Eigen::MatrixXd& normAdj = getNormalizedAdjacency(rdkMol);
    const Eigen::VectorXd& eigenvector = getDominantEigenvector(normAdj);
    if (eigenvector.size() == 0) return 0.0;
    return eigenvector.cwiseAbs().mean();
}

NormalizedAdjacencyMatrixRank::NormalizedAdjacencyMatrixRank() : EigenDescriptor("NormAdjRank", "Rank of normalized adjacency matrix") {}
std::variant<double, int, std::string> NormalizedAdjacencyMatrixRank::calculate(const Molecule& mol) const {
     if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    if (rdkMol->getNumAtoms() == 0) return 0;
    const Eigen::MatrixXd& normAdj = getNormalizedAdjacency(rdkMol);
    Eigen::FullPivLU<Eigen::MatrixXd> lu(normAdj);
    return static_cast<int>(lu.rank());
}

NormalizedAdjacencyEntropy::NormalizedAdjacencyEntropy() : EigenDescriptor("NormAdjEntropy", "Entropy based on abs normalized adjacency eigenvalues") {}
std::variant<double, int, std::string> NormalizedAdjacencyEntropy::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    if (rdkMol->getNumAtoms() == 0) return 0.0;
    const Eigen::MatrixXd& normAdj = getNormalizedAdjacency(rdkMol);
    const Eigen::VectorXd& eigenvalues = getEigenvalues(normAdj);
    if (eigenvalues.size() == 0) return 0.0;

    Eigen::VectorXd abs_eigs = eigenvalues.cwiseAbs();
    double sum_abs_eigs = abs_eigs.sum();
    if (sum_abs_eigs < 1e-10) return 0.0;

    double entropy = 0.0;
    for(int i=0; i<abs_eigs.size(); ++i) {
        if (abs_eigs(i) > 1e-10) {
             double p = abs_eigs(i) / sum_abs_eigs;
             entropy -= p * std::log2(p);
        }
    }
    return entropy;
}


// Signless Laplacian based
SignlessLaplacianSpectralRadius::SignlessLaplacianSpectralRadius() : EigenDescriptor("SignLapSpecRad", "Spectral radius of signless Laplacian") {}
std::variant<double, int, std::string> SignlessLaplacianSpectralRadius::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    if (rdkMol->getNumAtoms() == 0) return 0.0;
    const Eigen::MatrixXd& signlessLap = getSignlessLaplacianMatrix(rdkMol);
    const Eigen::VectorXd& eigenvalues = getEigenvalues(signlessLap);
    if (eigenvalues.size() == 0) return 0.0;
    return eigenvalues.maxCoeff();
}

SignlessLaplacianSmallestEigenvalue::SignlessLaplacianSmallestEigenvalue() : EigenDescriptor("SignLapMinEig", "Smallest eigenvalue of signless Laplacian") {}
std::variant<double, int, std::string> SignlessLaplacianSmallestEigenvalue::calculate(const Molecule& mol) const {
     if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    if (rdkMol->getNumAtoms() == 0) return 0.0;
    const Eigen::MatrixXd& signlessLap = getSignlessLaplacianMatrix(rdkMol);
    const Eigen::VectorXd& eigenvalues = getEigenvalues(signlessLap);
    if (eigenvalues.size() == 0) return 0.0;
    return eigenvalues(0);
}

SignlessLaplacianEnergy::SignlessLaplacianEnergy() : EigenDescriptor("SignLapEnergy", "Energy of signless Laplacian (sum |lambda_i|)") {}
std::variant<double, int, std::string> SignlessLaplacianEnergy::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    if (rdkMol->getNumAtoms() == 0) return 0.0;
    const Eigen::MatrixXd& signlessLap = getSignlessLaplacianMatrix(rdkMol);
    const Eigen::VectorXd& eigenvalues = getEigenvalues(signlessLap);
    if (eigenvalues.size() == 0) return 0.0;
    return eigenvalues.sum(); // Eigenvalues are non-negative
}

SignlessLaplacianTrace::SignlessLaplacianTrace() : EigenDescriptor("SignLapTrace", "Trace of signless Laplacian (sum of degrees)") {}
std::variant<double, int, std::string> SignlessLaplacianTrace::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    if (rdkMol->getNumAtoms() == 0) return 0.0;
    const Eigen::MatrixXd& signlessLap = getSignlessLaplacianMatrix(rdkMol);
    return signlessLap.trace();
}

SignlessLaplacianDeterminant::SignlessLaplacianDeterminant() : EigenDescriptor("SignLapDet", "Determinant of signless Laplacian") {}
std::variant<double, int, std::string> SignlessLaplacianDeterminant::calculate(const Molecule& mol) const {
     if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    if (rdkMol->getNumAtoms() == 0) return 1.0;
    const Eigen::MatrixXd& signlessLap = getSignlessLaplacianMatrix(rdkMol);
    return signlessLap.determinant();
}

SignlessLaplacianZeroEigenvalues::SignlessLaplacianZeroEigenvalues() : EigenDescriptor("SignLapZeroEig", "Number of zero eigenvalues of signless Laplacian") {}
std::variant<double, int, std::string> SignlessLaplacianZeroEigenvalues::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    if (rdkMol->getNumAtoms() == 0) return 0;
    const Eigen::MatrixXd& signlessLap = getSignlessLaplacianMatrix(rdkMol);
    const Eigen::VectorXd& eigenvalues = getEigenvalues(signlessLap);
    if (eigenvalues.size() == 0) return 0;
    return static_cast<int>((eigenvalues.array().abs() < 1e-10).count());
}

SignlessLaplacianPositiveEigenvalues::SignlessLaplacianPositiveEigenvalues() : EigenDescriptor("SignLapPosEig", "Number of positive eigenvalues of signless Laplacian") {}
std::variant<double, int, std::string> SignlessLaplacianPositiveEigenvalues::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    if (rdkMol->getNumAtoms() == 0) return 0;
    const Eigen::MatrixXd& signlessLap = getSignlessLaplacianMatrix(rdkMol);
    const Eigen::VectorXd& eigenvalues = getEigenvalues(signlessLap);
    if (eigenvalues.size() == 0) return 0;
    return static_cast<int>((eigenvalues.array() > 1e-10).count());
}

SignlessLaplacianNegativeEigenvalues::SignlessLaplacianNegativeEigenvalues() : EigenDescriptor("SignLapNegEig", "Number of negative eigenvalues of signless Laplacian (Should be 0)") {}
std::variant<double, int, std::string> SignlessLaplacianNegativeEigenvalues::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    return 0; // Signless Laplacian is positive semi-definite
}

SignlessLaplacianEigenvalueVariance::SignlessLaplacianEigenvalueVariance() : EigenDescriptor("SignLapEigVar", "Signless Laplacian eigenvalue variance") {}
std::variant<double, int, std::string> SignlessLaplacianEigenvalueVariance::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    if (rdkMol->getNumAtoms() == 0) return 0.0;
    const Eigen::MatrixXd& signlessLap = getSignlessLaplacianMatrix(rdkMol);
    const Eigen::VectorXd& eigenvalues = getEigenvalues(signlessLap);
    if (eigenvalues.size() == 0) return 0.0;
    return computeVariance(eigenvalues);
}

SignlessLaplacianEigenvalueSkewness::SignlessLaplacianEigenvalueSkewness() : EigenDescriptor("SignLapEigSkew", "Signless Laplacian eigenvalue skewness") {}
std::variant<double, int, std::string> SignlessLaplacianEigenvalueSkewness::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    if (rdkMol->getNumAtoms() < 3) return 0.0;
    const Eigen::MatrixXd& signlessLap = getSignlessLaplacianMatrix(rdkMol);
    const Eigen::VectorXd& eigenvalues = getEigenvalues(signlessLap);
    if (eigenvalues.size() < 3) return 0.0;
    return computeSkewness(eigenvalues);
}


// Weighted Adjacency based
WeightedAdjacencySpectralRadius::WeightedAdjacencySpectralRadius() : EigenDescriptor("WtAdjSpecRad", "Spectral radius of weighted adjacency matrix") {}
std::variant<double, int, std::string> WeightedAdjacencySpectralRadius::calculate(const Molecule& mol) const {
     if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    if (rdkMol->getNumAtoms() == 0) return 0.0;
    const Eigen::MatrixXd& wtAdj = getWeightedAdjacencyMatrix(rdkMol);
    const Eigen::VectorXd& eigenvalues = getEigenvalues(wtAdj); // WtAdj is symmetric
    if (eigenvalues.size() == 0) return 0.0;
    return eigenvalues.cwiseAbs().maxCoeff();
}


// Miscellaneous / Complex
MeanFirstPassageTime::MeanFirstPassageTime() : EigenDescriptor("MeanFirstPassage", "Mean first passage time (Requires pseudo-inverse)") {}
std::variant<double, int, std::string> MeanFirstPassageTime::calculate([[maybe_unused]] const Molecule& mol) const { 
    return "Not Implemented"; 
}

CommuteTimeDistance::CommuteTimeDistance() : EigenDescriptor("CommuteTime", "Commute time distance (Requires pseudo-inverse)") {}
std::variant<double, int, std::string> CommuteTimeDistance::calculate([[maybe_unused]] const Molecule& mol) const { 
    return "Not Implemented"; 
}

KirchhoffIndexVariance::KirchhoffIndexVariance() : EigenDescriptor("KirchhoffVar", "Kirchhoff index variance (Requires subgraphs/minors)") {}
std::variant<double, int, std::string> KirchhoffIndexVariance::calculate([[maybe_unused]] const Molecule& mol) const { 
    return "Not Implemented"; 
}

EffectiveGraphResistanceDistribution::EffectiveGraphResistanceDistribution() : EigenDescriptor("EffGraphResDist", "Effective graph resistance distribution (Requires pseudo-inverse)") {}
std::variant<double, int, std::string> EffectiveGraphResistanceDistribution::calculate([[maybe_unused]] const Molecule& mol) const { 
    return "Not Implemented"; 
}

LocalClusteringCoefficientDistribution::LocalClusteringCoefficientDistribution() : EigenDescriptor("LocalClustCoefDist", "Variance of local clustering coefficients") {}
std::variant<double, int, std::string> LocalClusteringCoefficientDistribution::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    int n = rdkMol->getNumAtoms();
    if (n < 3) return 0.0;

    const Eigen::MatrixXd& adj = getAdjacencyMatrix(rdkMol);
    const Eigen::MatrixXd& adj3 = getAdjPow3(rdkMol); // Use cached A^3
    if (adj3.rows() != n || adj.rows() != n) return "Error: Matrix size mismatch (Clustering)"; // Safety

    Eigen::VectorXd degrees = adj.rowwise().sum();
    if (degrees.size() != n) return "Error: Degree vector size mismatch (Clustering)"; // Safety

    Eigen::VectorXd coefficients(n);
    for (int i = 0; i < n; ++i) {
        double deg = degrees(i);
        if (deg < 2) {
            coefficients(i) = 0.0;
        } else {
            double triangles = adj3(i, i) / 2.0; // Closed 3-walks originating at i
            double max_possible = deg * (deg - 1.0) / 2.0;
            coefficients(i) = (max_possible > 1e-10) ? (triangles / max_possible) : 0.0;
        }
    }
    return computeVariance(coefficients);
}

GraphRobustnessIndex::GraphRobustnessIndex() : EigenDescriptor("GraphRobustIdx", "Graph robustness (Alg. Connectivity / Max Degree)") {}
std::variant<double, int, std::string> GraphRobustnessIndex::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    if (rdkMol->getNumAtoms() <= 1) return 0.0;

    auto alg_conn_val = LaplacianAlgebraicConnectivity().calculate(mol);
    auto max_deg_val = AdjacencyMaxDegree().calculate(mol);

    if (!std::holds_alternative<double>(alg_conn_val)) return alg_conn_val; // Propagate error
    if (!std::holds_alternative<int>(max_deg_val)) return max_deg_val; // Propagate error

    double alg_conn = std::get<double>(alg_conn_val);
    double max_deg = static_cast<double>(std::get<int>(max_deg_val));

    return (max_deg > 1e-10) ? (alg_conn / max_deg) : 0.0;
}

NormalizedEstradaIndex::NormalizedEstradaIndex() : EigenDescriptor("NormEstradaIdx", "Estrada index from normalized adjacency matrix") {}
std::variant<double, int, std::string> NormalizedEstradaIndex::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    if (rdkMol->getNumAtoms() == 0) return 0.0;
    const Eigen::MatrixXd& normAdj = getNormalizedAdjacency(rdkMol);
    const Eigen::VectorXd& eigenvalues = getEigenvalues(normAdj);
    if (eigenvalues.size() == 0) return 0.0;
    return eigenvalues.array().exp().sum();
}

GraphBipartivityIndex::GraphBipartivityIndex() : EigenDescriptor("GraphBipart", "Graph bipartivity index (sum neg_eig^2 / sum eig^2)") {}
std::variant<double, int, std::string> GraphBipartivityIndex::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    if (rdkMol->getNumAtoms() == 0) return 1.0; // Empty graph is bipartite
    const Eigen::MatrixXd& adjacency = getAdjacencyMatrix(rdkMol);
    const Eigen::VectorXd& eigenvalues = getEigenvalues(adjacency);
    if (eigenvalues.size() == 0) return 1.0; // No eigenvalues? Assume bipartite?

    double sum_sq = eigenvalues.squaredNorm();
    if (sum_sq < 1e-10) return 1.0; // Graph with no edges is bipartite

    // Check symmetry: lambda_min approx -lambda_max
    if (std::abs(eigenvalues.minCoeff() + eigenvalues.maxCoeff()) < 1e-9) {
        return 1.0; // Approximately symmetric spectrum -> bipartite
    } else {
        // Calculate index = 1 - |sum(lambda_i * exp(lambda_i))| / sum(|lambda_i| * exp(|lambda_i|))
        // Alternative: sum_sq_neg / sum_sq (as implemented before)
        double sum_sq_neg = 0.0;
        for(int i=0; i<eigenvalues.size(); ++i) {
            if (eigenvalues(i) < -1e-10) {
                sum_sq_neg += eigenvalues(i) * eigenvalues(i);
            }
        }
       return sum_sq_neg / sum_sq;
    }
}

SpanningTreeEntropy::SpanningTreeEntropy() : EigenDescriptor("SpanTreeEntropy", "Spanning tree entropy (Requires Laplacian eigenvalues)") {}
std::variant<double, int, std::string> SpanningTreeEntropy::calculate([[maybe_unused]] const Molecule& mol) const { 
    return "Not Implemented"; 
}

// --- Implementations for newly added descriptor classes ---

// DistMat based
DistMatSpecRad::DistMatSpecRad() : EigenDescriptor("DistMatSpecRad", "Spectral radius of Distance matrix") {}
std::variant<double, int, std::string> DistMatSpecRad::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    if (rdkMol->getNumAtoms() == 0) return 0.0;
    const Eigen::MatrixXd& distMat = getDistanceMatrix(rdkMol);
    if ((distMat.array() == std::numeric_limits<double>::infinity()).any()) {
        return "Error: Molecule has disconnected components (DistMat)";
    }
    const Eigen::VectorXd& eigenvalues = getEigenvalues(distMat);
    if (eigenvalues.size() == 0) return 0.0;
    return eigenvalues.cwiseAbs().maxCoeff();
}

DistMatMinEig::DistMatMinEig() : EigenDescriptor("DistMatMinEig", "Smallest eigenvalue of Distance matrix") {}
std::variant<double, int, std::string> DistMatMinEig::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    if (rdkMol->getNumAtoms() == 0) return 0.0;
    const Eigen::MatrixXd& distMat = getDistanceMatrix(rdkMol);
    if ((distMat.array() == std::numeric_limits<double>::infinity()).any()) {
        return "Error: Molecule has disconnected components (DistMat)";
    }
    const Eigen::VectorXd& eigenvalues = getEigenvalues(distMat);
    if (eigenvalues.size() == 0) return 0.0;
    return eigenvalues.minCoeff();
}

DistMatEigSum::DistMatEigSum() : EigenDescriptor("DistMatEigSum", "Sum of Distance matrix eigenvalues (Trace)") {}
std::variant<double, int, std::string> DistMatEigSum::calculate(const Molecule& mol) const {
     if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    if (rdkMol->getNumAtoms() == 0) return 0.0;
    const Eigen::MatrixXd& distMat = getDistanceMatrix(rdkMol);
     if ((distMat.array() == std::numeric_limits<double>::infinity()).any()) {
         return "Error: Molecule has disconnected components (DistMat)";
     }
    return 0.0; // Trace of distance matrix is 0
}

DistMatEigVar::DistMatEigVar() : EigenDescriptor("DistMatEigVar", "Variance of Distance matrix eigenvalues") {}
std::variant<double, int, std::string> DistMatEigVar::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    if (rdkMol->getNumAtoms() == 0) return 0.0;
    const Eigen::MatrixXd& distMat = getDistanceMatrix(rdkMol);
     if ((distMat.array() == std::numeric_limits<double>::infinity()).any()) {
         return "Error: Molecule has disconnected components (DistMat)";
     }
    const Eigen::VectorXd& eigenvalues = getEigenvalues(distMat);
    if (eigenvalues.size() == 0) return 0.0;
    return computeVariance(eigenvalues);
}

DistMatEigSkew::DistMatEigSkew() : EigenDescriptor("DistMatEigSkew", "Skewness of Distance matrix eigenvalues") {}
std::variant<double, int, std::string> DistMatEigSkew::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    if (rdkMol->getNumAtoms() < 3) return 0.0;
    const Eigen::MatrixXd& distMat = getDistanceMatrix(rdkMol);
     if ((distMat.array() == std::numeric_limits<double>::infinity()).any()) {
         return "Error: Molecule has disconnected components (DistMat)";
     }
    const Eigen::VectorXd& eigenvalues = getEigenvalues(distMat);
    if (eigenvalues.size() < 3) return 0.0;
    return computeSkewness(eigenvalues);
}

DistMatEigKurt::DistMatEigKurt() : EigenDescriptor("DistMatEigKurt", "Kurtosis of Distance matrix eigenvalues") {}
std::variant<double, int, std::string> DistMatEigKurt::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    if (rdkMol->getNumAtoms() < 4) return 0.0;
    const Eigen::MatrixXd& distMat = getDistanceMatrix(rdkMol);
     if ((distMat.array() == std::numeric_limits<double>::infinity()).any()) {
         return "Error: Molecule has disconnected components (DistMat)";
     }
    const Eigen::VectorXd& eigenvalues = getEigenvalues(distMat);
    if (eigenvalues.size() < 4) return 0.0;
    return computeKurtosis(eigenvalues);
}

DistMatEnergy::DistMatEnergy() : EigenDescriptor("DistMatEnergy", "Energy of Distance matrix (sum |lambda_i|)") {}
std::variant<double, int, std::string> DistMatEnergy::calculate(const Molecule& mol) const {
     if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    if (rdkMol->getNumAtoms() == 0) return 0.0;
    const Eigen::MatrixXd& distMat = getDistanceMatrix(rdkMol);
     if ((distMat.array() == std::numeric_limits<double>::infinity()).any()) {
         return "Error: Molecule has disconnected components (DistMat)";
     }
    const Eigen::VectorXd& eigenvalues = getEigenvalues(distMat);
    if (eigenvalues.size() == 0) return 0.0;
    return eigenvalues.cwiseAbs().sum();
}

DistMatFrobNorm::DistMatFrobNorm() : EigenDescriptor("DistMatFrobNorm", "Frobenius norm of Distance matrix") {}
std::variant<double, int, std::string> DistMatFrobNorm::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    if (rdkMol->getNumAtoms() == 0) return 0.0;
    const Eigen::MatrixXd& distMat = getDistanceMatrix(rdkMol);
    if ((distMat.array() == std::numeric_limits<double>::infinity()).any()) {
         return "Error: Molecule has disconnected components (DistMat)";
     }
    return distMat.norm();
}

DistMatRank::DistMatRank() : EigenDescriptor("DistMatRank", "Rank of Distance matrix") {}
std::variant<double, int, std::string> DistMatRank::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    if (rdkMol->getNumAtoms() == 0) return 0;
    const Eigen::MatrixXd& distMat = getDistanceMatrix(rdkMol);
    if ((distMat.array() == std::numeric_limits<double>::infinity()).any()) {
         return "Error: Molecule has disconnected components (DistMat)";
     }
    Eigen::FullPivLU<Eigen::MatrixXd> lu(distMat);
    return static_cast<int>(lu.rank());
}

DistMatDet::DistMatDet() : EigenDescriptor("DistMatDet", "Determinant of Distance matrix") {}
std::variant<double, int, std::string> DistMatDet::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    if (rdkMol->getNumAtoms() == 0) return 1.0;
    const Eigen::MatrixXd& distMat = getDistanceMatrix(rdkMol);
    if ((distMat.array() == std::numeric_limits<double>::infinity()).any()) {
         return "Error: Molecule has disconnected components (DistMat)";
     }
    return distMat.determinant();
}


DistMatConditionNumber::DistMatConditionNumber() : EigenDescriptor("DistMatConditionNumber", "Condition number of Distance matrix") {}
std::variant<double, int, std::string> DistMatConditionNumber::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    if (rdkMol->getNumAtoms() == 0) return 0.0;
    const Eigen::MatrixXd& distMat = getDistanceMatrix(rdkMol);
    if ((distMat.array() == std::numeric_limits<double>::infinity()).any()) {
         return "Error: Molecule has disconnected components (DistMat)";
     }
    const Eigen::VectorXd& singularValues = getSingularValues(distMat);
    if (singularValues.size() == 0) return 0.0;
    double max_sv = singularValues(0);
    double min_sv = 0.0;
     for (int i = singularValues.size() - 1; i >= 0; --i) {
         if (singularValues(i) > 1e-10) {
             min_sv = singularValues(i);
             break;
         }
     }
    return (min_sv > 1e-10) ? (max_sv / min_sv) : std::numeric_limits<double>::infinity();
}

DistMatSingValEntropy::DistMatSingValEntropy() : EigenDescriptor("DistMatSingValEntropy", "Entropy of normalized singular values of Distance matrix") {}
std::variant<double, int, std::string> DistMatSingValEntropy::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    if (rdkMol->getNumAtoms() == 0) return 0.0;
    const Eigen::MatrixXd& distMat = getDistanceMatrix(rdkMol);
    if ((distMat.array() == std::numeric_limits<double>::infinity()).any()) {
         return "Error: Molecule has disconnected components (DistMat)";
     }
    const Eigen::VectorXd& singularValues = getSingularValues(distMat);
    if (singularValues.size() == 0) return 0.0;
    double sum_sv = singularValues.sum();
    if (sum_sv < 1e-10) return 0.0;

    double entropy = 0.0;
    for(int i=0; i<singularValues.size(); ++i) {
        if (singularValues(i) > 1e-10) {
            double p = singularValues(i) / sum_sv;
            entropy -= p * std::log2(p);
        }
    }
    return entropy;
}

DistMatSpectralEffRank::DistMatSpectralEffRank() : EigenDescriptor("DistMatSpectralEffRank", "Spectral effective rank of Distance matrix") {}
std::variant<double, int, std::string> DistMatSpectralEffRank::calculate(const Molecule& mol) const {
    auto entropy_val = DistMatSingValEntropy().calculate(mol);
    if (std::holds_alternative<double>(entropy_val)) {
        return std::pow(2.0, std::get<double>(entropy_val));
    } else {
        return entropy_val; // Propagate error
    }
}

// WtNumAdj based
WtNumAdjSpecRad::WtNumAdjSpecRad() : EigenDescriptor("WtNumAdjSpecRad", "Spectral radius of AtomicNumber Weighted Adjacency") {}
std::variant<double, int, std::string> WtNumAdjSpecRad::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    if (rdkMol->getNumAtoms() == 0) return 0.0;
    const Eigen::MatrixXd& wtAdj = getAtomicNumberWeightedAdjacency(rdkMol);
    const Eigen::VectorXd& eigenvalues = getEigenvalues(wtAdj);
    if (eigenvalues.size() == 0) return 0.0;
    return eigenvalues.cwiseAbs().maxCoeff();
}

WtNumAdjMinEig::WtNumAdjMinEig() : EigenDescriptor("WtNumAdjMinEig", "Smallest eigenvalue of AtomicNumber Weighted Adjacency") {}
std::variant<double, int, std::string> WtNumAdjMinEig::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    if (rdkMol->getNumAtoms() == 0) return 0.0;
    const Eigen::MatrixXd& wtAdj = getAtomicNumberWeightedAdjacency(rdkMol);
    const Eigen::VectorXd& eigenvalues = getEigenvalues(wtAdj);
    if (eigenvalues.size() == 0) return 0.0;
    return eigenvalues.minCoeff();
}

WtNumAdjEigSum::WtNumAdjEigSum() : EigenDescriptor("WtNumAdjEigSum", "Sum of eigenvalues of AtomicNumber Weighted Adjacency (Trace)") {}
std::variant<double, int, std::string> WtNumAdjEigSum::calculate(const Molecule& mol) const {
     if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    if (rdkMol->getNumAtoms() == 0) return 0.0;
    const Eigen::MatrixXd& wtAdj = getAtomicNumberWeightedAdjacency(rdkMol);
    return wtAdj.trace(); // Trace is 0 for simple graphs
}

WtNumAdjEigVar::WtNumAdjEigVar() : EigenDescriptor("WtNumAdjEigVar", "Variance of eigenvalues of AtomicNumber Weighted Adjacency") {}
std::variant<double, int, std::string> WtNumAdjEigVar::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    if (rdkMol->getNumAtoms() == 0) return 0.0;
    const Eigen::MatrixXd& wtAdj = getAtomicNumberWeightedAdjacency(rdkMol);
    const Eigen::VectorXd& eigenvalues = getEigenvalues(wtAdj);
    if (eigenvalues.size() == 0) return 0.0;
    return computeVariance(eigenvalues);
}

WtNumAdjEnergy::WtNumAdjEnergy() : EigenDescriptor("WtNumAdjEnergy", "Energy of AtomicNumber Weighted Adjacency (sum |lambda_i|)") {}
std::variant<double, int, std::string> WtNumAdjEnergy::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    if (rdkMol->getNumAtoms() == 0) return 0.0;
    const Eigen::MatrixXd& wtAdj = getAtomicNumberWeightedAdjacency(rdkMol);
    const Eigen::VectorXd& eigenvalues = getEigenvalues(wtAdj);
    if (eigenvalues.size() == 0) return 0.0;
    return eigenvalues.cwiseAbs().sum();
}

// WtNumLap based
WtNumLapSpecRad::WtNumLapSpecRad() : EigenDescriptor("WtNumLapSpecRad", "Spectral radius of AtomicNumber Weighted Laplacian") {}
std::variant<double, int, std::string> WtNumLapSpecRad::calculate(const Molecule& mol) const {
     if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    if (rdkMol->getNumAtoms() == 0) return 0.0;
    const Eigen::MatrixXd& wtLap = getAtomicNumberWeightedLaplacian(rdkMol);
    const Eigen::VectorXd& eigenvalues = getEigenvalues(wtLap);
    if (eigenvalues.size() == 0) return 0.0;
    return eigenvalues.maxCoeff();
}

WtNumLapAlgConn::WtNumLapAlgConn() : EigenDescriptor("WtNumLapAlgConn", "Algebraic connectivity of AtomicNumber Weighted Laplacian") {}
std::variant<double, int, std::string> WtNumLapAlgConn::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    if (rdkMol->getNumAtoms() <= 1) return 0.0;
    const Eigen::MatrixXd& wtLap = getAtomicNumberWeightedLaplacian(rdkMol);
    const Eigen::VectorXd& eigenvalues = getEigenvalues(wtLap);
    if (eigenvalues.size() < 2) return 0.0;
    return eigenvalues(1);
}

WtNumLapEigSum::WtNumLapEigSum() : EigenDescriptor("WtNumLapEigSum", "Sum of eigenvalues of AtomicNumber Weighted Laplacian (Trace)") {}
std::variant<double, int, std::string> WtNumLapEigSum::calculate(const Molecule& mol) const {
     if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    if (rdkMol->getNumAtoms() == 0) return 0.0;
    const Eigen::MatrixXd& wtLap = getAtomicNumberWeightedLaplacian(rdkMol);
    return wtLap.trace();
}

WtNumLapEigVar::WtNumLapEigVar() : EigenDescriptor("WtNumLapEigVar", "Variance of eigenvalues of AtomicNumber Weighted Laplacian") {}
std::variant<double, int, std::string> WtNumLapEigVar::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    if (rdkMol->getNumAtoms() == 0) return 0.0;
    const Eigen::MatrixXd& wtLap = getAtomicNumberWeightedLaplacian(rdkMol);
    const Eigen::VectorXd& eigenvalues = getEigenvalues(wtLap);
    if (eigenvalues.size() == 0) return 0.0;
    return computeVariance(eigenvalues);
}

WtNumLapEnergy::WtNumLapEnergy() : EigenDescriptor("WtNumLapEnergy", "Energy of AtomicNumber Weighted Laplacian (sum |lambda_i - mean_wt_deg|)") {}
std::variant<double, int, std::string> WtNumLapEnergy::calculate(const Molecule& mol) const {
     if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    int n = rdkMol->getNumAtoms();
    if (n == 0) return 0.0;
    const Eigen::MatrixXd& wtLap = getAtomicNumberWeightedLaplacian(rdkMol);
    const Eigen::VectorXd& eigenvalues = getEigenvalues(wtLap);
    if (eigenvalues.size() == 0) return 0.0;
    double mean_wt_deg = wtLap.trace() / n;
    return (eigenvalues.array() - mean_wt_deg).abs().sum();
}

// AdjPow based
AdjPow2SpecRad::AdjPow2SpecRad() : EigenDescriptor("AdjPow2SpecRad", "Spectral radius of A^2") {}
std::variant<double, int, std::string> AdjPow2SpecRad::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    if (rdkMol->getNumAtoms() == 0) return 0.0;
    const Eigen::MatrixXd& adj2 = getAdjPow2(rdkMol);
    const Eigen::VectorXd& eigenvalues = getEigenvalues(adj2);
    if (eigenvalues.size() == 0) return 0.0;
    return eigenvalues.maxCoeff();
}

AdjPow2Energy::AdjPow2Energy() : EigenDescriptor("AdjPow2Energy", "Energy of A^2 (sum |lambda_i|)") {}
std::variant<double, int, std::string> AdjPow2Energy::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    if (rdkMol->getNumAtoms() == 0) return 0.0;
    const Eigen::MatrixXd& adj2 = getAdjPow2(rdkMol);
    const Eigen::VectorXd& eigenvalues = getEigenvalues(adj2);
    if (eigenvalues.size() == 0) return 0.0;
    return eigenvalues.sum(); // Since eigenvalues >= 0
}

AdjPow3SpecRad::AdjPow3SpecRad() : EigenDescriptor("AdjPow3SpecRad", "Spectral radius of A^3") {}
std::variant<double, int, std::string> AdjPow3SpecRad::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    if (rdkMol->getNumAtoms() == 0) return 0.0;
    const Eigen::MatrixXd& adj3 = getAdjPow3(rdkMol);
    const Eigen::VectorXd& eigenvalues = getEigenvalues(adj3);
    if (eigenvalues.size() == 0) return 0.0;
    return eigenvalues.cwiseAbs().maxCoeff();
}

AdjPow3Energy::AdjPow3Energy() : EigenDescriptor("AdjPow3Energy", "Energy of A^3 (sum |lambda_i|)") {}
std::variant<double, int, std::string> AdjPow3Energy::calculate(const Molecule& mol) const {
     if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    if (rdkMol->getNumAtoms() == 0) return 0.0;
    const Eigen::MatrixXd& adj3 = getAdjPow3(rdkMol);
    const Eigen::VectorXd& eigenvalues = getEigenvalues(adj3);
    if (eigenvalues.size() == 0) return 0.0;
    return eigenvalues.cwiseAbs().sum();
}

AdjPow4SpecRad::AdjPow4SpecRad() : EigenDescriptor("AdjPow4SpecRad", "Spectral radius of A^4") {}
std::variant<double, int, std::string> AdjPow4SpecRad::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    if (rdkMol->getNumAtoms() == 0) return 0.0;
    const Eigen::MatrixXd& adj4 = getAdjPow4(rdkMol);
    const Eigen::VectorXd& eigenvalues = getEigenvalues(adj4);
    if (eigenvalues.size() == 0) return 0.0;
    return eigenvalues.maxCoeff();
}

AdjPow4Energy::AdjPow4Energy() : EigenDescriptor("AdjPow4Energy", "Energy of A^4 (sum |lambda_i|)") {}
std::variant<double, int, std::string> AdjPow4Energy::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    if (rdkMol->getNumAtoms() == 0) return 0.0;
    const Eigen::MatrixXd& adj4 = getAdjPow4(rdkMol);
    const Eigen::VectorXd& eigenvalues = getEigenvalues(adj4);
    if (eigenvalues.size() == 0) return 0.0;
    return eigenvalues.sum(); // Since eigenvalues >= 0
}

// LapPow based
LapPow2SpecRad::LapPow2SpecRad() : EigenDescriptor("LapPow2SpecRad", "Spectral radius of L^2") {}
std::variant<double, int, std::string> LapPow2SpecRad::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    if (rdkMol->getNumAtoms() == 0) return 0.0;
    const Eigen::MatrixXd& lap2 = getLapPow2(rdkMol);
    const Eigen::VectorXd& eigenvalues = getEigenvalues(lap2);
    if (eigenvalues.size() == 0) return 0.0;
    return eigenvalues.maxCoeff();
}

LapPow2AlgConn::LapPow2AlgConn() : EigenDescriptor("LapPow2AlgConn", "Second smallest eigenvalue of L^2") {}
std::variant<double, int, std::string> LapPow2AlgConn::calculate(const Molecule& mol) const {
     if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    if (rdkMol->getNumAtoms() <= 1) return 0.0;
    const Eigen::MatrixXd& lap2 = getLapPow2(rdkMol);
    const Eigen::VectorXd& eigenvalues = getEigenvalues(lap2);
    if (eigenvalues.size() < 2) return 0.0;
    return eigenvalues(1);
}

LapPow2Energy::LapPow2Energy() : EigenDescriptor("LapPow2Energy", "Energy of L^2 (sum |lambda_i|)") {}
std::variant<double, int, std::string> LapPow2Energy::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    if (rdkMol->getNumAtoms() == 0) return 0.0;
    const Eigen::MatrixXd& lap2 = getLapPow2(rdkMol);
    const Eigen::VectorXd& eigenvalues = getEigenvalues(lap2);
    if (eigenvalues.size() == 0) return 0.0;
    return eigenvalues.sum(); // Since eigenvalues >= 0
}

LapPow2Trace::LapPow2Trace() : EigenDescriptor("LapPow2Trace", "Trace of L^2") {}
std::variant<double, int, std::string> LapPow2Trace::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    if (rdkMol->getNumAtoms() == 0) return 0.0;
    const Eigen::MatrixXd& lap2 = getLapPow2(rdkMol);
    return lap2.trace();
}

// Matrix Exponential based
DistMatEstradaIdx::DistMatEstradaIdx() : EigenDescriptor("DistMatEstradaIdx", "Estrada index of Distance matrix") {}
std::variant<double, int, std::string> DistMatEstradaIdx::calculate(const Molecule& mol) const {
     if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    if (rdkMol->getNumAtoms() == 0) return 0.0;
    const Eigen::MatrixXd& distMat = getDistanceMatrix(rdkMol);
    if ((distMat.array() == std::numeric_limits<double>::infinity()).any()) {
         return "Error: Molecule has disconnected components (DistMat)";
     }
    const Eigen::VectorXd& eigenvalues = getEigenvalues(distMat);
    if (eigenvalues.size() == 0) return 0.0;
    return eigenvalues.array().exp().sum();
}

WtNumAdjEstradaIdx::WtNumAdjEstradaIdx() : EigenDescriptor("WtNumAdjEstradaIdx", "Estrada index of AtomicNumber Weighted Adjacency") {}
std::variant<double, int, std::string> WtNumAdjEstradaIdx::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    if (rdkMol->getNumAtoms() == 0) return 0.0;
    const Eigen::MatrixXd& wtAdj = getAtomicNumberWeightedAdjacency(rdkMol);
    const Eigen::VectorXd& eigenvalues = getEigenvalues(wtAdj);
    if (eigenvalues.size() == 0) return 0.0;
    return eigenvalues.array().exp().sum();
}

WtNumLapEstradaIdx::WtNumLapEstradaIdx() : EigenDescriptor("WtNumLapEstradaIdx", "Estrada index of AtomicNumber Weighted Laplacian") {}
std::variant<double, int, std::string> WtNumLapEstradaIdx::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    if (rdkMol->getNumAtoms() == 0) return 0.0;
    const Eigen::MatrixXd& wtLap = getAtomicNumberWeightedLaplacian(rdkMol);
    const Eigen::VectorXd& eigenvalues = getEigenvalues(wtLap);
    if (eigenvalues.size() == 0) return 0.0;
    return eigenvalues.array().exp().sum();
}


// Eigenvector Statistics
AdjDomEigVecVar::AdjDomEigVecVar() : EigenDescriptor("AdjDomEigVecVar", "Variance of dominant adjacency eigenvector components") {}
std::variant<double, int, std::string> AdjDomEigVecVar::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    if (rdkMol->getNumAtoms() == 0) return 0.0;
    const Eigen::MatrixXd& adjacency = getAdjacencyMatrix(rdkMol);
    const Eigen::VectorXd& eigenvector = getDominantEigenvector(adjacency);
    if (eigenvector.size() == 0) return 0.0;
    return computeVariance(eigenvector);
}

AdjDomEigVecSkew::AdjDomEigVecSkew() : EigenDescriptor("AdjDomEigVecSkew", "Skewness of dominant adjacency eigenvector components") {}
std::variant<double, int, std::string> AdjDomEigVecSkew::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
     if (rdkMol->getNumAtoms() < 3) return 0.0;
    const Eigen::MatrixXd& adjacency = getAdjacencyMatrix(rdkMol);
    const Eigen::VectorXd& eigenvector = getDominantEigenvector(adjacency);
    if (eigenvector.size() < 3) return 0.0;
    return computeSkewness(eigenvector);
}

AdjDomEigVecKurt::AdjDomEigVecKurt() : EigenDescriptor("AdjDomEigVecKurt", "Kurtosis of dominant adjacency eigenvector components") {}
std::variant<double, int, std::string> AdjDomEigVecKurt::calculate(const Molecule& mol) const {
     if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
     if (rdkMol->getNumAtoms() < 4) return 0.0;
    const Eigen::MatrixXd& adjacency = getAdjacencyMatrix(rdkMol);
    const Eigen::VectorXd& eigenvector = getDominantEigenvector(adjacency);
    if (eigenvector.size() < 4) return 0.0;
    return computeKurtosis(eigenvector);
}

NormLapDomEigVecVar::NormLapDomEigVecVar() : EigenDescriptor("NormLapDomEigVecVar", "Variance of dominant normLaplacian eigenvector components") {}
std::variant<double, int, std::string> NormLapDomEigVecVar::calculate(const Molecule& mol) const {
     if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    if (rdkMol->getNumAtoms() == 0) return 0.0;
    const Eigen::MatrixXd& normLap = getNormalizedLaplacianMatrix(rdkMol);
    const Eigen::VectorXd& eigenvector = getDominantEigenvector(normLap);
    if (eigenvector.size() == 0) return 0.0;
    return computeVariance(eigenvector);
}

SignLapDomEigVecVar::SignLapDomEigVecVar() : EigenDescriptor("SignLapDomEigVecVar", "Variance of dominant signlessLaplacian eigenvector components") {}
std::variant<double, int, std::string> SignLapDomEigVecVar::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    if (rdkMol->getNumAtoms() == 0) return 0.0;
    const Eigen::MatrixXd& signLap = getSignlessLaplacianMatrix(rdkMol);
    const Eigen::VectorXd& eigenvector = getDominantEigenvector(signLap);
    if (eigenvector.size() == 0) return 0.0;
    return computeVariance(eigenvector);
}


// Miscellaneous Eigen-based
SumFormanRicci::SumFormanRicci() : EigenDescriptor("SumFormanRicci", "Sum of Forman-Ricci curvatures (Approximate)") {}
std::variant<double, int, std::string> SumFormanRicci::calculate([[maybe_unused]] const Molecule& mol) const { 
    return "Not Implemented"; 
}

LogDetLaplacian::LogDetLaplacian() : EigenDescriptor("LogDetLaplacian", "Log of pseudo-determinant of Laplacian") {}
std::variant<double, int, std::string> LogDetLaplacian::calculate(const Molecule& mol) const {
    // log(pseudo-det) = sum(log(lambda_i)) for non-zero lambda_i
     if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    int n = rdkMol->getNumAtoms();
    if (n <= 1) return 0.0; // log(pseudo-det) is 0 for N=1 (1 spanning tree), undefined for N=0

    auto zeroEigResult = LaplacianZeroEigenvalues().calculate(mol);
    if (std::holds_alternative<int>(zeroEigResult)) {
        if (std::get<int>(zeroEigResult) != 1) {
             return "Error: Disconnected graph (LogDetLaplacian)";
        }
    } else {
         return zeroEigResult; // Propagate error
    }

    const Eigen::MatrixXd& laplacian = getLaplacianMatrix(rdkMol);
    const Eigen::VectorXd& eigenvalues = getEigenvalues(laplacian); // Assumes sorted
    if (eigenvalues.size() < n) return "Error: Eigenvalue calculation failed (LogDetLaplacian)"; // Safety

    double log_sum = 0.0;
    int non_zero_count = 0;
    for (int i = 1; i < eigenvalues.size(); ++i) { // Skip lambda_0 = 0
        if (eigenvalues(i) > 1e-10) {
            log_sum += std::log(eigenvalues(i));
            non_zero_count++;
        } else {
             // This shouldn't happen for a connected graph
             return "Error: Unexpected non-positive eigenvalue > 0 (LogDetLaplacian)";
        }
    }
    // For a connected graph, there should be exactly n-1 non-zero eigenvalues
    if (non_zero_count != n - 1) {
       // This indicates a potential issue, possibly numerical instability or unexpected graph structure
       // For robustness, we proceed but could also return an error/warning here.
    }
    return log_sum;
}

AvgDegreeDistance::AvgDegreeDistance() : EigenDescriptor("AvgDegreeDistance", "Average degree distance (Sum(deg_i*deg_j*dist_ij))") {}
std::variant<double, int, std::string> AvgDegreeDistance::calculate(const Molecule& mol) const {
     if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    int n = rdkMol->getNumAtoms();
    if (n <= 1) return 0.0;

    const Eigen::MatrixXd& degreeMat = getDegreeMatrix(rdkMol);
    if (degreeMat.rows() != n) return "Error: Degree matrix size mismatch (AvgDegDist)"; // Safety
    Eigen::VectorXd degrees = degreeMat.diagonal();
    if (degrees.size() != n) return "Error: Degree vector size mismatch (AvgDegDist)"; // Safety

    const double* distMatRaw = RDKit::MolOps::getDistanceMat(*rdkMol, false, false, true);

    double sum_deg_dist = 0.0;
    int count = 0;
    bool connected = true;

    for(int i=0; i<n; ++i) {
        for(int j=i+1; j<n; ++j) {
             if (i * n + j >= n * n) return "Error: Distance matrix index out of bounds (AvgDegDist)"; // Safety
            double d = distMatRaw[i*n + j];
            if (d < 0) {
                connected = false;
                break;
            }
            sum_deg_dist += degrees(i) * degrees(j) * d;
            count++;
        }
        if (!connected) break;
    }

     if (!connected) {
         return "Error: Molecule has disconnected components (AvgDegDist)";
     }
     return (count > 0) ? (sum_deg_dist / count) : 0.0;
}

VarianceDegreeDistance::VarianceDegreeDistance() : EigenDescriptor("VarianceDegreeDistance", "Variance of degree distance (deg_i*deg_j*dist_ij)") {}
std::variant<double, int, std::string> VarianceDegreeDistance::calculate(const Molecule& mol) const {
    if (!mol.isValid()) return "Error: Invalid molecule";
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit molecule is null";
    int n = rdkMol->getNumAtoms();
    if (n <= 1) return 0.0;

    const Eigen::MatrixXd& degreeMat = getDegreeMatrix(rdkMol);
     if (degreeMat.rows() != n) return "Error: Degree matrix size mismatch (VarDegDist)"; // Safety
    Eigen::VectorXd degrees = degreeMat.diagonal();
    if (degrees.size() != n) return "Error: Degree vector size mismatch (VarDegDist)"; // Safety

    const double* distMatRaw = RDKit::MolOps::getDistanceMat(*rdkMol, false, false, true);

    std::vector<double> deg_dist_vec;
    bool connected = true;

    for(int i=0; i<n; ++i) {
        for(int j=i+1; j<n; ++j) {
            if (i * n + j >= n * n) return "Error: Distance matrix index out of bounds (VarDegDist)"; // Safety
            double d = distMatRaw[i*n + j];
            if (d < 0) {
                connected = false;
                break;
            }
            deg_dist_vec.push_back(degrees(i) * degrees(j) * d);
        }
        if (!connected) break;
    }

     if (!connected) {
         return "Error: Molecule has disconnected components (VarDegDist)";
     }

     if(deg_dist_vec.empty()) return 0.0;

     Eigen::Map<Eigen::VectorXd> deg_dists_eigen(deg_dist_vec.data(), deg_dist_vec.size());

     return computeVariance(deg_dists_eigen);
}

} // namespace descriptors
} // namespace desfact