#pragma once

#include "descriptors.hpp"
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <GraphMol/Atom.h>
#include <GraphMol/ROMol.h>
#include <variant>
#include <string>
#include <cmath>

namespace desfact {
namespace descriptors {

// Base class for Eigen-based descriptors
class EigenDescriptor : public Descriptor {
protected:
    // Cache management constants
    static constexpr size_t MAX_CACHE_ENTRIES = 1000;  // Maximum number of entries in each cache
    static constexpr size_t CACHE_CLEANUP_THRESHOLD = 800;  // When to trigger cleanup
    
    // Thread-local caches for matrices and eigenvalues
    struct CacheEntry {
        std::chrono::steady_clock::time_point lastAccess;
        size_t approximateSize;  // Rough estimate of memory usage in bytes
    };
    
    // Matrix caches with pointer as key (molecule-specific)
    struct MatrixCacheEntry : CacheEntry {
        Eigen::MatrixXd matrix;
    };
    
    // Eigenvalue cache with matrix hash as key
    struct EigenvalueCacheEntry : CacheEntry {
        Eigen::VectorXd eigenvalues;
    };
    
    // Thread-local caches
    static thread_local std::unordered_map<const RDKit::ROMol*, MatrixCacheEntry> adjacencyCache;
    static thread_local std::unordered_map<const RDKit::ROMol*, MatrixCacheEntry> degreeCache;
    static thread_local std::unordered_map<const RDKit::ROMol*, MatrixCacheEntry> laplacianCache;
    static thread_local std::unordered_map<const RDKit::ROMol*, MatrixCacheEntry> normalizedLaplacianCache;
    static thread_local std::unordered_map<const RDKit::ROMol*, MatrixCacheEntry> signlessLaplacianCache;
    static thread_local std::unordered_map<const RDKit::ROMol*, MatrixCacheEntry> weightedAdjacencyCache;
    static thread_local std::unordered_map<const RDKit::ROMol*, MatrixCacheEntry> distanceMatrixCache;
    
    // Hash function for matrices (for eigenvalue cache)
    struct MatrixHash {
        std::size_t operator()(const Eigen::MatrixXd& matrix) const {
            // Simple hash based on matrix dimensions and a sample of values
            std::size_t seed = matrix.rows() * 73 + matrix.cols();
            // Sample a few values from the matrix for the hash
            for (int i = 0; i < std::min(3, (int)matrix.rows()); ++i) {
                for (int j = 0; j < std::min(3, (int)matrix.cols()); ++j) {
                    seed ^= std::hash<double>{}(matrix(i,j)) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
                }
            }
            return seed;
        }
    };
    
    struct MatrixEqual {
        bool operator()(const Eigen::MatrixXd& a, const Eigen::MatrixXd& b) const {
            if (a.rows() != b.rows() || a.cols() != b.cols()) return false;
            // Check a sample of values for equality
            for (int i = 0; i < std::min(5, (int)a.rows()); ++i) {
                for (int j = 0; j < std::min(5, (int)a.cols()); ++j) {
                    if (std::abs(a(i,j) - b(i,j)) > 1e-10) return false;
                }
            }
            return true;
        }
    };
    
    // Eigenvalue cache with matrix hash
    static thread_local std::unordered_map<Eigen::MatrixXd, EigenvalueCacheEntry, MatrixHash, MatrixEqual> eigenvalueCache;
    
    // Cache management functions
    template<typename CacheType>
    static void cleanupCache(CacheType& cache, size_t keepEntries) {
        if (cache.size() <= keepEntries) return;
        
        // Sort entries by last access time
        std::vector<std::pair<typename CacheType::key_type, std::chrono::steady_clock::time_point>> entries;
        entries.reserve(cache.size());
        
        for (const auto& entry : cache) {
            entries.emplace_back(entry.first, entry.second.lastAccess);
        }
        
        // Sort by access time (oldest first)
        std::sort(entries.begin(), entries.end(), 
            [](const auto& a, const auto& b) { return a.second < b.second; });
        
        // Remove oldest entries
        size_t removeCount = cache.size() - keepEntries;
        for (size_t i = 0; i < removeCount; ++i) {
            cache.erase(entries[i].first);
        }
    }
    
    // Check all caches and clean up if needed
    static void checkAndCleanupCaches() {
        if (adjacencyCache.size() > CACHE_CLEANUP_THRESHOLD) {
            cleanupCache(adjacencyCache, MAX_CACHE_ENTRIES / 2);
        }
        if (degreeCache.size() > CACHE_CLEANUP_THRESHOLD) {
            cleanupCache(degreeCache, MAX_CACHE_ENTRIES / 2);
        }
        if (laplacianCache.size() > CACHE_CLEANUP_THRESHOLD) {
            cleanupCache(laplacianCache, MAX_CACHE_ENTRIES / 2);
        }
        if (normalizedLaplacianCache.size() > CACHE_CLEANUP_THRESHOLD) {
            cleanupCache(normalizedLaplacianCache, MAX_CACHE_ENTRIES / 2);
        }
        if (signlessLaplacianCache.size() > CACHE_CLEANUP_THRESHOLD) {
            cleanupCache(signlessLaplacianCache, MAX_CACHE_ENTRIES / 2);
        }
        if (weightedAdjacencyCache.size() > CACHE_CLEANUP_THRESHOLD) {
            cleanupCache(weightedAdjacencyCache, MAX_CACHE_ENTRIES / 2);
        }
        if (distanceMatrixCache.size() > CACHE_CLEANUP_THRESHOLD) {
            cleanupCache(distanceMatrixCache, MAX_CACHE_ENTRIES / 2);
        }
        if (eigenvalueCache.size() > CACHE_CLEANUP_THRESHOLD) {
            cleanupCache(eigenvalueCache, MAX_CACHE_ENTRIES / 2);
        }
    }
    
    // Get cached or build adjacency matrix
    Eigen::MatrixXd getAdjacencyMatrix(const RDKit::ROMol* mol) const {
        auto it = adjacencyCache.find(mol);
        if (it != adjacencyCache.end()) {
            // Update last access time
            it->second.lastAccess = std::chrono::steady_clock::now();
            return it->second.matrix;
        }
        
        // Build the matrix
        Eigen::MatrixXd adj = buildAdjacencyMatrix(mol);
        
        // Check and clean caches if needed
        checkAndCleanupCaches();
        
        // Cache the result
        size_t approxSize = adj.rows() * adj.cols() * sizeof(double);
        adjacencyCache[mol] = {std::chrono::steady_clock::now(), approxSize, adj};
        
        return adj;
    }
    
    // Get cached or build degree matrix
    Eigen::MatrixXd getDegreeMatrix(const RDKit::ROMol* mol) const {
        auto it = degreeCache.find(mol);
        if (it != degreeCache.end()) {
            it->second.lastAccess = std::chrono::steady_clock::now();
            return it->second.matrix;
        }
        
        // Get adjacency matrix (possibly from cache)
        Eigen::MatrixXd adj = getAdjacencyMatrix(mol);
        Eigen::MatrixXd degree = buildDegreeMatrix(adj);
        
        // Check and clean caches if needed
        checkAndCleanupCaches();
        
        // Cache the result
        size_t approxSize = degree.rows() * degree.cols() * sizeof(double);
        degreeCache[mol] = {std::chrono::steady_clock::now(), approxSize, degree};
        
        return degree;
    }
    
    // Get cached or build Laplacian matrix
    Eigen::MatrixXd getLaplacianMatrix(const RDKit::ROMol* mol) const {
        auto it = laplacianCache.find(mol);
        if (it != laplacianCache.end()) {
            it->second.lastAccess = std::chrono::steady_clock::now();
            return it->second.matrix;
        }
        
        // Get adjacency and degree matrices (possibly from cache)
        Eigen::MatrixXd adj = getAdjacencyMatrix(mol);
        Eigen::MatrixXd degree = getDegreeMatrix(mol);
        Eigen::MatrixXd laplacian = buildLaplacianMatrix(degree, adj);
        
        // Check and clean caches if needed
        checkAndCleanupCaches();
        
        // Cache the result
        size_t approxSize = laplacian.rows() * laplacian.cols() * sizeof(double);
        laplacianCache[mol] = {std::chrono::steady_clock::now(), approxSize, laplacian};
        
        return laplacian;
    }
    
    // Get cached or build normalized Laplacian matrix
    Eigen::MatrixXd getNormalizedLaplacianMatrix(const RDKit::ROMol* mol) const {
        auto it = normalizedLaplacianCache.find(mol);
        if (it != normalizedLaplacianCache.end()) {
            it->second.lastAccess = std::chrono::steady_clock::now();
            return it->second.matrix;
        }
        
        // Get adjacency and degree matrices (possibly from cache)
        Eigen::MatrixXd adj = getAdjacencyMatrix(mol);
        Eigen::MatrixXd degree = getDegreeMatrix(mol);
        Eigen::MatrixXd normalizedLaplacian = buildNormalizedLaplacian(degree, adj);
        
        // Check and clean caches if needed
        checkAndCleanupCaches();
        
        // Cache the result
        size_t approxSize = normalizedLaplacian.rows() * normalizedLaplacian.cols() * sizeof(double);
        normalizedLaplacianCache[mol] = {std::chrono::steady_clock::now(), approxSize, normalizedLaplacian};
        
        return normalizedLaplacian;
    }
    
    // Get cached or build signless Laplacian matrix
    Eigen::MatrixXd getSignlessLaplacianMatrix(const RDKit::ROMol* mol) const {
        auto it = signlessLaplacianCache.find(mol);
        if (it != signlessLaplacianCache.end()) {
            it->second.lastAccess = std::chrono::steady_clock::now();
            return it->second.matrix;
        }
        
        // Get adjacency and degree matrices (possibly from cache)
        Eigen::MatrixXd adj = getAdjacencyMatrix(mol);
        Eigen::MatrixXd degree = getDegreeMatrix(mol);
        Eigen::MatrixXd signlessLaplacian = buildSignlessLaplacian(degree, adj);
        
        // Check and clean caches if needed
        checkAndCleanupCaches();
        
        // Cache the result
        size_t approxSize = signlessLaplacian.rows() * signlessLaplacian.cols() * sizeof(double);
        signlessLaplacianCache[mol] = {std::chrono::steady_clock::now(), approxSize, signlessLaplacian};
        
        return signlessLaplacian;
    }
    
    // Get cached or build weighted adjacency matrix
    Eigen::MatrixXd getWeightedAdjacencyMatrix(const RDKit::ROMol* mol) const {
        auto it = weightedAdjacencyCache.find(mol);
        if (it != weightedAdjacencyCache.end()) {
            it->second.lastAccess = std::chrono::steady_clock::now();
            return it->second.matrix;
        }
        
        Eigen::MatrixXd weightedAdj = buildWeightedAdjacency(mol);
        
        // Check and clean caches if needed
        checkAndCleanupCaches();
        
        // Cache the result
        size_t approxSize = weightedAdj.rows() * weightedAdj.cols() * sizeof(double);
        weightedAdjacencyCache[mol] = {std::chrono::steady_clock::now(), approxSize, weightedAdj};
        
        return weightedAdj;
    }
    
    // Get cached or build distance matrix
    Eigen::MatrixXd getDistanceMatrix(const RDKit::ROMol* mol) const {
        auto it = distanceMatrixCache.find(mol);
        if (it != distanceMatrixCache.end()) {
            it->second.lastAccess = std::chrono::steady_clock::now();
            return it->second.matrix;
        }
        
        Eigen::MatrixXd distMatrix = buildDistanceMatrix(mol);
        
        // Check and clean caches if needed
        checkAndCleanupCaches();
        
        // Cache the result
        size_t approxSize = distMatrix.rows() * distMatrix.cols() * sizeof(double);
        distanceMatrixCache[mol] = {std::chrono::steady_clock::now(), approxSize, distMatrix};
        
        return distMatrix;
    }
    
    // Get cached or compute eigenvalues
    Eigen::VectorXd getEigenvalues(const Eigen::MatrixXd& matrix) const {
        auto it = eigenvalueCache.find(matrix);
        if (it != eigenvalueCache.end()) {
            it->second.lastAccess = std::chrono::steady_clock::now();
            return it->second.eigenvalues;
        }
        
        Eigen::VectorXd evals = computeEigenvalues(matrix);
        
        // Check and clean caches if needed
        checkAndCleanupCaches();
        
        // Cache the result
        size_t approxSize = evals.size() * sizeof(double);
        eigenvalueCache[matrix] = {std::chrono::steady_clock::now(), approxSize, evals};
        
        return evals;
    }
    
    // Original matrix building methods (now used by the caching layer)
    Eigen::MatrixXd buildAdjacencyMatrix(const RDKit::ROMol* mol) const;
    Eigen::MatrixXd buildDegreeMatrix(const Eigen::MatrixXd& adjacency) const;
    Eigen::MatrixXd buildLaplacianMatrix(const Eigen::MatrixXd& degree, const Eigen::MatrixXd& adjacency) const;
    Eigen::MatrixXd buildNormalizedLaplacian(const Eigen::MatrixXd& degree, const Eigen::MatrixXd& adjacency) const;
    Eigen::MatrixXd buildSignlessLaplacian(const Eigen::MatrixXd& degree, const Eigen::MatrixXd& adjacency) const;
    Eigen::MatrixXd buildWeightedAdjacency(const RDKit::ROMol* mol) const;
    Eigen::MatrixXd buildDistanceMatrix(const RDKit::ROMol* mol) const;
    Eigen::MatrixXd buildAtomicNumberWeightedAdjacency(const RDKit::ROMol* mol) const;
    Eigen::MatrixXd buildAtomicNumberWeightedLaplacian(const RDKit::ROMol* mol) const;
    Eigen::MatrixXd buildNormalizedAdjacency(const Eigen::MatrixXd& degree, const Eigen::MatrixXd& adjacency) const;
    
    // Compute eigenvalues (now used by the caching layer)
    Eigen::VectorXd computeEigenvalues(const Eigen::MatrixXd& matrix) const;
    Eigen::VectorXd computeSingularValues(const Eigen::MatrixXd& matrix) const;
    
    // Statistical helpers
    double computeVariance(const Eigen::VectorXd& values) const;
    double computeSkewness(const Eigen::VectorXd& values) const;
    double computeKurtosis(const Eigen::VectorXd& values) const;
    double computeEntropy(const Eigen::VectorXd& values) const;
    
    // Get dominant eigenvector
    Eigen::VectorXd getDominantEigenvector(const Eigen::MatrixXd& matrix) const;
    
public:
    EigenDescriptor(const std::string& name, const std::string& description);
    
    // Static method to clear all caches (call on program exit or SIGINT)
    static void clearAllCaches();
};

// Adjacency matrix descriptors
class AdjacencyNonZeroEntries : public EigenDescriptor {
public:
    AdjacencyNonZeroEntries();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class AdjacencyMatrixTrace : public EigenDescriptor {
public:
    AdjacencyMatrixTrace();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class AdjacencyFrobeniusNorm : public EigenDescriptor {
public:
    AdjacencyFrobeniusNorm();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class AdjacencySpectralRadius : public EigenDescriptor {
public:
    AdjacencySpectralRadius();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class AdjacencySmallestEigenvalue : public EigenDescriptor {
public:
    AdjacencySmallestEigenvalue();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class AdjacencySumEigenvalues : public EigenDescriptor {
public:
    AdjacencySumEigenvalues();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class AdjacencySumSquaresEigenvalues : public EigenDescriptor {
public:
    AdjacencySumSquaresEigenvalues();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class AdjacencyEigenvalueVariance : public EigenDescriptor {
public:
    AdjacencyEigenvalueVariance();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class AdjacencyEigenvalueSkewness : public EigenDescriptor {
public:
    AdjacencyEigenvalueSkewness();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class AdjacencyEigenvalueKurtosis : public EigenDescriptor {
public:
    AdjacencyEigenvalueKurtosis();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class AdjacencyPositiveEigenvalues : public EigenDescriptor {
public:
    AdjacencyPositiveEigenvalues();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class AdjacencyNegativeEigenvalues : public EigenDescriptor {
public:
    AdjacencyNegativeEigenvalues();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class AdjacencyZeroEigenvalues : public EigenDescriptor {
public:
    AdjacencyZeroEigenvalues();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class AdjacencyMaxDegree : public EigenDescriptor {
public:
    AdjacencyMaxDegree();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class AdjacencyMinDegree : public EigenDescriptor {
public:
    AdjacencyMinDegree();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class AdjacencyMeanDegree : public EigenDescriptor {
public:
    AdjacencyMeanDegree();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class AdjacencyDegreeVariance : public EigenDescriptor {
public:
    AdjacencyDegreeVariance();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class AdjacencyDegreeStdDev : public EigenDescriptor {
public:
    AdjacencyDegreeStdDev();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class AdjacencyMatrixRank : public EigenDescriptor {
public:
    AdjacencyMatrixRank();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class GraphEnergy : public EigenDescriptor {
public:
    GraphEnergy();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

// Laplacian matrix descriptors
class LaplacianSpectralRadius : public EigenDescriptor {
public:
    LaplacianSpectralRadius();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class LaplacianAlgebraicConnectivity : public EigenDescriptor {
public:
    LaplacianAlgebraicConnectivity();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class LaplacianZeroEigenvalues : public EigenDescriptor {
public:
    LaplacianZeroEigenvalues();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class LaplacianEnergy : public EigenDescriptor {
public:
    LaplacianEnergy();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class LaplacianMatrixTrace : public EigenDescriptor {
public:
    LaplacianMatrixTrace();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class LaplacianMatrixDeterminant : public EigenDescriptor {
public:
    LaplacianMatrixDeterminant();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class LaplacianTotalEffectiveResistance : public EigenDescriptor {
public:
    LaplacianTotalEffectiveResistance();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class LaplacianKirchhoffIndex : public EigenDescriptor {
public:
    LaplacianKirchhoffIndex();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class LaplacianEigenvalueVariance : public EigenDescriptor {
public:
    LaplacianEigenvalueVariance();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class LaplacianEigenvalueSkewness : public EigenDescriptor {
public:
    LaplacianEigenvalueSkewness();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

// Normalized Laplacian descriptors
class NormalizedLaplacianSpectralRadius : public EigenDescriptor {
public:
    NormalizedLaplacianSpectralRadius();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class NormalizedLaplacianSmallestNonzero : public EigenDescriptor {
public:
    NormalizedLaplacianSmallestNonzero();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class NormalizedLaplacianLargestEigenvalue : public EigenDescriptor {
public:
    NormalizedLaplacianLargestEigenvalue();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class NormalizedLaplacianEnergy : public EigenDescriptor {
public:
    NormalizedLaplacianEnergy();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class NormalizedLaplacianTrace : public EigenDescriptor {
public:
    NormalizedLaplacianTrace();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

// Degree matrix descriptors
class DegreeMatrixMaxDegree : public EigenDescriptor {
public:
    DegreeMatrixMaxDegree();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class DegreeMatrixMinDegree : public EigenDescriptor {
public:
    DegreeMatrixMinDegree();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class DegreeMatrixAvgDegree : public EigenDescriptor {
public:
    DegreeMatrixAvgDegree();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class DegreeMatrixVariance : public EigenDescriptor {
public:
    DegreeMatrixVariance();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class DegreeMatrixSkewness : public EigenDescriptor {
public:
    DegreeMatrixSkewness();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class DegreeMatrixEntropy : public EigenDescriptor {
public:
    DegreeMatrixEntropy();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

// Other linear algebra descriptors
class GraphIrregularity : public EigenDescriptor {
public:
    GraphIrregularity();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class WienerIndex : public EigenDescriptor {
public:
    WienerIndex();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

// class EstradaIndex : public EigenDescriptor {
// public:
//     EstradaIndex();
//     std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
// }; // Commented out due to zero or near-zero variance

class EstradaIndex : public EigenDescriptor {
public:
    EstradaIndex();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class NumberSpanningTrees : public EigenDescriptor {
public:
    NumberSpanningTrees();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class GraphEccentricity : public EigenDescriptor {
public:
    GraphEccentricity();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class SpectralGap : public EigenDescriptor {
public:
    SpectralGap();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class TraceMatrixPower2 : public EigenDescriptor {
public:
    TraceMatrixPower2();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class NumberTriangles : public EigenDescriptor {
public:
    NumberTriangles();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class GraphDiameter : public EigenDescriptor {
public:
    GraphDiameter();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

// Higher-order adjacency powers descriptors
class NumberOf2Walks : public EigenDescriptor {
public:
    NumberOf2Walks();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class NumberOf3Walks : public EigenDescriptor {
public:
    NumberOf3Walks();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class NumberOf4Walks : public EigenDescriptor {
public:
    NumberOf4Walks();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class MeanClosed3WalksPerNode : public EigenDescriptor {
public:
    MeanClosed3WalksPerNode();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class MeanClosed4WalksPerNode : public EigenDescriptor {
public:
    MeanClosed4WalksPerNode();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class Walk2Energy : public EigenDescriptor {
public:
    Walk2Energy();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class Walk3Energy : public EigenDescriptor {
public:
    Walk3Energy();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class Walk4Energy : public EigenDescriptor {
public:
    Walk4Energy();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class GraphIrregularityWalkCount : public EigenDescriptor {
public:
    GraphIrregularityWalkCount();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class TrianglesToPathsRatio : public EigenDescriptor {
public:
    TrianglesToPathsRatio();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

// Singular value decomposition (SVD) of A descriptors
class MaxSingularValue : public EigenDescriptor {
public:
    MaxSingularValue();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class MinNonZeroSingularValue : public EigenDescriptor {
public:
    MinNonZeroSingularValue();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class ConditionNumber : public EigenDescriptor {
public:
    ConditionNumber();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class SumSingularValues : public EigenDescriptor {
public:
    SumSingularValues();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class FrobeniusNormSingularValues : public EigenDescriptor {
public:
    FrobeniusNormSingularValues();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class SingularValueEntropy : public EigenDescriptor {
public:
    SingularValueEntropy();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class SingularValueVariance : public EigenDescriptor {
public:
    SingularValueVariance();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class SingularValueSkewness : public EigenDescriptor {
public:
    SingularValueSkewness();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class SpectralEffectiveRank : public EigenDescriptor {
public:
    SpectralEffectiveRank();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class NuclearNorm : public EigenDescriptor {
public:
    NuclearNorm();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

// Normalized adjacency matrix descriptors
class NormalizedAdjacencySpectralRadius : public EigenDescriptor {
public:
    NormalizedAdjacencySpectralRadius();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class NormalizedEigenvalueGap : public EigenDescriptor {
public:
    NormalizedEigenvalueGap();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class SumNormalizedEigenvalues : public EigenDescriptor {
public:
    SumNormalizedEigenvalues();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class VarianceNormalizedEigenvalues : public EigenDescriptor {
public:
    VarianceNormalizedEigenvalues();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class CountNormalizedEigenvaluesAboveHalf : public EigenDescriptor {
public:
    CountNormalizedEigenvaluesAboveHalf();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class NormalizedEnergy : public EigenDescriptor {
public:
    NormalizedEnergy();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class LargestNormalizedEigenvectorCentrality : public EigenDescriptor {
public:
    LargestNormalizedEigenvectorCentrality();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class AverageNormalizedEigenvectorCentrality : public EigenDescriptor {
public:
    AverageNormalizedEigenvectorCentrality();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class NormalizedAdjacencyMatrixRank : public EigenDescriptor {
public:
    NormalizedAdjacencyMatrixRank();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class NormalizedAdjacencyEntropy : public EigenDescriptor {
public:
    NormalizedAdjacencyEntropy();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

// Signless Laplacian matrix descriptors
class SignlessLaplacianSpectralRadius : public EigenDescriptor {
public:
    SignlessLaplacianSpectralRadius();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class SignlessLaplacianSmallestEigenvalue : public EigenDescriptor {
public:
    SignlessLaplacianSmallestEigenvalue();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class SignlessLaplacianEnergy : public EigenDescriptor {
public:
    SignlessLaplacianEnergy();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class SignlessLaplacianTrace : public EigenDescriptor {
public:
    SignlessLaplacianTrace();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class SignlessLaplacianDeterminant : public EigenDescriptor {
public:
    SignlessLaplacianDeterminant();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class SignlessLaplacianZeroEigenvalues : public EigenDescriptor {
public:
    SignlessLaplacianZeroEigenvalues();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class SignlessLaplacianPositiveEigenvalues : public EigenDescriptor {
public:
    SignlessLaplacianPositiveEigenvalues();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class SignlessLaplacianNegativeEigenvalues : public EigenDescriptor {
public:
    SignlessLaplacianNegativeEigenvalues();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class SignlessLaplacianEigenvalueVariance : public EigenDescriptor {
public:
    SignlessLaplacianEigenvalueVariance();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

// class SignlessLaplacianEigenvalueSkewness : public EigenDescriptor {
// public:
//     SignlessLaplacianEigenvalueSkewness();
//     std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
// }; // Commented out due to zero or near-zero variance

class SignlessLaplacianEigenvalueSkewness : public EigenDescriptor {
public:
    SignlessLaplacianEigenvalueSkewness();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

// Miscellaneous graph matrices descriptors
class WeightedAdjacencySpectralRadius : public EigenDescriptor {
public:
    WeightedAdjacencySpectralRadius();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class MeanFirstPassageTime : public EigenDescriptor {
public:
    MeanFirstPassageTime();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class CommuteTimeDistance : public EigenDescriptor {
public:
    CommuteTimeDistance();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class KirchhoffIndexVariance : public EigenDescriptor {
public:
    KirchhoffIndexVariance();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class EffectiveGraphResistanceDistribution : public EigenDescriptor {
public:
    EffectiveGraphResistanceDistribution();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class LocalClusteringCoefficientDistribution : public EigenDescriptor {
public:
    LocalClusteringCoefficientDistribution();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class GraphRobustnessIndex : public EigenDescriptor {
public:
    GraphRobustnessIndex();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class NormalizedEstradaIndex : public EigenDescriptor {
public:
    NormalizedEstradaIndex();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class GraphBipartivityIndex : public EigenDescriptor {
public:
    GraphBipartivityIndex();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class SpanningTreeEntropy : public EigenDescriptor {
public:
    SpanningTreeEntropy();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

// Distance Matrix Based (DistMat)
class DistMatSpecRad : public EigenDescriptor {
public: DistMatSpecRad(); std::variant<double, int, std::string> calculate(const Molecule& mol) const override; };
class DistMatMinEig : public EigenDescriptor {
public: DistMatMinEig(); std::variant<double, int, std::string> calculate(const Molecule& mol) const override; };
class DistMatEigSum : public EigenDescriptor {
public: DistMatEigSum(); std::variant<double, int, std::string> calculate(const Molecule& mol) const override; };
class DistMatEigVar : public EigenDescriptor {
public: DistMatEigVar(); std::variant<double, int, std::string> calculate(const Molecule& mol) const override; };
class DistMatEigSkew : public EigenDescriptor {
public: DistMatEigSkew(); std::variant<double, int, std::string> calculate(const Molecule& mol) const override; };
class DistMatEigKurt : public EigenDescriptor {
public: DistMatEigKurt(); std::variant<double, int, std::string> calculate(const Molecule& mol) const override; };
class DistMatEnergy : public EigenDescriptor {
public: DistMatEnergy(); std::variant<double, int, std::string> calculate(const Molecule& mol) const override; };
class DistMatFrobNorm : public EigenDescriptor {
public: DistMatFrobNorm(); std::variant<double, int, std::string> calculate(const Molecule& mol) const override; };
class DistMatRank : public EigenDescriptor {
public: DistMatRank(); std::variant<double, int, std::string> calculate(const Molecule& mol) const override; };
class DistMatDet : public EigenDescriptor {
public: DistMatDet(); std::variant<double, int, std::string> calculate(const Molecule& mol) const override; };
class DistMatConditionNumber : public EigenDescriptor {
public: DistMatConditionNumber(); std::variant<double, int, std::string> calculate(const Molecule& mol) const override; };
class DistMatSingValEntropy : public EigenDescriptor {
public: DistMatSingValEntropy(); std::variant<double, int, std::string> calculate(const Molecule& mol) const override; };
class DistMatSpectralEffRank : public EigenDescriptor {
public: DistMatSpectralEffRank(); std::variant<double, int, std::string> calculate(const Molecule& mol) const override; };

// Atomic Number Weighted Adjacency (WtNumAdj)
class WtNumAdjSpecRad : public EigenDescriptor {
public: WtNumAdjSpecRad(); std::variant<double, int, std::string> calculate(const Molecule& mol) const override; };
class WtNumAdjMinEig : public EigenDescriptor {
public: WtNumAdjMinEig(); std::variant<double, int, std::string> calculate(const Molecule& mol) const override; };
class WtNumAdjEigSum : public EigenDescriptor {
public: WtNumAdjEigSum(); std::variant<double, int, std::string> calculate(const Molecule& mol) const override; };
class WtNumAdjEigVar : public EigenDescriptor {
public: WtNumAdjEigVar(); std::variant<double, int, std::string> calculate(const Molecule& mol) const override; };
class WtNumAdjEnergy : public EigenDescriptor {
public: WtNumAdjEnergy(); std::variant<double, int, std::string> calculate(const Molecule& mol) const override; };

// Atomic Number Weighted Laplacian (WtNumLap)
class WtNumLapSpecRad : public EigenDescriptor {
public: WtNumLapSpecRad(); std::variant<double, int, std::string> calculate(const Molecule& mol) const override; };
class WtNumLapAlgConn : public EigenDescriptor {
public: WtNumLapAlgConn(); std::variant<double, int, std::string> calculate(const Molecule& mol) const override; };
class WtNumLapEigSum : public EigenDescriptor {
public: WtNumLapEigSum(); std::variant<double, int, std::string> calculate(const Molecule& mol) const override; };
class WtNumLapEigVar : public EigenDescriptor {
public: WtNumLapEigVar(); std::variant<double, int, std::string> calculate(const Molecule& mol) const override; };
class WtNumLapEnergy : public EigenDescriptor {
public: WtNumLapEnergy(); std::variant<double, int, std::string> calculate(const Molecule& mol) const override; };

// Adjacency Matrix Powers (A^k) Spectral/Energy
class AdjPow2SpecRad : public EigenDescriptor {
public: AdjPow2SpecRad(); std::variant<double, int, std::string> calculate(const Molecule& mol) const override; };
class AdjPow2Energy : public EigenDescriptor {
public: AdjPow2Energy(); std::variant<double, int, std::string> calculate(const Molecule& mol) const override; };
class AdjPow3SpecRad : public EigenDescriptor {
public: AdjPow3SpecRad(); std::variant<double, int, std::string> calculate(const Molecule& mol) const override; };
class AdjPow3Energy : public EigenDescriptor {
public: AdjPow3Energy(); std::variant<double, int, std::string> calculate(const Molecule& mol) const override; };
class AdjPow4SpecRad : public EigenDescriptor {
public: AdjPow4SpecRad(); std::variant<double, int, std::string> calculate(const Molecule& mol) const override; };
class AdjPow4Energy : public EigenDescriptor {
public: AdjPow4Energy(); std::variant<double, int, std::string> calculate(const Molecule& mol) const override; };

// Laplacian Matrix Powers (L^k) Spectral/Energy/Trace
class LapPow2SpecRad : public EigenDescriptor {
public: LapPow2SpecRad(); std::variant<double, int, std::string> calculate(const Molecule& mol) const override; };
class LapPow2AlgConn : public EigenDescriptor {
public: LapPow2AlgConn(); std::variant<double, int, std::string> calculate(const Molecule& mol) const override; };
class LapPow2Energy : public EigenDescriptor {
public: LapPow2Energy(); std::variant<double, int, std::string> calculate(const Molecule& mol) const override; };
class LapPow2Trace : public EigenDescriptor {
public: LapPow2Trace(); std::variant<double, int, std::string> calculate(const Molecule& mol) const override; };

// Matrix Exponential based (Estrada Index for other matrices)
// class LapEstradaIdx : public EigenDescriptor {
// public: LapEstradaIdx(); std::variant<double, int, std::string> calculate(const Molecule& mol) const override; }; // Commented out due to zero or near-zero variance
// class SignLapEstradaIdx : public EigenDescriptor {
// public: SignLapEstradaIdx(); std::variant<double, int, std::string> calculate(const Molecule& mol) const override; }; // Commented out due to zero or near-zero variance
class DistMatEstradaIdx : public EigenDescriptor {
public: DistMatEstradaIdx(); std::variant<double, int, std::string> calculate(const Molecule& mol) const override; };
class WtNumAdjEstradaIdx : public EigenDescriptor {
public: WtNumAdjEstradaIdx(); std::variant<double, int, std::string> calculate(const Molecule& mol) const override; };
class WtNumLapEstradaIdx : public EigenDescriptor {
public: WtNumLapEstradaIdx(); std::variant<double, int, std::string> calculate(const Molecule& mol) const override; };

// Eigenvector Statistics
class AdjDomEigVecVar : public EigenDescriptor {
public: AdjDomEigVecVar(); std::variant<double, int, std::string> calculate(const Molecule& mol) const override; };
class AdjDomEigVecSkew : public EigenDescriptor {
public: AdjDomEigVecSkew(); std::variant<double, int, std::string> calculate(const Molecule& mol) const override; };
class AdjDomEigVecKurt : public EigenDescriptor {
public: AdjDomEigVecKurt(); std::variant<double, int, std::string> calculate(const Molecule& mol) const override; };
class NormLapDomEigVecVar : public EigenDescriptor {
public: NormLapDomEigVecVar(); std::variant<double, int, std::string> calculate(const Molecule& mol) const override; };
class SignLapDomEigVecVar : public EigenDescriptor {
public: SignLapDomEigVecVar(); std::variant<double, int, std::string> calculate(const Molecule& mol) const override; };

// Miscellaneous Eigen-based
class SumFormanRicci : public EigenDescriptor {
public: SumFormanRicci(); std::variant<double, int, std::string> calculate(const Molecule& mol) const override; };
class LogDetLaplacian : public EigenDescriptor {
public: LogDetLaplacian(); std::variant<double, int, std::string> calculate(const Molecule& mol) const override; };
class AvgDegreeDistance : public EigenDescriptor {
public: AvgDegreeDistance(); std::variant<double, int, std::string> calculate(const Molecule& mol) const override; };
class VarianceDegreeDistance : public EigenDescriptor {
public: VarianceDegreeDistance(); std::variant<double, int, std::string> calculate(const Molecule& mol) const override; };

} // namespace descriptors
} // namespace desfact 