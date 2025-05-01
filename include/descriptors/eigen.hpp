#pragma once

#include "descriptors.hpp"
#include "descriptors/eigen_detail.hpp"
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <GraphMol/Atom.h>
#include <GraphMol/ROMol.h>
#include <variant>
#include <string>
#include <cmath>
#include <chrono> // For cache access time
#include <unordered_map> // For caches
#include <vector> // For sorting cache entries
#include <algorithm> // For std::sort, std::min

namespace desfact {
namespace descriptors {

// Base class for Eigen-based descriptors
class EigenDescriptor : public Descriptor {
protected:
    // --- Caching Structures ---
    // Maximum cache size constants
    static constexpr size_t MAX_CACHE_ENTRIES = 500;  // Adjusted max entries per cache
    static constexpr size_t CACHE_CLEANUP_THRESHOLD = 400; // Trigger cleanup earlier

    // Base cache entry struct
    struct CacheEntry {
        std::chrono::steady_clock::time_point lastAccess;
    };

    // Matrix cache entry: Stores the computed matrix
    struct MatrixCacheEntry : CacheEntry {
        Eigen::MatrixXd matrix;
    };

    // Eigenvalue cache entry: Stores computed eigenvalues
    struct EigenvalueCacheEntry : CacheEntry {
        Eigen::VectorXd eigenvalues;
    };

    // --- Thread-Local Caches ---
    // Use thread_local for thread safety without explicit locking for reads/writes *within* a thread.
    // Keyed by the const RDKit::ROMol pointer (assumes molecule object lifetime exceeds cache usage within a thread's task)
    static thread_local std::unordered_map<const RDKit::ROMol*, MatrixCacheEntry> adjacencyCache;
    static thread_local std::unordered_map<const RDKit::ROMol*, MatrixCacheEntry> degreeCache;
    static thread_local std::unordered_map<const RDKit::ROMol*, MatrixCacheEntry> laplacianCache;
    static thread_local std::unordered_map<const RDKit::ROMol*, MatrixCacheEntry> normalizedLaplacianCache;
    static thread_local std::unordered_map<const RDKit::ROMol*, MatrixCacheEntry> signlessLaplacianCache;
    static thread_local std::unordered_map<const RDKit::ROMol*, MatrixCacheEntry> weightedAdjacencyCache;
    static thread_local std::unordered_map<const RDKit::ROMol*, MatrixCacheEntry> distanceMatrixCache;
    static thread_local std::unordered_map<const RDKit::ROMol*, MatrixCacheEntry> atomicNumWeightedAdjCache;
    static thread_local std::unordered_map<const RDKit::ROMol*, MatrixCacheEntry> atomicNumWeightedLapCache;
    static thread_local std::unordered_map<const RDKit::ROMol*, MatrixCacheEntry> normalizedAdjCache;
    static thread_local std::unordered_map<const RDKit::ROMol*, MatrixCacheEntry> adjPow2Cache;
    static thread_local std::unordered_map<const RDKit::ROMol*, MatrixCacheEntry> adjPow3Cache;
    static thread_local std::unordered_map<const RDKit::ROMol*, MatrixCacheEntry> adjPow4Cache;
    static thread_local std::unordered_map<const RDKit::ROMol*, MatrixCacheEntry> lapPow2Cache;

    // --- Eigenvalue Cache ---
    // Hash function for Eigen matrices (simple version based on size and corner values)
    struct MatrixHash {
        std::size_t operator()(const Eigen::MatrixXd& matrix) const noexcept {
            std::size_t seed = 0;
            // Combine hash of rows and cols
            std::hash<Eigen::Index> indexHasher;
            seed ^= indexHasher(matrix.rows()) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
            seed ^= indexHasher(matrix.cols()) + 0x9e3779b9 + (seed << 6) + (seed >> 2);

            // Combine hash of a few corner/center elements for content check
            std::hash<double> doubleHasher;
            int rows = matrix.rows();
            int cols = matrix.cols();
            if (rows > 0 && cols > 0) {
                 seed ^= doubleHasher(matrix(0, 0)) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
                 if (rows > 1 && cols > 1)
                    seed ^= doubleHasher(matrix(rows-1, cols-1)) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
                 if (rows > 2 && cols > 2)
                     seed ^= doubleHasher(matrix(rows/2, cols/2)) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
            }
            return seed;
        }
    };

    // Equality check for Eigen matrices (needed for hash map)
    struct MatrixEqual {
        bool operator()(const Eigen::MatrixXd& a, const Eigen::MatrixXd& b) const noexcept {
            // Check dimensions first
            if (a.rows() != b.rows() || a.cols() != b.cols()) return false;
            // Use Eigen's built-in approximate comparison
            return a.isApprox(b, 1e-9); // Tolerance for floating point comparison
        }
    };
    // Cache for eigenvalues, keyed by the matrix itself using the custom hash and equality functors
    static thread_local std::unordered_map<Eigen::MatrixXd, EigenvalueCacheEntry, MatrixHash, MatrixEqual> eigenvalueCache;
    // Cache for singular values
    static thread_local std::unordered_map<Eigen::MatrixXd, EigenvalueCacheEntry, MatrixHash, MatrixEqual> singularValueCache; // Using EigenvalueCacheEntry as structure is same
    // Cache for dominant eigenvectors
    static thread_local std::unordered_map<Eigen::MatrixXd, EigenvalueCacheEntry, MatrixHash, MatrixEqual> dominantEigenvectorCache; // Use EigenvalueCacheEntry for VectorXd


    // --- Cache Management ---
    template<typename CacheType>
    static void cleanupCache(CacheType& cache, size_t keepEntries) {
        if (cache.size() <= keepEntries) return;

        // Create a vector of pairs: <key, lastAccessTime>
        std::vector<std::pair<typename CacheType::key_type, std::chrono::steady_clock::time_point>> entries;
        entries.reserve(cache.size());
        for (const auto& pair : cache) {
            entries.emplace_back(pair.first, pair.second.lastAccess);
        }

        // Sort by access time (oldest first)
        std::sort(entries.begin(), entries.end(),
                  [](const auto& a, const auto& b) { return a.second < b.second; });

        // Remove the oldest entries until the cache size is reduced to 'keepEntries'
        size_t removeCount = cache.size() - keepEntries;
        for (size_t i = 0; i < removeCount; ++i) {
            cache.erase(entries[i].first);
        }
    }

    // Call cleanup on all caches if any exceed threshold
    static void checkAndCleanupCaches() {
        // Check matrix caches
        if (adjacencyCache.size() > CACHE_CLEANUP_THRESHOLD) cleanupCache(adjacencyCache, MAX_CACHE_ENTRIES / 2);
        if (degreeCache.size() > CACHE_CLEANUP_THRESHOLD) cleanupCache(degreeCache, MAX_CACHE_ENTRIES / 2);
        if (laplacianCache.size() > CACHE_CLEANUP_THRESHOLD) cleanupCache(laplacianCache, MAX_CACHE_ENTRIES / 2);
        if (normalizedLaplacianCache.size() > CACHE_CLEANUP_THRESHOLD) cleanupCache(normalizedLaplacianCache, MAX_CACHE_ENTRIES / 2);
        if (signlessLaplacianCache.size() > CACHE_CLEANUP_THRESHOLD) cleanupCache(signlessLaplacianCache, MAX_CACHE_ENTRIES / 2);
        if (weightedAdjacencyCache.size() > CACHE_CLEANUP_THRESHOLD) cleanupCache(weightedAdjacencyCache, MAX_CACHE_ENTRIES / 2);
        if (distanceMatrixCache.size() > CACHE_CLEANUP_THRESHOLD) cleanupCache(distanceMatrixCache, MAX_CACHE_ENTRIES / 2);
        if (atomicNumWeightedAdjCache.size() > CACHE_CLEANUP_THRESHOLD) cleanupCache(atomicNumWeightedAdjCache, MAX_CACHE_ENTRIES / 2);
        if (atomicNumWeightedLapCache.size() > CACHE_CLEANUP_THRESHOLD) cleanupCache(atomicNumWeightedLapCache, MAX_CACHE_ENTRIES / 2);
        if (normalizedAdjCache.size() > CACHE_CLEANUP_THRESHOLD) cleanupCache(normalizedAdjCache, MAX_CACHE_ENTRIES / 2);
        if (adjPow2Cache.size() > CACHE_CLEANUP_THRESHOLD) cleanupCache(adjPow2Cache, MAX_CACHE_ENTRIES / 2);
        if (adjPow3Cache.size() > CACHE_CLEANUP_THRESHOLD) cleanupCache(adjPow3Cache, MAX_CACHE_ENTRIES / 2);
        if (adjPow4Cache.size() > CACHE_CLEANUP_THRESHOLD) cleanupCache(adjPow4Cache, MAX_CACHE_ENTRIES / 2);
        if (lapPow2Cache.size() > CACHE_CLEANUP_THRESHOLD) cleanupCache(lapPow2Cache, MAX_CACHE_ENTRIES / 2);


        // Check vector caches
        if (eigenvalueCache.size() > CACHE_CLEANUP_THRESHOLD) cleanupCache(eigenvalueCache, MAX_CACHE_ENTRIES / 2);
        if (singularValueCache.size() > CACHE_CLEANUP_THRESHOLD) cleanupCache(singularValueCache, MAX_CACHE_ENTRIES / 2);
        if (dominantEigenvectorCache.size() > CACHE_CLEANUP_THRESHOLD) cleanupCache(dominantEigenvectorCache, MAX_CACHE_ENTRIES / 2);
    }

    // --- Cached Getters ---
    // Template for getting matrix from cache or building it
    template<typename CacheMap, typename BuildFunc>
    const Eigen::MatrixXd& getOrBuildMatrix(const RDKit::ROMol* mol, CacheMap& cache, BuildFunc builder) const {
        if (!mol) {
            static const Eigen::MatrixXd empty_matrix(0, 0);
            return empty_matrix;
        }

        auto it = cache.find(mol);
        if (it != cache.end()) {
            if (!detail::isMatrixSizeValid(it->second.matrix)) {
                cache.erase(it); // Remove invalid cached matrix
            } else {
                it->second.lastAccess = std::chrono::steady_clock::now();
                return it->second.matrix;
            }
        }

        checkAndCleanupCaches(); // Check before building potentially large matrix
        Eigen::MatrixXd matrix = builder(mol);
        
        // Validate the built matrix
        if (!detail::isMatrixSizeValid(matrix)) {
            static const Eigen::MatrixXd empty_matrix(0, 0);
            return empty_matrix;
        }

        cache[mol] = MatrixCacheEntry{{std::chrono::steady_clock::now()}, std::move(matrix)};
        return cache[mol].matrix;
    }

    // Template for getting vector (Eigenvalues, SingularValues, Eigenvectors) from cache or computing
    template<typename CacheMap, typename ComputeFunc>
    const Eigen::VectorXd& getOrComputeVector(const Eigen::MatrixXd& matrix, CacheMap& cache, ComputeFunc computer) const {
        auto it = cache.find(matrix);
        if (it != cache.end()) {
            it->second.lastAccess = std::chrono::steady_clock::now();
            return it->second.eigenvalues;
        }

        checkAndCleanupCaches();
        Eigen::VectorXd vec = computer(matrix);
        cache[matrix] = EigenvalueCacheEntry{{std::chrono::steady_clock::now()}, std::move(vec)};
        return cache[matrix].eigenvalues;
    }

    // Specific getters using the templates
    const Eigen::MatrixXd& getAdjacencyMatrix(const RDKit::ROMol* mol) const {
        return getOrBuildMatrix(mol, adjacencyCache, [this](const RDKit::ROMol* m){ return buildAdjacencyMatrix(m); });
    }
    const Eigen::MatrixXd& getDegreeMatrix(const RDKit::ROMol* mol) const {
        return getOrBuildMatrix(mol, degreeCache, [this](const RDKit::ROMol* m){
            Eigen::MatrixXd adj = getAdjacencyMatrix(m); // Get potentially cached adj matrix
            return buildDegreeMatrix(adj);
        });
    }
     const Eigen::MatrixXd& getLaplacianMatrix(const RDKit::ROMol* mol) const {
         return getOrBuildMatrix(mol, laplacianCache, [this](const RDKit::ROMol* m){
             Eigen::MatrixXd adj = getAdjacencyMatrix(m);
             Eigen::MatrixXd deg = getDegreeMatrix(m); // Uses cache internally
             return buildLaplacianMatrix(deg, adj);
         });
     }
     const Eigen::MatrixXd& getNormalizedLaplacianMatrix(const RDKit::ROMol* mol) const {
         return getOrBuildMatrix(mol, normalizedLaplacianCache, [this](const RDKit::ROMol* m){
             Eigen::MatrixXd adj = getAdjacencyMatrix(m);
             Eigen::MatrixXd deg = getDegreeMatrix(m);
             return buildNormalizedLaplacian(deg, adj);
         });
     }
     const Eigen::MatrixXd& getSignlessLaplacianMatrix(const RDKit::ROMol* mol) const {
          return getOrBuildMatrix(mol, signlessLaplacianCache, [this](const RDKit::ROMol* m){
              Eigen::MatrixXd adj = getAdjacencyMatrix(m);
              Eigen::MatrixXd deg = getDegreeMatrix(m);
              return buildSignlessLaplacian(deg, adj);
          });
      }
     const Eigen::MatrixXd& getWeightedAdjacencyMatrix(const RDKit::ROMol* mol) const {
          return getOrBuildMatrix(mol, weightedAdjacencyCache, [this](const RDKit::ROMol* m){
              return buildWeightedAdjacency(m);
          });
      }
     const Eigen::MatrixXd& getDistanceMatrix(const RDKit::ROMol* mol) const {
          return getOrBuildMatrix(mol, distanceMatrixCache, [this](const RDKit::ROMol* m){
              return buildDistanceMatrix(m);
          });
      }
     const Eigen::MatrixXd& getAtomicNumberWeightedAdjacency(const RDKit::ROMol* mol) const {
          return getOrBuildMatrix(mol, atomicNumWeightedAdjCache, [this](const RDKit::ROMol* m){
              return buildAtomicNumberWeightedAdjacency(m);
          });
      }
      const Eigen::MatrixXd& getAtomicNumberWeightedLaplacian(const RDKit::ROMol* mol) const {
          return getOrBuildMatrix(mol, atomicNumWeightedLapCache, [this](const RDKit::ROMol* m){
              return buildAtomicNumberWeightedLaplacian(m);
          });
      }
      const Eigen::MatrixXd& getNormalizedAdjacency(const RDKit::ROMol* mol) const {
           return getOrBuildMatrix(mol, normalizedAdjCache, [this](const RDKit::ROMol* m){
               Eigen::MatrixXd adj = getAdjacencyMatrix(m);
               Eigen::MatrixXd deg = getDegreeMatrix(m);
               return buildNormalizedAdjacency(deg, adj);
           });
       }
       const Eigen::MatrixXd& getAdjPow2(const RDKit::ROMol* mol) const {
            return getOrBuildMatrix(mol, adjPow2Cache, [this](const RDKit::ROMol* m){
                const Eigen::MatrixXd& adj = getAdjacencyMatrix(m);
                return detail::safeMatrixMultiply(adj, adj);
            });
        }
        const Eigen::MatrixXd& getAdjPow3(const RDKit::ROMol* mol) const {
             return getOrBuildMatrix(mol, adjPow3Cache, [this](const RDKit::ROMol* m){
                 const Eigen::MatrixXd& adj2 = getAdjPow2(m);
                 const Eigen::MatrixXd& adj = getAdjacencyMatrix(m);
                 return detail::safeMatrixMultiply(adj2, adj);
             });
         }
        const Eigen::MatrixXd& getAdjPow4(const RDKit::ROMol* mol) const {
             return getOrBuildMatrix(mol, adjPow4Cache, [this](const RDKit::ROMol* m){
                 const Eigen::MatrixXd& adj2 = getAdjPow2(m);
                 return detail::safeMatrixMultiply(adj2, adj2);
             });
         }
        const Eigen::MatrixXd& getLapPow2(const RDKit::ROMol* mol) const {
             return getOrBuildMatrix(mol, lapPow2Cache, [this](const RDKit::ROMol* m){
                 Eigen::MatrixXd lap = getLaplacianMatrix(m); // Use cached L
                 return lap * lap;
             });
         }


    // Get cached or compute eigenvalues
    const Eigen::VectorXd& getEigenvalues(const Eigen::MatrixXd& matrix) const {
        return getOrComputeVector(matrix, eigenvalueCache, [this](const Eigen::MatrixXd& m){ return computeEigenvalues(m); });
    }
    // Get cached or compute singular values
    const Eigen::VectorXd& getSingularValues(const Eigen::MatrixXd& matrix) const {
        return getOrComputeVector(matrix, singularValueCache, [this](const Eigen::MatrixXd& m){ return computeSingularValues(m); });
    }
     // Get cached or compute dominant eigenvector
    const Eigen::VectorXd& getDominantEigenvector(const Eigen::MatrixXd& matrix) const {
        return getOrComputeVector(matrix, dominantEigenvectorCache, [this](const Eigen::MatrixXd& m){ return computeDominantEigenvector(m); });
    }


    // --- Original Matrix Building Methods (now called by caching layer if needed) ---
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

    // --- Compute Eigenvalues/vectors/SVD (now called by caching layer if needed) ---
    Eigen::VectorXd computeEigenvalues(const Eigen::MatrixXd& matrix) const;
    Eigen::VectorXd computeSingularValues(const Eigen::MatrixXd& matrix) const;
    Eigen::VectorXd computeDominantEigenvector(const Eigen::MatrixXd& matrix) const; // Renamed from getDominantEigenvector


    // --- Statistical Helpers ---
    double computeVariance(const Eigen::VectorXd& values) const;
    double computeSkewness(const Eigen::VectorXd& values) const;
    double computeKurtosis(const Eigen::VectorXd& values) const;
    double computeEntropy(const Eigen::VectorXd& values) const; // Add definition in .cpp if used

public:
    EigenDescriptor(const std::string& name, const std::string& description);
    ~EigenDescriptor() override = default; // Ensure virtual destructor

    // Static method to clear all thread-local caches (e.g., call at end of thread execution)
    static void clearThreadLocalCaches();
};

// --- Descriptor Classes (Inherit from EigenDescriptor) ---
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