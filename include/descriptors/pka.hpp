// include/descriptors/pka.hpp
#pragma once

#include "descriptors.hpp"
#include <GraphMol/GraphMol.h>
#include <GraphMol/Atom.h>
#include <GraphMol/SmilesParse/SmartsWrite.h> // Include necessary headers
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <string>
#include <variant>
#include <vector>          // Needed for RandomForestModel
#include <unordered_map>
#include <memory> // For unique_ptr

// Forward declare RDKit types
namespace RDKit {
    class ROMol;
    class RWMol; // Needed for SmartsToMol potentially
}

namespace desfact {
namespace descriptors {

// Moved definitions from pka.cpp
struct TreeNode {
    int feature;
    float threshold;
    int left_child;
    int right_child;
    float value;
};

struct RandomForestModel {
    int n_estimators;
    int n_features;
    std::vector<std::vector<TreeNode>> trees;
    std::unordered_map<std::string, int> vocabulary;
    std::vector<float> idf;
    int analyzer_type; // 0=char, 1=word
    int ngram_min;
    int ngram_max;
    
    // Keep predict definition here if simple, or move implementation to cpp
    float predict(const std::vector<float>& features) const {
        float sum = 0.0f;
        for (const auto& tree : trees) {
            int node_id = 0;
            while (true) {
                const TreeNode& node = tree[node_id];
                if (node.left_child == -1) { // Leaf node
                    sum += node.value;
                    break;
                }
                // Ensure features access is safe
                if (node.feature >= 0 && static_cast<size_t>(node.feature) < features.size()) {
                     if (features[node.feature] <= node.threshold) {
                         node_id = node.left_child;
                     } else {
                         node_id = node.right_child;
                     }
                } else {
                    // Handle error: invalid feature index
                    // For simplicity, let's break (consider logging/throwing)
                     // globalLogger.error("Invalid feature index in RF predict: " + std::to_string(node.feature));
                     break;
                }
            }
        }
        return (n_estimators > 0) ? (sum / n_estimators) : 0.0f; // Avoid division by zero
    }
};

// Common utility functions for PKA descriptors
namespace pka_utils {
    // Removed placeholder declarations, implementations are within PkaDescriptor or anonymous namespace
}

// Base class for pKa related descriptors (optional, but can group helpers)
class PkaDescriptor : public Descriptor {
protected:
    // Helper functions (to be defined in pka.cpp)
    static std::string classifyMoleculeType(const RDKit::ROMol& mol);
    static double estimateAcidPKa(const RDKit::ROMol& mol);
    static double estimateBasePKa(const RDKit::ROMol& mol);

    // Cache for SMARTS patterns (potentially static thread_local in cpp)
    // Use unique_ptr for automatic memory management
    struct PatternCache {
         std::unordered_map<std::string, std::unique_ptr<RDKit::ROMol>> patterns;
         // ~PatternCache(); // Destructor no longer needed with unique_ptr
         RDKit::ROMol* getPattern(const std::string& smarts); // Return raw pointer, cache owns it
    };
    // Make cache thread_local for thread safety
    static thread_local PatternCache s_patternCache;

    // Cache for ML models
    class ModelCache {
    public:
        const RandomForestModel* getModel(const std::string& modelPath);
        void clear() { models.clear(); }
    private:
        std::unordered_map<std::string, std::unique_ptr<RandomForestModel>> models;
    };
    static thread_local ModelCache s_modelCache;

public:
    PkaDescriptor(const std::string& name, const std::string& description);
    ~PkaDescriptor() override = default; // Ensure virtual destructor

    // Static method to clear thread-local caches
    static void clearThreadLocalCaches();
};

// Descriptor for acidic pKa
class AcidicPKADescriptor : public Descriptor {
public:
    AcidicPKADescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

// Descriptor for basic pKa
class BasicPKADescriptor : public Descriptor {
public:
    BasicPKADescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

// Descriptor for compound classification
class CompoundClassDescriptor : public Descriptor {
public:
    CompoundClassDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

// --- Specific pKa Descriptor Classes ---

class AcidityTypeDescriptor : public PkaDescriptor {
public:
    AcidityTypeDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class PkaEstimate1Descriptor : public PkaDescriptor {
public:
    PkaEstimate1Descriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class PkaEstimate2Descriptor : public PkaDescriptor {
public:
    PkaEstimate2Descriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

// --- Classification Boolean Descriptors ---

class IsAcidDescriptor : public PkaDescriptor {
public:
    IsAcidDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class IsBaseDescriptor : public PkaDescriptor {
public:
    IsBaseDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class IsNeutralDescriptor : public PkaDescriptor {
public:
    IsNeutralDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class IsZwitterionDescriptor : public PkaDescriptor {
public:
    IsZwitterionDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

} // namespace descriptors
} // namespace desfact