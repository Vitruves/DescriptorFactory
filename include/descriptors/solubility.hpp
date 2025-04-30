// ... existing code ...
#pragma once

#include "descriptors.hpp"
#include <GraphMol/GraphMol.h>
#include <string>
#include <variant>
#include <vector>
#include <unordered_map>
#include <memory>

namespace desfact {
namespace descriptors {

struct SolubilityTreeNode {
    int feature;
    float threshold;
    int left_child;
    int right_child;
    float value;
};

struct SolubilityRandomForestModel {
    int n_estimators;
    int n_features;
    std::vector<std::vector<SolubilityTreeNode>> trees;
    std::unordered_map<std::string, int> vocabulary;
    std::vector<float> idf;
    int analyzer_type;
    int ngram_min;
    int ngram_max;
    float predict(const std::vector<float>& features) const;
};

class SolubilityDescriptor : public Descriptor {
protected:
    class ModelCache {
    public:
        const SolubilityRandomForestModel* getModel(const std::string& modelPath);
        void clear() { models.clear(); }
    private:
        std::unordered_map<std::string, std::unique_ptr<SolubilityRandomForestModel>> models;
    };
    static thread_local ModelCache s_modelCache;
public:
    SolubilityDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
    static void clearThreadLocalCaches();
};

} // namespace descriptors
} // namespace desfact