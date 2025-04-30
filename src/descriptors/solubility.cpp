#include "descriptors/solubility.hpp"
#include <GraphMol/GraphMol.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <fstream>
#include <limits>
#include <memory>
#include "utils.hpp"

namespace desfact {
namespace descriptors {

thread_local SolubilityDescriptor::ModelCache SolubilityDescriptor::s_modelCache;

static SolubilityRandomForestModel loadSolubilityModel(const std::string& filename) {
    std::ifstream file(filename, std::ios::binary);
    if (!file) throw std::runtime_error("Cannot open solubility model file: " + filename);

    SolubilityRandomForestModel model;
    char magic[6] = {0};
    file.read(magic, 5);
    if (std::string(magic) != "RFBIN") throw std::runtime_error("Invalid solubility model file format");

    uint16_t version;
    file.read(reinterpret_cast<char*>(&version), sizeof(version));
    if (version != 1) throw std::runtime_error("Unsupported solubility model version");

    file.read(reinterpret_cast<char*>(&model.n_estimators), sizeof(model.n_estimators));
    file.read(reinterpret_cast<char*>(&model.n_features), sizeof(model.n_features));

    uint32_t vocab_size;
    file.read(reinterpret_cast<char*>(&vocab_size), sizeof(vocab_size));
    for (uint32_t i = 0; i < vocab_size; i++) {
        uint32_t word_len;
        file.read(reinterpret_cast<char*>(&word_len), sizeof(word_len));
        std::string word(word_len, ' ');
        file.read(&word[0], word_len);
        uint32_t idx;
        file.read(reinterpret_cast<char*>(&idx), sizeof(idx));
        model.vocabulary[word] = idx;
    }

    uint32_t idf_size;
    file.read(reinterpret_cast<char*>(&idf_size), sizeof(idf_size));
    model.idf.resize(idf_size);
    for (uint32_t i = 0; i < idf_size; i++) {
        file.read(reinterpret_cast<char*>(&model.idf[i]), sizeof(float));
    }

    uint8_t analyzer_type_byte;
    file.read(reinterpret_cast<char*>(&analyzer_type_byte), sizeof(analyzer_type_byte));
    model.analyzer_type = static_cast<int>(analyzer_type_byte);
    file.read(reinterpret_cast<char*>(&model.ngram_min), sizeof(model.ngram_min));
    file.read(reinterpret_cast<char*>(&model.ngram_max), sizeof(model.ngram_max));

    model.trees.resize(model.n_estimators);
    for (int tree_idx = 0; tree_idx < model.n_estimators; tree_idx++) {
        uint32_t node_count;
        file.read(reinterpret_cast<char*>(&node_count), sizeof(node_count));
        model.trees[tree_idx].resize(node_count);
        for (uint32_t node_idx = 0; node_idx < node_count; node_idx++) {
            SolubilityTreeNode& node = model.trees[tree_idx][node_idx];
            file.read(reinterpret_cast<char*>(&node.feature), sizeof(node.feature));
            file.read(reinterpret_cast<char*>(&node.threshold), sizeof(node.threshold));
            file.read(reinterpret_cast<char*>(&node.left_child), sizeof(node.left_child));
            file.read(reinterpret_cast<char*>(&node.right_child), sizeof(node.right_child));
            file.read(reinterpret_cast<char*>(&node.value), sizeof(node.value));
        }
    }
    return model;
}

float SolubilityRandomForestModel::predict(const std::vector<float>& features) const {
    float sum = 0.0f;
    for (const auto& tree : trees) {
        int node_id = 0;
        while (true) {
            const SolubilityTreeNode& node = tree[node_id];
            if (node.left_child == -1) {
                sum += node.value;
                break;
            }
            if (node.feature >= 0 && static_cast<size_t>(node.feature) < features.size()) {
                if (features[node.feature] <= node.threshold) {
                    node_id = node.left_child;
                } else {
                    node_id = node.right_child;
                }
            } else {
                break;
            }
        }
    }
    return (n_estimators > 0) ? (sum / n_estimators) : 0.0f;
}

static std::vector<float> createSolubilityTfidfFeatures(const std::string& smiles, const SolubilityRandomForestModel& model) {
    std::vector<float> features(model.n_features, 0.0f);
    if (model.analyzer_type == 0) {
        for (int n = model.ngram_min; n <= model.ngram_max; n++) {
            for (size_t i = 0; i + n <= smiles.length(); i++) {
                std::string ngram = smiles.substr(i, n);
                auto it = model.vocabulary.find(ngram);
                if (it != model.vocabulary.end()) {
                    int idx = it->second;
                    if (idx < model.n_features) features[idx] += 1.0f;
                }
            }
        }
    }
    for (int i = 0; i < model.n_features; i++) {
        if (features[i] > 0 && i < static_cast<int>(model.idf.size())) {
            features[i] *= model.idf[i];
        }
    }
    return features;
}

const SolubilityRandomForestModel* SolubilityDescriptor::ModelCache::getModel(const std::string& modelPath) {
    auto it = models.find(modelPath);
    if (it != models.end()) return it->second.get();
    try {
        auto model = std::make_unique<SolubilityRandomForestModel>(loadSolubilityModel(modelPath));
        SolubilityRandomForestModel* modelPtr = model.get();
        models[modelPath] = std::move(model);
        globalLogger.info("Loaded solubility model: " + modelPath);
        return modelPtr;
    } catch (const std::exception& e) {
        globalLogger.error("Error loading solubility model: " + std::string(e.what()));
        return nullptr;
    } catch (...) {
        globalLogger.error("Unknown error loading solubility model: " + modelPath);
        return nullptr;
    }
}

SolubilityDescriptor::SolubilityDescriptor()
    : Descriptor("Solubility", "Predicted aqueous solubility (logS, mol/L) using ML model") {}

std::variant<double, int, std::string> SolubilityDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return std::numeric_limits<double>::quiet_NaN();
    const RDKit::ROMol* rdmol = mol.getMolecule().get();
    if (!rdmol) return std::numeric_limits<double>::quiet_NaN();
    try {
        std::string smiles = RDKit::MolToSmiles(*rdmol);
        const SolubilityRandomForestModel* model = s_modelCache.getModel("models/solubilities.bin");
        if (!model) return std::numeric_limits<double>::quiet_NaN();
        std::vector<float> features = createSolubilityTfidfFeatures(smiles, *model);
        double prediction = model->predict(features);
        return prediction;
    } catch (const std::exception& e) {
        globalLogger.error("Solubility prediction error: " + std::string(e.what()));
        return std::numeric_limits<double>::quiet_NaN();
    } catch (...) {
        globalLogger.error("Unknown solubility prediction error");
        return std::numeric_limits<double>::quiet_NaN();
    }
}

void SolubilityDescriptor::clearThreadLocalCaches() {
    s_modelCache.clear();
}

} // namespace descriptors
} // namespace desfact
