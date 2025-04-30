#include "descriptors/pka.hpp"
#include <GraphMol/GraphMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <limits> // For std::numeric_limits
#include <vector>
#include <fstream>
#include <unordered_map>
#include <cmath>
#include <memory>
#include "utils.hpp" // For globalLogger
#include <GraphMol/SmilesParse/SmilesWrite.h>

namespace desfact {
namespace descriptors {

// Define the thread-local cache instance
thread_local PkaDescriptor::PatternCache PkaDescriptor::s_patternCache;
// Add model cache storage
thread_local PkaDescriptor::ModelCache PkaDescriptor::s_modelCache;

// Load model implementation
RandomForestModel loadModel(const std::string& filename) {
    std::ifstream file(filename, std::ios::binary);
    if (!file) {
        throw std::runtime_error("Cannot open model file: " + filename);
    }
    
    RandomForestModel model;
    
    // Read header
    char magic[6] = {0};
    file.read(magic, 5);
    if (std::string(magic) != "RFBIN") {
        throw std::runtime_error("Invalid model file format");
    }
    
    uint16_t version;
    file.read(reinterpret_cast<char*>(&version), sizeof(version));
    if (version != 1) {
        throw std::runtime_error("Unsupported model version");
    }
    
    file.read(reinterpret_cast<char*>(&model.n_estimators), sizeof(model.n_estimators));
    file.read(reinterpret_cast<char*>(&model.n_features), sizeof(model.n_features));
    
    // Read vocabulary
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
    
    // Read IDF values
    uint32_t idf_size;
    file.read(reinterpret_cast<char*>(&idf_size), sizeof(idf_size));
    model.idf.resize(idf_size);
    
    for (uint32_t i = 0; i < idf_size; i++) {
        file.read(reinterpret_cast<char*>(&model.idf[i]), sizeof(float));
    }
    
    // Read vectorizer parameters
    uint8_t analyzer_type_byte; 
    file.read(reinterpret_cast<char*>(&analyzer_type_byte), sizeof(analyzer_type_byte));
    model.analyzer_type = static_cast<int>(analyzer_type_byte); // Convert to int for the struct

    file.read(reinterpret_cast<char*>(&model.ngram_min), sizeof(model.ngram_min));
    file.read(reinterpret_cast<char*>(&model.ngram_max), sizeof(model.ngram_max));
    
    // Read trees
    model.trees.resize(model.n_estimators);
    
    for (int tree_idx = 0; tree_idx < model.n_estimators; tree_idx++) {
        uint32_t node_count;
        file.read(reinterpret_cast<char*>(&node_count), sizeof(node_count));
        
        model.trees[tree_idx].resize(node_count);
        
        for (uint32_t node_idx = 0; node_idx < node_count; node_idx++) {
            TreeNode& node = model.trees[tree_idx][node_idx];
            file.read(reinterpret_cast<char*>(&node.feature), sizeof(node.feature));
            file.read(reinterpret_cast<char*>(&node.threshold), sizeof(node.threshold));
            file.read(reinterpret_cast<char*>(&node.left_child), sizeof(node.left_child));
            file.read(reinterpret_cast<char*>(&node.right_child), sizeof(node.right_child));
            file.read(reinterpret_cast<char*>(&node.value), sizeof(node.value));
        }
    }
    
    return model;
}

// Create TF-IDF features from SMILES
std::vector<float> createTfidfFeatures(const std::string& smiles, const RandomForestModel& model) {
    std::vector<float> features(model.n_features, 0.0f);
    
    if (model.analyzer_type == 0) { // char analyzer
        // Extract character n-grams
        for (int n = model.ngram_min; n <= model.ngram_max; n++) {
            for (size_t i = 0; i <= smiles.length() - n; i++) {
                std::string ngram = smiles.substr(i, n);
                auto it = model.vocabulary.find(ngram);
                if (it != model.vocabulary.end()) {
                    int idx = it->second;
                    if (idx < model.n_features) {
                        features[idx] += 1.0f; // TF (simple count for now)
                    }
                }
            }
        }
    } else { // word analyzer
        // Implement if needed - for SMILES, char analyzer is typically more appropriate
        globalLogger.warning("Word analyzer not implemented for SMILES TF-IDF");
    }
    
    // Apply IDF weights
    for (int i = 0; i < model.n_features; i++) {
        if (features[i] > 0 && i < static_cast<int>(model.idf.size())) {
            features[i] *= model.idf[i]; // TF-IDF = TF * IDF
        }
    }
    
    return features;
}

// Get or load model from cache
const RandomForestModel* PkaDescriptor::ModelCache::getModel(const std::string& modelPath) {
    auto it = models.find(modelPath);
    if (it != models.end()) {
        return it->second.get();
    }

    try {
        auto model = std::make_unique<RandomForestModel>(loadModel(modelPath));
        RandomForestModel* modelPtr = model.get();
        models[modelPath] = std::move(model);
        globalLogger.info("Successfully loaded model from: " + modelPath);
        return modelPtr;
    } catch (const std::exception& e) {
        globalLogger.error("Error loading model from " + modelPath + ": " + std::string(e.what()));
        return nullptr;
    } catch (...) {
        globalLogger.error("Unknown error loading model from: " + modelPath);
        return nullptr;
    }
}

// Get or create SMARTS pattern, managing ownership
RDKit::ROMol* PkaDescriptor::PatternCache::getPattern(const std::string& smarts) {
    auto it = patterns.find(smarts);
    if (it != patterns.end()) {
        return it->second.get();
    }

    RDKit::ROMol* patternMol = nullptr;
    try {
        patternMol = RDKit::SmartsToMol(smarts);
        if (!patternMol) {
            globalLogger.warning("Failed to parse SMARTS pattern (nullptr): " + smarts);
            return nullptr;
        }
        patterns[smarts] = std::unique_ptr<RDKit::ROMol>(patternMol);
        return patternMol;
    } catch (const std::exception& e) {
        globalLogger.warning("SMARTS parsing exception for '" + smarts + "': " + std::string(e.what()));
        if (patternMol) delete patternMol; // Clean up if allocation happened before exception
        return nullptr;
    } catch (...) {
         globalLogger.warning("Unknown SMARTS parsing exception for '" + smarts + "'");
         if (patternMol) delete patternMol;
         return nullptr;
    }
}

// Base PkaDescriptor constructor
PkaDescriptor::PkaDescriptor(const std::string& name, const std::string& description)
    : Descriptor(name, description) {}

// Helper: Classify molecule type based on SMARTS patterns
std::string PkaDescriptor::classifyMoleculeType(const RDKit::ROMol& mol) {
    struct PatternInfo {
        std::string name;
        std::string smarts;
    };

    // Acid patterns
    const std::vector<PatternInfo> acidPatterns = {
        {"carboxylic_acid", "[CX3](=O)[OX2H1]"},
        {"phenol", "[c][OX2H1]"},
        {"sulfonic_acid", "[SX4](=O)(=O)[OX2H1]"},
        {"sulfinic_acid", "[S+2](=O)[O-]"}, // Check SMARTS validity
        {"phosphonic_acid", "[PX4](=O)([OX2H1])[OX2H1]"},
        {"phosphinic_acid", "[PX3](=O)[OX2H1]"},
        {"arsonic_acid", "[AsX4](=O)([OX2H1])[OX2H1]"}, // As might require valence check
        {"tetrazole", "c1[n][n][n][nH]1"},
        {"sulfonamide", "[NX3;H1,H2][S](=O)=O"},
        {"imide", "[NX3;H1](C=O)(C=O)"},
        {"thiol", "[SX2H1]"},
        {"alphaCarbonyl", "[CX4H1,CX4H2][CX3](=O)"},
        {"enol", "[CX3]=[CX3][OX2H1]"},
        {"boronic_acid", "[BX3]([OX2H1])[OX2H1]"}
    };

    // Base patterns
    const std::vector<PatternInfo> basePatterns = {
        {"primaryAmine", "[NX3;H2;!$(NC=[O,S,N]);!$(NS=O)]"},
        {"secondaryAmine", "[NX3;H1;!$(NC=[O,S,N]);!$(NS=O)]"},
        {"tertiaryAmine", "[NX3;H0;!$(NC=[O,S,N]);!$(NS=O)]"},
        {"aniline", "[c][NX3;H1,H2]"},
        {"pyridine", "[n;$(n:c:c)]"},
        {"pyridineLike", "[nX2r6]"},
        {"imidazoleLike", "[nX2r5]"},
        {"amidineGuanidine", "[NX3][CX3](=[NX2])[#6,#7]"}, // Check validity
        {"imine", "[CX3]=[NX2]"}
    };

    bool isAcid = false;
    for (const auto& p_info : acidPatterns) {
        RDKit::ROMol* pattern = s_patternCache.getPattern(p_info.smarts);
        if (pattern && RDKit::SubstructMatch(mol, *pattern).size() > 0) {
            isAcid = true;
            break;
        }
    }

    bool isBase = false;
    for (const auto& p_info : basePatterns) {
        RDKit::ROMol* pattern = s_patternCache.getPattern(p_info.smarts);
        if (pattern && RDKit::SubstructMatch(mol, *pattern).size() > 0) {
            isBase = true;
            break;
        }
    }

    if (isAcid && isBase) return "zwitterion";
    if (isAcid) return "acid";
    if (isBase) return "base";
    return "neutral";
}

// Helper: Estimate pKa for acidic groups using ML model or rule-based fallback
double PkaDescriptor::estimateAcidPKa(const RDKit::ROMol& mol) {
    // Try ML model first
    std::string modelPath = "models/acids.bin";
    globalLogger.info("Attempting to load acid model from: " + modelPath);
    
    const RandomForestModel* acidModel = s_modelCache.getModel(modelPath);
    if (acidModel) {
        globalLogger.info("Successfully loaded acid model!");
        try {
            std::string smarts = RDKit::MolToSmarts(mol);
            globalLogger.info("Processing acid SMARTS: " + smarts);
            std::vector<float> features = createTfidfFeatures(smarts, *acidModel);
            double prediction = acidModel->predict(features);
            globalLogger.info("ML acid prediction: " + std::to_string(prediction));
            return prediction;
        } catch (const std::exception& e) {
            globalLogger.warning("Error using acid ML model: " + std::string(e.what()) + ". Falling back to rule-based.");
        }
    } else {
        globalLogger.info("Acid model not found - using rule-based fallback method");
    }
    
    // Fallback to rule-based method
    const std::vector<std::pair<std::string, double>> acidRules = {
        {"[CX3](=O)[OX2H1]", 4.0}, // Carboxylic acid
        {"[c][OX2H1]", 10.0},      // Phenol
        {"[SX4](=O)(=O)[OX2H1]", 1.0}, // Sulfonic acid (approx)
        {"c1[n][n][n][nH]1", 5.0},    // Tetrazole (approx)
        {"[NX3;H1,H2][S](=O)=O", 9.0}, // Sulfonamide (approx)
        {"[NX3;H1](C=O)(C=O)", 9.5},    // Imide (approx)
        {"[SX2H1]", 8.5},             // Thiol
        {"[CX4H1,CX4H2][CX3](=O)", 19.0} // Alpha-carbonyl C-H (very weak acid)
    };

    double lowestPKa = std::numeric_limits<double>::infinity();

    for (const auto& rule : acidRules) {
        RDKit::ROMol* pattern = s_patternCache.getPattern(rule.first);
        if (pattern && RDKit::SubstructMatch(mol, *pattern).size() > 0) {
            lowestPKa = std::min(lowestPKa, rule.second);
        }
    }

    double finalPka = (lowestPKa == std::numeric_limits<double>::infinity()) ? std::numeric_limits<double>::quiet_NaN() : lowestPKa;
    globalLogger.info("Rule-based acid pKa: " + std::to_string(finalPka));
    return finalPka;
}

// Helper: Estimate pKa for basic groups using ML model or rule-based fallback
double PkaDescriptor::estimateBasePKa(const RDKit::ROMol& mol) {
    // Try ML model first
    const RandomForestModel* baseModel = s_modelCache.getModel("models/bases.bin");
    if (baseModel) {
        try {
            std::string smarts = RDKit::MolToSmarts(mol);
            std::vector<float> features = createTfidfFeatures(smarts, *baseModel);
            return baseModel->predict(features);
        } catch (const std::exception& e) {
            globalLogger.warning("Error using base ML model: " + std::string(e.what()) + ". Falling back to rule-based.");
        }
    }
    
    // Fallback to rule-based method
    const std::vector<std::pair<std::string, double>> baseRulesPKb = {
        {"[NX3;H2;!$(NC=[O,S,N]);!$(NS=O)]", 3.5},  // Primary aliphatic amine (pKb ~3.5 -> pKa ~10.5)
        {"[NX3;H1;!$(NC=[O,S,N]);!$(NS=O)]", 3.0},  // Secondary aliphatic amine (pKb ~3.0 -> pKa ~11.0)
        {"[NX3;H0;!$(NC=[O,S,N]);!$(NS=O)]", 4.0},  // Tertiary aliphatic amine (pKb ~4.0 -> pKa ~10.0)
        {"[c][NX3;H1,H2]", 9.5},                   // Aniline (pKb ~9.5 -> pKa ~4.5)
        {"[n;$(n:c:c)]", 8.8},                     // Pyridine (pKb ~8.8 -> pKa ~5.2)
        {"[nX2r5]", 7.0},                          // Imidazole-like (pKb ~7.0 -> pKa ~7.0)
        {"[NX3][CX3](=[NX2])[#6,#7]", 0.5}        // Guanidine/Amidine (very strong base, pKb ~0.5 -> pKa ~13.5)
    };

    double highestPKa = -std::numeric_limits<double>::infinity();

    for (const auto& rule : baseRulesPKb) {
        RDKit::ROMol* pattern = s_patternCache.getPattern(rule.first);
        if (pattern && RDKit::SubstructMatch(mol, *pattern).size() > 0) {
            double pKa_conj_acid = 14.0 - rule.second;
            highestPKa = std::max(highestPKa, pKa_conj_acid);
        }
    }

    return (highestPKa == -std::numeric_limits<double>::infinity()) ? std::numeric_limits<double>::quiet_NaN() : highestPKa;
}

// --- AcidityTypeDescriptor ---
AcidityTypeDescriptor::AcidityTypeDescriptor()
    : PkaDescriptor("AcidityType", "Classifies molecule as acid, base, zwitterion, or neutral based on SMARTS") {}

std::variant<double, int, std::string> AcidityTypeDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) {
        return std::string("Error: Invalid Molecule");
    }
    const RDKit::ROMol* rdmol = mol.getMolecule().get();
    if (!rdmol) return std::string("Error: Internal Null");

    try {
        return classifyMoleculeType(*rdmol);
    } catch (const std::exception& e) {
         globalLogger.error("Error classifying molecule: " + std::string(e.what()));
         return std::string("Error: Classification");
    } catch (...) {
         globalLogger.error("Unknown error during molecule classification.");
         return std::string("Error: Unknown Classification");
    }
}

// --- PkaEstimate1Descriptor ---
PkaEstimate1Descriptor::PkaEstimate1Descriptor()
    : PkaDescriptor("pKa1", "Estimated most acidic pKa (or basic pKa if only basic groups present)") {}

std::variant<double, int, std::string> PkaEstimate1Descriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    const RDKit::ROMol* rdmol = mol.getMolecule().get();
    if (!rdmol) return std::numeric_limits<double>::quiet_NaN();

    try {
        std::string type = classifyMoleculeType(*rdmol);
        double pka = std::numeric_limits<double>::quiet_NaN();

        if (type == "acid" || type == "zwitterion") {
            pka = estimateAcidPKa(*rdmol);
        } else if (type == "base") {
            pka = estimateBasePKa(*rdmol);
        }

        return pka;

    } catch (const std::exception& e) {
         globalLogger.error("Error estimating pKa1: " + std::string(e.what()));
         return std::string("Error: pKa1 Calc");
    } catch (...) {
         globalLogger.error("Unknown error during pKa1 estimation.");
         return std::string("Error: Unknown pKa1");
    }
}

// --- PkaEstimate2Descriptor ---
PkaEstimate2Descriptor::PkaEstimate2Descriptor()
    : PkaDescriptor("pKa2", "Estimated most basic pKa (only for zwitterions, otherwise NaN)") {}

std::variant<double, int, std::string> PkaEstimate2Descriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    const RDKit::ROMol* rdmol = mol.getMolecule().get();
    if (!rdmol) return std::numeric_limits<double>::quiet_NaN();

    try {
        std::string type = classifyMoleculeType(*rdmol);
        double pka = std::numeric_limits<double>::quiet_NaN();

        if (type == "zwitterion") {
            pka = estimateBasePKa(*rdmol);
        }

        return pka;

    } catch (const std::exception& e) {
         globalLogger.error("Error estimating pKa2: " + std::string(e.what()));
         return std::string("Error: pKa2 Calc");
    } catch (...) {
         globalLogger.error("Unknown error during pKa2 estimation.");
         return std::string("Error: Unknown pKa2");
    }
}

// --- IsAcidDescriptor ---
IsAcidDescriptor::IsAcidDescriptor()
    : PkaDescriptor("isAcid", "Returns 1 if molecule is classified as acidic, 0 otherwise") {}

std::variant<double, int, std::string> IsAcidDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) {
        return 0;
    }
    const RDKit::ROMol* rdmol = mol.getMolecule().get();
    if (!rdmol) return 0;

    try {
        std::string type = classifyMoleculeType(*rdmol);
        return (type == "acid" || type == "zwitterion") ? 1 : 0;
    } catch (...) {
         globalLogger.error("Unknown error during isAcid classification.");
        return 0;
    }
}

// --- IsBaseDescriptor ---
IsBaseDescriptor::IsBaseDescriptor()
    : PkaDescriptor("isBase", "Returns 1 if molecule is classified as basic, 0 otherwise") {}

std::variant<double, int, std::string> IsBaseDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) {
        return 0;
    }
    const RDKit::ROMol* rdmol = mol.getMolecule().get();
    if (!rdmol) return 0;

    try {
        std::string type = classifyMoleculeType(*rdmol);
        return (type == "base" || type == "zwitterion") ? 1 : 0;
    } catch (...) {
         globalLogger.error("Unknown error during isBase classification.");
        return 0;
    }
}

// --- IsNeutralDescriptor ---
IsNeutralDescriptor::IsNeutralDescriptor()
    : PkaDescriptor("isNeutral", "Returns 1 if molecule is classified as neutral, 0 otherwise") {}

std::variant<double, int, std::string> IsNeutralDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) {
        return 0;
    }
    const RDKit::ROMol* rdmol = mol.getMolecule().get();
    if (!rdmol) return 0;

    try {
        std::string type = classifyMoleculeType(*rdmol);
        return (type == "neutral") ? 1 : 0;
    } catch (...) {
         globalLogger.error("Unknown error during isNeutral classification.");
        return 0;
    }
}

// --- IsZwitterionDescriptor ---
IsZwitterionDescriptor::IsZwitterionDescriptor()
    : PkaDescriptor("isZwitterion", "Returns 1 if molecule is classified as zwitterionic, 0 otherwise") {}

std::variant<double, int, std::string> IsZwitterionDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) {
        return 0;
    }
    const RDKit::ROMol* rdmol = mol.getMolecule().get();
    if (!rdmol) return 0;

    try {
        std::string type = classifyMoleculeType(*rdmol);
        return (type == "zwitterion") ? 1 : 0;
    } catch (...) {
         globalLogger.error("Unknown error during isZwitterion classification.");
        return 0;
    }
}

// Implement clearThreadLocalCaches method
void PkaDescriptor::clearThreadLocalCaches() {
    s_patternCache.patterns.clear();
    s_modelCache.clear();
}

} // namespace descriptors
} // namespace desfact