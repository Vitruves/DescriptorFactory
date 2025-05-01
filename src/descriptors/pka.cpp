#include "descriptors/pka.hpp"
#include <GraphMol/GraphMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <limits>
#include <vector>
#include <fstream>
#include <unordered_map>
#include <cmath>
#include <memory>
#include <mutex>
#include "utils.hpp"
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <filesystem>

namespace desfact {
namespace descriptors {

// Define the thread-local cache instance
thread_local PkaDescriptor::PatternCache PkaDescriptor::s_patternCache;
// Add model cache storage with mutex protection
thread_local PkaDescriptor::ModelCache PkaDescriptor::s_modelCache;
std::mutex globalModelMutex;

// Check if model file exists to avoid crashes when trying to load missing files
bool modelFileExists(const std::string& filename) {
    std::filesystem::path modelPath(filename);
    return std::filesystem::exists(modelPath) && std::filesystem::is_regular_file(modelPath);
}

// Load model implementation with better error handling
RandomForestModel loadModel(const std::string& filename) {
    if (!modelFileExists(filename)) {
        throw std::runtime_error("Model file does not exist: " + filename);
    }
    
    std::ifstream file(filename, std::ios::binary);
    if (!file) {
        throw std::runtime_error("Cannot open model file: " + filename);
    }
    
    RandomForestModel model;
    
    try {
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
        
        // Validate estimators and features count to prevent crashes
        if (model.n_estimators <= 0 || model.n_estimators > 10000) {
            throw std::runtime_error("Invalid number of estimators: " + std::to_string(model.n_estimators));
        }
        
        if (model.n_features <= 0 || model.n_features > 1000000) {
            throw std::runtime_error("Invalid number of features: " + std::to_string(model.n_features));
        }
        
        // Read vocabulary
        uint32_t vocab_size;
        file.read(reinterpret_cast<char*>(&vocab_size), sizeof(vocab_size));
        
        // Validate vocabulary size
        if (vocab_size > 10000000) {
            throw std::runtime_error("Vocabulary size too large: " + std::to_string(vocab_size));
        }
        
        for (uint32_t i = 0; i < vocab_size; i++) {
            uint32_t word_len;
            file.read(reinterpret_cast<char*>(&word_len), sizeof(word_len));
            
            // Validate word length
            if (word_len > 1000) {
                throw std::runtime_error("Word length too large: " + std::to_string(word_len));
            }
            
            std::string word(word_len, ' ');
            file.read(&word[0], word_len);
            
            uint32_t idx;
            file.read(reinterpret_cast<char*>(&idx), sizeof(idx));
            
            // Validate index
            if (idx >= static_cast<uint32_t>(model.n_features)) {
                throw std::runtime_error("Feature index out of bounds: " + std::to_string(idx));
            }
            
            model.vocabulary[word] = idx;
        }
        
        // Read IDF values
        uint32_t idf_size;
        file.read(reinterpret_cast<char*>(&idf_size), sizeof(idf_size));
        
        // Validate IDF size
        if (idf_size > 1000000) {
            throw std::runtime_error("IDF size too large: " + std::to_string(idf_size));
        }
        
        model.idf.resize(idf_size);
        
        for (uint32_t i = 0; i < idf_size; i++) {
            file.read(reinterpret_cast<char*>(&model.idf[i]), sizeof(float));
        }
        
        // Read vectorizer parameters
        uint8_t analyzer_type_byte; 
        file.read(reinterpret_cast<char*>(&analyzer_type_byte), sizeof(analyzer_type_byte));
        model.analyzer_type = static_cast<int>(analyzer_type_byte);
    
        file.read(reinterpret_cast<char*>(&model.ngram_min), sizeof(model.ngram_min));
        file.read(reinterpret_cast<char*>(&model.ngram_max), sizeof(model.ngram_max));
        
        // Validate ngram range
        if (model.ngram_min < 1 || model.ngram_min > 10 || 
            model.ngram_max < model.ngram_min || model.ngram_max > 10) {
            throw std::runtime_error("Invalid ngram range: " + std::to_string(model.ngram_min) + 
                                   "-" + std::to_string(model.ngram_max));
        }
        
        // Read trees
        model.trees.resize(model.n_estimators);
        
        for (int tree_idx = 0; tree_idx < model.n_estimators; tree_idx++) {
            uint32_t node_count;
            file.read(reinterpret_cast<char*>(&node_count), sizeof(node_count));
            
            // Validate node count
            if (node_count > 1000000) {
                throw std::runtime_error("Tree node count too large: " + std::to_string(node_count));
            }
            
            model.trees[tree_idx].resize(node_count);
            
            for (uint32_t node_idx = 0; node_idx < node_count; node_idx++) {
                TreeNode& node = model.trees[tree_idx][node_idx];
                file.read(reinterpret_cast<char*>(&node.feature), sizeof(node.feature));
                file.read(reinterpret_cast<char*>(&node.threshold), sizeof(node.threshold));
                file.read(reinterpret_cast<char*>(&node.left_child), sizeof(node.left_child));
                file.read(reinterpret_cast<char*>(&node.right_child), sizeof(node.right_child));
                file.read(reinterpret_cast<char*>(&node.value), sizeof(node.value));
                
                // Validate node references
                if ((node.left_child != -1 && (node.left_child < 0 || node.left_child >= static_cast<int>(node_count))) ||
                    (node.right_child != -1 && (node.right_child < 0 || node.right_child >= static_cast<int>(node_count)))) {
                    throw std::runtime_error("Invalid tree node reference");
                }
                
                // Validate feature index
                if (node.left_child != -1 && (node.feature < 0 || node.feature >= model.n_features)) {
                    throw std::runtime_error("Invalid feature index in tree node");
                }
            }
        }

        // In the loadModel function, after the vocabulary section, add a maximum feature size check
        if (model.n_features > 1000000) {
            throw std::runtime_error("Feature dimension too large: " + std::to_string(model.n_features));
        }
    } catch (const std::exception& e) {
        throw std::runtime_error("Error reading model file: " + std::string(e.what()));
    }
    
    return model;
}

// Get or load model from cache with thread safety
const RandomForestModel* PkaDescriptor::ModelCache::getModel(const std::string& modelPath) {
    // Check if the model already exists in the cache
    {
        std::lock_guard<std::mutex> lock(globalModelMutex);
        auto it = models.find(modelPath);
        if (it != models.end()) {
            return it->second.get();
        }
    }
    
    // If we don't have a cached model, check if the file exists first
    if (!modelFileExists(modelPath)) {
        globalLogger.error("Model file does not exist: " + modelPath);
        return nullptr;
    }
    
    try {
        // Load the model outside the lock to minimize contention
        auto model = std::make_unique<RandomForestModel>(loadModel(modelPath));
        
        // Store the model in the cache with the lock
        std::lock_guard<std::mutex> lock(globalModelMutex);
        auto it = models.find(modelPath); // Check again in case another thread loaded it
        if (it != models.end()) {
            return it->second.get();
        }
        
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

// Modify all code that creates TF-IDF features to include size limits
std::vector<float> createTfidfFeatures(const std::string& smiles, const RandomForestModel& model) {
    // Add size limit for input string
    if (smiles.empty() || model.n_features <= 0 || smiles.length() > 20000) {
        return std::vector<float>(model.n_features, 0.0f); // Return zeroed vector
    }
    
    // Limit maximum features to avoid Eigen crashes
    if (model.n_features > 100000) {
        globalLogger.warning("Feature dimension too large: " + std::to_string(model.n_features));
        return std::vector<float>(model.n_features, 0.0f);
    }
    
    std::vector<float> features(model.n_features, 0.0f);
    
    if (model.analyzer_type == 0) { // char analyzer
        // Extract character n-grams
        for (int n = model.ngram_min; n <= model.ngram_max; n++) {
            if (n <= 0 || static_cast<size_t>(n) > smiles.length()) continue;
            
            for (size_t i = 0; i + n <= smiles.length(); i++) {
                std::string ngram = smiles.substr(i, n);
                auto it = model.vocabulary.find(ngram);
                if (it != model.vocabulary.end()) {
                    int idx = it->second;
                    if (idx >= 0 && idx < model.n_features) {
                        features[idx] += 1.0f;
                    }
                }
            }
        }
    }
    
    // Apply IDF weights with bounds checking
    for (int i = 0; i < model.n_features && i < static_cast<int>(model.idf.size()); i++) {
        if (features[i] > 0) {
            features[i] *= model.idf[i];
        }
    }
    
    return features;
}

// Safer RandomForestModel prediction
float RandomForestModel::predict(const std::vector<float>& features) const {
    // Validate inputs
    if (trees.empty() || n_estimators <= 0) {
        return 0.0f;
    }
    
    // Ensure feature vector has correct size
    if (features.empty() || static_cast<int>(features.size()) != n_features) {
        return 0.0f;
    }
    
    float sum = 0.0f;
    int valid_trees = 0;
    
    // Limit the number of trees used for very large feature spaces to avoid memory issues
    int max_trees_to_use = std::min(n_estimators, 500);
    
    for (int i = 0; i < max_trees_to_use && i < static_cast<int>(trees.size()); i++) {
        const auto& tree = trees[i];
        if (tree.empty()) continue;
        
        try {
            int node_id = 0;
            int max_nodes = static_cast<int>(tree.size());
            int steps = 0;
            const int max_steps = 1000; // Prevent infinite loops
            bool valid_prediction = false;
            
            while (node_id >= 0 && node_id < max_nodes && steps < max_steps) {
                const TreeNode& node = tree[node_id];
                
                // Leaf node - get value and exit
                if (node.left_child == -1) {
                    sum += node.value;
                    valid_prediction = true;
                    break;
                }
                
                // Check feature bounds
                if (node.feature < 0 || static_cast<size_t>(node.feature) >= features.size()) {
                    break; // Skip this tree
                }
                
                // Navigate tree safely
                if (features[node.feature] <= node.threshold) {
                    if (node.left_child < 0 || node.left_child >= max_nodes) {
                        break; // Invalid node reference
                    }
                    node_id = node.left_child;
                } else {
                    if (node.right_child < 0 || node.right_child >= max_nodes) {
                        break; // Invalid node reference
                    }
                    node_id = node.right_child;
                }
                
                steps++;
            }
            
            if (valid_prediction) {
                valid_trees++;
            }
        } catch (...) {
            // Skip any tree that causes an exception
            continue;
        }
    }
    
    // Avoid division by zero
    return valid_trees > 0 ? (sum / valid_trees) : 0.0f;
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

// Update the clearThreadLocalCaches method to use public interface
void PkaDescriptor::clearThreadLocalCaches() {
    s_patternCache.patterns.clear();
    s_modelCache.clear();
}

// Add to the existing code, more robust SMARTS conversion with max size limit
std::string safeMolToSmarts(const RDKit::ROMol& mol) {
    try {
        // Only attempt SMARTS conversion for reasonable molecule sizes
        if (mol.getNumAtoms() > 200) {
            globalLogger.warning("Molecule too large for SMARTS conversion: " + 
                                 std::to_string(mol.getNumAtoms()) + " atoms");
            return "";
        }
        
        std::string smarts = RDKit::MolToSmarts(mol);
        // Also check if the SMARTS string is extremely large
        if (smarts.length() > 10000) {
            globalLogger.warning("Generated SMARTS too large: " + 
                                std::to_string(smarts.length()) + " characters");
            return "";
        }
        return smarts;
    } catch (const std::exception& e) {
        globalLogger.error("Error converting molecule to SMARTS: " + std::string(e.what()));
        return "";
    } catch (...) {
        globalLogger.error("Unknown error converting molecule to SMARTS");
        return "";
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

    // Acid patterns with corrected SMARTS
    const std::vector<PatternInfo> acidPatterns = {
        {"carboxylic_acid", "[CX3](=O)[OX2H1]"},
        {"phenol", "[c][OX2H1]"},
        {"sulfonic_acid", "[SX4](=O)(=O)[OX2H1]"},
        {"sulfinic_acid", "[SX3](=O)[OX2H1]"}, // Fixed SMARTS
        {"phosphonic_acid", "[PX4](=O)([OX2H1])[OX2H1]"},
        {"phosphinic_acid", "[PX3](=O)[OX2H1]"},
        {"arsonic_acid", "[As](=O)([OX2H1])[OX2H1]"}, // Simplified SMARTS
        {"tetrazole", "c1[nH][n][n]n1"}, // More specific tetrazole
        {"sulfonamide", "[NX3;H1,H2][SX4](=O)(=O)"},
        {"imide", "[NX3;H1]([CX3]=O)[CX3]=O"},
        {"thiol", "[SX2H1]"},
        {"alphaCarbonyl", "[CX4;H1,H2][CX3](=O)"},
        {"enol", "[CX3]=[CX3][OX2H1]"},
        {"boronic_acid", "[BX3]([OX2H1])[OX2H1]"}
    };

    // Base patterns with refined SMARTS
    const std::vector<PatternInfo> basePatterns = {
        {"primaryAmine", "[NX3;H2;!$(NC=[O,S,N]);!$(NS=O)]"},
        {"secondaryAmine", "[NX3;H1;!$(NC=[O,S,N]);!$(NS=O)]"},
        {"tertiaryAmine", "[NX3;H0;!$(NC=[O,S,N]);!$(NS=O)]"},
        {"aniline", "[c][NX3;H1,H2]"},
        {"pyridine", "[nX2]1:c:c:c:c:c1"}, // More specific pyridine
        {"pyridineLike", "[nX2r6]"},
        {"imidazoleLike", "[nX2r5]"},
        {"amidineGuanidine", "[NX3][CX3]=[NX2]"}, // Simplified SMARTS
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
    if (!mol.getNumAtoms()) {
        return std::numeric_limits<double>::quiet_NaN(); // Empty molecule
    }
    
    // For very large molecules, skip ML model and use rule-based method
    if (mol.getNumAtoms() > 30) { // Even lower threshold from 50 to 30
        globalLogger.info("Large molecule detected, skipping ML model for acid pKa");
        // Fall through to rule-based approach
    } else {
        const RandomForestModel* acidModel = s_modelCache.getModel("models/acids.bin");
        if (acidModel) {
            try {
                std::string smarts = safeMolToSmarts(mol);
                if (smarts.empty()) {
                    globalLogger.warning("Empty SMARTS generated for molecule");
                    // Fall through to rule-based approach
                } else if (smarts.length() > 5000) { // Add max SMARTS length check
                    globalLogger.warning("SMARTS too long for ML model: " + std::to_string(smarts.length()));
                    // Fall through to rule-based approach
                } else {
                    std::vector<float> features = createTfidfFeatures(smarts, *acidModel);
                    if (!features.empty() && features.size() == static_cast<size_t>(acidModel->n_features)) {
                        double prediction = acidModel->predict(features);
                        return prediction;
                    }
                }
            } catch (const std::exception& e) {
                globalLogger.warning("Error using acid ML model: " + std::string(e.what()) + 
                                    ". Falling back to rule-based.");
            } catch (...) {
                globalLogger.warning("Unknown error in acid ML model. Falling back to rule-based.");
            }
        }
    }
    
    // Fallback to rule-based method
    const std::vector<std::pair<std::string, double>> acidRules = {
        {"[CX3](=O)[OX2H1]", 4.0}, // Carboxylic acid
        {"[c][OX2H1]", 10.0},      // Phenol
        {"[SX4](=O)(=O)[OX2H1]", 1.0}, // Sulfonic acid
        {"c1[nH][n][n]n1", 5.0},    // Tetrazole
        {"[NX3;H1,H2][SX4](=O)(=O)", 9.0}, // Sulfonamide
        {"[NX3;H1]([CX3]=O)[CX3]=O", 9.5},  // Imide
        {"[SX2H1]", 8.5},             // Thiol
        {"[CX4;H1,H2][CX3](=O)", 19.0} // Alpha-carbonyl C-H
    };

    double lowestPKa = std::numeric_limits<double>::infinity();

    for (const auto& rule : acidRules) {
        RDKit::ROMol* pattern = s_patternCache.getPattern(rule.first);
        if (pattern && RDKit::SubstructMatch(mol, *pattern).size() > 0) {
            lowestPKa = std::min(lowestPKa, rule.second);
        }
    }

    return (lowestPKa == std::numeric_limits<double>::infinity()) ? 
        std::numeric_limits<double>::quiet_NaN() : lowestPKa;
}

// Helper: Estimate pKa for basic groups using ML model or rule-based fallback
double PkaDescriptor::estimateBasePKa(const RDKit::ROMol& mol) {
    // For very large molecules, skip ML model and use rule-based method
    if (mol.getNumAtoms() > 30) { // Even lower threshold from 50 to 30
        globalLogger.info("Large molecule detected, skipping ML model for base pKa");
        // Fall through to rule-based approach
    } else {
        const RandomForestModel* baseModel = s_modelCache.getModel("models/bases.bin");
        if (baseModel) {
            try {
                std::string smarts = safeMolToSmarts(mol);
                if (!smarts.empty()) {
                    std::vector<float> features = createTfidfFeatures(smarts, *baseModel);
                    if (!features.empty()) {
                        double prediction = baseModel->predict(features);
                        return prediction;
                    }
                }
            } catch (const std::exception& e) {
                globalLogger.warning("Error using base ML model: " + std::string(e.what()) + 
                                    ". Falling back to rule-based.");
            } catch (...) {
                globalLogger.warning("Unknown error in base ML model. Falling back to rule-based.");
            }
        }
    }
    
    // Fallback to rule-based method
    const std::vector<std::pair<std::string, double>> baseRulesPKb = {
        {"[NX3;H2;!$(NC=[O,S,N]);!$(NS=O)]", 3.5},  // Primary aliphatic amine
        {"[NX3;H1;!$(NC=[O,S,N]);!$(NS=O)]", 3.0},  // Secondary aliphatic amine
        {"[NX3;H0;!$(NC=[O,S,N]);!$(NS=O)]", 4.0},  // Tertiary aliphatic amine
        {"[c][NX3;H1,H2]", 9.5},                    // Aniline
        {"[nX2]1:c:c:c:c:c1", 8.8},                 // Pyridine
        {"[nX2r5]", 7.0},                           // Imidazole-like
        {"[NX3][CX3]=[NX2]", 0.5}                   // Guanidine/Amidine
    };

    double highestPKa = -std::numeric_limits<double>::infinity();

    for (const auto& rule : baseRulesPKb) {
        RDKit::ROMol* pattern = s_patternCache.getPattern(rule.first);
        if (pattern && RDKit::SubstructMatch(mol, *pattern).size() > 0) {
            double pKa_conj_acid = 14.0 - rule.second;
            highestPKa = std::max(highestPKa, pKa_conj_acid);
        }
    }

    return (highestPKa == -std::numeric_limits<double>::infinity()) ? 
        std::numeric_limits<double>::quiet_NaN() : highestPKa;
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
         return std::numeric_limits<double>::quiet_NaN();
    } catch (...) {
         globalLogger.error("Unknown error during pKa1 estimation.");
         return std::numeric_limits<double>::quiet_NaN();
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
         return std::numeric_limits<double>::quiet_NaN();
    } catch (...) {
         globalLogger.error("Unknown error during pKa2 estimation.");
         return std::numeric_limits<double>::quiet_NaN();
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
        return 0;
    }
}

} // namespace descriptors
} // namespace desfact