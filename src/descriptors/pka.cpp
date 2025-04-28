#include "descriptors/pka.hpp"
#include <GraphMol/GraphMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <limits> // For std::numeric_limits
#include <vector>
#include "utils.hpp" // For globalLogger

namespace desfact {
namespace descriptors {

// Define the thread-local cache instance
thread_local PkaDescriptor::PatternCache PkaDescriptor::s_patternCache;

// Destructor for PatternCache to clean up ROMol pointers
PkaDescriptor::PatternCache::~PatternCache() {
    // unique_ptr handles deletion automatically
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

// Helper: Estimate pKa for acidic groups (simplified rule-based)
double PkaDescriptor::estimateAcidPKa(const RDKit::ROMol& mol) {
     // Define acidic SMARTS patterns with approximate pKa values
    const std::vector<std::pair<std::string, double>> acidRules = {
        {"[CX3](=O)[OX2H1]", 4.0}, // Carboxylic acid
        {"[c][OX2H1]", 10.0},      // Phenol
        {"[SX4](=O)(=O)[OX2H1]", 1.0}, // Sulfonic acid (approx)
        {"c1[n][n][n][nH]1", 5.0},    // Tetrazole (approx)
        {"[NX3;H1,H2][S](=O)=O", 9.0}, // Sulfonamide (approx)
        {"[NX3;H1](C=O)(C=O)", 9.5},    // Imide (approx)
        {"[SX2H1]", 8.5},             // Thiol
         {"[CX4H1,CX4H2][CX3](=O)", 19.0} // Alpha-carbonyl C-H (very weak acid)
        // Add more rules as needed
    };

    double lowestPKa = std::numeric_limits<double>::infinity(); // Initialize high

    for (const auto& rule : acidRules) {
        RDKit::ROMol* pattern = s_patternCache.getPattern(rule.first);
        if (pattern && RDKit::SubstructMatch(mol, *pattern).size() > 0) {
            lowestPKa = std::min(lowestPKa, rule.second);
        }
    }

    // Return lowest pKa found, or NaN if no acidic group matched
     return (lowestPKa == std::numeric_limits<double>::infinity()) ? std::numeric_limits<double>::quiet_NaN() : lowestPKa;
}

// Helper: Estimate pKa for basic groups (simplified rule-based)
double PkaDescriptor::estimateBasePKa(const RDKit::ROMol& mol) {
     // Define basic SMARTS patterns with approximate pKb values (convert to pKa later)
     // pKa = 14 - pKb (approx for conjugate acid)
    const std::vector<std::pair<std::string, double>> baseRulesPKb = {
         {"[NX3;H2;!$(NC=[O,S,N]);!$(NS=O)]", 3.5},  // Primary aliphatic amine (pKb ~3.5 -> pKa ~10.5)
         {"[NX3;H1;!$(NC=[O,S,N]);!$(NS=O)]", 3.0},  // Secondary aliphatic amine (pKb ~3.0 -> pKa ~11.0)
         {"[NX3;H0;!$(NC=[O,S,N]);!$(NS=O)]", 4.0},  // Tertiary aliphatic amine (pKb ~4.0 -> pKa ~10.0)
         {"[c][NX3;H1,H2]", 9.5},                   // Aniline (pKb ~9.5 -> pKa ~4.5)
         {"[n;$(n:c:c)]", 8.8},                     // Pyridine (pKb ~8.8 -> pKa ~5.2)
         {"[nX2r5]", 7.0},                          // Imidazole-like (pKb ~7.0 -> pKa ~7.0)
         {"[NX3][CX3](=[NX2])[#6,#7]", 0.5}        // Guanidine/Amidine (very strong base, pKb ~0.5 -> pKa ~13.5)
         // Add more rules as needed
    };

    double highestPKa = -std::numeric_limits<double>::infinity(); // Initialize low

    for (const auto& rule : baseRulesPKb) {
        RDKit::ROMol* pattern = s_patternCache.getPattern(rule.first);
        if (pattern && RDKit::SubstructMatch(mol, *pattern).size() > 0) {
             double pKa_conj_acid = 14.0 - rule.second; // Approximate conversion
            highestPKa = std::max(highestPKa, pKa_conj_acid);
        }
    }

     // Return highest pKa found, or NaN if no basic group matched
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
        return std::numeric_limits<double>::quiet_NaN(); // Use NaN for double errors
    }
    const RDKit::ROMol* rdmol = mol.getMolecule().get();
    if (!rdmol) return std::numeric_limits<double>::quiet_NaN();

    try {
        std::string type = classifyMoleculeType(*rdmol);
        double pka = std::numeric_limits<double>::quiet_NaN();

        if (type == "acid" || type == "zwitterion") {
            pka = estimateAcidPKa(*rdmol);
        } else if (type == "base") {
            // For bases, pKa1 conventionally represents the pKa of the conjugate acid
            pka = estimateBasePKa(*rdmol);
        }
        // For "neutral", pka remains NaN

        return pka; // Returns NaN if no relevant group found or type is neutral

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
        return std::numeric_limits<double>::quiet_NaN(); // Use NaN for double errors
    }
    const RDKit::ROMol* rdmol = mol.getMolecule().get();
    if (!rdmol) return std::numeric_limits<double>::quiet_NaN();

    try {
        std::string type = classifyMoleculeType(*rdmol);
        double pka = std::numeric_limits<double>::quiet_NaN();

        if (type == "zwitterion") {
            // For zwitterions, pKa2 represents the pKa of the basic group's conjugate acid
            pka = estimateBasePKa(*rdmol);
        }
        // For "acid", "base", "neutral", pKa2 is NaN

        return pka; // Returns NaN if not zwitterion or no basic group found

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
        return 0; // Return 0 for invalid molecules
    }
    const RDKit::ROMol* rdmol = mol.getMolecule().get();
    if (!rdmol) return 0;

    try {
        std::string type = classifyMoleculeType(*rdmol);
        return (type == "acid" || type == "zwitterion") ? 1 : 0;
    } catch (...) {
         // Log error? For a boolean descriptor, returning 0 might be safest on error.
         globalLogger.error("Unknown error during isAcid classification.");
        return 0;
    }
}

// --- IsBaseDescriptor ---
IsBaseDescriptor::IsBaseDescriptor()
    : PkaDescriptor("isBase", "Returns 1 if molecule is classified as basic, 0 otherwise") {}

std::variant<double, int, std::string> IsBaseDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) {
        return 0; // Return 0 for invalid molecules
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
        return 0; // Return 0 for invalid molecules (could argue otherwise)
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
        return 0; // Return 0 for invalid molecules
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

} // namespace descriptors
} // namespace desfact