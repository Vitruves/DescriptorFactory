// include/descriptors/pka.hpp
#pragma once

#include "descriptors.hpp"
#include <GraphMol/GraphMol.h>
#include <GraphMol/Atom.h>
#include <string>
#include <variant>
#include <unordered_map>

// Forward declare RDKit types
namespace RDKit {
    class ROMol;
}

namespace desfact {
namespace descriptors {

// Common utility functions for PKA descriptors
namespace pka_utils {
    std::string classifyMolecule(const Molecule& mol);
    double predictPKAFromModel(const Molecule& mol, const std::string& modelPath);
}

// Base class for pKa related descriptors (optional, but can group helpers)
class PkaDescriptor : public Descriptor {
protected:
    // Helper functions (to be defined in pka.cpp)
    static std::string classifyMoleculeType(const RDKit::ROMol& mol);
    static double estimateAcidPKa(const RDKit::ROMol& mol);
    static double estimateBasePKa(const RDKit::ROMol& mol);
    
    // Cache for SMARTS patterns (potentially static thread_local in cpp)
    struct PatternCache {
         std::unordered_map<std::string, std::unique_ptr<RDKit::ROMol>> patterns;
         ~PatternCache(); // Destructor to clean up ROMol pointers
         RDKit::ROMol* getPattern(const std::string& smarts);
    };
    static thread_local PatternCache s_patternCache;


public:
    PkaDescriptor(const std::string& name, const std::string& description);
    virtual ~PkaDescriptor() = default;
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