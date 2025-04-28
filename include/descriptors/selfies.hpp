#pragma once

#include "descriptors.hpp"
#include <vector>
#include <string>

namespace desfact {
namespace descriptors {

// Base class for SELFIES-based descriptors
class SelfiesDescriptor : public Descriptor {
protected:
    static std::vector<std::string> tokenizeSelfies(const std::string& selfiesStr);
    static bool isBranchToken(const std::string& token);
    static bool isRingToken(const std::string& token);
    static bool isAtomToken(const std::string& token);
    static bool isAromaticToken(const std::string& token);
    static bool isChargeToken(const std::string& token);
    static bool isHeteroatomToken(const std::string& token);
    static bool isStereoToken(const std::string& token);
    static bool isUnpairedElectronToken(const std::string& token);
    static int getRingSizeFromToken(const std::string& token);
    static double calculateTokenEntropy(const std::vector<std::string>& tokens);
    static int calculateBranchingDepth(const std::vector<std::string>& tokens);
    static double calculateSequenceSymmetry(const std::vector<std::string>& tokens);
    static int calculateLevenshteinDistance(const std::string& s1, const std::string& s2);

public:
    SelfiesDescriptor(const std::string& name, const std::string& description);
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
    virtual std::variant<double, int, std::string> calculateFromSelfies(const std::string& selfiesStr) const = 0;
};

class SelfiesTokenCountDescriptor : public SelfiesDescriptor {
public:
    SelfiesTokenCountDescriptor();
    std::variant<double, int, std::string> calculateFromSelfies(const std::string& selfiesStr) const override;
};

class SelfiesTokenDistributionDescriptor : public SelfiesDescriptor {
public:
    SelfiesTokenDistributionDescriptor();
    std::variant<double, int, std::string> calculateFromSelfies(const std::string& selfiesStr) const override;
};

class SelfiesAverageComplexityDescriptor : public SelfiesDescriptor {
public:
    SelfiesAverageComplexityDescriptor();
    std::variant<double, int, std::string> calculateFromSelfies(const std::string& selfiesStr) const override;
};

class SelfiesBranchingDepthDescriptor : public SelfiesDescriptor {
public:
    SelfiesBranchingDepthDescriptor();
    std::variant<double, int, std::string> calculateFromSelfies(const std::string& selfiesStr) const override;
};

class SelfiesMaxBranchingFanoutDescriptor : public SelfiesDescriptor {
public:
    SelfiesMaxBranchingFanoutDescriptor();
    std::variant<double, int, std::string> calculateFromSelfies(const std::string& selfiesStr) const override;
};

class SelfiesEmptyBranchesDescriptor : public SelfiesDescriptor {
public:
    SelfiesEmptyBranchesDescriptor();
    std::variant<double, int, std::string> calculateFromSelfies(const std::string& selfiesStr) const override;
};

class SelfiesRingTokenFrequencyDescriptor : public SelfiesDescriptor {
public:
    SelfiesRingTokenFrequencyDescriptor();
    std::variant<double, int, std::string> calculateFromSelfies(const std::string& selfiesStr) const override;
};

class SelfiesTokenEntropyDescriptor : public SelfiesDescriptor {
public:
    SelfiesTokenEntropyDescriptor();
    std::variant<double, int, std::string> calculateFromSelfies(const std::string& selfiesStr) const override;
};

class SelfiesStringLengthDescriptor : public SelfiesDescriptor {
public:
    SelfiesStringLengthDescriptor();
    std::variant<double, int, std::string> calculateFromSelfies(const std::string& selfiesStr) const override;
};

class SelfiesBranchTokenRatioDescriptor : public SelfiesDescriptor {
public:
    SelfiesBranchTokenRatioDescriptor();
    std::variant<double, int, std::string> calculateFromSelfies(const std::string& selfiesStr) const override;
};

class SelfiesRingTokenRatioDescriptor : public SelfiesDescriptor {
public:
    SelfiesRingTokenRatioDescriptor();
    std::variant<double, int, std::string> calculateFromSelfies(const std::string& selfiesStr) const override;
};

class SelfiesBigramCountDescriptor : public SelfiesDescriptor {
public:
    SelfiesBigramCountDescriptor();
    std::variant<double, int, std::string> calculateFromSelfies(const std::string& selfiesStr) const override;
};

class SelfiesTrigramCountDescriptor : public SelfiesDescriptor {
public:
    SelfiesTrigramCountDescriptor();
    std::variant<double, int, std::string> calculateFromSelfies(const std::string& selfiesStr) const override;
};

class SelfiesTokenRepetitionRateDescriptor : public SelfiesDescriptor {
public:
    SelfiesTokenRepetitionRateDescriptor();
    std::variant<double, int, std::string> calculateFromSelfies(const std::string& selfiesStr) const override;
};

class SelfiesAromaticTokenRatioDescriptor : public SelfiesDescriptor {
public:
    SelfiesAromaticTokenRatioDescriptor();
    std::variant<double, int, std::string> calculateFromSelfies(const std::string& selfiesStr) const override;
};

class SelfiesChargeTokenRatioDescriptor : public SelfiesDescriptor {
public:
    SelfiesChargeTokenRatioDescriptor();
    std::variant<double, int, std::string> calculateFromSelfies(const std::string& selfiesStr) const override;
};

class SelfiesAliphaticChainLengthDescriptor : public SelfiesDescriptor {
public:
    SelfiesAliphaticChainLengthDescriptor();
    std::variant<double, int, std::string> calculateFromSelfies(const std::string& selfiesStr) const override;
};

class SelfiesHeteroatomRatioDescriptor : public SelfiesDescriptor {
public:
    SelfiesHeteroatomRatioDescriptor();
    std::variant<double, int, std::string> calculateFromSelfies(const std::string& selfiesStr) const override;
};

class SelfiesMaxRingSizeDescriptor : public SelfiesDescriptor {
public:
    SelfiesMaxRingSizeDescriptor();
    std::variant<double, int, std::string> calculateFromSelfies(const std::string& selfiesStr) const override;
};

class SelfiesStereoTokenRatioDescriptor : public SelfiesDescriptor {
public:
    SelfiesStereoTokenRatioDescriptor();
    std::variant<double, int, std::string> calculateFromSelfies(const std::string& selfiesStr) const override;
};

class SelfiesGrammarDepthComplexityDescriptor : public SelfiesDescriptor {
public:
    SelfiesGrammarDepthComplexityDescriptor();
    std::variant<double, int, std::string> calculateFromSelfies(const std::string& selfiesStr) const override;
};

class SelfiesVocabularyDiversityDescriptor : public SelfiesDescriptor {
public:
    SelfiesVocabularyDiversityDescriptor();
    std::variant<double, int, std::string> calculateFromSelfies(const std::string& selfiesStr) const override;
};

class SelfiesBondComplexityDescriptor : public SelfiesDescriptor {
public:
    SelfiesBondComplexityDescriptor();
    std::variant<double, int, std::string> calculateFromSelfies(const std::string& selfiesStr) const override;
};

class SelfiesSequenceSymmetryDescriptor : public SelfiesDescriptor {
public:
    SelfiesSequenceSymmetryDescriptor();
    std::variant<double, int, std::string> calculateFromSelfies(const std::string& selfiesStr) const override;
};

class SelfiesTokenLengthVarianceDescriptor : public SelfiesDescriptor {
public:
    SelfiesTokenLengthVarianceDescriptor();
    std::variant<double, int, std::string> calculateFromSelfies(const std::string& selfiesStr) const override;
};

class SelfiesFormalChargeEntropyDescriptor : public SelfiesDescriptor {
public:
    SelfiesFormalChargeEntropyDescriptor();
    std::variant<double, int, std::string> calculateFromSelfies(const std::string& selfiesStr) const override;
};

class SelfiesElementDiversityDescriptor : public SelfiesDescriptor {
public:
    SelfiesElementDiversityDescriptor();
    std::variant<double, int, std::string> calculateFromSelfies(const std::string& selfiesStr) const override;
};

class SelfiesLevenshteinCompressionDescriptor : public SelfiesDescriptor {
public:
    SelfiesLevenshteinCompressionDescriptor();
    std::variant<double, int, std::string> calculateFromSelfies(const std::string& selfiesStr) const override;
};

class SelfiesLoopComplexityDescriptor : public SelfiesDescriptor {
public:
    SelfiesLoopComplexityDescriptor();
    std::variant<double, int, std::string> calculateFromSelfies(const std::string& selfiesStr) const override;
};

class SelfiesUnpairedElectronRatioDescriptor : public SelfiesDescriptor {
public:
    SelfiesUnpairedElectronRatioDescriptor();
    std::variant<double, int, std::string> calculateFromSelfies(const std::string& selfiesStr) const override;
};

} // namespace descriptors
} // namespace desfact
