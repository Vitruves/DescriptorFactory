// include/descriptors/morgan.hpp
#pragma once

#include "descriptors.hpp" // Base Descriptor class
#include <GraphMol/Fingerprints/MorganFingerprints.h> // Morgan specifics
#include <DataStructs/ExplicitBitVect.h> // ExplicitBitVect
#include <GraphMol/PartialCharges/GasteigerCharges.h> // Needed for 

namespace desfact {
namespace descriptors {

// --- Morgan Fingerprint Based Descriptors ---

class MorganBitDensityRatioDescriptor : public Descriptor {
public:
    MorganBitDensityRatioDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
private:
    unsigned int radius = 2; // Default radius for Morgan
    unsigned int nBits = 1024; // User specified
};

class MorganFragmentUniquenessScoreDescriptor : public Descriptor {
public:
    MorganFragmentUniquenessScoreDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
private:
    unsigned int radius = 2;
    unsigned int nBits = 1024;
};

class MorganRadiusInformationRatioDescriptor : public Descriptor {
public:
    MorganRadiusInformationRatioDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
private:
    unsigned int radius1 = 1;
    unsigned int radius2 = 2;
    unsigned int nBits = 1024;
};

class MorganBitClusteringCoefficientDescriptor : public Descriptor {
public:
    MorganBitClusteringCoefficientDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
private:
    unsigned int radius = 2;
    unsigned int nBits = 1024;
};

class MorganFrequentFragmentEntropyDescriptor : public Descriptor {
public:
    MorganFrequentFragmentEntropyDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
private:
    unsigned int radius = 2;
    unsigned int nBits = 1024;
};

class MorganFingerprintAsymmetryIndexDescriptor : public Descriptor {
public:
    MorganFingerprintAsymmetryIndexDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
private:
    unsigned int radius = 2;
    unsigned int nBits = 1024;
};

class MorganBitTransitionRateDescriptor : public Descriptor {
public:
    MorganBitTransitionRateDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
private:
    unsigned int radius = 2;
    unsigned int nBits = 1024;
};

class MorganFragmentSizeDistributionSkewnessDescriptor : public Descriptor {
public:
    MorganFragmentSizeDistributionSkewnessDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
private:
    unsigned int radius = 2;
    unsigned int nBits = 1024;
};

class MorganLongestCommonSubstructureScoreDescriptor : public Descriptor {
public:
    MorganLongestCommonSubstructureScoreDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
private:
    unsigned int radius = 2;
    unsigned int nBits = 1024;
};

class MorganPharmacophorePatternDensityDescriptor : public Descriptor {
public:
    MorganPharmacophorePatternDensityDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
private:
    unsigned int radius = 2;
    unsigned int nBits = 1024;
    // Note: Requires predefined pharmacophore SMARTS/SMIRKS
};

class MorganFragmentDiversityScoreDescriptor : public Descriptor {
public:
    MorganFragmentDiversityScoreDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
private:
    unsigned int radius = 2;
    unsigned int nBits = 1024;
};

class MorganBitCorrelationCoefficientDescriptor : public Descriptor {
public:
    MorganBitCorrelationCoefficientDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
private:
    unsigned int radius = 2;
    unsigned int nBits = 1024;
};

class MorganPatternRecurrenceFrequencyDescriptor : public Descriptor {
public:
    MorganPatternRecurrenceFrequencyDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
private:
    unsigned int radius = 2;
    unsigned int nBits = 1024;
};

class MorganBitSimilarityToReferenceSetDescriptor : public Descriptor {
public:
    MorganBitSimilarityToReferenceSetDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
private:
    unsigned int radius = 2;
    unsigned int nBits = 1024;
};

class MorganFragmentComplexityScoreDescriptor : public Descriptor {
public:
    MorganFragmentComplexityScoreDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
private:
    unsigned int radius = 2;
    unsigned int nBits = 1024;
};

class MorganInformationContentDensityDescriptor : public Descriptor {
public:
    MorganInformationContentDensityDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
private:
    unsigned int radius = 2;
    unsigned int nBits = 1024;
};

class MorganEnvironmentVariabilityIndexDescriptor : public Descriptor {
public:
    MorganEnvironmentVariabilityIndexDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
private:
    unsigned int radius = 2; // Radius for comparison environments
    unsigned int nBits = 1024;
};

class MorganBitPositionImportanceScoreDescriptor : public Descriptor {
public:
    MorganBitPositionImportanceScoreDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
private:
    unsigned int radius = 2;
    unsigned int nBits = 1024;
};

class MorganRingSystemRepresentationScoreDescriptor : public Descriptor {
public:
    MorganRingSystemRepresentationScoreDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
private:
    unsigned int radius = 2;
    unsigned int nBits = 1024;
};

class MorganBitProbabilityDistributionEntropyDescriptor : public Descriptor {
public:
    MorganBitProbabilityDistributionEntropyDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
private:
    unsigned int radius = 2;
    unsigned int nBits = 1024;
};

class MorganFragmentElectronegativitySpectrumDescriptor : public Descriptor {
public:
    MorganFragmentElectronegativitySpectrumDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
private:
    unsigned int radius = 2;
    unsigned int nBits = 1024;
};

class MorganBitPositionCorrelationMatrixDescriptor : public Descriptor {
public:
    MorganBitPositionCorrelationMatrixDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
private:
    unsigned int radius = 2;
    unsigned int nBits = 1024;
};

class MorganFragmentConnectivityPatternDescriptor : public Descriptor {
public:
    MorganFragmentConnectivityPatternDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
private:
    unsigned int radius = 2;
    unsigned int nBits = 1024;
};

class MorganBitOccurrenceFrequencySkewnessDescriptor : public Descriptor {
public:
    MorganBitOccurrenceFrequencySkewnessDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
private:
    unsigned int radius = 2;
    unsigned int nBits = 1024;
};

class MorganSubstructureHeterogeneityIndexDescriptor : public Descriptor {
public:
    MorganSubstructureHeterogeneityIndexDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
private:
    unsigned int radius = 2;
    unsigned int nBits = 1024;
};

class MorganBitPolarityDistributionDescriptor : public Descriptor {
public:
    MorganBitPolarityDistributionDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
private:
    unsigned int radius = 2;
    unsigned int nBits = 1024;
};

class MorganFragmentSimilarityNetworkDescriptor : public Descriptor {
public:
    MorganFragmentSimilarityNetworkDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
private:
    unsigned int radius = 2;
    unsigned int nBits = 1024;
};

class MorganBitPositionInformationGainDescriptor : public Descriptor {
public:
    MorganBitPositionInformationGainDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
private:
    unsigned int radius = 2;
    unsigned int nBits = 1024;
};

class MorganSubstructureDiversityGradientDescriptor : public Descriptor {
public:
    MorganSubstructureDiversityGradientDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
private:
    unsigned int radius1 = 1;
    unsigned int radius2 = 2;
    unsigned int nBits = 1024;
};

class MorganFingerprintSymmetryScoreDescriptor : public Descriptor {
public:
    MorganFingerprintSymmetryScoreDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
private:
    unsigned int radius = 2;
    unsigned int nBits = 1024;
};


} // namespace descriptors
} // namespace desfact