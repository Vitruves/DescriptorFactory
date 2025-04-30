// src/descriptors/morgan.cpp
#include "descriptors/morgan.hpp"
#include "utils.hpp" // Includes common headers and logger
#include "descriptors/sum.hpp" // For getAtomElectronegativity
#include <GraphMol/GraphMol.h>
#include <GraphMol/MolOps.h>
#include <RDGeneral/export.h>
#include <GraphMol/SanitException.h>
#include <DataManip/MetricMatrixCalc/MetricMatrixCalc.h>
#include <GraphMol/RingInfo.h> // Needed for 
#include <GraphMol/PartialCharges/GasteigerCharges.h> // Needed for 
#include <GraphMol/SmilesParse/SmilesParse.h> // For MolToSmiles if needed for fragments
#include <GraphMol/Substruct/SubstructMatch.h> // For pharmacophores if implemented
#include <DataStructs/ExplicitBitVect.h>
#include <DataStructs/BitOps.h> // For TanimotoSimilarity
#include <DataStructs/SparseIntVect.h> // Added for getFingerprint
#include <GraphMol/Fingerprints/MorganFingerprints.h>

#include <numeric>
#include <cmath>
#include <vector>
#include <string>
#include <set>
#include <map>
#include <memory> // For unique_ptr
#include <queue> // For queue
#include <limits> // For infinity()

namespace desfact {
namespace descriptors {

// --- Anonymous Namespace for Helpers ---
namespace {

    const std::string REFERENCE_SMILES = "CCN(CC)C(=O)C1CN(C)C2CC3=CNC4=CC=CC(=C34)C2=C1";
    
    // FORWARD DECLARATIONS
    std::unique_ptr<ExplicitBitVect> getMorganFingerprintBV(
        const RDKit::ROMol& mol, unsigned int radius, unsigned int nBits);
    
    std::unique_ptr<std::map<uint32_t,int>> getMorganFingerprintCounts(
        const RDKit::ROMol& mol, unsigned int radius);
    
    // Cache for reference fingerprints to avoid regenerating
    RDKit::ROMol* getReferenceMolecule() {
        static RDKit::ROMol* referenceMol = nullptr;
        if (!referenceMol) {
            try {
                referenceMol = RDKit::SmilesToMol(REFERENCE_SMILES);
                if (!referenceMol) {
                    globalLogger.error("Failed to create reference molecule");
                }
            } catch (const std::exception& e) {
                globalLogger.error("Exception creating reference: " + std::string(e.what()));
                referenceMol = nullptr;
            }
        }
        return referenceMol;
    }
    
    // Get reference fingerprint as bit vector
    std::unique_ptr<ExplicitBitVect> getReferenceFingerprintBV(unsigned int radius, unsigned int nBits) {
        static std::map<std::pair<unsigned int, unsigned int>, std::unique_ptr<ExplicitBitVect>> cache;
        
        auto key = std::make_pair(radius, nBits);
        if (cache.find(key) == cache.end()) {
            RDKit::ROMol* refMol = getReferenceMolecule();
            if (refMol) {
                cache[key] = getMorganFingerprintBV(*refMol, radius, nBits);
            }
        }
        
        if (cache.find(key) != cache.end()) {
            // Return a copy of the cached fingerprint
            return std::make_unique<ExplicitBitVect>(*cache[key]);
        }
        return nullptr;
    }
    
    // Get reference fingerprint counts
    std::unique_ptr<std::map<uint32_t,int>> getReferenceFingerprintCounts(unsigned int radius) {
        static std::map<unsigned int, std::unique_ptr<std::map<uint32_t,int>>> cache;
        
        if (cache.find(radius) == cache.end()) {
            RDKit::ROMol* refMol = getReferenceMolecule();
            if (refMol) {
                cache[radius] = getMorganFingerprintCounts(*refMol, radius);
            }
        }
        
        if (cache.find(radius) != cache.end()) {
            // Return a copy
            return std::make_unique<std::map<uint32_t,int>>(*cache[radius]);
        }
        return nullptr;
    }

    // Helper to get Morgan fingerprint as Bit Vector
    std::unique_ptr<ExplicitBitVect> getMorganFingerprintBV(
        const RDKit::ROMol& mol, unsigned int radius, unsigned int nBits)
    {
        try {
             ExplicitBitVect* fp = RDKit::MorganFingerprints::getFingerprintAsBitVect(
                 mol, radius, nBits, 
                 nullptr,  // atomIds
                 nullptr,  // ignoreAtoms
                 false,    // useChirality 
                 true,     // useBondTypes
                 false,    // useFeatures
                 nullptr,  // atomsSettingBits
                 false     // includeRedundantEnvironments - MUST BE bool, not nullptr
             );
             return std::unique_ptr<ExplicitBitVect>(fp);
        } catch (const std::exception& e) {
            globalLogger.debug("Morgan FP (BV) calculation failed: " + std::string(e.what()));
            return nullptr;
        } catch (...) {
            globalLogger.debug("Morgan FP (BV) calculation failed due to unknown error.");
            return nullptr;
        }
    }

    // Helper to get Morgan fingerprint counts
    std::unique_ptr<std::map<uint32_t,int>> getMorganFingerprintCounts(
        const RDKit::ROMol& mol, unsigned int radius)
    {
        auto info = std::make_unique<std::map<uint32_t,int>>();
        try {
            // Match exact parameter order from RDKit documentation:
            // mol, radius, invariants, fromAtoms, useChirality, useBondTypes, useCounts,  
            // onlyNonzeroInvariants, atomsSettingBits, includeRedundantEnvironments
            RDKit::SparseIntVect<uint32_t>* fp =
                RDKit::MorganFingerprints::getFingerprint(
                    mol, radius,
                    nullptr,  // invariants
                    nullptr,  // fromAtoms
                    false,    // useChirality
                    true,     // useBondTypes
                    true,     // useCounts
                    false,    // onlyNonzeroInvariants 
                    nullptr,  // atomsSettingBits
                    false     // includeRedundantEnvironments
                );
            if (fp) {
                for (auto [bit, cnt] : fp->getNonzeroElements()) {
                    (*info)[bit] = cnt;
                }
                delete fp;
            }
            return info;
        } catch (...) {
            globalLogger.debug("Morgan FP (Counts) calculation failed");
            return nullptr;
        }
    }

    // Helper to get Morgan fingerprint bit info
    std::unique_ptr<RDKit::MorganFingerprints::BitInfoMap> getMorganFingerprintBitInfo(
        const RDKit::ROMol& mol, unsigned int radius, unsigned int nBits)
    {
        auto bitInfo = std::make_unique<RDKit::MorganFingerprints::BitInfoMap>();
        try {
            // IMPORTANT: getFingerprintAsBitVect takes EXACTLY 10 args
            // The bitInfo map is passed as the 9th parameter (atomsSettingBits)
            ExplicitBitVect* fp = RDKit::MorganFingerprints::getFingerprintAsBitVect(
                mol, radius, nBits,
                nullptr,  // fromAtoms
                nullptr,  // ignoreAtoms
                false,    // useChirality
                true,     // useBondTypes
                false,    // useFeatures
                bitInfo.get(), // atomsSettingBits
                false     // includeRedundantEnvironments - MUST BE bool, not nullptr
            );
            if (fp) {
                delete fp;
                return bitInfo;
            } else {
                return nullptr;
            }
        } catch (...) {
            globalLogger.debug("Morgan FP (BitInfo) calculation failed");
            return nullptr;
        }
    }

    // Get Gasteiger partial charges
    bool getGasteigerCharges(const RDKit::ROMol* rdkMol, std::vector<double>& charges) {
        charges.assign(rdkMol->getNumAtoms(), 0.0);
        try {
            // Gasteiger charges require a writeable molecule, so we copy
            RDKit::RWMol nonConstMol(*rdkMol);
            RDKit::MolOps::sanitizeMol(nonConstMol); // Ensure sanitization
            RDKit::computeGasteigerCharges(nonConstMol, 12, true); // 12 iterations, throw on error

            for (unsigned int i = 0; i < rdkMol->getNumAtoms(); ++i) {
                 if (nonConstMol.getAtomWithIdx(i)->hasProp(RDKit::common_properties::_GasteigerCharge)) {
                    charges[i] = nonConstMol.getAtomWithIdx(i)->getProp<double>(RDKit::common_properties::_GasteigerCharge);
                 } else {
                     globalLogger.warning("Gasteiger charge calculation succeeded but charge property missing for atom " + std::to_string(i));
                     charges[i] = 0.0; // Assign default
                 }
            }
            return true;
        } catch (const std::exception& e) {
            globalLogger.debug("Failed to compute Gasteiger charges: " + std::string(e.what()));
            return false;
        } catch (...) {
            globalLogger.debug("Failed to compute Gasteiger charges due to unknown error.");
            return false;
        }
    }

    // Helper for Shannon Entropy calculation
    double calculateShannonEntropyInternal(const std::vector<int>& counts) {
        double totalCount = std::accumulate(counts.begin(), counts.end(), 0.0);
        if (totalCount <= 1.0) return 0.0;
        double entropy = 0.0;
        for (int count : counts) {
            if (count > 0) {
                double p = static_cast<double>(count) / totalCount;
                entropy -= p * std::log2(p);
            }
        }
        return entropy;
    }
    // Overload for maps (value type int)
    template <typename Key>
    double calculateShannonEntropyInternal(const std::map<Key, int>& counts) {
        double totalCount = 0;
        for(const auto& pair : counts) {
            totalCount += pair.second;
        }
        if (totalCount <= 1.0) return 0.0;
        double entropy = 0.0;
        for (const auto& pair : counts) {
            if (pair.second > 0) {
                double p = static_cast<double>(pair.second) / totalCount;
                entropy -= p * std::log2(p);
            }
        }
        return entropy;
    }
     // Overload for maps (value type double)
    template <typename Key>
    double calculateShannonEntropyInternal(const std::map<Key, double>& counts) {
        double totalCount = 0;
        for(const auto& pair : counts) {
            totalCount += pair.second;
        }
        if (totalCount <= 1.0) return 0.0;
        double entropy = 0.0;
        for (const auto& pair : counts) {
            if (pair.second > 1e-9) { // Avoid log(0) for floating point
                double p = pair.second / totalCount;
                entropy -= p * std::log2(p);
            }
        }
        return entropy;
    }


    // Helper for Skewness calculation
    double calculateSkewnessInternal(const std::vector<double>& values) {
        size_t n = values.size();
        if (n < 3) return 0.0;
        double sum = std::accumulate(values.begin(), values.end(), 0.0);
        double mean = sum / n;
        double m3 = 0.0;
        double m2 = 0.0;
        for (double val : values) {
            double diff = val - mean;
            m2 += diff * diff;
            m3 += diff * diff * diff;
        }
        m2 /= n;
        m3 /= n;
        if (m2 < 1e-9) return 0.0;
        double stddev = std::sqrt(m2);
        return m3 / (stddev * stddev * stddev);
    }

    // Tanimoto similarity for ExplicitBitVect
    double calculateTanimotoSimilarity(const ExplicitBitVect& bv1, const ExplicitBitVect& bv2) {
        // Use RDDataManip namespace instead of RDKit::DataStructs
        return RDDataManip::TanimotoSimilarityMetric(bv1, bv2, 0);  // 0 for dim parameter as it's unused
    }

} // end anonymous namespace

// --- Implementations ---

// MorganBitDensityRatio
MorganBitDensityRatioDescriptor::MorganBitDensityRatioDescriptor()
    : Descriptor("MorganBitDensityRatio", "Ratio of set bits to total bits in Morgan fingerprint (r=2, n=1024)") {}

std::variant<double, int, std::string> MorganBitDensityRatioDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit Mol is null";

    auto fp = getMorganFingerprintBV(*rdkMol, radius, nBits);
    if (!fp) return "Error: FP calculation failed";

    if (nBits == 0) return 0.0; // Avoid division by zero
    return static_cast<double>(fp->getNumOnBits()) / nBits;
}

// MorganFragmentUniquenessScore
MorganFragmentUniquenessScoreDescriptor::MorganFragmentUniquenessScoreDescriptor()
    : Descriptor("MorganFragUniqueness", "Score based on uniqueness of fragments (r=2, n=1024) vs reference (unimplemented)") {}

std::variant<double, int, std::string> MorganFragmentUniquenessScoreDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit Mol is null";

    // Get fingerprint counts for current molecule
    auto counts = getMorganFingerprintCounts(*rdkMol, radius);
    if (!counts) return "Error: FP calculation failed";
    
    // Get reference molecule fingerprint counts
    auto refCounts = getReferenceFingerprintCounts(radius);
    if (!refCounts) return "Error: Reference FP calculation failed";
    
    // Calculate uniqueness as 1 - Tanimoto similarity between counts
    double sumMin = 0.0, sumMol = 0.0, sumRef = 0.0;
    std::set<uint32_t> allBits;
    
    for (const auto& [bit, count] : *counts) {
        allBits.insert(bit);
        sumMol += count;
    }
    
    for (const auto& [bit, count] : *refCounts) {
        allBits.insert(bit);
        sumRef += count;
    }
    
    for (uint32_t bit : allBits) {
        int molCount = counts->count(bit) ? (*counts)[bit] : 0;
        int refCount = refCounts->count(bit) ? (*refCounts)[bit] : 0;
        sumMin += std::min(molCount, refCount);
    }
    
    if (sumMol + sumRef - sumMin <= 0) return 1.0;
    double similarity = sumMin / (sumMol + sumRef - sumMin);
    
    // Return uniqueness (1 - similarity)
    return 1.0 - similarity;
}

// MorganRadiusInformationRatio
MorganRadiusInformationRatioDescriptor::MorganRadiusInformationRatioDescriptor()
    : Descriptor("MorganRadiusInfoRatio", "Ratio of information content (entropy) r=2/r=1 (n=1024)") {}

std::variant<double, int, std::string> MorganRadiusInformationRatioDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit Mol is null";

    auto counts1 = getMorganFingerprintCounts(*rdkMol, radius1);
    auto counts2 = getMorganFingerprintCounts(*rdkMol, radius2);

    if (!counts1 || !counts2) return "Error: FP calculation failed";

    double entropy1 = calculateShannonEntropyInternal(*counts1);
    double entropy2 = calculateShannonEntropyInternal(*counts2);

    if (std::abs(entropy1) < 1e-9) { // Avoid division by zero
         return (std::abs(entropy2) < 1e-9) ? 1.0 : std::numeric_limits<double>::infinity(); // Or some large number
    }
    return entropy2 / entropy1;
}

// MorganBitClusteringCoefficient
MorganBitClusteringCoefficientDescriptor::MorganBitClusteringCoefficientDescriptor()
    : Descriptor("MorganBitClustering", "Measure of bit clustering in Morgan FP (r=2, n=1024)") {}

std::variant<double, int, std::string> MorganBitClusteringCoefficientDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit Mol is null";

    auto fp = getMorganFingerprintBV(*rdkMol, radius, nBits);
    if (!fp) return "Error: FP calculation failed";

    long long totalSetNeighbors = 0;
    long long numSetBits = fp->getNumOnBits();
    if (numSetBits <= 1) return 0.0; // Clustering requires at least 2 set bits

    for (unsigned int i = 0; i < nBits; ++i) {
        if (fp->getBit(i)) {
            if (i > 0 && fp->getBit(i - 1)) {
                totalSetNeighbors++;
            }
            if (i < nBits - 1 && fp->getBit(i + 1)) {
                totalSetNeighbors++;
            }
        }
    }
    // Max possible set neighbors is 2 * numSetBits (if all clumped)
    // Normalize? Let's return avg number of set neighbors per set bit.
    return static_cast<double>(totalSetNeighbors) / numSetBits;
}

// MorganFrequentFragmentEntropy
MorganFrequentFragmentEntropyDescriptor::MorganFrequentFragmentEntropyDescriptor()
    : Descriptor("MorganFragEntropy", "Entropy of Morgan fragment counts (r=2)") {}

std::variant<double, int, std::string> MorganFrequentFragmentEntropyDescriptor::calculate(const Molecule& mol) const {
     if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit Mol is null";

    auto counts = getMorganFingerprintCounts(*rdkMol, radius);
    if (!counts) return "Error: FP calculation failed";

    return calculateShannonEntropyInternal(*counts);
}

// MorganFingerprintAsymmetryIndex
MorganFingerprintAsymmetryIndexDescriptor::MorganFingerprintAsymmetryIndexDescriptor()
    : Descriptor("MorganFPAsymmetry", "Asymmetry measure of set bits across fingerprint (r=2, n=1024)") {}

std::variant<double, int, std::string> MorganFingerprintAsymmetryIndexDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit Mol is null";

    auto fp = getMorganFingerprintBV(*rdkMol, radius, nBits);
    if (!fp) return "Error: FP calculation failed";

    unsigned int halfPoint = nBits / 2;
    long long countFirstHalf = 0;
    long long countSecondHalf = 0;

    for (unsigned int i = 0; i < halfPoint; ++i) {
        if (fp->getBit(i)) countFirstHalf++;
    }
    for (unsigned int i = halfPoint; i < nBits; ++i) {
         if (fp->getBit(i)) countSecondHalf++;
    }

    long long totalBits = fp->getNumOnBits();
    if (totalBits == 0) return 0.0;

    // Calculate asymmetry: absolute difference normalized by total
    return static_cast<double>(std::abs(countFirstHalf - countSecondHalf)) / totalBits;
}

// MorganBitTransitionRate
MorganBitTransitionRateDescriptor::MorganBitTransitionRateDescriptor()
    : Descriptor("MorganBitTransitionRate", "Rate of 0->1 and 1->0 transitions in Morgan FP (r=2, n=1024)") {}

std::variant<double, int, std::string> MorganBitTransitionRateDescriptor::calculate(const Molecule& mol) const {
     if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit Mol is null";

    auto fp = getMorganFingerprintBV(*rdkMol, radius, nBits);
    if (!fp) return "Error: FP calculation failed";

    if (nBits <= 1) return 0.0;

    int transitions = 0;
    for (unsigned int i = 0; i < nBits - 1; ++i) {
        if (fp->getBit(i) != fp->getBit(i + 1)) {
            transitions++;
        }
    }

    // Rate is transitions / (total possible transitions)
    return static_cast<double>(transitions) / (nBits - 1);
}

// MorganFragmentSizeDistributionSkewness
MorganFragmentSizeDistributionSkewnessDescriptor::MorganFragmentSizeDistributionSkewnessDescriptor()
    : Descriptor("MorganFragSizeSkew", "Skewness of fragment size distribution (r=2, n=1024)") {}

std::variant<double, int, std::string> MorganFragmentSizeDistributionSkewnessDescriptor::calculate(const Molecule& mol) const {
     if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit Mol is null";

    auto bitInfo = getMorganFingerprintBitInfo(*rdkMol, radius, nBits);
    if (!bitInfo) return "Error: FP calculation failed";

    std::vector<double> fragmentSizes;
    std::set<unsigned int> processedAtomIndices; // Avoid double counting atoms if they appear in multiple fragments

    for(const auto& pair : *bitInfo) { // pair: <bitId, vector<pair<atomIdx, radius>>>
        const auto& atomRadiusVec = pair.second;
        if (atomRadiusVec.empty()) continue;

        unsigned int centerAtomIdx = atomRadiusVec[0].first;

        unsigned int fragmentRadius = atomRadiusVec[0].second;
        fragmentSizes.push_back(static_cast<double>(fragmentRadius));
    }

    if (fragmentSizes.size() < 3) return 0.0; // Skewness undefined

    return calculateSkewnessInternal(fragmentSizes);
}

// MorganLongestCommonSubstructureScore
MorganLongestCommonSubstructureScoreDescriptor::MorganLongestCommonSubstructureScoreDescriptor()
    : Descriptor("MorganLCSScore", "LCS score based on Morgan fragments (r=2, n=1024) (unimplemented)") {}

std::variant<double, int, std::string> MorganLongestCommonSubstructureScoreDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit Mol is null";

    auto fp = getMorganFingerprintBV(*rdkMol, radius, nBits);
    auto refFp = getReferenceFingerprintBV(radius, nBits);
    
    if (!fp || !refFp) return 0.0;
    
    // Find longest sequence of common bits as a simple measure of common substructure
    int maxCommonLength = 0;
    int currentLength = 0;
    
    for (unsigned int i = 0; i < nBits; ++i) {
        if (fp->getBit(i) && refFp->getBit(i)) {
            currentLength++;
            maxCommonLength = std::max(maxCommonLength, currentLength);
        } else {
            currentLength = 0;
        }
    }
    
    // Normalize by total possible length
    return static_cast<double>(maxCommonLength) / nBits;
}

// MorganPharmacophorePatternDensity
MorganPharmacophorePatternDensityDescriptor::MorganPharmacophorePatternDensityDescriptor()
    : Descriptor("MorganPharmacophoreDensity", "Density of fragments matching pharmacophore patterns (r=2, n=1024) (unimplemented)") {}

std::variant<double, int, std::string> MorganPharmacophorePatternDensityDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit Mol is null";
    
    // Basic pharmacophore bit "detector" based on atom properties
    int pharmacophoreBits = 0;
    auto bitInfo = getMorganFingerprintBitInfo(*rdkMol, radius, nBits);
    if (!bitInfo) return 0.0;
    
    // Get reference fingerprint for normalization
    auto refFp = getReferenceFingerprintBV(radius, nBits);
    if (!refFp) return 0.0;
    
    for (const auto& pair : *bitInfo) {
        const auto& atomRadiusVec = pair.second;
        if (atomRadiusVec.empty()) continue;
        
        unsigned int centerAtomIdx = atomRadiusVec[0].first;
        const RDKit::Atom* atom = rdkMol->getAtomWithIdx(centerAtomIdx);
        
        // Check for pharmacophore features - basic check for common types
        bool isPharmacophore = false;
        
        // H-bond acceptor
        if (atom->getAtomicNum() == 7 || atom->getAtomicNum() == 8) {
            isPharmacophore = true;
        }
        // H-bond donor
        else if ((atom->getAtomicNum() == 7 || atom->getAtomicNum() == 8) && 
                atom->getTotalNumHs() > 0) {
            isPharmacophore = true;
        }
        // Basic check for aromatic
        else if (atom->getIsAromatic()) {
            isPharmacophore = true;
        }
        
        if (isPharmacophore) {
            pharmacophoreBits++;
        }
    }
    
    // Normalize by reference fingerprint set bits for comparison
    int refBitsOn = refFp->getNumOnBits();
    if (refBitsOn == 0) return 0.0;
    
    return static_cast<double>(pharmacophoreBits) / refBitsOn;
}

// MorganFragmentDiversityScore
MorganFragmentDiversityScoreDescriptor::MorganFragmentDiversityScoreDescriptor()
    : Descriptor("MorganFragDiversity", "Diversity score (Shannon entropy) of Morgan fragments (r=2)") {}

std::variant<double, int, std::string> MorganFragmentDiversityScoreDescriptor::calculate(const Molecule& mol) const {
     if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit Mol is null";

    auto counts = getMorganFingerprintCounts(*rdkMol, radius);
    if (!counts) return "Error: FP calculation failed";

    // Shannon entropy of the fragment counts directly measures diversity
    return calculateShannonEntropyInternal(*counts);
}

// MorganBitCorrelationCoefficient
MorganBitCorrelationCoefficientDescriptor::MorganBitCorrelationCoefficientDescriptor()
    : Descriptor("MorganBitCorrelation", "Correlation between bits (r=2, n=1024) (requires dataset - unimplemented)") {}

std::variant<double, int, std::string> MorganBitCorrelationCoefficientDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit Mol is null";
    
    auto fp = getMorganFingerprintBV(*rdkMol, radius, nBits);
    auto refFp = getReferenceFingerprintBV(radius, nBits);
    
    if (!fp || !refFp) return 0.0;
    
    // Calculate correlation coefficient between fingerprints
    int n11 = 0, n10 = 0, n01 = 0, n00 = 0;
    
    for (unsigned int i = 0; i < nBits; ++i) {
        bool molBit = fp->getBit(i);
        bool refBit = refFp->getBit(i);
        
        if (molBit && refBit) n11++;
        else if (molBit && !refBit) n10++;
        else if (!molBit && refBit) n01++;
        else n00++;
    }
    
    // Compute Matthews correlation coefficient
    double numerator = (n11 * n00) - (n10 * n01);
    double denominator = std::sqrt(static_cast<double>((n11 + n10) * (n11 + n01) * (n00 + n10) * (n00 + n01)));
    
    if (std::abs(denominator) < 1e-9) return 0.0;
    
    return numerator / denominator;
}

// MorganPatternRecurrenceFrequency
MorganPatternRecurrenceFrequencyDescriptor::MorganPatternRecurrenceFrequencyDescriptor()
    : Descriptor("MorganPatternRecurrence", "Frequency of recurring fragments (r=2) (Avg Count > 1)") {}

std::variant<double, int, std::string> MorganPatternRecurrenceFrequencyDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit Mol is null";

    auto counts = getMorganFingerprintCounts(*rdkMol, radius);
    if (!counts || counts->empty()) return 0.0;

    double totalFragments = 0;
    double recurringFragmentsCount = 0; // Count identifiers that appear more than once
    for(const auto& pair : *counts) {
        totalFragments += pair.second;
        if (pair.second > 1) {
            recurringFragmentsCount++;
        }
    }

    if (counts->size() == 0) return 0.0; // Check size instead of empty() after potentially modifying
    // Return fraction of unique fragments that recur
    return recurringFragmentsCount / counts->size();
}

// MorganBitSimilarityToReferenceSet
MorganBitSimilarityToReferenceSetDescriptor::MorganBitSimilarityToReferenceSetDescriptor()
    : Descriptor("MorganSimilarityToRef", "Similarity to reference set (r=2, n=1024) (unimplemented)") {}

std::variant<double, int, std::string> MorganBitSimilarityToReferenceSetDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit Mol is null";
    
    // Get bit vector fingerprint for current molecule
    auto fp = getMorganFingerprintBV(*rdkMol, radius, nBits);
    if (!fp) return "Error: FP calculation failed";
    
    // Get reference molecule fingerprint
    auto refFp = getReferenceFingerprintBV(radius, nBits);
    if (!refFp) return "Error: Reference FP calculation failed";
    
    // Calculate Tanimoto similarity to reference
    return calculateTanimotoSimilarity(*fp, *refFp);
}

// MorganFragmentComplexityScore
MorganFragmentComplexityScoreDescriptor::MorganFragmentComplexityScoreDescriptor()
    : Descriptor("MorganFragComplexity", "Avg complexity (atom count) of fragments (r=2, n=1024)") {}

std::variant<double, int, std::string> MorganFragmentComplexityScoreDescriptor::calculate(const Molecule& mol) const {
     if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit Mol is null";

    auto bitInfo = getMorganFingerprintBitInfo(*rdkMol, radius, nBits);
    if (!bitInfo || bitInfo->empty()) return 0.0;

    double totalComplexity = 0;
    std::set<unsigned int> uniqueFragmentIDs; // Use bit IDs as proxy for fragment identifiers

    // Estimate complexity based on the radius of the environment generating the bit
    for(const auto& pair : *bitInfo) {
         uniqueFragmentIDs.insert(pair.first); // Use bit ID as fragment ID proxy
         const auto& atomRadiusVec = pair.second;
         if (!atomRadiusVec.empty()) {
             // Use max radius associated with this bit as complexity score
             unsigned int maxR = 0;
             for (const auto& arPair : atomRadiusVec) {
                 if (arPair.second > maxR) maxR = arPair.second;
             }
             totalComplexity += static_cast<double>(maxR + 1); // Include center atom (+1)
         }
    }


    if (uniqueFragmentIDs.empty()) return 0.0;
    // Average complexity over unique fragments (approximated by bit IDs)
    return totalComplexity / uniqueFragmentIDs.size();
}

// MorganInformationContentDensity
MorganInformationContentDensityDescriptor::MorganInformationContentDensityDescriptor()
    : Descriptor("MorganInfoDensity", "Information content (bit entropy) per set bit (r=2, n=1024)") {}

std::variant<double, int, std::string> MorganInformationContentDensityDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit Mol is null";

    auto fp = getMorganFingerprintBV(*rdkMol, radius, nBits);
    if (!fp) return "Error: FP calculation failed";

    long long numSetBits = fp->getNumOnBits();
    if (numSetBits == 0 || nBits == 0) return 0.0;

    double density = static_cast<double>(numSetBits) / nBits;
    if (density <= 1e-9 || density >= 1.0 - 1e-9) return 0.0; // Check bounds more carefully for log2

    // Shannon entropy of the bit vector (treating it as Bernoulli distribution)
    double bitEntropy = - (density * std::log2(density) + (1.0 - density) * std::log2(1.0 - density));
    double totalInformation = bitEntropy * nBits; // Total potential information

    // Density: information per set bit
    return totalInformation / numSetBits;
}

// MorganEnvironmentVariabilityIndex
MorganEnvironmentVariabilityIndexDescriptor::MorganEnvironmentVariabilityIndexDescriptor()
    : Descriptor("MorganEnvVariability", "Variability (avg Tanimoto dist) of atom environments (r=2, n=1024)") {}

std::variant<double, int, std::string> MorganEnvironmentVariabilityIndexDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol || rdkMol->getNumAtoms() < 2) return 0.0;

    std::vector<std::unique_ptr<ExplicitBitVect>> atomEnvFps;
    for (unsigned int i = 0; i < rdkMol->getNumAtoms(); ++i) {
        std::vector<std::uint32_t> atomIdVec = {i};
        try {
            ExplicitBitVect* fp = RDKit::MorganFingerprints::getFingerprintAsBitVect(
                *rdkMol, radius, nBits,
                &atomIdVec,  // fromAtoms
                nullptr,     // ignoreAtoms
                false,       // useChirality
                true,        // useBondTypes
                false,       // useFeatures
                nullptr,     // atomsSettingBits
                false        // includeRedundantEnvironments - MUST BE bool, not nullptr
            );
            if (fp) atomEnvFps.push_back(std::unique_ptr<ExplicitBitVect>(fp));
        } catch(...) { /* ignore atoms causing issues */ }
    }

    if (atomEnvFps.size() < 2) return 0.0;

    double totalDistance = 0;
    int pairCount = 0;
    for (size_t i = 0; i < atomEnvFps.size(); ++i) {
        for (size_t j = i + 1; j < atomEnvFps.size(); ++j) {
            // Explicitly check pointers before dereferencing
            if(atomEnvFps[i] && atomEnvFps[j]) {
                // Use corrected Tanimoto similarity calculation
                totalDistance += (1.0 - calculateTanimotoSimilarity(*atomEnvFps[i], *atomEnvFps[j]));
                pairCount++;
            }
        }
    }

    if (pairCount == 0) return 0.0;
    return totalDistance / pairCount; // Average Tanimoto distance
}

// MorganBitPositionImportanceScore
MorganBitPositionImportanceScoreDescriptor::MorganBitPositionImportanceScoreDescriptor()
    : Descriptor("MorganBitImportance", "Score based on bit importance (r=2, n=1024) (unimplemented)") {}

std::variant<double, int, std::string> MorganBitPositionImportanceScoreDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit Mol is null";
    
    auto fp = getMorganFingerprintBV(*rdkMol, radius, nBits);
    auto refFp = getReferenceFingerprintBV(radius, nBits);
    
    if (!fp || !refFp) return 0.0;
    
    // Define "important" bits as those set in the reference
    double importanceScore = 0.0;
    int refBitsOn = refFp->getNumOnBits();
    
    if (refBitsOn == 0) return 0.0;
    
    // Count bits that match the reference
    int matchingBits = 0;
    for (unsigned int i = 0; i < nBits; ++i) {
        if (refFp->getBit(i) && fp->getBit(i)) {
            matchingBits++;
        }
    }
    
    // Return proportion of reference bits captured
    return static_cast<double>(matchingBits) / refBitsOn;
}

// MorganRingSystemRepresentationScore
MorganRingSystemRepresentationScoreDescriptor::MorganRingSystemRepresentationScoreDescriptor()
    : Descriptor("MorganRingSysScore", "Score based on ring system representation (r=2, n=1024)") {}

std::variant<double, int, std::string> MorganRingSystemRepresentationScoreDescriptor::calculate(const Molecule& mol) const {
     if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit Mol is null";
    RDKit::MolOps::symmetrizeSSSR(*const_cast<RDKit::RWMol*>(static_cast<const RDKit::RWMol*>(rdkMol))); // Ensure SSSR
    const RDKit::RingInfo* ringInfo = rdkMol->getRingInfo();
    if (!ringInfo || !ringInfo->isInitialized() || ringInfo->numRings() == 0) return 0.0;


    auto bitInfo = getMorganFingerprintBitInfo(*rdkMol, radius, nBits);
    if (!bitInfo) return "Error: FP calculation failed";

    int ringAtomBits = 0;
    std::set<unsigned int> uniqueBitsConsidered;

    for (const auto& pair : *bitInfo) {
        // We use the *unfolded* bit ID from bitInfo map directly as fragment identifier proxy
        unsigned int unfoldedBitId = pair.first;
        unsigned int bitId = unfoldedBitId % nBits; // Calculate folded bit ID for checking if set later

        // Check if the *folded* bit is actually set in the final fingerprint
        // Need the BV for this check
        // Let's recalculate BV here - less efficient but ensures correctness
        auto fp_check = getMorganFingerprintBV(*rdkMol, radius, nBits);
        if (!fp_check || !fp_check->getBit(bitId)) continue; // Skip if the folded bit isn't set


         if (uniqueBitsConsidered.count(bitId)) continue; // Only count each final *set* bit once

         const auto& atomRadiusVec = pair.second;
         if (!atomRadiusVec.empty()) {
             unsigned int centerAtomIdx = atomRadiusVec[0].first;
             // Check if the center atom of the fragment is in a ring
             if (rdkMol->getRingInfo()->numAtomRings(centerAtomIdx) > 0) {
                 ringAtomBits++;
                 uniqueBitsConsidered.insert(bitId);
             }
         }
    }

     long long totalSetBits = uniqueBitsConsidered.size(); // Count the bits we actually considered
     if (totalSetBits == 0) return 0.0;

     // Return fraction of set bits originating from ring atoms
     return static_cast<double>(ringAtomBits) / totalSetBits;
}

// MorganBitProbabilityDistributionEntropy
MorganBitProbabilityDistributionEntropyDescriptor::MorganBitProbabilityDistributionEntropyDescriptor()
    : Descriptor("MorganBitProbEntropy", "Entropy of bit probability distribution (r=2, n=1024) (requires dataset - unimplemented)") {}

std::variant<double, int, std::string> MorganBitProbabilityDistributionEntropyDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit Mol is null";
    
    auto fp = getMorganFingerprintBV(*rdkMol, radius, nBits);
    auto refFp = getReferenceFingerprintBV(radius, nBits);
    
    if (!fp || !refFp) return 0.0;
    
    // Calculate entropy using reference as the "expected" probability
    // Create probability distribution from fingerprint bit frequencies
    std::map<int, double> bitProbs; // 0, 1, 2, 3 for bit patterns: (00,01,10,11)
    
    double totalCount = 0;
    
    for (unsigned int i = 0; i < nBits; i++) {
        bool molBit = fp->getBit(i);
        bool refBit = refFp->getBit(i);
        
        int pattern = (molBit ? 2 : 0) + (refBit ? 1 : 0); // 00=0, 01=1, 10=2, 11=3
        bitProbs[pattern] += 1.0;
        totalCount += 1.0;
    }
    
    // Calculate Shannon entropy
    double entropy = 0.0;
    for (const auto& pair : bitProbs) {
        double prob = pair.second / totalCount;
        if (prob > 0) {
            entropy -= prob * std::log2(prob);
        }
    }
    
    return entropy;
}

// MorganFragmentElectronegativitySpectrum
MorganFragmentElectronegativitySpectrumDescriptor::MorganFragmentElectronegativitySpectrumDescriptor()
    : Descriptor("MorganFragENSpectrum", "Avg EN of atoms in fragments (r=2, n=1024)") {}

std::variant<double, int, std::string> MorganFragmentElectronegativitySpectrumDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit Mol is null";

    auto bitInfo = getMorganFingerprintBitInfo(*rdkMol, radius, nBits);
     if (!bitInfo || bitInfo->empty()) return 0.0;

     double totalAvgEN = 0;
     int fragmentCount = 0;
     std::set<unsigned int> processedFragmentIDs; // Process each fragment ID once

     for (const auto& pair : *bitInfo) {
         unsigned int fragID = pair.first; // Use unfolded bit ID as fragment proxy
         if (processedFragmentIDs.count(fragID)) continue;

         const auto& atomRadiusVec = pair.second;
         if (atomRadiusVec.empty()) continue;

         std::set<unsigned int> atomsInFragment;
         unsigned int centerAtomIdx = atomRadiusVec[0].first;
         unsigned int maxRadius = atomRadiusVec[0].second;

         // BFS from center up to maxRadius to get atoms
         std::queue<std::pair<unsigned int, unsigned int>> q;
         q.push({centerAtomIdx, 0});
         std::set<unsigned int> visitedInFrag;
         visitedInFrag.insert(centerAtomIdx);
         atomsInFragment.insert(centerAtomIdx);

         while(!q.empty()) {
             auto current = q.front();
             q.pop();
             unsigned int currentAtomIdx = current.first;
             unsigned int currentR = current.second;

             if (currentR < maxRadius) {
                 for(const auto& nbr : rdkMol->atomNeighbors(rdkMol->getAtomWithIdx(currentAtomIdx))) {
                     unsigned int nbrIdx = nbr->getIdx();
                     if(visitedInFrag.find(nbrIdx) == visitedInFrag.end()) {
                         visitedInFrag.insert(nbrIdx);
                         atomsInFragment.insert(nbrIdx);
                         q.push({nbrIdx, currentR + 1});
                     }
                 }
             }
         }


         if (!atomsInFragment.empty()) {
             double currentFragENSum = 0;
             for (unsigned int atomIdx : atomsInFragment) {
                 // Assuming SumDescriptor is accessible or EN helper is local
                 currentFragENSum += SumDescriptor::getAtomElectronegativity(rdkMol->getAtomWithIdx(atomIdx));
             }
             totalAvgEN += currentFragENSum / atomsInFragment.size();
             fragmentCount++;
             processedFragmentIDs.insert(fragID);
         }
     }

     if (fragmentCount == 0) return 0.0;
     return totalAvgEN / fragmentCount; // Average of average fragment ENs
}

// MorganBitPositionCorrelationMatrix
MorganBitPositionCorrelationMatrixDescriptor::MorganBitPositionCorrelationMatrixDescriptor()
    : Descriptor("MorganBitCorrMatrix", "Correlation matrix of bits (r=2, n=1024) (requires dataset - unimplemented)") {}

std::variant<double, int, std::string> MorganBitPositionCorrelationMatrixDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit Mol is null";
    
    auto fp = getMorganFingerprintBV(*rdkMol, radius, nBits);
    auto refFp = getReferenceFingerprintBV(radius, nBits);
    
    if (!fp || !refFp) return 0.0;
    
    // Compute average correlation across neighboring bit positions
    double avgCorrelation = 0.0;
    int count = 0;
    
    for (unsigned int i = 0; i < nBits - 1; ++i) {
        bool pos1Mol = fp->getBit(i);
        bool pos1Ref = refFp->getBit(i);
        bool pos2Mol = fp->getBit(i+1);
        bool pos2Ref = refFp->getBit(i+1);
        
        // Simple correlation: 1 if both match, -1 if opposite, 0 if mixed
        double correlation = 0.0;
        if ((pos1Mol == pos1Ref) && (pos2Mol == pos2Ref)) {
            correlation = 1.0;
        } else if ((pos1Mol != pos1Ref) && (pos2Mol != pos2Ref)) {
            correlation = 1.0;
        } else {
            correlation = -1.0;
        }
        
        avgCorrelation += correlation;
        count++;
    }
    
    if (count == 0) return 0.0;
    
    return avgCorrelation / count;
}

// MorganFragmentConnectivityPattern
MorganFragmentConnectivityPatternDescriptor::MorganFragmentConnectivityPatternDescriptor()
    : Descriptor("MorganFragConnectivity", "Connectivity pattern of fragments (r=2, n=1024) (unimplemented)") {}

std::variant<double, int, std::string> MorganFragmentConnectivityPatternDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit Mol is null";
    
    auto fp = getMorganFingerprintBV(*rdkMol, radius, nBits);
    auto refFp = getReferenceFingerprintBV(radius, nBits);
    
    if (!fp || !refFp) return 0.0;
    
    // Measure connectivity through bit transition patterns
    int molTransitions = 0, refTransitions = 0, commonTransitions = 0;
    
    for (unsigned int i = 0; i < nBits - 1; ++i) {
        bool molTransition = fp->getBit(i) != fp->getBit(i+1);
        bool refTransition = refFp->getBit(i) != refFp->getBit(i+1);
        
        if (molTransition) molTransitions++;
        if (refTransition) refTransitions++;
        if (molTransition && refTransition) commonTransitions++;
    }
    
    if (refTransitions == 0) return 0.0;
    
    // Return similarity of transition patterns to reference
    return static_cast<double>(commonTransitions) / refTransitions;
}

// MorganBitOccurrenceFrequencySkewness
MorganBitOccurrenceFrequencySkewnessDescriptor::MorganBitOccurrenceFrequencySkewnessDescriptor()
    : Descriptor("MorganBitFreqSkew", "Skewness of bit frequency (r=2, n=1024) (requires dataset - unimplemented)") {}

std::variant<double, int, std::string> MorganBitOccurrenceFrequencySkewnessDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit Mol is null";
    
    auto fp = getMorganFingerprintBV(*rdkMol, radius, nBits);
    auto refFp = getReferenceFingerprintBV(radius, nBits);
    
    if (!fp || !refFp) return 0.0;
    
    // Calculate bit frequency difference normalized against reference
    std::vector<double> relativeFrequencies;
    
    // Count bits in 8 evenly-spaced regions
    int numRegions = 8;
    std::vector<int> molCounts(numRegions, 0);
    std::vector<int> refCounts(numRegions, 0);
    
    int bitsPerRegion = nBits / numRegions;
    
    for (int r = 0; r < numRegions; r++) {
        int start = r * bitsPerRegion;
        int end = (r == numRegions - 1) ? nBits : (r + 1) * bitsPerRegion;
        
        for (int i = start; i < end; i++) {
            if (fp->getBit(i)) molCounts[r]++;
            if (refFp->getBit(i)) refCounts[r]++;
        }
        
        // Avoid division by zero
        if (refCounts[r] > 0) {
            relativeFrequencies.push_back(static_cast<double>(molCounts[r]) / refCounts[r]);
        } else if (molCounts[r] == 0) {
            relativeFrequencies.push_back(1.0); // Both 0, consider equal
        } else {
            relativeFrequencies.push_back(2.0); // Infinite relative frequency, cap at 2
        }
    }
    
    // Calculate skewness of relative frequencies
    if (relativeFrequencies.size() < 3) return 0.0;
    
    return calculateSkewnessInternal(relativeFrequencies);
}

// MorganSubstructureHeterogeneityIndex
MorganSubstructureHeterogeneityIndexDescriptor::MorganSubstructureHeterogeneityIndexDescriptor()
    : Descriptor("MorganSubstructHeterogeneity", "Heterogeneity (entropy) of fragment counts (r=2)") {}

std::variant<double, int, std::string> MorganSubstructureHeterogeneityIndexDescriptor::calculate(const Molecule& mol) const {
     if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit Mol is null";

    // Same logic as MorganFragmentDiversityScore (using Shannon entropy)
    auto counts = getMorganFingerprintCounts(*rdkMol, radius);
    if (!counts) return "Error: FP calculation failed";

    return calculateShannonEntropyInternal(*counts);
}

// MorganBitPolarityDistribution
MorganBitPolarityDistributionDescriptor::MorganBitPolarityDistributionDescriptor()
    : Descriptor("MorganBitPolarityDist", "Distribution of polar bits (r=2, n=1024)") {}

std::variant<double, int, std::string> MorganBitPolarityDistributionDescriptor::calculate(const Molecule& mol) const {
     if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit Mol is null";

    auto bitInfo = getMorganFingerprintBitInfo(*rdkMol, radius, nBits);
     if (!bitInfo) return "Error: FP calculation failed";

     // Need Gasteiger charges to assess polarity
     std::vector<double> charges;
     if (!getGasteigerCharges(rdkMol, charges)) {
         globalLogger.debug("Failed to get Gasteiger charges for MorganBitPolarityDist.");
         return "Error: Gasteiger failed";
     }

     int polarBits = 0;
     std::set<unsigned int> uniqueBitsConsidered;
     double chargeThreshold = 0.1; // Threshold to consider an atom polar

     for (const auto& pair : *bitInfo) {
        unsigned int bitId = pair.first % nBits;
         // Check if the *folded* bit is actually set in the final fingerprint
        auto fp_check = getMorganFingerprintBV(*rdkMol, radius, nBits); // Recheck needed here
        if (!fp_check || !fp_check->getBit(bitId)) continue;

        if (uniqueBitsConsidered.count(bitId)) continue;

        const auto& atomRadiusVec = pair.second;
        if (!atomRadiusVec.empty()) {
            unsigned int centerAtomIdx = atomRadiusVec[0].first;
             // Check if the center atom has significant partial charge
            if (centerAtomIdx < charges.size() && std::abs(charges[centerAtomIdx]) > chargeThreshold) {
                 polarBits++;
                 uniqueBitsConsidered.insert(bitId);
            }
             // Alternative: Check if *any* atom in the fragment is polar? More complex.
        }
     }

     long long totalSetBits = uniqueBitsConsidered.size();
     if (totalSetBits == 0) return 0.0;

     // Return fraction of set bits deemed polar (based on center atom charge)
     return static_cast<double>(polarBits) / totalSetBits;
}

// MorganFragmentSimilarityNetwork
MorganFragmentSimilarityNetworkDescriptor::MorganFragmentSimilarityNetworkDescriptor()
    : Descriptor("MorganFragSimNetwork", "Network properties of fragments (r=2, n=1024) (unimplemented)") {}

std::variant<double, int, std::string> MorganFragmentSimilarityNetworkDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit Mol is null";
    
    // Compare fingerprint networks by examining bit patterns
    auto fp = getMorganFingerprintBV(*rdkMol, radius, nBits);
    auto refFp = getReferenceFingerprintBV(radius, nBits);
    
    if (!fp || !refFp) return 0.0;
    
    // Count clusters (adjacent set bits) in both fingerprints
    int molClusters = 0, refClusters = 0;
    bool inMolCluster = false, inRefCluster = false;
    
    for (unsigned int i = 0; i < nBits; ++i) {
        bool molBit = fp->getBit(i);
        bool refBit = refFp->getBit(i);
        
        // Track molecule clusters
        if (molBit && !inMolCluster) {
            inMolCluster = true;
            molClusters++;
        } else if (!molBit) {
            inMolCluster = false;
        }
        
        // Track reference clusters
        if (refBit && !inRefCluster) {
            inRefCluster = true;
            refClusters++;
        } else if (!refBit) {
            inRefCluster = false;
        }
    }
    
    if (refClusters == 0) return 1.0; // No reference clusters
    
    // Return similarity based on cluster count comparison
    return 1.0 - std::abs(static_cast<double>(molClusters - refClusters)) / refClusters;
}

// MorganBitPositionInformationGain
MorganBitPositionInformationGainDescriptor::MorganBitPositionInformationGainDescriptor()
    : Descriptor("MorganBitInfoGain", "Information gain from bits (r=2, n=1024) (requires context - unimplemented)") {}

std::variant<double, int, std::string> MorganBitPositionInformationGainDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit Mol is null";
    
    auto fp = getMorganFingerprintBV(*rdkMol, radius, nBits);
    auto refFp = getReferenceFingerprintBV(radius, nBits);
    
    if (!fp || !refFp) return 0.0;
    
    // Calculate information gain based on KL divergence to reference
    double infoGain = 0.0;
    
    // Count number of set bits in different regions of the fingerprint
    constexpr int numRegions = 16;
    const int bitsPerRegion = nBits / numRegions;
    
    for (int region = 0; region < numRegions; region++) {
        int start = region * bitsPerRegion;
        int end = (region == numRegions - 1) ? nBits : (region + 1) * bitsPerRegion;
        
        double molProb = 0.0;
        double refProb = 0.0;
        
        for (int i = start; i < end; i++) {
            if (fp->getBit(i)) molProb += 1.0;
            if (refFp->getBit(i)) refProb += 1.0;
        }
        
        // Normalize to probabilities
        molProb = molProb / (end - start) + 1e-6; // Add small constant to avoid log(0)
        refProb = refProb / (end - start) + 1e-6;
        
        // KL divergence term: p * log(p/q)
        infoGain += molProb * std::log2(molProb / refProb);
    }
    
    return std::abs(infoGain); // Return absolute value of information gain
}

// MorganSubstructureDiversityGradient
MorganSubstructureDiversityGradientDescriptor::MorganSubstructureDiversityGradientDescriptor()
    : Descriptor("MorganSubstructDivGradient", "Gradient of fragment diversity r=2 vs r=1") {}

std::variant<double, int, std::string> MorganSubstructureDiversityGradientDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit Mol is null";

    auto counts1 = getMorganFingerprintCounts(*rdkMol, radius1);
    auto counts2 = getMorganFingerprintCounts(*rdkMol, radius2);

    if (!counts1 || !counts2) return "Error: FP calculation failed";

    double entropy1 = calculateShannonEntropyInternal(*counts1);
    double entropy2 = calculateShannonEntropyInternal(*counts2);

    // Gradient = Difference in entropy (normalized by radius difference?)
    // Simple difference: entropy2 - entropy1
    return entropy2 - entropy1;
}

// MorganFingerprintSymmetryScore
MorganFingerprintSymmetryScoreDescriptor::MorganFingerprintSymmetryScoreDescriptor()
    : Descriptor("MorganFPSymmetry", "Symmetry score of Morgan fingerprint (r=2, n=1024)") {}

std::variant<double, int, std::string> MorganFingerprintSymmetryScoreDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return "Error: RDKit Mol is null";

    auto fp = getMorganFingerprintBV(*rdkMol, radius, nBits);
    if (!fp) return "Error: FP calculation failed";

    int matchingBits = 0;
    int comparableBits = nBits / 2; // Number of pairs to compare

    for (unsigned int i = 0; i < comparableBits; ++i) {
        if (fp->getBit(i) == fp->getBit(nBits - 1 - i)) {
            matchingBits++;
        }
    }

    if (comparableBits == 0) return 1.0; // Perfectly symmetric if no bits to compare

    // Return fraction of matching symmetric pairs
    return static_cast<double>(matchingBits) / comparableBits;
}

} // namespace descriptors
} // namespace desfact