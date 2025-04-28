#pragma once

#include "descriptors.hpp" // Base class
#include <string>
#include <variant>

// Forward declare Molecule in the global ::desfact namespace, not locally
namespace desfact { class Molecule; }

namespace desfact {
namespace descriptors {

// --- Topological and Structural ---
class LongestContinuousSp2FragmentLength : public Descriptor {
public:
    LongestContinuousSp2FragmentLength();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

class ChainBranchDepthMaximum : public Descriptor {
public:
    ChainBranchDepthMaximum();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

class ChainBranchingFrequency : public Descriptor {
public:
    ChainBranchingFrequency();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

class LongestAlternatingSingleDoubleBondPath : public Descriptor {
public:
    LongestAlternatingSingleDoubleBondPath();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

class AverageRingStrainProxy : public Descriptor {
public:
    AverageRingStrainProxy();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

class FractionBackboneCarbon : public Descriptor {
public:
    FractionBackboneCarbon();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

class TopologicalAsymmetryIndex : public Descriptor {
public:
    TopologicalAsymmetryIndex();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

class FractionConjugatedBondsBackbone : public Descriptor {
public:
    FractionConjugatedBondsBackbone();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

class MaxChainNoHetero : public Descriptor {
public:
    MaxChainNoHetero();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

class MeanRingSizeWeightedDegree : public Descriptor {
public:
    MeanRingSizeWeightedDegree();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

// --- Electronic / Electronegativity Pattern ---
class PolarClusterCount : public Descriptor {
public:
    PolarClusterCount();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

class MaxElectronegativityGradientChain : public Descriptor {
public:
    MaxElectronegativityGradientChain();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

class FractionPolarizableHeavy : public Descriptor {
public:
    FractionPolarizableHeavy();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

class NormalizedNetDipoleHeuristic : public Descriptor {
public:
    NormalizedNetDipoleHeuristic();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

class CumulativeENDifferenceBackbone : public Descriptor {
public:
    CumulativeENDifferenceBackbone();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

class HeavyAtomLocalChargeSkewness : public Descriptor {
public:
    HeavyAtomLocalChargeSkewness();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

class ElectronegativeEndChainCount : public Descriptor {
public:
    ElectronegativeEndChainCount();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

class FractionBondsLargeENGap : public Descriptor {
public:
    FractionBondsLargeENGap();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

class SymmetryElectronegativeDistribution : public Descriptor {
public:
    SymmetryElectronegativeDistribution();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

class FractionHeteroatomsTerminal : public Descriptor {
public:
    FractionHeteroatomsTerminal();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

// --- Conjugation / Resonance-Related ---
class ConjugatedSystemSizeMax : public Descriptor {
public:
    ConjugatedSystemSizeMax();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

class FractionAromaticAtomsConnectedChains : public Descriptor {
public:
    FractionAromaticAtomsConnectedChains();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

class MeanDistanceAromaticSystems : public Descriptor {
public:
    MeanDistanceAromaticSystems();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

class LengthLargestFullyConjugatedRingSystem : public Descriptor {
public:
    LengthLargestFullyConjugatedRingSystem();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

class RatioConjugatedCarbonsTotalCarbons : public Descriptor {
public:
    RatioConjugatedCarbonsTotalCarbons();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

class NormalizedConjugationPathwayDensity : public Descriptor {
public:
    NormalizedConjugationPathwayDensity();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

// --- Geometry Approximation / Exposure / Shape ---
class EstimatedMolecularAspectRatio : public Descriptor {
public:
    EstimatedMolecularAspectRatio();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

class TerminalHeavyAtomDispersion : public Descriptor {
public:
    TerminalHeavyAtomDispersion();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

class ConnectivityCorePeripheryGradient : public Descriptor {
public:
    ConnectivityCorePeripheryGradient();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

class HeavyAtomPlanarityHeuristic : public Descriptor {
public:
    HeavyAtomPlanarityHeuristic();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

class RingCoreChainPeripherySizeRatio : public Descriptor {
public:
    RingCoreChainPeripherySizeRatio();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

class RingChainAlternationFrequency : public Descriptor {
public:
    RingChainAlternationFrequency();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

// --- Charge and Ionization Site Heuristics ---
class PotentialDeprotonationSiteCount : public Descriptor {
public:
    PotentialDeprotonationSiteCount();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

class PotentialProtonationSiteCount : public Descriptor {
public:
    PotentialProtonationSiteCount();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

class FractionNonTerminalChargedSites : public Descriptor {
public:
    FractionNonTerminalChargedSites();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

class MeanTopologicalDistanceIonizableSites : public Descriptor {
public:
    MeanTopologicalDistanceIonizableSites();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

class RingIonizableSiteDensity : public Descriptor {
public:
    RingIonizableSiteDensity();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

class AromaticProtonDonorAcceptorRatioHeuristic : public Descriptor {
public:
    AromaticProtonDonorAcceptorRatioHeuristic();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

// --- Hydrophobicity / Polarity Features ---
class HydrophilicAtomClusterSizeMean : public Descriptor {
public:
    HydrophilicAtomClusterSizeMean();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

class HydrophobicIslandSizeMax : public Descriptor {
public:
    HydrophobicIslandSizeMax();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

class HydrophobicPolarSurfaceInterfaceEstimation : public Descriptor {
public:
    HydrophobicPolarSurfaceInterfaceEstimation();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

class NormalizedHydrophobicPolarRatioEstimate : public Descriptor {
public:
    NormalizedHydrophobicPolarRatioEstimate();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

class HydrophobicAtomPathLengthMean : public Descriptor {
public:
    HydrophobicAtomPathLengthMean();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

// Add 5 remaining descriptors to reach 50
class HeteroatomConnectivityIndex : public Descriptor {
public:
    HeteroatomConnectivityIndex();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

class RingFusionDensity : public Descriptor {
public:
    RingFusionDensity();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

class AverageBondPolarity : public Descriptor {
public:
    AverageBondPolarity();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

class StericHindranceProxy : public Descriptor {
public:
    StericHindranceProxy();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

class MolecularFlexibilityProxy : public Descriptor {
public:
    MolecularFlexibilityProxy();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};


} // namespace descriptors
} // namespace desfact
