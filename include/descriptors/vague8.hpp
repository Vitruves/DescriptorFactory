#pragma once

#include "descriptors/sum.hpp"
#include "descriptors/fractional.hpp"
#include "descriptors/strings.hpp" // Include if any descriptor needs string manipulation base
#include <GraphMol/GraphMol.h>
#include <GraphMol/Atom.h>
#include <GraphMol/Bond.h>

namespace desfact {
namespace descriptors {

// --- Vague8 Descriptors ---

// TopologicalChargeDistributionSkewness
class TopologicalChargeDistributionSkewnessDescriptor : public Descriptor {
public:
    TopologicalChargeDistributionSkewnessDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

// ElementNeighborhoodDiversity
class ElementNeighborhoodDiversityDescriptor : public Descriptor {
public:
    ElementNeighborhoodDiversityDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

// BondOrderAlternationPattern
class BondOrderAlternationPatternDescriptor : public Descriptor {
public:
    BondOrderAlternationPatternDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

// FunctionalGroupConnectivityMatrix -> AvgShortestPathBetweenDiffFuncGroups
class AvgShortestPathBetweenDiffFuncGroupsDescriptor : public Descriptor {
public:
    AvgShortestPathBetweenDiffFuncGroupsDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

// HeteroatomClusterTopology -> AvgHeteroClusterSize
class AvgHeteroClusterSizeDescriptor : public Descriptor {
public:
    AvgHeteroClusterSizeDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

// RingSubstitutionPatternCode -> EntropyOfRingSubstCountDist
class EntropyOfRingSubstCountDistDescriptor : public Descriptor {
public:
    EntropyOfRingSubstCountDistDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};


// FunctionalGroupDistanceHistogram -> EntropyOfFuncGroupDistances
class EntropyOfFuncGroupDistancesDescriptor : public Descriptor {
public:
    EntropyOfFuncGroupDistancesDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

// ChainBranchingRecursivePattern -> BranchPointsInChains / ChainAtoms
class ChainBranchingPatternDescriptor : public FractionalDescriptor { // Fractional seems appropriate
public:
    ChainBranchingPatternDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

// ElectronegativeAtomNeighborhoodScore -> Sum(CountENNeighborsRadius2)
class ElectronegativeAtomNeighborhoodScoreDescriptor : public SumDescriptor { // Sum seems appropriate
public:
    ElectronegativeAtomNeighborhoodScoreDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

// TopologicalChargeSeparationIndex -> MeanDistBetweenFormalCharges
class TopologicalChargeSeparationIndexDescriptor : public Descriptor {
public:
    TopologicalChargeSeparationIndexDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

// RingSystemConnectivityPattern -> InterRingSystemBonds / NumRingSystems
class RingSystemConnectivityPatternDescriptor : public Descriptor {
public:
    RingSystemConnectivityPatternDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

// HeteroatomPositionEntropy -> EntropyOfHeteroatomIndicesInCanonSmiles
class HeteroatomPositionEntropyDescriptor : public StringDescriptor { // StringDescriptor is suitable
public:
    HeteroatomPositionEntropyDescriptor();
    // Override calculateFromSmiles instead of calculate
    std::variant<double, int, std::string> calculateFromSmiles(const std::string& smiles) const override;
private:
    // calculate is implemented by the base StringDescriptor
    using StringDescriptor::calculate;
};

// BondOrderTransitionFrequency -> EntropyOfAdjacentBondTypePairs
class BondOrderTransitionFrequencyDescriptor : public Descriptor {
public:
    BondOrderTransitionFrequencyDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

// AtomNeighborhoodElectronegativityGradient -> Avg(AtomEN - AvgNeighborEN)
class AtomNeighborhoodElectronegativityGradientDescriptor : public SumDescriptor { // Sum/Avg pattern fits SumDescriptor
public:
    AtomNeighborhoodElectronegativityGradientDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

// ChainLengthDistributionEntropy
class ChainLengthDistributionEntropyDescriptor : public Descriptor {
public:
    ChainLengthDistributionEntropyDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

// RingSubstitutionSymmetry -> Avg(DistinctSubstituentTypesPerRing)
class RingSubstitutionSymmetryDescriptor : public Descriptor {
public:
    RingSubstitutionSymmetryDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

// HeteroatomSequencePatterns -> EntropyOfHeteroatomTriplets (A-B-C)
class HeteroatomSequencePatternsDescriptor : public Descriptor {
public:
    HeteroatomSequencePatternsDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

// FunctionalGroupIsolationTopology -> Avg(MinDistToOtherFuncGroup)
class FunctionalGroupIsolationTopologyDescriptor : public Descriptor {
public:
    FunctionalGroupIsolationTopologyDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

// HydrophobicPatchConnectivity -> NumHydrophobicPatches / NumHydrophobicAtoms
class HydrophobicPatchConnectivityDescriptor : public FractionalDescriptor { // Fractional fits
public:
    HydrophobicPatchConnectivityDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

// BondTopologicalEnvironmentFingerprint -> CountUniqueBondEnvFingerprints
class BondTopologicalEnvironmentFingerprintDescriptor : public Descriptor {
public:
    BondTopologicalEnvironmentFingerprintDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

// ChiralCenterTopologicalDistribution -> VarianceOfDistancesBetweenChiralCenters
class ChiralCenterTopologicalDistributionDescriptor : public Descriptor {
public:
    ChiralCenterTopologicalDistributionDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

// RingFusionPatternCode -> NumSharedRingBonds / NumRings
class RingFusionPatternCodeDescriptor : public Descriptor {
public:
    RingFusionPatternCodeDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

// ElectronegativityTopologicalMoment -> Sum(EN_i * dist(i, 0)) / Sum(EN_i)
class ElectronegativityTopologicalMomentDescriptor : public SumDescriptor { // Sum fits
public:
    ElectronegativityTopologicalMomentDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

// AtomicRadiiVarianceInNeighborhoods -> Avg(Variance(NeighborRadii))
class AtomicRadiiVarianceInNeighborhoodsDescriptor : public SumDescriptor { // Sum fits
public:
    AtomicRadiiVarianceInNeighborhoodsDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

// HeteroatomBondingPatternCode -> CountUniqueHeteroatomBondingCodes
class HeteroatomBondingPatternCodeDescriptor : public Descriptor {
public:
    HeteroatomBondingPatternCodeDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

// SubstructureFrequencySpectrum -> EntropyOfSubstructureFrequencies
class SubstructureFrequencySpectrumDescriptor : public Descriptor {
public:
    SubstructureFrequencySpectrumDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

// FormalChargeNeighborhoodPattern -> CountUniqueFormalChargeNeighborPatterns(R=1)
class FormalChargeNeighborhoodPatternDescriptor : public Descriptor {
public:
    FormalChargeNeighborhoodPatternDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};


} // namespace descriptors
} // namespace desfact
