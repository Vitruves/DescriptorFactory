#pragma once

#include "descriptors.hpp"
#include "sum.hpp"
#include "fractional.hpp"

namespace desfact {
namespace descriptors {

// Forward declarations if needed internally

// 1. Structural and Connectivity-Based
class AtomicConnectivityImbalanceDescriptor : public SumDescriptor {
public:
    AtomicConnectivityImbalanceDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class RingBridgeRatioDescriptor : public FractionalDescriptor {
public:
    RingBridgeRatioDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class SubstitutionPatternComplexityDescriptor : public Descriptor {
public:
    SubstitutionPatternComplexityDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class OpenChainSaturationRatioDescriptor : public FractionalDescriptor {
public:
    OpenChainSaturationRatioDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

// 2. Atomic Neighborhood Patterns
class HeteroatomNeighborhoodDiversityDescriptor : public SumDescriptor {
public:
    HeteroatomNeighborhoodDiversityDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class CarbonNeighborhoodUniformityDescriptor : public FractionalDescriptor {
public:
    CarbonNeighborhoodUniformityDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class PolarAtomNeighborhoodRatioDescriptor : public FractionalDescriptor {
public:
    PolarAtomNeighborhoodRatioDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class AverageAtomDegreeRangeDescriptor : public Descriptor {
public:
    AverageAtomDegreeRangeDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

// 3. Ring System Specific
class RingJunctionComplexityDescriptor : public Descriptor {
public:
    RingJunctionComplexityDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class NonFusedRingDensityDescriptor : public FractionalDescriptor {
public:
    NonFusedRingDensityDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class RingChainAttachmentDensityDescriptor : public FractionalDescriptor {
public:
    RingChainAttachmentDensityDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class RingTerminalSubstituentRatioDescriptor : public FractionalDescriptor {
public:
    RingTerminalSubstituentRatioDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class RingSaturationBalanceDescriptor : public FractionalDescriptor {
public:
    RingSaturationBalanceDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

// 4. Electronic Influence
class PolarizabilityGradientDescriptor : public SumDescriptor {
public:
    PolarizabilityGradientDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class ElectronWithdrawingAtomDensityDescriptor : public FractionalDescriptor {
public:
    ElectronWithdrawingAtomDensityDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class ElectronDonatingAtomDensityDescriptor : public FractionalDescriptor {
public:
    ElectronDonatingAtomDensityDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class ElectronegativityGradientDensityDescriptor : public FractionalDescriptor {
public:
    ElectronegativityGradientDensityDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class PeripheralElectronRichAtomRatioDescriptor : public FractionalDescriptor {
public:
    PeripheralElectronRichAtomRatioDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

// 5. Functional Group and Substitution Patterns
class FunctionalGroupIsolationIndexDescriptor : public FractionalDescriptor {
public:
    FunctionalGroupIsolationIndexDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class HydroxylGroupDispersionDescriptor : public Descriptor {
public:
    HydroxylGroupDispersionDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class AlkylChainDiversityDescriptor : public Descriptor {
public:
    AlkylChainDiversityDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class SubstituentPositionVarianceDescriptor : public Descriptor {
public:
    SubstituentPositionVarianceDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class TerminalFunctionalGroupClusteringDescriptor : public FractionalDescriptor {
public:
    TerminalFunctionalGroupClusteringDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

// 6. Bond-type and Hybridization
class ChainSaturationVariabilityDescriptor : public Descriptor {
public:
    ChainSaturationVariabilityDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class TripleBondTerminalRatioDescriptor : public FractionalDescriptor {
public:
    TripleBondTerminalRatioDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class AdjacentHybridizationTransitionDescriptor : public FractionalDescriptor {
public:
    AdjacentHybridizationTransitionDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class CyclicHybridizationHomogeneityDescriptor : public SumDescriptor {
public:
    CyclicHybridizationHomogeneityDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class InternalChainUnsaturationDensityDescriptor : public FractionalDescriptor {
public:
    InternalChainUnsaturationDensityDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

// 7. Hydrogen Bonding Patterns
class HydrogenBondDonorClusteringDescriptor : public FractionalDescriptor {
public:
    HydrogenBondDonorClusteringDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class AcceptorDonorRatioImbalanceDescriptor : public FractionalDescriptor {
public:
    AcceptorDonorRatioImbalanceDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class PeripheralDonorAcceptorBalanceDescriptor : public FractionalDescriptor {
public:
    PeripheralDonorAcceptorBalanceDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class IntraringHBondPotentialDescriptor : public FractionalDescriptor {
public:
    IntraringHBondPotentialDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class ChainEndHydrogenBondDensityDescriptor : public FractionalDescriptor {
public:
    ChainEndHydrogenBondDensityDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

// 8. Formal Charge Distribution
class FormalChargeNeighborhoodVarianceDescriptor : public SumDescriptor {
public:
    FormalChargeNeighborhoodVarianceDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class OppositeChargeNeighborRatioDescriptor : public FractionalDescriptor {
public:
    OppositeChargeNeighborRatioDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class ChargeGradientAlongChainsDescriptor : public Descriptor {
public:
    ChargeGradientAlongChainsDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class PeripheralChargeNeutralityDescriptor : public FractionalDescriptor {
public:
    PeripheralChargeNeutralityDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

// 9. Aromaticity & Conjugation Patterns
class InterRingConjugationRatioDescriptor : public FractionalDescriptor {
public:
    InterRingConjugationRatioDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class AromaticChainDensityDescriptor : public FractionalDescriptor {
public:
    AromaticChainDensityDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class ConjugationLengthVarianceDescriptor : public Descriptor {
public:
    ConjugationLengthVarianceDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class TerminalAromaticSubstitutionDescriptor : public FractionalDescriptor {
public:
    TerminalAromaticSubstitutionDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class CyclicVsChainAromaticRatioDescriptor : public FractionalDescriptor {
public:
    CyclicVsChainAromaticRatioDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

// 10. Structural Diversity and Complexity
class UniqueElementPairRatioDescriptor : public FractionalDescriptor {
public:
    UniqueElementPairRatioDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class HeavyAtomDegreeDiversityDescriptor : public Descriptor {
public:
    HeavyAtomDegreeDiversityDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class AtomPathDiversityDescriptor : public Descriptor {
public:
    AtomPathDiversityDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class InternalAtomComplexityDescriptor : public Descriptor {
public:
    InternalAtomComplexityDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class HeteroatomBondOrderVariabilityDescriptor : public SumDescriptor {
public:
    HeteroatomBondOrderVariabilityDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

} // namespace descriptors
} // namespace desfact
