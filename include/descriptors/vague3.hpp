#pragma once

#include "descriptors/sum.hpp"
#include "descriptors/fractional.hpp"
#include <GraphMol/GraphMol.h>
#include <GraphMol/Atom.h>
#include <GraphMol/Bond.h>


namespace desfact {
namespace descriptors {

// Connectivity & Topological Shape
class AtomCentralityVarianceDescriptor : public Descriptor {
public:
    AtomCentralityVarianceDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class PeripheryCoreRatioDescriptor : public FractionalDescriptor {
public:
    PeripheryCoreRatioDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class MinimumSpanningTreeDepthDescriptor : public Descriptor {
public:
    MinimumSpanningTreeDepthDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class TerminalBranchRatioDescriptor : public FractionalDescriptor {
public:
    TerminalBranchRatioDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class AverageAtomPathRedundancyDescriptor : public SumDescriptor {
public:
    AverageAtomPathRedundancyDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

// Atom-type Distribution
class CarbonClusterDensityDescriptor : public FractionalDescriptor {
public:
    CarbonClusterDensityDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class HeteroatomClusterDensityDescriptor : public FractionalDescriptor {
public:
    HeteroatomClusterDensityDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class TerminalHeavyAtomDiversityDescriptor : public Descriptor {
public:
    TerminalHeavyAtomDiversityDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class CentralHeteroatomRatioDescriptor : public FractionalDescriptor {
public:
    CentralHeteroatomRatioDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class ChainEndElementDiversityDescriptor : public Descriptor {
public:
    ChainEndElementDiversityDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

// Electronic and Polar Characteristics
class HighENNeighborDensityDescriptor : public FractionalDescriptor {
public:
    HighENNeighborDensityDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class ElectronDeficientCarbonRatioDescriptor : public FractionalDescriptor {
public:
    ElectronDeficientCarbonRatioDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class PolarizableAtomDispersionDescriptor : public Descriptor {
public:
    PolarizableAtomDispersionDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class DipoleMomentProxyDescriptor : public FractionalDescriptor {
public:
    DipoleMomentProxyDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

// Structural Patterns
class LongChainFragmentDensityDescriptor : public FractionalDescriptor {
public:
    LongChainFragmentDensityDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class ShortChainDensityDescriptor : public FractionalDescriptor {
public:
    ShortChainDensityDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class SubstitutionDensityPerRingDescriptor : public Descriptor {
public:
    SubstitutionDensityPerRingDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class LinearVsBranchedCarbonRatioDescriptor : public Descriptor {
public:
    LinearVsBranchedCarbonRatioDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class RingToBranchConnectivityRatioDescriptor : public FractionalDescriptor {
public:
    RingToBranchConnectivityRatioDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

// Hydrogen Bond Networks
class IsolatedHBondSiteDensityDescriptor : public FractionalDescriptor {
public:
    IsolatedHBondSiteDensityDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class FunctionalGroupSeparationIndexDescriptor : public Descriptor {
public:
    FunctionalGroupSeparationIndexDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class PeripheralHBondDensityDescriptor : public FractionalDescriptor {
public:
    PeripheralHBondDensityDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

// Aromatic and Conjugation Properties
class AromaticCoreRatioDescriptor : public FractionalDescriptor {
public:
    AromaticCoreRatioDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class ConjugationGapRatioDescriptor : public FractionalDescriptor {
public:
    ConjugationGapRatioDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class NonAromaticRingSubstitutionDescriptor : public FractionalDescriptor {
public:
    NonAromaticRingSubstitutionDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class AromaticToNonAromaticRingRatioDescriptor : public Descriptor {
public:
    AromaticToNonAromaticRingRatioDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class ConjugationTerminalDensityDescriptor : public FractionalDescriptor {
public:
    ConjugationTerminalDensityDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

// Charge Distribution Characteristics
class ChargePairDensityDescriptor : public Descriptor {
public:
    ChargePairDensityDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class FormalChargeAccessibilityDescriptor : public FractionalDescriptor {
public:
    FormalChargeAccessibilityDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class ChargeSeparationIndexDescriptor : public Descriptor {
public:
    ChargeSeparationIndexDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class NeutralAtomRatioDescriptor : public FractionalDescriptor {
public:
    NeutralAtomRatioDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class FunctionalGroupChargeBalanceDescriptor : public FractionalDescriptor {
public:
    FunctionalGroupChargeBalanceDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

// Spatial Configuration Proxies
class StereocenterDensityDescriptor : public FractionalDescriptor {
public:
    StereocenterDensityDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class DoubleBondConfigurationDensityDescriptor : public FractionalDescriptor {
public:
    DoubleBondConfigurationDensityDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class RingStereogenicAtomDensityDescriptor : public FractionalDescriptor {
public:
    RingStereogenicAtomDensityDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class TerminalStereocenterRatioDescriptor : public FractionalDescriptor {
public:
    TerminalStereocenterRatioDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class AdjacentStereocenterDensityDescriptor : public FractionalDescriptor {
public:
    AdjacentStereocenterDensityDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

// Elemental Composition Diversity
class GroupPeriodicDiversityDescriptor : public Descriptor {
public:
    GroupPeriodicDiversityDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class PeriodDiversityIndexDescriptor : public Descriptor {
public:
    PeriodDiversityIndexDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class AtomicMassVarianceDescriptor : public Descriptor {
public:
    AtomicMassVarianceDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class RareElementDensityDescriptor : public FractionalDescriptor {
public:
    RareElementDensityDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class AlkaliAlkalineEarthRatioDescriptor : public FractionalDescriptor {
public:
    AlkaliAlkalineEarthRatioDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

// Bonding Environment
class UnsaturationClustersDescriptor : public FractionalDescriptor {
public:
    UnsaturationClustersDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class RingSpanningBondDensityDescriptor : public FractionalDescriptor {
public:
    RingSpanningBondDensityDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class TerminalUnsaturatedBondRatioDescriptor : public FractionalDescriptor {
public:
    TerminalUnsaturatedBondRatioDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class HeavyAtomBondOrderVarianceDescriptor : public Descriptor {
public:
    HeavyAtomBondOrderVarianceDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class HeteroatomBondingDiversityDescriptor : public Descriptor {
public:
    HeteroatomBondingDiversityDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

// Charge Distribution Characteristics - add with other similar descriptors
class LocalizedChargeClustersDescriptor : public Descriptor {
public:
    LocalizedChargeClustersDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

} // namespace descriptors
} // namespace desfact
