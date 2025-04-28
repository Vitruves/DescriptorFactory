#pragma once

#include "descriptors.hpp" // Base class
#include <string>
#include <variant>

// Forward declare Molecule in the global ::desfact namespace, not locally
namespace desfact { class Molecule; }

namespace desfact {
namespace descriptors {

// --- Novel 2-D Descriptors (Group 6) ---

// 1. DegreeLogWeightedSum - Σ deg(i) · ln [deg(i)+1] over all heavy atoms
class DegreeLogWeightedSum : public Descriptor {
public:
    DegreeLogWeightedSum();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

// 2. BondOrderSkewness - Fisher skewness of the bond-order distribution (single = 1, double = 2, …)
class BondOrderSkewness : public Descriptor {
public:
    BondOrderSkewness();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

// 3. PathLengthGini - Gini coefficient of all heavy-atom shortest-path lengths
class PathLengthGini : public Descriptor {
public:
    PathLengthGini();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

// 4. LongestHomoelementPath - Length of the longest path whose atoms share an identical element symbol
class LongestHomoelementPath : public Descriptor {
public:
    LongestHomoelementPath();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

// 5. PeripheralAtomTypeDiversity - Shannon entropy of element types among terminal heavy atoms (degree = 1)
class PeripheralAtomTypeDiversity : public Descriptor {
public:
    PeripheralAtomTypeDiversity();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

// 6. HeteroRingEdgeRatio - Heteroatoms directly attached to rings but not ring members ÷ total heteroatoms
class HeteroRingEdgeRatio : public Descriptor {
public:
    HeteroRingEdgeRatio();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

// 7. BranchPointDensity - Branching atoms (degree > 2) per heavy atom
class BranchPointDensity : public Descriptor {
public:
    BranchPointDensity();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

// 8. MeanBranchSeparation - Mean shortest-path distance between all pairs of branch points
class MeanBranchSeparation : public Descriptor {
public:
    MeanBranchSeparation();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

// 9. CarbonChainBranchingIndex - Branch points located on the longest carbon chain ÷ chain length
class CarbonChainBranchingIndex : public Descriptor {
public:
    CarbonChainBranchingIndex();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

// 10. RingSystemCount - Number of distinct fused-ring systems (SSSR clusters)
class RingSystemCount : public Descriptor {
public:
    RingSystemCount();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

// 11. MeanRingPerimeterDegree - Average atomic degree among all atoms that belong to rings
class MeanRingPerimeterDegree : public Descriptor {
public:
    MeanRingPerimeterDegree();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

// 12. AtomTypeNeighborDiversity - For every element type, count distinct neighbor element types; mean over elements
class AtomTypeNeighborDiversity : public Descriptor {
public:
    AtomTypeNeighborDiversity();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

// 13. BondTypeAlternationRatio - Bonds that are part of an alternating single/double pattern of length ≥ 3 ÷ total bonds
class BondTypeAlternationRatio : public Descriptor {
public:
    BondTypeAlternationRatio();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

// 14. AtomValenceSpread - Max formal valence − Min formal valence among heavy atoms
class AtomValenceSpread : public Descriptor {
public:
    AtomValenceSpread();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

// 15. RingHydrogenRatio - Implicit H on ring atoms ÷ total implicit H in the molecule
class RingHydrogenRatio : public Descriptor {
public:
    RingHydrogenRatio();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

// 16. MaxConsecutiveSingleBonds - Length of the longest path consisting solely of single bonds
class MaxConsecutiveSingleBonds : public Descriptor {
public:
    MaxConsecutiveSingleBonds();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

// 17. HeteroatomPathFraction - Fraction of shortest paths (length ≤ 4) whose endpoints are heteroatoms
class HeteroatomPathFraction : public Descriptor {
public:
    HeteroatomPathFraction();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

// 18. MeanResonanceBondOrder - Mean bond order across conjugated (resonance-capable) bonds
class MeanResonanceBondOrder : public Descriptor {
public:
    MeanResonanceBondOrder();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

// 19. CarbonIsotopeCount - Count of carbon atoms carrying an explicit isotope label in SMILES
class CarbonIsotopeCount : public Descriptor {
public:
    CarbonIsotopeCount();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

// 20. ChargeBalanceSkew - Skewness of formal charge distribution
class ChargeBalanceSkew : public Descriptor {
public:
    ChargeBalanceSkew();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

// 21. ElementRunLengthVariance - Variance of consecutive same-element run lengths in canonical SMILES
class ElementRunLengthVariance : public Descriptor {
public:
    ElementRunLengthVariance();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

// 22. SSSRtoRingAtomRatio - Number of SSSR rings ÷ total number of ring atoms
class SSSRtoRingAtomRatio : public Descriptor {
public:
    SSSRtoRingAtomRatio();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

// 23. MaxRingDistance - Largest graph distance between any two atoms within the same ring
class MaxRingDistance : public Descriptor {
public:
    MaxRingDistance();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

// 24. MeanBondStereoStates - Average stereo code per bond (none = 0, / or \ = 1, E/Z = 2)
class MeanBondStereoStates : public Descriptor {
public:
    MeanBondStereoStates();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

// 25. TerminalDoubleBondCount - Count of atoms in double bonds that have degree = 1
class TerminalDoubleBondCount : public Descriptor {
public:
    TerminalDoubleBondCount();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

// 26. AromaticSubstituentDiversity - Shannon entropy of element types directly attached to aromatic atoms but outside rings
class AromaticSubstituentDiversity : public Descriptor {
public:
    AromaticSubstituentDiversity();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

// 27. HalogenNeighborHybridRatio - Halogen atoms bound to sp³ carbon ÷ total halogen atoms
class HalogenNeighborHybridRatio : public Descriptor {
public:
    HalogenNeighborHybridRatio();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

// 28. MeanAtomBetweennessCentrality - Average betweenness centrality of heavy atoms (un-weighted graph)
class MeanAtomBetweennessCentrality : public Descriptor {
public:
    MeanAtomBetweennessCentrality();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

// 29. BondOrderEntropy - Shannon entropy of bond-order frequencies
class BondOrderEntropy : public Descriptor {
public:
    BondOrderEntropy();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

// 30. HybridizationSymmetryIndex - Symmetry measure of hybridization distribution
class HybridizationSymmetryIndex : public Descriptor {
public:
    HybridizationSymmetryIndex();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

// 31. RingSizeMedian - Median ring size among all SSSR rings (0 if no rings)
class RingSizeMedian : public Descriptor {
public:
    RingSizeMedian();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

// 32. LongestPathBondOrderProduct - Product of bond orders along the molecule's longest simple path
class LongestPathBondOrderProduct : public Descriptor {
public:
    LongestPathBondOrderProduct();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

// 33. HeteroatomDegreeVariance - Variance of degrees among heteroatoms
class HeteroatomDegreeVariance : public Descriptor {
public:
    HeteroatomDegreeVariance();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

// 34. NonCarbonAtomAdjacencyCount - Number of bonds in which neither atom is carbon
class NonCarbonAtomAdjacencyCount : public Descriptor {
public:
    NonCarbonAtomAdjacencyCount();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

// 35. MaxBondENGap - Maximum absolute Pauling electronegativity difference across any single bond
class MaxBondENGap : public Descriptor {
public:
    MaxBondENGap();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

// 36. RingIndexSum - Σ (1 / ring size) over all SSSR rings
class RingIndexSum : public Descriptor {
public:
    RingIndexSum();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

// 37. BranchDepthAverage - Mean shortest-path distance from each branch point to the nearest terminal atom
class BranchDepthAverage : public Descriptor {
public:
    BranchDepthAverage();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

// 38. ConsecutiveHeteroBondFraction - Bonds connecting two heteroatoms ÷ total bonds
class ConsecutiveHeteroBondFraction : public Descriptor {
public:
    ConsecutiveHeteroBondFraction();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

// 39. HeavyAtomHydrogenRatio - (implicit H + explicit H) ÷ heavy atoms
class HeavyAtomHydrogenRatio : public Descriptor {
public:
    HeavyAtomHydrogenRatio();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

// 40. MeanAtomicMassPerDegree - Mean (atomic mass ÷ degree) over heavy atoms
class MeanAtomicMassPerDegree : public Descriptor {
public:
    MeanAtomicMassPerDegree();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

// 41. RareElementCount - Count of atoms with atomic number > 17
class RareElementCount : public Descriptor {
public:
    RareElementCount();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

// 42. TopoDistanceSkewness - Fisher skewness of all-pairs shortest-path length distribution
class TopoDistanceSkewness : public Descriptor {
public:
    TopoDistanceSkewness();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

// 43. PercentageIsolatedRingSystems - Isolated (non-fused) ring systems ÷ total ring systems (0 when no rings)
class PercentageIsolatedRingSystems : public Descriptor {
public:
    PercentageIsolatedRingSystems();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

// 44. ValenceMismatchBondCount - Bonds whose two atoms differ in allowed valence by ≥ 3
class ValenceMismatchBondCount : public Descriptor {
public:
    ValenceMismatchBondCount();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

// 45. HeteroAtomSequenceLengthMax - Longest run of consecutive heteroatom symbols in canonical SMILES
class HeteroAtomSequenceLengthMax : public Descriptor {
public:
    HeteroAtomSequenceLengthMax();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

// 46. ChiralAtomSkew - (R centers − S centers) ÷ total chiral centers (0 if none)
class ChiralAtomSkew : public Descriptor {
public:
    ChiralAtomSkew();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

// 47. RingAtomUnsaturationRatio - Unsaturated ring atoms (sp²/sp) ÷ ring atoms
class RingAtomUnsaturationRatio : public Descriptor {
public:
    RingAtomUnsaturationRatio();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

// 48. TerminalHeteroBondOrderSum - Sum of bond orders for bonds involving terminal heteroatoms
class TerminalHeteroBondOrderSum : public Descriptor {
public:
    TerminalHeteroBondOrderSum();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

// 49. AromaticNonAromaticBondCount - Count of bonds that connect an aromatic atom to a non-aromatic heavy atom
class AromaticNonAromaticBondCount : public Descriptor {
public:
    AromaticNonAromaticBondCount();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

// 50. RingAtomChargeVariance - Variance of formal charges restricted to ring atoms (0 if no rings or no charges)
class RingAtomChargeVariance : public Descriptor {
public:
    RingAtomChargeVariance();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

} // namespace descriptors
} // namespace desfact