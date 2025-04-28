#pragma once

#include "descriptors.hpp" // Base class
#include <string>
#include <variant>

// Forward declare Molecule in the global ::desfact namespace, not locally
namespace desfact { class Molecule; }

namespace desfact {
namespace descriptors {

// --- Novel 2-D Descriptors (Group 7) ---

// 1. HeteroatomFractionOfRings - Fraction of ring atoms that are heteroatoms
class HeteroatomFractionOfRings : public Descriptor {
public:
    HeteroatomFractionOfRings();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

// 2. AromaticRingRatio - Ratio of aromatic rings to total rings
class AromaticRingRatio : public Descriptor {
public:
    AromaticRingRatio();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

// 3. RingComplexityIndex - Measure of ring complexity based on size and connectivity
class RingComplexityIndex : public Descriptor {
public:
    RingComplexityIndex();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

// 4. ElectronegativeAtomDensity - Density of electronegative atoms (N, O, F, Cl, Br, I)
class ElectronegativeAtomDensity : public Descriptor {
public:
    ElectronegativeAtomDensity();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

// 5. ChargeAsymmetryIndex - Measure of charge distribution asymmetry
class ChargeAsymmetryIndex : public Descriptor {
public:
    ChargeAsymmetryIndex();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

// 6. TopologicalPolarSurfaceEfficiency - TPSA normalized by molecular weight
class TopologicalPolarSurfaceEfficiency : public Descriptor {
public:
    TopologicalPolarSurfaceEfficiency();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

// 7. HeterocyclicRingCount - Count of rings containing at least one heteroatom
class HeterocyclicRingCount : public Descriptor {
public:
    HeterocyclicRingCount();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

// 8. CarbocyclicRingCount - Count of rings containing only carbon atoms
class CarbocyclicRingCount : public Descriptor {
public:
    CarbocyclicRingCount();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

// 9. FusedRingSystems - Count of fused ring systems
class FusedRingSystems : public Descriptor {
public:
    FusedRingSystems();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

// 10. AverageFusedRingSize - Average size of fused ring systems
class AverageFusedRingSize : public Descriptor {
public:
    AverageFusedRingSize();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

// 11. HeterocyclicFusedRingRatio - Ratio of heterocyclic rings in fused systems
class HeterocyclicFusedRingRatio : public Descriptor {
public:
    HeterocyclicFusedRingRatio();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

// 12. BridgedRingCount - Count of bridged rings
class BridgedRingCount : public Descriptor {
public:
    BridgedRingCount();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

// 13. SpiroRingCount - Count of spiro rings
class SpiroRingCount : public Descriptor {
public:
    SpiroRingCount();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

// 14. MacrocyclicRingCount - Count of rings with size >= 12
class MacrocyclicRingCount : public Descriptor {
public:
    MacrocyclicRingCount();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

// 15. HydrogenBondAcceptorDensity - Density of H-bond acceptors
class HydrogenBondAcceptorDensity : public Descriptor {
public:
    HydrogenBondAcceptorDensity();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

// 16. HydrogenBondDonorDensity - Density of H-bond donors
class HydrogenBondDonorDensity : public Descriptor {
public:
    HydrogenBondDonorDensity();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

// 17. RotatableBondDensity - Density of rotatable bonds
class RotatableBondDensity : public Descriptor {
public:
    RotatableBondDensity();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

// 18. FunctionalGroupDiversity - Shannon entropy of functional group types
class FunctionalGroupDiversity : public Descriptor {
public:
    FunctionalGroupDiversity();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

// 19. AtomTypeEntropy - Shannon entropy of atom types
class AtomTypeEntropy : public Descriptor {
public:
    AtomTypeEntropy();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

// 20. BondTypeEntropy - Shannon entropy of bond types
class BondTypeEntropy : public Descriptor {
public:
    BondTypeEntropy();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

// 21. RingTypeEntropy - Shannon entropy of ring types
class RingTypeEntropy : public Descriptor {
public:
    RingTypeEntropy();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

// 22. ElectronegativityVariance - Variance of atom electronegativity values
class ElectronegativityVariance : public Descriptor {
public:
    ElectronegativityVariance();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

} // namespace descriptors
} // namespace desfact