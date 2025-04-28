#pragma once

#include "descriptors.hpp"
#include <GraphMol/GraphMol.h>
#include <GraphMol/Atom.h>
#include <GraphMol/Bond.h>
#include <GraphMol/MolOps.h>
#include <functional>

namespace desfact {
namespace descriptors {

// Base class for sum descriptors
class SumDescriptor : public Descriptor {
public:
    SumDescriptor(const std::string& name, const std::string& description)
        : Descriptor(name, description) {}
    
    // Helper methods for sum calculations
    static double calcAtomicSum(const Molecule& mol, 
                               const std::function<double(const RDKit::Atom*)>& valueFunction,
                               bool normalize = false);
    
    static double calcBondSum(const Molecule& mol,
                             const std::function<double(const RDKit::Bond*)>& valueFunction,
                             bool normalize = false);
                            
    // Property getters
    static double getAtomElectronegativity(const RDKit::Atom* atom);
    static double getAtomCovalentRadius(const RDKit::Atom* atom);
    static double getAtomVdWVolume(const RDKit::Atom* atom);
    static double getAtomPolarizability(const RDKit::Atom* atom);
    static double getAtomIonizationEnergy(const RDKit::Atom* atom);
    static double getAtomElectronAffinity(const RDKit::Atom* atom);
    
    // Helper predicates
    static bool isHeteroatom(const RDKit::Atom* atom);
    static bool isHalogen(const RDKit::Atom* atom);
    static bool isMetalAtom(const RDKit::Atom* atom);
    static bool isPolarBond(const RDKit::Bond* bond);
};

// Atomic property sums
class SumENDescriptor : public SumDescriptor {
public:
    SumENDescriptor() : SumDescriptor("SumEN", "Sum of atom electronegativities") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class SumCovalentRadiiDescriptor : public SumDescriptor {
public:
    SumCovalentRadiiDescriptor() : SumDescriptor("SumCovalentRadii", "Sum of atomic covalent radii") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class SumVdWVolumeDescriptor : public SumDescriptor {
public:
    SumVdWVolumeDescriptor() : SumDescriptor("SumVdWVolume", "Sum of Van der Waals atomic volumes") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class SumPolarizabilityDescriptor : public SumDescriptor {
public:
    SumPolarizabilityDescriptor() : SumDescriptor("SumPolarizability", "Sum of atomic polarizabilities") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class SumIonizationEnergyDescriptor : public SumDescriptor {
public:
    SumIonizationEnergyDescriptor() : SumDescriptor("SumIonizationEnergy", "Sum of first ionization energies") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class SumElectronAffinityDescriptor : public SumDescriptor {
public:
    SumElectronAffinityDescriptor() : SumDescriptor("SumElectronAffinity", "Sum of electron affinities") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

// Atom Class Sums
class SumAtomsDescriptor : public SumDescriptor {
public:
    SumAtomsDescriptor() : SumDescriptor("SumAtoms", "Total number of atoms") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class SumHeavyAtomsDescriptor : public SumDescriptor {
public:
    SumHeavyAtomsDescriptor() : SumDescriptor("SumHeavyAtoms", "Total number of atoms with Z > 1") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class SumHeteroatomsDescriptor : public SumDescriptor {
public:
    SumHeteroatomsDescriptor() : SumDescriptor("SumHeteroatoms", "Total number of non-C, non-H atoms") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class SumHalogensDescriptor : public SumDescriptor {
public:
    SumHalogensDescriptor() : SumDescriptor("SumHalogens", "Count of F, Cl, Br, I atoms") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class SumChargedAtomsDescriptor : public SumDescriptor {
public:
    SumChargedAtomsDescriptor() : SumDescriptor("SumChargedAtoms", "Count of atoms with formal charge ≠ 0") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

// Bond and Hybridization Sums
class SumBondsDescriptor : public SumDescriptor {
public:
    SumBondsDescriptor() : SumDescriptor("SumBonds", "Total number of bonds") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class SumDoubleBondsDescriptor : public SumDescriptor {
public:
    SumDoubleBondsDescriptor() : SumDescriptor("SumDoubleBonds", "Count of double bonds") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class SumTripleBondsDescriptor : public SumDescriptor {
public:
    SumTripleBondsDescriptor() : SumDescriptor("SumTripleBonds", "Count of triple bonds") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class SumAromaticBondsDescriptor : public SumDescriptor {
public:
    SumAromaticBondsDescriptor() : SumDescriptor("SumAromaticBonds", "Count of aromatic bonds") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class SumPolarBondsDescriptor : public SumDescriptor {
public:
    SumPolarBondsDescriptor() : SumDescriptor("SumPolarBonds", "Count of bonds with significant EN difference") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class SumUnpolarBondsDescriptor : public SumDescriptor {
public:
    SumUnpolarBondsDescriptor() : SumDescriptor("SumUnpolarBonds", "Count of bonds with small or zero EN difference") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class SumSp3BondsDescriptor : public SumDescriptor {
public:
    SumSp3BondsDescriptor() : SumDescriptor("SumSp3Bonds", "Count of sp3-hybridized central atom bonds") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class SumSp2BondsDescriptor : public SumDescriptor {
public:
    SumSp2BondsDescriptor() : SumDescriptor("SumSp2Bonds", "Count of sp2-hybridized central atom bonds") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

// Structural Sums
class SumRingsDescriptor : public SumDescriptor {
public:
    SumRingsDescriptor() : SumDescriptor("SumRings", "Total number of rings") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class SumAromaticRingsDescriptor : public SumDescriptor {
public:
    SumAromaticRingsDescriptor() : SumDescriptor("SumAromaticRings", "Total number of aromatic rings") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class SumRotatableBondsDescriptor : public SumDescriptor {
public:
    SumRotatableBondsDescriptor() : SumDescriptor("SumRotatableBonds", "Total number of rotatable bonds") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class SumBridgeAtomsDescriptor : public SumDescriptor {
public:
    SumBridgeAtomsDescriptor() : SumDescriptor("SumBridgeAtoms", "Number of atoms connecting multiple rings") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class SumRingAtomsDescriptor : public SumDescriptor {
public:
    SumRingAtomsDescriptor() : SumDescriptor("SumRingAtoms", "Number of atoms part of any ring") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

// Physicochemical Contribution Sums
class SumENMWDescriptor : public SumDescriptor {
public:
    SumENMWDescriptor() : SumDescriptor("SumENMW", "Sum of (electronegativity × atomic weight)") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class SumPolMWDescriptor : public SumDescriptor {
public:
    SumPolMWDescriptor() : SumDescriptor("SumPolMW", "Sum of (polarizability × atomic weight)") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class SumENRcovDescriptor : public SumDescriptor {
public:
    SumENRcovDescriptor() : SumDescriptor("SumENRcov", "Sum of (electronegativity × covalent radius)") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class SumENPolDescriptor : public SumDescriptor {
public:
    SumENPolDescriptor() : SumDescriptor("SumENPol", "Sum of (electronegativity × polarizability)") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

// Additional Hybridization Sums
class SumSp3AtomsDescriptor : public SumDescriptor {
public:
    SumSp3AtomsDescriptor() : SumDescriptor("SumSp3", "Total count of sp3-hybridized atoms") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class SumSp2AtomsDescriptor : public SumDescriptor {
public:
    SumSp2AtomsDescriptor() : SumDescriptor("SumSp2", "Total count of sp2-hybridized atoms") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

// Additional Pharmacophore-related Sums
class SumHBDonorsDescriptor : public SumDescriptor {
public:
    SumHBDonorsDescriptor() : SumDescriptor("SumHBDonors", "Total number of hydrogen bond donors") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class SumHBAcceptorsDescriptor : public SumDescriptor {
public:
    SumHBAcceptorsDescriptor() : SumDescriptor("SumHBAcceptors", "Total number of hydrogen bond acceptors") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

// Additional Functional Group Sums
class SumAmideGroupsDescriptor : public SumDescriptor {
public:
    SumAmideGroupsDescriptor() : SumDescriptor("SumAmideGroups", "Count of amide functional groups") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class SumCarboxylGroupsDescriptor : public SumDescriptor {
public:
    SumCarboxylGroupsDescriptor() : SumDescriptor("SumCarboxylGroups", "Count of carboxyl functional groups") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

// Additional Sum descriptors from the requested list
class SumElectronegativityRingDescriptor : public SumDescriptor {
public:
    SumElectronegativityRingDescriptor() : SumDescriptor("SumElectronegativityRing", "EN sum over atoms in rings") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class SumElectronegativityBondedDescriptor : public SumDescriptor {
public:
    SumElectronegativityBondedDescriptor() : SumDescriptor("SumElectronegativityBonded", "EN sum for bonded atoms only") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class SumAtomicRadiusHeavyDescriptor : public SumDescriptor {
public:
    SumAtomicRadiusHeavyDescriptor() : SumDescriptor("SumAtomicRadiusHeavy", "Covalent radius sum for heavy atoms") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class SumPolzENRatioDescriptor : public SumDescriptor {
public:
    SumPolzENRatioDescriptor() : SumDescriptor("SumPolzENRatio", "Sum of polarizability / EN") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class SumIEENRatioDescriptor : public SumDescriptor {
public:
    SumIEENRatioDescriptor() : SumDescriptor("SumIEENRatio", "Sum of IE / EN") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class SumEAMWENDescriptor : public SumDescriptor {
public:
    SumEAMWENDescriptor() : SumDescriptor("SumEAMWEN", "Sum of EA × MW × EN") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class SumSp2CAtomsDescriptor : public SumDescriptor {
public:
    SumSp2CAtomsDescriptor() : SumDescriptor("SumSp2CAtoms", "Count of sp2-hybridized carbon atoms") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class SumSpAtomsENDescriptor : public SumDescriptor {
public:
    SumSpAtomsENDescriptor() : SumDescriptor("SumSpAtomsEN", "Sum of sp-hybridized atoms × EN") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class SumFormalChargeAbsDescriptor : public SumDescriptor {
public:
    SumFormalChargeAbsDescriptor() : SumDescriptor("SumFormalChargeAbs", "Sum of absolute formal charges") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class SumFormalChargeHeavyDescriptor : public SumDescriptor {
public:
    SumFormalChargeHeavyDescriptor() : SumDescriptor("SumFormalChargeHeavy", "Formal charge sum over heavy atoms") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class SumFormalChargePositiveDescriptor : public SumDescriptor {
public:
    SumFormalChargePositiveDescriptor() : SumDescriptor("SumFormalChargePositive", "Sum of positive formal charges") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class SumFormalChargeNegativeDescriptor : public SumDescriptor {
public:
    SumFormalChargeNegativeDescriptor() : SumDescriptor("SumFormalChargeNegative", "Sum of negative formal charges") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class SumENMWRatioDescriptor : public SumDescriptor {
public:
    SumENMWRatioDescriptor() : SumDescriptor("SumENMWRatio", "Sum of EN / MW") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class SumPolzRingDescriptor : public SumDescriptor {
public:
    SumPolzRingDescriptor() : SumDescriptor("SumPolzRing", "Polarizability sum for atoms in ring systems") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class SumIEBondedDescriptor : public SumDescriptor {
public:
    SumIEBondedDescriptor() : SumDescriptor("SumIEBonded", "IE sum for atoms with ≥1 bond") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class SumEAHeavyBondedDescriptor : public SumDescriptor {
public:
    SumEAHeavyBondedDescriptor() : SumDescriptor("SumEAHeavyBonded", "EA sum over heavy atoms with ≥1 bond") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class SumSp3ENWeightedDescriptor : public SumDescriptor {
public:
    SumSp3ENWeightedDescriptor() : SumDescriptor("SumSp3ENWeighted", "Sum of EN × sp3 atoms") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class SumHeteroMWENDescriptor : public SumDescriptor {
public:
    SumHeteroMWENDescriptor() : SumDescriptor("SumHeteroMWEN", "Sum of MW × EN for heteroatoms") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class SumENSp2HeavyDescriptor : public SumDescriptor {
public:
    SumENSp2HeavyDescriptor() : SumDescriptor("SumENSp2Heavy", "EN of heavy sp2 atoms") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class SumENRadicalAtomsDescriptor : public SumDescriptor {
public:
    SumENRadicalAtomsDescriptor() : SumDescriptor("SumENRadicalAtoms", "EN of radical atoms") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class SumOxHeteroMWDescriptor : public SumDescriptor {
public:
    SumOxHeteroMWDescriptor() : SumDescriptor("SumOxHeteroMW", "Sum of ox state × MW for heteroatoms") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class SumENRingHeavyDescriptor : public SumDescriptor {
public:
    SumENRingHeavyDescriptor() : SumDescriptor("SumENRingHeavy", "EN of heavy atoms in rings") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class SumPolENBondedDescriptor : public SumDescriptor {
public:
    SumPolENBondedDescriptor() : SumDescriptor("SumPolENBonded", "Sum of polarizability × EN for bonded atoms") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

} // namespace descriptors
} // namespace desfact 