#pragma once

#include "descriptors.hpp"
#include <GraphMol/Atom.h>
#include <GraphMol/Bond.h>
#include <GraphMol/ROMol.h>
#include <functional>
#include <string>
#include <cmath>

namespace desfact {
namespace descriptors {

// Base class for all fractional descriptors
class FractionalDescriptor : public Descriptor {
public:
    FractionalDescriptor(const std::string& name, const std::string& description)
        : Descriptor(name, description) {}
    
    // Helper methods for fractional descriptor calculations
    static double calcAtomicFraction(const Molecule& mol, const std::function<bool(const RDKit::Atom*)>& predicate);
    static double calcBondFraction(const Molecule& mol, const std::function<bool(const RDKit::Bond*)>& predicate);
    
    // Utility for element-based checks
    static bool isElement(const RDKit::Atom* atom, int atomicNum);
    
    // Common atom property checks
    static bool isHalogen(const RDKit::Atom* atom);
    static bool isHeteroatom(const RDKit::Atom* atom);
    static bool isPolarAtom(const RDKit::Atom* atom);
    static bool isMetalAtom(const RDKit::Atom* atom);
    static bool isMetalloidAtom(const RDKit::Atom* atom);
    
    // Electronegativity utilities
    static double getAtomElectronegativity(const RDKit::Atom* atom);
    static double getMoleculeAverageElectronegativity(const RDKit::ROMol* mol);
    
    // Add this line:
    static double getAtomPolarizability(const RDKit::Atom* atom);
};

// Atomic Fraction Descriptors (elemental contribution to MW)
class FcCDescriptor : public FractionalDescriptor {
public:
    FcCDescriptor() : FractionalDescriptor("FcC", "Fractional contribution of Carbon to molecular weight") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class FcFDescriptor : public FractionalDescriptor {
public:
    FcFDescriptor() : FractionalDescriptor("FcF", "Fluorine contribution to molecular weight") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class FcODescriptor : public FractionalDescriptor {
public:
    FcODescriptor() : FractionalDescriptor("FcO", "Oxygen contribution to molecular weight") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class FcClDescriptor : public FractionalDescriptor {
public:
    FcClDescriptor() : FractionalDescriptor("FcCl", "Chlorine contribution to molecular weight") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class FcBrDescriptor : public FractionalDescriptor {
public:
    FcBrDescriptor() : FractionalDescriptor("FcBr", "Bromine contribution to molecular weight") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class FcIDescriptor : public FractionalDescriptor {
public:
    FcIDescriptor() : FractionalDescriptor("FcI", "Iodine contribution to molecular weight") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class FcSDescriptor : public FractionalDescriptor {
public:
    FcSDescriptor() : FractionalDescriptor("FcS", "Sulfur contribution to molecular weight") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class FcNDescriptor : public FractionalDescriptor {
public:
    FcNDescriptor() : FractionalDescriptor("FcN", "Nitrogen contribution to molecular weight") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class FcHDescriptor : public FractionalDescriptor {
public:
    FcHDescriptor() : FractionalDescriptor("FcH", "Hydrogen contribution to molecular weight") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class FcHaloDescriptor : public FractionalDescriptor {
public:
    FcHaloDescriptor() : FractionalDescriptor("FcHalo", "Combined halogen contribution: F, Cl, Br, I") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class FcHeteroDescriptor : public FractionalDescriptor {
public:
    FcHeteroDescriptor() : FractionalDescriptor("FcHetero", "Combined heteroatom contribution: all non-C, non-H atoms") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class FcPolarDescriptor : public FractionalDescriptor {
public:
    FcPolarDescriptor() : FractionalDescriptor("FcPolar", "Fraction of atoms with electronegativity > C") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class FcApolarDescriptor : public FractionalDescriptor {
public:
    FcApolarDescriptor() : FractionalDescriptor("FcApolar", "Fraction of atoms with electronegativity ≤ C") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

// Hybridization-Based Bond Fractions
class FcCSp3Descriptor : public FractionalDescriptor {
public:
    FcCSp3Descriptor() : FractionalDescriptor("FcCSp3", "Proportion of carbon-carbon sp3 bonds") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class FcCSp2Descriptor : public FractionalDescriptor {
public:
    FcCSp2Descriptor() : FractionalDescriptor("FcCSp2", "Proportion of carbon-carbon sp2 bonds") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

// Bond Polarity Descriptors
class FcUnpolDescriptor : public FractionalDescriptor {
public:
    FcUnpolDescriptor() : FractionalDescriptor("FcUnpol", "Proportion of bonds between atoms of equal electronegativity") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class FcPolDescriptor : public FractionalDescriptor {
public:
    FcPolDescriptor() : FractionalDescriptor("FcPol", "Proportion of polar bonds (atoms with different electronegativity)") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

// Electronegativity-Based Descriptors
class FcSumPolAtDescriptor : public FractionalDescriptor {
public:
    FcSumPolAtDescriptor() : FractionalDescriptor("FcSumPolAt", "Average electronegativity per atom") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class FcSumPolMWDescriptor : public FractionalDescriptor {
public:
    FcSumPolMWDescriptor() : FractionalDescriptor("FcSumPolMW", "Average electronegativity normalized by molecular weight") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

// Bond Valence Descriptors
class FcBondAtDescriptor : public FractionalDescriptor {
public:
    FcBondAtDescriptor() : FractionalDescriptor("FcBondAt", "Total bond valence (sum of bond orders) divided by atom count") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class FcBondNDescriptor : public FractionalDescriptor {
public:
    FcBondNDescriptor() : FractionalDescriptor("FcBondN", "Fraction of non-single bonds on nitrogen atoms") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class FcBondODescriptor : public FractionalDescriptor {
public:
    FcBondODescriptor() : FractionalDescriptor("FcBondO", "Fraction of non-single bonds on oxygen atoms") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class FcBondCDescriptor : public FractionalDescriptor {
public:
    FcBondCDescriptor() : FractionalDescriptor("FcBondC", "Fraction of non-single bonds on carbon atoms") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class FcBondSDescriptor : public FractionalDescriptor {
public:
    FcBondSDescriptor() : FractionalDescriptor("FcBondS", "Fraction of non-single bonds on sulfur atoms") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class FcBondPDescriptor : public FractionalDescriptor {
public:
    FcBondPDescriptor() : FractionalDescriptor("FcBondP", "Fraction of non-single bonds on phosphorus atoms") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

// Additional normalized descriptors
class FcHBDonorsDescriptor : public FractionalDescriptor {
public:
    FcHBDonorsDescriptor() : FractionalDescriptor("FcHBDonors", "Fraction of atoms that are H-bond donors") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class FcHBAcceptorsDescriptor : public FractionalDescriptor {
public:
    FcHBAcceptorsDescriptor() : FractionalDescriptor("FcHBAcceptors", "Fraction of atoms that are H-bond acceptors") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class FcAromaticAtomsDescriptor : public FractionalDescriptor {
public:
    FcAromaticAtomsDescriptor() : FractionalDescriptor("FcAromaticAtoms", "Fraction of atoms part of aromatic systems") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class FcRingAtomsDescriptor : public FractionalDescriptor {
public:
    FcRingAtomsDescriptor() : FractionalDescriptor("FcRingAtoms", "Fraction of atoms included in rings") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class FcBridgeAtomsDescriptor : public FractionalDescriptor {
public:
    FcBridgeAtomsDescriptor() : FractionalDescriptor("FcBridgeAtoms", "Fraction of atoms acting as ring junctions or bridges") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class FcChargedAtomsDescriptor : public FractionalDescriptor {
public:
    FcChargedAtomsDescriptor() : FractionalDescriptor("FcChargedAtoms", "Fraction of atoms carrying formal charge") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class FcHeavyAtomsDescriptor : public FractionalDescriptor {
public:
    FcHeavyAtomsDescriptor() : FractionalDescriptor("FcHeavyAtoms", "Fraction of atoms with atomic number > 1") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class FcMetalsDescriptor : public FractionalDescriptor {
public:
    FcMetalsDescriptor() : FractionalDescriptor("FcMetals", "Fraction of metal atoms") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

// Electronegativity-Based Fractions
class FcENAboveAvgDescriptor : public FractionalDescriptor {
public:
    FcENAboveAvgDescriptor() : FractionalDescriptor("FcENAboveAvg", "Fraction of atoms with electronegativity above average") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class FcENBelowAvgDescriptor : public FractionalDescriptor {
public:
    FcENBelowAvgDescriptor() : FractionalDescriptor("FcENBelowAvg", "Fraction of atoms with electronegativity below average") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class FcENHighDescriptor : public FractionalDescriptor {
public:
    FcENHighDescriptor() : FractionalDescriptor("FcENHigh", "Fraction of atoms with electronegativity > 3.5") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class FcENLowDescriptor : public FractionalDescriptor {
public:
    FcENLowDescriptor() : FractionalDescriptor("FcENLow", "Fraction of atoms with electronegativity < 2.0") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

// Atomic Radius
class FcSmallRDescriptor : public FractionalDescriptor {
public:
    FcSmallRDescriptor() : FractionalDescriptor("FcSmallR", "Fraction of atoms with covalent radius < 0.77 Å") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class FcLargeRDescriptor : public FractionalDescriptor {
public:
    FcLargeRDescriptor() : FractionalDescriptor("FcLargeR", "Fraction of atoms with covalent radius > 1.4 Å") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

// Polarizability
class FcLowPolzDescriptor : public FractionalDescriptor {
public:
    FcLowPolzDescriptor() : FractionalDescriptor("FcLowPolz", "Fraction of atoms with atomic polarizability < 5 Å³") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class FcHighPolzDescriptor : public FractionalDescriptor {
public:
    FcHighPolzDescriptor() : FractionalDescriptor("FcHighPolz", "Fraction of atoms with polarizability > 15 Å³") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

// Ionization Energy
// class FcLowIEDescriptor : public FractionalDescriptor {
// public:
//     FcLowIEDescriptor() : FractionalDescriptor("FcLowIE", "Fraction of atoms with ionization energy < 7 eV") {}
//     std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
// };

// Electron Affinity
class FcHighEADescriptor : public FractionalDescriptor {
public:
    FcHighEADescriptor() : FractionalDescriptor("FcHighEA", "Fraction of atoms with electron affinity > 2 eV") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class FcLowEADescriptor : public FractionalDescriptor {
public:
    FcLowEADescriptor() : FractionalDescriptor("FcLowEA", "Fraction with EA < 0.5 eV") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

// Van der Waals Volume
class FcSmallVdWDescriptor : public FractionalDescriptor {
public:
    FcSmallVdWDescriptor() : FractionalDescriptor("FcSmallVdW", "Fraction of atoms with VdW volume < 15 Å³") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class FcLargeVdWDescriptor : public FractionalDescriptor {
public:
    FcLargeVdWDescriptor() : FractionalDescriptor("FcLargeVdW", "Fraction of atoms with VdW volume > 25 Å³") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

// Electronegativity Contribution
// class FcENMWDescriptor : public FractionalDescriptor {
// public:
//     FcENMWDescriptor() : FractionalDescriptor("FcENMW", "Sum(atom electronegativity × atomic MW) / total MW") {}
//     std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
// };

class FcENBondedDescriptor : public FractionalDescriptor {
public:
    FcENBondedDescriptor() : FractionalDescriptor("FcENBonded", "Fraction of bonds where both atoms have EN > 3.0") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

// Radius-Based
class FcVdWMWDescriptor : public FractionalDescriptor {
public:
    FcVdWMWDescriptor() : FractionalDescriptor("FcVdWMW", "Sum(atom VdW volume × MW) / total MW") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class FcRcovMWDescriptor : public FractionalDescriptor {
public:
    FcRcovMWDescriptor() : FractionalDescriptor("FcRcovMW", "Sum(atom covalent radius × MW) / total MW") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

// Valence / Oxidation State
class FcHighOxStateDescriptor : public FractionalDescriptor {
public:
    FcHighOxStateDescriptor() : FractionalDescriptor("FcHighOxState", "Fraction of atoms with oxidation state ≥ +4") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

// Metal / Metalloid / Nonmetal Classification
class FcMetalloidDescriptor : public FractionalDescriptor {
public:
    FcMetalloidDescriptor() : FractionalDescriptor("FcMetalloid", "Fraction of metalloids") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

// Combined descriptors
class FcHETpolDescriptor : public FractionalDescriptor {
public:
    FcHETpolDescriptor() : FractionalDescriptor("FcHETpol", "Fraction of heteroatoms that are also polar") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class FcHALpolDescriptor : public FractionalDescriptor {
public:
    FcHALpolDescriptor() : FractionalDescriptor("FcHALpol", "Fraction of halogen atoms with polar bonds") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class FcHeavyPolDescriptor : public FractionalDescriptor {
public:
    FcHeavyPolDescriptor() : FractionalDescriptor("FcHeavyPol", "Fraction of heavy atoms (Z > 10) in polar bonds") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

// Additional Fractional Descriptors
class FcSp3HeavyAtomsDescriptor : public FractionalDescriptor {
public:
    FcSp3HeavyAtomsDescriptor() : FractionalDescriptor("FcSp3HeavyAtoms", "sp3 heavy atoms / total heavy atoms") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class FcSp2ENAboveAvgDescriptor : public FractionalDescriptor {
public:
    FcSp2ENAboveAvgDescriptor() : FractionalDescriptor("FcSp2ENAboveAvg", "sp2 atoms with EN > avg EN / total sp2 atoms") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

// class FcOxStateOddDescriptor : public FractionalDescriptor {
// public:
//     FcOxStateOddDescriptor() : FractionalDescriptor("FcOxStateOdd", "Atoms with odd ox state / total atoms") {}
//     std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
// };

class FcIEOddDescriptor : public FractionalDescriptor {
public:
    FcIEOddDescriptor() : FractionalDescriptor("FcIEOdd", "Atoms with odd IE / total atoms") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class FcEvenValenceAtomsDescriptor : public FractionalDescriptor {
public:
    FcEvenValenceAtomsDescriptor() : FractionalDescriptor("FcEvenValenceAtoms", "Even valence electron count / total atoms") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class FcNonHeteroHighEADescriptor : public FractionalDescriptor {
public:
    FcNonHeteroHighEADescriptor() : FractionalDescriptor("FcNonHeteroHighEA", "Non-heteroatoms with EA > 2 eV / total atoms") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class FcHeavyLowIEDescriptor : public FractionalDescriptor {
public:
    FcHeavyLowIEDescriptor() : FractionalDescriptor("FcHeavyLowIE", "Heavy atoms with IE < 7 / total heavy atoms") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class FcHeteroLowPolzDescriptor : public FractionalDescriptor {
public:
    FcHeteroLowPolzDescriptor() : FractionalDescriptor("FcHeteroLowPolz", "Heteroatoms with polarizability < 5 / total heteroatoms") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class FcOxENAboveThresholdDescriptor : public FractionalDescriptor {
public:
    FcOxENAboveThresholdDescriptor() : FractionalDescriptor("FcOxENAboveThreshold", "Atoms with (ox × EN) > 15 / total atoms") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class FcFormalChargeNonZeroDescriptor : public FractionalDescriptor {
public:
    FcFormalChargeNonZeroDescriptor() : FractionalDescriptor("FcFormalChargeNonZero", "Non-zero formal charge atoms / total atoms") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class FcFormalChargePositiveDescriptor : public FractionalDescriptor {
public:
    FcFormalChargePositiveDescriptor() : FractionalDescriptor("FcFormalChargePositive", "Atoms with formal charge > 0 / total atoms") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class FcFormalChargeNegativeDescriptor : public FractionalDescriptor {
public:
    FcFormalChargeNegativeDescriptor() : FractionalDescriptor("FcFormalChargeNegative", "Atoms with formal charge < 0 / total atoms") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

// class FcRadiusMWRatioAbove1Descriptor : public FractionalDescriptor {
// public:
//     FcRadiusMWRatioAbove1Descriptor() : FractionalDescriptor("FcRadiusMWRatioAbove1", "Atoms where radius/MW > 1 / total atoms") {}
//     std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
// };

class FcENMWAbove2Descriptor : public FractionalDescriptor {
public:
    FcENMWAbove2Descriptor() : FractionalDescriptor("FcENMWAbove2", "Atoms with EN × MW > 2 / total atoms") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class FcGroup16AtomsDescriptor : public FractionalDescriptor {
public:
    FcGroup16AtomsDescriptor() : FractionalDescriptor("FcGroup16Atoms", "Atoms in group 16 / total atoms") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class FcSpAtomsHighENDescriptor : public FractionalDescriptor {
public:
    FcSpAtomsHighENDescriptor() : FractionalDescriptor("FcSpAtomsHighEN", "sp-hybridized atoms with EN > 2.5 / total sp atoms") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class FcRingOxHighDescriptor : public FractionalDescriptor {
public:
    FcRingOxHighDescriptor() : FractionalDescriptor("FcRingOxHigh", "Atoms in ring with ox state > 2 / total ring atoms") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

// class FcLowRadicalENDescriptor : public FractionalDescriptor {
// public:
//     FcLowRadicalENDescriptor() : FractionalDescriptor("FcLowRadicalEN", "Radical atoms with EN < 2 / total radical atoms") {}
//     std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
// };

class FcHeavyFormalChargeDescriptor : public FractionalDescriptor {
public:
    FcHeavyFormalChargeDescriptor() : FractionalDescriptor("FcHeavyFormalCharge", "Heavy atoms with non-zero formal charge / heavy atom count") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class FcSp3PolarizabilityAbove10Descriptor : public FractionalDescriptor {
public:
    FcSp3PolarizabilityAbove10Descriptor() : FractionalDescriptor("FcSp3PolarizabilityAbove10", "sp3 atoms with polarizability > 10 / sp3 atom count") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class FcHeavyOxNegativeDescriptor : public FractionalDescriptor {
public:
    FcHeavyOxNegativeDescriptor() : FractionalDescriptor("FcHeavyOxNegative", "Heavy atoms with ox state < 0 / total heavy atoms") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

// class FcGroup17OxAbove1Descriptor : public FractionalDescriptor {
// public:
//     FcGroup17OxAbove1Descriptor() : FractionalDescriptor("FcGroup17OxAbove1", "Group 17 atoms with ox state > 1 / total group 17 atoms") {}
//     std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
// };

class FcENAboveMoleculeAvgDescriptor : public FractionalDescriptor {
public:
    FcENAboveMoleculeAvgDescriptor() : FractionalDescriptor("FcENAboveMoleculeAvg", "Atoms with EN > molecular average / total atoms") {}
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

}  // namespace descriptors
}  // namespace desfact