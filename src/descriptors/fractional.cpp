#include "descriptors/fractional.hpp"
#include "descriptors/element_properties.hpp"
#include <GraphMol/GraphMol.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/Atom.h>
#include <GraphMol/Bond.h>
#include <GraphMol/AtomIterators.h>
#include <GraphMol/BondIterators.h>
#include <GraphMol/Descriptors/MolDescriptors.h>
#include <GraphMol/Descriptors/Crippen.h>
#include <GraphMol/Descriptors/MolSurf.h>
#include <GraphMol/PeriodicTable.h>
#include <GraphMol/RingInfo.h> // Needed for ring ops
#include <GraphMol/Chirality.h> // Needed for chirality
#include <unordered_map>
#include <unordered_set>

namespace desfact {
namespace descriptors {

namespace {
    // ... other helper functions ...

    // Function to get the atomic mass of an atom
    double getAtomMass(const RDKit::Atom* atom) {
        if (!atom) return 0.0;
        return RDKit::PeriodicTable::getTable()->getAtomicWeight(atom->getAtomicNum());
    }

    // Function to get the average atomic mass of a molecule
    double getMoleculeAverageAtomicMass(const RDKit::ROMol* mol) {
        if (!mol || mol->getNumAtoms() == 0) return 0.0;
        double totalMass = 0.0;
        for (const RDKit::Atom* atom : mol->atoms()) {
            totalMass += getAtomMass(atom);
        }
        return totalMass / mol->getNumAtoms();
    }

    // Function to get the first ionization energy of an atom (in eV)
    double getAtomIonizationEnergy(const RDKit::Atom* atom) {
        if (!atom) return 0.0;
        return ElementProperties::getProperty(ElementProperties::ionizationEnergy, atom->getAtomicNum(), 10.0); // Default 10 eV
    }
    
    // Function to get the electron affinity of an atom (in eV)
    double getAtomElectronAffinity(const RDKit::Atom* atom) {
        if (!atom) return 0.0;
        return ElementProperties::getProperty(ElementProperties::electronAffinity, atom->getAtomicNum(), 1.0); // Default 1 eV
    }

    // Function to get the covalent radius of an atom (in Angstrom)
    double getAtomCovalentRadius(const RDKit::Atom* atom) {
        if (!atom) return 0.0;
        return ElementProperties::getProperty(ElementProperties::covalentRadius, atom->getAtomicNum(), 1.0); // Default 1 A
    }
    
    // Function to get the Van der Waals radius of an atom (in Angstrom)
    double getAtomVdWRadius(const RDKit::Atom* atom) {
        if (!atom) return 0.0;
        return ElementProperties::getProperty(ElementProperties::vanDerWaalsRadius, atom->getAtomicNum(), 1.6); // Default 1.6 A
    }

    // Helper to initialize ring info if needed
    void ensureRingInfo(const RDKit::ROMol& mol) {
        if (!mol.getRingInfo() || !mol.getRingInfo()->isInitialized()) {
            RDKit::MolOps::findSSSR(const_cast<RDKit::ROMol&>(mol));
        }
    }

} // End anonymous namespace

// FractionalDescriptor static methods
double FractionalDescriptor::calcAtomicFraction(const Molecule& mol, 
                                               const std::function<bool(const RDKit::Atom*)>& predicate) {
    if (!mol.isValid() || !mol.getMolecule()) {
        return 0.0; // Default to 0.0 if invalid
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    unsigned int matchCount = 0;
    unsigned int totalAtoms = rdkMol->getNumAtoms();
    
    if (totalAtoms == 0) return 0.0;
    
    for (const RDKit::Atom* atom : rdkMol->atoms()) {
        if (predicate(atom)) {
            matchCount++;
        }
    }
    
    return static_cast<double>(matchCount) / totalAtoms;
}

double FractionalDescriptor::calcBondFraction(const Molecule& mol, 
                                             const std::function<bool(const RDKit::Bond*)>& predicate) {
    if (!mol.isValid() || !mol.getMolecule()) {
        return 0.0; // Default to 0.0 if invalid
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    unsigned int matchCount = 0;
    unsigned int totalBonds = rdkMol->getNumBonds();
    
    if (totalBonds == 0) return 0.0;
    
    for (const RDKit::Bond* bond : rdkMol->bonds()) {
        if (predicate(bond)) {
            matchCount++;
        }
    }
    
    return static_cast<double>(matchCount) / totalBonds;
}

bool FractionalDescriptor::isElement(const RDKit::Atom* atom, int atomicNum) {
    return atom && atom->getAtomicNum() == atomicNum;
}

bool FractionalDescriptor::isHalogen(const RDKit::Atom* atom) {
    if (!atom) return false;
    int atomicNum = atom->getAtomicNum();
    return atomicNum == 9 || atomicNum == 17 || atomicNum == 35 || atomicNum == 53 || atomicNum == 85;
}

bool FractionalDescriptor::isHeteroatom(const RDKit::Atom* atom) {
    if (!atom) return false;
    int atomicNum = atom->getAtomicNum();
    return atomicNum != 6 && atomicNum != 1;  // Not C or H
}

bool FractionalDescriptor::isPolarAtom(const RDKit::Atom* atom) {
    if (!atom) return false;
    double en = getAtomElectronegativity(atom);
    return en > 2.55;  // Higher than carbon
}

bool FractionalDescriptor::isMetalAtom(const RDKit::Atom* atom) {
    if (!atom) return false;
    return ElementProperties::metals.find(atom->getAtomicNum()) != ElementProperties::metals.end();
}

bool FractionalDescriptor::isMetalloidAtom(const RDKit::Atom* atom) {
    if (!atom) return false;
    return ElementProperties::metalloids.find(atom->getAtomicNum()) != ElementProperties::metalloids.end();
}

double FractionalDescriptor::getAtomElectronegativity(const RDKit::Atom* atom) {
    if (!atom) return 0.0;
    return ElementProperties::getProperty(ElementProperties::electronegativity, atom->getAtomicNum(), 2.2);
}

double FractionalDescriptor::getMoleculeAverageElectronegativity(const RDKit::ROMol* mol) {
    if (!mol || mol->getNumAtoms() == 0) return 0.0;
    
    double sum = 0.0;
    for (const RDKit::Atom* atom : mol->atoms()) {
        sum += getAtomElectronegativity(atom);
    }
    
    return sum / mol->getNumAtoms();
}

double FractionalDescriptor::getAtomPolarizability(const RDKit::Atom* atom) {
    // Simple polarizability lookup based on atomic number
    static const std::unordered_map<int, double> polarizabilityMap = {
        {1, 0.667},   // H
        {6, 1.76},    // C
        {7, 1.10},    // N
        {8, 0.802},   // O
        {9, 0.557},   // F
        {15, 3.63},   // P
        {16, 2.90},   // S
        {17, 2.18},   // Cl
        {35, 3.05},   // Br
        {53, 5.35}    // I
    };
    
    int atomicNum = atom->getAtomicNum();
    auto it = polarizabilityMap.find(atomicNum);
    if (it != polarizabilityMap.end()) {
        return it->second;
    }
    
    // Default value for other atom types
    return 2.0;
}

// Implementation of descriptor calculation methods
std::variant<double, int, std::string> FcCDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) {
        return 0.0;
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    double totalMW = RDKit::Descriptors::calcExactMW(*rdkMol);
    
    if (totalMW <= 0.0) return 0.0;
    
    double carbonMW = 0.0;
    for (const RDKit::Atom* atom : rdkMol->atoms()) {
        if (atom->getAtomicNum() == 6) {  // Carbon
            carbonMW += RDKit::PeriodicTable::getTable()->getAtomicWeight(atom->getAtomicNum());
        }
    }
    
    return carbonMW / totalMW;
}

std::variant<double, int, std::string> FcFDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) {
        return 0.0;
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    double totalMW = RDKit::Descriptors::calcExactMW(*rdkMol);
    
    if (totalMW <= 0.0) return 0.0;
    
    double elementMW = 0.0;
    for (const RDKit::Atom* atom : rdkMol->atoms()) {
        if (atom->getAtomicNum() == 9) {  // Fluorine
            elementMW += RDKit::PeriodicTable::getTable()->getAtomicWeight(atom->getAtomicNum());
        }
    }
    
    return elementMW / totalMW;
}

std::variant<double, int, std::string> FcODescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) {
        return 0.0;
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    double totalMW = RDKit::Descriptors::calcExactMW(*rdkMol);
    
    if (totalMW <= 0.0) return 0.0;
    
    double elementMW = 0.0;
    for (const RDKit::Atom* atom : rdkMol->atoms()) {
        if (atom->getAtomicNum() == 8) {  // Oxygen
            elementMW += RDKit::PeriodicTable::getTable()->getAtomicWeight(atom->getAtomicNum());
        }
    }
    
    return elementMW / totalMW;
}

std::variant<double, int, std::string> FcClDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) {
        return 0.0;
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    double totalMW = RDKit::Descriptors::calcExactMW(*rdkMol);
    
    if (totalMW <= 0.0) return 0.0;
    
    double elementMW = 0.0;
    for (const RDKit::Atom* atom : rdkMol->atoms()) {
        if (atom->getAtomicNum() == 17) {  // Chlorine
            elementMW += RDKit::PeriodicTable::getTable()->getAtomicWeight(atom->getAtomicNum());
        }
    }
    
    return elementMW / totalMW;
}

std::variant<double, int, std::string> FcBrDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) {
        return 0.0;
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    double totalMW = RDKit::Descriptors::calcExactMW(*rdkMol);
    
    if (totalMW <= 0.0) return 0.0;
    
    double elementMW = 0.0;
    for (const RDKit::Atom* atom : rdkMol->atoms()) {
        if (atom->getAtomicNum() == 35) {  // Bromine
            elementMW += RDKit::PeriodicTable::getTable()->getAtomicWeight(atom->getAtomicNum());
        }
    }
    
    return elementMW / totalMW;
}

std::variant<double, int, std::string> FcIDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) {
        return 0.0;
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    double totalMW = RDKit::Descriptors::calcExactMW(*rdkMol);
    
    if (totalMW <= 0.0) return 0.0;
    
    double elementMW = 0.0;
    for (const RDKit::Atom* atom : rdkMol->atoms()) {
        if (atom->getAtomicNum() == 53) {  // Iodine
            elementMW += RDKit::PeriodicTable::getTable()->getAtomicWeight(atom->getAtomicNum());
        }
    }
    
    return elementMW / totalMW;
}

std::variant<double, int, std::string> FcSDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) {
        return 0.0;
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    double totalMW = RDKit::Descriptors::calcExactMW(*rdkMol);
    
    if (totalMW <= 0.0) return 0.0;
    
    double elementMW = 0.0;
    for (const RDKit::Atom* atom : rdkMol->atoms()) {
        if (atom->getAtomicNum() == 16) {  // Sulfur
            elementMW += RDKit::PeriodicTable::getTable()->getAtomicWeight(atom->getAtomicNum());
        }
    }
    
    return elementMW / totalMW;
}

std::variant<double, int, std::string> FcNDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) {
        return 0.0;
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    double totalMW = RDKit::Descriptors::calcExactMW(*rdkMol);
    
    if (totalMW <= 0.0) return 0.0;
    
    double elementMW = 0.0;
    for (const RDKit::Atom* atom : rdkMol->atoms()) {
        if (atom->getAtomicNum() == 7) {  // Nitrogen
            elementMW += RDKit::PeriodicTable::getTable()->getAtomicWeight(atom->getAtomicNum());
        }
    }
    
    return elementMW / totalMW;
}

std::variant<double, int, std::string> FcHDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return 0.0;
    
    double totalMW = 0.0;
    double hydrogenMW = 0.0;
    const double hAtomicWeight = 1.00794; // Precise hydrogen atomic weight
    
    // Count explicit hydrogens
    for (const RDKit::Atom* atom : rdkMol->atoms()) {
        double atomicWeight = RDKit::PeriodicTable::getTable()->getAtomicWeight(atom->getAtomicNum());
        totalMW += atomicWeight;
        
        if (atom->getAtomicNum() == 1) {
            hydrogenMW += atomicWeight;
        }
    }
    
    // Count implicit hydrogens
    for (const RDKit::Atom* atom : rdkMol->atoms()) {
        if (atom->getAtomicNum() > 1) { // Heavy atoms
            // Add implicit hydrogens
            hydrogenMW += atom->getTotalNumHs() * hAtomicWeight;
            totalMW += atom->getTotalNumHs() * hAtomicWeight;
        }
    }
    
    if (totalMW < 1e-6) return 0.0;
    
    return hydrogenMW / totalMW;
}

std::variant<double, int, std::string> FcHaloDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) {
        return 0.0;
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    double totalMW = RDKit::Descriptors::calcExactMW(*rdkMol);
    
    if (totalMW <= 0.0) return 0.0;
    
    double elementMW = 0.0;
    for (const RDKit::Atom* atom : rdkMol->atoms()) {
        int atomicNum = atom->getAtomicNum();
        if (atomicNum == 9 || atomicNum == 17 || atomicNum == 35 || atomicNum == 53) {  // F, Cl, Br, I
            elementMW += RDKit::PeriodicTable::getTable()->getAtomicWeight(atomicNum);
        }
    }
    
    return elementMW / totalMW;
}

std::variant<double, int, std::string> FcHeteroDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) {
        return 0.0;
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    double totalMW = RDKit::Descriptors::calcExactMW(*rdkMol);
    
    if (totalMW <= 0.0) return 0.0;
    
    double elementMW = 0.0;
    for (const RDKit::Atom* atom : rdkMol->atoms()) {
        int atomicNum = atom->getAtomicNum();
        if (atomicNum != 6 && atomicNum != 1) {  // Not C or H
            elementMW += RDKit::PeriodicTable::getTable()->getAtomicWeight(atomicNum);
        }
    }
    
    return elementMW / totalMW;
}

std::variant<double, int, std::string> FcPolarDescriptor::calculate(const Molecule& mol) const {
    return calcAtomicFraction(mol, [](const RDKit::Atom* atom) {
        return isPolarAtom(atom);
    });
}

std::variant<double, int, std::string> FcApolarDescriptor::calculate(const Molecule& mol) const {
    return calcAtomicFraction(mol, [](const RDKit::Atom* atom) {
        return !isPolarAtom(atom);
    });
}

std::variant<double, int, std::string> FcCSp3Descriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) {
        return 0.0;
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    unsigned int ccBondCount = 0;
    unsigned int ccSp3BondCount = 0;
    
    for (const RDKit::Bond* bond : rdkMol->bonds()) {
        const RDKit::Atom* beginAtom = bond->getBeginAtom();
        const RDKit::Atom* endAtom = bond->getEndAtom();
        
        if (beginAtom->getAtomicNum() == 6 && endAtom->getAtomicNum() == 6) {
            ccBondCount++;
            
            if (beginAtom->getHybridization() == RDKit::Atom::SP3 &&
                endAtom->getHybridization() == RDKit::Atom::SP3) {
                ccSp3BondCount++;
            }
        }
    }
    
    return (ccBondCount > 0) ? static_cast<double>(ccSp3BondCount) / ccBondCount : 0.0;
}

std::variant<double, int, std::string> FcCSp2Descriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) {
        return 0.0;
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    unsigned int ccBondCount = 0;
    unsigned int ccSp2BondCount = 0;
    
    for (const RDKit::Bond* bond : rdkMol->bonds()) {
        const RDKit::Atom* beginAtom = bond->getBeginAtom();
        const RDKit::Atom* endAtom = bond->getEndAtom();
        
        if (beginAtom->getAtomicNum() == 6 && endAtom->getAtomicNum() == 6) {
            ccBondCount++;
            
            if (beginAtom->getHybridization() == RDKit::Atom::SP2 &&
                endAtom->getHybridization() == RDKit::Atom::SP2) {
                ccSp2BondCount++;
            }
        }
    }
    
    return (ccBondCount > 0) ? static_cast<double>(ccSp2BondCount) / ccBondCount : 0.0;
}

std::variant<double, int, std::string> FcUnpolDescriptor::calculate(const Molecule& mol) const {
    return calcBondFraction(mol, [this](const RDKit::Bond* bond) {
        const RDKit::Atom* a1 = bond->getBeginAtom();
        const RDKit::Atom* a2 = bond->getEndAtom();
        double en1, en2;
        
        en1 = getAtomElectronegativity(a1);
        en2 = getAtomElectronegativity(a2);
        
        return std::abs(en1 - en2) < 0.2;  // Threshold for "equal" electronegativity
    });
}

std::variant<double, int, std::string> FcPolDescriptor::calculate(const Molecule& mol) const {
    return calcBondFraction(mol, [this](const RDKit::Bond* bond) {
        const RDKit::Atom* a1 = bond->getBeginAtom();
        const RDKit::Atom* a2 = bond->getEndAtom();
        double en1, en2;
        
        en1 = getAtomElectronegativity(a1);
        en2 = getAtomElectronegativity(a2);
        
        return std::abs(en1 - en2) >= 0.2;  // Threshold for "different" electronegativity
    });
}

std::variant<double, int, std::string> FcSumPolAtDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) {
        return 0.0;
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (rdkMol->getNumAtoms() == 0) return 0.0;
    
    double sumEN = 0.0;
    for (const RDKit::Atom* atom : rdkMol->atoms()) {
        sumEN += getAtomElectronegativity(atom);
    }
    
    return sumEN / rdkMol->getNumAtoms();
}

std::variant<double, int, std::string> FcSumPolMWDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) {
        return 0.0;
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    double totalMW = RDKit::Descriptors::calcExactMW(*rdkMol);
    
    if (totalMW <= 0.0) return 0.0;
    
    double sumEN = 0.0;
    for (const RDKit::Atom* atom : rdkMol->atoms()) {
        sumEN += getAtomElectronegativity(atom);
    }
    
    return sumEN / totalMW;
}

std::variant<double, int, std::string> FcBondAtDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) {
        return 0.0;
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (rdkMol->getNumAtoms() == 0) return 0.0;
    
    double totalBondOrder = 0.0;
    for (const RDKit::Bond* bond : rdkMol->bonds()) {
        totalBondOrder += bond->getBondTypeAsDouble();
    }
    
    return totalBondOrder / rdkMol->getNumAtoms();
}

std::variant<double, int, std::string> FcBondNDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) {
        return 0.0;
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    
    // Count non-single bonds on Nitrogen atoms
    int nBondCount = 0;
    int nNonSingleBondCount = 0;
    
    for (const RDKit::Atom* atom : rdkMol->atoms()) {
        if (atom->getAtomicNum() == 7) {  // Nitrogen
            RDKit::ROMol::ADJ_ITER nbrIdx, endNbrs;
            boost::tie(nbrIdx, endNbrs) = rdkMol->getAtomNeighbors(atom);
            
            while (nbrIdx != endNbrs) {
                const RDKit::Bond* bond = rdkMol->getBondBetweenAtoms(atom->getIdx(), *nbrIdx);
                nBondCount++;
                
                if (bond->getBondType() != RDKit::Bond::SINGLE) {
                    nNonSingleBondCount++;
                }
                
                ++nbrIdx;
            }
        }
    }
    
    return (nBondCount > 0) ? static_cast<double>(nNonSingleBondCount) / nBondCount : 0.0;
}

std::variant<double, int, std::string> FcBondODescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) {
        return 0.0;
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    
    // Count non-single bonds on Oxygen atoms
    int oBondCount = 0;
    int oNonSingleBondCount = 0;
    
    for (const RDKit::Atom* atom : rdkMol->atoms()) {
        if (atom->getAtomicNum() == 8) {  // Oxygen
            RDKit::ROMol::ADJ_ITER nbrIdx, endNbrs;
            boost::tie(nbrIdx, endNbrs) = rdkMol->getAtomNeighbors(atom);
            
            while (nbrIdx != endNbrs) {
                const RDKit::Bond* bond = rdkMol->getBondBetweenAtoms(atom->getIdx(), *nbrIdx);
                oBondCount++;
                
                if (bond->getBondType() != RDKit::Bond::SINGLE) {
                    oNonSingleBondCount++;
                }
                
                ++nbrIdx;
            }
        }
    }
    
    return (oBondCount > 0) ? static_cast<double>(oNonSingleBondCount) / oBondCount : 0.0;
}

std::variant<double, int, std::string> FcBondCDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) {
        return 0.0;
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    
    // Count non-single bonds on Carbon atoms
    int cBondCount = 0;
    int cNonSingleBondCount = 0;
    
    for (const RDKit::Atom* atom : rdkMol->atoms()) {
        if (atom->getAtomicNum() == 6) {  // Carbon
            RDKit::ROMol::ADJ_ITER nbrIdx, endNbrs;
            boost::tie(nbrIdx, endNbrs) = rdkMol->getAtomNeighbors(atom);
            
            while (nbrIdx != endNbrs) {
                const RDKit::Bond* bond = rdkMol->getBondBetweenAtoms(atom->getIdx(), *nbrIdx);
                cBondCount++;
                
                if (bond->getBondType() != RDKit::Bond::SINGLE) {
                    cNonSingleBondCount++;
                }
                
                ++nbrIdx;
            }
        }
    }
    
    return (cBondCount > 0) ? static_cast<double>(cNonSingleBondCount) / cBondCount : 0.0;
}

std::variant<double, int, std::string> FcBondSDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) {
        return 0.0;
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    
    // Count non-single bonds on Sulfur atoms
    int sBondCount = 0;
    int sNonSingleBondCount = 0;
    
    for (const RDKit::Atom* atom : rdkMol->atoms()) {
        if (atom->getAtomicNum() == 16) {  // Sulfur
            RDKit::ROMol::ADJ_ITER nbrIdx, endNbrs;
            boost::tie(nbrIdx, endNbrs) = rdkMol->getAtomNeighbors(atom);
            
            while (nbrIdx != endNbrs) {
                const RDKit::Bond* bond = rdkMol->getBondBetweenAtoms(atom->getIdx(), *nbrIdx);
                sBondCount++;
                
                if (bond->getBondType() != RDKit::Bond::SINGLE) {
                    sNonSingleBondCount++;
                }
                
                ++nbrIdx;
            }
        }
    }
    
    return (sBondCount > 0) ? static_cast<double>(sNonSingleBondCount) / sBondCount : 0.0;
}

std::variant<double, int, std::string> FcBondPDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) {
        return 0.0;
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    
    // Count non-single bonds on Phosphorus atoms
    int pBondCount = 0;
    int pNonSingleBondCount = 0;
    
    for (const RDKit::Atom* atom : rdkMol->atoms()) {
        if (atom->getAtomicNum() == 15) {  // Phosphorus
            RDKit::ROMol::ADJ_ITER nbrIdx, endNbrs;
            boost::tie(nbrIdx, endNbrs) = rdkMol->getAtomNeighbors(atom);
            
            while (nbrIdx != endNbrs) {
                const RDKit::Bond* bond = rdkMol->getBondBetweenAtoms(atom->getIdx(), *nbrIdx);
                pBondCount++;
                
                if (bond->getBondType() != RDKit::Bond::SINGLE) {
                    pNonSingleBondCount++;
                }
                
                ++nbrIdx;
            }
        }
    }
    
    return (pBondCount > 0) ? static_cast<double>(pNonSingleBondCount) / pBondCount : 0.0;
}

std::variant<double, int, std::string> FcHBDonorsDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) {
        return 0.0;
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (rdkMol->getNumAtoms() == 0) return 0.0;
    
    int hbdCount = RDKit::Descriptors::calcNumHBD(*rdkMol);
    return static_cast<double>(hbdCount) / rdkMol->getNumAtoms();
}

std::variant<double, int, std::string> FcHBAcceptorsDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) {
        return 0.0;
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (rdkMol->getNumAtoms() == 0) return 0.0;
    
    int hbaCount = RDKit::Descriptors::calcNumHBA(*rdkMol);
    return static_cast<double>(hbaCount) / rdkMol->getNumAtoms();
}

std::variant<double, int, std::string> FcAromaticAtomsDescriptor::calculate(const Molecule& mol) const {
    return calcAtomicFraction(mol, [](const RDKit::Atom* atom) {
        return atom->getIsAromatic();
    });
}

std::variant<double, int, std::string> FcRingAtomsDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) {
        return 0.0;
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (rdkMol->getNumAtoms() == 0) return 0.0;
    
    RDKit::RingInfo* ringInfo = rdkMol->getRingInfo();
    if (!ringInfo->isInitialized()) {
        RDKit::MolOps::findSSSR(*rdkMol);
    }
    
    int ringAtomCount = 0;
    for (const RDKit::Atom* atom : rdkMol->atoms()) {
        if (ringInfo->numAtomRings(atom->getIdx()) > 0) {
            ringAtomCount++;
        }
    }
    
    return static_cast<double>(ringAtomCount) / rdkMol->getNumAtoms();
}

std::variant<double, int, std::string> FcBridgeAtomsDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) {
        return 0.0;
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (rdkMol->getNumAtoms() == 0) return 0.0;
    
    RDKit::RingInfo* ringInfo = rdkMol->getRingInfo();
    if (!ringInfo->isInitialized()) {
        RDKit::MolOps::findSSSR(*rdkMol);
    }
    
    int bridgeAtomCount = 0;
    for (const RDKit::Atom* atom : rdkMol->atoms()) {
        if (ringInfo->numAtomRings(atom->getIdx()) >= 2) {
            bridgeAtomCount++;
        }
    }
    
    return static_cast<double>(bridgeAtomCount) / rdkMol->getNumAtoms();
}

std::variant<double, int, std::string> FcChargedAtomsDescriptor::calculate(const Molecule& mol) const {
    return calcAtomicFraction(mol, [](const RDKit::Atom* atom) {
        return atom->getFormalCharge() != 0;
    });
}

std::variant<double, int, std::string> FcHeavyAtomsDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return 0.0;
    
    int totalAtoms = 0;
    int heavyAtoms = 0;
    
    // Count both explicit atoms and implicit hydrogens
    for (const RDKit::Atom* atom : rdkMol->atoms()) {
        totalAtoms++;
        
        if (atom->getAtomicNum() > 1) {
            heavyAtoms++;
        }
        
        // Add implicit hydrogens to totalAtoms count
        if (atom->getAtomicNum() > 1) {
            totalAtoms += atom->getTotalNumHs();
        }
    }
    
    if (totalAtoms == 0) return 0.0;
    
    return static_cast<double>(heavyAtoms) / totalAtoms;
}

std::variant<double, int, std::string> FcMetalsDescriptor::calculate(const Molecule& mol) const {
    return calcAtomicFraction(mol, [this](const RDKit::Atom* atom) {
        return isMetalAtom(atom);
    });
}

std::variant<double, int, std::string> FcENAboveAvgDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) {
        return 0.0;
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (rdkMol->getNumAtoms() == 0) return 0.0;
    
    double avgEN = getMoleculeAverageElectronegativity(rdkMol);
    
    return calcAtomicFraction(mol, [avgEN, this](const RDKit::Atom* atom) {
        return getAtomElectronegativity(atom) > avgEN;
    });
}

std::variant<double, int, std::string> FcENBelowAvgDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) {
        return 0.0;
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (rdkMol->getNumAtoms() == 0) return 0.0;
    
    double avgEN = getMoleculeAverageElectronegativity(rdkMol);
    
    return calcAtomicFraction(mol, [avgEN, this](const RDKit::Atom* atom) {
        return getAtomElectronegativity(atom) < avgEN;
    });
}

std::variant<double, int, std::string> FcENHighDescriptor::calculate(const Molecule& mol) const {
    return calcAtomicFraction(mol, [this](const RDKit::Atom* atom) {
        return getAtomElectronegativity(atom) > 3.5;  // High electronegativity (e.g., O, F)
    });
}

std::variant<double, int, std::string> FcENLowDescriptor::calculate(const Molecule& mol) const {
    return calcAtomicFraction(mol, [this](const RDKit::Atom* atom) {
        return getAtomElectronegativity(atom) < 2.3;  // Adjusted threshold from 2.0 (captures H, P)
    });
}

std::variant<double, int, std::string> FcSmallRDescriptor::calculate(const Molecule& mol) const {
    return calcAtomicFraction(mol, [](const RDKit::Atom* atom) {
        double radius = ElementProperties::getProperty(ElementProperties::covalentRadius, atom->getAtomicNum(), 1.0);
        return radius < 0.77;  // Small radius (e.g., H, F)
    });
}

std::variant<double, int, std::string> FcLargeRDescriptor::calculate(const Molecule& mol) const {
    return calcAtomicFraction(mol, [](const RDKit::Atom* atom) {
        double radius = ElementProperties::getProperty(ElementProperties::covalentRadius, atom->getAtomicNum(), 1.0);
        return radius > 1.1;  // Adjusted threshold from 1.4 (captures P, Br, I)
    });
}

std::variant<double, int, std::string> FcLowPolzDescriptor::calculate(const Molecule& mol) const {
    return calcAtomicFraction(mol, [](const RDKit::Atom* atom) {
        double polarizability = ElementProperties::getProperty(ElementProperties::atomicPolarizability, atom->getAtomicNum(), 1.0);
        return polarizability < 5.0;  // Keep threshold 5.0 (excludes I)
    });
}

std::variant<double, int, std::string> FcHighPolzDescriptor::calculate(const Molecule& mol) const {
    return calcAtomicFraction(mol, [](const RDKit::Atom* atom) {
        double polarizability = ElementProperties::getProperty(ElementProperties::atomicPolarizability, atom->getAtomicNum(), 1.0);
        return polarizability > 3.5;  // Adjusted threshold from 15.0 (captures P, I)
    });
}

std::variant<double, int, std::string> FcHighEADescriptor::calculate(const Molecule& mol) const {
    // Use the static helper from the base class or ElementProperties namespace
    return calcAtomicFraction(mol, [](const RDKit::Atom* atom) {
        return ElementProperties::getProperty(ElementProperties::electronAffinity, atom->getAtomicNum(), 1.0) > 2.0;
    });
}

std::variant<double, int, std::string> FcLowEADescriptor::calculate(const Molecule& mol) const {
    return calcAtomicFraction(mol, [](const RDKit::Atom* atom) {
        double ea = ElementProperties::getProperty(ElementProperties::electronAffinity, atom->getAtomicNum(), 1.0);
        return ea < 0.5;  // Low electron affinity
    });
}

std::variant<double, int, std::string> FcSmallVdWDescriptor::calculate(const Molecule& mol) const {
    return calcAtomicFraction(mol, [](const RDKit::Atom* atom) {
        double vdw = ElementProperties::getProperty(ElementProperties::vanDerWaalsVolume, atom->getAtomicNum(), 20.0);
        return vdw < 15.0;  // Small VdW volume
    });
}

std::variant<double, int, std::string> FcLargeVdWDescriptor::calculate(const Molecule& mol) const {
    return calcAtomicFraction(mol, [](const RDKit::Atom* atom) {
        double vdw = ElementProperties::getProperty(ElementProperties::vanDerWaalsVolume, atom->getAtomicNum(), 20.0);
        return vdw > 25.0;  // Large VdW volume
    });
}

std::variant<double, int, std::string> FcENBondedDescriptor::calculate(const Molecule& mol) const {
    // calcBondFraction already returns 0.0 for invalid mol
    return calcBondFraction(mol, [this](const RDKit::Bond* bond) {
        double en1 = getAtomElectronegativity(bond->getBeginAtom());
        double en2 = getAtomElectronegativity(bond->getEndAtom());
        return en1 > 3.0 && en2 > 3.0;  // Both atoms have high electronegativity
    });
}

std::variant<double, int, std::string> FcVdWMWDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) {
        return 0.0;
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    double totalMW = RDKit::Descriptors::calcExactMW(*rdkMol);
    
    if (totalMW <= 0.0) return 0.0;
    
    double vdwMwSum = 0.0;
    for (const RDKit::Atom* atom : rdkMol->atoms()) {
        double vdw = ElementProperties::getProperty(ElementProperties::vanDerWaalsVolume, atom->getAtomicNum(), 20.0);
        double atomicWeight = RDKit::PeriodicTable::getTable()->getAtomicWeight(atom->getAtomicNum());
        vdwMwSum += vdw * atomicWeight;
    }
    
    return vdwMwSum / totalMW;
}

std::variant<double, int, std::string> FcRcovMWDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) {
        return 0.0;
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    double totalMW = RDKit::Descriptors::calcExactMW(*rdkMol);
    
    if (totalMW <= 0.0) return 0.0;
    
    double rcovMwSum = 0.0;
    for (const RDKit::Atom* atom : rdkMol->atoms()) {
        double rcov = ElementProperties::getProperty(ElementProperties::covalentRadius, atom->getAtomicNum(), 1.0);
        double atomicWeight = RDKit::PeriodicTable::getTable()->getAtomicWeight(atom->getAtomicNum());
        rcovMwSum += rcov * atomicWeight;
    }
    
    return rcovMwSum / totalMW;
}

std::variant<double, int, std::string> FcHighOxStateDescriptor::calculate(const Molecule& mol) const {
    // calcAtomicFraction already returns 0.0 for invalid mol
    return calcAtomicFraction(mol, [](const RDKit::Atom* atom) {
        int oxState = 0;
        try {
            oxState = atom->getProp<int>(RDKit::common_properties::_CIPRank); // Using CIPRank as a proxy, RDKit default ox state not readily available
        } catch (...) { /* ignore */ }
        return oxState >= 4;
    });
}

std::variant<double, int, std::string> FcMetalloidDescriptor::calculate(const Molecule& mol) const {
    // calcAtomicFraction already returns 0.0 for invalid mol
    return calcAtomicFraction(mol, [this](const RDKit::Atom* atom) {
        return isMetalloidAtom(atom);
    });
}

std::variant<double, int, std::string> FcHETpolDescriptor::calculate(const Molecule& mol) const {
    // calcAtomicFraction already returns 0.0 for invalid mol
    return calcAtomicFraction(mol, [this](const RDKit::Atom* atom) {
        return isHeteroatom(atom) && isPolarAtom(atom);
    });
}

std::variant<double, int, std::string> FcHALpolDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) {
        return 0.0;
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    unsigned int halogenCount = 0;
    unsigned int polarHalogenCount = 0;
    
    for (const RDKit::Atom* atom : rdkMol->atoms()) {
        if (isHalogen(atom)) {
            halogenCount++;
            
            // Check if this halogen has a polar bond
            RDKit::ROMol::ADJ_ITER nbrIdx, endNbrs;
            boost::tie(nbrIdx, endNbrs) = rdkMol->getAtomNeighbors(atom);
            
            while (nbrIdx != endNbrs) {
                const RDKit::Atom* neighbor = rdkMol->getAtomWithIdx(*nbrIdx);
                double en1 = getAtomElectronegativity(atom);
                double en2 = getAtomElectronegativity(neighbor);
                
                if (std::abs(en1 - en2) >= 0.5) {  // Consider a more significant threshold for polarity
                    polarHalogenCount++;
                    break;  // Count the atom once if any of its bonds is polar
                }
                
                ++nbrIdx;
            }
        }
    }
    
    return (halogenCount > 0) ? static_cast<double>(polarHalogenCount) / halogenCount : 0.0;
}

std::variant<double, int, std::string> FcHeavyPolDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) {
        return 0.0;
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    unsigned int heavyAtomCount = 0;
    unsigned int polarHeavyAtomCount = 0;
    
    for (const RDKit::Atom* atom : rdkMol->atoms()) {
        if (atom->getAtomicNum() > 10) {  // Heavy atoms (Z > 10)
            heavyAtomCount++;
            
            // Check if this heavy atom has a polar bond
            bool hasPolarBond = false;
            RDKit::ROMol::ADJ_ITER nbrIdx, endNbrs;
            boost::tie(nbrIdx, endNbrs) = rdkMol->getAtomNeighbors(atom);
            
            while (nbrIdx != endNbrs) {
                const RDKit::Atom* neighbor = rdkMol->getAtomWithIdx(*nbrIdx);
                double en1 = getAtomElectronegativity(atom);
                double en2 = getAtomElectronegativity(neighbor);
                
                if (std::abs(en1 - en2) >= 0.5) {
                    hasPolarBond = true;
                    break;
                }
                
                ++nbrIdx;
            }
            
            if (hasPolarBond) {
                polarHeavyAtomCount++;
            }
        }
    }
    
    return (heavyAtomCount > 0) ? static_cast<double>(polarHeavyAtomCount) / heavyAtomCount : 0.0;
}

// FcSp3PolarizabilityAbove10 - sp3 atoms with polarizability > 10 / sp3 atom count
std::variant<double, int, std::string> FcSp3PolarizabilityAbove10Descriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) {
        return 0.0;
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    int sp3Count = 0;
    int sp3HighPolzCount = 0;
    
    for (const RDKit::Atom* atom : rdkMol->atoms()) {
        if (atom->getHybridization() == RDKit::Atom::SP3) {
            sp3Count++;
            // Use ElementProperties map
            double polz = ElementProperties::getProperty(ElementProperties::atomicPolarizability, atom->getAtomicNum(), 1.0);
            
            if (polz > 3.5) { // Adjusted threshold from 10.0 (captures P, I)
                sp3HighPolzCount++;
            }
        }
    }
    
    return (sp3Count > 0) ? static_cast<double>(sp3HighPolzCount) / sp3Count : 0.0;
}

// FcHeavyOxNegative - Heavy atoms with ox state < 0 / total heavy atoms
std::variant<double, int, std::string> FcHeavyOxNegativeDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) {
        return 0.0;
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    int heavyCount = 0;
    int heavyNegOxCount = 0;
    
    for (const RDKit::Atom* atom : rdkMol->atoms()) {
        if (atom->getAtomicNum() > 1) {  // Heavy atom
            heavyCount++;
            
            // Approximate oxidation state
            int atomicNum = atom->getAtomicNum();
            int valence = atom->getTotalValence();
            int formalCharge = atom->getFormalCharge();
            int oxidationState = 0;
            
            // Rough oxidation state estimation
            if (atomicNum == 7) oxidationState = -3 + valence + formalCharge; // N
            else if (atomicNum == 8) oxidationState = -2 + valence + formalCharge; // O
            else if (atomicNum == 16) oxidationState = -2 + valence + formalCharge; // S
            else if (atomicNum == 9) oxidationState = -1 + valence + formalCharge; // F
            else if (atomicNum == 17) oxidationState = -1 + valence + formalCharge; // Cl
            else if (atomicNum == 35) oxidationState = -1 + valence + formalCharge; // Br
            else if (atomicNum == 53) oxidationState = -1 + valence + formalCharge; // I
            else if (atomicNum == 6) oxidationState = 4 - valence + formalCharge; // C
            else if (atomicNum == 15) oxidationState = 5 - valence + formalCharge; // P
            
            if (oxidationState < 0) {
                heavyNegOxCount++;
            }
        }
    }
    
    return (heavyCount > 0) ? static_cast<double>(heavyNegOxCount) / heavyCount : 0.0;
}


// FcENAboveMoleculeAvg - Atoms with EN > molecular average / total atoms
std::variant<double, int, std::string> FcENAboveMoleculeAvgDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) {
        return 0.0;
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    double avgEN = getMoleculeAverageElectronegativity(rdkMol);
    int atomCount = rdkMol->getNumAtoms();
    int highENCount = 0;
    
    for (const RDKit::Atom* atom : rdkMol->atoms()) {
        double en = getAtomElectronegativity(atom);
        
        if (en > avgEN) {
            highENCount++;
        }
    }
    
    return (atomCount > 0) ? static_cast<double>(highENCount) / atomCount : 0.0;
}

// FcSp2ENAboveAvgDescriptor
std::variant<double, int, std::string> FcSp2ENAboveAvgDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) {
        return 0.0;
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    double avgEN = getMoleculeAverageElectronegativity(rdkMol);
    int sp2Count = 0;
    int sp2HighENCount = 0;
    
    for (const RDKit::Atom* atom : rdkMol->atoms()) {
        if (atom->getHybridization() == RDKit::Atom::SP2) {
            sp2Count++;
            if (getAtomElectronegativity(atom) > avgEN) {
                sp2HighENCount++;
            }
        }
    }
    
    return (sp2Count > 0) ? static_cast<double>(sp2HighENCount) / sp2Count : 0.0;
}

// FcENMWAbove2Descriptor
std::variant<double, int, std::string> FcENMWAbove2Descriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) {
        return 0.0;
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    int atomCount = rdkMol->getNumAtoms();
    int highENMWCount = 0;
    
    for (const RDKit::Atom* atom : rdkMol->atoms()) {
        double en = getAtomElectronegativity(atom);
        double mw = RDKit::PeriodicTable::getTable()->getAtomicWeight(atom->getAtomicNum());
        
        if (en * mw > 40.0) { // Adjusted threshold from 2.0 (captures N, O, F, P, S, Cl, Br, I)
            highENMWCount++;
        }
    }
    
    return (atomCount > 0) ? static_cast<double>(highENMWCount) / atomCount : 0.0;
}

// FcFormalChargePositiveDescriptor
std::variant<double, int, std::string> FcFormalChargePositiveDescriptor::calculate(const Molecule& mol) const {
    return calcAtomicFraction(mol, [](const RDKit::Atom* atom) {
        return atom->getFormalCharge() > 0;
    });
}


// FcGroup16AtomsDescriptor
std::variant<double, int, std::string> FcGroup16AtomsDescriptor::calculate(const Molecule& mol) const {
    return calcAtomicFraction(mol, [](const RDKit::Atom* atom) {
        int n = atom->getAtomicNum();
        return (n == 8 || n == 16 || n == 34 || n == 52); // Group 16
    });
}


// FcIEOddDescriptor
std::variant<double, int, std::string> FcIEOddDescriptor::calculate(const Molecule& mol) const {
    return calcAtomicFraction(mol, [this](const RDKit::Atom* atom) {
        double ie = getAtomElectronegativity(atom);
        return static_cast<int>(ie) % 2 == 1;  // Rough approximation for odd IE
    });
}


// FcSpAtomsHighENDescriptor
std::variant<double, int, std::string> FcSpAtomsHighENDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0; // Changed NaN to 0.0
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return 0.0; // Changed NaN to 0.0
    int spAtomCount = 0;
    int spHighENCount = 0;
    for (const auto& atom : rdkMol->atoms()) {
        if (atom->getHybridization() == RDKit::Atom::HybridizationType::SP) {
            spAtomCount++;
            if (getAtomElectronegativity(atom) > 2.5) {
                spHighENCount++;
            }
        }
    }
    return (spAtomCount > 0) ? static_cast<double>(spHighENCount) / spAtomCount : 0.0;
}

// FcFormalChargeNegativeDescriptor
std::variant<double, int, std::string> FcFormalChargeNegativeDescriptor::calculate(const Molecule& mol) const {
    return calcAtomicFraction(mol, [](const RDKit::Atom* atom) {
        return atom->getFormalCharge() < 0;
    });
}

// FcSp3HeavyAtomsDescriptor
std::variant<double, int, std::string> FcSp3HeavyAtomsDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) {
        return 0.0;
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    int heavyCount = 0;
    int sp3HeavyCount = 0;
    
    for (const RDKit::Atom* atom : rdkMol->atoms()) {
        if (atom->getAtomicNum() > 1) {  // Heavy atom
            heavyCount++;
            if (atom->getHybridization() == RDKit::Atom::SP3) {
                sp3HeavyCount++;
            }
        }
    }
    
    return (heavyCount > 0) ? static_cast<double>(sp3HeavyCount) / heavyCount : 0.0;
}

// FcHeteroLowPolzDescriptor
std::variant<double, int, std::string> FcHeteroLowPolzDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0; // Changed NaN to 0.0
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return 0.0; // Changed NaN to 0.0
    int heteroAtomCount = 0;
    int heteroLowPolzCount = 0;
    for (const auto& atom : rdkMol->atoms()) {
        if (isHeteroatom(atom)) {
            heteroAtomCount++;
            if (getAtomPolarizability(atom) < 5.0) {
                heteroLowPolzCount++;
            }
        }
    }
    return (heteroAtomCount > 0) ? static_cast<double>(heteroLowPolzCount) / heteroAtomCount : 0.0;
}

// FcRingOxHighDescriptor
std::variant<double, int, std::string> FcRingOxHighDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) {
        return 0.0;
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    RDKit::RingInfo* ringInfo = rdkMol->getRingInfo();
    
    if (!ringInfo->isInitialized()) {
        RDKit::MolOps::findSSSR(*rdkMol);
    }
    
    int ringAtomCount = 0;
    int ringHighOxCount = 0;
    
    for (const RDKit::Atom* atom : rdkMol->atoms()) {
        if (ringInfo->numAtomRings(atom->getIdx()) > 0) {  // Atom is in a ring
            ringAtomCount++;
            
            // Approximate oxidation state
            int atomicNum = atom->getAtomicNum();
            int valence = atom->getTotalValence();
            int formalCharge = atom->getFormalCharge();
            int oxidationState = 0;
            
            // Rough oxidation state estimation
            if (atomicNum == 7) oxidationState = -3 + valence + formalCharge; // N
            else if (atomicNum == 8) oxidationState = -2 + valence + formalCharge; // O
            else if (atomicNum == 16) oxidationState = -2 + valence + formalCharge; // S
            else if (atomicNum == 9) oxidationState = -1 + valence + formalCharge; // F
            else if (atomicNum == 17) oxidationState = -1 + valence + formalCharge; // Cl
            else if (atomicNum == 35) oxidationState = -1 + valence + formalCharge; // Br
            else if (atomicNum == 53) oxidationState = -1 + valence + formalCharge; // I
            else if (atomicNum == 6) oxidationState = 4 - valence + formalCharge; // C
            else if (atomicNum == 15) oxidationState = 5 - valence + formalCharge; // P
            
            if (oxidationState > 2) {
                ringHighOxCount++;
            }
        }
    }
    
    return (ringAtomCount > 0) ? static_cast<double>(ringHighOxCount) / ringAtomCount : 0.0;
}

// FcOxENAboveThresholdDescriptor
std::variant<double, int, std::string> FcOxENAboveThresholdDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) {
        return 0.0;
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    int atomCount = rdkMol->getNumAtoms();
    int highOxENCount = 0;
    
    for (const RDKit::Atom* atom : rdkMol->atoms()) {
        // Approximate oxidation state
        int atomicNum = atom->getAtomicNum();
        int valence = atom->getTotalValence();
        int formalCharge = atom->getFormalCharge();
        int oxidationState = 0;
        
        // Rough oxidation state estimation
        if (atomicNum == 7) oxidationState = -3 + valence + formalCharge; // N
        else if (atomicNum == 8) oxidationState = -2 + valence + formalCharge; // O
        else if (atomicNum == 16) oxidationState = -2 + valence + formalCharge; // S
        else if (atomicNum == 9) oxidationState = -1 + valence + formalCharge; // F
        else if (atomicNum == 17) oxidationState = -1 + valence + formalCharge; // Cl
        else if (atomicNum == 35) oxidationState = -1 + valence + formalCharge; // Br
        else if (atomicNum == 53) oxidationState = -1 + valence + formalCharge; // I
        else if (atomicNum == 6) oxidationState = 4 - valence + formalCharge; // C
        else if (atomicNum == 15) oxidationState = 5 - valence + formalCharge; // P
        
        double en = getAtomElectronegativity(atom);
        
        if (std::abs(oxidationState) * en > 10.0) { // Adjusted threshold from 15.0
            highOxENCount++;
        }
    }
    
    return (atomCount > 0) ? static_cast<double>(highOxENCount) / atomCount : 0.0;
}

// FcFormalChargeNonZeroDescriptor
std::variant<double, int, std::string> FcFormalChargeNonZeroDescriptor::calculate(const Molecule& mol) const {
    return calcAtomicFraction(mol, [](const RDKit::Atom* atom) {
        return atom->getFormalCharge() != 0;
    });
}

// FcEvenValenceAtomsDescriptor
std::variant<double, int, std::string> FcEvenValenceAtomsDescriptor::calculate(const Molecule& mol) const {
    return calcAtomicFraction(mol, [](const RDKit::Atom* atom) {
        int valence = atom->getTotalValence();
        return valence % 2 == 0;  // Even valence
    });
}

// FcNonHeteroHighEADescriptor
std::variant<double, int, std::string> FcNonHeteroHighEADescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) {
        return 0.0;
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    int atomCount = rdkMol->getNumAtoms();
    int nonHeteroHighEACount = 0;
    
    for (const RDKit::Atom* atom : rdkMol->atoms()) {
        if (!isHeteroatom(atom)) {  // Not a heteroatom (C or H)
            // Use ElementProperties map
            int atomicNum = atom->getAtomicNum();
            double ea = ElementProperties::getProperty(ElementProperties::electronAffinity, atomicNum, 1.0);
            
            if (ea > 1.0) { // Adjusted threshold from 2.0 (captures C)
                nonHeteroHighEACount++;
            }
        }
    }
    
    return (atomCount > 0) ? static_cast<double>(nonHeteroHighEACount) / atomCount : 0.0;
}

// FcHeavyFormalChargeDescriptor
std::variant<double, int, std::string> FcHeavyFormalChargeDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) {
        return 0.0;
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    int heavyCount = 0;
    int heavyChargedCount = 0;
    
    for (const RDKit::Atom* atom : rdkMol->atoms()) {
        if (atom->getAtomicNum() > 1) {  // Heavy atom
            heavyCount++;
            if (atom->getFormalCharge() != 0) {
                heavyChargedCount++;
            }
        }
    }
    
    return (heavyCount > 0) ? static_cast<double>(heavyChargedCount) / heavyCount : 0.0;
}

// --- New Fractional Descriptor Implementations ---

std::variant<double, int, std::string> FcAtomMassAboveAvg::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (rdkMol->getNumAtoms() == 0) return 0.0;
    double avgMass = getMoleculeAverageAtomicMass(rdkMol);
    return calcAtomicFraction(mol, [avgMass](const RDKit::Atom* atom) {
        return getAtomMass(atom) > avgMass;
    });
}

std::variant<double, int, std::string> FcAtomMassBelowAvg::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (rdkMol->getNumAtoms() == 0) return 0.0;
    double avgMass = getMoleculeAverageAtomicMass(rdkMol);
    return calcAtomicFraction(mol, [avgMass](const RDKit::Atom* atom) {
        return getAtomMass(atom) < avgMass;
    });
}

std::variant<double, int, std::string> FcHighIE::calculate(const Molecule& mol) const {
    return calcAtomicFraction(mol, [](const RDKit::Atom* atom) {
        return getAtomIonizationEnergy(atom) > 13.0;
    });
}

std::variant<double, int, std::string> FcLowIE::calculate(const Molecule& mol) const {
    return calcAtomicFraction(mol, [](const RDKit::Atom* atom) {
        return getAtomIonizationEnergy(atom) < 9.0;
    });
}

std::variant<double, int, std::string> FcMediumPolz::calculate(const Molecule& mol) const {
    return calcAtomicFraction(mol, [](const RDKit::Atom* atom) {
        double polz = getAtomPolarizability(atom);
        return polz >= 1.0 && polz <= 3.0;
    });
}

std::variant<double, int, std::string> FcNearZeroEA::calculate(const Molecule& mol) const {
    return calcAtomicFraction(mol, [](const RDKit::Atom* atom) {
        double ea = getAtomElectronAffinity(atom);
        return ea >= -0.5 && ea <= 0.5;
    });
}

std::variant<double, int, std::string> FcValence1::calculate(const Molecule& mol) const {
    return calcAtomicFraction(mol, [](const RDKit::Atom* atom) {
        return atom->getTotalValence() == 1;
    });
}

std::variant<double, int, std::string> FcValence2::calculate(const Molecule& mol) const {
    return calcAtomicFraction(mol, [](const RDKit::Atom* atom) {
        return atom->getTotalValence() == 2;
    });
}

std::variant<double, int, std::string> FcValence3::calculate(const Molecule& mol) const {
    return calcAtomicFraction(mol, [](const RDKit::Atom* atom) {
        return atom->getTotalValence() == 3;
    });
}

std::variant<double, int, std::string> FcValence4::calculate(const Molecule& mol) const {
    return calcAtomicFraction(mol, [](const RDKit::Atom* atom) {
        return atom->getTotalValence() == 4;
    });
}

std::variant<double, int, std::string> FcValence5Plus::calculate(const Molecule& mol) const {
    return calcAtomicFraction(mol, [](const RDKit::Atom* atom) {
        return atom->getTotalValence() >= 5;
    });
}

std::variant<double, int, std::string> FcDegree1::calculate(const Molecule& mol) const {
    return calcAtomicFraction(mol, [](const RDKit::Atom* atom) {
        return atom->getDegree() == 1;
    });
}

std::variant<double, int, std::string> FcDegree2::calculate(const Molecule& mol) const {
    return calcAtomicFraction(mol, [](const RDKit::Atom* atom) {
        return atom->getDegree() == 2;
    });
}

std::variant<double, int, std::string> FcDegree3::calculate(const Molecule& mol) const {
    return calcAtomicFraction(mol, [](const RDKit::Atom* atom) {
        return atom->getDegree() == 3;
    });
}

std::variant<double, int, std::string> FcDegree4::calculate(const Molecule& mol) const {
    return calcAtomicFraction(mol, [](const RDKit::Atom* atom) {
        return atom->getDegree() == 4;
    });
}

std::variant<double, int, std::string> FcSingleBonds::calculate(const Molecule& mol) const {
    return calcBondFraction(mol, [](const RDKit::Bond* bond) {
        return bond->getBondType() == RDKit::Bond::SINGLE;
    });
}

std::variant<double, int, std::string> FcDoubleBonds::calculate(const Molecule& mol) const {
    return calcBondFraction(mol, [](const RDKit::Bond* bond) {
        return bond->getBondType() == RDKit::Bond::DOUBLE;
    });
}

std::variant<double, int, std::string> FcTripleBonds::calculate(const Molecule& mol) const {
    return calcBondFraction(mol, [](const RDKit::Bond* bond) {
        return bond->getBondType() == RDKit::Bond::TRIPLE;
    });
}

std::variant<double, int, std::string> FcAromaticBonds::calculate(const Molecule& mol) const {
    return calcBondFraction(mol, [](const RDKit::Bond* bond) {
        return bond->getIsAromatic();
    });
}

std::variant<double, int, std::string> FcBondPolarityHigh::calculate(const Molecule& mol) const {
    return calcBondFraction(mol, [](const RDKit::Bond* bond) {
        double en1 = getAtomElectronegativity(bond->getBeginAtom());
        double en2 = getAtomElectronegativity(bond->getEndAtom());
        return std::abs(en1 - en2) > 1.0;
    });
}

std::variant<double, int, std::string> FcBondPolarityMedium::calculate(const Molecule& mol) const {
    return calcBondFraction(mol, [](const RDKit::Bond* bond) {
        double en1 = getAtomElectronegativity(bond->getBeginAtom());
        double en2 = getAtomElectronegativity(bond->getEndAtom());
        double diff = std::abs(en1 - en2);
        return diff >= 0.5 && diff <= 1.0;
    });
}

std::variant<double, int, std::string> FcBondPolarityLow::calculate(const Molecule& mol) const {
    return calcBondFraction(mol, [](const RDKit::Bond* bond) {
        double en1 = getAtomElectronegativity(bond->getBeginAtom());
        double en2 = getAtomElectronegativity(bond->getEndAtom());
        return std::abs(en1 - en2) < 0.5;
    });
}

std::variant<double, int, std::string> FcBondHighLowPolz::calculate(const Molecule& mol) const {
    return calcBondFraction(mol, [](const RDKit::Bond* bond) {
        double polz1 = getAtomPolarizability(bond->getBeginAtom());
        double polz2 = getAtomPolarizability(bond->getEndAtom());
        bool high1 = polz1 > 3.5;
        bool low1 = polz1 < 1.0;
        bool high2 = polz2 > 3.5;
        bool low2 = polz2 < 1.0;
        return (high1 && low2) || (low1 && high2);
    });
}

std::variant<double, int, std::string> FcBondHighLowIE::calculate(const Molecule& mol) const {
    return calcBondFraction(mol, [](const RDKit::Bond* bond) {
        double ie1 = getAtomIonizationEnergy(bond->getBeginAtom());
        double ie2 = getAtomIonizationEnergy(bond->getEndAtom());
        bool high1 = ie1 > 13.0;
        bool low1 = ie1 < 9.0;
        bool high2 = ie2 > 13.0;
        bool low2 = ie2 < 9.0;
        return (high1 && low2) || (low1 && high2);
    });
}

std::variant<double, int, std::string> FcPotentialChiralAtoms::calculate(const Molecule& mol) const {
     if (!mol.isValid() || !mol.getMolecule()) return 0.0;
     // Create a mutable copy
     RDKit::ROMol molCopy = *mol.getMolecule(); 
     
     if (molCopy.getNumAtoms() == 0) return 0.0;

     // Force assignment on the copy
     RDKit::MolOps::assignStereochemistry(molCopy, true, true, true); 
     int potentialChiralCount = 0;
     // Iterate over the copy's atoms
     for (const RDKit::Atom* atom : molCopy.atoms()) { 
         if (atom->hasProp(RDKit::common_properties::_ChiralityPossible)) {
             potentialChiralCount++;
         }
     }
     return static_cast<double>(potentialChiralCount) / molCopy.getNumAtoms();
}


std::variant<double, int, std::string> FcSpecifiedChiralAtoms::calculate(const Molecule& mol) const {
    return calcAtomicFraction(mol, [](const RDKit::Atom* atom) {
        RDKit::Atom::ChiralType tag = atom->getChiralTag();
        return tag != RDKit::Atom::CHI_UNSPECIFIED && tag != RDKit::Atom::CHI_OTHER;
    });
}

std::variant<double, int, std::string> FcSpiroAtoms::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (rdkMol->getNumAtoms() == 0) return 0.0;
    ensureRingInfo(*rdkMol);
    const auto* ringInfo = rdkMol->getRingInfo();
    int spiroCount = 0;

    // Reuse logic from SpiroAtomCount descriptor in counts.cpp
     for(const auto& atom : rdkMol->atoms()){
         if(ringInfo->isAtomInRingOfSize(atom->getIdx(), 3) ||
            ringInfo->isAtomInRingOfSize(atom->getIdx(), 4) ||
            ringInfo->isAtomInRingOfSize(atom->getIdx(), 5) ||
            ringInfo->isAtomInRingOfSize(atom->getIdx(), 6) ||
            ringInfo->isAtomInRingOfSize(atom->getIdx(), 7) ||
            ringInfo->isAtomInRingOfSize(atom->getIdx(), 8))
         {
              std::vector<int> atomRings;
              for(size_t ring_idx = 0; ring_idx < ringInfo->atomRings().size(); ++ring_idx) {
                  const auto& ring = ringInfo->atomRings()[ring_idx];
                  if(std::find(ring.begin(), ring.end(), atom->getIdx()) != ring.end()) {
                      atomRings.push_back(ring_idx);
                  }
              }

              if (atomRings.size() >= 2) {
                  bool isSpiro = true;
                  for(size_t i = 0; i < atomRings.size(); ++i) {
                      for(size_t j = i + 1; j < atomRings.size(); ++j) {
                          const auto& ring1 = ringInfo->atomRings()[atomRings[i]];
                          const auto& ring2 = ringInfo->atomRings()[atomRings[j]];
                          int sharedAtoms = 0;
                          for(int atom_idx1 : ring1) {
                              for(int atom_idx2 : ring2) {
                                  if (atom_idx1 == atom_idx2) sharedAtoms++;
                              }
                          }
                          if (sharedAtoms > 1) {
                              isSpiro = false;
                              break;
                          }
                      }
                      if (!isSpiro) break;
                  }
                  if(isSpiro) spiroCount++;
              }
         }
     }
    return static_cast<double>(spiroCount) / rdkMol->getNumAtoms();
}

std::variant<double, int, std::string> FcCNBondFraction::calculate(const Molecule& mol) const {
    return calcBondFraction(mol, [](const RDKit::Bond* bond) {
        int z1 = bond->getBeginAtom()->getAtomicNum();
        int z2 = bond->getEndAtom()->getAtomicNum();
        return (z1 == 6 && z2 == 7) || (z1 == 7 && z2 == 6);
    });
}

std::variant<double, int, std::string> FcCOBondFraction::calculate(const Molecule& mol) const {
    return calcBondFraction(mol, [](const RDKit::Bond* bond) {
        int z1 = bond->getBeginAtom()->getAtomicNum();
        int z2 = bond->getEndAtom()->getAtomicNum();
        return (z1 == 6 && z2 == 8) || (z1 == 8 && z2 == 6);
    });
}

std::variant<double, int, std::string> FcCSBondFraction::calculate(const Molecule& mol) const {
    return calcBondFraction(mol, [](const RDKit::Bond* bond) {
        int z1 = bond->getBeginAtom()->getAtomicNum();
        int z2 = bond->getEndAtom()->getAtomicNum();
        return (z1 == 6 && z2 == 16) || (z1 == 16 && z2 == 6);
    });
}

}  // namespace descriptors
}  // namespace desfact