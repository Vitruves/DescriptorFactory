// ... existing code ...
#include "descriptors/counts.hpp"
#include <GraphMol/GraphMol.h>
#include <GraphMol/Atom.h>
#include <GraphMol/Bond.h>
#include <GraphMol/AtomIterators.h>
#include <GraphMol/BondIterators.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/Descriptors/MolDescriptors.h>
#include <set>
#include <map>
#include <vector>
#include <queue>
#include <cmath>
#include <GraphMol/RingInfo.h>
#include <GraphMol/Chirality.h>
#include <GraphMol/Chirality.h>

namespace desfact {
namespace descriptors {

namespace {
const std::set<int> halogens = {9, 17, 35, 53, 85};
const std::set<int> acidicElements = {8, 16, 7}; // O, S, N
const std::set<int> basicElements = {7}; // N
const std::set<int> electronegative = {7, 8, 9, 16, 17, 35, 53}; // > C
double getEN(int atomicNum) {
    static const std::map<int, double> en = {
        {1, 2.20}, {6, 2.55}, {7, 3.04}, {8, 3.44}, {9, 3.98},
        {16, 2.58}, {17, 3.16}, {35, 2.96}, {53, 2.66}
    };
    auto it = en.find(atomicNum);
    return it != en.end() ? it->second : 0.0;
}
bool isENAboveC(int atomicNum) { return getEN(atomicNum) > getEN(6); }
bool isAcidicAtom(const RDKit::Atom* atom) {
    int num = atom->getAtomicNum();
    if (num == 8 || num == 16) return atom->getDegree() > 0;
    if (num == 7) return atom->getDegree() > 1;
    return false;
}
bool isBasicAtom(const RDKit::Atom* atom) {
    int num = atom->getAtomicNum();
    if (num == 7) return atom->getDegree() > 1;
    return false;
}
std::vector<const RDKit::Atom*> getAtomsByPredicate(const RDKit::ROMol& mol, std::function<bool(const RDKit::Atom*)> pred) {
    std::vector<const RDKit::Atom*> result;
    for (const auto& atom : mol.atoms()) {
        if (pred(atom)) result.push_back(atom);
    }
    return result;
}
int countENAtomsAtDistance(const RDKit::ROMol& mol, const std::vector<const RDKit::Atom*>& startAtoms, int distance) {
    std::set<int> found;
    for (const auto* atom : startAtoms) {
        std::set<int> visited;
        std::queue<std::pair<const RDKit::Atom*, int>> q;
        q.push({atom, 0});
        visited.insert(atom->getIdx());
        while (!q.empty()) {
            auto [curr, dist] = q.front();
            q.pop();
            if (dist == distance && isENAboveC(curr->getAtomicNum())) found.insert(curr->getIdx());
            if (dist >= distance) continue;
            for (const auto& nbr : mol.atomNeighbors(curr)) {
                int idx = nbr->getIdx();
                if (!visited.count(idx)) {
                    visited.insert(idx);
                    q.push({nbr, dist + 1});
                }
            }
        }
    }
    return found.size();
}
} // namespace

HalogenCount::HalogenCount() : Descriptor("HalogenCount", "Number of halogen atoms (F, Cl, Br, I, At)") {}
std::variant<double, int, std::string> HalogenCount::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0;
    int count = 0;
    for (const auto& atom : mol.getMolecule()->atoms()) if (halogens.count(atom->getAtomicNum())) count++;
    return count;
}

CarbonCount::CarbonCount() : Descriptor("CarbonCount", "Number of carbon atoms") {}
std::variant<double, int, std::string> CarbonCount::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0;
    int count = 0;
    for (const auto& atom : mol.getMolecule()->atoms()) if (atom->getAtomicNum() == 6) count++;
    return count;
}

HydrogenCount::HydrogenCount() : Descriptor("HydrogenCount", "Number of hydrogen atoms (explicit + implicit)") {}
std::variant<double, int, std::string> HydrogenCount::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0;
    int count = 0;
    for (const auto& atom : mol.getMolecule()->atoms()) count += atom->getTotalNumHs();
    return count;
}

OxygenCount::OxygenCount() : Descriptor("OxygenCount", "Number of oxygen atoms") {}
std::variant<double, int, std::string> OxygenCount::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0;
    int count = 0;
    for (const auto& atom : mol.getMolecule()->atoms()) if (atom->getAtomicNum() == 8) count++;
    return count;
}

NitrogenCount::NitrogenCount() : Descriptor("NitrogenCount", "Number of nitrogen atoms") {}
std::variant<double, int, std::string> NitrogenCount::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0;
    int count = 0;
    for (const auto& atom : mol.getMolecule()->atoms()) if (atom->getAtomicNum() == 7) count++;
    return count;
}

SulfurCount::SulfurCount() : Descriptor("SulfurCount", "Number of sulfur atoms") {}
std::variant<double, int, std::string> SulfurCount::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0;
    int count = 0;
    for (const auto& atom : mol.getMolecule()->atoms()) if (atom->getAtomicNum() == 16) count++;
    return count;
}

AromaticAtomCount::AromaticAtomCount() : Descriptor("AromaticAtomCount", "Number of aromatic atoms") {}
std::variant<double, int, std::string> AromaticAtomCount::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0;
    int count = 0;
    for (const auto& atom : mol.getMolecule()->atoms()) if (atom->getIsAromatic()) count++;
    return count;
}

NonAromaticAtomCount::NonAromaticAtomCount() : Descriptor("NonAromaticAtomCount", "Number of non-aromatic atoms") {}
std::variant<double, int, std::string> NonAromaticAtomCount::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0;
    int count = 0;
    for (const auto& atom : mol.getMolecule()->atoms()) if (!atom->getIsAromatic()) count++;
    return count;
}

DoubleBondCount::DoubleBondCount() : Descriptor("DoubleBondCount", "Number of double bonds") {}
std::variant<double, int, std::string> DoubleBondCount::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0;
    int count = 0;
    for (const auto& bond : mol.getMolecule()->bonds()) if (bond->getBondType() == RDKit::Bond::DOUBLE) count++;
    return count;
}

TripleBondCount::TripleBondCount() : Descriptor("TripleBondCount", "Number of triple bonds") {}
std::variant<double, int, std::string> TripleBondCount::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0;
    int count = 0;
    for (const auto& bond : mol.getMolecule()->bonds()) if (bond->getBondType() == RDKit::Bond::TRIPLE) count++;
    return count;
}

SingleBondCount::SingleBondCount() : Descriptor("SingleBondCount", "Number of single bonds") {}
std::variant<double, int, std::string> SingleBondCount::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0;
    int count = 0;
    for (const auto& bond : mol.getMolecule()->bonds()) if (bond->getBondType() == RDKit::Bond::SINGLE) count++;
    return count;
}

AcidicFunctionCount::AcidicFunctionCount() : Descriptor("AcidicFunctionCount", "Number of acidic functional groups (heuristic)") {}
std::variant<double, int, std::string> AcidicFunctionCount::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0;
    int count = 0;
    for (const auto& atom : mol.getMolecule()->atoms()) if (isAcidicAtom(atom)) count++;
    return count;
}

BasicFunctionCount::BasicFunctionCount() : Descriptor("BasicFunctionCount", "Number of basic functional groups (heuristic)") {}
std::variant<double, int, std::string> BasicFunctionCount::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0;
    int count = 0;
    for (const auto& atom : mol.getMolecule()->atoms()) if (isBasicAtom(atom)) count++;
    return count;
}

ENAtoms2BondsFromAcidic::ENAtoms2BondsFromAcidic() : Descriptor("ENAtoms2BondsFromAcidic", "Number of electronegative atoms (>C) two bonds from acidic function") {}
std::variant<double, int, std::string> ENAtoms2BondsFromAcidic::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0;
    auto acids = getAtomsByPredicate(*mol.getMolecule(), isAcidicAtom);
    return countENAtomsAtDistance(*mol.getMolecule(), acids, 2);
}

ENAtoms2BondsFromBasic::ENAtoms2BondsFromBasic() : Descriptor("ENAtoms2BondsFromBasic", "Number of electronegative atoms (>C) two bonds from basic function") {}
std::variant<double, int, std::string> ENAtoms2BondsFromBasic::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0;
    auto basics = getAtomsByPredicate(*mol.getMolecule(), isBasicAtom);
    return countENAtomsAtDistance(*mol.getMolecule(), basics, 2);
}

ENAtoms3BondsFromAcidic::ENAtoms3BondsFromAcidic() : Descriptor("ENAtoms3BondsFromAcidic", "Number of electronegative atoms (>C) three bonds from acidic function") {}
std::variant<double, int, std::string> ENAtoms3BondsFromAcidic::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0;
    auto acids = getAtomsByPredicate(*mol.getMolecule(), isAcidicAtom);
    return countENAtomsAtDistance(*mol.getMolecule(), acids, 3);
}

ENAtoms3BondsFromBasic::ENAtoms3BondsFromBasic() : Descriptor("ENAtoms3BondsFromBasic", "Number of electronegative atoms (>C) three bonds from basic function") {}
std::variant<double, int, std::string> ENAtoms3BondsFromBasic::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0;
    auto basics = getAtomsByPredicate(*mol.getMolecule(), isBasicAtom);
    return countENAtomsAtDistance(*mol.getMolecule(), basics, 3);
}

LongestCSequenceSmiles::LongestCSequenceSmiles() : Descriptor("LongestCSequenceSmiles", "Longest uninterrupted sequence of 'c' or 'C' in SMILES") {}
std::variant<double, int, std::string> LongestCSequenceSmiles::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid()) return 0;
    const std::string& smiles = mol.getSmiles();
    int maxLen = 0, curLen = 0;
    for (char ch : smiles) {
        if (ch == 'c' || ch == 'C') {
            curLen++;
            if (curLen > maxLen) maxLen = curLen;
        } else {
            curLen = 0;
        }
    }
    return maxLen;
}

UppercaseCountSmiles::UppercaseCountSmiles() : Descriptor("UppercaseCountSmiles", "Count of uppercase letters in SMILES") {}
std::variant<double, int, std::string> UppercaseCountSmiles::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid()) return 0;
    const std::string& smiles = mol.getSmiles();
    int count = 0;
    for (char ch : smiles) if (ch >= 'A' && ch <= 'Z') count++;
    return count;
}

LowercaseCountSmiles::LowercaseCountSmiles() : Descriptor("LowercaseCountSmiles", "Count of lowercase letters in SMILES") {}
std::variant<double, int, std::string> LowercaseCountSmiles::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid()) return 0;
    const std::string& smiles = mol.getSmiles();
    int count = 0;
    for (char ch : smiles) if (ch >= 'a' && ch <= 'z') count++;
    return count;
}

AtomVolumeSum::AtomVolumeSum() : Descriptor("AtomVolumeSum", "Sum of atomic van der Waals volumes") {}
std::variant<double, int, std::string> AtomVolumeSum::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    static const std::map<int, double> vdw = {
        {1, 7.24}, {6, 20.58}, {7, 15.60}, {8, 14.71}, {9, 13.31},
        {15, 24.43}, {16, 24.43}, {17, 22.45}, {35, 27.07}, {53, 35.19}
    };
    double sum = 0.0;
    for (const auto& atom : mol.getMolecule()->atoms()) {
        int z = atom->getAtomicNum();
        auto it = vdw.find(z);
        sum += (it != vdw.end() ? it->second : 20.0);
    }
    return sum;
}

HeavyAtomCount::HeavyAtomCount() : Descriptor("HeavyAtomCount", "Number of heavy atoms (non-hydrogen)") {}
std::variant<double, int, std::string> HeavyAtomCount::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0;
    int count = 0;
    for (const auto& atom : mol.getMolecule()->atoms()) if (atom->getAtomicNum() > 1) count++;
    return count;
}

HeteroatomCount::HeteroatomCount() : Descriptor("HeteroatomCount", "Number of heteroatoms (not C or H)") {}
std::variant<double, int, std::string> HeteroatomCount::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0;
    int count = 0;
    for (const auto& atom : mol.getMolecule()->atoms()) {
        int z = atom->getAtomicNum();
        if (z > 1 && z != 6) count++;
    }
    return count;
}

RingAtomCount::RingAtomCount() : Descriptor("RingAtomCount", "Number of atoms in rings") {}
std::variant<double, int, std::string> RingAtomCount::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0;
    const auto* ringInfo = mol.getMolecule()->getRingInfo();
    int count = 0;
    for (const auto& atom : mol.getMolecule()->atoms())
        if (ringInfo->numAtomRings(atom->getIdx()) > 0) count++;
    return count;
}

RingCount::RingCount() : Descriptor("RingCount", "Number of rings (SSSR)") {}
std::variant<double, int, std::string> RingCount::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0;
    return static_cast<int>(mol.getMolecule()->getRingInfo()->numRings());
}

ChiralCenterCount::ChiralCenterCount() : Descriptor("ChiralCenterCount", "Number of chiral centers") {}
std::variant<double, int, std::string> ChiralCenterCount::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0;
    
    // Alternative approach using Atom properties
    int count = 0;
    for (const auto& atom : mol.getMolecule()->atoms()) {
        if (atom->getChiralTag() != RDKit::Atom::CHI_UNSPECIFIED && 
            atom->getChiralTag() != RDKit::Atom::CHI_OTHER) {
            count++;
        }
    }
    return count;
}

FormalChargeCount::FormalChargeCount() : Descriptor("FormalChargeCount", "Number of atoms with nonzero formal charge") {}
std::variant<double, int, std::string> FormalChargeCount::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0;
    int count = 0;
    for (const auto& atom : mol.getMolecule()->atoms()) if (atom->getFormalCharge() != 0) count++;
    return count;
}

PositiveChargeCount::PositiveChargeCount() : Descriptor("PositiveChargeCount", "Number of atoms with positive formal charge") {}
std::variant<double, int, std::string> PositiveChargeCount::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0;
    int count = 0;
    for (const auto& atom : mol.getMolecule()->atoms()) if (atom->getFormalCharge() > 0) count++;
    return count;
}

NegativeChargeCount::NegativeChargeCount() : Descriptor("NegativeChargeCount", "Number of atoms with negative formal charge") {}
std::variant<double, int, std::string> NegativeChargeCount::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0;
    int count = 0;
    for (const auto& atom : mol.getMolecule()->atoms()) if (atom->getFormalCharge() < 0) count++;
    return count;
}

RotatableBondCount::RotatableBondCount() : Descriptor("RotatableBondCount", "Number of rotatable bonds") {}
std::variant<double, int, std::string> RotatableBondCount::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0;
    return static_cast<int>(RDKit::Descriptors::calcNumRotatableBonds(*mol.getMolecule()));
}

BridgeheadAtomCount::BridgeheadAtomCount() : Descriptor("BridgeheadAtomCount", "Number of bridgehead atoms (proxy: atoms in >1 ring)") {}
std::variant<double, int, std::string> BridgeheadAtomCount::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0;
    const auto* ringInfo = mol.getMolecule()->getRingInfo();
    int nAtoms = mol.getMolecule()->getNumAtoms();
    int count = 0;
    for (int i = 0; i < nAtoms; ++i)
        if (ringInfo->numAtomRings(i) > 1) count++;
    return count;
}

} // namespace descriptors
} // namespace desfact