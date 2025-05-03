// ... existing code ...
#include "descriptors/counts.hpp"
#include <GraphMol/GraphMol.h>
#include <GraphMol/Atom.h>
#include <GraphMol/Bond.h>
#include <GraphMol/AtomIterators.h>
#include <GraphMol/BondIterators.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/Descriptors/MolDescriptors.h>
#include <GraphMol/SmilesParse/SmartsWrite.h> // For SMARTS patterns
#include <GraphMol/SmilesParse/SmilesParse.h> // For SmartsToMol
#include <GraphMol/PeriodicTable.h>
#include <GraphMol/RingInfo.h>
#include <GraphMol/Chirality.h>
#include <set>
#include <map>
#include <vector>
#include <queue>
#include <cmath>
#include "descriptors/element_properties.hpp"

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

// Helper for SMARTS matching count
int countSubstructMatches(const RDKit::ROMol& mol, RDKit::ROMol* patternMol) {
    if (!patternMol) return 0;
    std::vector<RDKit::MatchVectType> matches;
    unsigned int count = RDKit::SubstructMatch(mol, *patternMol, matches);
    return static_cast<int>(count);
}

// Helper to initialize ring info if needed
void ensureRingInfo(const RDKit::ROMol& mol) {
    if (!mol.getRingInfo() || !mol.getRingInfo()->isInitialized()) {
        RDKit::MolOps::findSSSR(const_cast<RDKit::ROMol&>(mol));
    }
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

// --- New Count Implementations ---

PhosphorusCount::PhosphorusCount() : Descriptor("PhosphorusCount", "Number of phosphorus atoms") {}
std::variant<double, int, std::string> PhosphorusCount::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0;
    int count = 0;
    for (const auto& atom : mol.getMolecule()->atoms()) if (atom->getAtomicNum() == 15) count++;
    return count;
}

SiliconCount::SiliconCount() : Descriptor("SiliconCount", "Number of silicon atoms") {}
std::variant<double, int, std::string> SiliconCount::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0;
    int count = 0;
    for (const auto& atom : mol.getMolecule()->atoms()) if (atom->getAtomicNum() == 14) count++;
    return count;
}

BoronCount::BoronCount() : Descriptor("BoronCount", "Number of boron atoms") {}
std::variant<double, int, std::string> BoronCount::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0;
    int count = 0;
    for (const auto& atom : mol.getMolecule()->atoms()) if (atom->getAtomicNum() == 5) count++;
    return count;
}

MetalAtomCount::MetalAtomCount() : Descriptor("MetalAtomCount", "Number of metal atoms") {}
std::variant<double, int, std::string> MetalAtomCount::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0;
    int count = 0;
    for (const auto& atom : mol.getMolecule()->atoms()) {
        if (desfact::descriptors::ElementProperties::metals.count(atom->getAtomicNum())) {
            count++;
        }
    }
    return count;
}

NonMetalAtomCount::NonMetalAtomCount() : Descriptor("NonMetalAtomCount", "Number of non-metal atoms (excluding metalloids)") {}
std::variant<double, int, std::string> NonMetalAtomCount::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0;
    int count = 0;
    for (const auto& atom : mol.getMolecule()->atoms()) {
        int atomicNum = atom->getAtomicNum();
        if (!desfact::descriptors::ElementProperties::metals.count(atomicNum) && !desfact::descriptors::ElementProperties::metalloids.count(atomicNum)) {
            count++;
        }
    }
    return count;
}

NumSpAtoms::NumSpAtoms() : Descriptor("NumSpAtoms", "Number of sp hybridized atoms") {}
std::variant<double, int, std::string> NumSpAtoms::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0;
    int count = 0;
    for (const auto& atom : mol.getMolecule()->atoms()) {
        if (atom->getHybridization() == RDKit::Atom::SP) count++;
    }
    return count;
}

NumSp2Atoms::NumSp2Atoms() : Descriptor("NumSp2Atoms", "Number of sp2 hybridized atoms") {}
std::variant<double, int, std::string> NumSp2Atoms::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0;
    int count = 0;
    for (const auto& atom : mol.getMolecule()->atoms()) {
        if (atom->getHybridization() == RDKit::Atom::SP2) count++;
    }
    return count;
}

NumSp3Atoms::NumSp3Atoms() : Descriptor("NumSp3Atoms", "Number of sp3 hybridized atoms") {}
std::variant<double, int, std::string> NumSp3Atoms::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0;
    int count = 0;
    for (const auto& atom : mol.getMolecule()->atoms()) {
        if (atom->getHybridization() == RDKit::Atom::SP3) count++;
    }
    return count;
}

RingCount3::RingCount3() : Descriptor("RingCount3", "Number of 3-membered rings") {}
std::variant<double, int, std::string> RingCount3::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0;
    ensureRingInfo(*mol.getMolecule());
    const auto* ringInfo = mol.getMolecule()->getRingInfo();
    int count = 0;
    for (const auto& ring : ringInfo->atomRings()) {
        if (ring.size() == 3) count++;
    }
    return count;
}

RingCount4::RingCount4() : Descriptor("RingCount4", "Number of 4-membered rings") {}
std::variant<double, int, std::string> RingCount4::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0;
    ensureRingInfo(*mol.getMolecule());
    const auto* ringInfo = mol.getMolecule()->getRingInfo();
    int count = 0;
    for (const auto& ring : ringInfo->atomRings()) {
        if (ring.size() == 4) count++;
    }
    return count;
}

RingCount5::RingCount5() : Descriptor("RingCount5", "Number of 5-membered rings") {}
std::variant<double, int, std::string> RingCount5::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0;
    ensureRingInfo(*mol.getMolecule());
    const auto* ringInfo = mol.getMolecule()->getRingInfo();
    int count = 0;
    for (const auto& ring : ringInfo->atomRings()) {
        if (ring.size() == 5) count++;
    }
    return count;
}

RingCount6::RingCount6() : Descriptor("RingCount6", "Number of 6-membered rings") {}
std::variant<double, int, std::string> RingCount6::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0;
    ensureRingInfo(*mol.getMolecule());
    const auto* ringInfo = mol.getMolecule()->getRingInfo();
    int count = 0;
    for (const auto& ring : ringInfo->atomRings()) {
        if (ring.size() == 6) count++;
    }
    return count;
}

RingCountLarge::RingCountLarge() : Descriptor("RingCountLarge", "Number of rings with > 6 atoms") {}
std::variant<double, int, std::string> RingCountLarge::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0;
    ensureRingInfo(*mol.getMolecule());
    const auto* ringInfo = mol.getMolecule()->getRingInfo();
    int count = 0;
    for (const auto& ring : ringInfo->atomRings()) {
        if (ring.size() > 6) count++;
    }
    return count;
}

FusedBondCount::FusedBondCount() : Descriptor("FusedBondCount", "Number of bonds shared by more than one ring") {}
std::variant<double, int, std::string> FusedBondCount::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0;
    ensureRingInfo(*mol.getMolecule());
    const auto* ringInfo = mol.getMolecule()->getRingInfo();
    int count = 0;
    for (const auto& bond : mol.getMolecule()->bonds()) {
        if (ringInfo->numBondRings(bond->getIdx()) > 1) count++;
    }
    return count;
}

SpiroAtomCount::SpiroAtomCount() : Descriptor("SpiroAtomCount", "Number of spiro atoms") {}
std::variant<double, int, std::string> SpiroAtomCount::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0;
    ensureRingInfo(*mol.getMolecule());
    const auto* ringInfo = mol.getMolecule()->getRingInfo();
    int count = 0;
    for(const auto& atom : mol.getMolecule()->atoms()){
        if(ringInfo->isAtomInRingOfSize(atom->getIdx(), 3) ||
           ringInfo->isAtomInRingOfSize(atom->getIdx(), 4) ||
           ringInfo->isAtomInRingOfSize(atom->getIdx(), 5) ||
           ringInfo->isAtomInRingOfSize(atom->getIdx(), 6) ||
           ringInfo->isAtomInRingOfSize(atom->getIdx(), 7) || // Check common sizes
           ringInfo->isAtomInRingOfSize(atom->getIdx(), 8)) // Adjust max size if needed
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
                 // Check if any pair of rings the atom belongs to share more than just this atom
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
                         if (sharedAtoms > 1) { // Share more than the potential spiro atom
                             isSpiro = false;
                             break;
                         }
                     }
                     if (!isSpiro) break;
                 }
                 if(isSpiro) count++;
             }
        }
    }
    return count;
}


MacrocycleRingCount::MacrocycleRingCount() : Descriptor("MacrocycleRingCount", "Number of macrocyclic rings (> 12 atoms)") {}
std::variant<double, int, std::string> MacrocycleRingCount::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0;
    ensureRingInfo(*mol.getMolecule());
    const auto* ringInfo = mol.getMolecule()->getRingInfo();
    int count = 0;
    for (const auto& ring : ringInfo->atomRings()) {
        if (ring.size() > 12) count++;
    }
    return count;
}

CNBondCount::CNBondCount() : Descriptor("CNBondCount", "Number of Carbon-Nitrogen bonds") {}
std::variant<double, int, std::string> CNBondCount::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0;
    int count = 0;
    for (const auto& bond : mol.getMolecule()->bonds()) {
        int z1 = bond->getBeginAtom()->getAtomicNum();
        int z2 = bond->getEndAtom()->getAtomicNum();
        if ((z1 == 6 && z2 == 7) || (z1 == 7 && z2 == 6)) count++;
    }
    return count;
}

COBondCount::COBondCount() : Descriptor("COBondCount", "Number of Carbon-Oxygen bonds") {}
std::variant<double, int, std::string> COBondCount::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0;
    int count = 0;
    for (const auto& bond : mol.getMolecule()->bonds()) {
        int z1 = bond->getBeginAtom()->getAtomicNum();
        int z2 = bond->getEndAtom()->getAtomicNum();
        if ((z1 == 6 && z2 == 8) || (z1 == 8 && z2 == 6)) count++;
    }
    return count;
}

CSBondCount::CSBondCount() : Descriptor("CSBondCount", "Number of Carbon-Sulfur bonds") {}
std::variant<double, int, std::string> CSBondCount::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0;
    int count = 0;
    for (const auto& bond : mol.getMolecule()->bonds()) {
        int z1 = bond->getBeginAtom()->getAtomicNum();
        int z2 = bond->getEndAtom()->getAtomicNum();
        if ((z1 == 6 && z2 == 16) || (z1 == 16 && z2 == 6)) count++;
    }
    return count;
}

CHaloBondCount::CHaloBondCount() : Descriptor("CHaloBondCount", "Number of Carbon-Halogen bonds") {}
std::variant<double, int, std::string> CHaloBondCount::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0;
    int count = 0;
    for (const auto& bond : mol.getMolecule()->bonds()) {
        int z1 = bond->getBeginAtom()->getAtomicNum();
        int z2 = bond->getEndAtom()->getAtomicNum();
        bool z1_is_C = (z1 == 6);
        bool z2_is_C = (z2 == 6);
        bool z1_is_Halo = halogens.count(z1);
        bool z2_is_Halo = halogens.count(z2);
        if ((z1_is_C && z2_is_Halo) || (z1_is_Halo && z2_is_C)) count++;
    }
    return count;
}

// Static SMARTS patterns initialization
RDKit::ROMol* AlcoholGroupCount::alcoholSmarts = RDKit::SmartsToMol("[#6!$([a])][OX2H1]");
RDKit::ROMol* EtherGroupCount::etherSmarts = RDKit::SmartsToMol("[OD2]([#6])[#6]");
RDKit::ROMol* AmineCount::amineSmarts = RDKit::SmartsToMol("[NX3;!$(N=O);!$(N#*);!$([N]=*)]"); // N bonded only to C/H, not double/triple bonded, not nitro/amide like
RDKit::ROMol* EsterGroupCount::esterSmarts = RDKit::SmartsToMol("[#6][CX3](=[OX1])[OX2H0][#6]");
RDKit::ROMol* KetoneGroupCount::ketoneSmarts = RDKit::SmartsToMol("[#6][CX3](=[OX1])[#6]");
RDKit::ROMol* AldehydeGroupCount::aldehydeSmarts = RDKit::SmartsToMol("[CX3H1](=O)[#6,#1]"); // Ensure C is bonded to C or H
RDKit::ROMol* NitroGroupCount::nitroSmarts = RDKit::SmartsToMol("[$([NX3](=O)=O),$([NX3+](=O)[O-])]");
RDKit::ROMol* SulfonylGroupCount::sulfonylSmarts = RDKit::SmartsToMol("[$(S(=O)(=O))]");
RDKit::ROMol* CyanoGroupCount::cyanoSmarts = RDKit::SmartsToMol("[CX1]#[NX1]");
RDKit::ROMol* PhenylGroupCount::phenylSmarts = RDKit::SmartsToMol("c1ccccc1");


AlcoholGroupCount::AlcoholGroupCount() : Descriptor("AlcoholGroupCount", "Number of alcohol groups (-OH on non-aromatic C)") {}
std::variant<double, int, std::string> AlcoholGroupCount::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0;
    return countSubstructMatches(*mol.getMolecule(), alcoholSmarts);
}

EtherGroupCount::EtherGroupCount() : Descriptor("EtherGroupCount", "Number of ether groups (C-O-C)") {}
std::variant<double, int, std::string> EtherGroupCount::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0;
    return countSubstructMatches(*mol.getMolecule(), etherSmarts);
}

AmineCount::AmineCount() : Descriptor("AmineCount", "Number of amine nitrogen atoms (excluding amides, nitros, etc.)") {}
std::variant<double, int, std::string> AmineCount::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0;
    return countSubstructMatches(*mol.getMolecule(), amineSmarts);
}

EsterGroupCount::EsterGroupCount() : Descriptor("EsterGroupCount", "Number of ester groups") {}
std::variant<double, int, std::string> EsterGroupCount::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0;
    return countSubstructMatches(*mol.getMolecule(), esterSmarts);
}

KetoneGroupCount::KetoneGroupCount() : Descriptor("KetoneGroupCount", "Number of ketone groups (C(=O)C)") {}
std::variant<double, int, std::string> KetoneGroupCount::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0;
    return countSubstructMatches(*mol.getMolecule(), ketoneSmarts);
}

AldehydeGroupCount::AldehydeGroupCount() : Descriptor("AldehydeGroupCount", "Number of aldehyde groups (C(=O)H)") {}
std::variant<double, int, std::string> AldehydeGroupCount::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0;
    return countSubstructMatches(*mol.getMolecule(), aldehydeSmarts);
}

NitroGroupCount::NitroGroupCount() : Descriptor("NitroGroupCount", "Number of nitro groups") {}
std::variant<double, int, std::string> NitroGroupCount::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0;
    return countSubstructMatches(*mol.getMolecule(), nitroSmarts);
}

SulfonylGroupCount::SulfonylGroupCount() : Descriptor("SulfonylGroupCount", "Number of sulfonyl groups (S(=O)=O)") {}
std::variant<double, int, std::string> SulfonylGroupCount::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0;
    return countSubstructMatches(*mol.getMolecule(), sulfonylSmarts);
}

CyanoGroupCount::CyanoGroupCount() : Descriptor("CyanoGroupCount", "Number of cyano groups (C#N)") {}
std::variant<double, int, std::string> CyanoGroupCount::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0;
    return countSubstructMatches(*mol.getMolecule(), cyanoSmarts);
}

PhenylGroupCount::PhenylGroupCount() : Descriptor("PhenylGroupCount", "Number of phenyl rings") {}
std::variant<double, int, std::string> PhenylGroupCount::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0;
    return countSubstructMatches(*mol.getMolecule(), phenylSmarts);
}


} // namespace descriptors
} // namespace desfact