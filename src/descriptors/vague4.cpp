#include "descriptors/vague4.hpp"
#include <GraphMol/GraphMol.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Descriptors/MolDescriptors.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/AtomIterators.h>
#include <GraphMol/BondIterators.h>
#include <GraphMol/RingInfo.h>
#include <GraphMol/PeriodicTable.h>
#include <GraphMol/MolTransforms/MolTransforms.h> // May need for positions
#include <numeric>
#include <unordered_set>
#include <unordered_map>
#include <queue>
#include <cmath>
#include <limits>
#include <set>
#include <vector>
#include <algorithm>
#include "utils.hpp"
#include <map> // <--- Make sure map is included
#include <set> // <--- Make sure set is included

namespace desfact {
namespace descriptors {

// Anonymous namespace for helper functions specific to vague4 descriptors
namespace {

    // --- Basic Atom/Bond Properties (adapted/copied from vague3) ---
    inline unsigned int getAtomDegree(const RDKit::Atom* atom) {
        return atom ? atom->getDegree() : 0;
    }
    inline bool isTerminalAtom(const RDKit::Atom* atom) {
        return atom ? getAtomDegree(atom) <= 1 : false; // Allow degree 0 for isolated atoms
    }
    inline bool isAromaticAtom(const RDKit::Atom* atom) {
        return atom && atom->getIsAromatic();
    }
    inline bool isAtomInRing(const RDKit::Atom* atom) {
        return atom && atom->getOwningMol().getRingInfo()->numAtomRings(atom->getIdx()) > 0;
    }
    inline bool isCarbon(const RDKit::Atom* atom) {
        return atom && atom->getAtomicNum() == 6;
    }
    inline bool isHeteroatom(const RDKit::Atom* atom) {
        return atom && atom->getAtomicNum() != 6 && atom->getAtomicNum() != 1;
    }
    inline double getBondOrder(const RDKit::Bond* bond) {
        if (!bond) return 0.0;
        switch (bond->getBondType()) {
            case RDKit::Bond::BondType::SINGLE: return 1.0;
            case RDKit::Bond::BondType::DOUBLE: return 2.0;
            case RDKit::Bond::BondType::TRIPLE: return 3.0;
            case RDKit::Bond::BondType::AROMATIC: return 1.5;
            default: return 1.0; // Treat others as single?
        }
    }
    inline bool isSaturatedBond(const RDKit::Bond* bond) {
        return bond && bond->getBondType() == RDKit::Bond::SINGLE && !bond->getIsAromatic();
    }
    inline bool isUnsaturatedBond(const RDKit::Bond* bond) {
         return bond && (bond->getBondType() == RDKit::Bond::BondType::DOUBLE ||
                         bond->getBondType() == RDKit::Bond::BondType::TRIPLE ||
                         bond->getBondType() == RDKit::Bond::BondType::AROMATIC); // Include aromatic? Yes based on name.
    }
    inline bool isTripleBond(const RDKit::Bond* bond) {
        return bond && bond->getBondType() == RDKit::Bond::BondType::TRIPLE;
    }
     inline bool isConjugatedBond(const RDKit::Bond* bond) {
        return bond && bond->getIsConjugated();
    }


    // --- Electronegativity & Polarizability ---
    // Pauling scale (from vague3)
    double getElectronegativity(const RDKit::Atom* atom) {
        if (!atom) return 0.0;
        static const std::unordered_map<int, double> electronegativityValues = {
            {1, 2.20}, {2, 0.00}, {3, 0.98}, {4, 1.57}, {5, 2.04}, {6, 2.55},
            {7, 3.04}, {8, 3.44}, {9, 3.98}, {10, 0.00}, {11, 0.93}, {12, 1.31},
            {13, 1.61}, {14, 1.90}, {15, 2.19}, {16, 2.58}, {17, 3.16}, {18, 0.00},
            {19, 0.82}, {20, 1.00}, {31, 1.81}, {32, 2.01}, {33, 2.18}, {34, 2.55},
            {35, 2.96}, {36, 3.00}, {37, 0.82}, {38, 0.95}, {49, 1.78}, {50, 1.96},
            {51, 2.05}, {52, 2.10}, {53, 2.66}, {54, 2.60}
            // Add more if needed, default 0.0
        };
        int atomicNum = atom->getAtomicNum();
        auto it = electronegativityValues.find(atomicNum);
        return it != electronegativityValues.end() ? it->second : 0.0;
    }

    // Atomic Polarizability (in Å^3, approximate values)
    double getPolarizability(const RDKit::Atom* atom) {
        if (!atom) return 0.0;
        // Source: CRC Handbook, Wikipedia, various computational chem resources (values can vary)
        static const std::unordered_map<int, double> polarizabilityValues = {
            {1, 0.667}, // H
            {6, 1.76},  // C
            {7, 1.10},  // N
            {8, 0.802}, // O
            {9, 0.557}, // F
            {14, 5.38}, // Si
            {15, 3.63}, // P
            {16, 2.90}, // S
            {17, 2.18}, // Cl
            {35, 3.05}, // Br
            {53, 5.35}  // I
            // Add more if needed, default 0.0
        };
        int atomicNum = atom->getAtomicNum();
        auto it = polarizabilityValues.find(atomicNum);
        return it != polarizabilityValues.end() ? it->second : 1.0; // Default to 1.0? or 0.0? Let's use 1.0 as a neutral guess.
    }

    // --- Math Helpers (adapted/copied from vague3) ---
    double calculateVariance(const std::vector<double>& values) {
        size_t count = values.size();
        if (count <= 1) return 0.0; // Variance requires >1 point

        double sum = std::accumulate(values.begin(), values.end(), 0.0);
        double mean = sum / count;
        double squaredDiffSum = 0.0;
        for (double value : values) {
            double diff = value - mean;
            squaredDiffSum += diff * diff;
        }
        return squaredDiffSum / count; // Population variance
    }

    double calculateDiversityIndex(const std::vector<int>& counts) {
        double totalCount = std::accumulate(counts.begin(), counts.end(), 0.0);
        if (totalCount <= 1.0) return 0.0; // No diversity with 0 or 1 item

        double diversity = 0.0;
        for (int count : counts) {
            if (count > 0) {
                double p = static_cast<double>(count) / totalCount;
                diversity -= p * std::log(p); // Natural log
            }
        }
        // Normalize by log(number of categories) for consistency? Optional.
        // Let's return raw Shannon index for now.
        return diversity > 0.0 ? diversity : 0.0; // Ensure non-negative
    }


    // --- Graph Traversal (adapted/copied from vague3) ---
    int getShortestPath(const RDKit::ROMol* mol, unsigned int atom1Idx, unsigned int atom2Idx) {
        if (!mol || atom1Idx >= mol->getNumAtoms() || atom2Idx >= mol->getNumAtoms()) return -1;
        if (atom1Idx == atom2Idx) return 0;

        std::vector<int> distances(mol->getNumAtoms(), -1);
        std::queue<unsigned int> queue;

        distances[atom1Idx] = 0;
        queue.push(atom1Idx);

        while (!queue.empty()) {
            unsigned int current = queue.front();
            queue.pop();

            if (current == atom2Idx) return distances[current]; // Found

            for (const auto& nbr : mol->atomNeighbors(mol->getAtomWithIdx(current))) {
                unsigned int nbrIdx = nbr->getIdx();
                if (distances[nbrIdx] == -1) {
                    distances[nbrIdx] = distances[current] + 1;
                    queue.push(nbrIdx);
                }
            }
        }
        return -1; // Not connected
    }

    // --- Chain/Ring Identification ---

    // Identify atoms NOT belonging to any SSSR ring
    std::vector<int> getChainAtomIndices(const RDKit::ROMol* mol) {
        std::vector<int> chainIndices;
        if (!mol) return chainIndices;
        const RDKit::RingInfo* ringInfo = mol->getRingInfo();
        if (!ringInfo || !ringInfo->isInitialized()) {
            RDKit::MolOps::symmetrizeSSSR(*const_cast<RDKit::ROMol*>(mol));
            ringInfo = mol->getRingInfo();
        }

        for (unsigned int i = 0; i < mol->getNumAtoms(); ++i) {
            if (!ringInfo->numAtomRings(i)) {
                chainIndices.push_back(i);
            }
        }
        return chainIndices;
    }

    // Find connected components among a subset of atoms (e.g., chain atoms)
    std::vector<std::vector<int>> findConnectedComponents(const RDKit::ROMol* mol, const std::vector<int>& atomIndices) {
        std::vector<std::vector<int>> components;
        if (!mol || atomIndices.empty()) return components;

        std::unordered_set<int> indexSet(atomIndices.begin(), atomIndices.end());
        std::unordered_set<int> visited;

        for (int startIdx : atomIndices) {
            if (visited.find(startIdx) == visited.end()) {
                std::vector<int> currentComponent;
                std::queue<int> queue;

                queue.push(startIdx);
                visited.insert(startIdx);

                while (!queue.empty()) {
                    int current = queue.front();
                    queue.pop();
                    currentComponent.push_back(current);

                    const RDKit::Atom* currentAtom = mol->getAtomWithIdx(current);
                    for (const auto& nbr : mol->atomNeighbors(currentAtom)) {
                        int nbrIdx = nbr->getIdx();
                        // Check if neighbor is in the allowed set and not visited
                        if (indexSet.count(nbrIdx) && visited.find(nbrIdx) == visited.end()) {
                            visited.insert(nbrIdx);
                            queue.push(nbrIdx);
                        }
                    }
                }
                if (!currentComponent.empty()) {
                    components.push_back(currentComponent);
                }
            }
        }
        return components;
    }

    // Identify open chain fragments (connected components of non-ring atoms)
    std::vector<std::vector<int>> getChainFragments(const RDKit::ROMol* mol) {
         if (!mol) return {};
         return findConnectedComponents(mol, getChainAtomIndices(mol));
    }


    // Get bonds belonging only to chain atoms
    std::vector<const RDKit::Bond*> getChainBonds(const RDKit::ROMol* mol, const std::vector<int>& chainAtomIndices) {
        std::vector<const RDKit::Bond*> chainBonds;
         if (!mol) return chainBonds;
        std::unordered_set<int> chainSet(chainAtomIndices.begin(), chainAtomIndices.end());
        for (const RDKit::Bond* bond : mol->bonds()) {
            if (chainSet.count(bond->getBeginAtomIdx()) && chainSet.count(bond->getEndAtomIdx())) {
                chainBonds.push_back(bond);
            }
        }
        return chainBonds;
    }

    // --- Hydrogen Bonding ---
    bool isHBondDonor(const RDKit::Atom* atom) {
        if (!atom) return false;
        int atomicNum = atom->getAtomicNum();
        // Basic definition: O, N, S with explicit H
        if (atomicNum == 7 || atomicNum == 8 || atomicNum == 16) {
            return atom->getNumExplicitHs() > 0 || atom->getNumImplicitHs() > 0;
        }
        return false;
    }

    bool isHBondAcceptor(const RDKit::Atom* atom) {
        if (!atom) return false;
        int atomicNum = atom->getAtomicNum();
         // Basic definition: N, O, F, maybe S with lone pairs
         // RDKit's Lipinski definition is more precise but requires computation.
         // Simple heuristic: N, O, F not positively charged?
        if ((atomicNum == 7 || atomicNum == 8 || atomicNum == 9) && atom->getFormalCharge() <= 0) {
             // Check if it has lone pairs implicitly (based on valence)
             unsigned int valence = atom->getExplicitValence();
             // Needs a proper check based on typical valences to infer lone pairs
             // Simplified: just check type for now
             return true;
        }
        return false;
    }


    // --- Functional Groups & SMARTS ---
    // Helper to get atoms matching SMARTS patterns
    std::vector<std::vector<int>> getSmartsMatches(const RDKit::ROMol* mol, const std::string& smarts) {
        std::vector<std::vector<int>> matches_indices;
        if (!mol || smarts.empty()) return matches_indices;

        RDKit::RWMol* query = nullptr;
        try {
            query = RDKit::SmartsToMol(smarts);
        } catch (const std::exception& e) {
             globalLogger.debug("SMARTS parse error for '" + smarts + "': " + e.what());
             return matches_indices; // Return empty on parse error
        } catch (...) {
             globalLogger.debug("Unknown SMARTS parse error for '" + smarts + "'");
             return matches_indices; // Return empty on parse error
        }


        if (query) {
            std::vector<RDKit::MatchVectType> matches_vect;
             try {
                 RDKit::SubstructMatch(*mol, *query, matches_vect);
             } catch (const std::exception& e) {
                 globalLogger.debug("SubstructMatch error for '" + smarts + "': " + e.what());
                 delete query;
                 return matches_indices;
             } catch (...) {
                 globalLogger.debug("Unknown SubstructMatch error for '" + smarts + "'");
                 delete query;
                 return matches_indices;
             }


            for(const auto& match : matches_vect) {
                std::vector<int> current_match_indices;
                for(const auto& pair : match) {
                    current_match_indices.push_back(pair.second);
                }
                if (!current_match_indices.empty()) {
                    matches_indices.push_back(current_match_indices);
                }
            }
            delete query;
        }
        return matches_indices;
    }

    // Define common functional groups - more can be added
    const std::vector<std::pair<std::string, std::string>>& getCommonFunctionalGroups() {
         static const std::vector<std::pair<std::string, std::string>> groups = {
            {"Hydroxyl", "[#6][OX2H1]"}, // Attached to Carbon
            {"Phenol", "[c][OX2H1]"},    // Attached to aromatic C
            {"Ether", "[#6][OX2][#6]"},
            {"Amine (Pri)", "[NX3;H2;!$(NC=O)]"},
            {"Amine (Sec)", "[NX3;H1;!$(NC=O)]"},
            {"Amine (Ter)", "[NX3;H0;!$(NC=O)]"},
            {"CarboxylicAcid", "[CX3](=O)[OX2H1]"},
            {"Ester", "[CX3](=O)[OX2][#6]"},
            {"Amide", "[CX3](=O)[NX3]"},
            {"Ketone", "[#6][CX3](=O)[#6]"},
            {"Aldehyde", "[CX3H1](=O)"},
            {"Nitrile", "[NX1]#[CX2]"},
            {"Nitro", "[NX3+](=O)[O-]"},
            {"Thiol", "[#6][SX2H1]"},
            {"Sulfide", "[#6][SX2][#6]"}
         };
         return groups;
    }

    // Find all atoms belonging to any defined functional group
    std::unordered_map<int, std::string> getFunctionalGroupMembership(const RDKit::ROMol* mol) {
         std::unordered_map<int, std::string> membership; // atom_idx -> group_name
         if (!mol) return membership;

         for (const auto& group_pair : getCommonFunctionalGroups()) {
             const std::string& name = group_pair.first;
             const std::string& smarts = group_pair.second;
             auto matches = getSmartsMatches(mol, smarts);
             for(const auto& match : matches) {
                 for(int idx : match) {
                      // Allow atom to be part of multiple groups? For now, first match wins.
                      if (membership.find(idx) == membership.end()) {
                           membership[idx] = name;
                      }
                 }
             }
         }
         return membership;
    }


    // --- Other Specific Helpers ---
    // Check if an atom is peripheral (degree <= 1)
    bool isPeripheralAtom(const RDKit::Atom* atom, const RDKit::ROMol& mol) {
         if (!atom) return false;
         // Consider isolated atoms (degree 0) as peripheral
         return atom->getDegree() <= 1;
    }

    // Get atoms belonging to substituents attached to a ring atom
    std::vector<const RDKit::Atom*> getRingSubstituentAtoms(const RDKit::Atom* ringAtom) {
         std::vector<const RDKit::Atom*> substituentAtoms;
         if (!ringAtom || !isAtomInRing(ringAtom)) return substituentAtoms;

         const RDKit::ROMol& mol = ringAtom->getOwningMol();
         for (const auto& neighbor : mol.atomNeighbors(ringAtom)) {
              if (!isAtomInRing(neighbor)) {
                   // Explore the substituent branch starting from this neighbor
                   std::queue<const RDKit::Atom*> q;
                   std::unordered_set<int> visited;
                   q.push(neighbor);
                   visited.insert(neighbor->getIdx());
                   substituentAtoms.push_back(neighbor);

                   while (!q.empty()) {
                        const RDKit::Atom* current = q.front();
                        q.pop();

                        for (const auto& subNeighbor : mol.atomNeighbors(current)) {
                             // Add if not the original ring atom and not visited and not in any ring
                             if (subNeighbor->getIdx() != ringAtom->getIdx() &&
                                 visited.find(subNeighbor->getIdx()) == visited.end() &&
                                 !isAtomInRing(subNeighbor)) {
                                  visited.insert(subNeighbor->getIdx());
                                  q.push(subNeighbor);
                                  substituentAtoms.push_back(subNeighbor);
                             }
                        }
                   }
              }
         }
         return substituentAtoms;
    }

     // Check if a substituent atom is terminal within its substituent branch
     bool isTerminalSubstituentAtom(const RDKit::Atom* substAtom, const RDKit::Atom* ringAttachmentAtom) {
          if (!substAtom || !ringAttachmentAtom) return false;
          const RDKit::ROMol& mol = substAtom->getOwningMol();
          int nonRingNeighbors = 0;
          for (const auto& neighbor : mol.atomNeighbors(substAtom)) {
               // Count neighbors that are not the atom connecting to the ring
               // and are also not part of any ring themselves
               if (neighbor->getIdx() != ringAttachmentAtom->getIdx() && !isAtomInRing(neighbor)) {
                    nonRingNeighbors++;
               }
          }
          // If the only neighbor (excluding the ring attachment) is none, it's terminal
          // Or handle the case where the substituent is just one atom.
          if (substAtom->getDegree() == 1) return true; // It only connects to the ring atom
          return nonRingNeighbors == 0; // Only connected back towards the ring
     }

     // Find atoms at the ends of chain fragments
     std::vector<const RDKit::Atom*> getChainEndAtoms(const RDKit::ROMol* mol, const std::vector<std::vector<int>>& chainFragments) {
          std::vector<const RDKit::Atom*> endAtoms;
           if (!mol) return endAtoms;
          std::unordered_set<int> chainAtomSet;
          for (const auto& frag : chainFragments) {
               for (int idx : frag) {
                    chainAtomSet.insert(idx);
               }
          }

          for (const auto& fragment : chainFragments) {
               if (fragment.size() == 1) {
                    // Single atom chain, considered an end
                    endAtoms.push_back(mol->getAtomWithIdx(fragment[0]));
               } else {
                    for (int atomIdx : fragment) {
                         const RDKit::Atom* atom = mol->getAtomWithIdx(atomIdx);
                         int chainNeighbors = 0;
                         for (const auto& nbr : mol->atomNeighbors(atom)) {
                              if (chainAtomSet.count(nbr->getIdx())) {
                                   chainNeighbors++;
                              }
                         }
                         if (chainNeighbors <= 1) {
                              // Atom connects to 1 or 0 other atoms within its chain fragment
                              endAtoms.push_back(atom);
                         }
                    }
               }
          }
          return endAtoms;
     }


} // anonymous namespace

// --- Implementations ---

// 1. Structural and Connectivity-Based
AtomicConnectivityImbalanceDescriptor::AtomicConnectivityImbalanceDescriptor()
    : SumDescriptor("AtomicConnectivityImbalance", "Average absolute difference in degrees between bonded heavy atoms") {}

std::variant<double, int, std::string> AtomicConnectivityImbalanceDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return std::numeric_limits<double>::quiet_NaN();
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return std::numeric_limits<double>::quiet_NaN();

    double totalDifference = 0.0;
    int bondCount = 0;

    for (const RDKit::Bond* bond : rdkMol->bonds()) {
        const RDKit::Atom* atom1 = bond->getBeginAtom();
        const RDKit::Atom* atom2 = bond->getEndAtom();
        // Only consider bonds between heavy atoms? User didn't specify, let's include H for now.
        // If excluding H: if (atom1->getAtomicNum() > 1 && atom2->getAtomicNum() > 1) { ... }
        totalDifference += std::abs(static_cast<double>(getAtomDegree(atom1)) - static_cast<double>(getAtomDegree(atom2)));
        bondCount++;
    }

    return (bondCount > 0) ? (totalDifference / bondCount) : 0.0;
}


RingBridgeRatioDescriptor::RingBridgeRatioDescriptor()
    : FractionalDescriptor("RingBridgeRatio", "Fraction of bonds bridging two different SSSR ring systems") {}

std::variant<double, int, std::string> RingBridgeRatioDescriptor::calculate(const Molecule& mol) const {
     if (!mol.isValid() || !mol.getMolecule()) return std::numeric_limits<double>::quiet_NaN();
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
     if (!rdkMol) return std::numeric_limits<double>::quiet_NaN();

    const RDKit::RingInfo* ringInfo = rdkMol->getRingInfo();
    if (!ringInfo || !ringInfo->isInitialized()) {
        RDKit::MolOps::symmetrizeSSSR(*const_cast<RDKit::ROMol*>(rdkMol));
        ringInfo = rdkMol->getRingInfo();
    }
     if (!ringInfo) return std::numeric_limits<double>::quiet_NaN();


    int bridgeBonds = 0;
    int totalBonds = rdkMol->getNumBonds();
    if (totalBonds == 0) return 0.0;

    for (const RDKit::Bond* bond : rdkMol->bonds()) {
        int idx1 = bond->getBeginAtomIdx();
        int idx2 = bond->getEndAtomIdx();

        // Check if both atoms are in rings, but the bond itself is not in any ring
        if (ringInfo->numAtomRings(idx1) > 0 &&
            ringInfo->numAtomRings(idx2) > 0 &&
            ringInfo->numBondRings(bond->getIdx()) == 0) {
             bridgeBonds++;
        }
    }

    return static_cast<double>(bridgeBonds) / totalBonds;
}


SubstitutionPatternComplexityDescriptor::SubstitutionPatternComplexityDescriptor()
    : Descriptor("SubstitutionPatternComplexity", "Number of unique substitution patterns on aromatic rings (heuristic)") {}

std::variant<double, int, std::string> SubstitutionPatternComplexityDescriptor::calculate(const Molecule& mol) const {
     if (!mol.isValid() || !mol.getMolecule()) return -1; // Return int count
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
     if (!rdkMol) return -1;

    const RDKit::RingInfo* ringInfo = rdkMol->getRingInfo();
    if (!ringInfo || !ringInfo->isInitialized()) {
        RDKit::MolOps::symmetrizeSSSR(*const_cast<RDKit::ROMol*>(rdkMol));
        ringInfo = rdkMol->getRingInfo();
    }
     if (!ringInfo) return -1;

    std::set<std::string> uniquePatterns;
    const auto& atomRings = ringInfo->atomRings(); // vector<vector<int>>: ringIdx -> [atomIndices]

    for (const auto& ringAtomIndices : atomRings) {
         bool isAromaticRing = true;
         for (int idx : ringAtomIndices) {
              if (!rdkMol->getAtomWithIdx(idx)->getIsAromatic()) {
                   isAromaticRing = false;
                   break;
              }
         }
         if (!isAromaticRing) continue;

         // Generate a canonical pattern string for this aromatic ring
         std::map<int, std::vector<int>> substitutionMap; // ring_atom_idx -> sorted list of neighbor atomic numbers (outside ring)
         std::unordered_set<int> ringIndexSet(ringAtomIndices.begin(), ringAtomIndices.end());

         for (int ringAtomIdx : ringAtomIndices) {
             const RDKit::Atom* ringAtom = rdkMol->getAtomWithIdx(ringAtomIdx);
             std::vector<int> substituents;
             for (const auto& neighbor : rdkMol->atomNeighbors(ringAtom)) {
                 if (ringIndexSet.find(neighbor->getIdx()) == ringIndexSet.end()) {
                     // Neighbor is outside the current ring
                     substituents.push_back(neighbor->getAtomicNum());
                 }
             }
             if (!substituents.empty()) {
                  std::sort(substituents.begin(), substituents.end());
                  substitutionMap[ringAtomIdx] = substituents;
             }
         }

         // Create a string representation - canonicalize based on ring indices?
         // Simpler: just create a sorted string of substitutions
         std::string patternString = "R" + std::to_string(ringAtomIndices.size()) + ":";
         std::vector<std::string> substStrings;
         for(int ringIdx : ringAtomIndices){ // Iterate in SSSR order
            if(substitutionMap.count(ringIdx)){
                std::string s = std::to_string(ringIdx) + "(";
                for(size_t i=0; i < substitutionMap[ringIdx].size(); ++i){
                    s += std::to_string(substitutionMap[ringIdx][i]) + (i == substitutionMap[ringIdx].size() - 1 ? "" : ",");
                }
                s += ")";
                substStrings.push_back(s);
            } else {
                 substStrings.push_back(std::to_string(ringIdx) + "()");
            }
         }
         // Needs a truly canonical representation, maybe based on graph isomorphism
         // of the substituted ring, but that's too complex here.
         // Let's just sort the individual substituent strings and concatenate.
         std::sort(substStrings.begin(), substStrings.end());
         for(const auto& s : substStrings) patternString += s + ";";

         uniquePatterns.insert(patternString);
    }

    return static_cast<int>(uniquePatterns.size());
}


OpenChainSaturationRatioDescriptor::OpenChainSaturationRatioDescriptor()
    : FractionalDescriptor("OpenChainSaturationRatio", "Fraction of saturated bonds in open chains") {}

std::variant<double, int, std::string> OpenChainSaturationRatioDescriptor::calculate(const Molecule& mol) const {
     if (!mol.isValid() || !mol.getMolecule()) return std::numeric_limits<double>::quiet_NaN();
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
     if (!rdkMol) return std::numeric_limits<double>::quiet_NaN();

    auto chainAtomIndices = getChainAtomIndices(rdkMol);
    if (chainAtomIndices.empty()) return 1.0; // No chains -> vacuously all saturated? Or 0.0? Let's say 1.0 if no chain bonds exist.

    auto chainBonds = getChainBonds(rdkMol, chainAtomIndices);
    if (chainBonds.empty()) return 1.0; // Or 0.0? Let's stick to 1.0 (no unsaturated bonds found)

    int saturatedChainBonds = 0;
    for (const RDKit::Bond* bond : chainBonds) {
        if (isSaturatedBond(bond)) {
            saturatedChainBonds++;
        }
    }

    return static_cast<double>(saturatedChainBonds) / chainBonds.size();
}

// 2. Atomic Neighborhood Patterns
HeteroatomNeighborhoodDiversityDescriptor::HeteroatomNeighborhoodDiversityDescriptor()
    : SumDescriptor("HeteroatomNeighborhoodDiversity", "Avg variance of atomic numbers of heteroatom neighbors per atom") {}

std::variant<double, int, std::string> HeteroatomNeighborhoodDiversityDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return std::numeric_limits<double>::quiet_NaN();
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
     if (!rdkMol) return std::numeric_limits<double>::quiet_NaN();

    double totalVariance = 0.0;
    int atomCount = 0;

    for (const RDKit::Atom* atom : rdkMol->atoms()) {
        std::vector<double> heteroNeighborAtomicNums;
        for (const auto& neighbor : rdkMol->atomNeighbors(atom)) {
            if (isHeteroatom(neighbor)) {
                heteroNeighborAtomicNums.push_back(static_cast<double>(neighbor->getAtomicNum()));
            }
        }

        if (heteroNeighborAtomicNums.size() > 1) { // Variance requires at least 2 neighbors
            totalVariance += calculateVariance(heteroNeighborAtomicNums);
            atomCount++;
        }
    }

    return (atomCount > 0) ? (totalVariance / atomCount) : 0.0;
}


CarbonNeighborhoodUniformityDescriptor::CarbonNeighborhoodUniformityDescriptor()
    : FractionalDescriptor("CarbonNeighborhoodUniformity", "Fraction of carbon atoms bonded exclusively to carbons") {}

std::variant<double, int, std::string> CarbonNeighborhoodUniformityDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return std::numeric_limits<double>::quiet_NaN();
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
     if (!rdkMol) return std::numeric_limits<double>::quiet_NaN();

    int totalCarbons = 0;
    int uniformCarbons = 0;

    for (const RDKit::Atom* atom : rdkMol->atoms()) {
        if (isCarbon(atom)) {
            totalCarbons++;
            bool onlyCarbonNeighbors = true;
            if (atom->getDegree() == 0) { // Isolated carbon
                onlyCarbonNeighbors = false; // Or true? Let's say false.
            } else {
                for (const auto& neighbor : rdkMol->atomNeighbors(atom)) {
                    if (!isCarbon(neighbor)) {
                        onlyCarbonNeighbors = false;
                        break;
                    }
                }
            }
            if (onlyCarbonNeighbors) {
                uniformCarbons++;
            }
        }
    }

    return (totalCarbons > 0) ? (static_cast<double>(uniformCarbons) / totalCarbons) : 0.0;
}


PolarAtomNeighborhoodRatioDescriptor::PolarAtomNeighborhoodRatioDescriptor()
    : FractionalDescriptor("PolarAtomNeighborhoodRatio", "Fraction of polar atoms adjacent only to other polar atoms") {}

std::variant<double, int, std::string> PolarAtomNeighborhoodRatioDescriptor::calculate(const Molecule& mol) const {
     if (!mol.isValid() || !mol.getMolecule()) return std::numeric_limits<double>::quiet_NaN();
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
     if (!rdkMol) return std::numeric_limits<double>::quiet_NaN();

    // Define polar atoms (e.g., N, O, S, F, Cl, Br, I)
    std::unordered_set<int> polarElements = {7, 8, 16, 9, 17, 35, 53};
    auto isPolarAtom = [&](const RDKit::Atom* atom) {
        return atom && polarElements.count(atom->getAtomicNum());
    };

    int totalPolarAtoms = 0;
    int polarOnlyNeighbors = 0;

    for (const RDKit::Atom* atom : rdkMol->atoms()) {
        if (isPolarAtom(atom)) {
            totalPolarAtoms++;
            bool onlyPolarNeighbors = true;
            if (atom->getDegree() == 0) {
                 onlyPolarNeighbors = false; // Isolated polar atom
            } else {
                 for (const auto& neighbor : rdkMol->atomNeighbors(atom)) {
                     if (!isPolarAtom(neighbor)) {
                          onlyPolarNeighbors = false;
                          break;
                     }
                 }
            }
            if (onlyPolarNeighbors) {
                polarOnlyNeighbors++;
            }
        }
    }

    return (totalPolarAtoms > 0) ? (static_cast<double>(polarOnlyNeighbors) / totalPolarAtoms) : 0.0;
}


AverageAtomDegreeRangeDescriptor::AverageAtomDegreeRangeDescriptor()
    : Descriptor("AtomDegreeRange", "Difference between max and min atom degrees") {} // Corrected name

std::variant<double, int, std::string> AverageAtomDegreeRangeDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return -1; // Return int
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
     if (!rdkMol) return -1;
    if (rdkMol->getNumAtoms() == 0) return 0;


    int minDegree = std::numeric_limits<int>::max();
    int maxDegree = 0; // Degree is non-negative

    for (const RDKit::Atom* atom : rdkMol->atoms()) {
         int degree = getAtomDegree(atom);
         if (degree < minDegree) minDegree = degree;
         if (degree > maxDegree) maxDegree = degree;
    }

     if (minDegree == std::numeric_limits<int>::max()) return 0; // Should only happen for 0 atoms, already handled

    return maxDegree - minDegree;
}

// 3. Ring System Specific
RingJunctionComplexityDescriptor::RingJunctionComplexityDescriptor()
    : Descriptor("RingJunctionComplexity", "Number of atoms connecting >= 2 distinct SSSR ring systems") {}

std::variant<double, int, std::string> RingJunctionComplexityDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return -1; // Return int count
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
     if (!rdkMol) return -1;

    const RDKit::RingInfo* ringInfo = rdkMol->getRingInfo();
    if (!ringInfo || !ringInfo->isInitialized()) {
        RDKit::MolOps::symmetrizeSSSR(*const_cast<RDKit::ROMol*>(rdkMol));
        ringInfo = rdkMol->getRingInfo();
    }
     if (!ringInfo) return -1;

    int junctionAtoms = 0;
    for (unsigned int i = 0; i < rdkMol->getNumAtoms(); ++i) {
        if (ringInfo->numAtomRings(i) >= 2) { // Belongs to at least two rings
            junctionAtoms++;
        }
    }

    return junctionAtoms;
}


NonFusedRingDensityDescriptor::NonFusedRingDensityDescriptor()
    : FractionalDescriptor("NonFusedRingDensity", "Fraction of rings that are non-fused (isolated)") {}

std::variant<double, int, std::string> NonFusedRingDensityDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return std::numeric_limits<double>::quiet_NaN();
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
     if (!rdkMol) return std::numeric_limits<double>::quiet_NaN();

    const RDKit::RingInfo* ringInfo = rdkMol->getRingInfo();
    if (!ringInfo || !ringInfo->isInitialized()) {
        RDKit::MolOps::symmetrizeSSSR(*const_cast<RDKit::ROMol*>(rdkMol));
        ringInfo = rdkMol->getRingInfo();
    }
     if (!ringInfo) return std::numeric_limits<double>::quiet_NaN();

    int totalRings = ringInfo->numRings();
    if (totalRings == 0) return 1.0; // Or 0.0? Let's say 1.0 (no fused rings found).

    int nonFusedRings = 0;
    const auto& atomRings = ringInfo->atomRings(); // ringIdx -> [atomIndices]

    for (const auto& ringAtomIndices : atomRings) {
        bool isFused = false;
        for (int atomIdx : ringAtomIndices) {
            if (ringInfo->numAtomRings(atomIdx) > 1) { // Atom belongs to more than one ring
                isFused = true;
                break;
            }
        }
        if (!isFused) {
            nonFusedRings++;
        }
    }

    return static_cast<double>(nonFusedRings) / totalRings;
}


RingChainAttachmentDensityDescriptor::RingChainAttachmentDensityDescriptor()
    : FractionalDescriptor("RingChainAttachmentDensity", "Fraction of chains attached to rings") {}

std::variant<double, int, std::string> RingChainAttachmentDensityDescriptor::calculate(const Molecule& mol) const {
     if (!mol.isValid() || !mol.getMolecule()) return std::numeric_limits<double>::quiet_NaN();
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
     if (!rdkMol) return std::numeric_limits<double>::quiet_NaN();

    auto chainFragments = getChainFragments(rdkMol);
    if (chainFragments.empty()) return 0.0; // No chains

    const RDKit::RingInfo* ringInfo = rdkMol->getRingInfo();
     if (!ringInfo || !ringInfo->isInitialized()) {
        RDKit::MolOps::symmetrizeSSSR(*const_cast<RDKit::ROMol*>(rdkMol));
        ringInfo = rdkMol->getRingInfo();
    }
      if (!ringInfo) return std::numeric_limits<double>::quiet_NaN(); // Should not happen if mol is valid


    int attachedChains = 0;
    for (const auto& fragment : chainFragments) {
        bool isAttached = false;
        for (int chainAtomIdx : fragment) {
            const RDKit::Atom* chainAtom = rdkMol->getAtomWithIdx(chainAtomIdx);
            for (const auto& neighbor : rdkMol->atomNeighbors(chainAtom)) {
                if (ringInfo->numAtomRings(neighbor->getIdx()) > 0) { // Neighbor is in a ring
                    isAttached = true;
                    break;
                }
            }
            if (isAttached) break;
        }
        if (isAttached) {
            attachedChains++;
        }
    }

    return static_cast<double>(attachedChains) / chainFragments.size();
}


RingTerminalSubstituentRatioDescriptor::RingTerminalSubstituentRatioDescriptor()
    : FractionalDescriptor("RingTerminalSubstituentRatio", "Fraction of ring substituents that are terminal atoms") {}

std::variant<double, int, std::string> RingTerminalSubstituentRatioDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return std::numeric_limits<double>::quiet_NaN();
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
     if (!rdkMol) return std::numeric_limits<double>::quiet_NaN();

    const RDKit::RingInfo* ringInfo = rdkMol->getRingInfo();
    if (!ringInfo || !ringInfo->isInitialized()) {
        RDKit::MolOps::symmetrizeSSSR(*const_cast<RDKit::ROMol*>(rdkMol));
        ringInfo = rdkMol->getRingInfo();
    }
     if (!ringInfo) return std::numeric_limits<double>::quiet_NaN();

    int totalSubstituents = 0;
    int terminalSubstituents = 0;

    for (unsigned int i = 0; i < rdkMol->getNumAtoms(); ++i) {
         const RDKit::Atom* atom = rdkMol->getAtomWithIdx(i);
         if (ringInfo->numAtomRings(i) > 0) { // If it's a ring atom
              for (const auto& neighbor : rdkMol->atomNeighbors(atom)) {
                   if (ringInfo->numAtomRings(neighbor->getIdx()) == 0) { // If neighbor is a substituent atom
                        totalSubstituents++;
                        // Check if this neighbor is terminal in the context of the whole molecule
                        if (isTerminalAtom(neighbor)) {
                             terminalSubstituents++;
                        }
                   }
              }
         }
    }


    return (totalSubstituents > 0) ? (static_cast<double>(terminalSubstituents) / totalSubstituents) : 0.0;
}


RingSaturationBalanceDescriptor::RingSaturationBalanceDescriptor()
    : FractionalDescriptor("RingSaturationBalance", "Ratio of saturated (sp3) to unsaturated ring atoms") {}

std::variant<double, int, std::string> RingSaturationBalanceDescriptor::calculate(const Molecule& mol) const {
     if (!mol.isValid() || !mol.getMolecule()) return std::numeric_limits<double>::quiet_NaN();
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
     if (!rdkMol) return std::numeric_limits<double>::quiet_NaN();

    const RDKit::RingInfo* ringInfo = rdkMol->getRingInfo();
     if (!ringInfo || !ringInfo->isInitialized()) {
        RDKit::MolOps::symmetrizeSSSR(*const_cast<RDKit::ROMol*>(rdkMol));
        ringInfo = rdkMol->getRingInfo();
    }
      if (!ringInfo) return std::numeric_limits<double>::quiet_NaN();


    int saturatedRingAtoms = 0;
    int unsaturatedRingAtoms = 0;
    std::unordered_set<int> processedAtoms; // Avoid double counting atoms in multiple rings

    for(const auto& ring : ringInfo->atomRings()){
        for(int idx : ring){
             if (processedAtoms.insert(idx).second) { // If atom was not yet processed
                 const RDKit::Atom* atom = rdkMol->getAtomWithIdx(idx);
                 if (atom->getHybridization() == RDKit::Atom::HybridizationType::SP3) {
                     saturatedRingAtoms++;
                 } else { // SP2, SP, Aromatic, etc.
                     unsaturatedRingAtoms++;
                 }
             }
        }
    }


    if (unsaturatedRingAtoms == 0) {
        // Denominator is zero
        if (saturatedRingAtoms == 0) {
             // 0 / 0 case: No ring atoms, or only non-sp3/non-unsaturated types? Define as 0.0
             return 0.0;
        } else {
            // N / 0 case: All ring atoms are saturated. Use large number for ratio.
             return 9999.9;
        }
    }

    return static_cast<double>(saturatedRingAtoms) / unsaturatedRingAtoms;
}


// 4. Electronic Influence
PolarizabilityGradientDescriptor::PolarizabilityGradientDescriptor()
    : SumDescriptor("PolarizabilityGradient", "Average absolute polarizability difference between directly bonded atoms") {}

std::variant<double, int, std::string> PolarizabilityGradientDescriptor::calculate(const Molecule& mol) const {
     if (!mol.isValid() || !mol.getMolecule()) return std::numeric_limits<double>::quiet_NaN();
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
     if (!rdkMol) return std::numeric_limits<double>::quiet_NaN();

    double totalDifference = 0.0;
    int bondCount = 0;

    for (const RDKit::Bond* bond : rdkMol->bonds()) {
        const RDKit::Atom* atom1 = bond->getBeginAtom();
        const RDKit::Atom* atom2 = bond->getEndAtom();
        // Consider only heavy atoms? Let's include H for now.
        totalDifference += std::abs(getPolarizability(atom1) - getPolarizability(atom2));
        bondCount++;
    }

    return (bondCount > 0) ? (totalDifference / bondCount) : 0.0;
}


ElectronWithdrawingAtomDensityDescriptor::ElectronWithdrawingAtomDensityDescriptor()
    : FractionalDescriptor("ElectronWithdrawingAtomDensity", "Fraction of atoms bonded directly to strong EWGs (NO2, CN, CF3)") {}

std::variant<double, int, std::string> ElectronWithdrawingAtomDensityDescriptor::calculate(const Molecule& mol) const {
     if (!mol.isValid() || !mol.getMolecule()) return std::numeric_limits<double>::quiet_NaN();
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
     if (!rdkMol) return std::numeric_limits<double>::quiet_NaN();
    int nAtoms = rdkMol->getNumAtoms();
    if (nAtoms == 0) return 0.0;

    std::vector<std::string> ewgSmarts = {
        "[NX3v5](=O)=O",        // Nitro: More specific N
        "[CX2v4]#[NX1v3]",      // Cyano: More specific C, N
        "[CX4v4](F)(F)F"        // Trifluoromethyl
    };
    // Alternate SMARTS:
    // "[N+](=O)[O-]" // Nitro
    // "C#N"         // Cyano
    // "C(F)(F)F"    // CF3

    std::unordered_set<int> attachedAtoms; // Atoms attached TO the EWG

    for (const std::string& smarts : ewgSmarts) {
        auto matches = getSmartsMatches(rdkMol, smarts);
        for (const auto& matchIndices : matches) {
            // Find the atom in the match that connects to the rest of the molecule
            for (int matchedAtomIdx : matchIndices) {
                const RDKit::Atom* matchedAtom = rdkMol->getAtomWithIdx(matchedAtomIdx);
                for (const auto& neighbor : rdkMol->atomNeighbors(matchedAtom)) {
                     bool neighborInMatch = false;
                     for(int k : matchIndices) if(neighbor->getIdx() == k) { neighborInMatch = true; break;}

                     if (!neighborInMatch) {
                         // This neighbor is the attachment point
                         attachedAtoms.insert(neighbor->getIdx());
                         // Assuming EWG attaches via one atom; might need refinement
                     }
                }
            }
        }
    }

    return static_cast<double>(attachedAtoms.size()) / nAtoms;
}


ElectronDonatingAtomDensityDescriptor::ElectronDonatingAtomDensityDescriptor()
    : FractionalDescriptor("ElectronDonatingAtomDensity", "Fraction of atoms bonded directly to strong EDGs (alkoxy, amino)") {}

std::variant<double, int, std::string> ElectronDonatingAtomDensityDescriptor::calculate(const Molecule& mol) const {
     if (!mol.isValid() || !mol.getMolecule()) return std::numeric_limits<double>::quiet_NaN();
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
     if (!rdkMol) return std::numeric_limits<double>::quiet_NaN();
    int nAtoms = rdkMol->getNumAtoms();
    if (nAtoms == 0) return 0.0;

    std::vector<std::string> edgSmarts = {
        "[O;X2;!$(O=*)][#6]",        // Alkoxy/Aryloxy (O bonded to C, not carbonyl)
        "[N;X3;!$(N=*)][#6]",        // Amino (N bonded to C, not carbonyl/nitro etc) - needs refinement maybe exclude amides explicitly
    };
    // Alternative SMARTS:
    // "[OX2][CX4]", // Alkoxy (simple)
    // "[NX3][CX4]", // Amino (simple)

    std::unordered_set<int> attachedAtoms; // Atoms attached TO the EDG

    for (const std::string& smarts : edgSmarts) {
        auto matches = getSmartsMatches(rdkMol, smarts);
        for (const auto& matchIndices : matches) {
            // Logic similar to EWG: find attachment point outside the matched group
            for (int matchedAtomIdx : matchIndices) {
                const RDKit::Atom* matchedAtom = rdkMol->getAtomWithIdx(matchedAtomIdx);
                 // Heuristic: Assume O or N is the donating atom connecting out
                 if (matchedAtom->getAtomicNum() == 8 || matchedAtom->getAtomicNum() == 7) {
                     for (const auto& neighbor : rdkMol->atomNeighbors(matchedAtom)) {
                          bool neighborInMatch = false;
                          for(int k : matchIndices) if(neighbor->getIdx() == k) { neighborInMatch = true; break;}

                          if (!neighborInMatch) {
                              attachedAtoms.insert(neighbor->getIdx());
                          }
                     }
                 }
            }
        }
    }

    return static_cast<double>(attachedAtoms.size()) / nAtoms;
}


ElectronegativityGradientDensityDescriptor::ElectronegativityGradientDensityDescriptor()
    : FractionalDescriptor("ElectronegativityGradientDensity", "Fraction of bonds with electronegativity difference > 1.2") {}

std::variant<double, int, std::string> ElectronegativityGradientDensityDescriptor::calculate(const Molecule& mol) const {
     if (!mol.isValid() || !mol.getMolecule()) return std::numeric_limits<double>::quiet_NaN();
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
     if (!rdkMol) return std::numeric_limits<double>::quiet_NaN();

    int highGradientBonds = 0;
    int totalBonds = rdkMol->getNumBonds();
    if (totalBonds == 0) return 0.0;

    for (const RDKit::Bond* bond : rdkMol->bonds()) {
        double en1 = getElectronegativity(bond->getBeginAtom());
        double en2 = getElectronegativity(bond->getEndAtom());
        if (std::abs(en1 - en2) > 1.2) {
            highGradientBonds++;
        }
    }

    return static_cast<double>(highGradientBonds) / totalBonds;
}


PeripheralElectronRichAtomRatioDescriptor::PeripheralElectronRichAtomRatioDescriptor()
    : FractionalDescriptor("PeripheralElectronRichAtomRatio", "Fraction of peripheral heavy atoms with lone electron pairs") {}

std::variant<double, int, std::string> PeripheralElectronRichAtomRatioDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return std::numeric_limits<double>::quiet_NaN();
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
     if (!rdkMol) return std::numeric_limits<double>::quiet_NaN();


    int peripheralHeavyAtoms = 0;
    int peripheralElectronRich = 0;

    for (const RDKit::Atom* atom : rdkMol->atoms()) {
         if (atom->getAtomicNum() > 1 && isPeripheralAtom(atom, *rdkMol)) {
            peripheralHeavyAtoms++;
            // Check for lone pairs (heuristic: N, O, S, Halogens with typical valence)
            int atomicNum = atom->getAtomicNum();
            unsigned int explicitValence = atom->getExplicitValence();
             int formalCharge = atom->getFormalCharge();
             unsigned int numLonePairs = 0; // RDKit doesn't store this directly AFAIK

             // Simplified check based on common elements and valence rules
             bool hasLonePair = false;
             if ((atomicNum == 7 && explicitValence <= 3) || // N
                 (atomicNum == 8 && explicitValence <= 2) || // O
                 (atomicNum == 16 && explicitValence <= 2) || // S
                 ((atomicNum == 9 || atomicNum == 17 || atomicNum == 35 || atomicNum == 53) && explicitValence <= 1)) { // Halogens
                  // This check isn't perfect, doesn't handle charges well.
                  // A better way might use atom->getNumRadicalElectrons() if available/relevant
                  // or check RDKit::Descriptors::NumLipinskiHBA?
                  // Let's use isHBondAcceptor as a proxy for now.
                  if(isHBondAcceptor(atom)){ // Proxy for lone pair availability on N, O, F
                       hasLonePair = true;
                  } else if (atomicNum == 16 || atomicNum == 17 || atomicNum == 35 || atomicNum == 53) { // S, Cl, Br, I might have lone pairs even if not acceptors
                      // Basic valence check
                      unsigned int defaultValence = RDKit::PeriodicTable::getTable()->getDefaultValence(atomicNum);
                      if (explicitValence < defaultValence && formalCharge <=0) hasLonePair = true; // Simplified guess
                  }
             }


            if (hasLonePair) {
                peripheralElectronRich++;
            }
        }
    }

    return (peripheralHeavyAtoms > 0) ? (static_cast<double>(peripheralElectronRich) / peripheralHeavyAtoms) : 0.0;
}


// 5. Functional Group and Substitution Patterns
FunctionalGroupIsolationIndexDescriptor::FunctionalGroupIsolationIndexDescriptor()
    : FractionalDescriptor("FunctionalGroupIsolationIndex", "Fraction of functional groups separated by >= 4 non-functional atoms from others") {}

std::variant<double, int, std::string> FunctionalGroupIsolationIndexDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return std::numeric_limits<double>::quiet_NaN();
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
     if (!rdkMol) return std::numeric_limits<double>::quiet_NaN();

    // 1. Find all functional groups and the atoms belonging to them
    std::vector<std::vector<int>> functionalGroupsAtoms;
    std::unordered_set<int> allFunctionalAtoms;
    for (const auto& group_pair : getCommonFunctionalGroups()) {
        auto matches = getSmartsMatches(rdkMol, group_pair.second);
        for (const auto& match : matches) {
            functionalGroupsAtoms.push_back(match);
            for (int idx : match) {
                allFunctionalAtoms.insert(idx);
            }
        }
    }

    if (functionalGroupsAtoms.size() <= 1) {
        return 1.0; // 0 or 1 group is considered isolated
    }

    // 2. For each group, check minimum distance to any other group
    int isolatedGroups = 0;
    for (size_t i = 0; i < functionalGroupsAtoms.size(); ++i) {
        int minDistance = std::numeric_limits<int>::max();
        bool isIsolated = true;

        for (size_t j = 0; j < functionalGroupsAtoms.size(); ++j) {
            if (i == j) continue;

            // Find min path distance between atoms of group i and group j
            int currentPairMinDist = std::numeric_limits<int>::max();
            for (int atomI : functionalGroupsAtoms[i]) {
                for (int atomJ : functionalGroupsAtoms[j]) {
                    int dist = getShortestPath(rdkMol, atomI, atomJ);
                    if (dist != -1 && dist < currentPairMinDist) {
                        currentPairMinDist = dist;
                    }
                }
            }

             // Check if path length (excluding endpoints) consists of non-functional atoms
             // This requires path reconstruction, which is complex.
             // Let's simplify: just check the shortest path length between the groups.
            if (currentPairMinDist != std::numeric_limits<int>::max() && currentPairMinDist < 4) { // If distance < 4 bonds
                isIsolated = false;
                break;
            }
        }
        if (isIsolated) {
            isolatedGroups++;
        }
    }

    return static_cast<double>(isolatedGroups) / functionalGroupsAtoms.size();
}


HydroxylGroupDispersionDescriptor::HydroxylGroupDispersionDescriptor()
    : Descriptor("HydroxylGroupDispersion", "Variance in distances among oxygen atoms of hydroxyl groups") {}

std::variant<double, int, std::string> HydroxylGroupDispersionDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return std::numeric_limits<double>::quiet_NaN();
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
     if (!rdkMol) return std::numeric_limits<double>::quiet_NaN();

    std::vector<int> hydroxylOxygens;
    // Use SMARTS to find hydroxyl O atoms ([#6][OX2H1], [c][OX2H1])
    std::string ohSmarts = "[OX2H1;$(O-[#6;!$(C=O)])]"; // O-H attached to non-carbonyl C
    auto matches = getSmartsMatches(rdkMol, ohSmarts);
    for(const auto& match : matches) {
        if (!match.empty()) hydroxylOxygens.push_back(match[0]); // Assume O is first atom in match
    }


    if (hydroxylOxygens.size() <= 1) {
        return 0.0; // No variance with 0 or 1 group
    }

    std::vector<double> distances;
    for (size_t i = 0; i < hydroxylOxygens.size(); ++i) {
        for (size_t j = i + 1; j < hydroxylOxygens.size(); ++j) {
            int dist = getShortestPath(rdkMol, hydroxylOxygens[i], hydroxylOxygens[j]);
            if (dist != -1) { // Check if connected
                distances.push_back(static_cast<double>(dist));
            }
        }
    }

     // If no paths found between groups (disconnected molecule?), variance is 0?
     if(distances.empty()) return 0.0;

    return calculateVariance(distances);
}


AlkylChainDiversityDescriptor::AlkylChainDiversityDescriptor()
    : Descriptor("AlkylChainDiversity", "Number of distinct alkyl chain lengths per molecule") {}

std::variant<double, int, std::string> AlkylChainDiversityDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return -1; // Return int count
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
     if (!rdkMol) return -1;

    // 1. Find all sp3 carbons not in rings
    std::vector<int> alkylCarbons;
     const RDKit::RingInfo* ringInfo = rdkMol->getRingInfo();
     if (!ringInfo || !ringInfo->isInitialized()) RDKit::MolOps::symmetrizeSSSR(*const_cast<RDKit::ROMol*>(rdkMol));
     ringInfo = rdkMol->getRingInfo();
      if (!ringInfo) return -1;


    for (const RDKit::Atom* atom : rdkMol->atoms()) {
        if (atom->getAtomicNum() == 6 &&
            atom->getHybridization() == RDKit::Atom::HybridizationType::SP3 &&
            !ringInfo->numAtomRings(atom->getIdx())) {
            alkylCarbons.push_back(atom->getIdx());
        }
    }

    if (alkylCarbons.empty()) return 0;

    // 2. Find connected components of these carbons
    auto alkylFragments = findConnectedComponents(rdkMol, alkylCarbons);

    // 3. Determine the length of each chain (longest path within the fragment)
    std::set<int> uniqueLengths;
    for (const auto& fragment : alkylFragments) {
        if (fragment.empty()) continue;
        if (fragment.size() == 1) {
             // Need to consider attachment, is it CH3? or internal CH2?
             // Let's count length by number of C atoms for simplicity.
            uniqueLengths.insert(1);
            continue;
        }

        int maxLength = 0;
        // Find longest path within the fragment (expensive)
        // Approximation: just use the number of atoms in the fragment?
        // Let's use fragment size for simplicity.
        uniqueLengths.insert(fragment.size());


        // Proper way: Find diameter of the subgraph induced by fragment atoms
        // for (size_t i = 0; i < fragment.size(); ++i) {
        //     for (size_t j = i + 1; j < fragment.size(); ++j) {
        //         int dist = getShortestPath(rdkMol, fragment[i], fragment[j]); // This uses full graph path, WRONG
        //         // Need shortest path *within* the fragment subgraph
        //         // maxLength = std::max(maxLength, dist);
        //     }
        // }
        // uniqueLengths.insert(maxLength + 1); // Length in atoms
    }

    return static_cast<int>(uniqueLengths.size());
}


SubstituentPositionVarianceDescriptor::SubstituentPositionVarianceDescriptor()
    : Descriptor("SubstituentPositionVariance", "Variance of attachment positions of substituents on aromatic rings") {}

std::variant<double, int, std::string> SubstituentPositionVarianceDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return std::numeric_limits<double>::quiet_NaN();
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
     if (!rdkMol) return std::numeric_limits<double>::quiet_NaN();


    const RDKit::RingInfo* ringInfo = rdkMol->getRingInfo();
     if (!ringInfo || !ringInfo->isInitialized()) RDKit::MolOps::symmetrizeSSSR(*const_cast<RDKit::ROMol*>(rdkMol));
     ringInfo = rdkMol->getRingInfo();
      if (!ringInfo) return std::numeric_limits<double>::quiet_NaN();


    std::vector<double> positions;
    const auto& atomRings = ringInfo->atomRings(); // ringIdx -> [atomIndices]

    for (const auto& ringAtomIndices : atomRings) {
         bool isAromatic = true;
         for(int idx : ringAtomIndices) if(!rdkMol->getAtomWithIdx(idx)->getIsAromatic()) isAromatic = false;
         if (!isAromatic) continue;

         std::unordered_set<int> ringIndexSet(ringAtomIndices.begin(), ringAtomIndices.end());
         int posCounter = 0; // Assign positions 0, 1, 2... around the ring

         for (int ringAtomIdx : ringAtomIndices) {
             const RDKit::Atom* ringAtom = rdkMol->getAtomWithIdx(ringAtomIdx);
             for (const auto& neighbor : rdkMol->atomNeighbors(ringAtom)) {
                 if (ringIndexSet.find(neighbor->getIdx()) == ringIndexSet.end()) {
                     // Found a substituent attached at this position
                     positions.push_back(static_cast<double>(posCounter));
                 }
             }
             posCounter++;
         }
    }

    if (positions.size() <= 1) {
        return 0.0; // No variance if 0 or 1 substituents
    }

    // Calculate variance of collected positions [0, 1, 2, ..., 1, 4, ...]
    return calculateVariance(positions);
}


TerminalFunctionalGroupClusteringDescriptor::TerminalFunctionalGroupClusteringDescriptor()
    : FractionalDescriptor("TerminalFunctionalGroupClustering", "Fraction of terminal functional groups clustered within 2 bonds of another") {}

std::variant<double, int, std::string> TerminalFunctionalGroupClusteringDescriptor::calculate(const Molecule& mol) const {
     if (!mol.isValid() || !mol.getMolecule()) return std::numeric_limits<double>::quiet_NaN();
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
     if (!rdkMol) return std::numeric_limits<double>::quiet_NaN();

    // 1. Find all functional groups and their atoms
    std::vector<std::vector<int>> functionalGroupsAtoms;
    std::vector<bool> isTerminalGroup; // Mark if group contains a terminal atom

    for (const auto& group_pair : getCommonFunctionalGroups()) {
        auto matches = getSmartsMatches(rdkMol, group_pair.second);
        for (const auto& match : matches) {
             functionalGroupsAtoms.push_back(match);
             bool terminal = false;
             for(int idx : match){
                  if(isTerminalAtom(rdkMol->getAtomWithIdx(idx))){
                       terminal = true;
                       break;
                  }
             }
             isTerminalGroup.push_back(terminal);
        }
    }

    // 2. Filter for terminal groups
     std::vector<std::vector<int>> terminalGroups;
     for(size_t i=0; i<functionalGroupsAtoms.size(); ++i){
          if(isTerminalGroup[i]) terminalGroups.push_back(functionalGroupsAtoms[i]);
     }

     if (terminalGroups.size() <= 1) return 0.0; // No clusters if 0 or 1 terminal group

    // 3. Check pairwise distances for terminal groups
    std::unordered_set<int> clusteredGroupIndices;
    for (size_t i = 0; i < terminalGroups.size(); ++i) {
        for (size_t j = i + 1; j < terminalGroups.size(); ++j) {
            // Find min distance between group i and group j
            int minDist = std::numeric_limits<int>::max();
            for (int atomI : terminalGroups[i]) {
                for (int atomJ : terminalGroups[j]) {
                    int dist = getShortestPath(rdkMol, atomI, atomJ);
                    if (dist != -1 && dist < minDist) {
                        minDist = dist;
                    }
                }
            }

            if (minDist != -1 && minDist <= 2) { // If distance <= 2 bonds
                clusteredGroupIndices.insert(i);
                clusteredGroupIndices.insert(j);
            }
        }
    }

    return static_cast<double>(clusteredGroupIndices.size()) / terminalGroups.size();
}


// 6. Bond-type and Hybridization
ChainSaturationVariabilityDescriptor::ChainSaturationVariabilityDescriptor()
    : Descriptor("ChainSaturationVariability", "Variance of saturated bond fraction across different chains") {}

std::variant<double, int, std::string> ChainSaturationVariabilityDescriptor::calculate(const Molecule& mol) const {
     if (!mol.isValid() || !mol.getMolecule()) return std::numeric_limits<double>::quiet_NaN();
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
     if (!rdkMol) return std::numeric_limits<double>::quiet_NaN();

    auto chainFragments = getChainFragments(rdkMol);
    if (chainFragments.size() <= 1) return 0.0; // No variance with 0 or 1 chain

    std::vector<double> saturationFractions;
    for (const auto& fragment : chainFragments) {
        std::vector<const RDKit::Bond*> fragmentBonds = getChainBonds(rdkMol, fragment);
        if (fragmentBonds.empty()) {
            // saturationFractions.push_back(1.0); // Chain with no bonds is fully saturated? Or skip? Skip for now.
            continue;
        }

        int saturatedCount = 0;
        for (const RDKit::Bond* bond : fragmentBonds) {
            if (isSaturatedBond(bond)) {
                saturatedCount++;
            }
        }
        saturationFractions.push_back(static_cast<double>(saturatedCount) / fragmentBonds.size());
    }

     if (saturationFractions.size() <= 1) return 0.0; // Could happen if some chains had no bonds

    return calculateVariance(saturationFractions);
}


TripleBondTerminalRatioDescriptor::TripleBondTerminalRatioDescriptor()
    : FractionalDescriptor("TripleBondTerminalRatio", "Fraction of triple bonds involving a terminal atom") {}

std::variant<double, int, std::string> TripleBondTerminalRatioDescriptor::calculate(const Molecule& mol) const {
     if (!mol.isValid() || !mol.getMolecule()) return std::numeric_limits<double>::quiet_NaN();
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
     if (!rdkMol) return std::numeric_limits<double>::quiet_NaN();

    int totalTripleBonds = 0;
    int terminalTripleBonds = 0;

    for (const RDKit::Bond* bond : rdkMol->bonds()) {
        if (isTripleBond(bond)) {
            totalTripleBonds++;
            if (isTerminalAtom(bond->getBeginAtom()) || isTerminalAtom(bond->getEndAtom())) {
                terminalTripleBonds++;
            }
        }
    }

    return (totalTripleBonds > 0) ? (static_cast<double>(terminalTripleBonds) / totalTripleBonds) : 0.0;
}


AdjacentHybridizationTransitionDescriptor::AdjacentHybridizationTransitionDescriptor()
    : FractionalDescriptor("AdjacentHybridizationTransition", "Fraction of bonds connecting atoms with different hybridization") {}

std::variant<double, int, std::string> AdjacentHybridizationTransitionDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return std::numeric_limits<double>::quiet_NaN();
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
     if (!rdkMol) return std::numeric_limits<double>::quiet_NaN();

    int transitionBonds = 0;
    int totalBonds = rdkMol->getNumBonds();
    if (totalBonds == 0) return 0.0;

    for (const RDKit::Bond* bond : rdkMol->bonds()) {
        if (bond->getBeginAtom()->getHybridization() != bond->getEndAtom()->getHybridization()) {
            transitionBonds++;
        }
    }

    return static_cast<double>(transitionBonds) / totalBonds;
}


CyclicHybridizationHomogeneityDescriptor::CyclicHybridizationHomogeneityDescriptor()
    : SumDescriptor("CyclicHybridizationHomogeneity", "Average variance of hybridization states within each ring") {}

std::variant<double, int, std::string> CyclicHybridizationHomogeneityDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return std::numeric_limits<double>::quiet_NaN();
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
     if (!rdkMol) return std::numeric_limits<double>::quiet_NaN();

    const RDKit::RingInfo* ringInfo = rdkMol->getRingInfo();
     if (!ringInfo || !ringInfo->isInitialized()) RDKit::MolOps::symmetrizeSSSR(*const_cast<RDKit::ROMol*>(rdkMol));
     ringInfo = rdkMol->getRingInfo();
      if (!ringInfo) return std::numeric_limits<double>::quiet_NaN();

    int numRings = ringInfo->numRings();
    if (numRings == 0) return 0.0; // Perfect homogeneity if no rings? Or NaN? Let's use 0.0.

    double totalVariance = 0.0;
    auto hybToDouble = [](RDKit::Atom::HybridizationType h) -> double {
        switch(h) {
            case RDKit::Atom::HybridizationType::SP: return 1.0;
            case RDKit::Atom::HybridizationType::SP2: return 2.0;
            case RDKit::Atom::HybridizationType::SP3: return 3.0;
            // How to treat AROMATIC? Assign it SP2 value?
            // case RDKit::Atom::HybridizationType::AROMATIC: return 2.0; // Handled by getHybridization returning SP2
            default: return 0.0; // Other/unspecified
        }
    };


    for (const auto& ringAtomIndices : ringInfo->atomRings()) {
        if (ringAtomIndices.size() <= 1) continue; // Need >1 atom for variance

        std::vector<double> hybridizations;
        for (int idx : ringAtomIndices) {
            hybridizations.push_back(hybToDouble(rdkMol->getAtomWithIdx(idx)->getHybridization()));
        }
        totalVariance += calculateVariance(hybridizations);
    }

    return totalVariance / numRings;
}


InternalChainUnsaturationDensityDescriptor::InternalChainUnsaturationDensityDescriptor()
    : FractionalDescriptor("InternalChainUnsaturationDensity", "Fraction of unsaturated bonds strictly internal within chains") {}

std::variant<double, int, std::string> InternalChainUnsaturationDensityDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return std::numeric_limits<double>::quiet_NaN();
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
     if (!rdkMol) return std::numeric_limits<double>::quiet_NaN();

    auto chainFragments = getChainFragments(rdkMol);
    if (chainFragments.empty()) return 0.0;

    int internalChainBonds = 0;
    int unsaturatedInternalChainBonds = 0;

    for (const auto& fragment : chainFragments) {
         if (fragment.size() <= 2) continue; // Need at least 3 atoms for an internal bond

         std::unordered_set<int> fragmentSet(fragment.begin(), fragment.end());
         std::unordered_set<int> fragmentEndAtoms;

         // Find end atoms within this fragment context
         for(int atomIdx : fragment) {
              const RDKit::Atom* atom = rdkMol->getAtomWithIdx(atomIdx);
              int fragmentNeighbors = 0;
              for(const auto& nbr : rdkMol->atomNeighbors(atom)) {
                   if(fragmentSet.count(nbr->getIdx())) fragmentNeighbors++;
              }
              if(fragmentNeighbors <= 1) fragmentEndAtoms.insert(atomIdx);
         }


         // Iterate bonds within the fragment
         for (const RDKit::Bond* bond : getChainBonds(rdkMol, fragment)) {
             int idx1 = bond->getBeginAtomIdx();
             int idx2 = bond->getEndAtomIdx();

             // Check if *neither* atom is an end atom of this fragment
             if (fragmentEndAtoms.find(idx1) == fragmentEndAtoms.end() &&
                 fragmentEndAtoms.find(idx2) == fragmentEndAtoms.end())
             {
                 internalChainBonds++;
                 if (isUnsaturatedBond(bond)) {
                     unsaturatedInternalChainBonds++;
                 }
             }
         }
    }

    return (internalChainBonds > 0) ? (static_cast<double>(unsaturatedInternalChainBonds) / internalChainBonds) : 0.0;
}


// 7. Hydrogen Bonding Patterns
HydrogenBondDonorClusteringDescriptor::HydrogenBondDonorClusteringDescriptor()
    : FractionalDescriptor("HydrogenBondDonorClustering", "Fraction of H-bond donors clustered within <= 2 bonds of another donor") {}

std::variant<double, int, std::string> HydrogenBondDonorClusteringDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return std::numeric_limits<double>::quiet_NaN();
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
     if (!rdkMol) return std::numeric_limits<double>::quiet_NaN();

    std::vector<int> donorIndices;
    for (const RDKit::Atom* atom : rdkMol->atoms()) {
        if (isHBondDonor(atom)) {
            donorIndices.push_back(atom->getIdx());
        }
    }

    if (donorIndices.size() <= 1) return 0.0; // No clusters possible

    std::unordered_set<int> clusteredDonors;
    for (size_t i = 0; i < donorIndices.size(); ++i) {
        for (size_t j = i + 1; j < donorIndices.size(); ++j) {
            int dist = getShortestPath(rdkMol, donorIndices[i], donorIndices[j]);
            if (dist != -1 && dist <= 2) {
                clusteredDonors.insert(donorIndices[i]);
                clusteredDonors.insert(donorIndices[j]);
            }
        }
    }

    return static_cast<double>(clusteredDonors.size()) / donorIndices.size();
}


AcceptorDonorRatioImbalanceDescriptor::AcceptorDonorRatioImbalanceDescriptor()
    : FractionalDescriptor("AcceptorDonorRatioImbalance", "Absolute difference between H-bond donors and acceptors normalized by heavy atoms") {}

std::variant<double, int, std::string> AcceptorDonorRatioImbalanceDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return std::numeric_limits<double>::quiet_NaN();
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
     if (!rdkMol) return std::numeric_limits<double>::quiet_NaN();

    int donors = 0;
    int acceptors = 0;
    int heavyAtoms = 0;

    for (const RDKit::Atom* atom : rdkMol->atoms()) {
        if (atom->getAtomicNum() > 1) heavyAtoms++;
        if (isHBondDonor(atom)) donors++;
        if (isHBondAcceptor(atom)) acceptors++;
    }

    if (heavyAtoms == 0) return 0.0;

    return static_cast<double>(std::abs(donors - acceptors)) / heavyAtoms;
}


PeripheralDonorAcceptorBalanceDescriptor::PeripheralDonorAcceptorBalanceDescriptor()
    : FractionalDescriptor("PeripheralDonorAcceptorBalance", "Ratio of peripheral H-bond donors to peripheral acceptors") {}

std::variant<double, int, std::string> PeripheralDonorAcceptorBalanceDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return std::numeric_limits<double>::quiet_NaN();
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
     if (!rdkMol) return std::numeric_limits<double>::quiet_NaN();

    int peripheralDonors = 0;
    int peripheralAcceptors = 0;

    for (const RDKit::Atom* atom : rdkMol->atoms()) {
         if (isPeripheralAtom(atom, *rdkMol)) {
              if (isHBondDonor(atom)) peripheralDonors++;
              if (isHBondAcceptor(atom)) peripheralAcceptors++;
         }
    }

    if (peripheralAcceptors == 0) {
        return (peripheralDonors > 0) ? 9999.9 : 1.0; // Avoid infinity, use 1.0 if both 0
    }

    return static_cast<double>(peripheralDonors) / peripheralAcceptors;
}


IntraringHBondPotentialDescriptor::IntraringHBondPotentialDescriptor()
    : FractionalDescriptor("IntraringHBondPotential", "Potential H-bond pairs within rings normalized by total ring atoms") {}

std::variant<double, int, std::string> IntraringHBondPotentialDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return std::numeric_limits<double>::quiet_NaN();
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
     if (!rdkMol) return std::numeric_limits<double>::quiet_NaN();

    const RDKit::RingInfo* ringInfo = rdkMol->getRingInfo();
     if (!ringInfo || !ringInfo->isInitialized()) RDKit::MolOps::symmetrizeSSSR(*const_cast<RDKit::ROMol*>(rdkMol));
     ringInfo = rdkMol->getRingInfo();
      if (!ringInfo) return std::numeric_limits<double>::quiet_NaN();

    int potentialPairs = 0;
    std::unordered_set<int> allRingAtomsSet;

    for (const auto& ringAtomIndices : ringInfo->atomRings()) {
        std::vector<int> ringDonors;
        std::vector<int> ringAcceptors;
        for (int idx : ringAtomIndices) {
             allRingAtomsSet.insert(idx); // Collect all unique ring atoms
            const RDKit::Atom* atom = rdkMol->getAtomWithIdx(idx);
            if (isHBondDonor(atom)) ringDonors.push_back(idx);
            if (isHBondAcceptor(atom)) ringAcceptors.push_back(idx);
        }

        // Count potential pairs within this ring (simple D*A count)
        // A more sophisticated version might check geometry/distance within the ring.
        potentialPairs += ringDonors.size() * ringAcceptors.size();
    }

    int totalRingAtoms = allRingAtomsSet.size();
    if (totalRingAtoms == 0) return 0.0;

    return static_cast<double>(potentialPairs) / totalRingAtoms;
}


ChainEndHydrogenBondDensityDescriptor::ChainEndHydrogenBondDensityDescriptor()
    : FractionalDescriptor("ChainEndHydrogenBondDensity", "Fraction of chain-terminal heavy atoms capable of H-bonding") {}

std::variant<double, int, std::string> ChainEndHydrogenBondDensityDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return std::numeric_limits<double>::quiet_NaN();
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
     if (!rdkMol) return std::numeric_limits<double>::quiet_NaN();

    auto chainFragments = getChainFragments(rdkMol);
    if (chainFragments.empty()) return 0.0;

    auto chainEndAtoms = getChainEndAtoms(rdkMol, chainFragments);
    if (chainEndAtoms.empty()) return 0.0;

    int hBondCapableEnds = 0;
    int totalHeavyEnds = 0;

    for(const RDKit::Atom* atom : chainEndAtoms){
         if(atom->getAtomicNum() > 1) { // Only consider heavy atoms
              totalHeavyEnds++;
              if (isHBondDonor(atom) || isHBondAcceptor(atom)) {
                   hBondCapableEnds++;
              }
         }
    }


    return (totalHeavyEnds > 0) ? (static_cast<double>(hBondCapableEnds) / totalHeavyEnds) : 0.0;
}


// 8. Formal Charge Distribution
FormalChargeNeighborhoodVarianceDescriptor::FormalChargeNeighborhoodVarianceDescriptor()
    : SumDescriptor("FormalChargeNeighborhoodVariance", "Average variance of formal charges among neighbors for each atom") {}

std::variant<double, int, std::string> FormalChargeNeighborhoodVarianceDescriptor::calculate(const Molecule& mol) const {
     if (!mol.isValid() || !mol.getMolecule()) return std::numeric_limits<double>::quiet_NaN();
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
     if (!rdkMol) return std::numeric_limits<double>::quiet_NaN();

    double totalVariance = 0.0;
    int atomCount = 0; // Count atoms with >= 2 neighbors

    for (const RDKit::Atom* atom : rdkMol->atoms()) {
        std::vector<double> neighborCharges;
         unsigned int degree = atom->getDegree();
         if (degree < 2) continue; // Variance requires at least 2 neighbors

        for (const auto& neighbor : rdkMol->atomNeighbors(atom)) {
            neighborCharges.push_back(static_cast<double>(neighbor->getFormalCharge()));
        }

        totalVariance += calculateVariance(neighborCharges);
        atomCount++;
    }

    return (atomCount > 0) ? (totalVariance / atomCount) : 0.0;
}


OppositeChargeNeighborRatioDescriptor::OppositeChargeNeighborRatioDescriptor()
    : FractionalDescriptor("OppositeChargeNeighborRatio", "Fraction of bonded charged atom pairs with opposite charges") {}

std::variant<double, int, std::string> OppositeChargeNeighborRatioDescriptor::calculate(const Molecule& mol) const {
     if (!mol.isValid() || !mol.getMolecule()) return std::numeric_limits<double>::quiet_NaN();
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
     if (!rdkMol) return std::numeric_limits<double>::quiet_NaN();

    int oppositeChargeBonds = 0;
    int totalChargedBonds = 0; // Bonds between two charged atoms

    for (const RDKit::Bond* bond : rdkMol->bonds()) {
        int charge1 = bond->getBeginAtom()->getFormalCharge();
        int charge2 = bond->getEndAtom()->getFormalCharge();

        if (charge1 != 0 && charge2 != 0) {
            totalChargedBonds++;
            if (charge1 * charge2 < 0) { // Check if signs are opposite
                oppositeChargeBonds++;
            }
        }
    }

    return (totalChargedBonds > 0) ? (static_cast<double>(oppositeChargeBonds) / totalChargedBonds) : 0.0;
}


ChargeGradientAlongChainsDescriptor::ChargeGradientAlongChainsDescriptor()
    : Descriptor("ChargeGradientAlongChains", "Count of formal charge sign changes along linear chain sequences") {}

std::variant<double, int, std::string> ChargeGradientAlongChainsDescriptor::calculate(const Molecule& mol) const {
     if (!mol.isValid() || !mol.getMolecule()) return -1; // Return int count
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
     if (!rdkMol) return -1;

    auto chainFragments = getChainFragments(rdkMol);
    if (chainFragments.empty()) return 0;

    int totalSignChanges = 0;

    for (const auto& fragment : chainFragments) {
         if (fragment.size() < 2) continue; // Need at least 2 atoms for a change

         // How to define "linear sequences"? Traverse all simple paths? Expensive.
         // Let's approximate: iterate through bonds within the chain.
         std::unordered_set<int> fragmentSet(fragment.begin(), fragment.end());
         for(const RDKit::Bond* bond : getChainBonds(rdkMol, fragment)) {
              int charge1 = bond->getBeginAtom()->getFormalCharge();
              int charge2 = bond->getEndAtom()->getFormalCharge();
              // Count change from +/- to 0, 0 to +/-, + to - or - to +
              if ((charge1 > 0 && charge2 <= 0) || (charge1 < 0 && charge2 >= 0) ||
                  (charge1 == 0 && charge2 != 0) )
              {
                   totalSignChanges++;
              }
         }
    }

    return totalSignChanges;
}

PeripheralChargeNeutralityDescriptor::PeripheralChargeNeutralityDescriptor()
    : FractionalDescriptor("PeripheralChargeNeutrality", "Fraction of peripheral heavy atoms with zero charge") {}

std::variant<double, int, std::string> PeripheralChargeNeutralityDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return std::numeric_limits<double>::quiet_NaN();
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
     if (!rdkMol) return std::numeric_limits<double>::quiet_NaN();

    int peripheralHeavyAtoms = 0;
    int neutralPeripheralHeavy = 0;

    for (const RDKit::Atom* atom : rdkMol->atoms()) {
         if (atom->getAtomicNum() > 1 && isPeripheralAtom(atom, *rdkMol)) {
            peripheralHeavyAtoms++;
            if (atom->getFormalCharge() == 0) {
                neutralPeripheralHeavy++;
            }
        }
    }

    return (peripheralHeavyAtoms > 0) ? (static_cast<double>(neutralPeripheralHeavy) / peripheralHeavyAtoms) : 0.0;
}


// 9. Aromaticity & Conjugation Patterns
InterRingConjugationRatioDescriptor::InterRingConjugationRatioDescriptor()
    : FractionalDescriptor("InterRingConjugationRatio", "Fraction of conjugated bonds connecting different aromatic rings") {}

std::variant<double, int, std::string> InterRingConjugationRatioDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return std::numeric_limits<double>::quiet_NaN();
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
     if (!rdkMol) return std::numeric_limits<double>::quiet_NaN();

    const RDKit::RingInfo* ringInfo = rdkMol->getRingInfo();
     if (!ringInfo || !ringInfo->isInitialized()) RDKit::MolOps::symmetrizeSSSR(*const_cast<RDKit::ROMol*>(rdkMol));
     ringInfo = rdkMol->getRingInfo();
      if (!ringInfo) return std::numeric_limits<double>::quiet_NaN();

    int conjugatedInterRingBonds = 0;
    int totalAromaticRingBonds = 0; // Bonds where both atoms are aromatic and in a ring

    for (const RDKit::Bond* bond : rdkMol->bonds()) {
         const RDKit::Atom* atom1 = bond->getBeginAtom();
         const RDKit::Atom* atom2 = bond->getEndAtom();
         int idx1 = atom1->getIdx();
         int idx2 = atom2->getIdx();

         bool atom1AroInRing = atom1->getIsAromatic() && ringInfo->numAtomRings(idx1) > 0;
         bool atom2AroInRing = atom2->getIsAromatic() && ringInfo->numAtomRings(idx2) > 0;

         if(atom1AroInRing && atom2AroInRing) {
              totalAromaticRingBonds++;

              // Check if bond connects different rings (is not in a ring itself)
              // AND is conjugated
              if (ringInfo->numBondRings(bond->getIdx()) == 0 && isConjugatedBond(bond)) {
                   conjugatedInterRingBonds++;
              }
         }
    }

    return (totalAromaticRingBonds > 0) ? (static_cast<double>(conjugatedInterRingBonds) / totalAromaticRingBonds) : 0.0;
}


AromaticChainDensityDescriptor::AromaticChainDensityDescriptor()
    : FractionalDescriptor("AromaticChainDensity", "Fraction of aromatic atoms having at least one non-ring neighbor") {}

std::variant<double, int, std::string> AromaticChainDensityDescriptor::calculate(const Molecule& mol) const {
     if (!mol.isValid() || !mol.getMolecule()) return std::numeric_limits<double>::quiet_NaN();
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
     if (!rdkMol) return std::numeric_limits<double>::quiet_NaN();

     const RDKit::RingInfo* ringInfo = rdkMol->getRingInfo();
     if (!ringInfo || !ringInfo->isInitialized()) RDKit::MolOps::symmetrizeSSSR(*const_cast<RDKit::ROMol*>(rdkMol));
     ringInfo = rdkMol->getRingInfo();
      if (!ringInfo) return std::numeric_limits<double>::quiet_NaN();


    int totalAromaticAtoms = 0;
    int aromaticAtomsWithChainNeighbor = 0;

    for (const RDKit::Atom* atom : rdkMol->atoms()) {
        if (isAromaticAtom(atom)) {
            totalAromaticAtoms++;
            bool hasChainNeighbor = false;
            for (const auto& neighbor : rdkMol->atomNeighbors(atom)) {
                if (ringInfo->numAtomRings(neighbor->getIdx()) == 0) { // Neighbor is not in any ring
                    hasChainNeighbor = true;
                    break;
                }
            }
            if (hasChainNeighbor) {
                aromaticAtomsWithChainNeighbor++;
            }
        }
    }

    return (totalAromaticAtoms > 0) ? (static_cast<double>(aromaticAtomsWithChainNeighbor) / totalAromaticAtoms) : 0.0;
}


ConjugationLengthVarianceDescriptor::ConjugationLengthVarianceDescriptor()
    : Descriptor("ConjugationLengthVariance", "Variance in lengths (number of bonds) of conjugated segments") {}

std::variant<double, int, std::string> ConjugationLengthVarianceDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return std::numeric_limits<double>::quiet_NaN();
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
     if (!rdkMol) return std::numeric_limits<double>::quiet_NaN();

    std::vector<double> segmentLengths; // Store length of each conjugated segment in bonds
    std::vector<bool> visitedBonds(rdkMol->getNumBonds(), false);

    for (unsigned int i = 0; i < rdkMol->getNumBonds(); ++i) {
        const RDKit::Bond* startBond = rdkMol->getBondWithIdx(i);
        if (!visitedBonds[i] && isConjugatedBond(startBond)) {
            int currentLength = 0;
            std::queue<unsigned int> q;
            q.push(i);
            visitedBonds[i] = true;

            while (!q.empty()) {
                unsigned int currentBondIdx = q.front();
                q.pop();
                currentLength++;
                const RDKit::Bond* currentBond = rdkMol->getBondWithIdx(currentBondIdx);
                const RDKit::Atom* atom1 = currentBond->getBeginAtom();
                const RDKit::Atom* atom2 = currentBond->getEndAtom();

                // Explore neighbors through atoms
                 for(const RDKit::Atom* atom : {atom1, atom2}){
                      for(const RDKit::Bond* nextBond : rdkMol->atomBonds(atom)){
                           unsigned int nextBondIdx = nextBond->getIdx();
                            if (nextBondIdx != currentBondIdx && !visitedBonds[nextBondIdx] && isConjugatedBond(nextBond)) {
                                 visitedBonds[nextBondIdx] = true;
                                 q.push(nextBondIdx);
                            }
                      }
                 }

            }
            segmentLengths.push_back(static_cast<double>(currentLength));
        }
    }

    if (segmentLengths.size() <= 1) return 0.0; // No variance with 0 or 1 segment

    return calculateVariance(segmentLengths);
}


TerminalAromaticSubstitutionDescriptor::TerminalAromaticSubstitutionDescriptor()
    : FractionalDescriptor("TerminalAromaticSubstitution", "Fraction of substituents on aromatic rings that are terminal") {}

std::variant<double, int, std::string> TerminalAromaticSubstitutionDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return std::numeric_limits<double>::quiet_NaN();
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
     if (!rdkMol) return std::numeric_limits<double>::quiet_NaN();

    const RDKit::RingInfo* ringInfo = rdkMol->getRingInfo();
     if (!ringInfo || !ringInfo->isInitialized()) RDKit::MolOps::symmetrizeSSSR(*const_cast<RDKit::ROMol*>(rdkMol));
     ringInfo = rdkMol->getRingInfo();
      if (!ringInfo) return std::numeric_limits<double>::quiet_NaN();

    int totalSubstituents = 0;
    int terminalSubstituents = 0;

    for (unsigned int i = 0; i < rdkMol->getNumAtoms(); ++i) {
         const RDKit::Atom* atom = rdkMol->getAtomWithIdx(i);
         if (ringInfo->numAtomRings(i) > 0) { // If it's a ring atom
              for (const auto& neighbor : rdkMol->atomNeighbors(atom)) {
                   if (ringInfo->numAtomRings(neighbor->getIdx()) == 0) { // If neighbor is a substituent atom
                        totalSubstituents++;
                        // Check if this neighbor is terminal in the context of the whole molecule
                        if (isTerminalAtom(neighbor)) {
                             terminalSubstituents++;
                        }
                   }
              }
         }
    }


    return (totalSubstituents > 0) ? (static_cast<double>(terminalSubstituents) / totalSubstituents) : 0.0;
}


CyclicVsChainAromaticRatioDescriptor::CyclicVsChainAromaticRatioDescriptor()
    : FractionalDescriptor("CyclicVsChainAromaticRatio", "Ratio of aromatic atoms in SSSR rings vs aromatic atoms not in SSSR rings") {}

std::variant<double, int, std::string> CyclicVsChainAromaticRatioDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return std::numeric_limits<double>::quiet_NaN();
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
     if (!rdkMol) return std::numeric_limits<double>::quiet_NaN();

    const RDKit::RingInfo* ringInfo = rdkMol->getRingInfo();
     if (!ringInfo || !ringInfo->isInitialized()) RDKit::MolOps::symmetrizeSSSR(*const_cast<RDKit::ROMol*>(rdkMol));
     ringInfo = rdkMol->getRingInfo();
      if (!ringInfo) return std::numeric_limits<double>::quiet_NaN();

    int aromaticInRing = 0;
    int aromaticNotInRing = 0;

    for (const RDKit::Atom* atom : rdkMol->atoms()) {
        if (isAromaticAtom(atom)) {
            if (ringInfo->numAtomRings(atom->getIdx()) > 0) {
                aromaticInRing++;
            } else {
                aromaticNotInRing++; // Should be rare with standard aromaticity models
            }
        }
    }

    if (aromaticNotInRing == 0) {
        // Denominator is zero
        if (aromaticInRing == 0) {
            // 0 / 0 case: No aromatic atoms, ratio is 0.0.
            return 0.0;
        } else {
            // N / 0 case (where N > 0): All aromatic atoms are cyclic. Use large number for ratio.
            return 9999.9;
        }
    }

    // Denominator is non-zero, calculate the ratio normally
    return static_cast<double>(aromaticInRing) / aromaticNotInRing;
}


// 10. Structural Diversity and Complexity
UniqueElementPairRatioDescriptor::UniqueElementPairRatioDescriptor()
    : FractionalDescriptor("UniqueElementPairRatio", "Fraction of unique bonded element pairs relative to total bonds") {}

std::variant<double, int, std::string> UniqueElementPairRatioDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return std::numeric_limits<double>::quiet_NaN();
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
     if (!rdkMol) return std::numeric_limits<double>::quiet_NaN();


    int totalBonds = rdkMol->getNumBonds();
    if (totalBonds == 0) return 0.0; // Or 1.0? No pairs, no unique pairs. Let's use 0.0.

    std::set<std::pair<int, int>> uniquePairs;
    for (const RDKit::Bond* bond : rdkMol->bonds()) {
        int el1 = bond->getBeginAtom()->getAtomicNum();
        int el2 = bond->getEndAtom()->getAtomicNum();
        // Ensure canonical order (smaller element first)
        if (el1 > el2) std::swap(el1, el2);
        uniquePairs.insert({el1, el2});
    }

    return static_cast<double>(uniquePairs.size()) / totalBonds;
}


HeavyAtomDegreeDiversityDescriptor::HeavyAtomDegreeDiversityDescriptor()
    : Descriptor("HeavyAtomDegreeDiversity", "Shannon diversity index of heavy atom degrees") {}

std::variant<double, int, std::string> HeavyAtomDegreeDiversityDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return std::numeric_limits<double>::quiet_NaN();
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
     if (!rdkMol) return std::numeric_limits<double>::quiet_NaN();

    std::map<int, int> degreeCounts; // degree -> count
    int heavyAtomCount = 0;
    for (const RDKit::Atom* atom : rdkMol->atoms()) {
        if (atom->getAtomicNum() > 1) {
             heavyAtomCount++;
            degreeCounts[getAtomDegree(atom)]++;
        }
    }

    if (heavyAtomCount == 0) return 0.0;

    std::vector<int> counts;
    for (const auto& pair : degreeCounts) {
        counts.push_back(pair.second);
    }

    return calculateDiversityIndex(counts);
}


AtomPathDiversityDescriptor::AtomPathDiversityDescriptor()
    : Descriptor("AtomPathDiversity", "Number of unique shortest path lengths between distinct heavy element types") {}

std::variant<double, int, std::string> AtomPathDiversityDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return -1; // Return int count
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
     if (!rdkMol) return -1;
    int nAtoms = rdkMol->getNumAtoms();
    if (nAtoms <= 1) return 0;

    std::vector<const RDKit::Atom*> heavyAtoms;
     for(const RDKit::Atom* atom : rdkMol->atoms()){
          if(atom->getAtomicNum() > 1) heavyAtoms.push_back(atom);
     }
     if (heavyAtoms.size() <= 1) return 0;


    std::set<int> uniquePathLengths;
    for (size_t i = 0; i < heavyAtoms.size(); ++i) {
        for (size_t j = i + 1; j < heavyAtoms.size(); ++j) {
             const RDKit::Atom* atom1 = heavyAtoms[i];
             const RDKit::Atom* atom2 = heavyAtoms[j];

             // Calculate path only if element types are distinct
             if (atom1->getAtomicNum() != atom2->getAtomicNum()) {
                 int dist = getShortestPath(rdkMol, atom1->getIdx(), atom2->getIdx());
                 if (dist > 0) { // Exclude disconnected pairs and same atom
                     uniquePathLengths.insert(dist);
                 }
             }
        }
    }

    return static_cast<int>(uniquePathLengths.size());
}


InternalAtomComplexityDescriptor::InternalAtomComplexityDescriptor()
    : Descriptor("InternalAtomComplexity", "Number of internal heavy atoms bonded to atoms of >= 3 different elements") {}

std::variant<double, int, std::string> InternalAtomComplexityDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return -1; // Return int count
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
     if (!rdkMol) return -1;

    int complexInternalAtoms = 0;
    for (const RDKit::Atom* atom : rdkMol->atoms()) {
         if (atom->getAtomicNum() > 1 && atom->getDegree() > 1) { // Internal heavy atom
             std::unordered_set<int> neighborElements;
             for (const auto& neighbor : rdkMol->atomNeighbors(atom)) {
                 neighborElements.insert(neighbor->getAtomicNum());
             }
             if (neighborElements.size() >= 3) {
                 complexInternalAtoms++;
             }
         }
    }
    return complexInternalAtoms;
}


HeteroatomBondOrderVariabilityDescriptor::HeteroatomBondOrderVariabilityDescriptor()
    : SumDescriptor("HeteroatomBondOrderVariability", "Average variance in bond orders around each heteroatom") {}

std::variant<double, int, std::string> HeteroatomBondOrderVariabilityDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return std::numeric_limits<double>::quiet_NaN();
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
     if (!rdkMol) return std::numeric_limits<double>::quiet_NaN();

    double totalVariance = 0.0;
    int heteroatomCount = 0;

    for (const RDKit::Atom* atom : rdkMol->atoms()) {
        if (isHeteroatom(atom)) {
             heteroatomCount++;
             std::vector<double> bondOrders;
             // Check degree before calculating variance
             unsigned int degree = atom->getDegree();
             if (degree <= 1) continue; // Variance requires >1 bond

             for (const RDKit::Bond* bond : rdkMol->atomBonds(atom)) {
                 bondOrders.push_back(getBondOrder(bond));
             }
             totalVariance += calculateVariance(bondOrders);
        }
    }

    // Avoid division by zero if no heteroatoms meet the criteria
    return (heteroatomCount > 0 && totalVariance > 0.0) ? (totalVariance / heteroatomCount) : 0.0;
}


} // namespace descriptors
} // namespace desfact
