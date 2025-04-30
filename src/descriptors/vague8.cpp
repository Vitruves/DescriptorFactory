#include "descriptors/vague8.hpp"
#include "utils.hpp" // Includes common headers and logger
#include <GraphMol/GraphMol.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmartsWrite.h> // For MolToSmarts
#include <GraphMol/Descriptors/MolDescriptors.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/AtomIterators.h>
#include <GraphMol/BondIterators.h>
#include <GraphMol/RingInfo.h>
#include <GraphMol/PeriodicTable.h>
#include <GraphMol/PartialCharges/GasteigerCharges.h>
#include <GraphMol/MolTransforms/MolTransforms.h> // Needed? Maybe not for topological
#include <GraphMol/ChemTransforms/ChemTransforms.h> // For MolFragmenter
#include <GraphMol/Fingerprints/MorganFingerprints.h>
#include <GraphMol/Fingerprints/Fingerprints.h> // TopologicalTorsion Fingerprint maybe?
#include <DataStructs/ExplicitBitVect.h> // Needed for Morgan FP result

#include <numeric>
#include <cmath>
#include <vector>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <set>
#include <algorithm>
#include <limits>
#include <sstream>

// Re-use helpers from vague3 by including its implementation detail?
// Or better, move common helpers to utils.hpp/cpp or a dedicated helpers file.
// For now, let's redefine necessary helpers or slightly adapt them here.
#include "descriptors/sum.hpp" // Access base class static helpers
#include "descriptors/fractional.hpp" // Access base class static helpers


namespace desfact {
namespace descriptors {

// --- Anonymous Namespace for Helpers ---
namespace {

    // Basic atom/bond checks (reuse or redefine slightly)
    inline bool isHeteroatomV8(const RDKit::Atom* atom) {
        return atom->getAtomicNum() != 6 && atom->getAtomicNum() != 1;
    }
    inline bool isAtomInRingV8(const RDKit::Atom* atom) {
        return atom->getOwningMol().getRingInfo()->numAtomRings(atom->getIdx()) > 0;
    }
    inline unsigned int getAtomDegreeV8(const RDKit::Atom* atom) {
        return atom->getDegree();
    }
    inline bool isCarbonV8(const RDKit::Atom* atom) {
        return atom->getAtomicNum() == 6;
    }
    inline bool isTerminalAtomV8(const RDKit::Atom* atom) {
        return getAtomDegreeV8(atom) == 1;
    }
    inline bool isChiralCenterV8(const RDKit::Atom* atom) {
        return atom->getChiralTag() != RDKit::Atom::ChiralType::CHI_UNSPECIFIED;
    }
     inline double getBondOrderV8(const RDKit::Bond* bond) {
        switch (bond->getBondType()) {
            case RDKit::Bond::BondType::SINGLE: return 1.0;
            case RDKit::Bond::BondType::DOUBLE: return 2.0;
            case RDKit::Bond::BondType::TRIPLE: return 3.0;
            case RDKit::Bond::BondType::AROMATIC: return 1.5;
            default: return 1.0; // Default for other types like DATIVE, etc.
        }
    }

    // Electronegativity (reuse from sum.cpp's helpers via base class)
    double getAtomENV8(const RDKit::Atom* atom) {
        return SumDescriptor::getAtomElectronegativity(atom);
    }

     // Atomic Radii (reuse from sum.cpp's helpers via base class)
    double getAtomCovalentRadiusV8(const RDKit::Atom* atom) {
        return SumDescriptor::getAtomCovalentRadius(atom);
    }

    // Shortest path (BFS implementation)
    int getShortestPathV8(const RDKit::ROMol* mol, unsigned int atom1Idx, unsigned int atom2Idx) {
        if (atom1Idx == atom2Idx) return 0;
        unsigned int nAtoms = mol->getNumAtoms();
        if (atom1Idx >= nAtoms || atom2Idx >= nAtoms) return -1; // Invalid indices

        std::vector<int> distances(nAtoms, -1);
        std::queue<unsigned int> queue;

        distances[atom1Idx] = 0;
        queue.push(atom1Idx);

        while (!queue.empty()) {
            unsigned int current = queue.front();
            queue.pop();

            if (current == atom2Idx) return distances[current];

            for (const auto& nbr : mol->atomNeighbors(mol->getAtomWithIdx(current))) {
                unsigned int nbrIdx = nbr->getIdx();
                if (distances[nbrIdx] == -1) {
                    distances[nbrIdx] = distances[current] + 1;
                    queue.push(nbrIdx);
                }
            }
        }
        return -1; // Path not found (disconnected graph)
    }

    // All-pairs shortest path (Floyd-Warshall) - compute only if needed by multiple logic points
    std::vector<std::vector<int>> getAllPairsShortestPath(const RDKit::ROMol* mol) {
        unsigned int nAtoms = mol->getNumAtoms();
        const int INF = std::numeric_limits<int>::max() / 2; // Avoid overflow
        std::vector<std::vector<int>> distMatrix(nAtoms, std::vector<int>(nAtoms, INF));

        for (unsigned int i = 0; i < nAtoms; ++i) distMatrix[i][i] = 0;

        for (const auto& bond : mol->bonds()) {
            unsigned int u = bond->getBeginAtomIdx();
            unsigned int v = bond->getEndAtomIdx();
            distMatrix[u][v] = 1;
            distMatrix[v][u] = 1;
        }

        for (unsigned int k = 0; k < nAtoms; ++k) {
            for (unsigned int i = 0; i < nAtoms; ++i) {
                for (unsigned int j = 0; j < nAtoms; ++j) {
                    if (distMatrix[i][k] != INF && distMatrix[k][j] != INF &&
                        distMatrix[i][k] + distMatrix[k][j] < distMatrix[i][j]) {
                        distMatrix[i][j] = distMatrix[i][k] + distMatrix[k][j];
                    }
                }
            }
        }
         // Convert INF back to -1 for consistency with single path function
        for (unsigned int i = 0; i < nAtoms; ++i) {
            for (unsigned int j = 0; j < nAtoms; ++j) {
                if (distMatrix[i][j] == INF) distMatrix[i][j] = -1;
            }
        }

        return distMatrix;
    }


    // Calculate Shannon diversity index
    double calculateShannonEntropy(const std::vector<int>& counts) {
        double totalCount = std::accumulate(counts.begin(), counts.end(), 0.0);
        if (totalCount <= 1.0) return 0.0; // Entropy is 0 if 0 or 1 item

        double entropy = 0.0;
        for (int count : counts) {
            if (count > 0) {
                double p = static_cast<double>(count) / totalCount;
                entropy -= p * std::log2(p); // Use log base 2 for bits, or ln for nats
            }
        }
        return entropy;
    }
     // Overload for maps (string keys)
    double calculateShannonEntropy(const std::unordered_map<std::string, int>& counts) {
        double totalCount = 0;
        for(const auto& pair : counts) {
            totalCount += pair.second;
        }
        if (totalCount <= 1.0) return 0.0;

        double entropy = 0.0;
        for (const auto& pair : counts) {
            if (pair.second > 0) {
                double p = static_cast<double>(pair.second) / totalCount;
                entropy -= p * std::log2(p);
            }
        }
        return entropy;
    }
    // Overload for maps (int keys)
    double calculateShannonEntropy(const std::unordered_map<int, int>& counts) {
         double totalCount = 0;
        for(const auto& pair : counts) {
            totalCount += pair.second;
        }
        if (totalCount <= 1.0) return 0.0;

        double entropy = 0.0;
        for (const auto& pair : counts) {
            if (pair.second > 0) {
                double p = static_cast<double>(pair.second) / totalCount;
                entropy -= p * std::log2(p);
            }
        }
        return entropy;
    }


    // Calculate Variance
    double calculateVarianceV8(const std::vector<double>& values) {
        size_t count = values.size();
        if (count <= 1) return 0.0;

        double sum = std::accumulate(values.begin(), values.end(), 0.0);
        double mean = sum / count;
        double squaredDiffSum = 0.0;
        for (double value : values) {
            double diff = value - mean;
            squaredDiffSum += diff * diff;
        }
        return squaredDiffSum / count; // Population variance
    }

    // Calculate Skewness
    double calculateSkewness(const std::vector<double>& values) {
        size_t n = values.size();
        if (n < 3) return 0.0; // Skewness undefined for < 3 points

        double sum = std::accumulate(values.begin(), values.end(), 0.0);
        double mean = sum / n;

        double m3 = 0.0; // Third central moment
        double m2 = 0.0; // Second central moment (variance numerator)
        for (double val : values) {
            double diff = val - mean;
            m2 += diff * diff;
            m3 += diff * diff * diff;
        }
        m2 /= n;
        m3 /= n;

        if (m2 < 1e-9) return 0.0; // Avoid division by zero if variance is negligible

        double stddev = std::sqrt(m2);
        return m3 / (stddev * stddev * stddev);
    }


    // Get Gasteiger partial charges
    bool getGasteigerCharges(const RDKit::ROMol* rdkMol, std::vector<double>& charges) {
        charges.assign(rdkMol->getNumAtoms(), 0.0);
        try {
            // Gasteiger charges require a writeable molecule, so we copy
            RDKit::RWMol nonConstMol(*rdkMol);
            RDKit::MolOps::sanitizeMol(nonConstMol); // Ensure sanitization
            RDKit::computeGasteigerCharges(nonConstMol, 12, true); // 12 iterations, throw on error

            for (unsigned int i = 0; i < rdkMol->getNumAtoms(); ++i) {
                 if (nonConstMol.getAtomWithIdx(i)->hasProp(RDKit::common_properties::_GasteigerCharge)) {
                    charges[i] = nonConstMol.getAtomWithIdx(i)->getProp<double>(RDKit::common_properties::_GasteigerCharge);
                 } else {
                     // This case shouldn't happen if computeGasteigerCharges succeeds
                     globalLogger.warning("Gasteiger charge calculation succeeded but charge property missing for atom " + std::to_string(i)); // Fixed: warn -> warning
                     charges[i] = 0.0; // Assign default
                 }
            }
            return true;
        } catch (const std::exception& e) {
            globalLogger.debug("Failed to compute Gasteiger charges: " + std::string(e.what()));
            return false;
        } catch (...) {
            globalLogger.debug("Failed to compute Gasteiger charges due to unknown error.");
            return false;
        }
    }

    // Define functional groups (example list, can be expanded)
    const std::vector<std::pair<std::string, std::string>>& getFunctionalGroupSMARTS() {
        static const std::vector<std::pair<std::string, std::string>> functionalGroupSmarts = {
            {"Alcohol", "[#6][OX2H]"}, // Aliphatic alcohol C-OH
            {"Phenol", "[#6;a][OX2H]"}, // Phenol Ar-OH
            {"Ether", "[OD2]([#6])[#6]"}, // Ether C-O-C
            {"Thiol", "[#6][SX2H]"}, // Thiol C-SH
            {"Sulfide", "[SD2]([#6])[#6]"}, // Sulfide C-S-C
            {"PrimaryAmine", "[NX3;H2;!$(NC=[O,S,N])][#6]"}, // Primary amine R-NH2
            {"SecondaryAmine", "[NX3;H1;!$(NC=[O,S,N])]([#6])[#6]"}, // Secondary amine R-NH-R'
            {"TertiaryAmine", "[NX3;H0;!$(NC=[O,S,N])]([#6])([#6])[#6]"}, // Tertiary amine R3N
            {"Aldehyde", "[CX3H1](=O)[#6]"}, // Aldehyde R-CHO
            {"Ketone", "[#6][CX3](=O)[#6]"}, // Ketone R-CO-R'
            {"CarboxylicAcid", "[CX3](=O)[OX2H]"}, // Carboxylic acid R-COOH
            {"Ester", "[#6][CX3](=O)[OX2][#6]"}, // Ester R-COO-R'
            {"Amide", "[CX3](=O)[NX3]"}, // Amide R-CONH2, R-CONHR', R-CONR'R''
            {"Nitro", "[NX3+](=O)[O-]"}, // Nitro R-NO2
            {"Halogen", "[F,Cl,Br,I]"} // Halogen (attached to Carbon)
        };
        return functionalGroupSmarts;
    }

    // Find functional group instances (returns map: group name -> list of atom indices)
    std::unordered_map<std::string, std::vector<std::vector<int>>> findFunctionalGroups(const RDKit::ROMol* mol) {
        std::unordered_map<std::string, std::vector<std::vector<int>>> groupInstances;
        const auto& smartsList = getFunctionalGroupSMARTS();

        for (const auto& pair : smartsList) {
            const std::string& name = pair.first;
            const std::string& smarts = pair.second;
            std::unique_ptr<RDKit::ROMol> query(RDKit::SmartsToMol(smarts));
            if (!query) {
                globalLogger.warning("Failed to parse SMARTS for functional group: " + name); // Fixed: warn -> warning
                continue;
            }

            std::vector<std::vector<std::pair<int, int>>> matches;
            unsigned int nMatches = RDKit::SubstructMatch(*mol, *query, matches);

            if (nMatches > 0) {
                std::vector<std::vector<int>>& instances = groupInstances[name];
                for (const auto& match : matches) {
                    std::vector<int> atomIndices;
                    for (const auto& atomPair : match) {
                        atomIndices.push_back(atomPair.second);
                    }
                    if (!atomIndices.empty()) {
                         // Sort indices to handle potential duplicates/order differences if needed
                         std::sort(atomIndices.begin(), atomIndices.end());
                         instances.push_back(atomIndices);
                    }
                }
                 // Optional: Remove duplicate matches if SMARTS could overlap identically
                 std::sort(instances.begin(), instances.end());
                 instances.erase(std::unique(instances.begin(), instances.end()), instances.end());
            }
        }
        return groupInstances;
    }


    // Identify ring systems (connected components of ring atoms)
    std::vector<std::vector<int>> getRingSystems(const RDKit::ROMol* mol) {
        const RDKit::RingInfo* ringInfo = mol->getRingInfo();
        if (!ringInfo->isInitialized() || ringInfo->numRings() == 0) {
            return {};
        }

        std::vector<std::vector<int>> ringSystems;
        std::vector<bool> visited(mol->getNumAtoms(), false);
        std::vector<int> ringAtomIndices;

        // Collect all atoms that are part of any ring
         for (unsigned int i = 0; i < mol->getNumAtoms(); ++i) {
             if (ringInfo->numAtomRings(i) > 0) {
                 ringAtomIndices.push_back(i);
             }
         }

        // Perform BFS/DFS only on ring atoms to find connected components
        for (int startAtomIdx : ringAtomIndices) {
            if (visited[startAtomIdx]) continue;

            std::vector<int> currentSystem;
            std::queue<int> queue;
            queue.push(startAtomIdx);
            visited[startAtomIdx] = true;

            while (!queue.empty()) {
                int u = queue.front();
                queue.pop();
                currentSystem.push_back(u);

                for (const auto& neighbor : mol->atomNeighbors(mol->getAtomWithIdx(u))) {
                    int v = neighbor->getIdx();
                    // Check if neighbor is also a ring atom and not visited
                    if (ringInfo->numAtomRings(v) > 0 && !visited[v]) {
                         // Check if they share a bond that is in a ring
                         const RDKit::Bond* bond = mol->getBondBetweenAtoms(u, v);
                         if (bond && ringInfo->numBondRings(bond->getIdx()) > 0) {
                            visited[v] = true;
                            queue.push(v);
                         }
                    }
                }
            }
            if (!currentSystem.empty()) {
                ringSystems.push_back(currentSystem);
            }
        }
        return ringSystems;
    }


     // Find chain fragments (non-ring atoms) - similar to vague3 helper
    std::vector<std::vector<int>> findChainFragmentsV8(const RDKit::ROMol* mol) {
        std::vector<std::vector<int>> chainFragments;
        std::vector<bool> visited(mol->getNumAtoms(), false);
        const RDKit::RingInfo* ringInfo = mol->getRingInfo();
        bool ringsInitialized = ringInfo && ringInfo->isInitialized();

        // Mark all ring atoms as visited
        for (unsigned int i = 0; i < mol->getNumAtoms(); ++i) {
            if (ringsInitialized && ringInfo->numAtomRings(i) > 0) {
                visited[i] = true;
            }
        }

        // Find chain fragments
        for (unsigned int i = 0; i < mol->getNumAtoms(); ++i) {
            if (visited[i]) continue;

            std::vector<int> fragment;
            std::queue<int> queue;
            queue.push(i);
            visited[i] = true;

            while (!queue.empty()) {
                int current = queue.front();
                queue.pop();
                fragment.push_back(current);

                for (const auto& nbr : mol->atomNeighbors(mol->getAtomWithIdx(current))) {
                    int nbrIdx = nbr->getIdx();
                    // Add neighbor if it's not visited AND not a ring atom
                    if (!visited[nbrIdx] && (!ringsInitialized || ringInfo->numAtomRings(nbrIdx) == 0)) {
                        queue.push(nbrIdx);
                        visited[nbrIdx] = true;
                    }
                }
            }
            if (!fragment.empty()) {
                 chainFragments.push_back(fragment);
            }
        }
        return chainFragments;
    }

    // Define electronegative atoms (adjust as needed)
    bool isElectronegativeV8(const RDKit::Atom* atom) {
        int atomicNum = atom->getAtomicNum();
        // Common EN atoms: O, N, F, Cl, Br, S, P, I
        return (atomicNum == 8 || atomicNum == 7 || atomicNum == 9 ||
                atomicNum == 17 || atomicNum == 35 || atomicNum == 16 ||
                atomicNum == 15 || atomicNum == 53);
    }

     // Define hydrophobic atoms (Carbon primarily, maybe others?)
    bool isHydrophobicV8(const RDKit::Atom* atom) {
        // Primarily carbon atoms not bonded to highly electronegative atoms?
        // Simplification: just non-polar carbons and maybe sulfur?
        int atomicNum = atom->getAtomicNum();
         if (atomicNum == 6) {
             // Check neighbors? A simple definition: Carbon.
             return true;
         }
         // Maybe non-charged Sulfur?
         // if (atomicNum == 16 && atom->getFormalCharge() == 0) return true;
         return false;
    }

     // Get neighbors within a certain radius
     std::vector<int> getNeighborsWithinRadius(const RDKit::ROMol* mol, int startAtomIdx, int radius) {
         std::unordered_set<int> neighbors;
         std::queue<std::pair<int, int>> queue; // {atomIdx, distance}
         std::unordered_set<int> visited;

         queue.push({startAtomIdx, 0});
         visited.insert(startAtomIdx);

         while (!queue.empty()) {
             auto currentPair = queue.front();
             queue.pop();
             int currentAtomIdx = currentPair.first;
             int currentDist = currentPair.second;

             if (currentDist > 0) { // Don't include the start atom itself
                 neighbors.insert(currentAtomIdx);
             }

             if (currentDist < radius) {
                 for (const auto& nbr : mol->atomNeighbors(mol->getAtomWithIdx(currentAtomIdx))) {
                     int nbrIdx = nbr->getIdx();
                     if (visited.find(nbrIdx) == visited.end()) {
                         visited.insert(nbrIdx);
                         queue.push({nbrIdx, currentDist + 1});
                     }
                 }
             }
         }
         return std::vector<int>(neighbors.begin(), neighbors.end());
     }

     // Define standard substructures (example list)
     const std::vector<std::string>& getStandardSubstructuresSMARTS() {
         static const std::vector<std::string> smarts = {
             "c1ccccc1", // Benzene
             "[CX3](=O)[OX2H]", // Carboxylic acid
             "[CX3](=O)[OX2][#6]", // Ester
             "[CX3](=O)[NX3]", // Amide
             "[#6][OX2H]", // Alcohol
             "[NX3;H2]", // Primary amine
             "[NX3;H1]", // Secondary amine
             "[NX3;H0]", // Tertiary amine
             "[SX2H]", // Thiol
             "[CX3]=[CX3]", // Alkene
             "[CX2]#[CX2]", // Alkyne
             "[F,Cl,Br,I]", // Halogen
             "[!#6;!#1;!H0]~[CH3]", // Heteroatom-Methyl
             "O=C-N", // Amide fragment
             "C-O-C", // Ether linkage
             "c:n:c", // Pyridine-like N
             "S", // Sulfur atom
             "P" // Phosphorus atom
             // Add more common fragments
         };
         return smarts;
     }

} // end anonymous namespace


// --- Descriptor Implementations ---

// TopologicalChargeDistributionSkewness
TopologicalChargeDistributionSkewnessDescriptor::TopologicalChargeDistributionSkewnessDescriptor()
    : Descriptor("TopChargeDistSkew", "Skewness of partial charges across topological distances from center") {}

std::variant<double, int, std::string> TopologicalChargeDistributionSkewnessDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol || rdkMol->getNumAtoms() == 0) return 0.0;

    std::vector<double> charges;
    if (!getGasteigerCharges(rdkMol, charges)) {
        // Return 0 if charges cannot be computed
        return 0.0;
    }

    if (rdkMol->getNumAtoms() < 3) return 0.0; // Skewness needs >= 3 points

    // Calculate topological distances from atom 0 (or geometric center if coords available, but let's use atom 0)
    std::vector<double> chargeDistValues;
    for (unsigned int i = 0; i < rdkMol->getNumAtoms(); ++i) {
        if (std::abs(charges[i]) > 1e-4) { // Only consider atoms with non-negligible charge
             int dist = getShortestPathV8(rdkMol, 0, i);
             if (dist >= 0) {
                 // Weight charge by distance? Or just look at distribution of charges at distances?
                 // Let's calculate skewness of the charges themselves, potentially weighted later.
                 // For now, just the skewness of the charge values.
                 chargeDistValues.push_back(charges[i]);
             }
        }
    }

    if (chargeDistValues.size() < 3) return 0.0;

    return calculateSkewness(chargeDistValues);
}


// ElementNeighborhoodDiversity
ElementNeighborhoodDiversityDescriptor::ElementNeighborhoodDiversityDescriptor()
    : Descriptor("ElemNeighborDiversity", "Avg Shannon entropy of element types in immediate neighborhood (radius 1)") {}

std::variant<double, int, std::string> ElementNeighborhoodDiversityDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol || rdkMol->getNumAtoms() == 0) return 0.0;

    double totalEntropy = 0.0;
    int atomCount = 0;

    for (const auto& atom : rdkMol->atoms()) {
        std::unordered_map<int, int> neighborElements;
        int degree = atom->getDegree();
        if (degree == 0) continue; // Skip isolated atoms

        atomCount++;
        for (const auto& nbr : rdkMol->atomNeighbors(atom)) {
            neighborElements[nbr->getAtomicNum()]++;
        }

        std::vector<int> counts;
        for(const auto& pair : neighborElements) {
            counts.push_back(pair.second);
        }
        totalEntropy += calculateShannonEntropy(counts);
    }

    if (atomCount == 0) return 0.0;
    return totalEntropy / atomCount;
}

// BondOrderAlternationPattern
BondOrderAlternationPatternDescriptor::BondOrderAlternationPatternDescriptor()
    : Descriptor("BondOrderAlternation", "Fraction of adjacent bond pairs with alternating order (e.g., S-D, D-S, A-A)") {}

std::variant<double, int, std::string> BondOrderAlternationPatternDescriptor::calculate(const Molecule& mol) const {
     if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol || rdkMol->getNumAtoms() < 3) return 0.0;

    int alternatingPairs = 0;
    int totalPairs = 0;

    for (const auto& atom : rdkMol->atoms()) {
        if (atom->getDegree() < 2) continue; // Need at least two bonds to form a pair

        std::vector<const RDKit::Bond*> bonds;
        for(const auto& bond : rdkMol->atomBonds(atom)) {
            bonds.push_back(bond);
        }

        // Iterate through pairs of bonds connected to this atom
        for (size_t i = 0; i < bonds.size(); ++i) {
            for (size_t j = i + 1; j < bonds.size(); ++j) {
                totalPairs++;
                RDKit::Bond::BondType type1 = bonds[i]->getBondType();
                RDKit::Bond::BondType type2 = bonds[j]->getBondType();

                bool isAlternating = false;
                if ((type1 == RDKit::Bond::SINGLE && type2 == RDKit::Bond::DOUBLE) ||
                    (type1 == RDKit::Bond::DOUBLE && type2 == RDKit::Bond::SINGLE) ||
                    (type1 == RDKit::Bond::AROMATIC && type2 == RDKit::Bond::AROMATIC) || // A-A can be alternating
                    (type1 == RDKit::Bond::SINGLE && type2 == RDKit::Bond::AROMATIC) || // Consider S-A / A-S?
                    (type1 == RDKit::Bond::AROMATIC && type2 == RDKit::Bond::SINGLE) ||
                     (type1 == RDKit::Bond::DOUBLE && type2 == RDKit::Bond::TRIPLE) || // D-T / T-D ?
                    (type1 == RDKit::Bond::TRIPLE && type2 == RDKit::Bond::DOUBLE) )
                 {
                    // Basic definition: Single/Double or Aromatic/Aromatic seems most common
                     isAlternating = true;
                }


                // Simplified check: Are bond orders different (excluding triple maybe)?
                 double order1 = getBondOrderV8(bonds[i]);
                 double order2 = getBondOrderV8(bonds[j]);
                 // Let's count if orders are 1 and 2, or 1.5 and 1.5
                 if ((order1 == 1.0 && order2 == 2.0) || (order1 == 2.0 && order2 == 1.0) ||
                     (order1 == 1.5 && order2 == 1.5) ){
                     alternatingPairs++;
                 }
                 // Alternative: Fraction of non-single bonds adjacent to single bonds?
            }
        }
    }

    if (totalPairs == 0) return 0.0;
    return static_cast<double>(alternatingPairs) / totalPairs;
}


// FunctionalGroupConnectivityMatrix -> AvgShortestPathBetweenDiffFuncGroups
AvgShortestPathBetweenDiffFuncGroupsDescriptor::AvgShortestPathBetweenDiffFuncGroupsDescriptor()
    : Descriptor("AvgPathDiffFuncGroups", "Average shortest path between atoms of different functional group types") {}

std::variant<double, int, std::string> AvgShortestPathBetweenDiffFuncGroupsDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol || rdkMol->getNumAtoms() < 2) return 0.0;

    auto funcGroups = findFunctionalGroups(rdkMol);
    if (funcGroups.size() < 2) return 0.0; // Need at least two different group types

    std::vector<std::string> groupNames;
    std::vector<std::vector<int>> groupAtomLists; // Flattened list of atoms per group type
     for(auto const& [name, instances] : funcGroups) {
         if (!instances.empty()) {
             groupNames.push_back(name);
             std::vector<int> atoms;
             for(const auto& instance : instances) {
                 atoms.insert(atoms.end(), instance.begin(), instance.end());
             }
             // Remove duplicate atoms within the same group type list
             std::sort(atoms.begin(), atoms.end());
             atoms.erase(std::unique(atoms.begin(), atoms.end()), atoms.end());
             groupAtomLists.push_back(atoms);
         }
     }

     if (groupAtomLists.size() < 2) return 0.0;

    double totalMinDistanceSum = 0;
    int pairCount = 0;
    auto distMatrix = getAllPairsShortestPath(rdkMol); // Precompute distances

    for (size_t i = 0; i < groupAtomLists.size(); ++i) {
        for (size_t j = i + 1; j < groupAtomLists.size(); ++j) {
            int minDistance = std::numeric_limits<int>::max();
            bool foundPath = false;

            for (int atomI : groupAtomLists[i]) {
                for (int atomJ : groupAtomLists[j]) {
                    int dist = (atomI < distMatrix.size() && atomJ < distMatrix[atomI].size()) ? distMatrix[atomI][atomJ] : -1;
                    if (dist != -1 && dist < minDistance) {
                        minDistance = dist;
                        foundPath = true;
                    }
                }
            }

            if (foundPath && minDistance > 0) { // Only consider groups that are connected
                totalMinDistanceSum += minDistance;
                pairCount++;
            }
        }
    }

    if (pairCount == 0) return 0.0;
    return totalMinDistanceSum / pairCount;
}

// HeteroatomClusterTopology -> AvgHeteroClusterSize
AvgHeteroClusterSizeDescriptor::AvgHeteroClusterSizeDescriptor()
    : Descriptor("AvgHeteroClusterSize", "Average size (number of atoms) of connected heteroatom clusters") {}

std::variant<double, int, std::string> AvgHeteroClusterSizeDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol || rdkMol->getNumAtoms() == 0) return 0.0;

    std::vector<std::vector<int>> clusters;
    std::vector<bool> visited(rdkMol->getNumAtoms(), false);
    double totalSize = 0;
    int clusterCount = 0;

    for (unsigned int i = 0; i < rdkMol->getNumAtoms(); ++i) {
        const RDKit::Atom* atom = rdkMol->getAtomWithIdx(i);
        if (visited[i] || !isHeteroatomV8(atom)) continue;

        std::vector<int> currentCluster;
        std::queue<int> queue;
        queue.push(i);
        visited[i] = true;

        while (!queue.empty()) {
            int u = queue.front();
            queue.pop();
            currentCluster.push_back(u);

            for (const auto& nbr : rdkMol->atomNeighbors(rdkMol->getAtomWithIdx(u))) {
                int v = nbr->getIdx();
                if (!visited[v] && isHeteroatomV8(nbr)) {
                    visited[v] = true;
                    queue.push(v);
                }
            }
        }

        if (!currentCluster.empty()) {
             // Consider minimum cluster size? Original vague3 used minSize=2. Let's average all.
            totalSize += currentCluster.size();
            clusterCount++;
        }
    }

    if (clusterCount == 0) return 0.0;
    return totalSize / clusterCount;
}

// RingSubstitutionPatternCode -> EntropyOfRingSubstCountDist
EntropyOfRingSubstCountDistDescriptor::EntropyOfRingSubstCountDistDescriptor()
    : Descriptor("EntropyRingSubstCount", "Entropy of the distribution of substitution counts per ring atom") {}

std::variant<double, int, std::string> EntropyOfRingSubstCountDistDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    const RDKit::RingInfo* ringInfo = rdkMol->getRingInfo();
    if (!ringInfo || !ringInfo->isInitialized() || ringInfo->numRings() == 0) return 0.0;

    std::unordered_map<int, int> substitutionCounts; // Map: count -> frequency
    int totalRingAtoms = 0;

    for (unsigned int i = 0; i < rdkMol->getNumAtoms(); ++i) {
        if (ringInfo->numAtomRings(i) > 0) {
            totalRingAtoms++;
            const RDKit::Atom* atom = rdkMol->getAtomWithIdx(i);
            int substituentCount = 0;
            for (const auto& nbr : rdkMol->atomNeighbors(atom)) {
                if (ringInfo->numAtomRings(nbr->getIdx()) == 0) {
                    substituentCount++;
                }
            }
            substitutionCounts[substituentCount]++;
        }
    }

    if (totalRingAtoms == 0) return 0.0;

    std::vector<int> counts;
     for(const auto& pair: substitutionCounts) {
         counts.push_back(pair.second);
     }

    return calculateShannonEntropy(counts);
}


// FunctionalGroupDistanceHistogram -> EntropyOfFuncGroupDistances
EntropyOfFuncGroupDistancesDescriptor::EntropyOfFuncGroupDistancesDescriptor()
    : Descriptor("EntropyFuncGroupDists", "Entropy of histogram of shortest path distances between functional group atoms") {}

std::variant<double, int, std::string> EntropyOfFuncGroupDistancesDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol || rdkMol->getNumAtoms() < 2) return 0.0;

    auto funcGroups = findFunctionalGroups(rdkMol);
    if (funcGroups.empty()) return 0.0;

    std::vector<int> allFuncGroupAtoms;
     for(auto const& [name, instances] : funcGroups) {
         for(const auto& instance : instances) {
            allFuncGroupAtoms.insert(allFuncGroupAtoms.end(), instance.begin(), instance.end());
         }
     }
    // Remove duplicates
    std::sort(allFuncGroupAtoms.begin(), allFuncGroupAtoms.end());
    allFuncGroupAtoms.erase(std::unique(allFuncGroupAtoms.begin(), allFuncGroupAtoms.end()), allFuncGroupAtoms.end());

    if (allFuncGroupAtoms.size() < 2) return 0.0;

    std::unordered_map<int, int> distanceCounts; // distance -> count
    auto distMatrix = getAllPairsShortestPath(rdkMol); // Precompute distances

    for (size_t i = 0; i < allFuncGroupAtoms.size(); ++i) {
        for (size_t j = i + 1; j < allFuncGroupAtoms.size(); ++j) {
            int atomI = allFuncGroupAtoms[i];
            int atomJ = allFuncGroupAtoms[j];
            int dist = (atomI < distMatrix.size() && atomJ < distMatrix[atomI].size()) ? distMatrix[atomI][atomJ] : -1;

            if (dist > 0) { // Only consider connected atoms, ignore distance 0
                distanceCounts[dist]++;
            }
        }
    }

    if (distanceCounts.empty()) return 0.0;

    std::vector<int> counts;
    for(const auto& pair : distanceCounts) {
        counts.push_back(pair.second);
    }

    return calculateShannonEntropy(counts);
}


// ChainBranchingRecursivePattern -> BranchPointsInChains / ChainAtoms
ChainBranchingPatternDescriptor::ChainBranchingPatternDescriptor()
    : FractionalDescriptor("ChainBranchingPattern", "Branch points (degree >= 3) in non-ring chains / total atoms in non-ring chains") {}

std::variant<double, int, std::string> ChainBranchingPatternDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol || rdkMol->getNumAtoms() == 0) return 0.0;

    auto chainFragments = findChainFragmentsV8(rdkMol);
    if (chainFragments.empty()) return 0.0;

    int branchPoints = 0;
    int totalChainAtoms = 0;

    for (const auto& fragment : chainFragments) {
        totalChainAtoms += fragment.size();
        for (int atomIdx : fragment) {
            const RDKit::Atom* atom = rdkMol->getAtomWithIdx(atomIdx);
             // Degree within the whole molecule, or just within the chain?
             // Let's use degree within the whole molecule to define branching.
            if (atom->getDegree() >= 3) {
                branchPoints++;
            }
        }
    }

    if (totalChainAtoms == 0) return 0.0;
    return static_cast<double>(branchPoints) / totalChainAtoms;
}


// ElectronegativeAtomNeighborhoodScore -> Sum(CountENNeighborsRadius2)
ElectronegativeAtomNeighborhoodScoreDescriptor::ElectronegativeAtomNeighborhoodScoreDescriptor()
    : SumDescriptor("ENAtomNeighborScore", "Sum over all atoms of count of electronegative neighbors within radius 2") {}

std::variant<double, int, std::string> ElectronegativeAtomNeighborhoodScoreDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0; // Return int
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol || rdkMol->getNumAtoms() == 0) return 0;

    int totalScore = 0;
     for (unsigned int i = 0; i < rdkMol->getNumAtoms(); ++i) {
         std::vector<int> neighbors = getNeighborsWithinRadius(rdkMol, i, 2);
         int count = 0;
         for (int nbrIdx : neighbors) {
             if (isElectronegativeV8(rdkMol->getAtomWithIdx(nbrIdx))) {
                 count++;
             }
         }
         totalScore += count;
     }

    return totalScore;
}

// TopologicalChargeSeparationIndex -> MeanDistBetweenFormalCharges
TopologicalChargeSeparationIndexDescriptor::TopologicalChargeSeparationIndexDescriptor()
    : Descriptor("TopoChargeSepIndex", "Mean shortest path distance between atoms with non-zero formal charge") {}

std::variant<double, int, std::string> TopologicalChargeSeparationIndexDescriptor::calculate(const Molecule& mol) const {
     if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol || rdkMol->getNumAtoms() < 2) return 0.0;

    std::vector<int> chargedAtoms;
    for (const auto& atom : rdkMol->atoms()) {
        if (atom->getFormalCharge() != 0) {
            chargedAtoms.push_back(atom->getIdx());
        }
    }

    if (chargedAtoms.size() < 2) return 0.0;

    double totalDistance = 0;
    int pairCount = 0;
    auto distMatrix = getAllPairsShortestPath(rdkMol); // Precompute

    for (size_t i = 0; i < chargedAtoms.size(); ++i) {
        for (size_t j = i + 1; j < chargedAtoms.size(); ++j) {
            int atomI = chargedAtoms[i];
            int atomJ = chargedAtoms[j];
             int dist = (atomI < distMatrix.size() && atomJ < distMatrix[atomI].size()) ? distMatrix[atomI][atomJ] : -1;

            if (dist > 0) { // Only connected pairs
                totalDistance += dist;
                pairCount++;
            }
        }
    }

    if (pairCount == 0) return 0.0; // All charged atoms might be in disconnected components
    return totalDistance / pairCount;
}


// RingSystemConnectivityPattern -> InterRingSystemBonds / NumRingSystems
RingSystemConnectivityPatternDescriptor::RingSystemConnectivityPatternDescriptor()
    : Descriptor("RingSysConnectPattern", "Bonds connecting different ring systems / Number of ring systems") {}

std::variant<double, int, std::string> RingSystemConnectivityPatternDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol || rdkMol->getNumAtoms() == 0) return 0.0;

    std::vector<std::vector<int>> ringSystems = getRingSystems(rdkMol);
    int numSystems = ringSystems.size();
    if (numSystems <= 1) return 0.0; // Need at least two systems for connections

    // Create a mapping from atom index to system index
    std::vector<int> atomToSystem(rdkMol->getNumAtoms(), -1);
    for (int sysIdx = 0; sysIdx < numSystems; ++sysIdx) {
        for (int atomIdx : ringSystems[sysIdx]) {
            atomToSystem[atomIdx] = sysIdx;
        }
    }

    int interSystemBonds = 0;
    for (const auto& bond : rdkMol->bonds()) {
        int atom1 = bond->getBeginAtomIdx();
        int atom2 = bond->getEndAtomIdx();
        int sys1 = atomToSystem[atom1];
        int sys2 = atomToSystem[atom2];

        // Check if both atoms are in ring systems AND they are in different systems
        if (sys1 != -1 && sys2 != -1 && sys1 != sys2) {
            interSystemBonds++;
        }
    }

    // The definition asks for bonds / systems
    return static_cast<double>(interSystemBonds) / numSystems;
}


// HeteroatomPositionEntropy -> EntropyOfHeteroatomIndicesInCanonSmiles
HeteroatomPositionEntropyDescriptor::HeteroatomPositionEntropyDescriptor()
    : StringDescriptor("HeteroatomPosEntropy", "Entropy of positional indices of heteroatoms in canonical SMILES") {}

std::variant<double, int, std::string> HeteroatomPositionEntropyDescriptor::calculateFromSmiles(const std::string& smiles) const {
    if (smiles.empty()) return 0.0;

    std::unordered_map<int, int> positionCounts; // Position -> Count (should always be 1)
    std::vector<int> heteroatomPositions;

    for (int i = 0; i < smiles.length(); ++i) {
        char c = smiles[i];
        // Basic check for common heteroatom element symbols (case-sensitive in SMILES)
        // More robust: parse SMILES to ROMol and check atomic numbers?
        // StringDescriptor base implies calculation should be from string only.
        // This is tricky without full SMILES parsing. Let's approximate.
        // Consider O, N, S, P, F, Cl, Br, I
         if (c == 'O' || c == 'N' || c == 'S' || c == 'P' || c == 'F' || c == 'I') {
             heteroatomPositions.push_back(i);
         } else if (c == 'C' && (i + 1) < smiles.length() && smiles[i+1] == 'l') { // Chlorine 'Cl'
              heteroatomPositions.push_back(i);
              i++; // Skip 'l'
         } else if (c == 'B' && (i + 1) < smiles.length() && smiles[i+1] == 'r') { // Bromine 'Br'
             heteroatomPositions.push_back(i);
             i++; // Skip 'r'
         }
         // This simple check misses bracket atoms like [nH], [O-], etc.
         // A truly robust implementation might require temporary Mol parsing here.
         // Let's stick to the simple version based on StringDescriptor principle.
    }

     if (heteroatomPositions.empty()) return 0.0;

     // Calculate entropy of the distribution of positions.
     // If all positions are unique, entropy might just reflect the number of heteroatoms.
     // Let's try entropy of the *intervals* between heteroatoms.
     if (heteroatomPositions.size() < 2) return 0.0; // Need at least two for intervals

     std::vector<int> intervals;
     for(size_t i = 1; i < heteroatomPositions.size(); ++i) {
         intervals.push_back(heteroatomPositions[i] - heteroatomPositions[i-1]);
     }

     // Count frequencies of intervals
     std::unordered_map<int, int> intervalCounts;
     for (int interval : intervals) {
         intervalCounts[interval]++;
     }

     if (intervalCounts.empty()) return 0.0;

     std::vector<int> counts;
     for(const auto& pair : intervalCounts) {
         counts.push_back(pair.second);
     }

    return calculateShannonEntropy(counts);
}


// BondOrderTransitionFrequency -> EntropyOfAdjacentBondTypePairs
BondOrderTransitionFrequencyDescriptor::BondOrderTransitionFrequencyDescriptor()
    : Descriptor("BondOrderTransFreq", "Entropy of frequencies of adjacent bond type pairs (A-B, B-C)") {}

std::variant<double, int, std::string> BondOrderTransitionFrequencyDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol || rdkMol->getNumAtoms() < 3) return 0.0;

    std::unordered_map<std::string, int> transitionCounts;

    for (const auto& atomB : rdkMol->atoms()) {
        if (atomB->getDegree() < 2) continue;

        std::vector<const RDKit::Bond*> bondsFromB;
         for(const auto& bond : rdkMol->atomBonds(atomB)) {
             bondsFromB.push_back(bond);
         }

        for (size_t i = 0; i < bondsFromB.size(); ++i) {
            for (size_t j = i + 1; j < bondsFromB.size(); ++j) {
                 const RDKit::Bond* bondAB = bondsFromB[i];
                 const RDKit::Bond* bondBC = bondsFromB[j];

                 // Get bond types as strings for map key
                 std::string type1 = std::to_string(static_cast<int>(bondAB->getBondType()));
                 std::string type2 = std::to_string(static_cast<int>(bondBC->getBondType()));

                 // Ensure canonical order for the pair (e.g., smaller number first)
                 std::string key = (type1 < type2) ? (type1 + "-" + type2) : (type2 + "-" + type1);
                 transitionCounts[key]++;
            }
        }
    }

    if (transitionCounts.empty()) return 0.0;

    return calculateShannonEntropy(transitionCounts);
}


// AtomNeighborhoodElectronegativityGradient -> Avg(AtomEN - AvgNeighborEN)
AtomNeighborhoodElectronegativityGradientDescriptor::AtomNeighborhoodElectronegativityGradientDescriptor()
    : SumDescriptor("AtomNeighborENGrad", "Average difference between atom EN and its neighbors' average EN") {}

std::variant<double, int, std::string> AtomNeighborhoodElectronegativityGradientDescriptor::calculate(const Molecule& mol) const {
     if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol || rdkMol->getNumAtoms() == 0) return 0.0;

    double totalGradientSum = 0.0;
    int atomCount = 0;

    for (const auto& atom : rdkMol->atoms()) {
        double atomEN = getAtomENV8(atom);
        double neighborENSum = 0.0;
        int degree = atom->getDegree();

        if (degree == 0) continue;
        atomCount++;

        for (const auto& nbr : rdkMol->atomNeighbors(atom)) {
            neighborENSum += getAtomENV8(nbr);
        }
        double avgNeighborEN = neighborENSum / degree;
        totalGradientSum += std::abs(atomEN - avgNeighborEN); // Use absolute difference? Or signed? Let's use absolute.
    }

    if (atomCount == 0) return 0.0;
    return totalGradientSum / atomCount;
}


// ChainLengthDistributionEntropy
ChainLengthDistributionEntropyDescriptor::ChainLengthDistributionEntropyDescriptor()
    : Descriptor("ChainLenDistEntropy", "Entropy of the distribution of non-ring chain lengths") {}

std::variant<double, int, std::string> ChainLengthDistributionEntropyDescriptor::calculate(const Molecule& mol) const {
     if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol || rdkMol->getNumAtoms() == 0) return 0.0;

    auto chainFragments = findChainFragmentsV8(rdkMol);
    if (chainFragments.empty()) return 0.0;

    std::unordered_map<int, int> lengthCounts; // length -> count
    for (const auto& fragment : chainFragments) {
        if (!fragment.empty()) {
             lengthCounts[fragment.size()]++;
        }
    }

    if (lengthCounts.empty()) return 0.0;

    std::vector<int> counts;
     for(const auto& pair : lengthCounts) {
         counts.push_back(pair.second);
     }

    return calculateShannonEntropy(counts);
}

// RingSubstitutionSymmetry -> Avg(DistinctSubstituentTypesPerRing)
RingSubstitutionSymmetryDescriptor::RingSubstitutionSymmetryDescriptor()
    : Descriptor("RingSubstSymmetry", "Inverse of average number of distinct substituent SMARTS per ring") {}

std::variant<double, int, std::string> RingSubstitutionSymmetryDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    const RDKit::RingInfo* ringInfo = rdkMol->getRingInfo();
     if (!ringInfo || !ringInfo->isInitialized()) {
         RDKit::RWMol nonConstMolCopy(*rdkMol);
         RDKit::MolOps::findSSSR(nonConstMolCopy);
         ringInfo = nonConstMolCopy.getRingInfo();
         if (!ringInfo || !ringInfo->isInitialized()) {
              globalLogger.debug("Ring info could not be initialized for RingSubstSymmetry.");
              return 0.0;
         }
     }
     if (ringInfo->numRings() == 0) return 0.0;

    double totalDistinctSubstCount = 0;
    int ringCount = ringInfo->numRings();
    RDKit::RWMol nonConstMol(*rdkMol); // Use a copy if modifications happen (though FP shouldn't modify)

    for (const auto& ringBondsIdx : ringInfo->bondRings()) {
        std::unordered_set<int> ringAtomsIdx;
        for(int bondIdx : ringBondsIdx) {
            const RDKit::Bond* bond = nonConstMol.getBondWithIdx(bondIdx);
            ringAtomsIdx.insert(bond->getBeginAtomIdx());
            ringAtomsIdx.insert(bond->getEndAtomIdx());
        }

        std::unordered_set<std::string> distinctSubstituentsForRing;

         try {
            for (int atomIdx : ringAtomsIdx) {
                const RDKit::Atom* atom = nonConstMol.getAtomWithIdx(atomIdx);
                for (const auto& nbr : nonConstMol.atomNeighbors(atom)) {
                    int nbrIdx = nbr->getIdx();
                    // If neighbor is NOT in the current ring system
                    if (ringAtomsIdx.find(nbrIdx) == ringAtomsIdx.end()) {
                         const RDKit::Bond* bondToSubst = nonConstMol.getBondBetweenAtoms(atomIdx, nbrIdx);
                         if (!bondToSubst) continue;

                         // --- Fingerprint Fix Start ---
                         // Use Morgan fingerprint (ECFP-like) radius 1 centered on the neighbor atom
                         const unsigned int radius = 1;
                         const unsigned int nBits = 2048; // Standard size
                         bool useChirality = false;
                         bool useBondTypes = true;

                         // Fix: Use uint32_t instead of int for atomIndicesToUse
                         std::vector<std::uint32_t> atomIndicesToUse = {static_cast<std::uint32_t>(nbrIdx)};

                         ExplicitBitVect* fp = RDKit::MorganFingerprints::getFingerprintAsBitVect(
                            nonConstMol, 
                            radius,
                            nBits,
                            &atomIndicesToUse, // Center on this atom
                            nullptr, // invariants
                            useChirality,
                            useBondTypes,
                            false // useCounts = false for BitVect
                        );
                        // --- Fingerprint Fix End ---

                        if(fp) {
                            std::string fpStr = fp->toString(); // Convert bit vector to string
                            distinctSubstituentsForRing.insert(fpStr);
                            delete fp; // Clean up the pointer
                        }
                    }
                }
            }
            totalDistinctSubstCount += distinctSubstituentsForRing.size();

         } catch(const std::exception& e) {
              globalLogger.debug("Exception processing substituents for a ring: " + std::string(e.what()));
         } catch(...) {
              globalLogger.debug("Unknown error processing substituents for a ring.");
         }
    }

    if (ringCount == 0) return 0.0;
    double avgDistinctSubsts = totalDistinctSubstCount / ringCount;
    return 1.0 / (1.0 + avgDistinctSubsts);
}


// HeteroatomSequencePatterns -> EntropyOfHeteroatomTriplets (A-B-C)
HeteroatomSequencePatternsDescriptor::HeteroatomSequencePatternsDescriptor()
    : Descriptor("HeteroatomSeqPatterns", "Entropy of element sequence triplets (A-B-C) involving at least one heteroatom") {}

std::variant<double, int, std::string> HeteroatomSequencePatternsDescriptor::calculate(const Molecule& mol) const {
     if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol || rdkMol->getNumAtoms() < 3) return 0.0;

    std::unordered_map<std::string, int> sequenceCounts;

     for (const auto& atomB : rdkMol->atoms()) {
         int atomBIdx = atomB->getIdx();
         if (atomB->getDegree() < 2) continue;

         std::vector<const RDKit::Atom*> neighbors;
          for(const auto& nbr : rdkMol->atomNeighbors(atomB)) {
              neighbors.push_back(nbr);
          }

         for (size_t i = 0; i < neighbors.size(); ++i) {
             for (size_t j = i + 1; j < neighbors.size(); ++j) {
                 const RDKit::Atom* atomA = neighbors[i];
                 const RDKit::Atom* atomC = neighbors[j];
                 int elemA = atomA->getAtomicNum();
                 int elemB = atomB->getAtomicNum();
                 int elemC = atomC->getAtomicNum();

                 // Check if at least one is heteroatom
                 if (isHeteroatomV8(atomA) || isHeteroatomV8(atomB) || isHeteroatomV8(atomC)) {
                     // Create canonical sequence string (e.g., sorted elements)
                     std::vector<int> elems = {elemA, elemB, elemC};
                     // Sort based on A-B-C path? No, just the triplet. Sort numerically.
                     std::sort(elems.begin(), elems.end());
                     std::string seqKey = std::to_string(elems[0]) + "-" + std::to_string(elems[1]) + "-" + std::to_string(elems[2]);
                     sequenceCounts[seqKey]++;
                 }
             }
         }
     }

    if (sequenceCounts.empty()) return 0.0;

    return calculateShannonEntropy(sequenceCounts);
}


// FunctionalGroupIsolationTopology -> Avg(MinDistToOtherFuncGroup)
FunctionalGroupIsolationTopologyDescriptor::FunctionalGroupIsolationTopologyDescriptor()
    : Descriptor("FuncGroupIsolationTopo", "Average minimum shortest path distance from each functional group instance to another instance") {}

std::variant<double, int, std::string> FunctionalGroupIsolationTopologyDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol || rdkMol->getNumAtoms() < 2) return 0.0;

     auto funcGroups = findFunctionalGroups(rdkMol);
     if (funcGroups.empty()) return 0.0;

     // Flatten into a list of all instances, storing their atoms
     std::vector<std::vector<int>> allInstancesAtoms;
     for(const auto& pair : funcGroups) {
         allInstancesAtoms.insert(allInstancesAtoms.end(), pair.second.begin(), pair.second.end());
     }

     if (allInstancesAtoms.size() < 2) return 0.0; // Need at least two instances total

     double totalMinDistSum = 0;
     int instanceCount = allInstancesAtoms.size();
     auto distMatrix = getAllPairsShortestPath(rdkMol); // Precompute

     for(size_t i = 0; i < allInstancesAtoms.size(); ++i) {
         int currentMinDist = std::numeric_limits<int>::max();
         bool foundOther = false;

         for (size_t j = 0; j < allInstancesAtoms.size(); ++j) {
             if (i == j) continue; // Don't compare an instance to itself

             // Find min distance between any atom in instance i and any atom in instance j
             int pairMinDist = std::numeric_limits<int>::max();
             bool foundPairPath = false;
             for (int atomI : allInstancesAtoms[i]) {
                 for (int atomJ : allInstancesAtoms[j]) {
                    int dist = (atomI < distMatrix.size() && atomJ < distMatrix[atomI].size()) ? distMatrix[atomI][atomJ] : -1;
                    if (dist != -1 && dist < pairMinDist) {
                        pairMinDist = dist;
                        foundPairPath = true;
                    }
                 }
             }

             if (foundPairPath && pairMinDist < currentMinDist) {
                  currentMinDist = pairMinDist;
                  foundOther = true;
             }
         }

         if (foundOther && currentMinDist > 0) { // Only add if another connected group was found
             totalMinDistSum += currentMinDist;
         } else {
             // If an instance is isolated or only connected to itself, how to count?
             // Decrement the count of instances we average over?
             // For now, let's assume isolated groups contribute 0 to the sum and reduce the count.
             instanceCount--;
         }
     }

     if (instanceCount <= 0) return 0.0;
     return totalMinDistSum / instanceCount;
}


// HydrophobicPatchConnectivity -> NumHydrophobicPatches / NumHydrophobicAtoms
HydrophobicPatchConnectivityDescriptor::HydrophobicPatchConnectivityDescriptor()
    : FractionalDescriptor("HydrophobicPatchConnect", "Number of connected hydrophobic patches / Total number of hydrophobic atoms") {}

std::variant<double, int, std::string> HydrophobicPatchConnectivityDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol || rdkMol->getNumAtoms() == 0) return 0.0;

    std::vector<int> hydrophobicAtoms;
     for(const auto& atom : rdkMol->atoms()) {
         if (isHydrophobicV8(atom)) {
             hydrophobicAtoms.push_back(atom->getIdx());
         }
     }

     if (hydrophobicAtoms.empty()) return 0.0;

     int patchCount = 0;
     std::vector<bool> visited(rdkMol->getNumAtoms(), false); // Use full size for indexing convenience

     for (int startAtomIdx : hydrophobicAtoms) {
         if (visited[startAtomIdx]) continue;

         patchCount++;
         std::queue<int> queue;
         queue.push(startAtomIdx);
         visited[startAtomIdx] = true;

         while (!queue.empty()) {
             int u = queue.front();
             queue.pop();

             for (const auto& nbr : rdkMol->atomNeighbors(rdkMol->getAtomWithIdx(u))) {
                 int v = nbr->getIdx();
                 // Add to patch if neighbor is hydrophobic and not visited
                 if (!visited[v] && isHydrophobicV8(nbr)) {
                      visited[v] = true;
                      queue.push(v);
                 }
             }
         }
     }

    // The descriptor asks for Patches / HydrophobicAtoms
    return static_cast<double>(patchCount) / hydrophobicAtoms.size();
}


// BondTopologicalEnvironmentFingerprint -> CountUniqueBondEnvFingerprints
BondTopologicalEnvironmentFingerprintDescriptor::BondTopologicalEnvironmentFingerprintDescriptor()
    : Descriptor("BondTopoEnvFP", "Count of unique bond topological environments (Morgan-like, radius 1)") {}

std::variant<double, int, std::string> BondTopologicalEnvironmentFingerprintDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0; // Return int count
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol || rdkMol->getNumBonds() == 0) return 0;

    std::set<std::vector<uint32_t>> uniqueFps; // Use set for automatic uniqueness

    // Morgan fingerprints are atom-centered. How to get bond-centered?
    // Option 1: Generate FP for each atom, combine FPs of bonded atoms?
    // Option 2: Generate FP centered on atom A including bond AB, centered on B including bond BA, combine?
    // Option 3: Use RDKit's AtomPair or TopologicalTorsion fingerprints which consider bonds?
    // Let's try a simplified Morgan-like approach for bonds:
    // For bond A-B, collect invariants for A, B, and their neighbors (excluding B from A's neighbors, A from B's neighbors).
    // Hash these invariants.

    std::set<std::string> uniqueBondEnvHashes;

    for (const auto& bond : rdkMol->bonds()) {
        std::vector<std::string> invariants;
        const RDKit::Atom* atomA = bond->getBeginAtom();
        const RDKit::Atom* atomB = bond->getEndAtom();
        int idxA = atomA->getIdx();
        int idxB = atomB->getIdx();

        // Bond invariant
        invariants.push_back("B" + std::to_string(static_cast<int>(bond->getBondType())));

        // Atom A invariants
        invariants.push_back("A" + std::to_string(atomA->getAtomicNum()));
        invariants.push_back("D" + std::to_string(atomA->getDegree()));
        invariants.push_back("C" + std::to_string(atomA->getFormalCharge()));
        invariants.push_back("R" + std::to_string(isAtomInRingV8(atomA)));
        // Atom A neighbors (excluding B)
        for(const auto& nbrA : rdkMol->atomNeighbors(atomA)) {
            if (nbrA->getIdx() != idxB) {
                 invariants.push_back("N" + std::to_string(nbrA->getAtomicNum()));
                 // Add bond type to neighbor?
                 const RDKit::Bond* bondToNbr = rdkMol->getBondBetweenAtoms(idxA, nbrA->getIdx());
                 if(bondToNbr) invariants.push_back("b" + std::to_string(static_cast<int>(bondToNbr->getBondType())));
            }
        }

         // Atom B invariants
        invariants.push_back("A" + std::to_string(atomB->getAtomicNum()));
        invariants.push_back("D" + std::to_string(atomB->getDegree()));
        invariants.push_back("C" + std::to_string(atomB->getFormalCharge()));
        invariants.push_back("R" + std::to_string(isAtomInRingV8(atomB)));
         // Atom B neighbors (excluding A)
        for(const auto& nbrB : rdkMol->atomNeighbors(atomB)) {
            if (nbrB->getIdx() != idxA) {
                 invariants.push_back("N" + std::to_string(nbrB->getAtomicNum()));
                 const RDKit::Bond* bondToNbr = rdkMol->getBondBetweenAtoms(idxB, nbrB->getIdx());
                 if(bondToNbr) invariants.push_back("b" + std::to_string(static_cast<int>(bondToNbr->getBondType())));
            }
        }

        // Create canonical string hash
        std::sort(invariants.begin(), invariants.end());
        std::stringstream ss;
        for(const auto& inv : invariants) ss << inv << ";";
        uniqueBondEnvHashes.insert(ss.str());
    }


    return static_cast<int>(uniqueBondEnvHashes.size());
}


// ChiralCenterTopologicalDistribution -> VarianceOfDistancesBetweenChiralCenters
ChiralCenterTopologicalDistributionDescriptor::ChiralCenterTopologicalDistributionDescriptor()
    : Descriptor("ChiralCenterTopoDist", "Variance of shortest path distances between chiral centers") {}

std::variant<double, int, std::string> ChiralCenterTopologicalDistributionDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol || rdkMol->getNumAtoms() < 2) return 0.0;

    std::vector<int> chiralCenters;
    for (const auto& atom : rdkMol->atoms()) {
        if (isChiralCenterV8(atom)) {
            chiralCenters.push_back(atom->getIdx());
        }
    }

    if (chiralCenters.size() < 2) return 0.0; // Variance requires at least 2 points

    std::vector<double> distances;
    auto distMatrix = getAllPairsShortestPath(rdkMol); // Precompute

    for (size_t i = 0; i < chiralCenters.size(); ++i) {
        for (size_t j = i + 1; j < chiralCenters.size(); ++j) {
            int atomI = chiralCenters[i];
            int atomJ = chiralCenters[j];
             int dist = (atomI < distMatrix.size() && atomJ < distMatrix[atomI].size()) ? distMatrix[atomI][atomJ] : -1;

            if (dist > 0) { // Only connected pairs
                distances.push_back(static_cast<double>(dist));
            }
        }
    }

    if (distances.size() < 1) return 0.0; // Need at least one distance pair for variance

    return calculateVarianceV8(distances);
}


// RingFusionPatternCode -> NumSharedRingBonds / NumRings
RingFusionPatternCodeDescriptor::RingFusionPatternCodeDescriptor()
    : Descriptor("RingFusionCode", "Number of bonds shared between rings / Total number of rings") {}

std::variant<double, int, std::string> RingFusionPatternCodeDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    const RDKit::RingInfo* ringInfo = rdkMol->getRingInfo();
    if (!ringInfo || !ringInfo->isInitialized() || ringInfo->numRings() == 0) return 0.0;

    int numRings = ringInfo->numRings();
    int sharedBonds = 0;

    for (unsigned int bondIdx = 0; bondIdx < rdkMol->getNumBonds(); ++bondIdx) {
        if (ringInfo->numBondRings(bondIdx) > 1) { // Bond is in more than one ring
            sharedBonds++;
        }
    }

    // Avoid division by zero (though numRings is > 0 here)
    return static_cast<double>(sharedBonds) / numRings;
}


// ElectronegativityTopologicalMoment -> Sum(EN_i * dist(i, 0)) / Sum(EN_i)
ElectronegativityTopologicalMomentDescriptor::ElectronegativityTopologicalMomentDescriptor()
    : SumDescriptor("ENTopoMoment", "Sum(EN_i * dist(i, atom0)) / Sum(EN_i)") {}

std::variant<double, int, std::string> ElectronegativityTopologicalMomentDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol || rdkMol->getNumAtoms() == 0) return 0.0;

    double weightedDistSum = 0.0;
    double totalEN = 0.0;
    auto distMatrix = getAllPairsShortestPath(rdkMol); // Get all distances from atom 0


    for (unsigned int i = 0; i < rdkMol->getNumAtoms(); ++i) {
        const RDKit::Atom* atom = rdkMol->getAtomWithIdx(i);
        double en = getAtomENV8(atom);
        if (en > 0) { // Only consider atoms with defined EN
             int dist = (0 < distMatrix.size() && i < distMatrix[0].size()) ? distMatrix[0][i] : -1;
             if (dist != -1) { // Only if connected to atom 0
                weightedDistSum += en * dist;
                totalEN += en;
             }
        }
    }

    if (totalEN < 1e-6) return 0.0; // Avoid division by zero if total EN is negligible
    return weightedDistSum / totalEN;
}


// AtomicRadiiVarianceInNeighborhoods -> Avg(Variance(NeighborRadii))
AtomicRadiiVarianceInNeighborhoodsDescriptor::AtomicRadiiVarianceInNeighborhoodsDescriptor()
    : SumDescriptor("AtomRadiiNeighborVar", "Average variance of covalent radii in immediate atom neighborhoods") {}

std::variant<double, int, std::string> AtomicRadiiVarianceInNeighborhoodsDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol || rdkMol->getNumAtoms() == 0) return 0.0;

    double totalVarianceSum = 0.0;
    int atomCount = 0;

     for (const auto& atom : rdkMol->atoms()) {
         int degree = atom->getDegree();
         if (degree < 2) continue; // Variance needs at least 2 neighbors

         atomCount++;
         std::vector<double> neighborRadii;
         for(const auto& nbr : rdkMol->atomNeighbors(atom)) {
             neighborRadii.push_back(getAtomCovalentRadiusV8(nbr));
         }

         totalVarianceSum += calculateVarianceV8(neighborRadii);
     }

    if (atomCount == 0) return 0.0;
    return totalVarianceSum / atomCount;
}


// HeteroatomBondingPatternCode -> CountUniqueHeteroatomBondingCodes
HeteroatomBondingPatternCodeDescriptor::HeteroatomBondingPatternCodeDescriptor()
    : Descriptor("HeteroatomBondPattern", "Count of unique bonding patterns for heteroatoms") {}

std::variant<double, int, std::string> HeteroatomBondingPatternCodeDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0; // Return int count
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol || rdkMol->getNumAtoms() == 0) return 0;

    std::set<std::string> uniquePatterns;

    for (const auto& atom : rdkMol->atoms()) {
        if (!isHeteroatomV8(atom)) continue;

        std::vector<std::string> bondCodes;
        for (const auto& bond : rdkMol->atomBonds(atom)) {
             const RDKit::Atom* neighbor = bond->getOtherAtom(atom);
             // Code: BondType-NeighborElement
             std::string code = std::to_string(static_cast<int>(bond->getBondType())) + "-" + std::to_string(neighbor->getAtomicNum());
             bondCodes.push_back(code);
        }
        // Add atom's own properties? Element, Charge?
        std::string atomCode = "A" + std::to_string(atom->getAtomicNum()) + "C" + std::to_string(atom->getFormalCharge());
        bondCodes.push_back(atomCode);


        // Sort codes to create canonical pattern string
        std::sort(bondCodes.begin(), bondCodes.end());
        std::stringstream ss;
        for(const auto& code : bondCodes) ss << code << ";";
        uniquePatterns.insert(ss.str());
    }

    return static_cast<int>(uniquePatterns.size());
}


// SubstructureFrequencySpectrum -> EntropyOfSubstructureFrequencies
SubstructureFrequencySpectrumDescriptor::SubstructureFrequencySpectrumDescriptor()
    : Descriptor("SubstructFreqEntropy", "Entropy of frequencies of predefined common substructures") {}

std::variant<double, int, std::string> SubstructureFrequencySpectrumDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol || rdkMol->getNumAtoms() == 0) return 0.0;

    const auto& smartsList = getStandardSubstructuresSMARTS();
    std::unordered_map<std::string, int> frequencyCounts;
    bool matchFound = false;

    for (const std::string& smarts : smartsList) {
        std::unique_ptr<RDKit::ROMol> query(RDKit::SmartsToMol(smarts));
        if (!query) continue; // Skip invalid SMARTS

         std::vector<std::vector<std::pair<int, int>>> matches;
         try {
            unsigned int nMatches = RDKit::SubstructMatch(*rdkMol, *query, matches, true, true); // uniquify=true, recursionPossible=true
            if (nMatches > 0) {
                frequencyCounts[smarts] = nMatches;
                matchFound = true;
            }
         } catch(...) {
             globalLogger.debug("Substructure match failed for SMARTS: " + smarts);
             // Continue with other SMARTS
         }
    }

    if (!matchFound || frequencyCounts.empty()) return 0.0;

    std::vector<int> counts;
    for(const auto& pair : frequencyCounts) {
        counts.push_back(pair.second);
    }

    return calculateShannonEntropy(counts);
}


// FormalChargeNeighborhoodPattern -> CountUniqueFormalChargeNeighborPatterns(R=1)
FormalChargeNeighborhoodPatternDescriptor::FormalChargeNeighborhoodPatternDescriptor()
    : Descriptor("FormalChargeNeighborPat", "Count of unique formal charge patterns in immediate neighborhoods") {}

std::variant<double, int, std::string> FormalChargeNeighborhoodPatternDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0; // Return int count
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol || rdkMol->getNumAtoms() == 0) return 0;

    std::set<std::string> uniquePatterns;

    for (const auto& atom : rdkMol->atoms()) {
        std::vector<int> neighborCharges;
        // Include central atom's charge? Yes.
        neighborCharges.push_back(atom->getFormalCharge());

        for (const auto& nbr : rdkMol->atomNeighbors(atom)) {
            neighborCharges.push_back(nbr->getFormalCharge());
        }

        // Sort charges for canonical pattern
        std::sort(neighborCharges.begin(), neighborCharges.end());

        std::stringstream ss;
         for(size_t i=0; i < neighborCharges.size(); ++i) {
             ss << neighborCharges[i] << (i == neighborCharges.size() - 1 ? "" : ",");
         }
         uniquePatterns.insert(ss.str());
    }

    return static_cast<int>(uniquePatterns.size());
}


} // namespace descriptors
} // namespace desfact
