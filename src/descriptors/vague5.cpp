#include "descriptors/vague5.hpp"
#include "utils.hpp" // For globalLogger
#include <GraphMol/GraphMol.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/Atom.h>
#include <GraphMol/Bond.h>
#include <GraphMol/AtomIterators.h>
#include <GraphMol/BondIterators.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/Descriptors/MolDescriptors.h>
#include <GraphMol/PartialCharges/GasteigerCharges.h>
#include <GraphMol/ChemTransforms/ChemTransforms.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/PeriodicTable.h>
#include <numeric>
#include <vector>
#include <cmath>
#include <map>
#include <set>
#include <deque>
#include <algorithm>
#include <limits>
#include <functional>
#include <GraphMol/RingInfo.h>

namespace desfact {
namespace descriptors {

// --- Helper Functions (Anonymous Namespace) ---
namespace {

// Re-introduce minimal needed element properties locally
namespace ElementPropertiesLocal {
    // Pauling electronegativity values
    static const std::map<int, double> electronegativity = {
        {1, 2.20},  // H
        {6, 2.55},  // C
        {7, 3.04},  // N
        {8, 3.44},  // O
        {9, 3.98},  // F
        {15, 2.19}, // P
        {16, 2.58}, // S
        {17, 3.16}, // Cl
        {35, 2.96}, // Br
        {53, 2.66}  // I
        // Add others if needed by placeholders later
    };
    // Helper to get property safely
    double getProperty(const std::map<int, double>& propertyMap, int atomicNum, double defaultValue) {
        auto it = propertyMap.find(atomicNum);
        return (it != propertyMap.end()) ? it->second : defaultValue;
    }
} // namespace ElementPropertiesLocal


// Placeholder for complex calculations - replace with actual logic
double placeholderCalculation(const RDKit::ROMol& mol, const std::string& descName) {
    // globalLogger.debug("Placeholder calculation for " + descName);
    // Simple placeholder: return hash of SMILES modulo 100
    try {
        std::string smiles = RDKit::MolToSmiles(mol);
        size_t hash_val = std::hash<std::string>{}(smiles);
        return static_cast<double>(hash_val % 100);
    } catch (...) {
        return std::numeric_limits<double>::quiet_NaN();
    }
}

int placeholderIntCalculation(const RDKit::ROMol& mol, const std::string& descName) {
    // globalLogger.debug("Placeholder calculation for " + descName);
    // Simple placeholder: return atom count
    try {
         return mol.getNumAtoms();
    } catch (...) {
        return -1;
    }
}

// Basic element check
bool isPolarAtom(const RDKit::Atom* atom) {
    int atomicNum = atom->getAtomicNum();
    return atomicNum == 7 || atomicNum == 8 || atomicNum == 9 || atomicNum == 15 || atomicNum == 16 || atomicNum == 17 || atomicNum == 35 || atomicNum == 53; // N, O, F, P, S, Cl, Br, I
}

bool isHydrophobicAtom(const RDKit::Atom* atom) {
    if (atom->getAtomicNum() != 6) return false; // Must be carbon
    // Use getAtomNeighbors iterator
    const RDKit::ROMol& mol = atom->getOwningMol();
    for (const auto& neighbor : mol.atomNeighbors(atom)) {
         if (neighbor->getAtomicNum() != 6 && neighbor->getAtomicNum() != 1) {
             return false; // Bonded to something other than C or H
         }
    }
    return true;
}

double getAtomENPauling(const RDKit::Atom* atom) {
    // Use local map
    return ElementPropertiesLocal::getProperty(ElementPropertiesLocal::electronegativity, atom->getAtomicNum(), 0.0);
}

// Get longest path in a graph represented by adjacency list
int getLongestPath(int startNode, int n, const std::vector<std::vector<int>>& adj) {
    std::vector<int> dist(n, -1);
    std::deque<int> q;

    dist[startNode] = 0;
    q.push_back(startNode);
    int maxDist = 0;
    int farthestNode = startNode;

    while (!q.empty()) {
        int u = q.front();
        q.pop_front();

        for (int v : adj[u]) {
            if (dist[v] == -1) {
                dist[v] = dist[u] + 1;
                if (dist[v] > maxDist) {
                    maxDist = dist[v];
                    farthestNode = v;
                }
                q.push_back(v);
            }
        }
    }

    std::fill(dist.begin(), dist.end(), -1);
    dist[farthestNode] = 0;
    q.push_back(farthestNode);
    maxDist = 0;

     while (!q.empty()) {
        int u = q.front();
        q.pop_front();

        for (int v : adj[u]) {
            if (dist[v] == -1) {
                dist[v] = dist[u] + 1;
                maxDist = std::max(maxDist, dist[v]);
                q.push_back(v);
            }
        }
    }
    return maxDist;
}

// Find connected components
void findComponentsDFS(int u, int componentId, const std::vector<std::vector<int>>& adj, std::vector<int>& componentMap, std::vector<bool>& visited) {
    visited[u] = true;
    componentMap[u] = componentId;
    for (int v : adj[u]) {
        if (!visited[v]) {
            findComponentsDFS(v, componentId, adj, componentMap, visited);
        }
    }
}

} // anonymous namespace


// --- Implementations ---

// 1. LongestContinuousSp2FragmentLength
LongestContinuousSp2FragmentLength::LongestContinuousSp2FragmentLength()
    : Descriptor("LongestContinuousSp2FragmentLength", "Length (number of atoms) of the longest continuous path of sp2 hybridized atoms") {}

std::variant<double, int, std::string> LongestContinuousSp2FragmentLength::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0;
    const RDKit::ROMol& rdkMol = *mol.getMolecule();

    std::vector<int> sp2_indices;
    std::map<int, int> atom_idx_to_sp2_idx; // Map original atom index to index within sp2 subset
    for (const auto& atom : rdkMol.atoms()) {
        if (atom->getHybridization() == RDKit::Atom::SP2) {
            atom_idx_to_sp2_idx[atom->getIdx()] = sp2_indices.size();
            sp2_indices.push_back(atom->getIdx());
        }
    }

    if (sp2_indices.empty()) return 0;

    int n_sp2 = sp2_indices.size();
    std::vector<std::vector<int>> adj(n_sp2);
    for (int i = 0; i < n_sp2; ++i) {
        const RDKit::Atom* atom1 = rdkMol.getAtomWithIdx(sp2_indices[i]);
        // Use getAtomNeighbors iterator
        for (const auto& neighbor : rdkMol.atomNeighbors(atom1)) {
            if (neighbor->getHybridization() == RDKit::Atom::SP2) {
                // Check if neighbor is in our sp2 list
                auto map_it = atom_idx_to_sp2_idx.find(neighbor->getIdx());
                if(map_it != atom_idx_to_sp2_idx.end()){
                    int neighbor_sp2_idx = map_it->second;
                    if (neighbor_sp2_idx > i) { // Avoid duplicate edges and self-loops
                        adj[i].push_back(neighbor_sp2_idx);
                        adj[neighbor_sp2_idx].push_back(i);
                    }
                }
            }
        }
    }

    int max_len = 0;
    std::vector<bool> visited(n_sp2, false);
    for (int i = 0; i < n_sp2; ++i) {
        if (!visited[i]) {
            // Find all nodes in the current component
            std::vector<int> current_component_nodes;
            std::vector<int> component_map(n_sp2, -1); // Map sp2_idx to index within component
            std::vector<std::vector<int>> component_adj;
            std::deque<int> q;

            q.push_back(i);
            visited[i] = true;
            int component_node_idx = 0;

            while (!q.empty()) {
                int u = q.front();
                q.pop_front();

                component_map[u] = component_node_idx++;
                current_component_nodes.push_back(u);
                component_adj.resize(component_node_idx); // Ensure size

                for (int v : adj[u]) {
                    if (!visited[v]) {
                        visited[v] = true;
                        q.push_back(v);
                    }
                }
            }

             // Build adjacency list for the current component only
            for(int node_sp2_idx_in_map : current_component_nodes) { // Use the index from the map key (original sp2 index)
                int u_comp_idx = component_map[node_sp2_idx_in_map]; // Get the component-local index
                 for(int neighbor_sp2_idx_in_adj : adj[node_sp2_idx_in_map]){ // Iterate neighbors using the original sp2 index
                     // Check if the neighbor is part of the current component using the component_map and visited array
                     if(component_map[neighbor_sp2_idx_in_adj] != -1){ // Ensure neighbor is in this component
                         int v_comp_idx = component_map[neighbor_sp2_idx_in_adj];
                         // Add edge only once using component-local indices
                         if(v_comp_idx > u_comp_idx){
                             component_adj[u_comp_idx].push_back(v_comp_idx);
                             component_adj[v_comp_idx].push_back(u_comp_idx);
                         }
                     }
                 }
            }

            if (!current_component_nodes.empty()) {
                 if (!component_adj.empty()) {
                    int component_longest_path = getLongestPath(0, component_adj.size(), component_adj);
                    max_len = std::max(max_len, component_longest_path + 1); // Path length is edges, size is atoms
                 } else if (current_component_nodes.size() == 1) {
                     max_len = std::max(max_len, 1); // Single sp2 atom counts as length 1
                 }
            }
        }
    }

    return max_len;
}

// 2. ChainBranchDepthMaximum
ChainBranchDepthMaximum::ChainBranchDepthMaximum()
    : Descriptor("ChainBranchDepthMaximum", "Maximum depth of branching from the main chain (heuristic)") {}

std::variant<double, int, std::string> ChainBranchDepthMaximum::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0;
    const RDKit::ROMol& rdkMol = *mol.getMolecule();
    // Placeholder - complex graph traversal needed. Using placeholder.
    return placeholderIntCalculation(rdkMol, getName());
}

// 3. ChainBranchingFrequency
ChainBranchingFrequency::ChainBranchingFrequency()
    : Descriptor("ChainBranchingFrequency", "Number of branch points (degree > 2, non-ring) per heavy atom") {}

std::variant<double, int, std::string> ChainBranchingFrequency::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    const RDKit::ROMol& rdkMol = *mol.getMolecule();
    int heavyAtomCount = rdkMol.getNumHeavyAtoms();
    if (heavyAtomCount == 0) return 0.0;

    RDKit::RingInfo* ringInfo = rdkMol.getRingInfo();
    if (!ringInfo->isInitialized()) {
         RDKit::MolOps::findSSSR(const_cast<RDKit::ROMol&>(rdkMol));
    }

    int branchPoints = 0;
    for (const auto& atom : rdkMol.atoms()) {
        if (atom->getAtomicNum() == 1) continue; // Skip hydrogens
        // Check if atom is in any ring using ringInfo
        if (atom->getDegree() > 2 && !(ringInfo->numAtomRings(atom->getIdx()) > 0)) {
            branchPoints++;
        }
    }
    return static_cast<double>(branchPoints) / heavyAtomCount;
}

// 4. LongestAlternatingSingleDoubleBondPath
LongestAlternatingSingleDoubleBondPath::LongestAlternatingSingleDoubleBondPath()
    : Descriptor("LongestAlternatingSingleDoubleBondPath", "Length of the longest path with alternating single/double bonds") {}

std::variant<double, int, std::string> LongestAlternatingSingleDoubleBondPath::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0;
    const RDKit::ROMol& rdkMol = *mol.getMolecule();
    // Placeholder - requires specific path finding algorithm. Using placeholder.
    return placeholderIntCalculation(rdkMol, getName());
}

// 5. AverageRingStrainProxy
AverageRingStrainProxy::AverageRingStrainProxy()
    : Descriptor("AverageRingStrainProxy", "Proxy for ring strain based on inverse average size of fused ring systems") {}

std::variant<double, int, std::string> AverageRingStrainProxy::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    const RDKit::ROMol& rdkMol = *mol.getMolecule();
     // Placeholder - Complex ring analysis needed. Using placeholder.
     return placeholderCalculation(rdkMol, getName());
}

// 6. FractionBackboneCarbon
FractionBackboneCarbon::FractionBackboneCarbon()
    : Descriptor("FractionBackboneCarbon", "Fraction of carbon atoms belonging to the molecular backbone (longest chain, heuristic)") {}

std::variant<double, int, std::string> FractionBackboneCarbon::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    const RDKit::ROMol& rdkMol = *mol.getMolecule();
    // Placeholder - Backbone identification needed. Using placeholder.
    return placeholderCalculation(rdkMol, getName());
}

// 7. TopologicalAsymmetryIndex
TopologicalAsymmetryIndex::TopologicalAsymmetryIndex()
    : Descriptor("TopologicalAsymmetryIndex", "Difference between the two longest paths starting from a central point (heuristic)") {}

std::variant<double, int, std::string> TopologicalAsymmetryIndex::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    const RDKit::ROMol& rdkMol = *mol.getMolecule();
    // Placeholder - Path finding needed. Using placeholder.
    return placeholderCalculation(rdkMol, getName());
}

// 8. FractionConjugatedBondsBackbone
FractionConjugatedBondsBackbone::FractionConjugatedBondsBackbone()
    : Descriptor("FractionConjugatedBondsBackbone", "Fraction of bonds in the backbone that are part of a conjugated system") {}

std::variant<double, int, std::string> FractionConjugatedBondsBackbone::calculate(const ::desfact::Molecule& mol) const {
     if (!mol.isValid() || !mol.getMolecule()) return 0.0;
     const RDKit::ROMol& rdkMol = *mol.getMolecule();
     // Placeholder - Requires backbone ID and conjugation analysis. Using placeholder.
     return placeholderCalculation(rdkMol, getName());
}

// 9. MaxChainNoHetero
MaxChainNoHetero::MaxChainNoHetero()
    : Descriptor("MaxChainNoHetero", "Length of the longest chain containing no heteroatoms") {}

std::variant<double, int, std::string> MaxChainNoHetero::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0;
    const RDKit::ROMol& rdkMol = *mol.getMolecule();
    // Placeholder - Requires path finding on carbon subgraph. Using placeholder.
    return placeholderIntCalculation(rdkMol, getName());
}

// 10. MeanRingSizeWeightedDegree
MeanRingSizeWeightedDegree::MeanRingSizeWeightedDegree()
    : Descriptor("MeanRingSizeWeightedDegree", "Mean ring size weighted by the degree of atoms participating in the ring") {}

std::variant<double, int, std::string> MeanRingSizeWeightedDegree::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    const RDKit::ROMol& rdkMol = *mol.getMolecule();
    // Placeholder - Requires detailed ring analysis. Using placeholder.
    return placeholderCalculation(rdkMol, getName());
}


// --- Electronic / Electronegativity Pattern ---

// 11. PolarClusterCount
PolarClusterCount::PolarClusterCount()
    : Descriptor("PolarClusterCount", "Number of connected clusters of 2 or more polar atoms (N, O, S, P, Halogens)") {}

std::variant<double, int, std::string> PolarClusterCount::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0;
    const RDKit::ROMol& rdkMol = *mol.getMolecule();

    std::vector<int> polar_indices;
    std::map<int, int> atom_idx_to_polar_idx;
    for (const auto& atom : rdkMol.atoms()) {
        if (isPolarAtom(atom)) {
            atom_idx_to_polar_idx[atom->getIdx()] = polar_indices.size();
            polar_indices.push_back(atom->getIdx());
        }
    }

    if (polar_indices.size() < 2) return 0;

    int n_polar = polar_indices.size();
    std::vector<std::vector<int>> adj(n_polar);
    for (int i = 0; i < n_polar; ++i) {
        const RDKit::Atom* atom1 = rdkMol.getAtomWithIdx(polar_indices[i]);
        // Use getAtomNeighbors iterator
        for (const auto& neighbor : rdkMol.atomNeighbors(atom1)) {
             if (isPolarAtom(neighbor)) {
                // Check if neighbor is in our polar list
                auto map_it = atom_idx_to_polar_idx.find(neighbor->getIdx());
                 if(map_it != atom_idx_to_polar_idx.end()){
                    int neighbor_polar_idx = map_it->second;
                    if (neighbor_polar_idx > i) { // Avoid duplicate edges
                        adj[i].push_back(neighbor_polar_idx);
                        adj[neighbor_polar_idx].push_back(i);
                    }
                 }
             }
        }
    }

    int clusterCount = 0;
    std::vector<bool> visited(n_polar, false);
    for (int i = 0; i < n_polar; ++i) {
        if (!visited[i]) {
            int componentSize = 0;
            std::deque<int> q;
            q.push_back(i);
            visited[i] = true;
            componentSize++;

            while(!q.empty()){
                int u = q.front();
                q.pop_front();
                for(int v : adj[u]){
                    if(!visited[v]){
                        visited[v] = true;
                        q.push_back(v);
                        componentSize++;
                    }
                }
            }
            if (componentSize >= 2) {
                clusterCount++;
            }
        }
    }
    return clusterCount;
}


// 12. MaxElectronegativityGradientChain
MaxElectronegativityGradientChain::MaxElectronegativityGradientChain()
    : Descriptor("MaxElectronegativityGradientChain", "Maximum difference in Pauling EN between adjacent atoms along any chain") {}

std::variant<double, int, std::string> MaxElectronegativityGradientChain::calculate(const ::desfact::Molecule& mol) const {
     if (!mol.isValid() || !mol.getMolecule()) return 0.0;
     const RDKit::ROMol& rdkMol = *mol.getMolecule();
     double maxGradient = 0.0;
      for (const auto& bond : rdkMol.bonds()) {
         double en1 = getAtomENPauling(bond->getBeginAtom());
         double en2 = getAtomENPauling(bond->getEndAtom());
         maxGradient = std::max(maxGradient, std::abs(en1 - en2));
      }
      return maxGradient;
}

// 13. FractionPolarizableHeavy
FractionPolarizableHeavy::FractionPolarizableHeavy()
    : Descriptor("FractionPolarizableHeavy", "Fraction of heavy atoms considered highly polarizable (e.g., S, P, Br, I)") {}

std::variant<double, int, std::string> FractionPolarizableHeavy::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    const RDKit::ROMol& rdkMol = *mol.getMolecule();
    int heavyAtomCount = rdkMol.getNumHeavyAtoms();
    if(heavyAtomCount == 0) return 0.0;

    int polarizableCount = 0;
    std::set<int> polarizableElements = {15, 16, 34, 35, 52, 53}; // P, S, Se, Br, Te, I
     for (const auto& atom : rdkMol.atoms()) {
         if (atom->getAtomicNum() > 1 && polarizableElements.count(atom->getAtomicNum())) {
             polarizableCount++;
         }
     }
    return static_cast<double>(polarizableCount) / heavyAtomCount;
}

// 14. NormalizedNetDipoleHeuristic
NormalizedNetDipoleHeuristic::NormalizedNetDipoleHeuristic()
    : Descriptor("NormalizedNetDipoleHeuristic", "Simple heuristic for net molecular dipole moment magnitude (normalized)") {}

std::variant<double, int, std::string> NormalizedNetDipoleHeuristic::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    const RDKit::ROMol& rdkMol = *mol.getMolecule();
    // Placeholder - Very complex, requires geometry approximation or Gasteiger charges maybe. Using placeholder.
    return placeholderCalculation(rdkMol, getName());
}

// 15. CumulativeENDifferenceBackbone
CumulativeENDifferenceBackbone::CumulativeENDifferenceBackbone()
    : Descriptor("CumulativeENDifferenceBackbone", "Sum of absolute EN differences across bonds in the presumed backbone") {}

std::variant<double, int, std::string> CumulativeENDifferenceBackbone::calculate(const ::desfact::Molecule& mol) const {
     if (!mol.isValid() || !mol.getMolecule()) return 0.0;
     const RDKit::ROMol& rdkMol = *mol.getMolecule();
     // Placeholder - Requires backbone identification. Using placeholder.
     return placeholderCalculation(rdkMol, getName());
}

// 16. HeavyAtomLocalChargeSkewness
HeavyAtomLocalChargeSkewness::HeavyAtomLocalChargeSkewness()
    : Descriptor("HeavyAtomLocalChargeSkewness", "Skewness of the distribution of formal charges on heavy atoms") {}

std::variant<double, int, std::string> HeavyAtomLocalChargeSkewness::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    const RDKit::ROMol& rdkMol = *mol.getMolecule();

    std::vector<double> charges;
    for (const auto& atom : rdkMol.atoms()) {
        if (atom->getAtomicNum() > 1) { // Heavy atoms only
            charges.push_back(static_cast<double>(atom->getFormalCharge()));
        }
    }

    if (charges.size() < 3) return 0.0; // Skewness requires at least 3 points

    double sum = std::accumulate(charges.begin(), charges.end(), 0.0);
    double mean = sum / charges.size();
    double variance_sum = 0.0;
    double skewness_sum = 0.0;

    for (double charge : charges) {
        variance_sum += (charge - mean) * (charge - mean);
        skewness_sum += (charge - mean) * (charge - mean) * (charge - mean);
    }

    double variance = variance_sum / charges.size();
    if (variance < 1e-9) return 0.0; // Avoid division by zero if all charges are the same

    double stddev = std::sqrt(variance);
    double skewness = skewness_sum / (charges.size() * stddev * stddev * stddev);

    return skewness;
}


// 17. ElectronegativeEndChainCount
ElectronegativeEndChainCount::ElectronegativeEndChainCount()
    : Descriptor("ElectronegativeEndChainCount", "Number of chain terminations ending with a highly electronegative atom (O, N, F, Cl)") {}

std::variant<double, int, std::string> ElectronegativeEndChainCount::calculate(const ::desfact::Molecule& mol) const {
     if (!mol.isValid() || !mol.getMolecule()) return 0;
     const RDKit::ROMol& rdkMol = *mol.getMolecule();
     int count = 0;
     std::set<int> enElements = {7, 8, 9, 17}; // N, O, F, Cl
     RDKit::RingInfo* ringInfo = rdkMol.getRingInfo(); // Get RingInfo once
      if (!ringInfo->isInitialized()) {
           RDKit::MolOps::findSSSR(const_cast<RDKit::ROMol&>(rdkMol));
      }

     for (const auto& atom : rdkMol.atoms()) {
          // Terminal atom: degree 1, not H, and not in a ring
          if (atom->getDegree() == 1 && atom->getAtomicNum() != 1) {
                 // Use RingInfo to check if atom is in any ring
                 if (!(ringInfo->numAtomRings(atom->getIdx()) > 0)) {
                    if (enElements.count(atom->getAtomicNum())) {
                        count++;
                    }
                 }
          }
     }
     return count;
}

// 18. FractionBondsLargeENGap
FractionBondsLargeENGap::FractionBondsLargeENGap()
    : Descriptor("FractionBondsLargeENGap", "Fraction of bonds connecting atoms with an EN difference > 1.5 Pauling units") {}

std::variant<double, int, std::string> FractionBondsLargeENGap::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    const RDKit::ROMol& rdkMol = *mol.getMolecule();
    int totalBonds = rdkMol.getNumBonds();
    if(totalBonds == 0) return 0.0;

    int largeGapBonds = 0;
    for (const auto& bond : rdkMol.bonds()) {
         double en1 = getAtomENPauling(bond->getBeginAtom());
         double en2 = getAtomENPauling(bond->getEndAtom());
         if (std::abs(en1 - en2) > 1.5) {
             largeGapBonds++;
         }
    }
    return static_cast<double>(largeGapBonds) / totalBonds;
}


// 19. SymmetryElectronegativeDistribution
SymmetryElectronegativeDistribution::SymmetryElectronegativeDistribution()
    : Descriptor("SymmetryElectronegativeDistribution", "Measure of symmetry in the distribution of electronegative atoms (heuristic)") {}

std::variant<double, int, std::string> SymmetryElectronegativeDistribution::calculate(const ::desfact::Molecule& mol) const {
     if (!mol.isValid() || !mol.getMolecule()) return 0.0;
     const RDKit::ROMol& rdkMol = *mol.getMolecule();
     // Placeholder - Complex heuristic needed. Using placeholder.
     return placeholderCalculation(rdkMol, getName());
}

// 20. FractionHeteroatomsTerminal
FractionHeteroatomsTerminal::FractionHeteroatomsTerminal()
    : Descriptor("FractionHeteroatomsTerminal", "Fraction of heteroatoms that are in terminal positions (degree 1)") {}

std::variant<double, int, std::string> FractionHeteroatomsTerminal::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    const RDKit::ROMol& rdkMol = *mol.getMolecule();

    int terminalHetero = 0;
    int totalHetero = 0;
    RDKit::RingInfo* ringInfo = rdkMol.getRingInfo(); // Get RingInfo once
     if (!ringInfo->isInitialized()) {
          RDKit::MolOps::findSSSR(const_cast<RDKit::ROMol&>(rdkMol));
     }

    for (const auto& atom : rdkMol.atoms()) {
        if (atom->getAtomicNum() != 1 && atom->getAtomicNum() != 6) { // Heteroatom
            totalHetero++;
            // Check if terminal (degree 1 and not in any ring)
             if (atom->getDegree() == 1) {
                  if (!(ringInfo->numAtomRings(atom->getIdx()) > 0)) {
                     terminalHetero++;
                 }
             }
        }
    }

    if (totalHetero == 0) return 0.0;
    return static_cast<double>(terminalHetero) / totalHetero;
}


// --- Conjugation / Resonance-Related ---

// 21. ConjugatedSystemSizeMax
ConjugatedSystemSizeMax::ConjugatedSystemSizeMax()
    : Descriptor("ConjugatedSystemSizeMax", "Number of atoms in the largest conjugated system") {}

std::variant<double, int, std::string> ConjugatedSystemSizeMax::calculate(const ::desfact::Molecule& mol) const {
     if (!mol.isValid() || !mol.getMolecule()) return 0;
     const RDKit::ROMol& rdkMol = *mol.getMolecule();
      // Placeholder - Requires identifying conjugated atoms and finding components. Using placeholder.
     return placeholderIntCalculation(rdkMol, getName());
}

// 22. FractionAromaticAtomsConnectedChains
FractionAromaticAtomsConnectedChains::FractionAromaticAtomsConnectedChains()
    : Descriptor("FractionAromaticAtomsConnectedChains", "Fraction of aromatic atoms that are connected to non-aromatic heavy atoms (chains)") {}

std::variant<double, int, std::string> FractionAromaticAtomsConnectedChains::calculate(const ::desfact::Molecule& mol) const {
     if (!mol.isValid() || !mol.getMolecule()) return 0.0;
     const RDKit::ROMol& rdkMol = *mol.getMolecule();
     int aromaticConnectedToChain = 0;
     int totalAromatic = 0;
      for(const auto& atom : rdkMol.atoms()){
          if(atom->getIsAromatic()){
              totalAromatic++;
              bool connectedToChain = false;
              // Use getAtomNeighbors iterator
               for (const auto& neighbor : rdkMol.atomNeighbors(atom)) {
                  if(neighbor->getAtomicNum() > 1 && !neighbor->getIsAromatic()){
                      connectedToChain = true;
                      break;
                  }
              }
              if(connectedToChain){
                  aromaticConnectedToChain++;
              }
          }
      }
      if(totalAromatic == 0) return 0.0;
      return static_cast<double>(aromaticConnectedToChain) / totalAromatic;
}

// 23. MeanDistanceAromaticSystems
MeanDistanceAromaticSystems::MeanDistanceAromaticSystems()
    : Descriptor("MeanDistanceAromaticSystems", "Mean topological distance between distinct aromatic systems (heuristic)") {}

std::variant<double, int, std::string> MeanDistanceAromaticSystems::calculate(const ::desfact::Molecule& mol) const {
     if (!mol.isValid() || !mol.getMolecule()) return 0.0;
     const RDKit::ROMol& rdkMol = *mol.getMolecule();
     // Placeholder - Requires identifying aromatic systems and path finding. Using placeholder.
     return placeholderCalculation(rdkMol, getName());
}

// 24. LengthLargestFullyConjugatedRingSystem
LengthLargestFullyConjugatedRingSystem::LengthLargestFullyConjugatedRingSystem()
    : Descriptor("LengthLargestFullyConjugatedRingSystem", "Number of atoms in the largest fully conjugated ring system") {}

std::variant<double, int, std::string> LengthLargestFullyConjugatedRingSystem::calculate(const ::desfact::Molecule& mol) const {
     if (!mol.isValid() || !mol.getMolecule()) return 0;
     const RDKit::ROMol& rdkMol = *mol.getMolecule();
     // Placeholder - Requires specific ring system analysis. Using placeholder.
     return placeholderIntCalculation(rdkMol, getName());
}

// 25. RatioConjugatedCarbonsTotalCarbons
RatioConjugatedCarbonsTotalCarbons::RatioConjugatedCarbonsTotalCarbons()
    : Descriptor("RatioConjugatedCarbonsTotalCarbons", "Ratio of sp2/sp carbon atoms to the total number of carbon atoms") {}

std::variant<double, int, std::string> RatioConjugatedCarbonsTotalCarbons::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    const RDKit::ROMol& rdkMol = *mol.getMolecule();
    int conjugatedCarbon = 0;
    int totalCarbon = 0;
    for(const auto& atom : rdkMol.atoms()){
        if(atom->getAtomicNum() == 6){
            totalCarbon++;
            RDKit::Atom::HybridizationType hyb = atom->getHybridization();
            if(hyb == RDKit::Atom::SP2 || hyb == RDKit::Atom::SP){
                conjugatedCarbon++;
            }
        }
    }
    if(totalCarbon == 0) return 0.0;
    return static_cast<double>(conjugatedCarbon) / totalCarbon;
}


// 26. NormalizedConjugationPathwayDensity
NormalizedConjugationPathwayDensity::NormalizedConjugationPathwayDensity()
    : Descriptor("NormalizedConjugationPathwayDensity", "Number of conjugated atoms divided by the total number of heavy atoms") {}

std::variant<double, int, std::string> NormalizedConjugationPathwayDensity::calculate(const ::desfact::Molecule& mol) const {
     if (!mol.isValid() || !mol.getMolecule()) return 0.0;
     const RDKit::ROMol& rdkMol = *mol.getMolecule();
     int heavyAtomCount = rdkMol.getNumHeavyAtoms();
     if(heavyAtomCount == 0) return 0.0;
     // Placeholder - Requires identifying all conjugated atoms. Using placeholder.
     int conjugatedAtoms = placeholderIntCalculation(rdkMol, "ConjugatedAtomCount"); // Hypothetical intermediate calc
     if (conjugatedAtoms < 0) conjugatedAtoms = 0; // Handle potential error from placeholder
     return static_cast<double>(conjugatedAtoms) / heavyAtomCount;
}


// --- Geometry Approximation / Exposure / Shape ---

// 27. EstimatedMolecularAspectRatio
EstimatedMolecularAspectRatio::EstimatedMolecularAspectRatio()
    : Descriptor("EstimatedMolecularAspectRatio", "Heuristic aspect ratio based on longest path vs. width estimate") {}

std::variant<double, int, std::string> EstimatedMolecularAspectRatio::calculate(const ::desfact::Molecule& mol) const {
     if (!mol.isValid() || !mol.getMolecule()) return 0.0;
     const RDKit::ROMol& rdkMol = *mol.getMolecule();
     // Placeholder - Requires complex path finding / shape estimation. Using placeholder.
     return placeholderCalculation(rdkMol, getName());
}

// 28. TerminalHeavyAtomDispersion
TerminalHeavyAtomDispersion::TerminalHeavyAtomDispersion()
    : Descriptor("TerminalHeavyAtomDispersion", "Measure of the spatial spread of terminal heavy atoms (heuristic)") {}

std::variant<double, int, std::string> TerminalHeavyAtomDispersion::calculate(const ::desfact::Molecule& mol) const {
     if (!mol.isValid() || !mol.getMolecule()) return 0.0;
     const RDKit::ROMol& rdkMol = *mol.getMolecule();
     // Placeholder - Requires distance calculations between terminal atoms. Using placeholder.
     return placeholderCalculation(rdkMol, getName());
}

// 29. ConnectivityCorePeripheryGradient
ConnectivityCorePeripheryGradient::ConnectivityCorePeripheryGradient()
    : Descriptor("ConnectivityCorePeripheryGradient", "Difference in average atom degree between core and peripheral atoms (heuristic)") {}

std::variant<double, int, std::string> ConnectivityCorePeripheryGradient::calculate(const ::desfact::Molecule& mol) const {
     if (!mol.isValid() || !mol.getMolecule()) return 0.0;
     const RDKit::ROMol& rdkMol = *mol.getMolecule();
     // Placeholder - Requires defining core/periphery. Using placeholder.
     return placeholderCalculation(rdkMol, getName());
}

// 30. HeavyAtomPlanarityHeuristic
HeavyAtomPlanarityHeuristic::HeavyAtomPlanarityHeuristic()
    : Descriptor("HeavyAtomPlanarityHeuristic", "Fraction of heavy atoms likely part of a planar system (sp2, aromatic)") {}

std::variant<double, int, std::string> HeavyAtomPlanarityHeuristic::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    const RDKit::ROMol& rdkMol = *mol.getMolecule();
    int heavyAtomCount = rdkMol.getNumHeavyAtoms();
    if(heavyAtomCount == 0) return 0.0;

    int planarHeavy = 0;
    for(const auto& atom : rdkMol.atoms()){
        if(atom->getAtomicNum() > 1){ // Heavy atom
             if(atom->getIsAromatic() || atom->getHybridization() == RDKit::Atom::SP2) {
                 planarHeavy++;
             }
        }
    }
    return static_cast<double>(planarHeavy) / heavyAtomCount;
}


// 31. RingCoreChainPeripherySizeRatio
RingCoreChainPeripherySizeRatio::RingCoreChainPeripherySizeRatio()
    : Descriptor("RingCoreChainPeripherySizeRatio", "Ratio of atoms in rings vs atoms in chains attached to rings (heuristic)") {}

std::variant<double, int, std::string> RingCoreChainPeripherySizeRatio::calculate(const ::desfact::Molecule& mol) const {
     if (!mol.isValid() || !mol.getMolecule()) return 0.0;
     const RDKit::ROMol& rdkMol = *mol.getMolecule();
     // Placeholder - Requires distinguishing ring/chain atoms. Using placeholder.
     return placeholderCalculation(rdkMol, getName());
}

// 32. RingChainAlternationFrequency
RingChainAlternationFrequency::RingChainAlternationFrequency()
    : Descriptor("RingChainAlternationFrequency", "Frequency of transitions between ring systems and chains along the molecule's 'spine' (heuristic)") {}

std::variant<double, int, std::string> RingChainAlternationFrequency::calculate(const ::desfact::Molecule& mol) const {
     if (!mol.isValid() || !mol.getMolecule()) return 0.0;
     const RDKit::ROMol& rdkMol = *mol.getMolecule();
     // Placeholder - Complex traversal needed. Using placeholder.
     return placeholderCalculation(rdkMol, getName());
}

// --- Charge and Ionization Site Heuristics ---

// 33. PotentialDeprotonationSiteCount
PotentialDeprotonationSiteCount::PotentialDeprotonationSiteCount()
    : Descriptor("PotentialDeprotonationSiteCount", "Count of potential acidic sites (O, N, S with attached H, or acidic C-H)") {}

std::variant<double, int, std::string> PotentialDeprotonationSiteCount::calculate(const ::desfact::Molecule& mol) const {
     if (!mol.isValid() || !mol.getMolecule()) return 0; // Return int error
     const RDKit::ROMol& rdkMol = *mol.getMolecule();
     try {
          // HBD should handle necessary calculations/checks
          // Remove sanitizeMol call
          unsigned int hbdCount = RDKit::Descriptors::calcNumHBD(rdkMol);
          return static_cast<int>(hbdCount); // Cast to int for the variant
     } catch (const std::exception& e) {
         globalLogger.debug("HBD Calc Error: " + std::string(e.what()) + " for " + mol.getOriginalSmiles());
          return std::string("Error: HBD Calc");
     } catch (...) {
          globalLogger.debug("Unknown HBD Calc Error for " + mol.getOriginalSmiles());
         return std::string("Error: HBD Calc Unknown");
     }
}

// 34. PotentialProtonationSiteCount
PotentialProtonationSiteCount::PotentialProtonationSiteCount()
    : Descriptor("PotentialProtonationSiteCount", "Count of potential basic sites (basic N mainly, possibly O)") {}

std::variant<double, int, std::string> PotentialProtonationSiteCount::calculate(const ::desfact::Molecule& mol) const {
     if (!mol.isValid() || !mol.getMolecule()) return 0; // Return int error
      const RDKit::ROMol& rdkMol = *mol.getMolecule();
     try {
         // HBA should handle necessary calculations/checks
         // Remove sanitizeMol call
          unsigned int hbaCount = RDKit::Descriptors::calcNumHBA(rdkMol);
          return static_cast<int>(hbaCount); // Cast to int for the variant
     } catch (const std::exception& e) {
         globalLogger.debug("HBA Calc Error: " + std::string(e.what()) + " for " + mol.getOriginalSmiles());
          return std::string("Error: HBA Calc");
     } catch (...) {
         globalLogger.debug("Unknown HBA Calc Error for " + mol.getOriginalSmiles());
          return std::string("Error: HBA Calc Unknown");
     }
}

// 35. FractionNonTerminalChargedSites
FractionNonTerminalChargedSites::FractionNonTerminalChargedSites()
    : Descriptor("FractionNonTerminalChargedSites", "Fraction of formally charged heavy atoms that are not in terminal positions") {}

std::variant<double, int, std::string> FractionNonTerminalChargedSites::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    const RDKit::ROMol& rdkMol = *mol.getMolecule();
    int nonTerminalCharged = 0;
    int totalCharged = 0;
    RDKit::RingInfo* ringInfo = rdkMol.getRingInfo(); // Get RingInfo once
     if (!ringInfo->isInitialized()) {
          RDKit::MolOps::findSSSR(const_cast<RDKit::ROMol&>(rdkMol));
     }

    for(const auto& atom : rdkMol.atoms()){
        if(atom->getAtomicNum() > 1 && atom->getFormalCharge() != 0){
            totalCharged++;
            // Check if non-terminal (degree > 1 or in a ring)
             if (atom->getDegree() > 1 || (ringInfo->numAtomRings(atom->getIdx()) > 0) ) {
                 nonTerminalCharged++;
            }
        }
    }

    if(totalCharged == 0) return 0.0;
    return static_cast<double>(nonTerminalCharged) / totalCharged;
}

// 36. MeanTopologicalDistanceIonizableSites
MeanTopologicalDistanceIonizableSites::MeanTopologicalDistanceIonizableSites()
    : Descriptor("MeanTopologicalDistanceIonizableSites", "Mean topological distance between potential acidic/basic sites") {}

std::variant<double, int, std::string> MeanTopologicalDistanceIonizableSites::calculate(const ::desfact::Molecule& mol) const {
     if (!mol.isValid() || !mol.getMolecule()) return 0.0;
     const RDKit::ROMol& rdkMol = *mol.getMolecule();
     // Placeholder - Requires identifying sites and distance matrix. Using placeholder.
     return placeholderCalculation(rdkMol, getName());
}

// 37. RingIonizableSiteDensity
RingIonizableSiteDensity::RingIonizableSiteDensity()
    : Descriptor("RingIonizableSiteDensity", "Number of potential ionizable sites within ring systems per ring atom") {}

std::variant<double, int, std::string> RingIonizableSiteDensity::calculate(const ::desfact::Molecule& mol) const {
     if (!mol.isValid() || !mol.getMolecule()) return 0.0;
     const RDKit::ROMol& rdkMol = *mol.getMolecule();
     // Placeholder - Requires identifying sites and ring atoms. Using placeholder.
     return placeholderCalculation(rdkMol, getName());
}

// 38. AromaticProtonDonorAcceptorRatioHeuristic
AromaticProtonDonorAcceptorRatioHeuristic::AromaticProtonDonorAcceptorRatioHeuristic()
    : Descriptor("AromaticProtonDonorAcceptorRatioHeuristic", "Ratio of H-bond donors to acceptors within aromatic systems") {}

std::variant<double, int, std::string> AromaticProtonDonorAcceptorRatioHeuristic::calculate(const ::desfact::Molecule& mol) const {
     if (!mol.isValid() || !mol.getMolecule()) return 0.0;
     const RDKit::ROMol& rdkMol = *mol.getMolecule();
     // Placeholder - Requires identifying aromatic donors/acceptors. Using placeholder.
     return placeholderCalculation(rdkMol, getName());
}


// --- Hydrophobicity / Polarity Features ---

// 39. HydrophilicAtomClusterSizeMean
HydrophilicAtomClusterSizeMean::HydrophilicAtomClusterSizeMean()
    : Descriptor("HydrophilicAtomClusterSizeMean", "Mean size of connected clusters of polar atoms (N, O, S, P, Halogens)") {}

std::variant<double, int, std::string> HydrophilicAtomClusterSizeMean::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    const RDKit::ROMol& rdkMol = *mol.getMolecule();

    std::vector<int> polar_indices;
    std::map<int, int> atom_idx_to_polar_idx;
    for (const auto& atom : rdkMol.atoms()) {
        if (isPolarAtom(atom)) {
             atom_idx_to_polar_idx[atom->getIdx()] = polar_indices.size();
            polar_indices.push_back(atom->getIdx());
        }
    }

    if (polar_indices.empty()) return 0.0;

    int n_polar = polar_indices.size();
    std::vector<std::vector<int>> adj(n_polar);
     for (int i = 0; i < n_polar; ++i) {
        const RDKit::Atom* atom1 = rdkMol.getAtomWithIdx(polar_indices[i]);
        // Use getAtomNeighbors iterator
        for (const auto& neighbor : rdkMol.atomNeighbors(atom1)) {
             if (isPolarAtom(neighbor)) {
                 // Check if neighbor is in our polar list
                 auto map_it = atom_idx_to_polar_idx.find(neighbor->getIdx());
                  if(map_it != atom_idx_to_polar_idx.end()){
                    int neighbor_polar_idx = map_it->second;
                    if (neighbor_polar_idx > i) {
                        adj[i].push_back(neighbor_polar_idx);
                        adj[neighbor_polar_idx].push_back(i);
                    }
                  }
             }
        }
    }

    double totalSizeSum = 0;
    int clusterCount = 0;
    std::vector<bool> visited(n_polar, false);
    for (int i = 0; i < n_polar; ++i) {
        if (!visited[i]) {
            int componentSize = 0;
            std::deque<int> q;
            q.push_back(i);
            visited[i] = true;
            componentSize++;
            clusterCount++;

             while(!q.empty()){
                int u = q.front();
                q.pop_front();
                for(int v : adj[u]){
                    if(!visited[v]){
                        visited[v] = true;
                        q.push_back(v);
                        componentSize++;
                    }
                }
            }
            totalSizeSum += componentSize;
        }
    }

    if (clusterCount == 0) return 0.0;
    return totalSizeSum / clusterCount;
}


// 40. HydrophobicIslandSizeMax
HydrophobicIslandSizeMax::HydrophobicIslandSizeMax()
    : Descriptor("HydrophobicIslandSizeMax", "Number of atoms in the largest connected fragment of hydrophobic atoms (C only bonded to C, H)") {}

std::variant<double, int, std::string> HydrophobicIslandSizeMax::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0;
    const RDKit::ROMol& rdkMol = *mol.getMolecule();

    std::vector<int> hydrophobic_indices;
    std::map<int, int> atom_idx_to_hydrophobic_idx;
    for (const auto& atom : rdkMol.atoms()) {
        if (isHydrophobicAtom(atom)) {
            atom_idx_to_hydrophobic_idx[atom->getIdx()] = hydrophobic_indices.size();
            hydrophobic_indices.push_back(atom->getIdx());
        }
    }

     if (hydrophobic_indices.empty()) return 0;

    int n_hydrophobic = hydrophobic_indices.size();
    std::vector<std::vector<int>> adj(n_hydrophobic);
    for (int i = 0; i < n_hydrophobic; ++i) {
        const RDKit::Atom* atom1 = rdkMol.getAtomWithIdx(hydrophobic_indices[i]);
        // Use getAtomNeighbors iterator
        for (const auto& neighbor : rdkMol.atomNeighbors(atom1)) {
             if (isHydrophobicAtom(neighbor)) { // Check if neighbor is also hydrophobic
                 // Check if neighbor is in our hydrophobic list
                 auto map_it = atom_idx_to_hydrophobic_idx.find(neighbor->getIdx());
                 if(map_it != atom_idx_to_hydrophobic_idx.end()){
                    int neighbor_hydrophobic_idx = map_it->second;
                    if (neighbor_hydrophobic_idx > i) {
                        adj[i].push_back(neighbor_hydrophobic_idx);
                        adj[neighbor_hydrophobic_idx].push_back(i);
                    }
                 }
             }
        }
    }

    int maxSize = 0;
    std::vector<bool> visited(n_hydrophobic, false);
    for (int i = 0; i < n_hydrophobic; ++i) {
        if (!visited[i]) {
            int componentSize = 0;
            std::deque<int> q;
            q.push_back(i);
            visited[i] = true;
            componentSize++;

            while(!q.empty()){
                int u = q.front();
                q.pop_front();
                for(int v : adj[u]){
                    if(!visited[v]){
                        visited[v] = true;
                        q.push_back(v);
                        componentSize++;
                    }
                }
            }
             maxSize = std::max(maxSize, componentSize);
        }
    }
    return maxSize;
}


// 41. HydrophobicPolarSurfaceInterfaceEstimation
HydrophobicPolarSurfaceInterfaceEstimation::HydrophobicPolarSurfaceInterfaceEstimation()
    : Descriptor("HydrophobicPolarSurfaceInterfaceEstimation", "Number of bonds between polar and non-polar heavy atoms") {}

std::variant<double, int, std::string> HydrophobicPolarSurfaceInterfaceEstimation::calculate(const ::desfact::Molecule& mol) const {
     if (!mol.isValid() || !mol.getMolecule()) return 0;
     const RDKit::ROMol& rdkMol = *mol.getMolecule();
     int interfaceBonds = 0;
     for (const auto& bond : rdkMol.bonds()) {
         const RDKit::Atom* atom1 = bond->getBeginAtom();
         const RDKit::Atom* atom2 = bond->getEndAtom();
         // Exclude bonds to hydrogen
         if (atom1->getAtomicNum() > 1 && atom2->getAtomicNum() > 1) {
             bool atom1_polar = isPolarAtom(atom1);
             bool atom2_polar = isPolarAtom(atom2);
             if (atom1_polar != atom2_polar) { // One is polar, the other is not
                 interfaceBonds++;
             }
         }
     }
     return interfaceBonds;
}

// 42. NormalizedHydrophobicPolarRatioEstimate
NormalizedHydrophobicPolarRatioEstimate::NormalizedHydrophobicPolarRatioEstimate()
    : Descriptor("NormalizedHydrophobicPolarRatioEstimate", "Ratio of hydrophobic heavy atoms to polar heavy atoms") {}

std::variant<double, int, std::string> NormalizedHydrophobicPolarRatioEstimate::calculate(const ::desfact::Molecule& mol) const {
     if (!mol.isValid() || !mol.getMolecule()) return 0.0;
     const RDKit::ROMol& rdkMol = *mol.getMolecule();
     int polarCount = 0;
     int hydrophobicCount = 0; // Using non-polar as proxy for hydrophobic here
     for(const auto& atom : rdkMol.atoms()){
         if(atom->getAtomicNum() > 1){ // Heavy atom
             if(isPolarAtom(atom)){
                 polarCount++;
             } else {
                 hydrophobicCount++;
             }
         }
     }
     if (polarCount == 0) {
          // Avoid division by zero. If there are hydrophobic atoms, return large number, else 0.
         return (hydrophobicCount > 0) ? 999.0 : 0.0;
     }
     return static_cast<double>(hydrophobicCount) / polarCount;
}

// 43. HydrophobicAtomPathLengthMean
HydrophobicAtomPathLengthMean::HydrophobicAtomPathLengthMean()
    : Descriptor("HydrophobicAtomPathLengthMean", "Mean topological distance between pairs of hydrophobic atoms") {}

std::variant<double, int, std::string> HydrophobicAtomPathLengthMean::calculate(const ::desfact::Molecule& mol) const {
     if (!mol.isValid() || !mol.getMolecule()) return 0.0;
     const RDKit::ROMol& rdkMol = *mol.getMolecule();
     // Placeholder - Requires identifying hydrophobic atoms and distance matrix. Using placeholder.
     return placeholderCalculation(rdkMol, getName());
}

// 44. HeteroatomConnectivityIndex
HeteroatomConnectivityIndex::HeteroatomConnectivityIndex()
    : Descriptor("HeteroatomConnectivityIndex", "Sum of connectivity values (1/sqrt(deg1*deg2)) for bonds involving at least one heteroatom") {}

std::variant<double, int, std::string> HeteroatomConnectivityIndex::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    const RDKit::ROMol& rdkMol = *mol.getMolecule();
    double index = 0.0;

    for(const auto& bond : rdkMol.bonds()){
        const RDKit::Atom* atom1 = bond->getBeginAtom();
        const RDKit::Atom* atom2 = bond->getEndAtom();
        int an1 = atom1->getAtomicNum();
        int an2 = atom2->getAtomicNum();

        if(an1 != 6 && an1 != 1 || an2 != 6 && an2 != 1) { // At least one heteroatom
             int deg1 = atom1->getDegree();
             int deg2 = atom2->getDegree();
             if (deg1 > 0 && deg2 > 0) {
                 index += 1.0 / std::sqrt(static_cast<double>(deg1 * deg2));
             }
        }
    }
    return index;
}


// 45. RingFusionDensity
RingFusionDensity::RingFusionDensity()
    : Descriptor("RingFusionDensity", "Number of ring atoms involved in >1 SSSR ring / total ring atoms") {}

std::variant<double, int, std::string> RingFusionDensity::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    const RDKit::ROMol& rdkMol = *mol.getMolecule();
    RDKit::RingInfo* ringInfo = rdkMol.getRingInfo();
    if (!ringInfo->isInitialized()) {
        RDKit::MolOps::findSSSR(const_cast<RDKit::ROMol&>(rdkMol));
    }

    int fusionAtoms = 0;
    std::set<int> allRingAtoms;

    for (const auto& ring : ringInfo->atomRings()) {
        for (int atomIdx : ring) {
            allRingAtoms.insert(atomIdx);
            if (ringInfo->numAtomRings(atomIdx) > 1) {
                // Use a set to count each fusion atom only once, even if in >2 rings
            }
        }
    }
     // Recalculate fusionAtoms using the set approach
     for (int atomIdx : allRingAtoms) {
         if (ringInfo->numAtomRings(atomIdx) > 1) {
             fusionAtoms++;
         }
     }


    if (allRingAtoms.empty()) return 0.0;

    return static_cast<double>(fusionAtoms) / allRingAtoms.size();
}


// 46. AverageBondPolarity
AverageBondPolarity::AverageBondPolarity()
    : Descriptor("AverageBondPolarity", "Average absolute difference in Pauling electronegativity across all bonds") {}

std::variant<double, int, std::string> AverageBondPolarity::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    const RDKit::ROMol& rdkMol = *mol.getMolecule();
    int numBonds = rdkMol.getNumBonds();
    if (numBonds == 0) return 0.0;

    double totalPolarity = 0.0;
    for(const auto& bond : rdkMol.bonds()){
        double en1 = getAtomENPauling(bond->getBeginAtom());
        double en2 = getAtomENPauling(bond->getEndAtom());
        totalPolarity += std::abs(en1 - en2);
    }

    return totalPolarity / numBonds;
}


// 47. StericHindranceProxy
StericHindranceProxy::StericHindranceProxy()
    : Descriptor("StericHindranceProxy", "Sum of branch point contributions (e.g., count of quaternary carbons)") {}

std::variant<double, int, std::string> StericHindranceProxy::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0;
    const RDKit::ROMol& rdkMol = *mol.getMolecule();
    int steric_proxy = 0;
     for (const auto& atom : rdkMol.atoms()) {
         if (atom->getAtomicNum() == 6 && atom->getDegree() == 4) { // Quaternary Carbon
             steric_proxy++;
         }
         // Could add other contributions here (e.g., degree 3 heavy atoms?)
     }
     return steric_proxy;
}

// 48. MolecularFlexibilityProxy
MolecularFlexibilityProxy::MolecularFlexibilityProxy()
    : Descriptor("MolecularFlexibilityProxy", "Number of rotatable bonds (RDKit implementation)") {}

std::variant<double, int, std::string> MolecularFlexibilityProxy::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0; // Return int error
     const RDKit::ROMol& rdkMol = *mol.getMolecule();
     try {
          // NumRotatableBonds should handle necessary calculations/checks
          // Remove sanitizeMol call
          unsigned int rotBonds = RDKit::Descriptors::calcNumRotatableBonds(rdkMol);
          return static_cast<int>(rotBonds); // Cast to int for the variant
     } catch (const std::exception& e) {
         globalLogger.debug("RotBond Calc Error: " + std::string(e.what()) + " for " + mol.getOriginalSmiles());
          return std::string("Error: RotBond Calc");
     } catch (...) {
          globalLogger.debug("Unknown RotBond Calc Error for " + mol.getOriginalSmiles());
          return std::string("Error: RotBond Calc Unknown");
     }
}


} // namespace descriptors
} // namespace desfact
