#include "descriptors/vague3.hpp"
#include <GraphMol/GraphMol.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Descriptors/MolDescriptors.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/AtomIterators.h>
#include <GraphMol/BondIterators.h>
#include <GraphMol/RingInfo.h>
#include <GraphMol/PeriodicTable.h>
#include <numeric>
#include <unordered_set>
#include <unordered_map>
#include <queue>
#include <cmath>
#include "utils.hpp"
#include <GraphMol/PartialCharges/GasteigerCharges.h>

namespace desfact {
namespace descriptors {

// Helper functions
namespace {
    // Calculate atom degree (number of connected atoms)
    inline unsigned int getAtomDegree(const RDKit::Atom* atom) {
        return atom->getDegree();
    }
    
    // Check if atom is a peripheral atom (degree ≤ 1)
    inline bool isPeripheralAtom(const RDKit::Atom* atom) {
        return getAtomDegree(atom) <= 1;
    }
    
    // Check if atom is a core atom (degree ≥ 3)
    inline bool isCoreAtom(const RDKit::Atom* atom) {
        return getAtomDegree(atom) >= 3;
    }
    
    // Check if atom is a carbon
    inline bool isCarbon(const RDKit::Atom* atom) {
        return atom->getAtomicNum() == 6;
    }
    
    // Check if atom is a heteroatom (not carbon or hydrogen)
    inline bool isHeteroatom(const RDKit::Atom* atom) {
        return atom->getAtomicNum() != 6 && atom->getAtomicNum() != 1;
    }
    
    // Check if atom is a terminal atom
    inline bool isTerminalAtom(const RDKit::Atom* atom) {
        return getAtomDegree(atom) == 1;
    }
    
    inline bool isAtomInRing(const RDKit::Atom* atom) {
        return atom->getIsAromatic() || atom->getOwningMol().getRingInfo()->numAtomRings(atom->getIdx()) > 0;
    }
    
    inline bool isAromaticAtom(const RDKit::Atom* atom) {
        return atom->getIsAromatic();
    }
    
    // Check if bond is aromatic
    inline bool isAromaticBond(const RDKit::Bond* bond) {
        return bond->getIsAromatic();
    }
    
    // Get electronegativity of an atom
    double getElectronegativity(const RDKit::Atom* atom) {
        static const std::unordered_map<int, double> electronegativityValues = {
            {1, 2.20}, {2, 0.00}, {3, 0.98}, {4, 1.57}, {5, 2.04}, {6, 2.55},
            {7, 3.04}, {8, 3.44}, {9, 3.98}, {10, 0.00}, {11, 0.93}, {12, 1.31},
            {13, 1.61}, {14, 1.90}, {15, 2.19}, {16, 2.58}, {17, 3.16}, {18, 0.00},
            {19, 0.82}, {20, 1.00}, {31, 1.81}, {32, 2.01}, {33, 2.18}, {34, 2.55},
            {35, 2.96}, {36, 3.00}, {37, 0.82}, {38, 0.95}, {49, 1.78}, {50, 1.96},
            {51, 2.05}, {52, 2.10}, {53, 2.66}, {54, 2.60}
        };
        int atomicNum = atom->getAtomicNum();
        auto it = electronegativityValues.find(atomicNum);
        return it != electronegativityValues.end() ? it->second : 0.0;
    }
    
    // Check if atom is highly electronegative (EN > 3.0)
    inline bool isHighlyElectronegative(const RDKit::Atom* atom) {
        return getElectronegativity(atom) > 3.0;
    }
    
    // Check if atom has high polarizability
    inline bool isHighlyPolarizable(const RDKit::Atom* atom) {
        // Atoms typically considered highly polarizable (halogens except F, S, P, etc.)
        std::vector<int> polarizableElements = {16, 15, 35, 53, 34, 52, 33, 51};
        return std::find(polarizableElements.begin(), polarizableElements.end(),
                       atom->getAtomicNum()) != polarizableElements.end();
    }
    
    // Check if atom is an H-bond donor
    inline bool isHBondDonor(const RDKit::Atom* atom) {
        int atomicNum = atom->getAtomicNum();
        if (atomicNum == 7 || atomicNum == 8 || atomicNum == 16) { // N, O, S
            for (const auto& nbr : atom->getOwningMol().atomNeighbors(atom)) {
                if (atom->getOwningMol().getAtomWithIdx(nbr->getIdx())->getAtomicNum() == 1) { // Has H attached
                    return true;
                }
            }
        }
        return false;
    }
    
    // Check if atom is an H-bond acceptor
    inline bool isHBondAcceptor(const RDKit::Atom* atom) {
        int atomicNum = atom->getAtomicNum();
        return (atomicNum == 7 || atomicNum == 8) && atom->getTotalValence() < 4; // N, O with lone pairs
    }
    
    // Check if atom is charged
    inline bool isCharged(const RDKit::Atom* atom) {
        return atom->getFormalCharge() != 0;
    }
    
    // Check if atom is a stereocenter
    inline bool isStereocenter(const RDKit::Atom* atom) {
        return atom->getChiralTag() != RDKit::Atom::ChiralType::CHI_UNSPECIFIED;
    }
    
    // Check if atom is alkali or alkaline earth
    inline bool isAlkaliOrAlkalineEarth(const RDKit::Atom* atom) {
        int atomicNum = atom->getAtomicNum();
        return (atomicNum == 3 || atomicNum == 11 || atomicNum == 19 || atomicNum == 37 || atomicNum == 55 || // Alkali
                atomicNum == 4 || atomicNum == 12 || atomicNum == 20 || atomicNum == 38 || atomicNum == 56);  // Alkaline earth
    }
    
    // Check if atom is from a rare group (group number ≥ 15, excluding common heteroatoms)
    inline bool isRareElement(const RDKit::Atom* atom) {
        int atomicNum = atom->getAtomicNum();
        if (atomicNum == 7 || atomicNum == 8 || atomicNum == 9 || atomicNum == 16 || atomicNum == 17) {
            return false; // Common heteroatoms (N, O, F, S, Cl)
        }
        
        // Simplified periodic table groups
        const std::unordered_map<int, int> elementGroups = {
            {5, 13}, {6, 14}, {7, 15}, {8, 16}, {9, 17}, {15, 15}, {16, 16}, {17, 17},
            {33, 15}, {34, 16}, {35, 17}, {51, 15}, {52, 16}, {53, 17}
        };
        
        auto it = elementGroups.find(atomicNum);
        int group = (it != elementGroups.end()) ? it->second : 0;
        
        return group >= 15;
    }
    
    // Check if bond is unsaturated (double or triple)
    inline bool isUnsaturatedBond(const RDKit::Bond* bond) {
        return bond->getBondType() == RDKit::Bond::BondType::DOUBLE ||
               bond->getBondType() == RDKit::Bond::BondType::TRIPLE;
    }
    
    // Get the shortest path between two atoms
    int getShortestPath(const RDKit::ROMol* mol, unsigned int atom1Idx, unsigned int atom2Idx) {
        // Manually compute distance using BFS
        std::vector<int> distances(mol->getNumAtoms(), -1);
        std::queue<unsigned int> queue;
        
        distances[atom1Idx] = 0;
        queue.push(atom1Idx);
        
        while (!queue.empty() && distances[atom2Idx] == -1) {
            unsigned int current = queue.front();
            queue.pop();
            
            for (const auto& nbr : mol->atomNeighbors(mol->getAtomWithIdx(current))) {
                unsigned int nbrIdx = nbr->getIdx();
                if (distances[nbrIdx] == -1) {
                    distances[nbrIdx] = distances[current] + 1;
                    queue.push(nbrIdx);
                }
            }
        }
        
        return distances[atom2Idx];
    }
    
    // Get the bond order as a numeric value
    double getBondOrder(const RDKit::Bond* bond) {
        switch (bond->getBondType()) {
            case RDKit::Bond::BondType::SINGLE: return 1.0;
            case RDKit::Bond::BondType::DOUBLE: return 2.0;
            case RDKit::Bond::BondType::TRIPLE: return 3.0;
            case RDKit::Bond::BondType::AROMATIC: return 1.5;
            default: return 1.0;
        }
    }
    
    // Calculate Shannon diversity index
    double calculateDiversityIndex(const std::vector<int>& counts) {
        double totalCount = std::accumulate(counts.begin(), counts.end(), 0);
        if (totalCount == 0) return 0.0;
        
        double diversity = 0.0;
        for (int count : counts) {
            if (count > 0) {
                double p = static_cast<double>(count) / totalCount;
                diversity -= p * std::log(p);
            }
        }
        return diversity;
    }
    
    // Calculate variance of a vector of values
    double calculateVariance(const std::vector<double>& values) {
        if (values.empty()) return 0.0;
        
        double sum = std::accumulate(values.begin(), values.end(), 0.0);
        double mean = sum / values.size();
        
        double squaredDiffSum = 0.0;
        for (double value : values) {
            double diff = value - mean;
            squaredDiffSum += diff * diff;
        }
        
        return squaredDiffSum / values.size();
    }
    
    // Find all carbon clusters (≥3 contiguous carbon atoms)
    std::vector<std::vector<unsigned int>> findCarbonClusters(const RDKit::ROMol* mol, int minSize = 3) {
        std::vector<std::vector<unsigned int>> clusters;
        std::vector<bool> visited(mol->getNumAtoms(), false);
        
        for (unsigned int i = 0; i < mol->getNumAtoms(); ++i) {
            if (visited[i] || mol->getAtomWithIdx(i)->getAtomicNum() != 6) continue;
            
            std::vector<unsigned int> cluster;
            std::queue<unsigned int> queue;
            queue.push(i);
            visited[i] = true;
            
            while (!queue.empty()) {
                unsigned int current = queue.front();
                queue.pop();
                cluster.push_back(current);
                
                for (const auto& nbr : mol->atomNeighbors(mol->getAtomWithIdx(current))) {
                    unsigned int nbrIdx = nbr->getIdx();
                    if (!visited[nbrIdx] && mol->getAtomWithIdx(nbrIdx)->getAtomicNum() == 6) {
                        queue.push(nbrIdx);
                        visited[nbrIdx] = true;
                    }
                }
            }
            
            if (cluster.size() >= static_cast<size_t>(minSize)) {
                clusters.push_back(cluster);
            }
        }
        
        return clusters;
    }
    
    // Find all heteroatom clusters (≥2 contiguous heteroatoms)
    std::vector<std::vector<unsigned int>> findHeteroatomClusters(const RDKit::ROMol* mol, int minSize = 2) {
        std::vector<std::vector<unsigned int>> clusters;
        std::vector<bool> visited(mol->getNumAtoms(), false);
        
        for (unsigned int i = 0; i < mol->getNumAtoms(); ++i) {
            if (visited[i] || !isHeteroatom(mol->getAtomWithIdx(i))) continue;
            
            std::vector<unsigned int> cluster;
            std::queue<unsigned int> queue;
            queue.push(i);
            visited[i] = true;
            
            while (!queue.empty()) {
                unsigned int current = queue.front();
                queue.pop();
                cluster.push_back(current);
                
                for (const auto& nbr : mol->atomNeighbors(mol->getAtomWithIdx(current))) {
                    unsigned int nbrIdx = nbr->getIdx();
                    if (!visited[nbrIdx] && isHeteroatom(mol->getAtomWithIdx(nbrIdx))) {
                        queue.push(nbrIdx);
                        visited[nbrIdx] = true;
                    }
                }
            }
            
            if (cluster.size() >= static_cast<size_t>(minSize)) {
                clusters.push_back(cluster);
            }
        }
        
        return clusters;
    }
    
    // Calculate closeness centrality for all atoms
    std::vector<double> calculateClosenessCentrality(const RDKit::ROMol* mol) {
        int nAtoms = mol->getNumAtoms();
        std::vector<double> centrality(nAtoms, 0.0);
        
        for (int i = 0; i < nAtoms; ++i) {
            double sumDistances = 0.0;
            int reachableAtoms = 0;
            
            for (int j = 0; j < nAtoms; ++j) {
                if (i == j) continue;
                
                int distance = getShortestPath(mol, i, j);
                if (distance > 0) {  // If atoms are connected
                    sumDistances += distance;
                    reachableAtoms++;
                }
            }
            
            if (reachableAtoms > 0) {
                centrality[i] = static_cast<double>(reachableAtoms) / sumDistances;
            }
        }
        
        return centrality;
    }
    
    // Get minimum spanning tree as adjacency list
    std::vector<std::vector<int>> getMinimumSpanningTree(const RDKit::ROMol* mol) {
        int nAtoms = mol->getNumAtoms();
        std::vector<std::vector<int>> adjacencyList(nAtoms);
        std::vector<bool> visited(nAtoms, false);
        std::vector<int> parent(nAtoms, -1);
        std::vector<double> key(nAtoms, std::numeric_limits<double>::max());
        
        key[0] = 0;
        
        for (int count = 0; count < nAtoms; ++count) {
            int u = -1;
            for (int v = 0; v < nAtoms; ++v) {
                if (!visited[v] && (u == -1 || key[v] < key[u])) {
                    u = v;
                }
            }
            
            if (u == -1) break; // Disconnected graph
            
            visited[u] = true;
            
            if (parent[u] != -1) {
                adjacencyList[u].push_back(parent[u]);
                adjacencyList[parent[u]].push_back(u);
            }
            
            for (const auto& nbr : mol->atomNeighbors(mol->getAtomWithIdx(u))) {
                int v = nbr->getIdx();
                double weight = 1.0; // All edges have equal weight for simplicity
                
                if (!visited[v] && weight < key[v]) {
                    parent[v] = u;
                    key[v] = weight;
                }
            }
        }
        
        return adjacencyList;
    }
    
    // Calculate the depth of a tree from its adjacency list
    int calculateTreeDepth(const std::vector<std::vector<int>>& adjacencyList) {
        if (adjacencyList.empty()) return 0;
        
        int nNodes = adjacencyList.size();
        std::vector<int> distances(nNodes, -1);
        std::queue<int> queue;
        
        // Start BFS from node 0
        distances[0] = 0;
        queue.push(0);
        
        while (!queue.empty()) {
            int u = queue.front();
            queue.pop();
            
            for (int v : adjacencyList[u]) {
                if (distances[v] == -1) {
                    distances[v] = distances[u] + 1;
                    queue.push(v);
                }
            }
        }
        
        // Find the maximum distance
        return *std::max_element(distances.begin(), distances.end());
    }
    
    // Identify chain fragments
    std::vector<std::vector<int>> findChainFragments(const RDKit::ROMol* mol) {
        std::vector<std::vector<int>> chainFragments;
        std::vector<bool> visited(mol->getNumAtoms(), false);
        
        // Mark all ring atoms as visited
        for (unsigned int i = 0; i < mol->getNumAtoms(); ++i) {
            if (isAtomInRing(mol->getAtomWithIdx(i))) {
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
                    if (!visited[nbrIdx] && !isAtomInRing(mol->getAtomWithIdx(nbrIdx))) {
                        queue.push(nbrIdx);
                        visited[nbrIdx] = true;
                    }
                }
            }
            
            chainFragments.push_back(fragment);
        }
        
        return chainFragments;
    }
    
    // Identify branches
    std::vector<std::vector<int>> findBranches(const RDKit::ROMol* mol) {
        std::vector<std::vector<int>> branches;
        std::vector<bool> visited(mol->getNumAtoms(), false);
        
        // Mark all core atoms (degree >= 3) as visited
        for (unsigned int i = 0; i < mol->getNumAtoms(); ++i) {
            if (isCoreAtom(mol->getAtomWithIdx(i))) {
                visited[i] = true;
            }
        }
        
        // Find branches starting from core atoms
        for (unsigned int i = 0; i < mol->getNumAtoms(); ++i) {
            if (!visited[i]) continue;
            
            const RDKit::Atom* atom = mol->getAtomWithIdx(i);
            
            for (const auto& nbr : mol->atomNeighbors(atom)) {
                int nbrIdx = nbr->getIdx();
                if (visited[nbrIdx]) continue;
                
                std::vector<int> branch;
                std::queue<int> queue;
                queue.push(nbrIdx);
                visited[nbrIdx] = true;
                
                while (!queue.empty()) {
                    int current = queue.front();
                    queue.pop();
                    branch.push_back(current);
                    
                    for (const auto& subNbr : mol->atomNeighbors(mol->getAtomWithIdx(current))) {
                        int subNbrIdx = subNbr->getIdx();
                        if (!visited[subNbrIdx]) {
                            queue.push(subNbrIdx);
                            visited[subNbrIdx] = true;
                        }
                    }
                }
                
                branches.push_back(branch);
            }
        }
        
        return branches;
    }
    
    // Is bond terminal (one of the atoms has degree 1)
    bool isTerminalBond(const RDKit::Bond* bond) {
        const RDKit::Atom* atom1 = bond->getBeginAtom();
        const RDKit::Atom* atom2 = bond->getEndAtom();
        return isTerminalAtom(atom1) || isTerminalAtom(atom2);
    }
    
    // Is bond part of a conjugated system
    bool isConjugatedBond(const RDKit::Bond* bond) {
        return bond->getIsConjugated();
    }
    
    // Count rings in molecule
    int countRings(const RDKit::ROMol* mol) {
        return mol->getRingInfo()->numRings();
    }
    
    // Count aromatic rings
    int countAromaticRings(const RDKit::ROMol* mol) {
        int aromaticRings = 0;
        const RDKit::RingInfo* ringInfo = mol->getRingInfo();
        
        for (unsigned int i = 0; i < ringInfo->numRings(); ++i) {
            const std::vector<int>& ringAtoms = ringInfo->atomRings()[i];
            bool isAromatic = true;
            
            for (int atomIdx : ringAtoms) {
                if (!mol->getAtomWithIdx(atomIdx)->getIsAromatic()) {
                    isAromatic = false;
                    break;
                }
            }
            
            if (isAromatic) {
                aromaticRings++;
            }
        }
        
        return aromaticRings;
    }
}

// Connectivity & Topological Shape Descriptors
AtomCentralityVarianceDescriptor::AtomCentralityVarianceDescriptor()
    : Descriptor("AtomCentralityVar", "Variance of atom closeness centrality in molecular graph") {}

std::variant<double, int, std::string> AtomCentralityVarianceDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) {
        return "Error: Invalid Molecule";
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    
    std::vector<double> centrality = calculateClosenessCentrality(rdkMol);
    double variance = calculateVariance(centrality);
    
    return variance;
}

PeripheryCoreRatioDescriptor::PeripheryCoreRatioDescriptor()
    : FractionalDescriptor("PeripheryCoreRatio", "Peripheral atoms (degree ≤1) / core atoms (degree ≥3)") {}

std::variant<double, int, std::string> PeripheryCoreRatioDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) {
        return "Error: Invalid Molecule";
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    
    int peripheralCount = 0;
    int coreCount = 0;
    
    for (const RDKit::Atom* atom : rdkMol->atoms()) {
        if (isPeripheralAtom(atom)) {
            peripheralCount++;
        } else if (isCoreAtom(atom)) {
            coreCount++;
        }
    }
    
    if (coreCount == 0) {
        return 0.0; // No core atoms
    }
    
    return static_cast<double>(peripheralCount) / coreCount;
}

MinimumSpanningTreeDepthDescriptor::MinimumSpanningTreeDepthDescriptor()
    : Descriptor("MSTreeDepth", "Depth of the molecular minimum spanning tree") {}

std::variant<double, int, std::string> MinimumSpanningTreeDepthDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) {
        return "Error: Invalid Molecule";
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    
    std::vector<std::vector<int>> mst = getMinimumSpanningTree(rdkMol);
    int depth = calculateTreeDepth(mst);
    
    return depth;
}

TerminalBranchRatioDescriptor::TerminalBranchRatioDescriptor()
    : FractionalDescriptor("TerminalBranchRatio", "Terminal branches / total branches") {}

std::variant<double, int, std::string> TerminalBranchRatioDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) {
        return "Error: Invalid Molecule";
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    
    std::vector<std::vector<int>> branches = findBranches(rdkMol);
    
    if (branches.empty()) {
        return 0.0; // No branches
    }
    
    int terminalBranches = 0;
    for (const auto& branch : branches) {
        bool isTerminal = false;
        for (int atomIdx : branch) {
            if (isTerminalAtom(rdkMol->getAtomWithIdx(atomIdx))) {
                isTerminal = true;
                break;
            }
        }
        if (isTerminal) {
            terminalBranches++;
        }
    }
    
    return static_cast<double>(terminalBranches) / branches.size();
}

AverageAtomPathRedundancyDescriptor::AverageAtomPathRedundancyDescriptor()
    : SumDescriptor("AvgAtomPathRedundancy", "Mean number of alternative shortest paths between atom pairs") {}

std::variant<double, int, std::string> AverageAtomPathRedundancyDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) {
        return "Error: Invalid Molecule";
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    int nAtoms = rdkMol->getNumAtoms();
    
    if (nAtoms <= 2) {
        return 0.0; // No path redundancy possible
    }
    
    // This is a simplified approximation using ring membership
    int ringCount = countRings(rdkMol);
    double avgRedundancy = static_cast<double>(ringCount) / nAtoms;
    
    return avgRedundancy;
}

// Atom-type Distribution Descriptors
CarbonClusterDensityDescriptor::CarbonClusterDensityDescriptor()
    : FractionalDescriptor("CarbonClusterDensity", "Clusters of ≥3 contiguous carbon atoms / total carbons") {}

std::variant<double, int, std::string> CarbonClusterDensityDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) {
        return "Error: Invalid Molecule";
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    
    int totalCarbons = 0;
    for (const RDKit::Atom* atom : rdkMol->atoms()) {
        if (isCarbon(atom)) {
            totalCarbons++;
        }
    }
    
    if (totalCarbons == 0) {
        return 0.0; // No carbon atoms
    }
    
    std::vector<std::vector<unsigned int>> carbonClusters = findCarbonClusters(rdkMol, 3);
    int clusteredCarbons = 0;
    
    for (const auto& cluster : carbonClusters) {
        clusteredCarbons += cluster.size();
    }
    
    return static_cast<double>(clusteredCarbons) / totalCarbons;
}

HeteroatomClusterDensityDescriptor::HeteroatomClusterDensityDescriptor()
    : FractionalDescriptor("HeteroatomClusterDensity", "Clusters of ≥2 contiguous heteroatoms / total heteroatoms") {}

std::variant<double, int, std::string> HeteroatomClusterDensityDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) {
        return "Error: Invalid Molecule";
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    
    int totalHeteroatoms = 0;
    for (const RDKit::Atom* atom : rdkMol->atoms()) {
        if (isHeteroatom(atom)) {
            totalHeteroatoms++;
        }
    }
    
    if (totalHeteroatoms == 0) {
        return 0.0; // No heteroatoms
    }
    
    std::vector<std::vector<unsigned int>> heteroatomClusters = findHeteroatomClusters(rdkMol, 2);
    int clusteredHeteroatoms = 0;
    
    for (const auto& cluster : heteroatomClusters) {
        clusteredHeteroatoms += cluster.size();
    }
    
    return static_cast<double>(clusteredHeteroatoms) / totalHeteroatoms;
}

TerminalHeavyAtomDiversityDescriptor::TerminalHeavyAtomDiversityDescriptor()
    : Descriptor("TerminalHeavyAtomDiversity", "Diversity of terminal atom types") {}

std::variant<double, int, std::string> TerminalHeavyAtomDiversityDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) {
        return "Error: Invalid Molecule";
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    
    std::unordered_map<int, int> elementCounts;
    
    for (const RDKit::Atom* atom : rdkMol->atoms()) {
        if (isTerminalAtom(atom) && atom->getAtomicNum() > 1) { // Heavy terminal atoms
            elementCounts[atom->getAtomicNum()]++;
        }
    }
    
    std::vector<int> counts;
    for (const auto& pair : elementCounts) {
        counts.push_back(pair.second);
    }
    
    double diversity = calculateDiversityIndex(counts);
    
    return diversity;
}

CentralHeteroatomRatioDescriptor::CentralHeteroatomRatioDescriptor()
    : FractionalDescriptor("CentralHeteroatomRatio", "Central heteroatoms (degree ≥3) / total heteroatoms") {}

std::variant<double, int, std::string> CentralHeteroatomRatioDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) {
        return "Error: Invalid Molecule";
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    
    int centralHeteroatoms = 0;
    int totalHeteroatoms = 0;
    
    for (const RDKit::Atom* atom : rdkMol->atoms()) {
        if (isHeteroatom(atom)) {
            totalHeteroatoms++;
            
            if (isCoreAtom(atom)) {
                centralHeteroatoms++;
            }
        }
    }
    
    if (totalHeteroatoms == 0) {
        return 0.0; // No heteroatoms
    }
    
    return static_cast<double>(centralHeteroatoms) / totalHeteroatoms;
}

ChainEndElementDiversityDescriptor::ChainEndElementDiversityDescriptor()
    : Descriptor("ChainEndElementDiversity", "Unique element types at chain ends") {}

std::variant<double, int, std::string> ChainEndElementDiversityDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) {
        return "Error: Invalid Molecule";
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    
    std::vector<std::vector<int>> chainFragments = findChainFragments(rdkMol);
    std::unordered_set<int> uniqueElements;
    
    for (const auto& fragment : chainFragments) {
        if (fragment.empty()) continue;
        
        // Find chain ends (atoms with only one connection within the chain)
        for (int atomIdx : fragment) {
            const RDKit::Atom* atom = rdkMol->getAtomWithIdx(atomIdx);
            int connections = 0;
            
            for (const auto& nbr : rdkMol->atomNeighbors(atom)) {
                int nbrIdx = nbr->getIdx();
                if (std::find(fragment.begin(), fragment.end(), nbrIdx) != fragment.end()) {
                    connections++;
                }
            }
            
            if (connections <= 1) {
                uniqueElements.insert(atom->getAtomicNum());
            }
        }
    }
    
    return static_cast<double>(uniqueElements.size());
}

// Electronic and Polar Characteristics Descriptors
HighENNeighborDensityDescriptor::HighENNeighborDensityDescriptor()
    : FractionalDescriptor("HighENNeighborDensity", "Atoms bonded to atoms with EN > 3.0 / total atoms") {}

std::variant<double, int, std::string> HighENNeighborDensityDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) {
        return "Error: Invalid Molecule";
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    
    int atomsWithHighENNeighbors = 0;
    int totalAtoms = rdkMol->getNumAtoms();
    
    for (const RDKit::Atom* atom : rdkMol->atoms()) {
        if (atom->getAtomicNum() == 1) continue; // Skip hydrogens
        
        bool hasHighENNeighbor = false;
        for (const auto& nbr : rdkMol->atomNeighbors(atom)) {
            if (isHighlyElectronegative(nbr)) {
                hasHighENNeighbor = true;
                break;
            }
        }
        
        if (hasHighENNeighbor) {
            atomsWithHighENNeighbors++;
        }
    }
    
    return static_cast<double>(atomsWithHighENNeighbors) / totalAtoms;
}

ElectronDeficientCarbonRatioDescriptor::ElectronDeficientCarbonRatioDescriptor()
    : FractionalDescriptor("ElectronDeficientCRatio", "Carbons bonded directly to electronegative substituents / total carbons") {}

std::variant<double, int, std::string> ElectronDeficientCarbonRatioDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) {
        return "Error: Invalid Molecule";
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    
    int electronDeficientCarbons = 0;
    int totalCarbons = 0;
    
    for (const RDKit::Atom* atom : rdkMol->atoms()) {
        if (isCarbon(atom)) {
            totalCarbons++;
            
            bool hasENSubstituent = false;
            for (const auto& nbr : rdkMol->atomNeighbors(atom)) {
                if (isHighlyElectronegative(nbr)) {
                    hasENSubstituent = true;
                    break;
                }
            }
            
            if (hasENSubstituent) {
                electronDeficientCarbons++;
            }
        }
    }
    
    if (totalCarbons == 0) {
        return 0.0;
    }
    
    return static_cast<double>(electronDeficientCarbons) / totalCarbons;
}

PolarizableAtomDispersionDescriptor::PolarizableAtomDispersionDescriptor()
    : Descriptor("PolarizableAtomDispersion", "Variance of distances between atoms with high polarizability") {}

std::variant<double, int, std::string> PolarizableAtomDispersionDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) {
        return "Error: Invalid Molecule";
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    
    std::vector<unsigned int> polarizableAtoms;
    for (unsigned int i = 0; i < rdkMol->getNumAtoms(); ++i) {
        if (isHighlyPolarizable(rdkMol->getAtomWithIdx(i))) {
            polarizableAtoms.push_back(i);
        }
    }
    
    if (polarizableAtoms.size() <= 1) {
        return 0.0; // Not enough polarizable atoms
    }
    
    std::vector<double> distances;
    for (size_t i = 0; i < polarizableAtoms.size(); ++i) {
        for (size_t j = i + 1; j < polarizableAtoms.size(); ++j) {
            int dist = getShortestPath(rdkMol, polarizableAtoms[i], polarizableAtoms[j]);
            distances.push_back(static_cast<double>(dist));
        }
    }
    
    double variance = calculateVariance(distances);
    
    return variance;
}

DipoleMomentProxyDescriptor::DipoleMomentProxyDescriptor()
    : FractionalDescriptor("DipoleMomentProxy", "Fraction of bonds between atoms with EN difference > 1.0") {}

std::variant<double, int, std::string> DipoleMomentProxyDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) {
        return 0.0; // Changed NaN to 0.0
    }
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return 0.0; // Changed NaN to 0.0

    int polarBonds = 0;
    int totalBonds = rdkMol->getNumBonds();

    // Explicitly handle the case of zero bonds to avoid potential division issues
    if (totalBonds == 0) {
        return 0.0;
    }

    for (const RDKit::Bond* bond : rdkMol->bonds()) {
        const RDKit::Atom* atom1 = bond->getBeginAtom();
        const RDKit::Atom* atom2 = bond->getEndAtom();

        double en1 = getElectronegativity(atom1);
        double en2 = getElectronegativity(atom2);
        double enDiff = std::abs(en1 - en2);

        if (enDiff > 1.0) {
            polarBonds++;
        }
    }

    // totalBonds is guaranteed > 0 here
    return static_cast<double>(polarBonds) / totalBonds;
}

// Structural Patterns Descriptors
LongChainFragmentDensityDescriptor::LongChainFragmentDensityDescriptor()
    : FractionalDescriptor("LongChainFragmentDensity", "Number of chain fragments ≥5 atoms / total chain fragments") {}

std::variant<double, int, std::string> LongChainFragmentDensityDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) {
        return "Error: Invalid Molecule";
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    
    std::vector<std::vector<int>> chainFragments = findChainFragments(rdkMol);
    
    if (chainFragments.empty()) {
        return 0.0; // No chain fragments
    }
    
    int longChainFragments = 0;
    
    for (const auto& fragment : chainFragments) {
        if (fragment.size() >= 5) {
            longChainFragments++;
        }
    }
    
    return static_cast<double>(longChainFragments) / chainFragments.size();
}

ShortChainDensityDescriptor::ShortChainDensityDescriptor()
    : FractionalDescriptor("ShortChainDensity", "Chain fragments ≤3 atoms / total chain fragments") {}

std::variant<double, int, std::string> ShortChainDensityDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) {
        return "Error: Invalid Molecule";
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    
    std::vector<std::vector<int>> chainFragments = findChainFragments(rdkMol);
    
    if (chainFragments.empty()) {
        return 0.0; // No chain fragments
    }
    
    int shortChainFragments = 0;
    
    for (const auto& fragment : chainFragments) {
        if (fragment.size() <= 3) {
            shortChainFragments++;
        }
    }
    
    return static_cast<double>(shortChainFragments) / chainFragments.size();
}

SubstitutionDensityPerRingDescriptor::SubstitutionDensityPerRingDescriptor()
    : Descriptor("SubstitutionDensityPerRing", "Mean substituents per ring atom") {}

std::variant<double, int, std::string> SubstitutionDensityPerRingDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) {
        return "Error: Invalid Molecule";
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    
    int ringAtoms = 0;
    int substituents = 0;
    
    for (const RDKit::Atom* atom : rdkMol->atoms()) {
        if (isAtomInRing(atom)) {
            ringAtoms++;
            
            for (const auto& nbr : rdkMol->atomNeighbors(atom)) {
                if (!isAtomInRing(nbr)) {
                    substituents++;
                }
            }
        }
    }
    
    if (ringAtoms == 0) {
        return 0.0; // No ring atoms
    }
    
    return static_cast<double>(substituents) / ringAtoms;
}

LinearVsBranchedCarbonRatioDescriptor::LinearVsBranchedCarbonRatioDescriptor()
    : Descriptor("LinearVsBranchedCRatio", "Linear-chain carbons / branched carbons") {}

std::variant<double, int, std::string> LinearVsBranchedCarbonRatioDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) {
        return 0.0; // Changed NaN to 0.0
    }
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return 0.0; // Changed NaN to 0.0

    int linearCarbons = 0;
    int branchedCarbons = 0;
    
    for (const RDKit::Atom* atom : rdkMol->atoms()) {
        if (isCarbon(atom) && !isAtomInRing(atom)) {
            unsigned int degree = getAtomDegree(atom);
            
            if (degree <= 2) {
                linearCarbons++;
            } else {
                branchedCarbons++;
            }
        }
    }
    
    if (branchedCarbons == 0) {
        return linearCarbons > 0 ? 9999.9 : 0.0; // Keep large num/0.0
    }
    
    return static_cast<double>(linearCarbons) / branchedCarbons;
}

RingToBranchConnectivityRatioDescriptor::RingToBranchConnectivityRatioDescriptor()
    : FractionalDescriptor("RingToBranchConnectivityRatio", "Ring atoms connected directly to branches / total ring atoms") {}

std::variant<double, int, std::string> RingToBranchConnectivityRatioDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) {
        return "Error: Invalid Molecule";
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    
    int ringAtoms = 0;
    int ringAtomsToBranches = 0;
    
    for (const RDKit::Atom* atom : rdkMol->atoms()) {
        if (isAtomInRing(atom)) {
            ringAtoms++;
            
            bool connectedToBranch = false;
            for (const auto& nbr : rdkMol->atomNeighbors(atom)) {
                if (!isAtomInRing(nbr)) {
                    connectedToBranch = true;
                    break;
                }
            }
            
            if (connectedToBranch) {
                ringAtomsToBranches++;
            }
        }
    }
    
    if (ringAtoms == 0) {
        return 0.0; // No ring atoms
    }
    
    return static_cast<double>(ringAtomsToBranches) / ringAtoms;
}

// Hydrogen Bond Networks
IsolatedHBondSiteDensityDescriptor::IsolatedHBondSiteDensityDescriptor()
    : FractionalDescriptor("IsolatedHBondSiteDensity", "H-bond capable atoms with no neighboring H-bond capable atoms") {}

std::variant<double, int, std::string> IsolatedHBondSiteDensityDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) {
        return "Error: Invalid Molecule";
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    
    int hbondCapableAtoms = 0;
    int isolatedHbondSites = 0;
    
    for (unsigned int i = 0; i < rdkMol->getNumAtoms(); ++i) {
        const RDKit::Atom* atom = rdkMol->getAtomWithIdx(i);
        
        if (isHBondDonor(atom) || isHBondAcceptor(atom)) {
            hbondCapableAtoms++;
            
            bool isIsolated = true;
            for (const auto& nbr : rdkMol->atomNeighbors(atom)) {
                if (isHBondDonor(nbr) || isHBondAcceptor(nbr)) {
                    isIsolated = false;
                    break;
                }
            }
            
            if (isIsolated) {
                isolatedHbondSites++;
            }
        }
    }
    
    if (hbondCapableAtoms == 0) {
        return 0.0; // No H-bond capable atoms
    }
    
    return static_cast<double>(isolatedHbondSites) / hbondCapableAtoms;
}

// Hydrogen Bond Networks (continued)
FunctionalGroupSeparationIndexDescriptor::FunctionalGroupSeparationIndexDescriptor()
    : Descriptor("FunctionalGroupSeparationIndex", "Mean shortest path between different functional group atoms") {}

std::variant<double, int, std::string> FunctionalGroupSeparationIndexDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) {
        return "Error: Invalid Molecule";
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    
    // Define functional groups using SMARTS patterns
    std::vector<std::string> functionalGroupSmarts = {
        "[OH]", // Alcohol/Phenol
        "[NH,NH2]", // Amine
        "[CX3](=O)[OX2H,OX1-]", // Carboxylic acid
        "[CX3](=O)[OX2]", // Ester
        "[CX3](=O)[NX3]", // Amide
        "[CX3](=O)", // Carbonyl
        "[NX3]", // Amine/Amide
        "[SX2]", // Thiol/Thioether
    };
    
    std::vector<std::vector<int>> functionalGroupAtoms;
    
    // Find atoms belonging to each functional group
    for (const std::string& smarts : functionalGroupSmarts) {
        try {
            RDKit::RWMol* query = RDKit::SmartsToMol(smarts);
            if (!query) continue;
            
            std::vector<std::vector<std::pair<int, int>>> matches;
            RDKit::SubstructMatch(*rdkMol, *query, matches);
            
            if (!matches.empty()) {
                std::vector<int> groupAtoms;
                for (const auto& match : matches) {
                    for (const auto& pair : match) {
                        groupAtoms.push_back(pair.second);
                    }
                }
                if (!groupAtoms.empty()) {
                    functionalGroupAtoms.push_back(groupAtoms);
                }
            }
            
            delete query;
        } catch (...) {
            // Skip if SMARTS parsing fails
            continue;
        }
    }
    
    if (functionalGroupAtoms.size() <= 1) {
        return 0.0; // No or only one functional group found
    }
    
    double totalDistance = 0.0;
    int pairCount = 0;
    
    // Calculate mean distance between all pairs of different functional groups
    for (size_t i = 0; i < functionalGroupAtoms.size(); ++i) {
        for (size_t j = i + 1; j < functionalGroupAtoms.size(); ++j) {
            double minDistance = std::numeric_limits<double>::max();
            
            // Find minimum distance between any atoms of the two functional groups
            for (int atomI : functionalGroupAtoms[i]) {
                for (int atomJ : functionalGroupAtoms[j]) {
                    int distance = getShortestPath(rdkMol, atomI, atomJ);
                    if (distance > 0 && distance < minDistance) {
                        minDistance = distance;
                    }
                }
            }
            
            if (minDistance < std::numeric_limits<double>::max()) {
                totalDistance += minDistance;
                pairCount++;
            }
        }
    }
    
    if (pairCount == 0) {
        return 0.0;
    }
    
    return totalDistance / pairCount;
}

PeripheralHBondDensityDescriptor::PeripheralHBondDensityDescriptor()
    : FractionalDescriptor("PeripheralHBondDensity", "Peripheral atoms capable of hydrogen bonding / total peripheral atoms") {}

std::variant<double, int, std::string> PeripheralHBondDensityDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) {
        return "Error: Invalid Molecule";
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    
    int peripheralAtoms = 0;
    int peripheralHBondAtoms = 0;
    
    for (const RDKit::Atom* atom : rdkMol->atoms()) {
        if (isPeripheralAtom(atom) && atom->getAtomicNum() > 1) { // Only count heavy atoms
            peripheralAtoms++;
            
            if (isHBondDonor(atom) || isHBondAcceptor(atom)) {
                peripheralHBondAtoms++;
            }
        }
    }
    
    if (peripheralAtoms == 0) {
        return 0.0;
    }
    
    return static_cast<double>(peripheralHBondAtoms) / peripheralAtoms;
}

// Aromatic and Conjugation Properties
AromaticCoreRatioDescriptor::AromaticCoreRatioDescriptor()
    : FractionalDescriptor("AromaticCoreRatio", "Aromatic atoms in ring cores / total aromatic atoms") {}

std::variant<double, int, std::string> AromaticCoreRatioDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) {
        return "Error: Invalid Molecule";
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    
    int aromaticAtoms = 0;
    int aromaticCoreAtoms = 0;
    
    for (const RDKit::Atom* atom : rdkMol->atoms()) {
        if (isAromaticAtom(atom)) {
            aromaticAtoms++;
            
            if (getAtomDegree(atom) >= 3) { // Core atoms have higher degree
                aromaticCoreAtoms++;
            }
        }
    }
    
    if (aromaticAtoms == 0) {
        return 0.0;
    }
    
    return static_cast<double>(aromaticCoreAtoms) / aromaticAtoms;
}

ConjugationGapRatioDescriptor::ConjugationGapRatioDescriptor()
    : FractionalDescriptor("ConjugationGapRatio", "Number of conjugation breaks / total conjugated fragments") {}

std::variant<double, int, std::string> ConjugationGapRatioDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) {
        return "Error: Invalid Molecule";
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    
    // Find conjugated fragments
    std::vector<bool> visited(rdkMol->getNumAtoms(), false);
    std::vector<std::vector<int>> conjugatedFragments;
    
    for (unsigned int i = 0; i < rdkMol->getNumAtoms(); ++i) {
        if (visited[i]) continue;
        
        const RDKit::Atom* atom = rdkMol->getAtomWithIdx(i);
        
        // Skip atoms that can't be conjugated
        if (atom->getIsAromatic() || 
            atom->getHybridization() == RDKit::Atom::HybridizationType::SP2 ||
            atom->getHybridization() == RDKit::Atom::HybridizationType::SP) {
            
            std::vector<int> fragment;
            std::queue<int> queue;
            queue.push(i);
            visited[i] = true;
            
            while (!queue.empty()) {
                int current = queue.front();
                queue.pop();
                fragment.push_back(current);
                
                const RDKit::Atom* currentAtom = rdkMol->getAtomWithIdx(current);
                
                for (const auto& nbr : rdkMol->atomNeighbors(currentAtom)) {
                    int nbrIdx = nbr->getIdx();
                    if (visited[nbrIdx]) continue;
                    
                    // Check if the bond is conjugated
                    const RDKit::Bond* bond = rdkMol->getBondBetweenAtoms(current, nbrIdx);
                    if (bond && (bond->getIsAromatic() || bond->getIsConjugated() ||
                                bond->getBondType() == RDKit::Bond::BondType::DOUBLE ||
                                bond->getBondType() == RDKit::Bond::BondType::TRIPLE)) {
                        queue.push(nbrIdx);
                        visited[nbrIdx] = true;
                    }
                }
            }
            
            if (!fragment.empty()) {
                conjugatedFragments.push_back(fragment);
            }
        }
    }
    
    if (conjugatedFragments.empty()) {
        return 0.0; // No conjugated fragments
    }
    
    // Count conjugation breaks (number of fragments - 1)
    // assuming a single fully conjugated molecule would have no breaks
    return static_cast<double>(conjugatedFragments.size() - 1) / conjugatedFragments.size();
}

NonAromaticRingSubstitutionDescriptor::NonAromaticRingSubstitutionDescriptor()
    : FractionalDescriptor("NonAromaticRingSubstitution", "Substituents on non-aromatic rings / total ring substituents") {}

std::variant<double, int, std::string> NonAromaticRingSubstitutionDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) {
        return "Error: Invalid Molecule";
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    
    int totalRingSubstituents = 0;
    int nonAromaticRingSubstituents = 0;
    
    for (const RDKit::Atom* atom : rdkMol->atoms()) {
        if (isAtomInRing(atom)) {
            bool isAromatic = isAromaticAtom(atom);
            
            for (const auto& nbr : rdkMol->atomNeighbors(atom)) {
                if (!isAtomInRing(nbr)) {
                    totalRingSubstituents++;
                    
                    if (!isAromatic) {
                        nonAromaticRingSubstituents++;
                    }
                }
            }
        }
    }
    
    if (totalRingSubstituents == 0) {
        return 0.0;
    }
    
    return static_cast<double>(nonAromaticRingSubstituents) / totalRingSubstituents;
}

AromaticToNonAromaticRingRatioDescriptor::AromaticToNonAromaticRingRatioDescriptor()
    : Descriptor("AromaticToNonAromaticRingRatio", "Aromatic rings / total rings") {}

std::variant<double, int, std::string> AromaticToNonAromaticRingRatioDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) {
        return "Error: Invalid Molecule";
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    
    int totalRings = countRings(rdkMol);
    if (totalRings == 0) {
        return 0.0;
    }
    
    int aromaticRings = countAromaticRings(rdkMol);
    
    return static_cast<double>(aromaticRings) / totalRings;
}

ConjugationTerminalDensityDescriptor::ConjugationTerminalDensityDescriptor()
    : FractionalDescriptor("ConjugationTerminalDensity", "Terminal atoms involved in conjugation / total conjugated atoms") {}

std::variant<double, int, std::string> ConjugationTerminalDensityDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) {
        return "Error: Invalid Molecule";
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    
    int conjugatedAtoms = 0;
    int terminalConjugatedAtoms = 0;
    
    for (const RDKit::Atom* atom : rdkMol->atoms()) {
        bool isConjugated = atom->getIsAromatic() || 
                           atom->getHybridization() == RDKit::Atom::HybridizationType::SP2 ||
                           atom->getHybridization() == RDKit::Atom::HybridizationType::SP;
        
        if (isConjugated) {
            conjugatedAtoms++;
            
            if (isTerminalAtom(atom)) {
                terminalConjugatedAtoms++;
            }
        }
    }
    
    if (conjugatedAtoms == 0) {
        return 0.0;
    }
    
    return static_cast<double>(terminalConjugatedAtoms) / conjugatedAtoms;
}

// Charge Distribution Characteristics
ChargePairDensityDescriptor::ChargePairDensityDescriptor()
    : Descriptor("ChargePairDensity", "Count of pairs of oppositely charged atoms separated by ≤4 bonds") {}

std::variant<double, int, std::string> ChargePairDensityDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) {
        return "Error: Invalid Molecule";
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    
    std::vector<int> positiveAtoms;
    std::vector<int> negativeAtoms;
    
    for (unsigned int i = 0; i < rdkMol->getNumAtoms(); ++i) {
        int charge = rdkMol->getAtomWithIdx(i)->getFormalCharge();
        if (charge > 0) {
            positiveAtoms.push_back(i);
        } else if (charge < 0) {
            negativeAtoms.push_back(i);
        }
    }
    
    int chargePairs = 0;
    
    for (int posIdx : positiveAtoms) {
        for (int negIdx : negativeAtoms) {
            int distance = getShortestPath(rdkMol, posIdx, negIdx);
            if (distance > 0 && distance <= 4) {
                chargePairs++;
            }
        }
    }
    
    return chargePairs;
}

FormalChargeAccessibilityDescriptor::FormalChargeAccessibilityDescriptor()
    : FractionalDescriptor("FormalChargeAccessibility", "Charged atoms at molecular periphery / total charged atoms") {}

std::variant<double, int, std::string> FormalChargeAccessibilityDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) {
        return "Error: Invalid Molecule";
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    
    int chargedAtoms = 0;
    int peripheralChargedAtoms = 0;
    
    for (const RDKit::Atom* atom : rdkMol->atoms()) {
        if (atom->getFormalCharge() != 0) {
            chargedAtoms++;
            
            if (isPeripheralAtom(atom)) {
                peripheralChargedAtoms++;
            }
        }
    }
    
    if (chargedAtoms == 0) {
        return 0.0;
    }
    
    return static_cast<double>(peripheralChargedAtoms) / chargedAtoms;
}

ChargeSeparationIndexDescriptor::ChargeSeparationIndexDescriptor()
    : Descriptor("ChargeSeparationIndex", "Mean graph distance between charged atoms") {}

std::variant<double, int, std::string> ChargeSeparationIndexDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) {
        return "Error: Invalid Molecule";
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    
    std::vector<int> chargedAtoms;
    
    for (unsigned int i = 0; i < rdkMol->getNumAtoms(); ++i) {
        if (rdkMol->getAtomWithIdx(i)->getFormalCharge() != 0) {
            chargedAtoms.push_back(i);
        }
    }
    
    if (chargedAtoms.size() <= 1) {
        return 0.0; // No or only one charged atom
    }
    
    double totalDistance = 0.0;
    int pairCount = 0;
    
    for (size_t i = 0; i < chargedAtoms.size(); ++i) {
        for (size_t j = i + 1; j < chargedAtoms.size(); ++j) {
            int distance = getShortestPath(rdkMol, chargedAtoms[i], chargedAtoms[j]);
            if (distance > 0) {
                totalDistance += distance;
                pairCount++;
            }
        }
    }
    
    if (pairCount == 0) {
        return 0.0;
    }
    
    return totalDistance / pairCount;
}

NeutralAtomRatioDescriptor::NeutralAtomRatioDescriptor()
    : FractionalDescriptor("NeutralAtomRatio", "Atoms with zero formal charge / total atoms") {}

std::variant<double, int, std::string> NeutralAtomRatioDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) {
        return "Error: Invalid Molecule";
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    
    int neutralAtoms = 0;
    int totalAtoms = rdkMol->getNumAtoms();
    
    for (const RDKit::Atom* atom : rdkMol->atoms()) {
        if (atom->getFormalCharge() == 0) {
            neutralAtoms++;
        }
    }
    
    return static_cast<double>(neutralAtoms) / totalAtoms;
}

FunctionalGroupChargeBalanceDescriptor::FunctionalGroupChargeBalanceDescriptor()
    : FractionalDescriptor("FunctionalGroupChargeBalance", "Charged functional groups / total functional groups") {}

std::variant<double, int, std::string> FunctionalGroupChargeBalanceDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) {
        return "Error: Invalid Molecule";
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    
    // Define functional groups using SMARTS patterns
    std::vector<std::string> functionalGroupSmarts = {
        "[OH]", // Alcohol/Phenol
        "[NH,NH2]", // Amine
        "[CX3](=O)[OX2H,OX1-]", // Carboxylic acid
        "[CX3](=O)[OX2]", // Ester
        "[CX3](=O)[NX3]", // Amide
        "[CX3](=O)", // Carbonyl
        "[NX3]", // Amine/Amide
        "[SX2]", // Thiol/Thioether
    };
    
    int totalFunctionalGroups = 0;
    int chargedFunctionalGroups = 0;
    
    // Find functional groups
    for (const std::string& smarts : functionalGroupSmarts) {
        try {
            RDKit::RWMol* query = RDKit::SmartsToMol(smarts);
            if (!query) continue;
            
            std::vector<std::vector<std::pair<int, int>>> matches;
            RDKit::SubstructMatch(*rdkMol, *query, matches);
            
            for (const auto& match : matches) {
                totalFunctionalGroups++;
                
                // Check if any atom in the functional group is charged
                bool hasCharge = false;
                for (const auto& pair : match) {
                    if (rdkMol->getAtomWithIdx(pair.second)->getFormalCharge() != 0) {
                        hasCharge = true;
                        break;
                    }
                }
                
                if (hasCharge) {
                    chargedFunctionalGroups++;
                }
            }
            
            delete query;
        } catch (...) {
            // Skip if SMARTS parsing fails
            continue;
        }
    }
    
    if (totalFunctionalGroups == 0) {
        return 0.0;
    }
    
    return static_cast<double>(chargedFunctionalGroups) / totalFunctionalGroups;
}

// Spatial Configuration Proxies
StereocenterDensityDescriptor::StereocenterDensityDescriptor()
    : FractionalDescriptor("StereocenterDensity", "Chiral centers / total heavy atoms") {}

std::variant<double, int, std::string> StereocenterDensityDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) {
        return "Error: Invalid Molecule";
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    
    int stereocenters = 0;
    int heavyAtoms = 0;
    
    for (const RDKit::Atom* atom : rdkMol->atoms()) {
        if (atom->getAtomicNum() > 1) { // Heavy atoms
            heavyAtoms++;
            
            if (isStereocenter(atom)) {
                stereocenters++;
            }
        }
    }
    
    if (heavyAtoms == 0) {
        return 0.0;
    }
    
    return static_cast<double>(stereocenters) / heavyAtoms;
}

DoubleBondConfigurationDensityDescriptor::DoubleBondConfigurationDensityDescriptor()
    : FractionalDescriptor("DoubleBondConfigurationDensity", "Double bonds with distinct substituents / total double bonds") {}

std::variant<double, int, std::string> DoubleBondConfigurationDensityDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) {
        return "Error: Invalid Molecule";
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    
    int configuredDoubleBonds = 0;
    int totalDoubleBonds = 0;
    
    for (const RDKit::Bond* bond : rdkMol->bonds()) {
        if (bond->getBondType() == RDKit::Bond::BondType::DOUBLE) {
            totalDoubleBonds++;
            
            // Check if this is an E/Z configured double bond
            if (bond->getStereo() == RDKit::Bond::BondStereo::STEREOE ||
                bond->getStereo() == RDKit::Bond::BondStereo::STEREOZ) {
                configuredDoubleBonds++;
            }
        }
    }
    
    if (totalDoubleBonds == 0) {
        return 0.0;
    }
    
    return static_cast<double>(configuredDoubleBonds) / totalDoubleBonds;
}

RingStereogenicAtomDensityDescriptor::RingStereogenicAtomDensityDescriptor()
    : FractionalDescriptor("RingStereogenicAtomDensity", "Chiral centers within rings / total ring atoms") {}

std::variant<double, int, std::string> RingStereogenicAtomDensityDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) {
        return "Error: Invalid Molecule";
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    
    int ringAtoms = 0;
    int ringStereogenicAtoms = 0;
    
    for (const RDKit::Atom* atom : rdkMol->atoms()) {
        if (isAtomInRing(atom)) {
            ringAtoms++;
            
            if (isStereocenter(atom)) {
                ringStereogenicAtoms++;
            }
        }
    }
    
    if (ringAtoms == 0) {
        return 0.0;
    }
    
    return static_cast<double>(ringStereogenicAtoms) / ringAtoms;
}

TerminalStereocenterRatioDescriptor::TerminalStereocenterRatioDescriptor()
    : FractionalDescriptor("TerminalStereocenterRatio", "Chiral terminal atoms / total terminal atoms") {}

std::variant<double, int, std::string> TerminalStereocenterRatioDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) {
        return "Error: Invalid Molecule";
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    
    int terminalAtoms = 0;
    int terminalStereocenters = 0;
    
    for (const RDKit::Atom* atom : rdkMol->atoms()) {
        if (isTerminalAtom(atom) && atom->getAtomicNum() > 1) { // Heavy terminal atoms
            terminalAtoms++;
            
            if (isStereocenter(atom)) {
                terminalStereocenters++;
            }
        }
    }
    
    if (terminalAtoms == 0) {
        return 0.0;
    }
    
    return static_cast<double>(terminalStereocenters) / terminalAtoms;
}

AdjacentStereocenterDensityDescriptor::AdjacentStereocenterDensityDescriptor()
    : FractionalDescriptor("AdjacentStereocenterDensity", "Pairs of directly bonded stereocenters / total stereocenters") {}

std::variant<double, int, std::string> AdjacentStereocenterDensityDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) {
        return "Error: Invalid Molecule";
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    
    int stereocenters = 0;
    int adjacentStereocenters = 0;
    
    for (const RDKit::Atom* atom : rdkMol->atoms()) {
        if (isStereocenter(atom)) {
            stereocenters++;
            
            for (const auto& nbr : rdkMol->atomNeighbors(atom)) {
                if (isStereocenter(nbr)) {
                    adjacentStereocenters++;
                }
            }
        }
    }
    
    if (stereocenters == 0) {
        return 0.0;
    }
    
    return static_cast<double>(adjacentStereocenters) / stereocenters;
}

// Implementations of Constructors for Elemental Composition Diversity Descriptors
GroupPeriodicDiversityDescriptor::GroupPeriodicDiversityDescriptor()
    : Descriptor("GroupPeriodicDiversity", "Diversity index across periodic groups") {}

PeriodDiversityIndexDescriptor::PeriodDiversityIndexDescriptor()
    : Descriptor("PeriodDiversityIndex", "Diversity index across periodic periods") {}

AtomicMassVarianceDescriptor::AtomicMassVarianceDescriptor()
    : Descriptor("AtomicMassVariance", "Variance of atomic masses of constituent atoms") {}

RareElementDensityDescriptor::RareElementDensityDescriptor()
    : FractionalDescriptor("RareElementDensity", "Fraction of atoms from uncommon periodic groups (group number ≥15, excluding common heteroatoms)") {}

AlkaliAlkalineEarthRatioDescriptor::AlkaliAlkalineEarthRatioDescriptor()
    : FractionalDescriptor("AlkaliAlkalineEarthRatio", "Alkali and alkaline earth atoms / total atoms") {}

// Implementations of Constructors for Bonding Environment Descriptors
UnsaturationClustersDescriptor::UnsaturationClustersDescriptor()
    : FractionalDescriptor("UnsaturationClusters", "Clusters of ≥2 adjacent unsaturated bonds (double or triple) / total unsaturated bonds") {}

RingSpanningBondDensityDescriptor::RingSpanningBondDensityDescriptor()
    : FractionalDescriptor("RingSpanningBondDensity", "Bonds connecting different rings / total bonds") {}

TerminalUnsaturatedBondRatioDescriptor::TerminalUnsaturatedBondRatioDescriptor()
    : FractionalDescriptor("TerminalUnsaturatedBondRatio", "Terminal unsaturated bonds / total unsaturated bonds") {}

HeavyAtomBondOrderVarianceDescriptor::HeavyAtomBondOrderVarianceDescriptor()
    : Descriptor("HeavyAtomBondOrderVariance", "Variance in bond orders per heavy atom") {}

HeteroatomBondingDiversityDescriptor::HeteroatomBondingDiversityDescriptor()
    : Descriptor("HeteroatomBondingDiversity", "Average distinct atom types bonded to each heteroatom") {}

// Now let's complete the implementation of the calculate methods for each missing descriptor

std::variant<double, int, std::string> GroupPeriodicDiversityDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) {
        return "Error: Invalid Molecule";
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    
    // Map of periodic groups to count
    std::unordered_map<int, int> groupCounts;
    
    // Define group mappings (simplified)
    static const std::unordered_map<int, int> elementGroups = {
        {1, 1}, {3, 1}, {11, 1}, {19, 1}, {37, 1}, {55, 1}, // Group 1 (H, Li, Na, K, Rb, Cs)
        {4, 2}, {12, 2}, {20, 2}, {38, 2}, {56, 2}, // Group 2 (Be, Mg, Ca, Sr, Ba)
        {5, 13}, {13, 13}, {31, 13}, {49, 13}, {81, 13}, // Group 13 (B, Al, Ga, In, Tl)
        {6, 14}, {14, 14}, {32, 14}, {50, 14}, {82, 14}, // Group 14 (C, Si, Ge, Sn, Pb)
        {7, 15}, {15, 15}, {33, 15}, {51, 15}, {83, 15}, // Group 15 (N, P, As, Sb, Bi)
        {8, 16}, {16, 16}, {34, 16}, {52, 16}, {84, 16}, // Group 16 (O, S, Se, Te, Po)
        {9, 17}, {17, 17}, {35, 17}, {53, 17}, {85, 17}, // Group 17 (F, Cl, Br, I, At)
        {2, 18}, {10, 18}, {18, 18}, {36, 18}, {54, 18}, {86, 18} // Group 18 (He, Ne, Ar, Kr, Xe, Rn)
    };
    
    for (const RDKit::Atom* atom : rdkMol->atoms()) {
        int atomicNum = atom->getAtomicNum();
        if (atomicNum > 1) { // Skip hydrogens
            auto it = elementGroups.find(atomicNum);
            int group = (it != elementGroups.end()) ? it->second : 0;
            
            if (group > 0) {
                groupCounts[group]++;
            }
        }
    }
    
    std::vector<int> counts;
    for (const auto& pair : groupCounts) {
        counts.push_back(pair.second);
    }
    
    return calculateDiversityIndex(counts);
}

std::variant<double, int, std::string> PeriodDiversityIndexDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) {
        return "Error: Invalid Molecule";
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    
    // Map of periodic periods to count
    std::unordered_map<int, int> periodCounts;
    
    // Define period mappings (simplified)
    static const std::unordered_map<int, int> elementPeriods = {
        {1, 1}, {2, 1}, // Period 1 (H, He)
        {3, 2}, {4, 2}, {5, 2}, {6, 2}, {7, 2}, {8, 2}, {9, 2}, {10, 2}, // Period 2
        {11, 3}, {12, 3}, {13, 3}, {14, 3}, {15, 3}, {16, 3}, {17, 3}, {18, 3}, // Period 3
        {19, 4}, {20, 4}, {31, 4}, {32, 4}, {33, 4}, {34, 4}, {35, 4}, {36, 4}, // Period 4 (main group)
        {37, 5}, {38, 5}, {49, 5}, {50, 5}, {51, 5}, {52, 5}, {53, 5}, {54, 5}, // Period 5 (main group)
        {55, 6}, {56, 6}, {81, 6}, {82, 6}, {83, 6}, {84, 6}, {85, 6}, {86, 6}  // Period 6 (main group)
    };
    
    for (const RDKit::Atom* atom : rdkMol->atoms()) {
        int atomicNum = atom->getAtomicNum();
        if (atomicNum > 1) { // Skip hydrogens
            auto it = elementPeriods.find(atomicNum);
            int period = (it != elementPeriods.end()) ? it->second : 0;
            
            if (period > 0) {
                periodCounts[period]++;
            }
        }
    }
    
    std::vector<int> counts;
    for (const auto& pair : periodCounts) {
        counts.push_back(pair.second);
    }
    
    return calculateDiversityIndex(counts);
}

std::variant<double, int, std::string> AtomicMassVarianceDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) {
        return "Error: Invalid Molecule";
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    
    std::vector<double> atomicMasses;
    
    for (const RDKit::Atom* atom : rdkMol->atoms()) {
        if (atom->getAtomicNum() > 1) { // Skip hydrogens
            atomicMasses.push_back(atom->getMass());
        }
    }
    
    if (atomicMasses.empty()) {
        return 0.0;
    }
    
    return calculateVariance(atomicMasses);
}

std::variant<double, int, std::string> RareElementDensityDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) {
        return "Error: Invalid Molecule";
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    
    int rareElements = 0;
    int totalHeavyAtoms = 0;
    
    for (const RDKit::Atom* atom : rdkMol->atoms()) {
        if (atom->getAtomicNum() > 1) { // Heavy atoms
            totalHeavyAtoms++;
            
            if (isRareElement(atom)) {
                rareElements++;
            }
        }
    }
    
    if (totalHeavyAtoms == 0) {
        return 0.0;
    }
    
    return static_cast<double>(rareElements) / totalHeavyAtoms;
}

std::variant<double, int, std::string> AlkaliAlkalineEarthRatioDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) {
        return "Error: Invalid Molecule";
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    
    int alkaliAlkalineEarthAtoms = 0;
    int totalHeavyAtoms = 0;
    
    for (const RDKit::Atom* atom : rdkMol->atoms()) {
        if (atom->getAtomicNum() > 1) { // Heavy atoms
            totalHeavyAtoms++;
            
            if (isAlkaliOrAlkalineEarth(atom)) {
                alkaliAlkalineEarthAtoms++;
            }
        }
    }
    
    if (totalHeavyAtoms == 0) {
        return 0.0;
    }
    
    return static_cast<double>(alkaliAlkalineEarthAtoms) / totalHeavyAtoms;
}

std::variant<double, int, std::string> UnsaturationClustersDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) {
        return "Error: Invalid Molecule";
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    
    int totalUnsaturatedBonds = 0;
    int clusteredUnsaturatedBonds = 0;
    
    // Count total unsaturated bonds
    for (const RDKit::Bond* bond : rdkMol->bonds()) {
        if (isUnsaturatedBond(bond)) {
            totalUnsaturatedBonds++;
        }
    }
    
    if (totalUnsaturatedBonds == 0) {
        return 0.0;
    }
    
    // Find clusters of unsaturated bonds
    std::vector<bool> visited(rdkMol->getNumBonds(), false);
    
    for (unsigned int i = 0; i < rdkMol->getNumBonds(); ++i) {
        if (visited[i] || !isUnsaturatedBond(rdkMol->getBondWithIdx(i))) {
            continue;
        }
        
        std::vector<unsigned int> cluster;
        std::queue<unsigned int> queue;
        queue.push(i);
        visited[i] = true;
        
        while (!queue.empty()) {
            unsigned int current = queue.front();
            queue.pop();
            cluster.push_back(current);
            
            const RDKit::Bond* bond = rdkMol->getBondWithIdx(current);
            unsigned int beginAtomIdx = bond->getBeginAtomIdx();
            unsigned int endAtomIdx = bond->getEndAtomIdx();
            
            // Check all bonds connected to either end of this bond
            for (const RDKit::Bond* nbr : rdkMol->bonds()) {
                if (nbr->getIdx() == current) continue;
                
                if ((nbr->getBeginAtomIdx() == beginAtomIdx || nbr->getEndAtomIdx() == beginAtomIdx ||
                     nbr->getBeginAtomIdx() == endAtomIdx || nbr->getEndAtomIdx() == endAtomIdx) &&
                    isUnsaturatedBond(nbr) && !visited[nbr->getIdx()]) {
                    queue.push(nbr->getIdx());
                    visited[nbr->getIdx()] = true;
                }
            }
        }
        
        if (cluster.size() >= 2) {
            clusteredUnsaturatedBonds += cluster.size();
        }
    }
    
    return static_cast<double>(clusteredUnsaturatedBonds) / totalUnsaturatedBonds;
}

std::variant<double, int, std::string> RingSpanningBondDensityDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) {
        return "Error: Invalid Molecule";
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    const RDKit::RingInfo* ringInfo = rdkMol->getRingInfo();
    
    int ringSpanningBonds = 0;
    int totalBonds = rdkMol->getNumBonds();
    
    if (totalBonds == 0) {
        return 0.0;
    }
    
    for (const RDKit::Bond* bond : rdkMol->bonds()) {
        unsigned int beginAtomIdx = bond->getBeginAtomIdx();
        unsigned int endAtomIdx = bond->getEndAtomIdx();
        
        bool beginInRing = ringInfo->numAtomRings(beginAtomIdx) > 0;
        bool endInRing = ringInfo->numAtomRings(endAtomIdx) > 0;
        
        if (beginInRing && endInRing) {
            bool differentRings = false;
            for (unsigned int i = 0; i < ringInfo->numRings(); ++i) {
                const std::vector<int>& ringAtoms = ringInfo->atomRings()[i];
                bool containsBegin = std::find(ringAtoms.begin(), ringAtoms.end(), beginAtomIdx) != ringAtoms.end();
                bool containsEnd = std::find(ringAtoms.begin(), ringAtoms.end(), endAtomIdx) != ringAtoms.end();
                
                if (containsBegin != containsEnd) {
                    differentRings = true;
                    break;
                }
            }
            
            if (differentRings) {
                ringSpanningBonds++;
            }
        }
    }
    
    return static_cast<double>(ringSpanningBonds) / totalBonds;
}

std::variant<double, int, std::string> TerminalUnsaturatedBondRatioDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) {
        return "Error: Invalid Molecule";
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    
    int terminalUnsaturatedBonds = 0;
    int totalUnsaturatedBonds = 0;
    
    for (const RDKit::Bond* bond : rdkMol->bonds()) {
        if (isUnsaturatedBond(bond)) {
            totalUnsaturatedBonds++;
            
            if (isTerminalBond(bond)) {
                terminalUnsaturatedBonds++;
            }
        }
    }
    
    if (totalUnsaturatedBonds == 0) {
        return 0.0;
    }
    
    return static_cast<double>(terminalUnsaturatedBonds) / totalUnsaturatedBonds;
}

std::variant<double, int, std::string> HeavyAtomBondOrderVarianceDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) {
        return "Error: Invalid Molecule";
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    
    std::vector<double> atomBondOrderVariances;
    
    for (const RDKit::Atom* atom : rdkMol->atoms()) {
        if (atom->getAtomicNum() > 1) { // Heavy atoms
            std::vector<double> bondOrders;
            
            for (const RDKit::Bond* bond : rdkMol->atomBonds(atom)) {
                bondOrders.push_back(getBondOrder(bond));
            }
            
            if (!bondOrders.empty()) {
                double variance = calculateVariance(bondOrders);
                atomBondOrderVariances.push_back(variance);
            }
        }
    }
    
    if (atomBondOrderVariances.empty()) {
        return 0.0;
    }
    
    // Return the average variance across all heavy atoms
    double sum = std::accumulate(atomBondOrderVariances.begin(), atomBondOrderVariances.end(), 0.0);
    return sum / atomBondOrderVariances.size();
}

std::variant<double, int, std::string> HeteroatomBondingDiversityDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) {
        return "Error: Invalid Molecule";
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    
    int totalHeteroatoms = 0;
    double totalDiversity = 0.0;
    
    for (const RDKit::Atom* atom : rdkMol->atoms()) {
        if (isHeteroatom(atom)) {
            totalHeteroatoms++;
            
            std::unordered_set<int> neighborAtomTypes;
            for (const auto& nbr : rdkMol->atomNeighbors(atom)) {
                neighborAtomTypes.insert(nbr->getAtomicNum());
            }
            
            totalDiversity += neighborAtomTypes.size();
        }
    }
    
    if (totalHeteroatoms == 0) {
        return 0.0;
    }
    
    return totalDiversity / totalHeteroatoms;
}

LocalizedChargeClustersDescriptor::LocalizedChargeClustersDescriptor()
    : Descriptor("LocalizedChargeClusters", "Fraction of atoms with significant partial charges") {}

std::variant<double, int, std::string> LocalizedChargeClustersDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return 0.0;
    
    // Compute partial charges if needed
    std::vector<double> charges;
    RDKit::RWMol nonConstMol(*rdkMol);
    RDKit::computeGasteigerCharges(nonConstMol);
    
    charges.resize(rdkMol->getNumAtoms());
    for (unsigned int i = 0; i < rdkMol->getNumAtoms(); ++i) {
        if (nonConstMol.getAtomWithIdx(i)->hasProp(RDKit::common_properties::_GasteigerCharge)) {
            charges[i] = nonConstMol.getAtomWithIdx(i)->getProp<double>(RDKit::common_properties::_GasteigerCharge);
        } else {
            charges[i] = 0.0;
        }
    }
    
    // Lower the threshold for charge clusters
    const double chargeThreshold = 0.05; // Lower from original 0.1
    
    // Find clusters of atoms with significant partial charges
    std::vector<bool> visited(rdkMol->getNumAtoms(), false);
    int clusterCount = 0;
    
    for (unsigned int i = 0; i < rdkMol->getNumAtoms(); ++i) {
        if (visited[i] || std::abs(charges[i]) < chargeThreshold) continue;
        
        // Start a new cluster
        std::queue<unsigned int> queue;
        queue.push(i);
        visited[i] = true;
        bool validCluster = false;
        
        while (!queue.empty()) {
            unsigned int current = queue.front();
            queue.pop();
            
            // Look at neighbors
            for (const auto& nbr : rdkMol->atomNeighbors(rdkMol->getAtomWithIdx(current))) {
                unsigned int nbrIdx = nbr->getIdx();
                if (!visited[nbrIdx] && std::abs(charges[nbrIdx]) >= chargeThreshold) {
                    if (charges[nbrIdx] * charges[current] < 0) { // Opposite charges
                        validCluster = true;
                    }
                    queue.push(nbrIdx);
                    visited[nbrIdx] = true;
                }
            }
        }
        
        if (validCluster) {
            clusterCount++;
        }
    }
    
    // Return the count normalized by heavy atom count
    double heavyAtomCount = 0;
    for (const RDKit::Atom* atom : rdkMol->atoms()) {
        if (atom->getAtomicNum() > 1) heavyAtomCount++;
    }
    
    if (heavyAtomCount < 1) return 0.0;
    return static_cast<double>(clusterCount) / heavyAtomCount;
}
}
}