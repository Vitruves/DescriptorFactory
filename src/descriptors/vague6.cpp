#include "descriptors/vague6.hpp"
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
#include <GraphMol/RingInfo.h>
#include <numeric>
#include <algorithm>
#include <cmath>
#include <queue>
#include <stack>
#include <map>
#include <set>
#include <algorithm>
#include <limits>
#include <functional>

namespace desfact {
namespace descriptors {

// --- Helper Functions (Anonymous Namespace) ---
namespace {
    // Helper functions for calculations
    
    // Calculate Shannon entropy from a frequency map
    double calculateShannonEntropy(const std::map<int, int>& frequencies, int totalCount) {
        double entropy = 0.0;
        for (const auto& pair : frequencies) {
            double p = static_cast<double>(pair.second) / totalCount;
            if (p > 0) entropy -= p * std::log2(p);
        }
        return entropy;
    }
    
    // Calculate Gini coefficient from a vector of values
    double calculateGiniCoefficient(std::vector<double> values) {
        if (values.empty()) return 0.0;
        std::sort(values.begin(), values.end());
        double n = static_cast<double>(values.size());
        double sum = 0.0;
        for (size_t i = 0; i < values.size(); ++i) {
            sum += (2.0 * (i + 1) - n - 1) * values[i];
        }
        double totalSum = std::accumulate(values.begin(), values.end(), 0.0);
        if (totalSum == 0.0) return 0.0;
        return sum / (n * totalSum);
    }
    
    // Calculate Fisher's skewness coefficient
    double calculateFisherSkewness(const std::vector<double>& values) {
        if (values.size() < 3) return 0.0;
        
        double mean = std::accumulate(values.begin(), values.end(), 0.0) / values.size();
        
        double variance = 0.0;
        double skewness = 0.0;
        for (double val : values) {
            double diff = val - mean;
            double diff2 = diff * diff;
            variance += diff2;
            skewness += diff * diff2;
        }
        
        variance /= values.size();
        if (variance < 1e-10) return 0.0;
        
        double stdDev = std::sqrt(variance);
        skewness /= (values.size() * std::pow(stdDev, 3));
        
        return skewness;
    }
    
    // Check if atom is a heteroatom (not C or H)
    bool isHeteroatom(const RDKit::Atom* atom) {
        int atomicNum = atom->getAtomicNum();
        return atomicNum > 1 && atomicNum != 6;
    }
    
    // Check if atom is a halogen
    bool isHalogen(const RDKit::Atom* atom) {
        int atomicNum = atom->getAtomicNum();
        return atomicNum == 9 || atomicNum == 17 || atomicNum == 35 || atomicNum == 53 || atomicNum == 85;
    }
    
    // Get Pauling electronegativity for an atom
    double getAtomENPauling(const RDKit::Atom* atom) {
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
        };
        
        auto it = electronegativity.find(atom->getAtomicNum());
        return (it != electronegativity.end()) ? it->second : 0.0;
    }
    
    // Compute shortest paths between all pairs of atoms
    std::vector<std::vector<int>> computeAllPairsShortestPaths(const RDKit::ROMol& mol) {
        int numAtoms = mol.getNumAtoms();
        std::vector<std::vector<int>> dist(numAtoms, std::vector<int>(numAtoms, std::numeric_limits<int>::max() / 2));
        
        // Initialize with direct connections
        for (const auto& bond : mol.bonds()) {
            int u = bond->getBeginAtomIdx();
            int v = bond->getEndAtomIdx();
            dist[u][v] = 1;
            dist[v][u] = 1;
        }
        
        // Self-loops are 0
        for (int i = 0; i < numAtoms; ++i) {
            dist[i][i] = 0;
        }
        
        // Floyd-Warshall algorithm
        for (int k = 0; k < numAtoms; ++k) {
            for (int i = 0; i < numAtoms; ++i) {
                for (int j = 0; j < numAtoms; ++j) {
                    if (dist[i][k] + dist[k][j] < dist[i][j]) {
                        dist[i][j] = dist[i][k] + dist[k][j];
                    }
                }
            }
        }
        
        return dist;
    }
}

// --- Descriptor Implementations ---

// 1. DegreeLogWeightedSum
DegreeLogWeightedSum::DegreeLogWeightedSum()
    : Descriptor("DegreeLogWeightedSum", "Σ deg(i) · ln [deg(i)+1] over all heavy atoms") {}

std::variant<double, int, std::string> DegreeLogWeightedSum::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    
    const RDKit::ROMol& rdkMol = *mol.getMolecule();
    double sum = 0.0;
    
    for (const auto& atom : rdkMol.atoms()) {
        if (atom->getAtomicNum() > 1) { // Heavy atoms only
            int degree = atom->getDegree();
            sum += degree * std::log(degree + 1.0);
        }
    }
    
    return sum;
}

// 2. BondOrderSkewness
BondOrderSkewness::BondOrderSkewness()
    : Descriptor("BondOrderSkewness", "Fisher skewness of the bond-order distribution") {}

std::variant<double, int, std::string> BondOrderSkewness::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    
    const RDKit::ROMol& rdkMol = *mol.getMolecule();
    std::vector<double> bondOrders;
    
    for (const auto& bond : rdkMol.bonds()) {
        double order = 0.0;
        if (bond->getBondType() == RDKit::Bond::SINGLE) order = 1.0;
        else if (bond->getBondType() == RDKit::Bond::DOUBLE) order = 2.0;
        else if (bond->getBondType() == RDKit::Bond::TRIPLE) order = 3.0;
        else if (bond->getBondType() == RDKit::Bond::AROMATIC) order = 1.5;
        bondOrders.push_back(order);
    }
    
    if (bondOrders.size() < 3) return 0.0;
    return calculateFisherSkewness(bondOrders);
}

// 3. PathLengthGini
PathLengthGini::PathLengthGini()
    : Descriptor("PathLengthGini", "Gini coefficient of all heavy-atom shortest-path lengths") {}

std::variant<double, int, std::string> PathLengthGini::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    
    const RDKit::ROMol& rdkMol = *mol.getMolecule();
    int numHeavyAtoms = 0;
    
    // Count heavy atoms
    for (const auto& atom : rdkMol.atoms()) {
        if (atom->getAtomicNum() > 1) numHeavyAtoms++;
    }
    
    if (numHeavyAtoms <= 1) return 0.0; // No paths to measure
    
    // Compute all shortest paths
    std::vector<std::vector<int>> shortestPaths = computeAllPairsShortestPaths(rdkMol);
    
    // Extract path lengths between heavy atoms
    std::vector<double> pathLengths;
    for (int i = 0; i < rdkMol.getNumAtoms(); ++i) {
        if (rdkMol.getAtomWithIdx(i)->getAtomicNum() <= 1) continue;
        
        for (int j = i + 1; j < rdkMol.getNumAtoms(); ++j) {
            if (rdkMol.getAtomWithIdx(j)->getAtomicNum() <= 1) continue;
            
            pathLengths.push_back(static_cast<double>(shortestPaths[i][j]));
        }
    }
    
    return calculateGiniCoefficient(pathLengths);
}

// 4. LongestHomoelementPath
LongestHomoelementPath::LongestHomoelementPath()
    : Descriptor("LongestHomoelementPath", "Length of the longest path whose atoms share an identical element symbol") {}

std::variant<double, int, std::string> LongestHomoelementPath::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0;
    
    const RDKit::ROMol& rdkMol = *mol.getMolecule();
    int numAtoms = rdkMol.getNumAtoms();
    if (numAtoms <= 1) return 1; // Single atom is its own path
    
    // Group atoms by element type
    std::map<int, std::vector<int>> elementGroups;
    for (int i = 0; i < numAtoms; ++i) {
        int atomicNum = rdkMol.getAtomWithIdx(i)->getAtomicNum();
        if (atomicNum > 1) { // Only consider heavy atoms
            elementGroups[atomicNum].push_back(i);
        }
    }
    
    int longestPath = 1; // Minimum path length is 1 (single atom)
    
    // For each element type, find the longest path
    for (const auto& group : elementGroups) {
        if (group.second.size() <= 1) continue; // Need at least 2 atoms for a path
        
        // Build adjacency list for atoms of this element
        std::vector<std::vector<int>> adj(numAtoms);
        for (const auto& bond : rdkMol.bonds()) {
            int u = bond->getBeginAtomIdx();
            int v = bond->getEndAtomIdx();
            
            if (rdkMol.getAtomWithIdx(u)->getAtomicNum() == group.first && 
                rdkMol.getAtomWithIdx(v)->getAtomicNum() == group.first) {
                adj[u].push_back(v);
                adj[v].push_back(u);
            }
        }
        
        // Find longest path using BFS from each atom in the group
        for (int startAtom : group.second) {
            std::vector<bool> visited(numAtoms, false);
            std::vector<int> distance(numAtoms, 0);
            std::queue<int> q;
            
            q.push(startAtom);
            visited[startAtom] = true;
            
            while (!q.empty()) {
                int u = q.front();
                q.pop();
                
                for (int v : adj[u]) {
                    if (!visited[v]) {
                        visited[v] = true;
                        distance[v] = distance[u] + 1;
                        q.push(v);
                        
                        longestPath = std::max(longestPath, distance[v] + 1); // +1 to count nodes not edges
                    }
                }
            }
        }
    }
    
    return longestPath;
}

// 5. PeripheralAtomTypeDiversity
PeripheralAtomTypeDiversity::PeripheralAtomTypeDiversity()
    : Descriptor("PeripheralAtomTypeDiversity", "Shannon entropy of element types among terminal heavy atoms") {}

std::variant<double, int, std::string> PeripheralAtomTypeDiversity::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    
    const RDKit::ROMol& rdkMol = *mol.getMolecule();
    std::map<int, int> terminalElementFreq; // Map of atomic number to frequency
    int totalTerminalAtoms = 0;
    
    for (const auto& atom : rdkMol.atoms()) {
        if (atom->getAtomicNum() > 1 && atom->getDegree() == 1) { // Terminal heavy atom
            terminalElementFreq[atom->getAtomicNum()]++;
            totalTerminalAtoms++;
        }
    }
    
    if (totalTerminalAtoms == 0) return 0.0;
    
    return calculateShannonEntropy(terminalElementFreq, totalTerminalAtoms);
}

// 6. HeteroRingEdgeRatio
HeteroRingEdgeRatio::HeteroRingEdgeRatio()
    : Descriptor("HeteroRingEdgeRatio", "Heteroatoms directly attached to rings but not ring members ÷ total heteroatoms") {}

std::variant<double, int, std::string> HeteroRingEdgeRatio::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    
    const RDKit::ROMol& rdkMol = *mol.getMolecule();
    RDKit::RingInfo* ringInfo = rdkMol.getRingInfo();
    if (!ringInfo->isInitialized()) {
        RDKit::MolOps::findSSSR(const_cast<RDKit::ROMol&>(rdkMol));
    }
    
    int heteroRingEdges = 0; // Heteroatoms directly attached to rings but not in rings
    int totalHeteroatoms = 0; // Total heteroatoms
    
    for (const auto& atom : rdkMol.atoms()) {
        if (isHeteroatom(atom)) {
            totalHeteroatoms++;
            
            // Check if atom is not in a ring but connected to a ring atom
            if (!ringInfo->numAtomRings(atom->getIdx())) {
                bool connectedToRing = false;
                for (const auto& nbr : rdkMol.atomNeighbors(atom)) {
                    if (ringInfo->numAtomRings(nbr->getIdx())) {
                        connectedToRing = true;
                        break;
                    }
                }
                if (connectedToRing) {
                    heteroRingEdges++;
                }
            }
        }
    }
    
    if (totalHeteroatoms == 0) return 0.0;
    return static_cast<double>(heteroRingEdges) / totalHeteroatoms;
}

// 7. BranchPointDensity
BranchPointDensity::BranchPointDensity()
    : Descriptor("BranchPointDensity", "Branching atoms (degree > 2) per heavy atom") {}

std::variant<double, int, std::string> BranchPointDensity::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    
    const RDKit::ROMol& rdkMol = *mol.getMolecule();
    int branchPoints = 0;
    int heavyAtoms = 0;
    
    for (const auto& atom : rdkMol.atoms()) {
        if (atom->getAtomicNum() > 1) { // Heavy atom
            heavyAtoms++;
            if (atom->getDegree() > 2) { // Branching atom
                branchPoints++;
            }
        }
    }
    
    if (heavyAtoms == 0) return 0.0;
    return static_cast<double>(branchPoints) / heavyAtoms;
}

// 8. MeanBranchSeparation
MeanBranchSeparation::MeanBranchSeparation()
    : Descriptor("MeanBranchSeparation", "Mean shortest-path distance between all pairs of branch points") {}

std::variant<double, int, std::string> MeanBranchSeparation::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    
    const RDKit::ROMol& rdkMol = *mol.getMolecule();
    
    // Identify branch points
    std::vector<int> branchPoints;
    for (const auto& atom : rdkMol.atoms()) {
        if (atom->getAtomicNum() > 1 && atom->getDegree() > 2) {
            branchPoints.push_back(atom->getIdx());
        }
    }
    
    if (branchPoints.size() < 2) return 0.0; // Need at least 2 branch points
    
    // Compute all shortest paths
    std::vector<std::vector<int>> shortestPaths = computeAllPairsShortestPaths(rdkMol);
    
    // Calculate mean distance between branch points
    double totalDistance = 0.0;
    int pairCount = 0;
    
    for (size_t i = 0; i < branchPoints.size(); ++i) {
        for (size_t j = i + 1; j < branchPoints.size(); ++j) {
            totalDistance += shortestPaths[branchPoints[i]][branchPoints[j]];
            pairCount++;
        }
    }
    
    if (pairCount == 0) return 0.0;
    return totalDistance / pairCount;
}

// 9. CarbonChainBranchingIndex
CarbonChainBranchingIndex::CarbonChainBranchingIndex()
    : Descriptor("CarbonChainBranchingIndex", "Branch points located on the longest carbon chain ÷ chain length") {}

std::variant<double, int, std::string> CarbonChainBranchingIndex::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    
    const RDKit::ROMol& rdkMol = *mol.getMolecule();
    int numAtoms = rdkMol.getNumAtoms();
    
    // Create adjacency list for carbon atoms only
    std::vector<std::vector<int>> carbonAdj(numAtoms);
    for (const auto& bond : rdkMol.bonds()) {
        int u = bond->getBeginAtomIdx();
        int v = bond->getEndAtomIdx();
        
        if (rdkMol.getAtomWithIdx(u)->getAtomicNum() == 6 && 
            rdkMol.getAtomWithIdx(v)->getAtomicNum() == 6) {
            carbonAdj[u].push_back(v);
            carbonAdj[v].push_back(u);
        }
    }
    
    // Find the longest carbon chain using BFS
    std::vector<int> longestChain;
    int maxLength = 0;
    
    for (int i = 0; i < numAtoms; ++i) {
        if (rdkMol.getAtomWithIdx(i)->getAtomicNum() != 6) continue;
        
        // BFS from this carbon atom
        std::vector<bool> visited(numAtoms, false);
        std::vector<int> distance(numAtoms, 0);
        std::vector<int> parent(numAtoms, -1);
        std::queue<int> q;
        
        q.push(i);
        visited[i] = true;
        
        int farthestNode = i;
        int maxDist = 0;
        
        while (!q.empty()) {
            int u = q.front();
            q.pop();
            
            for (int v : carbonAdj[u]) {
                if (!visited[v]) {
                    visited[v] = true;
                    distance[v] = distance[u] + 1;
                    parent[v] = u;
                    q.push(v);
                    
                    if (distance[v] > maxDist) {
                        maxDist = distance[v];
                        farthestNode = v;
                    }
                }
            }
        }
        
        // Reconstruct the path from i to farthestNode
        if (maxDist > maxLength) {
            maxLength = maxDist;
            longestChain.clear();
            
            int current = farthestNode;
            while (current != -1) {
                longestChain.push_back(current);
                current = parent[current];
            }
        }
    }
    
    if (longestChain.empty()) return 0.0;
    
    // Count branch points on the longest chain
    int branchPoints = 0;
    std::set<int> chainAtoms(longestChain.begin(), longestChain.end());
    
    for (int atom : longestChain) {
        int branchCount = 0;
        for (const auto& nbr : rdkMol.atomNeighbors(rdkMol.getAtomWithIdx(atom))) {
            if (nbr->getAtomicNum() == 6 && chainAtoms.find(nbr->getIdx()) == chainAtoms.end()) {
                branchCount++;
            }
        }
        if (branchCount > 0) branchPoints++;
    }
    
    return static_cast<double>(branchPoints) / longestChain.size();
}

// 10. RingSystemCount
RingSystemCount::RingSystemCount()
    : Descriptor("RingSystemCount", "Number of distinct fused-ring systems (SSSR clusters)") {}

std::variant<double, int, std::string> RingSystemCount::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0;
    
    const RDKit::ROMol& rdkMol = *mol.getMolecule();
    RDKit::RingInfo* ringInfo = rdkMol.getRingInfo();
    if (!ringInfo->isInitialized()) {
        RDKit::MolOps::findSSSR(const_cast<RDKit::ROMol&>(rdkMol));
    }
    
    const auto& rings = ringInfo->atomRings();
    if (rings.empty()) return 0;
    
    // Build a graph where nodes are rings and edges represent ring fusion
    int numRings = rings.size();
    std::vector<std::vector<int>> ringAdjList(numRings);
    
    for (int i = 0; i < numRings; ++i) {
        const auto& ringI = rings[i];
        std::set<int> atomsInRingI(ringI.begin(), ringI.end());
        
        for (int j = i + 1; j < numRings; ++j) {
            const auto& ringJ = rings[j];
            
            // Check if rings share any atoms (fusion)
            bool fused = false;
            for (int atom : ringJ) {
                if (atomsInRingI.find(atom) != atomsInRingI.end()) {
                    fused = true;
                    break;
                }
            }
            
            if (fused) {
                ringAdjList[i].push_back(j);
                ringAdjList[j].push_back(i);
            }
        }
    }
    
    // Count connected components (ring systems) using DFS
    std::vector<bool> visited(numRings, false);
    int systemCount = 0;
    
    for (int i = 0; i < numRings; ++i) {
        if (!visited[i]) {
            systemCount++;
            
            // DFS to mark all rings in this system
            std::stack<int> stack;
            stack.push(i);
            visited[i] = true;
            
            while (!stack.empty()) {
                int ring = stack.top();
                stack.pop();
                
                for (int adjRing : ringAdjList[ring]) {
                    if (!visited[adjRing]) {
                        visited[adjRing] = true;
                        stack.push(adjRing);
                    }
                }
            }
        }
    }
    
    return systemCount;
}

// 11. MeanRingPerimeterDegree
MeanRingPerimeterDegree::MeanRingPerimeterDegree()
    : Descriptor("MeanRingPerimeterDegree", "Average atomic degree among all atoms that belong to rings") {}

std::variant<double, int, std::string> MeanRingPerimeterDegree::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    
    const RDKit::ROMol& rdkMol = *mol.getMolecule();
    RDKit::RingInfo* ringInfo = rdkMol.getRingInfo();
    if (!ringInfo->isInitialized()) {
        RDKit::MolOps::findSSSR(const_cast<RDKit::ROMol&>(rdkMol));
    }
    
    double totalDegree = 0.0;
    int ringAtomCount = 0;
    std::set<int> ringAtoms;
    
    // Collect all atoms that are part of at least one ring
    for (const auto& ring : ringInfo->atomRings()) {
        for (int atomIdx : ring) {
            ringAtoms.insert(atomIdx);
        }
    }
    
    // Calculate the mean degree of ring atoms
    for (int atomIdx : ringAtoms) {
        const RDKit::Atom* atom = rdkMol.getAtomWithIdx(atomIdx);
        totalDegree += atom->getDegree();
        ringAtomCount++;
    }
    
    if (ringAtomCount == 0) return 0.0;
    return totalDegree / ringAtomCount;
}

// 12. AtomTypeNeighborDiversity
AtomTypeNeighborDiversity::AtomTypeNeighborDiversity()
    : Descriptor("AtomTypeNeighborDiversity", "For every element type, count distinct neighbor element types; mean over elements") {}

std::variant<double, int, std::string> AtomTypeNeighborDiversity::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    
    const RDKit::ROMol& rdkMol = *mol.getMolecule();
    
    // Map of element type -> set of neighbor element types
    std::map<int, std::set<int>> elementNeighborTypes;
    
    // For each atom, record its neighbors' element types
    for (const auto& atom : rdkMol.atoms()) {
        int atomicNum = atom->getAtomicNum();
        if (atomicNum <= 1) continue; // Skip H atoms
        
        for (const auto& nbr : rdkMol.atomNeighbors(atom)) {
            int nbrAtomicNum = nbr->getAtomicNum();
            if (nbrAtomicNum > 1) { // Only consider heavy atom neighbors
                elementNeighborTypes[atomicNum].insert(nbrAtomicNum);
            }
        }
    }
    
    // Calculate mean diversity
    double totalDiversity = 0.0;
    int elementCount = 0;
    
    for (const auto& pair : elementNeighborTypes) {
        totalDiversity += pair.second.size();
        elementCount++;
    }
    
    if (elementCount == 0) return 0.0;
    return totalDiversity / elementCount;
}

// 13. BondTypeAlternationRatio
BondTypeAlternationRatio::BondTypeAlternationRatio()
    : Descriptor("BondTypeAlternationRatio", "Bonds that are part of an alternating single/double pattern of length ≥ 3 ÷ total bonds") {}

std::variant<double, int, std::string> BondTypeAlternationRatio::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    
    const RDKit::ROMol& rdkMol = *mol.getMolecule();
    int totalBonds = rdkMol.getNumBonds();
    if (totalBonds == 0) return 0.0;
    
    // Build adjacency list with bond types
    int numAtoms = rdkMol.getNumAtoms();
    std::vector<std::vector<std::pair<int, RDKit::Bond::BondType>>> adj(numAtoms);
    
    for (const auto& bond : rdkMol.bonds()) {
        int u = bond->getBeginAtomIdx();
        int v = bond->getEndAtomIdx();
        RDKit::Bond::BondType bondType = bond->getBondType();
        
        adj[u].push_back({v, bondType});
        adj[v].push_back({u, bondType});
    }
    
    // Count bonds in alternating patterns
    std::set<std::pair<int, int>> alternatingBonds; // Set of (u,v) pairs to avoid duplicates
    
    for (int startAtom = 0; startAtom < numAtoms; ++startAtom) {
        // Try DFS from each atom
        std::vector<bool> visited(numAtoms, false);
        std::vector<RDKit::Bond::BondType> pathBondTypes;
        std::vector<std::pair<int, int>> pathBonds;
        
        std::function<void(int, RDKit::Bond::BondType, int)> dfs = [&](int atom, RDKit::Bond::BondType prevBondType, int prevAtom) {
            visited[atom] = true;
            
            for (const auto& [neighbor, bondType] : adj[atom]) {
                if (neighbor == prevAtom) continue; // Don't go back
                
                // Check if this bond continues the alternating pattern
                bool isAlternating = false;
                if (prevBondType == RDKit::Bond::SINGLE && bondType == RDKit::Bond::DOUBLE) {
                    isAlternating = true;
                } else if (prevBondType == RDKit::Bond::DOUBLE && bondType == RDKit::Bond::SINGLE) {
                    isAlternating = true;
                }
                
                if (isAlternating) {
                    // Add this bond to the current path
                    pathBondTypes.push_back(bondType);
                    pathBonds.push_back({std::min(atom, neighbor), std::max(atom, neighbor)});
                    
                    // If path length is at least 3, mark all bonds in the path as alternating
                    if (pathBondTypes.size() >= 3) {
                        for (const auto& bond : pathBonds) {
                            alternatingBonds.insert(bond);
                        }
                    }
                    
                    if (!visited[neighbor]) {
                        dfs(neighbor, bondType, atom);
                    }
                    
                    // Backtrack
                    pathBondTypes.pop_back();
                    pathBonds.pop_back();
                } else {
                    // Start a new path with this bond
                    std::vector<RDKit::Bond::BondType> newPath = {bondType};
                    std::vector<std::pair<int, int>> newBonds = {{std::min(atom, neighbor), std::max(atom, neighbor)}};
                    
                    if (!visited[neighbor]) {
                        dfs(neighbor, bondType, atom);
                    }
                }
            }
            
            visited[atom] = false; // Allow revisiting in different paths
        };
        
        // Start DFS with a dummy previous bond type
        dfs(startAtom, RDKit::Bond::UNSPECIFIED, -1);
    }
    
    return static_cast<double>(alternatingBonds.size()) / totalBonds;
}

// 14. AtomValenceSpread
AtomValenceSpread::AtomValenceSpread()
    : Descriptor("AtomValenceSpread", "Max formal valence − Min formal valence among heavy atoms") {}

std::variant<double, int, std::string> AtomValenceSpread::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0;
    
    const RDKit::ROMol& rdkMol = *mol.getMolecule();
    int minValence = std::numeric_limits<int>::max();
    int maxValence = std::numeric_limits<int>::min();
    bool foundHeavyAtom = false;
    
    for (const auto& atom : rdkMol.atoms()) {
        if (atom->getAtomicNum() > 1) { // Heavy atom
            foundHeavyAtom = true;
            int valence = atom->getTotalValence();
            minValence = std::min(minValence, valence);
            maxValence = std::max(maxValence, valence);
        }
    }
    
    if (!foundHeavyAtom) return 0;
    return maxValence - minValence;
}

// 15. RingHydrogenRatio
RingHydrogenRatio::RingHydrogenRatio()
    : Descriptor("RingHydrogenRatio", "Implicit H on ring atoms ÷ total implicit H in the molecule") {}

std::variant<double, int, std::string> RingHydrogenRatio::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    
    const RDKit::ROMol& rdkMol = *mol.getMolecule();
    RDKit::RingInfo* ringInfo = rdkMol.getRingInfo();
    if (!ringInfo->isInitialized()) {
        RDKit::MolOps::findSSSR(const_cast<RDKit::ROMol&>(rdkMol));
    }
    
    int ringImplicitH = 0;
    int totalImplicitH = 0;
    
    for (const auto& atom : rdkMol.atoms()) {
        int implicitH = atom->getImplicitValence();
        totalImplicitH += implicitH;
        
        if (ringInfo->numAtomRings(atom->getIdx()) > 0) {
            ringImplicitH += implicitH;
        }
    }
    
    if (totalImplicitH == 0) return 0.0;
    return static_cast<double>(ringImplicitH) / totalImplicitH;
}

// 16. MaxConsecutiveSingleBonds
MaxConsecutiveSingleBonds::MaxConsecutiveSingleBonds()
    : Descriptor("MaxConsecutiveSingleBonds", "Length of the longest path consisting solely of single bonds") {}

std::variant<double, int, std::string> MaxConsecutiveSingleBonds::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0;
    
    const RDKit::ROMol& rdkMol = *mol.getMolecule();
    int numAtoms = rdkMol.getNumAtoms();
    if (numAtoms <= 1) return 0;
    
    // Build adjacency list for single bonds only
    std::vector<std::vector<int>> singleBondAdj(numAtoms);
    for (const auto& bond : rdkMol.bonds()) {
        if (bond->getBondType() == RDKit::Bond::SINGLE) {
            int u = bond->getBeginAtomIdx();
            int v = bond->getEndAtomIdx();
            singleBondAdj[u].push_back(v);
            singleBondAdj[v].push_back(u);
        }
    }
    
    // Find longest path using BFS from each atom
    int maxPathLength = 0;
    
    for (int startAtom = 0; startAtom < numAtoms; ++startAtom) {
        std::vector<bool> visited(numAtoms, false);
        std::vector<int> distance(numAtoms, 0);
        std::queue<int> q;
        
        q.push(startAtom);
        visited[startAtom] = true;
        
        while (!q.empty()) {
            int u = q.front();
            q.pop();
            
            for (int v : singleBondAdj[u]) {
                if (!visited[v]) {
                    visited[v] = true;
                    distance[v] = distance[u] + 1;
                    q.push(v);
                    
                    maxPathLength = std::max(maxPathLength, distance[v]);
                }
            }
        }
    }
    
    return maxPathLength + 1; // Convert from bond count to atom count in the path
}

// 17. HeteroatomPathFraction
HeteroatomPathFraction::HeteroatomPathFraction()
    : Descriptor("HeteroatomPathFraction", "Fraction of shortest paths (length ≤ 4) whose endpoints are heteroatoms") {}

std::variant<double, int, std::string> HeteroatomPathFraction::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    
    const RDKit::ROMol& rdkMol = *mol.getMolecule();
    int numAtoms = rdkMol.getNumAtoms();
    
    // Compute all shortest paths
    std::vector<std::vector<int>> shortestPaths = computeAllPairsShortestPaths(rdkMol);
    
    int heteroEndpointPaths = 0;
    int totalPaths = 0;
    
    // Count paths with length <= 4 and heteroatom endpoints
    for (int i = 0; i < numAtoms; ++i) {
        for (int j = i + 1; j < numAtoms; ++j) { // Only count each path once
            if (shortestPaths[i][j] <= 4 && shortestPaths[i][j] > 0) { // Path exists and has length <= 4
                totalPaths++;
                
                const RDKit::Atom* atomI = rdkMol.getAtomWithIdx(i);
                const RDKit::Atom* atomJ = rdkMol.getAtomWithIdx(j);
                
                if (isHeteroatom(atomI) && isHeteroatom(atomJ)) {
                    heteroEndpointPaths++;
                }
            }
        }
    }
    
    if (totalPaths == 0) return 0.0;
    return static_cast<double>(heteroEndpointPaths) / totalPaths;
}

// 18. MeanResonanceBondOrder
MeanResonanceBondOrder::MeanResonanceBondOrder()
    : Descriptor("MeanResonanceBondOrder", "Mean bond order across conjugated (resonance-capable) bonds") {}

std::variant<double, int, std::string> MeanResonanceBondOrder::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    
    const RDKit::ROMol& rdkMol = *mol.getMolecule();
    
    // Identify conjugated bonds (part of resonance systems)
    // For this implementation, we'll consider aromatic bonds and bonds between sp2 hybridized atoms
    double totalBondOrder = 0.0;
    int conjugatedBondCount = 0;
    
    for (const auto& bond : rdkMol.bonds()) {
        bool isConjugated = false;
        
        // Aromatic bonds are conjugated
        if (bond->getBondType() == RDKit::Bond::AROMATIC) {
            isConjugated = true;
        }
        // Bonds between sp2 hybridized atoms can be part of conjugated systems
        else if (bond->getBeginAtom()->getHybridization() == RDKit::Atom::SP2 && 
                 bond->getEndAtom()->getHybridization() == RDKit::Atom::SP2) {
            isConjugated = true;
        }
        
        if (isConjugated) {
            double bondOrder = 0.0;
            if (bond->getBondType() == RDKit::Bond::SINGLE) bondOrder = 1.0;
            else if (bond->getBondType() == RDKit::Bond::DOUBLE) bondOrder = 2.0;
            else if (bond->getBondType() == RDKit::Bond::TRIPLE) bondOrder = 3.0;
            else if (bond->getBondType() == RDKit::Bond::AROMATIC) bondOrder = 1.5;
            
            totalBondOrder += bondOrder;
            conjugatedBondCount++;
        }
    }
    
    if (conjugatedBondCount == 0) return 0.0;
    return totalBondOrder / conjugatedBondCount;
}

// 19. CarbonIsotopeCount
CarbonIsotopeCount::CarbonIsotopeCount()
    : Descriptor("CarbonIsotopeCount", "Count of carbon atoms carrying an explicit isotope label in SMILES") {}

std::variant<double, int, std::string> CarbonIsotopeCount::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0;
    
    const RDKit::ROMol& rdkMol = *mol.getMolecule();
    int isotopeCount = 0;
    
    for (const auto& atom : rdkMol.atoms()) {
        if (atom->getAtomicNum() == 6 && atom->getIsotope() > 0) { // Carbon with explicit isotope
            isotopeCount++;
        }
    }
    
    return isotopeCount;
}

// 20. ChargeBalanceSkew
ChargeBalanceSkew::ChargeBalanceSkew()
    : Descriptor("ChargeBalanceSkew", "Skewness of formal charge distribution") {}

std::variant<double, int, std::string> ChargeBalanceSkew::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    
    const RDKit::ROMol& rdkMol = *mol.getMolecule();
    std::vector<double> charges;
    
    for (const auto& atom : rdkMol.atoms()) {
        if (atom->getAtomicNum() > 1) { // Only consider heavy atoms
            charges.push_back(static_cast<double>(atom->getFormalCharge()));
        }
    }
    
    if (charges.size() < 3) return 0.0; // Need at least 3 values for skewness
    return calculateFisherSkewness(charges);
}

// 21. ElementRunLengthVariance
ElementRunLengthVariance::ElementRunLengthVariance()
    : Descriptor("ElementRunLengthVariance", "Variance of consecutive same-element run lengths in canonical SMILES") {}

std::variant<double, int, std::string> ElementRunLengthVariance::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    
    // Get canonical SMILES
    std::string smiles;
    try {
        smiles = RDKit::MolToSmiles(*mol.getMolecule(), true); // canonical=true
    } catch (...) {
        return 0.0;
    }
    
    // Parse SMILES to find runs of same element
    std::vector<int> runLengths;
    char currentElement = '\0';
    int currentRunLength = 0;
    
    for (char c : smiles) {
        // Only consider uppercase letters as start of element symbols
        if (std::isupper(c)) {
            if (currentElement != '\0') {
                runLengths.push_back(currentRunLength);
            }
            currentElement = c;
            currentRunLength = 1;
        }
        // For lowercase letters, check if they're part of the current element (e.g., 'Cl')
        else if (std::islower(c) && currentElement != '\0') {
            // Still the same element, just continue the run
        }
        // For digits, brackets, etc., they're not part of element symbols
        else if (!std::isalpha(c)) {
            if (currentElement != '\0') {
                runLengths.push_back(currentRunLength);
                currentElement = '\0';
                currentRunLength = 0;
            }
        }
    }
    
    // Add the last run if there is one
    if (currentRunLength > 0) {
        runLengths.push_back(currentRunLength);
    }
    
    // Calculate variance of run lengths
    if (runLengths.size() < 2) return 0.0;
    
    double mean = std::accumulate(runLengths.begin(), runLengths.end(), 0.0) / runLengths.size();
    double variance = 0.0;
    
    for (int length : runLengths) {
        double diff = length - mean;
        variance += diff * diff;
    }
    
    variance /= runLengths.size();
    return variance;
}

// 22. SSSRtoRingAtomRatio
SSSRtoRingAtomRatio::SSSRtoRingAtomRatio()
    : Descriptor("SSSRtoRingAtomRatio", "Number of SSSR rings ÷ total number of ring atoms") {}

std::variant<double, int, std::string> SSSRtoRingAtomRatio::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    
    const RDKit::ROMol& rdkMol = *mol.getMolecule();
    RDKit::RingInfo* ringInfo = rdkMol.getRingInfo();
    if (!ringInfo->isInitialized()) {
        RDKit::MolOps::findSSSR(const_cast<RDKit::ROMol&>(rdkMol));
    }
    
    const auto& rings = ringInfo->atomRings();
    if (rings.empty()) return 0.0;
    
    // Count unique ring atoms
    std::set<int> uniqueRingAtoms;
    for (const auto& ring : rings) {
        uniqueRingAtoms.insert(ring.begin(), ring.end());
    }
    
    if (uniqueRingAtoms.empty()) return 0.0;
    return static_cast<double>(rings.size()) / uniqueRingAtoms.size();
}

// 23. MaxRingDistance
MaxRingDistance::MaxRingDistance()
    : Descriptor("MaxRingDistance", "Largest graph distance between any two atoms within the same ring") {}

std::variant<double, int, std::string> MaxRingDistance::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0;
    
    const RDKit::ROMol& rdkMol = *mol.getMolecule();
    RDKit::RingInfo* ringInfo = rdkMol.getRingInfo();
    if (!ringInfo->isInitialized()) {
        RDKit::MolOps::findSSSR(const_cast<RDKit::ROMol&>(rdkMol));
    }
    
    const auto& rings = ringInfo->atomRings();
    if (rings.empty()) return 0;
    
    // Compute all shortest paths
    std::vector<std::vector<int>> shortestPaths = computeAllPairsShortestPaths(rdkMol);
    
    int maxDistance = 0;
    
    // For each ring, find the maximum distance between any two atoms in the ring
    for (const auto& ring : rings) {
        for (size_t i = 0; i < ring.size(); ++i) {
            for (size_t j = i + 1; j < ring.size(); ++j) {
                int distance = shortestPaths[ring[i]][ring[j]];
                maxDistance = std::max(maxDistance, distance);
            }
        }
    }
    
    return maxDistance;
}

// 24. MeanBondStereoStates
MeanBondStereoStates::MeanBondStereoStates()
    : Descriptor("MeanBondStereoStates", "Average stereo code per bond (none = 0, / or \\ = 1, E/Z = 2)") {}

std::variant<double, int, std::string> MeanBondStereoStates::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    
    const RDKit::ROMol& rdkMol = *mol.getMolecule();
    int totalBonds = rdkMol.getNumBonds();
    if (totalBonds == 0) return 0.0;
    
    double totalStereoCode = 0.0;
    
    for (const auto& bond : rdkMol.bonds()) {
        int stereoCode = 0;
        
        // Check for stereochemistry
        if (bond->getStereo() == RDKit::Bond::STEREOANY) {
            stereoCode = 1; // Wiggly bond
        }
        else if (bond->getStereo() == RDKit::Bond::STEREOE ||
                 bond->getStereo() == RDKit::Bond::STEREOZ) {
            stereoCode = 2; // E/Z stereochemistry
        }
        else if (bond->getStereo() == RDKit::Bond::STEREOCIS ||
                 bond->getStereo() == RDKit::Bond::STEREOTRANS) {
            stereoCode = 2; // Cis/Trans stereochemistry
        }
        else if (bond->getBondDir() == RDKit::Bond::BEGINWEDGE ||
                 bond->getBondDir() == RDKit::Bond::BEGINDASH) {
            stereoCode = 1; // Wedge or dash bond
        }
        
        totalStereoCode += stereoCode;
    }
    
    return totalStereoCode / totalBonds;
}

// 25. TerminalDoubleBondCount
TerminalDoubleBondCount::TerminalDoubleBondCount()
    : Descriptor("TerminalDoubleBondCount", "Count of atoms in double bonds that have degree = 1") {}

std::variant<double, int, std::string> TerminalDoubleBondCount::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0;
    
    const RDKit::ROMol& rdkMol = *mol.getMolecule();
    int count = 0;
    
    for (const auto& bond : rdkMol.bonds()) {
        if (bond->getBondType() == RDKit::Bond::DOUBLE) {
            const RDKit::Atom* beginAtom = bond->getBeginAtom();
            const RDKit::Atom* endAtom = bond->getEndAtom();
            
            // Check if either atom has degree 1
            if (beginAtom->getDegree() == 1 || endAtom->getDegree() == 1) {
                count++;
            }
        }
    }
    
    return count;
}

// 26. AromaticSubstituentDiversity
AromaticSubstituentDiversity::AromaticSubstituentDiversity()
    : Descriptor("AromaticSubstituentDiversity", "Shannon entropy of element types directly attached to aromatic atoms but outside rings") {}

std::variant<double, int, std::string> AromaticSubstituentDiversity::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    
    const RDKit::ROMol& rdkMol = *mol.getMolecule();
    RDKit::RingInfo* ringInfo = rdkMol.getRingInfo();
    if (!ringInfo->isInitialized()) {
        RDKit::MolOps::findSSSR(const_cast<RDKit::ROMol&>(rdkMol));
    }
    
    // Map of element type -> frequency for substituents
    std::map<int, int> substituentElementFreq;
    int totalSubstituents = 0;
    
    for (const auto& atom : rdkMol.atoms()) {
        // Check if atom is aromatic and in a ring
        if (atom->getIsAromatic() && ringInfo->numAtomRings(atom->getIdx()) > 0) {
            // Look at its neighbors
            for (const auto& nbr : rdkMol.atomNeighbors(atom)) {
                // If neighbor is not in a ring, it's a substituent
                if (ringInfo->numAtomRings(nbr->getIdx()) == 0) {
                    substituentElementFreq[nbr->getAtomicNum()]++;
                    totalSubstituents++;
                }
            }
        }
    }
    
    if (totalSubstituents == 0) return 0.0;
    return calculateShannonEntropy(substituentElementFreq, totalSubstituents);
}

// 27. HalogenNeighborHybridRatio
HalogenNeighborHybridRatio::HalogenNeighborHybridRatio()
    : Descriptor("HalogenNeighborHybridRatio", "Halogen atoms bound to sp³ carbon ÷ total halogen atoms") {}

std::variant<double, int, std::string> HalogenNeighborHybridRatio::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    
    const RDKit::ROMol& rdkMol = *mol.getMolecule();
    int halogensSp3 = 0;
    int totalHalogens = 0;
    
    for (const auto& atom : rdkMol.atoms()) {
        if (isHalogen(atom)) {
            totalHalogens++;
            
            // Check if connected to sp3 carbon
            for (const auto& nbr : rdkMol.atomNeighbors(atom)) {
                if (nbr->getAtomicNum() == 6 && nbr->getHybridization() == RDKit::Atom::SP3) {
                    halogensSp3++;
                    break;
                }
            }
        }
    }
    
    if (totalHalogens == 0) return 0.0;
    return static_cast<double>(halogensSp3) / totalHalogens;
}

// 28. MeanAtomBetweennessCentrality
MeanAtomBetweennessCentrality::MeanAtomBetweennessCentrality()
    : Descriptor("MeanAtomBetweennessCentrality", "Average betweenness centrality of heavy atoms (un-weighted graph)") {}

std::variant<double, int, std::string> MeanAtomBetweennessCentrality::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    
    const RDKit::ROMol& rdkMol = *mol.getMolecule();
    int numAtoms = rdkMol.getNumAtoms();
    if (numAtoms <= 2) return 0.0; // Betweenness is only meaningful for graphs with at least 3 nodes
    
    // Count heavy atoms
    int numHeavyAtoms = 0;
    for (const auto& atom : rdkMol.atoms()) {
        if (atom->getAtomicNum() > 1) numHeavyAtoms++;
    }
    
    if (numHeavyAtoms <= 2) return 0.0;
    
    // Build adjacency list
    std::vector<std::vector<int>> adj(numAtoms);
    for (const auto& bond : rdkMol.bonds()) {
        int u = bond->getBeginAtomIdx();
        int v = bond->getEndAtomIdx();
        adj[u].push_back(v);
        adj[v].push_back(u);
    }
    
    // Calculate betweenness centrality for each atom
    std::vector<double> betweenness(numAtoms, 0.0);
    
    // For each source vertex
    for (int s = 0; s < numAtoms; ++s) {
        if (rdkMol.getAtomWithIdx(s)->getAtomicNum() <= 1) continue; // Skip hydrogen atoms
        
        // BFS to find shortest paths and count them
        std::vector<std::vector<int>> pred(numAtoms); // predecessors on shortest paths from s
        std::vector<int> dist(numAtoms, -1); // distances from source s
        std::vector<int> sigma(numAtoms, 0); // number of shortest paths from s to each vertex
        std::queue<int> q;
        
        dist[s] = 0;
        sigma[s] = 1;
        q.push(s);
        
        while (!q.empty()) {
            int v = q.front();
            q.pop();
            
            for (int w : adj[v]) {
                // Path discovery
                if (dist[w] < 0) {
                    dist[w] = dist[v] + 1;
                    q.push(w);
                }
                
                // Path counting
                if (dist[w] == dist[v] + 1) {
                    sigma[w] += sigma[v];
                    pred[w].push_back(v);
                }
            }
        }
        
        // Accumulation phase
        std::vector<double> delta(numAtoms, 0.0);
        std::vector<int> stack; // vertices in order of non-increasing distance from s
        
        for (int w = 0; w < numAtoms; ++w) {
            if (w != s && dist[w] > 0) {
                stack.push_back(w);
            }
        }
        
        std::sort(stack.begin(), stack.end(), [&dist](int a, int b) {
            return dist[a] > dist[b];
        });
        
        for (int w : stack) {
            for (int v : pred[w]) {
                delta[v] += (static_cast<double>(sigma[v]) / sigma[w]) * (1.0 + delta[w]);
            }
            if (w != s) {
                betweenness[w] += delta[w];
            }
        }
    }
    
    // Calculate mean betweenness of heavy atoms
    double totalBetweenness = 0.0;
    int heavyAtomCount = 0;
    
    for (int i = 0; i < numAtoms; ++i) {
        if (rdkMol.getAtomWithIdx(i)->getAtomicNum() > 1) { // Heavy atom
            totalBetweenness += betweenness[i];
            heavyAtomCount++;
        }
    }
    
    if (heavyAtomCount == 0) return 0.0;
    return totalBetweenness / heavyAtomCount;
}

// 29. BondOrderEntropy
BondOrderEntropy::BondOrderEntropy()
    : Descriptor("BondOrderEntropy", "Shannon entropy of bond-order frequencies") {}

std::variant<double, int, std::string> BondOrderEntropy::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    
    const RDKit::ROMol& rdkMol = *mol.getMolecule();
    int totalBonds = rdkMol.getNumBonds();
    if (totalBonds == 0) return 0.0;
    
    // Count bond orders
    std::map<int, int> bondOrderFreq;
    
    for (const auto& bond : rdkMol.bonds()) {
        int bondOrder = 0;
        if (bond->getBondType() == RDKit::Bond::SINGLE) bondOrder = 1;
        else if (bond->getBondType() == RDKit::Bond::DOUBLE) bondOrder = 2;
        else if (bond->getBondType() == RDKit::Bond::TRIPLE) bondOrder = 3;
        else if (bond->getBondType() == RDKit::Bond::AROMATIC) bondOrder = 4; // Use 4 for aromatic to distinguish from others
        
        bondOrderFreq[bondOrder]++;
    }
    
    // Calculate Shannon entropy
    double entropy = 0.0;
    for (const auto& pair : bondOrderFreq) {
        double p = static_cast<double>(pair.second) / totalBonds;
        entropy -= p * std::log2(p);
    }
    
    return entropy;
}

// 30. HybridizationSymmetryIndex
HybridizationSymmetryIndex::HybridizationSymmetryIndex()
    : Descriptor("HybridizationSymmetryIndex", "Symmetry measure of hybridization distribution") {}

std::variant<double, int, std::string> HybridizationSymmetryIndex::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    
    const RDKit::ROMol& rdkMol = *mol.getMolecule();
    
    // Count hybridization types
    std::map<RDKit::Atom::HybridizationType, int> hybridizationFreq;
    int totalHeavyAtoms = 0;
    
    for (const auto& atom : rdkMol.atoms()) {
        if (atom->getAtomicNum() > 1) { // Only consider heavy atoms
            hybridizationFreq[atom->getHybridization()]++;
            totalHeavyAtoms++;
        }
    }
    
    if (totalHeavyAtoms <= 1) return 1.0; // Perfect symmetry with only one atom
    
    // Calculate symmetry index (1 - Gini coefficient of hybridization distribution)
    std::vector<double> frequencies;
    for (const auto& pair : hybridizationFreq) {
        frequencies.push_back(static_cast<double>(pair.second));
    }
    
    double gini = calculateGiniCoefficient(frequencies);
    return 1.0 - gini; // Higher value means more symmetric distribution
}

// 31. RingSizeMedian
RingSizeMedian::RingSizeMedian()
    : Descriptor("RingSizeMedian", "Median ring size among all SSSR rings (0 if no rings)") {}

std::variant<double, int, std::string> RingSizeMedian::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0;
    
    const RDKit::ROMol& rdkMol = *mol.getMolecule();
    RDKit::RingInfo* ringInfo = rdkMol.getRingInfo();
    if (!ringInfo->isInitialized()) {
        RDKit::MolOps::findSSSR(const_cast<RDKit::ROMol&>(rdkMol));
    }
    
    const auto& rings = ringInfo->atomRings();
    if (rings.empty()) return 0;
    
    // Collect ring sizes
    std::vector<int> ringSizes;
    for (const auto& ring : rings) {
        ringSizes.push_back(ring.size());
    }
    
    // Calculate median
    std::sort(ringSizes.begin(), ringSizes.end());
    if (ringSizes.size() % 2 == 0) {
        // Even number of rings, average the middle two
        int mid1 = ringSizes[ringSizes.size() / 2 - 1];
        int mid2 = ringSizes[ringSizes.size() / 2];
        return (mid1 + mid2) / 2;
    } else {
        // Odd number of rings, return the middle one
        return ringSizes[ringSizes.size() / 2];
    }
}

// 32. LongestPathBondOrderProduct
LongestPathBondOrderProduct::LongestPathBondOrderProduct()
    : Descriptor("LongestPathBondOrderProduct", "Product of bond orders along the molecule's longest simple path") {}

std::variant<double, int, std::string> LongestPathBondOrderProduct::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    
    const RDKit::ROMol& rdkMol = *mol.getMolecule();
    int numAtoms = rdkMol.getNumAtoms();
    if (numAtoms <= 1) return 0.0;
    
    // Build adjacency list with bond indices
    std::vector<std::vector<std::pair<int, int>>> adj(numAtoms); // (neighbor, bond index)
    for (const auto& bond : rdkMol.bonds()) {
        int u = bond->getBeginAtomIdx();
        int v = bond->getEndAtomIdx();
        int bondIdx = bond->getIdx();
        
        adj[u].push_back({v, bondIdx});
        adj[v].push_back({u, bondIdx});
    }
    
    // Find longest path using BFS from each atom
    std::vector<int> longestPath;
    int maxPathLength = 0;
    
    for (int startAtom = 0; startAtom < numAtoms; ++startAtom) {
        std::vector<bool> visited(numAtoms, false);
        std::vector<int> distance(numAtoms, 0);
        std::vector<int> parent(numAtoms, -1);
        std::vector<int> parentBond(numAtoms, -1);
        std::queue<int> q;
        
        q.push(startAtom);
        visited[startAtom] = true;
        
        int farthestNode = startAtom;
        int maxDist = 0;
        
        while (!q.empty()) {
            int u = q.front();
            q.pop();
            
            for (const auto& [v, bondIdx] : adj[u]) {
                if (!visited[v]) {
                    visited[v] = true;
                    distance[v] = distance[u] + 1;
                    parent[v] = u;
                    parentBond[v] = bondIdx;
                    q.push(v);
                    
                    if (distance[v] > maxDist) {
                        maxDist = distance[v];
                        farthestNode = v;
                    }
                }
            }
        }
        
        // If this path is longer than the current longest path, update it
        if (maxDist > maxPathLength) {
            maxPathLength = maxDist;
            
            // Reconstruct the path
            std::vector<int> path;
            std::vector<int> bonds;
            int current = farthestNode;
            
            while (current != -1) {
                path.push_back(current);
                if (parentBond[current] != -1) {
                    bonds.push_back(parentBond[current]);
                }
                current = parent[current];
            }
            
            longestPath = path;
        }
    }
    
    if (longestPath.size() <= 1) return 0.0;
    
    // Calculate product of bond orders along the path
    double product = 1.0;
    for (int i = 0; i < longestPath.size() - 1; ++i) {
        int u = longestPath[i];
        int v = longestPath[i + 1];
        
        // Find the bond between u and v
        const RDKit::Bond* bond = rdkMol.getBondBetweenAtoms(u, v);
        if (!bond) continue;
        
        double bondOrder = 0.0;
        if (bond->getBondType() == RDKit::Bond::SINGLE) bondOrder = 1.0;
        else if (bond->getBondType() == RDKit::Bond::DOUBLE) bondOrder = 2.0;
        else if (bond->getBondType() == RDKit::Bond::TRIPLE) bondOrder = 3.0;
        else if (bond->getBondType() == RDKit::Bond::AROMATIC) bondOrder = 1.5;
        
        product *= bondOrder;
    }
    
    return product;
}

// 33. HeteroatomDegreeVariance
HeteroatomDegreeVariance::HeteroatomDegreeVariance()
    : Descriptor("HeteroatomDegreeVariance", "Variance of degrees among heteroatoms") {}

std::variant<double, int, std::string> HeteroatomDegreeVariance::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    
    const RDKit::ROMol& rdkMol = *mol.getMolecule();
    std::vector<int> heteroDegrees;
    
    for (const auto& atom : rdkMol.atoms()) {
        if (isHeteroatom(atom)) {
            heteroDegrees.push_back(atom->getDegree());
        }
    }
    
    if (heteroDegrees.size() < 2) return 0.0; // Need at least 2 heteroatoms for variance
    
    // Calculate variance
    double mean = std::accumulate(heteroDegrees.begin(), heteroDegrees.end(), 0.0) / heteroDegrees.size();
    double variance = 0.0;
    
    for (int degree : heteroDegrees) {
        double diff = degree - mean;
        variance += diff * diff;
    }
    
    variance /= heteroDegrees.size();
    return variance;
}

// 34. NonCarbonAtomAdjacencyCount
NonCarbonAtomAdjacencyCount::NonCarbonAtomAdjacencyCount()
    : Descriptor("NonCarbonAtomAdjacencyCount", "Number of bonds in which neither atom is carbon") {}

std::variant<double, int, std::string> NonCarbonAtomAdjacencyCount::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0;
    
    const RDKit::ROMol& rdkMol = *mol.getMolecule();
    int count = 0;
    
    for (const auto& bond : rdkMol.bonds()) {
        const RDKit::Atom* beginAtom = bond->getBeginAtom();
        const RDKit::Atom* endAtom = bond->getEndAtom();
        
        if (beginAtom->getAtomicNum() != 6 && endAtom->getAtomicNum() != 6 && 
            beginAtom->getAtomicNum() > 1 && endAtom->getAtomicNum() > 1) {
            // Neither atom is carbon, and both are heavy atoms
            count++;
        }
    }
    
    return count;
}

// 35. MaxBondENGap
MaxBondENGap::MaxBondENGap()
    : Descriptor("MaxBondENGap", "Maximum absolute Pauling electronegativity difference across any single bond") {}

std::variant<double, int, std::string> MaxBondENGap::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    
    const RDKit::ROMol& rdkMol = *mol.getMolecule();
    double maxENDiff = 0.0;
    
    for (const auto& bond : rdkMol.bonds()) {
        double en1 = getAtomENPauling(bond->getBeginAtom());
        double en2 = getAtomENPauling(bond->getEndAtom());
        double enDiff = std::abs(en1 - en2);
        
        maxENDiff = std::max(maxENDiff, enDiff);
    }
    
    return maxENDiff;
}

// 36. RingIndexSum
RingIndexSum::RingIndexSum()
    : Descriptor("RingIndexSum", "Σ (1 / ring size) over all SSSR rings") {}

std::variant<double, int, std::string> RingIndexSum::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    
    const RDKit::ROMol& rdkMol = *mol.getMolecule();
    RDKit::RingInfo* ringInfo = rdkMol.getRingInfo();
    if (!ringInfo->isInitialized()) {
        RDKit::MolOps::findSSSR(const_cast<RDKit::ROMol&>(rdkMol));
    }
    
    const auto& rings = ringInfo->atomRings();
    if (rings.empty()) return 0.0;
    
    double sum = 0.0;
    for (const auto& ring : rings) {
        sum += 1.0 / ring.size();
    }
    
    return sum;
}

// 37. BranchDepthAverage
BranchDepthAverage::BranchDepthAverage()
    : Descriptor("BranchDepthAverage", "Mean shortest-path distance from each branch point to the nearest terminal atom") {}

std::variant<double, int, std::string> BranchDepthAverage::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    
    const RDKit::ROMol& rdkMol = *mol.getMolecule();
    int numAtoms = rdkMol.getNumAtoms();
    
    // Identify branch points and terminal atoms
    std::vector<int> branchPoints;
    std::vector<int> terminalAtoms;
    
    for (const auto& atom : rdkMol.atoms()) {
        if (atom->getAtomicNum() <= 1) continue; // Skip hydrogen atoms
        
        int degree = atom->getDegree();
        if (degree > 2) {
            branchPoints.push_back(atom->getIdx());
        } else if (degree == 1) {
            terminalAtoms.push_back(atom->getIdx());
        }
    }
    
    if (branchPoints.empty() || terminalAtoms.empty()) return 0.0;
    
    // Compute all shortest paths
    std::vector<std::vector<int>> shortestPaths = computeAllPairsShortestPaths(rdkMol);
    
    // For each branch point, find distance to nearest terminal atom
    double totalDepth = 0.0;
    
    for (int branchPoint : branchPoints) {
        int minDistance = std::numeric_limits<int>::max();
        
        for (int terminalAtom : terminalAtoms) {
            minDistance = std::min(minDistance, shortestPaths[branchPoint][terminalAtom]);
        }
        
        totalDepth += minDistance;
    }
    
    return totalDepth / branchPoints.size();
}

// 38. ConsecutiveHeteroBondFraction
ConsecutiveHeteroBondFraction::ConsecutiveHeteroBondFraction()
    : Descriptor("ConsecutiveHeteroBondFraction", "Bonds connecting two heteroatoms ÷ total bonds") {}

std::variant<double, int, std::string> ConsecutiveHeteroBondFraction::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    
    const RDKit::ROMol& rdkMol = *mol.getMolecule();
    int totalBonds = rdkMol.getNumBonds();
    if (totalBonds == 0) return 0.0;
    
    int heteroBonds = 0;
    
    for (const auto& bond : rdkMol.bonds()) {
        const RDKit::Atom* beginAtom = bond->getBeginAtom();
        const RDKit::Atom* endAtom = bond->getEndAtom();
        
        if (isHeteroatom(beginAtom) && isHeteroatom(endAtom)) {
            heteroBonds++;
        }
    }
    
    return static_cast<double>(heteroBonds) / totalBonds;
}

// 39. HeavyAtomHydrogenRatio
HeavyAtomHydrogenRatio::HeavyAtomHydrogenRatio()
    : Descriptor("HeavyAtomHydrogenRatio", "(implicit H + explicit H) ÷ heavy atoms") {}

std::variant<double, int, std::string> HeavyAtomHydrogenRatio::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    
    const RDKit::ROMol& rdkMol = *mol.getMolecule();
    int heavyAtoms = 0;
    int totalHydrogens = 0;
    
    for (const auto& atom : rdkMol.atoms()) {
        if (atom->getAtomicNum() > 1) { // Heavy atom
            heavyAtoms++;
            totalHydrogens += atom->getImplicitValence(); // Implicit hydrogens
            totalHydrogens += atom->getNumExplicitHs(); // Explicit hydrogens
        } else if (atom->getAtomicNum() == 1) { // Explicit hydrogen atom
            totalHydrogens++;
        }
    }
    
    if (heavyAtoms == 0) return 0.0;
    return static_cast<double>(totalHydrogens) / heavyAtoms;
}

// 40. MeanAtomicMassPerDegree
MeanAtomicMassPerDegree::MeanAtomicMassPerDegree()
    : Descriptor("MeanAtomicMassPerDegree", "Mean (atomic mass ÷ degree) over heavy atoms") {}

std::variant<double, int, std::string> MeanAtomicMassPerDegree::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    
    const RDKit::ROMol& rdkMol = *mol.getMolecule();
    double totalRatio = 0.0;
    int heavyAtomCount = 0;
    
    for (const auto& atom : rdkMol.atoms()) {
        if (atom->getAtomicNum() > 1) { // Heavy atom
            int degree = atom->getDegree();
            if (degree > 0) { // Avoid division by zero
                double atomicMass = atom->getMass();
                totalRatio += atomicMass / degree;
                heavyAtomCount++;
            }
        }
    }
    
    if (heavyAtomCount == 0) return 0.0;
    return totalRatio / heavyAtomCount;
}

// 41. RareElementCount
RareElementCount::RareElementCount()
    : Descriptor("RareElementCount", "Count of atoms with atomic number > 17") {}

std::variant<double, int, std::string> RareElementCount::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0;
    
    const RDKit::ROMol& rdkMol = *mol.getMolecule();
    int count = 0;
    
    for (const auto& atom : rdkMol.atoms()) {
        if (atom->getAtomicNum() > 17) { // Elements heavier than Cl
            count++;
        }
    }
    
    return count;
}

// 42. TopoDistanceSkewness
TopoDistanceSkewness::TopoDistanceSkewness()
    : Descriptor("TopoDistanceSkewness", "Fisher skewness of all-pairs shortest-path length distribution") {}

std::variant<double, int, std::string> TopoDistanceSkewness::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    
    const RDKit::ROMol& rdkMol = *mol.getMolecule();
    int numAtoms = rdkMol.getNumAtoms();
    if (numAtoms <= 2) return 0.0; // Need at least 3 atoms for skewness
    
    // Compute all shortest paths
    std::vector<std::vector<int>> shortestPaths = computeAllPairsShortestPaths(rdkMol);
    
    // Collect all path lengths
    std::vector<double> pathLengths;
    for (int i = 0; i < numAtoms; ++i) {
        for (int j = i + 1; j < numAtoms; ++j) { // Only count each path once
            if (shortestPaths[i][j] < std::numeric_limits<int>::max() / 2) { // Valid path
                pathLengths.push_back(static_cast<double>(shortestPaths[i][j]));
            }
        }
    }
    
    if (pathLengths.size() < 3) return 0.0; // Need at least 3 paths for skewness
    return calculateFisherSkewness(pathLengths);
}

// 43. PercentageIsolatedRingSystems
PercentageIsolatedRingSystems::PercentageIsolatedRingSystems()
    : Descriptor("PercentageIsolatedRingSystems", "Isolated (non-fused) ring systems ÷ total ring systems (0 when no rings)") {}

std::variant<double, int, std::string> PercentageIsolatedRingSystems::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    
    const RDKit::ROMol& rdkMol = *mol.getMolecule();
    RDKit::RingInfo* ringInfo = rdkMol.getRingInfo();
    if (!ringInfo->isInitialized()) {
        RDKit::MolOps::findSSSR(const_cast<RDKit::ROMol&>(rdkMol));
    }
    
    const auto& rings = ringInfo->atomRings();
    if (rings.empty()) return 0.0;
    
    // Build a graph where nodes are rings and edges represent ring fusion
    int numRings = rings.size();
    std::vector<std::vector<int>> ringAdjList(numRings);
    
    for (int i = 0; i < numRings; ++i) {
        const auto& ringI = rings[i];
        std::set<int> atomsInRingI(ringI.begin(), ringI.end());
        
        for (int j = i + 1; j < numRings; ++j) {
            const auto& ringJ = rings[j];
            
            // Check if rings share any atoms (fusion)
            bool fused = false;
            for (int atom : ringJ) {
                if (atomsInRingI.find(atom) != atomsInRingI.end()) {
                    fused = true;
                    break;
                }
            }
            
            if (fused) {
                ringAdjList[i].push_back(j);
                ringAdjList[j].push_back(i);
            }
        }
    }
    
    // Count connected components (ring systems) and isolated rings
    std::vector<bool> visited(numRings, false);
    int totalRingSystems = 0;
    int isolatedRingSystems = 0;
    
    for (int i = 0; i < numRings; ++i) {
        if (!visited[i]) {
            totalRingSystems++;
            bool isIsolated = true;
            
            // DFS to mark all rings in this system
            std::stack<int> stack;
            stack.push(i);
            visited[i] = true;
            
            while (!stack.empty()) {
                int ring = stack.top();
                stack.pop();
                
                for (int adjRing : ringAdjList[ring]) {
                    isIsolated = false; // If there are adjacent rings, it's not isolated
                    if (!visited[adjRing]) {
                        visited[adjRing] = true;
                        stack.push(adjRing);
                    }
                }
            }
            
            if (isIsolated) {
                isolatedRingSystems++;
            }
        }
    }
    
    return static_cast<double>(isolatedRingSystems) / totalRingSystems;
}

// 44. ValenceMismatchBondCount
ValenceMismatchBondCount::ValenceMismatchBondCount()
    : Descriptor("ValenceMismatchBondCount", "Bonds whose two atoms differ in allowed valence by ≥ 3") {}

std::variant<double, int, std::string> ValenceMismatchBondCount::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0;
    
    const RDKit::ROMol& rdkMol = *mol.getMolecule();
    int count = 0;
    
    for (const auto& bond : rdkMol.bonds()) {
        const RDKit::Atom* beginAtom = bond->getBeginAtom();
        const RDKit::Atom* endAtom = bond->getEndAtom();
        
        // Get the default valence for each element
        int beginValence = RDKit::PeriodicTable::getTable()->getDefaultValence(beginAtom->getAtomicNum());
        int endValence = RDKit::PeriodicTable::getTable()->getDefaultValence(endAtom->getAtomicNum());
        
        // Check if valence difference is >= 3
        if (std::abs(beginValence - endValence) >= 3) {
            count++;
        }
    }
    
    return count;
}

// 45. HeteroAtomSequenceLengthMax
HeteroAtomSequenceLengthMax::HeteroAtomSequenceLengthMax()
    : Descriptor("HeteroAtomSequenceLengthMax", "Longest run of consecutive heteroatom symbols in canonical SMILES") {}

std::variant<double, int, std::string> HeteroAtomSequenceLengthMax::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0;
    
    // Get canonical SMILES
    std::string smiles;
    try {
        smiles = RDKit::MolToSmiles(*mol.getMolecule(), true); // canonical=true
    } catch (...) {
        return 0;
    }
    
    // Parse SMILES to find runs of heteroatoms
    int maxRun = 0;
    int currentRun = 0;
    bool inHeteroRun = false;
    
    for (size_t i = 0; i < smiles.size(); ++i) {
        char c = smiles[i];
        
        // Check if this is the start of an element symbol (uppercase letter)
        if (std::isupper(c)) {
            // Check if it's a heteroatom (not C)
            bool isHetero = (c != 'C');
            
            // Handle two-letter elements (e.g., 'Cl', 'Br')
            if (i + 1 < smiles.size() && std::islower(smiles[i + 1])) {
                // Still a heteroatom if it's not 'C' followed by something
                isHetero = (c != 'C');
            }
            
            if (isHetero) {
                if (!inHeteroRun) {
                    inHeteroRun = true;
                    currentRun = 1;
                } else {
                    currentRun++;
                }
                maxRun = std::max(maxRun, currentRun);
            } else {
                inHeteroRun = false;
                currentRun = 0;
            }
        }
        // For non-element characters (digits, brackets, etc.), don't reset the run
        // as they might be part of the element description in SMILES
    }
    
    return maxRun;
}

// 46. ChiralAtomSkew
ChiralAtomSkew::ChiralAtomSkew()
    : Descriptor("ChiralAtomSkew", "(R centers − S centers) ÷ total chiral centers (0 if none)") {}

std::variant<double, int, std::string> ChiralAtomSkew::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    
    const RDKit::ROMol& rdkMol = *mol.getMolecule();
    int rCenters = 0;
    int sCenters = 0;
    
    for (const auto& atom : rdkMol.atoms()) {
        if (atom->getChiralTag() == RDKit::Atom::CHI_TETRAHEDRAL_CW) {
            rCenters++;
        } else if (atom->getChiralTag() == RDKit::Atom::CHI_TETRAHEDRAL_CCW) {
            sCenters++;
        }
    }
    
    int totalChiralCenters = rCenters + sCenters;
    if (totalChiralCenters == 0) return 0.0;
    
    return static_cast<double>(rCenters - sCenters) / totalChiralCenters;
}

// 47. RingAtomUnsaturationRatio
RingAtomUnsaturationRatio::RingAtomUnsaturationRatio()
    : Descriptor("RingAtomUnsaturationRatio", "Unsaturated ring atoms (sp²/sp) ÷ ring atoms") {}

std::variant<double, int, std::string> RingAtomUnsaturationRatio::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    
    const RDKit::ROMol& rdkMol = *mol.getMolecule();
    RDKit::RingInfo* ringInfo = rdkMol.getRingInfo();
    if (!ringInfo->isInitialized()) {
        RDKit::MolOps::findSSSR(const_cast<RDKit::ROMol&>(rdkMol));
    }
    
    int unsaturatedRingAtoms = 0;
    int totalRingAtoms = 0;
    
    for (const auto& atom : rdkMol.atoms()) {
        if (ringInfo->numAtomRings(atom->getIdx()) > 0) {
            totalRingAtoms++;
            
            // Check if atom is unsaturated (sp2 or sp hybridization)
            if (atom->getHybridization() == RDKit::Atom::SP2 || 
                atom->getHybridization() == RDKit::Atom::SP) {
                unsaturatedRingAtoms++;
            }
        }
    }
    
    if (totalRingAtoms == 0) return 0.0;
    return static_cast<double>(unsaturatedRingAtoms) / totalRingAtoms;
}

// 48. TerminalHeteroBondOrderSum
TerminalHeteroBondOrderSum::TerminalHeteroBondOrderSum()
    : Descriptor("TerminalHeteroBondOrderSum", "Sum of bond orders for bonds involving terminal heteroatoms") {}

std::variant<double, int, std::string> TerminalHeteroBondOrderSum::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    
    const RDKit::ROMol& rdkMol = *mol.getMolecule();
    double sum = 0.0;
    
    for (const auto& atom : rdkMol.atoms()) {
        // Check if atom is a terminal heteroatom
        if (isHeteroatom(atom) && atom->getDegree() == 1) {
            // Find the bond to this terminal atom
            for (const auto& bond : rdkMol.atomBonds(atom)) {
                double bondOrder = 0.0;
                if (bond->getBondType() == RDKit::Bond::SINGLE) bondOrder = 1.0;
                else if (bond->getBondType() == RDKit::Bond::DOUBLE) bondOrder = 2.0;
                else if (bond->getBondType() == RDKit::Bond::TRIPLE) bondOrder = 3.0;
                else if (bond->getBondType() == RDKit::Bond::AROMATIC) bondOrder = 1.5;
                
                sum += bondOrder;
            }
        }
    }
    
    return sum;
}

// 49. AromaticNonAromaticBondCount
AromaticNonAromaticBondCount::AromaticNonAromaticBondCount()
    : Descriptor("AromaticNonAromaticBondCount", "Count of bonds that connect an aromatic atom to a non-aromatic heavy atom") {}

std::variant<double, int, std::string> AromaticNonAromaticBondCount::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0;
    
    const RDKit::ROMol& rdkMol = *mol.getMolecule();
    int count = 0;
    
    for (const auto& bond : rdkMol.bonds()) {
        const RDKit::Atom* beginAtom = bond->getBeginAtom();
        const RDKit::Atom* endAtom = bond->getEndAtom();
        
        // Check if one atom is aromatic and the other is not, and both are heavy atoms
        if (beginAtom->getAtomicNum() > 1 && endAtom->getAtomicNum() > 1) {
            if ((beginAtom->getIsAromatic() && !endAtom->getIsAromatic()) ||
                (!beginAtom->getIsAromatic() && endAtom->getIsAromatic())) {
                count++;
            }
        }
    }
    
    return count;
}

// 50. RingAtomChargeVariance
RingAtomChargeVariance::RingAtomChargeVariance()
    : Descriptor("RingAtomChargeVariance", "Variance of formal charges restricted to ring atoms (0 if no rings or no charges)") {}

std::variant<double, int, std::string> RingAtomChargeVariance::calculate(const ::desfact::Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    
    const RDKit::ROMol& rdkMol = *mol.getMolecule();
    RDKit::RingInfo* ringInfo = rdkMol.getRingInfo();
    if (!ringInfo->isInitialized()) {
        RDKit::MolOps::findSSSR(const_cast<RDKit::ROMol&>(rdkMol));
    }
    
    std::vector<double> ringAtomCharges;
    
    for (const auto& atom : rdkMol.atoms()) {
        if (ringInfo->numAtomRings(atom->getIdx()) > 0) {
            ringAtomCharges.push_back(static_cast<double>(atom->getFormalCharge()));
        }
    }
    
    if (ringAtomCharges.empty()) return 0.0;
    
    // Calculate variance
    double mean = std::accumulate(ringAtomCharges.begin(), ringAtomCharges.end(), 0.0) / ringAtomCharges.size();
    double variance = 0.0;
    
    for (double charge : ringAtomCharges) {
        double diff = charge - mean;
        variance += diff * diff;
    }
    
    variance /= ringAtomCharges.size();
    return variance;
}

} // namespace descriptors
} // namespace desfact