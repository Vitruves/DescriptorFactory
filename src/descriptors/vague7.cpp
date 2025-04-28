#include "descriptors/vague7.hpp"
#include <GraphMol/GraphMol.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/Descriptors/MolDescriptors.h>
#include <GraphMol/Descriptors/Crippen.h>
#include <GraphMol/Descriptors/MolSurf.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <cmath>
#include <numeric>
#include <map>
#include <set>

namespace desfact {
namespace descriptors {

// 1. HeteroatomFractionOfRings
HeteroatomFractionOfRings::HeteroatomFractionOfRings()
    : Descriptor("HeteroatomFractionOfRings", "Fraction of ring atoms that are heteroatoms")
{
}

std::variant<double, int, std::string> HeteroatomFractionOfRings::calculate(const ::desfact::Molecule& mol) const
{
    if (!mol.isValid() || !mol.getMolecule()) {
        return 0.0;
    }
    
    const auto rdmol = mol.getMolecule().get();
    
    // Get ring information
    std::vector<std::vector<int>> sssr;
    RDKit::MolOps::findSSSR(*rdmol, sssr);
    
    if (sssr.empty()) {
        return 0.0; // No rings
    }
    
    int ringAtomCount = 0;
    int heteroRingAtomCount = 0;
    
    std::set<int> ringAtoms;
    for (const auto& ring : sssr) {
        for (int atomIdx : ring) {
            if (ringAtoms.insert(atomIdx).second) {
                ringAtomCount++;
                const RDKit::Atom* atom = rdmol->getAtomWithIdx(atomIdx);
                if (atom->getAtomicNum() != 6) { // Not carbon
                    heteroRingAtomCount++;
                }
            }
        }
    }
    
    return ringAtomCount > 0 ? static_cast<double>(heteroRingAtomCount) / ringAtomCount : 0.0;
}

// 2. AromaticRingRatio
AromaticRingRatio::AromaticRingRatio()
    : Descriptor("AromaticRingRatio", "Ratio of aromatic rings to total rings")
{
}

std::variant<double, int, std::string> AromaticRingRatio::calculate(const ::desfact::Molecule& mol) const
{
    if (!mol.isValid() || !mol.getMolecule()) {
        return 0.0;
    }
    
    const auto rdmol = mol.getMolecule().get();
    
    // Get ring information
    std::vector<std::vector<int>> sssr;
    RDKit::MolOps::findSSSR(*rdmol, sssr);
    
    if (sssr.empty()) {
        return 0.0; // No rings
    }
    
    int aromaticRingCount = 0;
    
    for (const auto& ring : sssr) {
        bool isAromatic = true;
        for (int atomIdx : ring) {
            const RDKit::Atom* atom = rdmol->getAtomWithIdx(atomIdx);
            if (!atom->getIsAromatic()) {
                isAromatic = false;
                break;
            }
        }
        if (isAromatic) {
            aromaticRingCount++;
        }
    }
    
    return static_cast<double>(aromaticRingCount) / sssr.size();
}

// 3. RingComplexityIndex
RingComplexityIndex::RingComplexityIndex()
    : Descriptor("RingComplexityIndex", "Measure of ring complexity based on size and connectivity")
{
}

std::variant<double, int, std::string> RingComplexityIndex::calculate(const ::desfact::Molecule& mol) const
{
    if (!mol.isValid() || !mol.getMolecule()) {
        return 0.0;
    }
    
    const auto rdmol = mol.getMolecule().get();
    
    // Get ring information
    std::vector<std::vector<int>> sssr;
    RDKit::MolOps::findSSSR(*rdmol, sssr);
    
    if (sssr.empty()) {
        return 0.0; // No rings
    }
    
    double complexity = 0.0;
    
    // Calculate complexity based on ring sizes and fusion points
    std::map<int, int> atomRingCount; // Count how many rings each atom is part of
    
    for (const auto& ring : sssr) {
        // Add complexity based on ring size (larger rings are more complex)
        complexity += std::log(ring.size());
        
        // Count ring membership for each atom
        for (int atomIdx : ring) {
            atomRingCount[atomIdx]++;
        }
    }
    
    // Add complexity for fusion points (atoms shared by multiple rings)
    for (const auto& [atomIdx, count] : atomRingCount) {
        if (count > 1) {
            complexity += (count - 1) * 0.5; // Each fusion point adds complexity
        }
    }
    
    return complexity;
}

// 4. ElectronegativeAtomDensity
ElectronegativeAtomDensity::ElectronegativeAtomDensity()
    : Descriptor("ElectronegativeAtomDensity", "Density of electronegative atoms (N, O, F, Cl, Br, I)")
{
}

std::variant<double, int, std::string> ElectronegativeAtomDensity::calculate(const ::desfact::Molecule& mol) const
{
    if (!mol.isValid() || !mol.getMolecule()) {
        return 0.0;
    }
    
    const auto rdmol = mol.getMolecule().get();
    int totalAtoms = rdmol->getNumHeavyAtoms();
    
    if (totalAtoms == 0) {
        return 0.0;
    }
    
    int electronegativeCount = 0;
    
    // Count electronegative atoms (N, O, F, Cl, Br, I)
    for (const RDKit::Atom* atom : rdmol->atoms()) {
        int atomicNum = atom->getAtomicNum();
        if (atomicNum == 7 || atomicNum == 8 || atomicNum == 9 || 
            atomicNum == 17 || atomicNum == 35 || atomicNum == 53) {
            electronegativeCount++;
        }
    }
    
    return static_cast<double>(electronegativeCount) / totalAtoms;
}

// 5. ChargeAsymmetryIndex
ChargeAsymmetryIndex::ChargeAsymmetryIndex()
    : Descriptor("ChargeAsymmetryIndex", "Measure of charge distribution asymmetry")
{
}

std::variant<double, int, std::string> ChargeAsymmetryIndex::calculate(const ::desfact::Molecule& mol) const
{
    if (!mol.isValid() || !mol.getMolecule()) {
        return 0.0;
    }
    
    const auto rdmol = mol.getMolecule().get();
    int totalAtoms = rdmol->getNumHeavyAtoms();
    
    if (totalAtoms == 0) {
        return 0.0;
    }
    
    int positiveCharges = 0;
    int negativeCharges = 0;
    
    // Count positive and negative formal charges
    for (const RDKit::Atom* atom : rdmol->atoms()) {
        int charge = atom->getFormalCharge();
        if (charge > 0) {
            positiveCharges += charge;
        } else if (charge < 0) {
            negativeCharges -= charge; // Make positive for calculation
        }
    }
    
    int totalCharges = positiveCharges + negativeCharges;
    if (totalCharges == 0) {
        return 0.0; // No charges
    }
    
    // Calculate asymmetry as the absolute difference between positive and negative charges
    // normalized by the total charge magnitude
    return std::abs(positiveCharges - negativeCharges) / static_cast<double>(totalCharges);
}

// 6. TopologicalPolarSurfaceEfficiency
TopologicalPolarSurfaceEfficiency::TopologicalPolarSurfaceEfficiency()
    : Descriptor("TopologicalPolarSurfaceEfficiency", "TPSA normalized by molecular weight")
{
}

std::variant<double, int, std::string> TopologicalPolarSurfaceEfficiency::calculate(const ::desfact::Molecule& mol) const
{
    if (!mol.isValid() || !mol.getMolecule()) {
        return 0.0;
    }
    
    const auto rdmol = mol.getMolecule().get();
    
    // Calculate TPSA
    double tpsa = RDKit::Descriptors::calcTPSA(*rdmol);
    
    // Calculate molecular weight
    double mw = RDKit::Descriptors::calcAMW(*rdmol);
    
    if (mw == 0) {
        return 0.0;
    }
    
    // Return TPSA per unit of molecular weight
    return tpsa / mw;
}

// 7. HeterocyclicRingCount
HeterocyclicRingCount::HeterocyclicRingCount()
    : Descriptor("HeterocyclicRingCount", "Count of rings containing at least one heteroatom")
{
}

std::variant<double, int, std::string> HeterocyclicRingCount::calculate(const ::desfact::Molecule& mol) const
{
    if (!mol.isValid() || !mol.getMolecule()) {
        return 0;
    }
    
    const auto rdmol = mol.getMolecule().get();
    
    // Get ring information
    std::vector<std::vector<int>> sssr;
    RDKit::MolOps::findSSSR(*rdmol, sssr);
    
    int heterocyclicCount = 0;
    
    for (const auto& ring : sssr) {
        bool hasHeteroatom = false;
        for (int atomIdx : ring) {
            const RDKit::Atom* atom = rdmol->getAtomWithIdx(atomIdx);
            if (atom->getAtomicNum() != 6) { // Not carbon
                hasHeteroatom = true;
                break;
            }
        }
        if (hasHeteroatom) {
            heterocyclicCount++;
        }
    }
    
    return heterocyclicCount;
}

// 8. CarbocyclicRingCount
CarbocyclicRingCount::CarbocyclicRingCount()
    : Descriptor("CarbocyclicRingCount", "Count of rings containing only carbon atoms")
{
}

std::variant<double, int, std::string> CarbocyclicRingCount::calculate(const ::desfact::Molecule& mol) const
{
    if (!mol.isValid() || !mol.getMolecule()) {
        return 0;
    }
    
    const auto rdmol = mol.getMolecule().get();
    
    // Get ring information
    std::vector<std::vector<int>> sssr;
    RDKit::MolOps::findSSSR(*rdmol, sssr);
    
    int carbocyclicCount = 0;
    
    for (const auto& ring : sssr) {
        bool allCarbon = true;
        for (int atomIdx : ring) {
            const RDKit::Atom* atom = rdmol->getAtomWithIdx(atomIdx);
            if (atom->getAtomicNum() != 6) { // Not carbon
                allCarbon = false;
                break;
            }
        }
        if (allCarbon) {
            carbocyclicCount++;
        }
    }
    
    return carbocyclicCount;
}

// 9. FusedRingSystems
FusedRingSystems::FusedRingSystems()
    : Descriptor("FusedRingSystems", "Count of fused ring systems")
{
}

std::variant<double, int, std::string> FusedRingSystems::calculate(const ::desfact::Molecule& mol) const
{
    if (!mol.isValid() || !mol.getMolecule()) {
        return 0;
    }
    
    const auto rdmol = mol.getMolecule().get();
    
    // Get ring information
    std::vector<std::vector<int>> sssr;
    RDKit::MolOps::findSSSR(*rdmol, sssr);
    
    if (sssr.empty()) {
        return 0; // No rings
    }
    
    // Create a graph where rings are nodes and edges connect fused rings
    std::vector<std::set<int>> ringGraph(sssr.size());
    
    // Identify fused rings (rings that share at least one atom)
    for (size_t i = 0; i < sssr.size(); ++i) {
        const auto& ring1 = sssr[i];
        std::set<int> ring1Atoms(ring1.begin(), ring1.end());
        
        for (size_t j = i + 1; j < sssr.size(); ++j) {
            const auto& ring2 = sssr[j];
            
            // Check if rings share any atoms
            bool fused = false;
            for (int atomIdx : ring2) {
                if (ring1Atoms.count(atomIdx) > 0) {
                    fused = true;
                    break;
                }
            }
            
            if (fused) {
                ringGraph[i].insert(j);
                ringGraph[j].insert(i);
            }
        }
    }
    
    // Count connected components (fused ring systems)
    std::vector<bool> visited(sssr.size(), false);
    int systemCount = 0;
    
    for (size_t i = 0; i < sssr.size(); ++i) {
        if (!visited[i]) {
            // Start a new system
            systemCount++;
            
            // DFS to mark all rings in this system
            std::vector<size_t> stack = {i};
            visited[i] = true;
            
            while (!stack.empty()) {
                size_t current = stack.back();
                stack.pop_back();
                
                for (int neighbor : ringGraph[current]) {
                    if (!visited[neighbor]) {
                        visited[neighbor] = true;
                        stack.push_back(neighbor);
                    }
                }
            }
        }
    }
    
    return systemCount;
}

// 10. AverageFusedRingSize
AverageFusedRingSize::AverageFusedRingSize()
    : Descriptor("AverageFusedRingSize", "Average size of fused ring systems")
{
}

std::variant<double, int, std::string> AverageFusedRingSize::calculate(const ::desfact::Molecule& mol) const
{
    if (!mol.isValid() || !mol.getMolecule()) {
        return 0.0;
    }
    
    const auto rdmol = mol.getMolecule().get();
    
    // Get ring information
    std::vector<std::vector<int>> sssr;
    RDKit::MolOps::findSSSR(*rdmol, sssr);
    
    if (sssr.empty()) {
        return 0.0; // No rings
    }
    
    // Create a graph where rings are nodes and edges connect fused rings
    std::vector<std::set<int>> ringGraph(sssr.size());
    
    // Identify fused rings (rings that share at least one atom)
    for (size_t i = 0; i < sssr.size(); ++i) {
        const auto& ring1 = sssr[i];
        std::set<int> ring1Atoms(ring1.begin(), ring1.end());
        
        for (size_t j = i + 1; j < sssr.size(); ++j) {
            const auto& ring2 = sssr[j];
            
            // Check if rings share any atoms
            bool fused = false;
            for (int atomIdx : ring2) {
                if (ring1Atoms.count(atomIdx) > 0) {
                    fused = true;
                    break;
                }
            }
            
            if (fused) {
                ringGraph[i].insert(j);
                ringGraph[j].insert(i);
            }
        }
    }
    
    // Find connected components (fused ring systems) and calculate sizes
    std::vector<bool> visited(sssr.size(), false);
    std::vector<int> systemSizes; // Number of rings in each system
    
    for (size_t i = 0; i < sssr.size(); ++i) {
        if (!visited[i]) {
            // Start a new system
            int systemSize = 0;
            
            // DFS to mark all rings in this system
            std::vector<size_t> stack = {i};
            visited[i] = true;
            systemSize++;
            
            while (!stack.empty()) {
                size_t current = stack.back();
                stack.pop_back();
                
                for (int neighbor : ringGraph[current]) {
                    if (!visited[neighbor]) {
                        visited[neighbor] = true;
                        systemSize++;
                        stack.push_back(neighbor);
                    }
                }
            }
            
            systemSizes.push_back(systemSize);
        }
    }
    
    if (systemSizes.empty()) {
        return 0.0;
    }
    
    // Calculate average system size
    double totalSize = std::accumulate(systemSizes.begin(), systemSizes.end(), 0.0);
    return totalSize / systemSizes.size();
}

// 11. HeterocyclicFusedRingRatio
HeterocyclicFusedRingRatio::HeterocyclicFusedRingRatio()
    : Descriptor("HeterocyclicFusedRingRatio", "Ratio of heterocyclic rings in fused systems")
{
}

std::variant<double, int, std::string> HeterocyclicFusedRingRatio::calculate(const ::desfact::Molecule& mol) const
{
    if (!mol.isValid() || !mol.getMolecule()) {
        return 0.0;
    }
    
    const auto rdmol = mol.getMolecule().get();
    
    // Get ring information
    std::vector<std::vector<int>> sssr;
    RDKit::MolOps::findSSSR(*rdmol, sssr);
    
    if (sssr.empty()) {
        return 0.0; // No rings
    }
    
    // Create a graph where rings are nodes and edges connect fused rings
    std::vector<std::set<int>> ringGraph(sssr.size());
    
    // Identify fused rings (rings that share at least one atom)
    for (size_t i = 0; i < sssr.size(); ++i) {
        const auto& ring1 = sssr[i];
        std::set<int> ring1Atoms(ring1.begin(), ring1.end());
        
        for (size_t j = i + 1; j < sssr.size(); ++j) {
            const auto& ring2 = sssr[j];
            
            // Check if rings share any atoms
            bool fused = false;
            for (int atomIdx : ring2) {
                if (ring1Atoms.count(atomIdx) > 0) {
                    fused = true;
                    break;
                }
            }
            
            if (fused) {
                ringGraph[i].insert(j);
                ringGraph[j].insert(i);
            }
        }
    }
    
    // Find connected components (fused ring systems) and count heterocyclic rings in each
    std::vector<bool> visited(sssr.size(), false);
    int totalHeterocyclicRings = 0;
    int totalRings = 0;
    
    for (size_t i = 0; i < sssr.size(); ++i) {
        if (!visited[i]) {
            // Start a new system
            std::vector<size_t> systemRings;
            
            // DFS to mark all rings in this system
            std::vector<size_t> stack = {i};
            visited[i] = true;
            systemRings.push_back(i);
            
            while (!stack.empty()) {
                size_t current = stack.back();
                stack.pop_back();
                
                for (int neighbor : ringGraph[current]) {
                    if (!visited[neighbor]) {
                        visited[neighbor] = true;
                        systemRings.push_back(neighbor);
                        stack.push_back(neighbor);
                    }
                }
            }
            
            // Count heterocyclic rings in this system
            int heterocyclicCount = 0;
            for (size_t ringIdx : systemRings) {
                const auto& ring = sssr[ringIdx];
                bool hasHeteroatom = false;
                for (int atomIdx : ring) {
                    const RDKit::Atom* atom = rdmol->getAtomWithIdx(atomIdx);
                    if (atom->getAtomicNum() != 6) { // Not carbon
                        hasHeteroatom = true;
                        break;
                    }
                }
                if (hasHeteroatom) {
                    heterocyclicCount++;
                }
            }
            
            totalHeterocyclicRings += heterocyclicCount;
            totalRings += systemRings.size();
        }
    }
    
    return totalRings > 0 ? static_cast<double>(totalHeterocyclicRings) / totalRings : 0.0;
}

// 12. BridgedRingCount
BridgedRingCount::BridgedRingCount()
    : Descriptor("BridgedRingCount", "Count of bridged rings")
{
}

std::variant<double, int, std::string> BridgedRingCount::calculate(const ::desfact::Molecule& mol) const
{
    if (!mol.isValid() || !mol.getMolecule()) {
        return 0;
    }
    
    const auto rdmol = mol.getMolecule().get();
    
    // Get ring information
    std::vector<std::vector<int>> sssr;
    RDKit::MolOps::findSSSR(*rdmol, sssr);
    
    if (sssr.empty()) {
        return 0; // No rings
    }
    
    int bridgedCount = 0;
    
    // A bridged ring has at least one atom that connects to 3+ other atoms in the ring
    for (const auto& ring : sssr) {
        std::map<int, int> atomConnections;
        
        // Count connections within the ring for each atom
        for (size_t i = 0; i < ring.size(); ++i) {
            int atom1 = ring[i];
            for (size_t j = i + 1; j < ring.size(); ++j) {
                int atom2 = ring[j];
                const RDKit::Bond* bond = rdmol->getBondBetweenAtoms(atom1, atom2);
                if (bond) {
                    atomConnections[atom1]++;
                    atomConnections[atom2]++;
                }
            }
        }
        
        // Check if any atom has 3+ connections (bridged)
        bool isBridged = false;
        for (const auto& [atomIdx, connections] : atomConnections) {
            if (connections >= 3) {
                isBridged = true;
                break;
            }
        }
        
        if (isBridged) {
            bridgedCount++;
        }
    }
    
    return bridgedCount;
}

// 13. SpiroRingCount
SpiroRingCount::SpiroRingCount()
    : Descriptor("SpiroRingCount", "Count of spiro rings")
{
}

std::variant<double, int, std::string> SpiroRingCount::calculate(const ::desfact::Molecule& mol) const
{
    if (!mol.isValid() || !mol.getMolecule()) {
        return 0;
    }
    
    const auto rdmol = mol.getMolecule().get();
    
    // Get ring information
    std::vector<std::vector<int>> sssr;
    RDKit::MolOps::findSSSR(*rdmol, sssr);
    
    if (sssr.size() < 2) {
        return 0; // Need at least 2 rings for spiro
    }
    
    int spiroCount = 0;
    
    // Check each pair of rings for spiro fusion (sharing exactly one atom)
    for (size_t i = 0; i < sssr.size(); ++i) {
        const auto& ring1 = sssr[i];
        std::set<int> ring1Atoms(ring1.begin(), ring1.end());
        
        for (size_t j = i + 1; j < sssr.size(); ++j) {
            const auto& ring2 = sssr[j];
            
            // Count shared atoms
            int sharedAtoms = 0;
            for (int atomIdx : ring2) {
                if (ring1Atoms.count(atomIdx) > 0) {
                    sharedAtoms++;
                }
            }
            
            // Spiro fusion has exactly one shared atom
            if (sharedAtoms == 1) {
                spiroCount++;
            }
        }
    }
    
    return spiroCount;
}

// 14. MacrocyclicRingCount
MacrocyclicRingCount::MacrocyclicRingCount()
    : Descriptor("MacrocyclicRingCount", "Count of rings with size >= 12")
{
}

std::variant<double, int, std::string> MacrocyclicRingCount::calculate(const ::desfact::Molecule& mol) const
{
    if (!mol.isValid() || !mol.getMolecule()) {
        return 0;
    }
    
    const auto rdmol = mol.getMolecule().get();
    
    // Get ring information
    std::vector<std::vector<int>> sssr;
    RDKit::MolOps::findSSSR(*rdmol, sssr);
    
    int macrocyclicCount = 0;
    
    // Count rings with size >= 12
    for (const auto& ring : sssr) {
        if (ring.size() >= 12) {
            macrocyclicCount++;
        }
    }
    
    return macrocyclicCount;
}

// 15. HydrogenBondAcceptorDensity
HydrogenBondAcceptorDensity::HydrogenBondAcceptorDensity()
    : Descriptor("HydrogenBondAcceptorDensity", "Density of H-bond acceptors")
{
}

std::variant<double, int, std::string> HydrogenBondAcceptorDensity::calculate(const ::desfact::Molecule& mol) const
{
    if (!mol.isValid() || !mol.getMolecule()) {
        return 0.0;
    }
    
    const auto rdmol = mol.getMolecule().get();
    int totalAtoms = rdmol->getNumHeavyAtoms();
    
    if (totalAtoms == 0) {
        return 0.0;
    }
    
    // Calculate number of H-bond acceptors
    int acceptorCount = 0;
    
    for (const RDKit::Atom* atom : rdmol->atoms()) {
        int atomicNum = atom->getAtomicNum();
        if (atomicNum == 7 || atomicNum == 8) { // N or O
            // Simple heuristic: N and O atoms with available lone pairs
            int explicitValence = atom->getExplicitValence();
            int implicitValence = atom->getImplicitValence();
            int totalValence = explicitValence + implicitValence;
            
            if ((atomicNum == 7 && totalValence < 3) || // N with lone pair
                (atomicNum == 8 && totalValence < 2)) {  // O with lone pair
                acceptorCount++;
            }
        }
    }
    
    return static_cast<double>(acceptorCount) / totalAtoms;
}

// 16. HydrogenBondDonorDensity
HydrogenBondDonorDensity::HydrogenBondDonorDensity()
    : Descriptor("HydrogenBondDonorDensity", "Density of H-bond donors")
{
}

std::variant<double, int, std::string> HydrogenBondDonorDensity::calculate(const ::desfact::Molecule& mol) const
{
    if (!mol.isValid() || !mol.getMolecule()) {
        return 0.0;
    }
    
    const auto rdmol = mol.getMolecule().get();
    int totalAtoms = rdmol->getNumHeavyAtoms();
    
    if (totalAtoms == 0) {
        return 0.0;
    }
    
    // Calculate number of H-bond donors
    int donorCount = 0;
    
    for (const RDKit::Atom* atom : rdmol->atoms()) {
        int atomicNum = atom->getAtomicNum();
        if (atomicNum == 7 || atomicNum == 8) { // N or O
            // Count explicit hydrogens
            int explicitHs = atom->getNumExplicitHs();
            // Count implicit hydrogens
            int implicitHs = atom->getNumImplicitHs();
            
            if (explicitHs + implicitHs > 0) {
                donorCount++;
            }
        }
    }
    
    return static_cast<double>(donorCount) / totalAtoms;
}

// 17. RotatableBondDensity
RotatableBondDensity::RotatableBondDensity()
    : Descriptor("RotatableBondDensity", "Density of rotatable bonds")
{
}

std::variant<double, int, std::string> RotatableBondDensity::calculate(const ::desfact::Molecule& mol) const
{
    if (!mol.isValid() || !mol.getMolecule()) {
        return 0.0;
    }
    
    const auto rdmol = mol.getMolecule().get();
    int totalBonds = rdmol->getNumBonds();
    
    if (totalBonds == 0) {
        return 0.0;
    }
    
    // Get ring information
    std::vector<std::vector<int>> sssr;
    RDKit::MolOps::findSSSR(*rdmol, sssr);
    
    int rotatableBondCount = 0;
    
    for (const RDKit::Bond* bond : rdmol->bonds()) {
        // Rotatable bonds are single bonds that are not in a ring and not terminal
        if (bond->getBondType() == RDKit::Bond::SINGLE) {
            // Check if atoms are in rings
            bool beginAtomInRing = false;
            bool endAtomInRing = false;
            for (const auto& ring : sssr) {
                if (std::find(ring.begin(), ring.end(), bond->getBeginAtomIdx()) != ring.end()) {
                    beginAtomInRing = true;
                }
                if (std::find(ring.begin(), ring.end(), bond->getEndAtomIdx()) != ring.end()) {
                    endAtomInRing = true;
                }
                if (beginAtomInRing && endAtomInRing) {
                    break;
                }
            }
            
            // If both atoms are in rings and they're in the same ring, the bond is in a ring
            bool bondInRing = false;
            if (beginAtomInRing && endAtomInRing) {
                for (const auto& ring : sssr) {
                    if (std::find(ring.begin(), ring.end(), bond->getBeginAtomIdx()) != ring.end() &&
                        std::find(ring.begin(), ring.end(), bond->getEndAtomIdx()) != ring.end()) {
                        bondInRing = true;
                        break;
                    }
                }
            }
            
            // Rotatable if not in ring and not terminal
            if (!bondInRing && 
                rdmol->getAtomWithIdx(bond->getBeginAtomIdx())->getDegree() > 1 && 
                rdmol->getAtomWithIdx(bond->getEndAtomIdx())->getDegree() > 1) {
                rotatableBondCount++;
            }
        }
    }
    
    return static_cast<double>(rotatableBondCount) / totalBonds;
}

// 18. FunctionalGroupDiversity
FunctionalGroupDiversity::FunctionalGroupDiversity()
    : Descriptor("FunctionalGroupDiversity", "Shannon entropy of functional group types")
{
}

std::variant<double, int, std::string> FunctionalGroupDiversity::calculate(const ::desfact::Molecule& mol) const
{
    if (!mol.isValid() || !mol.getMolecule()) {
        return 0.0;
    }
    
    const auto rdmol = mol.getMolecule().get();
    
    // Define common functional group SMARTS patterns
    std::vector<std::pair<std::string, std::string>> functionalGroups = {
        {"Alcohol", "[OX2H]"},
        {"Aldehyde", "[CX3H1](=O)[#6]"},
        {"Ketone", "[#6][CX3](=O)[#6]"},
        {"Carboxylic Acid", "[CX3](=O)[OX2H1]"},
        {"Ester", "[#6][CX3](=O)[OX2][#6]"},
        {"Amide", "[NX3][CX3](=[OX1])[#6]"},
        {"Amine", "[NX3;H2,H1,H0;!$(NC=O)]"},
        {"Nitro", "[$([NX3](=O)=O),$([NX3+](=O)[O-])]"},
        {"Nitrile", "[NX1]#[CX2]"},
        {"Thiol", "[SX2H]"},
        {"Halogen", "[F,Cl,Br,I]"},
        {"Phosphate", "[$(P(=[OX1])([$([OX2H]),$([OX1-]),$([OX2]P)])([$([OX2H]),$([OX1-]),$([OX2]P)])[$([OX2H]),$([OX1-]),$([OX2]P)]),$([P+]([OX1-])([$([OX2H]),$([OX1-]),$([OX2]P)])([$([OX2H]),$([OX1-]),$([OX2]P)])[$([OX2H]),$([OX1-]),$([OX2]P)])]"},
        {"Sulfate", "[$([SX4](=[OX1])(=[OX1])([$([OX2H]),$([OX1-]),$([OX2]S)])([$([OX2H]),$([OX1-]),$([OX2]S)]),$([SX4+2]([OX1-])([OX1-])([$([OX2H]),$([OX1-]),$([OX2]S)])([$([OX2H]),$([OX1-]),$([OX2]S)]))]"}
    };
    
    std::map<std::string, int> groupCounts;
    int totalGroups = 0;
    
    // Count occurrences of each functional group
    for (const auto& [name, smarts] : functionalGroups) {
        RDKit::ROMol* pattern = nullptr;
        try {
            pattern = RDKit::SmartsToMol(smarts);
            if (pattern) {
                std::vector<std::vector<std::pair<int, int>>> matches;
                int matchCount = RDKit::SubstructMatch(*rdmol, *pattern, matches);
                if (matchCount > 0) {
                    groupCounts[name] = matchCount;
                    totalGroups += matchCount;
                }
                delete pattern;
            }
        } catch (...) {
            if (pattern) delete pattern;
        }
    }
    
    if (totalGroups == 0) {
        return 0.0; // No functional groups found
    }
    
    // Calculate Shannon entropy
    double entropy = 0.0;
    for (const auto& [name, count] : groupCounts) {
        double probability = static_cast<double>(count) / totalGroups;
        entropy -= probability * std::log2(probability);
    }
    
    return entropy;
}

// 19. AtomTypeEntropy
AtomTypeEntropy::AtomTypeEntropy()
    : Descriptor("AtomTypeEntropy", "Shannon entropy of atom types")
{
}

std::variant<double, int, std::string> AtomTypeEntropy::calculate(const ::desfact::Molecule& mol) const
{
    if (!mol.isValid() || !mol.getMolecule()) {
        return 0.0;
    }
    
    const auto rdmol = mol.getMolecule().get();
    int totalAtoms = rdmol->getNumHeavyAtoms();
    
    if (totalAtoms == 0) {
        return 0.0;
    }
    
    std::map<int, int> atomTypeCounts; // Map of atomic number to count
    
    for (const RDKit::Atom* atom : rdmol->atoms()) {
        if (atom->getAtomicNum() > 1) { // Skip hydrogens
            atomTypeCounts[atom->getAtomicNum()]++;
        }
    }
    
    // Calculate Shannon entropy
    double entropy = 0.0;
    for (const auto& [atomicNum, count] : atomTypeCounts) {
        double probability = static_cast<double>(count) / totalAtoms;
        entropy -= probability * std::log2(probability);
    }
    
    return entropy;
}

// 20. BondTypeEntropy
BondTypeEntropy::BondTypeEntropy()
    : Descriptor("BondTypeEntropy", "Shannon entropy of bond types")
{
}

std::variant<double, int, std::string> BondTypeEntropy::calculate(const ::desfact::Molecule& mol) const
{
    if (!mol.isValid() || !mol.getMolecule()) {
        return 0.0;
    }
    
    const auto rdmol = mol.getMolecule().get();
    int totalBonds = rdmol->getNumBonds();
    
    if (totalBonds == 0) {
        return 0.0;
    }
    
    std::map<RDKit::Bond::BondType, int> bondTypeCounts;
    
    for (const RDKit::Bond* bond : rdmol->bonds()) {
        bondTypeCounts[bond->getBondType()]++;
    }
    
    // Calculate Shannon entropy
    double entropy = 0.0;
    for (const auto& [bondType, count] : bondTypeCounts) {
        double probability = static_cast<double>(count) / totalBonds;
        entropy -= probability * std::log2(probability);
    }
    
    return entropy;
}

// 21. RingTypeEntropy
RingTypeEntropy::RingTypeEntropy()
    : Descriptor("RingTypeEntropy", "Shannon entropy of ring types")
{
}

std::variant<double, int, std::string> RingTypeEntropy::calculate(const ::desfact::Molecule& mol) const
{
    if (!mol.isValid() || !mol.getMolecule()) {
        return 0.0;
    }
    
    const auto rdmol = mol.getMolecule().get();
    
    // Get ring information
    std::vector<std::vector<int>> sssr;
    RDKit::MolOps::findSSSR(*rdmol, sssr);
    
    if (sssr.empty()) {
        return 0.0; // No rings
    }
    
    // Categorize rings by size and aromaticity
    std::map<std::pair<int, bool>, int> ringTypeCounts; // (size, isAromatic) -> count
    
    for (const auto& ring : sssr) {
        int ringSize = ring.size();
        
        // Check if ring is aromatic
        bool isAromatic = true;
        for (int atomIdx : ring) {
            const RDKit::Atom* atom = rdmol->getAtomWithIdx(atomIdx);
            if (!atom->getIsAromatic()) {
                isAromatic = false;
                break;
            }
        }
        
        ringTypeCounts[{ringSize, isAromatic}]++;
    }
    
    // Calculate Shannon entropy
    double entropy = 0.0;
    for (const auto& [ringType, count] : ringTypeCounts) {
        double probability = static_cast<double>(count) / sssr.size();
        entropy -= probability * std::log2(probability);
    }
    
    return entropy;
}

// 22. ElectronegativityVariance
ElectronegativityVariance::ElectronegativityVariance()
    : Descriptor("ElectronegativityVariance", "Variance of atom electronegativity values")
{
}

std::variant<double, int, std::string> ElectronegativityVariance::calculate(const ::desfact::Molecule& mol) const
{
    if (!mol.isValid() || !mol.getMolecule()) {
        return 0.0;
    }
    
    const auto rdmol = mol.getMolecule().get();
    int totalAtoms = rdmol->getNumHeavyAtoms();
    
    if (totalAtoms == 0) {
        return 0.0;
    }
    
    // Pauling electronegativity values for common elements
    std::map<int, double> electronegativity = {
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
    
    std::vector<double> values;
    double sum = 0.0;
    
    // Collect electronegativity values
    for (const RDKit::Atom* atom : rdmol->atoms()) {
        int atomicNum = atom->getAtomicNum();
        if (electronegativity.find(atomicNum) != electronegativity.end()) {
            double value = electronegativity[atomicNum];
            values.push_back(value);
            sum += value;
        }
    }
    
    if (values.empty()) {
        return 0.0;
    }
    
    // Calculate mean
    double mean = sum / values.size();
    
    // Calculate variance
    double variance = 0.0;
    for (double value : values) {
        double diff = value - mean;
        variance += diff * diff;
    }
    variance /= values.size();
    
    return variance;
}

} // namespace descriptors
} // namespace desfact
