#include "descriptors/electronic.hpp"
#include <GraphMol/GraphMol.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/Atom.h>
#include <GraphMol/Bond.h>
#include <GraphMol/AtomIterators.h>
#include <GraphMol/BondIterators.h>
#include <GraphMol/PeriodicTable.h>
#include <GraphMol/PartialCharges/GasteigerCharges.h>
#include <GraphMol/RingInfo.h>
#include <vector>
#include <numeric>
#include <cmath>
#include <limits>
#include <queue>          // Keep queue
#include <unordered_map> // Keep unordered_map
#include "utils.hpp"      // Keep utils

namespace desfact {
namespace descriptors {

namespace { // Anonymous namespace for helpers

    // --- Helper Functions ---

    // --- Basic Atom/Bond Properties (from vague3/vague4) ---
    inline unsigned int getAtomDegree(const RDKit::Atom* atom) {
        return atom ? atom->getDegree() : 0;
    }
    inline bool isTerminalAtom(const RDKit::Atom* atom) {
        return atom ? getAtomDegree(atom) <= 1 : false;
    }
    inline bool isAromaticAtom(const RDKit::Atom* atom) {
        return atom && atom->getIsAromatic();
    }
    inline bool isAtomInRing(const RDKit::Atom* atom) {
         // Ensure RingInfo is initialized if needed, though should be by Molecule class ideally
         if (!atom) return false;
         const RDKit::ROMol& mol = atom->getOwningMol();
         const RDKit::RingInfo* ringInfo = mol.getRingInfo();
         if (!ringInfo || !ringInfo->isInitialized()) {
             // This is a fallback, ideally the molecule passed should have SSSR computed
             RDKit::MolOps::symmetrizeSSSR(*const_cast<RDKit::ROMol*>(&mol));
             ringInfo = mol.getRingInfo();
         }
         return ringInfo && ringInfo->numAtomRings(atom->getIdx()) > 0;
    }

    // Check if atom is an H-bond donor (Moved from vague3.cpp)
    inline bool isHBondDonor(const RDKit::Atom* atom) {
        if (!atom) return false;
        int atomicNum = atom->getAtomicNum();
        // Basic definition: O, N, S with explicit H or implicit H if applicable
        if (atomicNum == 7 || atomicNum == 8 || atomicNum == 16) {
             // RDKit's NumExplicitHs() + NumImplicitHs() gives total H count
            return atom->getTotalNumHs() > 0;
        }
        return false;
    }

    // Check if atom is an H-bond acceptor (Moved from vague3.cpp)
    inline bool isHBondAcceptor(const RDKit::Atom* atom) {
        if (!atom) return false;
        int atomicNum = atom->getAtomicNum();
         // Basic definition: N, O, F not positively charged?
         // RDKit's Lipinski definition is more precise. Let's stick to a simple heuristic.
        if ((atomicNum == 7 || atomicNum == 8 || atomicNum == 9) && atom->getFormalCharge() <= 0) {
            // Additional checks could be added (e.g., hybridization)
             return true;
        }
        return false;
    }


    // Get VdW Radius (consistent access)
    double getVdwRadius(const RDKit::Atom* atom) {
        if (!atom) return 0.0;
        return RDKit::PeriodicTable::getTable()->getRvdw(atom->getAtomicNum());
    }

    // Get Electronegativity (Pauling, consistent access)
    double getElectronegativity(const RDKit::Atom* atom) {
         if (!atom) return 0.0;
        // Using a simplified map for common organic elements + H
        static const std::unordered_map<int, double> electronegativityValues = {
            {1, 2.20}, {5, 2.04}, {6, 2.55}, {7, 3.04}, {8, 3.44}, {9, 3.98},
            {14, 1.90},{15, 2.19}, {16, 2.58}, {17, 3.16}, {35, 2.96}, {53, 2.66}
            // Add others if needed, default 2.2 for H if not found? Or 0? Let's use 0 as default for unknown.
        };
        int atomicNum = atom->getAtomicNum();
        auto it = electronegativityValues.find(atomicNum);
        return it != electronegativityValues.end() ? it->second : 0.0;
    }

    // Get Covalent Radius
     double getCovalentRadius(const RDKit::Atom* atom) {
        if (!atom) return 0.0;
        return RDKit::PeriodicTable::getTable()->getRcovalent(atom->getAtomicNum());
     }


    // Compute Gasteiger Charges - potentially expensive, called internally
    bool computeAndGetGasteigerCharges(const RDKit::ROMol* mol, std::vector<double>& charges) {
        if (!mol) return false;
        try {
            // RDKit's computeGasteigerCharges modifies the molecule by adding properties
            // To avoid modifying the const molecule, create a temporary copy
            RDKit::RWMol nonConstMol(*mol);
            RDKit::computeGasteigerCharges(nonConstMol); // Default iterations=12

            charges.resize(mol->getNumAtoms());
            for (unsigned int i = 0; i < mol->getNumAtoms(); ++i) {
                 // Check if the property exists before trying to get it
                 if (nonConstMol.getAtomWithIdx(i)->hasProp(RDKit::common_properties::_GasteigerCharge)) {
                    charges[i] = nonConstMol.getAtomWithIdx(i)->getProp<double>(RDKit::common_properties::_GasteigerCharge);
                 } else {
                    // Handle cases where charge calculation might fail for an atom
                     globalLogger.debug("Gasteiger charge not computed for atom index " + std::to_string(i));
                    charges[i] = 0.0; // Assign default?
                 }
            }
            return true;
        } catch (const std::exception& e) {
            globalLogger.error("Gasteiger charge calculation failed: " + std::string(e.what()));
            return false;
        } catch (...) {
            globalLogger.error("Unknown error during Gasteiger charge calculation.");
            return false;
        }
    }

    // Calculate basic statistics
    struct StatsResult {
        double sum = 0.0;
        double mean = 0.0;
        double stddev = 0.0;
        double minVal = std::numeric_limits<double>::quiet_NaN();
        double maxVal = std::numeric_limits<double>::quiet_NaN();
        double range = std::numeric_limits<double>::quiet_NaN();
        double skewness = std::numeric_limits<double>::quiet_NaN();
        double kurtosis = std::numeric_limits<double>::quiet_NaN();
    };

    StatsResult calculateStatistics(const std::vector<double>& data) {
        StatsResult result;
        std::vector<double> validData;
        
        // Filter out NaN and Inf values
        for (double val : data) {
            if (!std::isnan(val) && !std::isinf(val)) {
                validData.push_back(val);
            }
        }
        
        size_t n = validData.size();
        if (n == 0) {
            // Initialize all fields to 0.0 instead of NaNs
            result.sum = 0.0;
            result.mean = 0.0;
            result.stddev = 0.0;
            result.minVal = 0.0;
            result.maxVal = 0.0;
            result.range = 0.0;
            result.skewness = 0.0;
            result.kurtosis = 0.0;
            return result;
        }

        // Proceed with valid data only
        result.sum = std::accumulate(validData.begin(), validData.end(), 0.0);
        result.mean = result.sum / n;
        result.minVal = *std::min_element(validData.begin(), validData.end());
        result.maxVal = *std::max_element(validData.begin(), validData.end());
        result.range = result.maxVal - result.minVal;

        if (n <= 1) { // Stddev, skew, kurtosis require more data
             result.stddev = 0.0;
             return result;
        }

        double variance_sum = 0.0;
        for (double val : validData) {
            variance_sum += std::pow(val - result.mean, 2);
        }
        double variance = variance_sum / n; // Population variance
        result.stddev = std::sqrt(variance);

        if (n <= 2 || result.stddev < 1e-9) { // Skewness, kurtosis require n>2 and variance > 0
            return result;
        }

        double skew_sum = 0.0;
        double kurt_sum = 0.0;
        for (double val : validData) {
            double term = (val - result.mean) / result.stddev;
            skew_sum += std::pow(term, 3);
            kurt_sum += std::pow(term, 4);
        }

        // Corrected sample skewness (adjust for population if needed)
        // double n_fl = static_cast<double>(n);
        // result.skewness = (n_fl / ((n_fl - 1.0) * (n_fl - 2.0))) * skew_sum; // Sample skewness
        result.skewness = skew_sum / n; // Population skewness

        // Excess kurtosis (adjust for population if needed)
        // result.kurtosis = ((n_fl * (n_fl + 1.0)) / ((n_fl - 1.0) * (n_fl - 2.0) * (n_fl - 3.0))) * kurt_sum
        //                  - (3.0 * std::pow(n_fl - 1.0, 2)) / ((n_fl - 2.0) * (n_fl - 3.0)); // Sample excess kurtosis
         result.kurtosis = (kurt_sum / n) - 3.0; // Population excess kurtosis


        return result;
    }

    // Identify atoms with lone pairs (heuristic) - Moved from vague3.cpp
    bool hasLonePair(const RDKit::Atom* atom) {
        if (!atom) return false;
        int atomicNum = atom->getAtomicNum();
        // Basic check: N, O, S, Halogens with less than typical covalent bonds and not positively charged
        if (atom->getFormalCharge() > 0) return false;

        if (atomicNum == 7) return atom->getExplicitValence() < 3; // N
        if (atomicNum == 8) return atom->getExplicitValence() < 2; // O
        if (atomicNum == 16) return atom->getExplicitValence() < 2; // S
        if (atomicNum == 9 || atomicNum == 17 || atomicNum == 35 || atomicNum == 53) return atom->getExplicitValence() < 1; // Halogens
        // Could refine for P, etc.
        return false;
    }

    // Identify atoms in pi systems
    bool isInPiSystem(const RDKit::Atom* atom) {
        if (!atom) return false;
        return atom->getIsAromatic() ||
               atom->getHybridization() == RDKit::Atom::HybridizationType::SP2 ||
               atom->getHybridization() == RDKit::Atom::HybridizationType::SP;
    }

     // BFS for shortest path (from vague3 helper)
     int getShortestPath(const RDKit::ROMol* mol, unsigned int atom1Idx, unsigned int atom2Idx) {
        if (!mol || atom1Idx >= mol->getNumAtoms() || atom2Idx >= mol->getNumAtoms()) return -1;
        if (atom1Idx == atom2Idx) return 0;

        std::vector<int> distances(mol->getNumAtoms(), -1);
        std::queue<unsigned int> queue; // Corrected: use std::queue

        distances[atom1Idx] = 0;
        queue.push(atom1Idx); // Corrected: use queue

        while (!queue.empty()) { // Corrected: use queue
            unsigned int current = queue.front(); // Corrected: use queue
            queue.pop(); // Corrected: use queue

            if (current == atom2Idx) return distances[current]; // Found

            for (const auto& nbr : mol->atomNeighbors(mol->getAtomWithIdx(current))) {
                unsigned int nbrIdx = nbr->getIdx();
                if (distances[nbrIdx] == -1) {
                    distances[nbrIdx] = distances[current] + 1;
                    queue.push(nbrIdx); // Corrected: use queue
                }
            }
        }
        return -1; // Not connected
    }

    // Identify branch points
    bool isBranchPoint(const RDKit::Atom* atom) {
        return atom && atom->getDegree() >= 3;
    }

     // Get heavy atoms
     std::vector<const RDKit::Atom*> getHeavyAtoms(const RDKit::ROMol* mol) {
         std::vector<const RDKit::Atom*> heavyAtoms;
         if (!mol) return heavyAtoms;
         for (const auto* atom : mol->atoms()) {
             if (atom->getAtomicNum() > 1) {
                 heavyAtoms.push_back(atom);
             }
         }
         return heavyAtoms;
     }


} // anonymous namespace

// --- Descriptor Implementations ---

// --- Category 1: Global Descriptors based on Statistics of (Radius & Charge) ---

SumRadiusChargeDescriptor::SumRadiusChargeDescriptor()
    : Descriptor("SumRadiusCharge", "Sum of (VdW Radius * Gasteiger Partial Charge)") {}

std::variant<double, int, std::string> SumRadiusChargeDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0; // Changed NaN to 0.0
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return 0.0; // Changed NaN to 0.0

    std::vector<double> charges;
    if (!computeAndGetGasteigerCharges(rdkMol, charges)) {
        return "Error: Gasteiger";
    }

    double sumVal = 0.0;
    for (const RDKit::Atom* atom : rdkMol->atoms()) {
        sumVal += getVdwRadius(atom) * charges[atom->getIdx()];
    }
    return std::isfinite(sumVal) ? sumVal : 0.0;
}


AvgRadiusChargeDescriptor::AvgRadiusChargeDescriptor()
    : Descriptor("AvgRadiusCharge", "Average (VdW Radius * Gasteiger Partial Charge)") {}

std::variant<double, int, std::string> AvgRadiusChargeDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0; // Changed NaN to 0.0
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
     if (!rdkMol || rdkMol->getNumAtoms() == 0) return 0.0;


    std::vector<double> charges;
    if (!computeAndGetGasteigerCharges(rdkMol, charges)) {
        return "Error: Gasteiger";
    }

    double sumVal = 0.0;
    for (const RDKit::Atom* atom : rdkMol->atoms()) {
        sumVal += getVdwRadius(atom) * charges[atom->getIdx()];
    }
    double result = sumVal / rdkMol->getNumAtoms();
    return std::isfinite(result) ? result : 0.0;
}


StdDevRadiusChargeDescriptor::StdDevRadiusChargeDescriptor()
    : Descriptor("StdDevRadiusCharge", "Std Dev of (VdW Radius * Gasteiger Partial Charge)") {}

std::variant<double, int, std::string> StdDevRadiusChargeDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol || rdkMol->getNumAtoms() <= 1) return 0.0;

    std::vector<double> charges;
    if (!computeAndGetGasteigerCharges(rdkMol, charges)) return 0.0;

    std::vector<double> values;
    for (const RDKit::Atom* atom : rdkMol->atoms()) {
        double val = getVdwRadius(atom) * charges[atom->getIdx()];
        if (std::isfinite(val)) {  // Both !isnan and !isinf
            values.push_back(val);
        }
    }
    
    if (values.size() <= 1) return 0.0;
    
    // Calculate mean
    double sum = std::accumulate(values.begin(), values.end(), 0.0);
    double mean = sum / values.size();
    
    // Calculate variance and stddev
    double variance = 0.0;
    for (double val : values) {
        variance += (val - mean) * (val - mean);
    }
    variance /= values.size();
    
    double stddev = std::sqrt(variance);
    
    // IMPORTANT: Final check to prevent NaN
    return std::isfinite(stddev) ? stddev : 0.0;
}


RangeRadiusChargeDescriptor::RangeRadiusChargeDescriptor()
    : Descriptor("RangeRadiusCharge", "Range of (VdW Radius * Gasteiger Partial Charge)") {}

std::variant<double, int, std::string> RangeRadiusChargeDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol || rdkMol->getNumAtoms() == 0) return 0.0;

    std::vector<double> charges;
    if (!computeAndGetGasteigerCharges(rdkMol, charges)) return 0.0;

    // Use only non-NaN, non-inf values
    std::vector<double> values;
    for (const RDKit::Atom* atom : rdkMol->atoms()) {
        if (atom->getAtomicNum() <= 1) continue; // Skip hydrogens
        
        double r = getVdwRadius(atom);
        double c = charges[atom->getIdx()];
        double val = r * c;
        
        if (!std::isnan(val) && !std::isinf(val)) {
            values.push_back(val);
        }
    }
    
    if (values.empty()) return 0.0;
    
    auto minmax = std::minmax_element(values.begin(), values.end());
    double range = *minmax.second - *minmax.first;
    
    // IMPORTANT: Catch any potential NaN/inf
    return std::isfinite(range) ? range : 0.0;
}


SkewRadiusChargeDescriptor::SkewRadiusChargeDescriptor()
    : Descriptor("SkewRadiusCharge", "Skewness of (VdW Radius * Gasteiger Partial Charge)") {}

std::variant<double, int, std::string> SkewRadiusChargeDescriptor::calculate(const Molecule& mol) const {
     if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
     if (!rdkMol || rdkMol->getNumAtoms() <= 2) return 0.0; // Skewness requires n > 2

    std::vector<double> charges;
    if (!computeAndGetGasteigerCharges(rdkMol, charges)) return 0.0; // Return 0 on Gasteiger error

    std::vector<double> values;
    for (const RDKit::Atom* atom : rdkMol->atoms()) {
        double val = getVdwRadius(atom) * charges[atom->getIdx()];
        if (std::isfinite(val)) {
            values.push_back(val);
        }
    }
    
    size_t n = values.size();
    if (n <= 2) return 0.0; // Skewness requires n > 2 valid points

    // Calculate mean
    double sum = std::accumulate(values.begin(), values.end(), 0.0);
    double mean = sum / n;
    
    // Calculate stddev
    double variance = 0.0;
    for (double val : values) {
        variance += (val - mean) * (val - mean);
    }
    variance /= n; // Population variance
    double stddev = std::sqrt(variance);

    if (stddev < 1e-9) return 0.0; // Avoid division by zero if stddev is effectively zero

    // Calculate skewness
    double skew_sum = 0.0;
    for (double val : values) {
        skew_sum += std::pow((val - mean) / stddev, 3);
    }
    double skewness = skew_sum / n; // Population skewness

    return std::isfinite(skewness) ? skewness : 0.0; // Final check
}


KurtosisRadiusChargeDescriptor::KurtosisRadiusChargeDescriptor()
    : Descriptor("KurtosisRadiusCharge", "Kurtosis of (VdW Radius * Gasteiger Partial Charge)") {}

std::variant<double, int, std::string> KurtosisRadiusChargeDescriptor::calculate(const Molecule& mol) const {
     if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol || rdkMol->getNumAtoms() <= 3) return 0.0; // Kurtosis often requires n > 3

    std::vector<double> charges;
    if (!computeAndGetGasteigerCharges(rdkMol, charges)) return 0.0; // Return 0 on Gasteiger error

    std::vector<double> values;
    for (const RDKit::Atom* atom : rdkMol->atoms()) {
        double val = getVdwRadius(atom) * charges[atom->getIdx()];
        if (std::isfinite(val)) {
            values.push_back(val);
        }
    }

    size_t n = values.size();
     // Technically needs n > 3 for stable kurtosis, but let's use n > 2 like skewness for consistency here.
     if (n <= 2) return 0.0;

    // Calculate mean
    double sum = std::accumulate(values.begin(), values.end(), 0.0);
    double mean = sum / n;
    
    // Calculate stddev
    double variance = 0.0;
    for (double val : values) {
        variance += (val - mean) * (val - mean);
    }
    variance /= n; // Population variance
    double stddev = std::sqrt(variance);

    if (stddev < 1e-9) return 0.0; // Avoid division by zero

    // Calculate kurtosis
    double kurt_sum = 0.0;
    for (double val : values) {
        kurt_sum += std::pow((val - mean) / stddev, 4);
    }
     // Population excess kurtosis
    double kurtosis = (kurt_sum / n) - 3.0; 

    return std::isfinite(kurtosis) ? kurtosis : 0.0; // Final check
}


SumRadiusPerAbsChargeDescriptor::SumRadiusPerAbsChargeDescriptor()
    : Descriptor("SumRadiusPerAbsCharge", "Sum of VdW Radius / (|Gasteiger Charge| + epsilon)") {}

std::variant<double, int, std::string> SumRadiusPerAbsChargeDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0; // Changed NaN to 0.0
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return 0.0; // Changed NaN to 0.0

    std::vector<double> charges;
    if (!computeAndGetGasteigerCharges(rdkMol, charges)) {
        return "Error: Gasteiger";
    }

    double sumVal = 0.0;
    double epsilon = 1e-6; // Small value to prevent division by zero

    for (const RDKit::Atom* atom : rdkMol->atoms()) {
        double absCharge = std::abs(charges[atom->getIdx()]);
        sumVal += getVdwRadius(atom) / (absCharge + epsilon);
    }
    return std::isfinite(sumVal) ? sumVal : 0.0;
}


AvgRadiusPositiveChargeDescriptor::AvgRadiusPositiveChargeDescriptor()
    : Descriptor("AvgRadiusPositiveCharge", "Average VdW Radius of positively charged atoms") {}

std::variant<double, int, std::string> AvgRadiusPositiveChargeDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0; // Changed NaN to 0.0
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return 0.0; // Changed NaN to 0.0

    std::vector<double> charges;
    if (!computeAndGetGasteigerCharges(rdkMol, charges)) {
        return "Error: Gasteiger";
    }

    double sumRadius = 0.0;
    int count = 0;
    for (const RDKit::Atom* atom : rdkMol->atoms()) {
        if (charges[atom->getIdx()] > 0) {
            sumRadius += getVdwRadius(atom);
            count++;
        }
    }
    return (count > 0) ? sumRadius / count : 0.0; // Return 0 if no positive atoms
}


AvgRadiusNegativeChargeDescriptor::AvgRadiusNegativeChargeDescriptor()
    : Descriptor("AvgRadiusNegativeCharge", "Average VdW Radius of negatively charged atoms") {}

std::variant<double, int, std::string> AvgRadiusNegativeChargeDescriptor::calculate(const Molecule& mol) const {
     if (!mol.isValid() || !mol.getMolecule()) return 0.0; // Changed NaN to 0.0
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
     if (!rdkMol) return 0.0; // Changed NaN to 0.0


    std::vector<double> charges;
    if (!computeAndGetGasteigerCharges(rdkMol, charges)) {
        return "Error: Gasteiger";
    }

    double sumRadius = 0.0;
    int count = 0;
    for (const RDKit::Atom* atom : rdkMol->atoms()) {
        if (charges[atom->getIdx()] < 0) {
            sumRadius += getVdwRadius(atom);
            count++;
        }
    }
    return (count > 0) ? sumRadius / count : 0.0; // Return 0 if no negative atoms
}


RatioSumRadiusChargedDescriptor::RatioSumRadiusChargedDescriptor()
    : Descriptor("RatioSumRadiusCharged", "Ratio: Sum VdW Radius (Pos Charge) / Sum VdW Radius (Neg Charge)") {}

std::variant<double, int, std::string> RatioSumRadiusChargedDescriptor::calculate(const Molecule& mol) const {
     if (!mol.isValid() || !mol.getMolecule()) return 0.0; // Changed NaN to 0.0
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
     if (!rdkMol) return 0.0; // Changed NaN to 0.0


    std::vector<double> charges;
    if (!computeAndGetGasteigerCharges(rdkMol, charges)) {
        return "Error: Gasteiger";
    }

    double sumRadiusPos = 0.0;
    double sumRadiusNeg = 0.0;

    for (const RDKit::Atom* atom : rdkMol->atoms()) {
         double charge = charges[atom->getIdx()];
         if (charge > 0) {
             sumRadiusPos += getVdwRadius(atom);
         } else if (charge < 0) {
             sumRadiusNeg += getVdwRadius(atom);
         }
    }

     if (sumRadiusNeg < 1e-9) { // Avoid division by zero
         return (sumRadiusPos > 1e-9) ? 9999.9 : 1.0; // Keep large num for Inf, 1.0 for 0/0
     }

    return sumRadiusPos / sumRadiusNeg;
}


WeightedAvgElectronegativityDescriptor::WeightedAvgElectronegativityDescriptor()
    : Descriptor("WeightedAvgElectronegativity", "Weighted Average Electronegativity (by VdW Radius)") {}

std::variant<double, int, std::string> WeightedAvgElectronegativityDescriptor::calculate(const Molecule& mol) const {
     if (!mol.isValid() || !mol.getMolecule()) return 0.0; // Changed NaN to 0.0
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
     if (!rdkMol || rdkMol->getNumAtoms() == 0) return 0.0;


    double weightedSumEN = 0.0;
    double sumWeights = 0.0;

    for (const RDKit::Atom* atom : rdkMol->atoms()) {
        double radius = getVdwRadius(atom);
        weightedSumEN += getElectronegativity(atom) * radius;
        sumWeights += radius;
    }

     if (sumWeights < 1e-9) return 0.0; // Avoid division by zero

    return weightedSumEN / sumWeights;
}


WeightedStdDevElectronegativityDescriptor::WeightedStdDevElectronegativityDescriptor()
    : Descriptor("WeightedStdDevElectronegativity", "Std Dev of Electronegativity (weighted by VdW Radius)") {}

std::variant<double, int, std::string> WeightedStdDevElectronegativityDescriptor::calculate(const Molecule& mol) const {
     if (!mol.isValid() || !mol.getMolecule()) return 0.0; // Changed NaN to 0.0
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
     if (!rdkMol || rdkMol->getNumAtoms() <= 1) return 0.0; // Need > 1 atom


    std::vector<double> enValues;
    std::vector<double> weights; // VdW radii
    enValues.reserve(rdkMol->getNumAtoms());
    weights.reserve(rdkMol->getNumAtoms());

    for (const RDKit::Atom* atom : rdkMol->atoms()) {
        enValues.push_back(getElectronegativity(atom));
        weights.push_back(getVdwRadius(atom));
    }

    // Calculate weighted mean
    double weightedSumEN = 0.0;
    double sumWeights = 0.0;
    for (size_t i = 0; i < enValues.size(); ++i) {
        weightedSumEN += enValues[i] * weights[i];
        sumWeights += weights[i];
    }

     if (sumWeights < 1e-9) return 0.0; // All weights zero?
    double weightedMean = weightedSumEN / sumWeights;

    // Calculate weighted variance
    double weightedVarianceSum = 0.0;
    for (size_t i = 0; i < enValues.size(); ++i) {
        weightedVarianceSum += weights[i] * std::pow(enValues[i] - weightedMean, 2);
    }

    // Using sum of weights as denominator for population variance estimate
    double weightedVariance = weightedVarianceSum / sumWeights;

    return std::isfinite(std::sqrt(weightedVariance)) ? std::sqrt(weightedVariance) : 0.0;
}


SumRadiusSqAbsChargeDescriptor::SumRadiusSqAbsChargeDescriptor()
    : Descriptor("SumRadiusSqAbsCharge", "Sum of VdW Radius^2 * abs(Gasteiger Charge)") {}

std::variant<double, int, std::string> SumRadiusSqAbsChargeDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0; // Changed NaN to 0.0
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return 0.0; // Changed NaN to 0.0

    std::vector<double> charges;
    if (!computeAndGetGasteigerCharges(rdkMol, charges)) {
        return "Error: Gasteiger";
    }

    double sumVal = 0.0;
    for (const RDKit::Atom* atom : rdkMol->atoms()) {
         double radius = getVdwRadius(atom);
         sumVal += (radius * radius) * std::abs(charges[atom->getIdx()]);
    }
    return std::isfinite(sumVal) ? sumVal : 0.0;
}


WeightedStdDevChargeDescriptor::WeightedStdDevChargeDescriptor()
    : Descriptor("WeightedStdDevCharge", "Std Dev of Gasteiger Charge (weighted by VdW Radius)") {}

std::variant<double, int, std::string> WeightedStdDevChargeDescriptor::calculate(const Molecule& mol) const {
     if (!mol.isValid() || !mol.getMolecule()) return 0.0; // Changed NaN to 0.0
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
     if (!rdkMol || rdkMol->getNumAtoms() <= 1) return 0.0; // Need > 1 atom


    std::vector<double> charges;
    if (!computeAndGetGasteigerCharges(rdkMol, charges)) return 0.0;

    // Collect valid charge-weight pairs
    std::vector<std::pair<double, double>> data;
    double totalWeight = 0.0;
    
    for (const RDKit::Atom* atom : rdkMol->atoms()) {
        double charge = charges[atom->getIdx()];
        double weight = getVdwRadius(atom);
        
        if (std::isfinite(charge) && std::isfinite(weight) && weight > 0) {
            data.emplace_back(charge, weight);
            totalWeight += weight;
        }
    }
    
    if (data.empty() || totalWeight <= 0) return 0.0;
    
    // Calculate weighted mean
    double weightedSum = 0.0;
    for (const auto& [charge, weight] : data) {
        weightedSum += charge * weight;
    }
    double weightedMean = weightedSum / totalWeight;
    
    // Calculate weighted variance
    double variance = 0.0;
    for (const auto& [charge, weight] : data) {
        variance += weight * (charge - weightedMean) * (charge - weightedMean);
    }
    variance /= totalWeight;
    
    double stddev = std::sqrt(variance);
    
    // Final safety check
    return std::isfinite(stddev) ? stddev : 0.0;
}


TotalAbsChargePerRadiusDescriptor::TotalAbsChargePerRadiusDescriptor()
    : Descriptor("TotalAbsChargePerRadius", "Total (Abs(Charge) / VdW Radius)") {}

std::variant<double, int, std::string> TotalAbsChargePerRadiusDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0; // Changed NaN to 0.0
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return 0.0; // Changed NaN to 0.0

    std::vector<double> charges;
    if (!computeAndGetGasteigerCharges(rdkMol, charges)) {
        return "Error: Gasteiger";
    }

    double sumVal = 0.0;
    double epsilon = 1e-6; // Prevent division by zero radius

    for (const RDKit::Atom* atom : rdkMol->atoms()) {
         double radius = getVdwRadius(atom);
         sumVal += std::abs(charges[atom->getIdx()]) / (radius + epsilon);
    }
    return std::isfinite(sumVal) ? sumVal : 0.0;
}


SumInvRadiusAbsChargeDescriptor::SumInvRadiusAbsChargeDescriptor()
    : Descriptor("SumInvRadiusAbsCharge", "Sum of (1 / VdW Radius * Abs(Gasteiger Charge))") {}

std::variant<double, int, std::string> SumInvRadiusAbsChargeDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0; // Changed NaN to 0.0
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return 0.0; // Changed NaN to 0.0

    std::vector<double> charges;
    if (!computeAndGetGasteigerCharges(rdkMol, charges)) {
        return "Error: Gasteiger";
    }

    double sumVal = 0.0;
    double epsilon = 1e-6; // Prevent division by zero radius

    for (const RDKit::Atom* atom : rdkMol->atoms()) {
         double radius = getVdwRadius(atom);
         sumVal += (1.0 / (radius + epsilon)) * std::abs(charges[atom->getIdx()]);
    }
    return std::isfinite(sumVal) ? sumVal : 0.0;
}


SumEnRadiusSqDescriptor::SumEnRadiusSqDescriptor()
    : Descriptor("SumEnRadiusSq", "Sum of Electronegativity * (VdW Radius)^2") {}

std::variant<double, int, std::string> SumEnRadiusSqDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0; // Changed NaN to 0.0
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return 0.0; // Changed NaN to 0.0

    double sumVal = 0.0;
    for (const RDKit::Atom* atom : rdkMol->atoms()) {
         double radius = getVdwRadius(atom);
         sumVal += getElectronegativity(atom) * (radius * radius);
    }
    return std::isfinite(sumVal) ? sumVal : 0.0;
}


RatioAvgEnRadiusRingChainDescriptor::RatioAvgEnRadiusRingChainDescriptor()
    : Descriptor("RatioAvgEnRadiusRingChain", "Ratio of Avg(EN*Radius) for ring atoms vs chain atoms") {}

std::variant<double, int, std::string> RatioAvgEnRadiusRingChainDescriptor::calculate(const Molecule& mol) const {
     if (!mol.isValid() || !mol.getMolecule()) return 0.0; // Changed NaN to 0.0
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
     if (!rdkMol) return 0.0; // Changed NaN to 0.0

    const RDKit::RingInfo* ringInfo = rdkMol->getRingInfo();
     if (!ringInfo || !ringInfo->isInitialized()) RDKit::MolOps::symmetrizeSSSR(*const_cast<RDKit::ROMol*>(rdkMol));
     ringInfo = rdkMol->getRingInfo();
      if (!ringInfo) return 0.0; // Changed NaN to 0.0


    double sumEnRadiusRing = 0.0;
    int countRing = 0;
    double sumEnRadiusChain = 0.0;
    int countChain = 0;

    for (const RDKit::Atom* atom : rdkMol->atoms()) {
         double en = getElectronegativity(atom);
         double radius = getVdwRadius(atom);
         double val = en * radius;
         if (ringInfo->numAtomRings(atom->getIdx()) > 0) {
             sumEnRadiusRing += val;
             countRing++;
         } else {
             sumEnRadiusChain += val;
             countChain++;
         }
    }

     if (countRing == 0 && countChain == 0) return 1.0; // No atoms? Define as 1?
     if (countChain == 0) return (countRing > 0) ? 9999.9 : 1.0; // Use large number for inf, 1.0 for 0/0
     if (countRing == 0) return 0.0; // 0 / N case

     double avgRing = sumEnRadiusRing / countRing;
     double avgChain = sumEnRadiusChain / countChain;

     if (std::abs(avgChain) < 1e-9) {
         return (std::abs(avgRing) > 1e-9) ? 9999.9 : 1.0; // Use large number for inf, 1.0 for 0/0
     }

    return avgRing / avgChain;
}


AvgRatioChargeRadiusDescriptor::AvgRatioChargeRadiusDescriptor()
    : Descriptor("AvgRatioChargeRadius", "Average Ratio (Gasteiger Charge / VdW Radius) over all atoms") {}

std::variant<double, int, std::string> AvgRatioChargeRadiusDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0; // Changed NaN to 0.0
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
     if (!rdkMol || rdkMol->getNumAtoms() == 0) return 0.0;


    std::vector<double> charges;
    if (!computeAndGetGasteigerCharges(rdkMol, charges)) {
        return "Error: Gasteiger";
    }

    double sumRatio = 0.0;
    double epsilon = 1e-6;

    for (const RDKit::Atom* atom : rdkMol->atoms()) {
         double radius = getVdwRadius(atom);
         sumRatio += charges[atom->getIdx()] / (radius + epsilon);
    }
    double result = sumRatio / rdkMol->getNumAtoms();
    return std::isfinite(result) ? result : 0.0;
}


RatioAvgRChargeToAvgRPlusChargeDescriptor::RatioAvgRChargeToAvgRPlusChargeDescriptor()
    : Descriptor("RatioAvgRChargeToAvgRPlusCharge", "Ratio of Avg (R * |q|) to Avg (R + |q|)") {}

std::variant<double, int, std::string> RatioAvgRChargeToAvgRPlusChargeDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0; // Changed NaN to 0.0
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
     if (!rdkMol || rdkMol->getNumAtoms() == 0) return 1.0; // Or 0.0? Ratio is undefined. Let's use 1.0

    std::vector<double> charges;
    if (!computeAndGetGasteigerCharges(rdkMol, charges)) {
        return "Error: Gasteiger";
    }

    double sumRq = 0.0;
    double sumRplusQ = 0.0;

    for (const RDKit::Atom* atom : rdkMol->atoms()) {
         double radius = getVdwRadius(atom);
         double absCharge = std::abs(charges[atom->getIdx()]);
         sumRq += radius * absCharge;
         sumRplusQ += radius + absCharge;
    }

     double avgRq = sumRq / rdkMol->getNumAtoms();
     double avgRplusQ = sumRplusQ / rdkMol->getNumAtoms();

      if (std::abs(avgRplusQ) < 1e-9) {
          return (std::abs(avgRq) > 1e-9) ? 9999.9 : 1.0; // Avoid division by zero, handle 0/0
      }

    return avgRq / avgRplusQ;
}


AvgRadiusLonePairAtomsDescriptor::AvgRadiusLonePairAtomsDescriptor()
    : Descriptor("AvgRadiusLonePairAtoms", "Avg VdW Radius for atoms with lone pairs (heuristic)") {}

std::variant<double, int, std::string> AvgRadiusLonePairAtomsDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0; // Changed NaN to 0.0
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return 0.0; // Changed NaN to 0.0

    double sumRadius = 0.0;
    int count = 0;
    for (const RDKit::Atom* atom : rdkMol->atoms()) {
        if (hasLonePair(atom)) {
            sumRadius += getVdwRadius(atom);
            count++;
        }
    }
    return (count > 0) ? sumRadius / count : 0.0; // Return 0 if no atoms with lone pairs found
}


AvgRadiusPiSystemAtomsDescriptor::AvgRadiusPiSystemAtomsDescriptor()
    : Descriptor("AvgRadiusPiSystemAtoms", "Avg VdW Radius for atoms in pi systems") {}

std::variant<double, int, std::string> AvgRadiusPiSystemAtomsDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0; // Changed NaN to 0.0
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return 0.0; // Changed NaN to 0.0

    double sumRadius = 0.0;
    int count = 0;
    for (const RDKit::Atom* atom : rdkMol->atoms()) {
        if (isInPiSystem(atom)) {
            sumRadius += getVdwRadius(atom);
            count++;
        }
    }
    return (count > 0) ? sumRadius / count : 0.0;
}


StdDevRatioChargeRadiusDescriptor::StdDevRatioChargeRadiusDescriptor()
    : Descriptor("StdDevRatioChargeRadius", "Std Dev of Ratios (|Gasteiger Charge| / VdW Radius)") {}

std::variant<double, int, std::string> StdDevRatioChargeRadiusDescriptor::calculate(const Molecule& mol) const {
     if (!mol.isValid() || !mol.getMolecule()) return 0.0; // Changed NaN to 0.0
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
     if (!rdkMol || rdkMol->getNumAtoms() == 0) return 0.0;

    std::vector<double> charges;
    if (!computeAndGetGasteigerCharges(rdkMol, charges)) {
        return "Error: Gasteiger";
    }

    std::vector<double> ratios;
    ratios.reserve(rdkMol->getNumAtoms());
    double epsilon = 1e-6;

    for (const RDKit::Atom* atom : rdkMol->atoms()) {
         double radius = getVdwRadius(atom);
         ratios.push_back(std::abs(charges[atom->getIdx()]) / (radius + epsilon));
    }

    return std::isfinite(calculateStatistics(ratios).stddev) ? calculateStatistics(ratios).stddev : 0.0;
}


SumRadiusPositiveChargeThreshDescriptor::SumRadiusPositiveChargeThreshDescriptor()
    : Descriptor("SumRadiusPositiveChargeThresh", "Sum of VdW Radius for atoms with Gasteiger Charge > 0.1") {}

std::variant<double, int, std::string> SumRadiusPositiveChargeThreshDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0; // Changed NaN to 0.0
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return 0.0; // Changed NaN to 0.0

    std::vector<double> charges;
    if (!computeAndGetGasteigerCharges(rdkMol, charges)) {
        return "Error: Gasteiger";
    }

    double sumRadius = 0.0;
    double threshold = 0.1;
    for (const RDKit::Atom* atom : rdkMol->atoms()) {
        if (charges[atom->getIdx()] > threshold) {
            sumRadius += getVdwRadius(atom);
        }
    }
    return std::isfinite(sumRadius) ? sumRadius : 0.0;
}


SumRadiusNegativeChargeThreshDescriptor::SumRadiusNegativeChargeThreshDescriptor()
    : Descriptor("SumRadiusNegativeChargeThresh", "Sum of VdW Radius for atoms with Gasteiger Charge < -0.1") {}

std::variant<double, int, std::string> SumRadiusNegativeChargeThreshDescriptor::calculate(const Molecule& mol) const {
     if (!mol.isValid() || !mol.getMolecule()) return 0.0; // Changed NaN to 0.0
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
     if (!rdkMol) return 0.0; // Changed NaN to 0.0


    std::vector<double> charges;
    if (!computeAndGetGasteigerCharges(rdkMol, charges)) {
        return "Error: Gasteiger";
    }

    double sumRadius = 0.0;
    double threshold = -0.1;
    for (const RDKit::Atom* atom : rdkMol->atoms()) {
        if (charges[atom->getIdx()] < threshold) {
            sumRadius += getVdwRadius(atom);
        }
    }
    return std::isfinite(sumRadius) ? sumRadius : 0.0;
}


// Topological Descriptors (Placeholders/Approximations)
WeightedPathCountDescriptor::WeightedPathCountDescriptor()
    : Descriptor("WeightedPathCount", "Placeholder: Weighted path count (R*|q|)") {}

std::variant<double, int, std::string> WeightedPathCountDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol || rdkMol->getNumHeavyAtoms() < 2) return 0.0;
    
    std::vector<double> charges;
    if (!computeAndGetGasteigerCharges(rdkMol, charges)) return 0.0;
    
    // Calculate R*|q| weights for each atom
    std::vector<double> weights(rdkMol->getNumAtoms());
    for (const RDKit::Atom* atom : getHeavyAtoms(rdkMol)) {
        weights[atom->getIdx()] = getVdwRadius(atom) * std::abs(charges[atom->getIdx()]);
    }
    
    // Count paths of length 2-4, weighted by endpoints
    double weightedPathCount = 0.0;
    for (size_t i = 0; i < rdkMol->getNumAtoms(); ++i) {
        for (size_t j = i+1; j < rdkMol->getNumAtoms(); ++j) {
            int dist = getShortestPath(rdkMol, i, j);
            if (dist >= 2 && dist <= 4) {
                weightedPathCount += weights[i] * weights[j];
            }
        }
    }
    
    return std::isfinite(weightedPathCount) ? weightedPathCount : 0.0;
}

AvgWeightedPathLengthDescriptor::AvgWeightedPathLengthDescriptor()
    : Descriptor("AvgWeightedPathLength", "Placeholder: Average weighted path length (R*|q|)") {}

std::variant<double, int, std::string> AvgWeightedPathLengthDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol || rdkMol->getNumHeavyAtoms() < 2) return 0.0;
    
    std::vector<double> charges;
    if (!computeAndGetGasteigerCharges(rdkMol, charges)) return 0.0;
    
    double sumWeightedPaths = 0.0;
    double sumWeights = 0.0;
    int pathCount = 0;
    
    auto heavyAtoms = getHeavyAtoms(rdkMol);
    for (size_t i = 0; i < heavyAtoms.size(); ++i) {
        for (size_t j = i+1; j < heavyAtoms.size(); ++j) {
            int dist = getShortestPath(rdkMol, heavyAtoms[i]->getIdx(), heavyAtoms[j]->getIdx());
            if (dist > 0) {
                double weight1 = getVdwRadius(heavyAtoms[i]) * std::abs(charges[heavyAtoms[i]->getIdx()]);
                double weight2 = getVdwRadius(heavyAtoms[j]) * std::abs(charges[heavyAtoms[j]->getIdx()]);
                double pathWeight = weight1 * weight2;
                
                sumWeightedPaths += dist * pathWeight;
                sumWeights += pathWeight;
                pathCount++;
            }
        }
    }
    
    if (pathCount == 0 || sumWeights < 1e-9) return 0.0;
    double result = sumWeightedPaths / sumWeights;
    return std::isfinite(result) ? result : 0.0;
}

BalabanLikeIndexRChargeDescriptor::BalabanLikeIndexRChargeDescriptor()
    : Descriptor("BalabanLikeIndexRCharge", "Placeholder: Balaban-like index (R+|q|)") {}

std::variant<double, int, std::string> BalabanLikeIndexRChargeDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol || rdkMol->getNumHeavyAtoms() < 3) return 0.0;
    
    std::vector<double> charges;
    if (!computeAndGetGasteigerCharges(rdkMol, charges)) return 0.0;
    
    // Compute distance matrix
    std::vector<std::vector<int>> distMatrix;
    auto heavyAtoms = getHeavyAtoms(rdkMol);
    size_t n = heavyAtoms.size();
    distMatrix.resize(n, std::vector<int>(n, 0));
    
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i+1; j < n; ++j) {
            int dist = getShortestPath(rdkMol, heavyAtoms[i]->getIdx(), heavyAtoms[j]->getIdx());
            if (dist <= 0) dist = 999; // Not connected
            distMatrix[i][j] = distMatrix[j][i] = dist;
        }
    }
    
    // Calculate sum of distance terms for each vertex
    std::vector<double> distSums(n, 0.0);
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            if (i != j && distMatrix[i][j] < 999) {
                distSums[i] += distMatrix[i][j];
            }
        }
    }
    
    // Calculate product of weights for bond contributions
    double sum = 0.0;
    for (const RDKit::Bond* bond : rdkMol->bonds()) {
        unsigned int idx1 = bond->getBeginAtomIdx();
        unsigned int idx2 = bond->getEndAtomIdx();
        
        // Find corresponding indices in heavyAtoms
        size_t i = 0, j = 0;
        bool found1 = false, found2 = false;
        for (size_t k = 0; k < heavyAtoms.size(); ++k) {
            if (heavyAtoms[k]->getIdx() == idx1) {
                i = k;
                found1 = true;
            }
            if (heavyAtoms[k]->getIdx() == idx2) {
                j = k;
                found2 = true;
            }
        }
        
        if (!found1 || !found2) continue; // Skip hydrogens
        
        double weight1 = getVdwRadius(heavyAtoms[i]) * std::abs(charges[idx1]);
        double weight2 = getVdwRadius(heavyAtoms[j]) * std::abs(charges[idx2]);
        
        if (distSums[i] > 0 && distSums[j] > 0) {
            sum += 1.0 / std::sqrt(distSums[i] * distSums[j] * weight1 * weight2);
        }
    }
    
    // Calculate cyclomatic number
    int cycloNum = rdkMol->getNumBonds() - rdkMol->getNumHeavyAtoms() + 1;
    if (cycloNum <= 0) return 0.0;
    
    return std::isfinite(static_cast<double>(rdkMol->getNumBonds()) * sum / cycloNum) ? static_cast<double>(rdkMol->getNumBonds()) * sum / cycloNum : 0.0;
}

WienerLikeIndexRChargeDescriptor::WienerLikeIndexRChargeDescriptor()
    : Descriptor("WienerLikeIndexRCharge", "Placeholder: Wiener-like index (R*|q|)") {}

std::variant<double, int, std::string> WienerLikeIndexRChargeDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol || rdkMol->getNumHeavyAtoms() < 2) return 0.0;
    
    std::vector<double> charges;
    if (!computeAndGetGasteigerCharges(rdkMol, charges)) return 0.0;
    
    double wienerSum = 0.0;
    auto heavyAtoms = getHeavyAtoms(rdkMol);
    
    for (size_t i = 0; i < heavyAtoms.size(); ++i) {
        for (size_t j = i+1; j < heavyAtoms.size(); ++j) {
            int dist = getShortestPath(rdkMol, heavyAtoms[i]->getIdx(), heavyAtoms[j]->getIdx());
            if (dist > 0) {
                double weight1 = getVdwRadius(heavyAtoms[i]) * std::abs(charges[heavyAtoms[i]->getIdx()]);
                double weight2 = getVdwRadius(heavyAtoms[j]) * std::abs(charges[heavyAtoms[j]->getIdx()]);
                
                wienerSum += dist * weight1 * weight2;
            }
        }
    }
    
    return std::isfinite(wienerSum) ? wienerSum : 0.0;
}

EccentricityMaxRadiusAtomChargeWeightedDescriptor::EccentricityMaxRadiusAtomChargeWeightedDescriptor()
    : Descriptor("EccentricityMaxRadiusAtomChargeWeighted", "Eccentricity of max R atom, weighted by its charge") {}

std::variant<double, int, std::string> EccentricityMaxRadiusAtomChargeWeightedDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0; // Changed NaN to 0.0
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
     if (!rdkMol || rdkMol->getNumHeavyAtoms() <= 1) return 0.0; // Needs heavy atoms for radius/eccentricity

    std::vector<double> charges;
    if (!computeAndGetGasteigerCharges(rdkMol, charges)) {
        return "Error: Gasteiger";
    }

    int maxRAtomIdx = -1;
    double maxR = -1.0;
    for (const RDKit::Atom* atom : getHeavyAtoms(rdkMol)) {
         double r = getVdwRadius(atom);
         if (r > maxR) {
             maxR = r;
             maxRAtomIdx = atom->getIdx();
         }
    }
    if (maxRAtomIdx == -1) return 0.0; // No heavy atoms found


    // Calculate eccentricity (max shortest path from this atom)
    int maxPath = 0;
     for (const RDKit::Atom* otherAtom : getHeavyAtoms(rdkMol)) {
         if (otherAtom->getIdx() == static_cast<unsigned int>(maxRAtomIdx)) continue;
         int dist = getShortestPath(rdkMol, maxRAtomIdx, otherAtom->getIdx());
         if (dist > maxPath) {
             maxPath = dist;
         }
     }

     return std::isfinite(static_cast<double>(maxPath) * charges[maxRAtomIdx]) ? static_cast<double>(maxPath) * charges[maxRAtomIdx] : 0.0; // Weight by charge
}

TopoDistMaxRMaxPosChargeDescriptor::TopoDistMaxRMaxPosChargeDescriptor()
    : Descriptor("TopoDistMaxRMaxPosCharge", "Topological distance between max R atom and max positive charge atom") {}

std::variant<double, int, std::string> TopoDistMaxRMaxPosChargeDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return -1; // Return distance as int or -1
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
     if (!rdkMol || rdkMol->getNumHeavyAtoms() <= 1) return 0; // Return 0 for int descriptor

    std::vector<double> charges;
    if (!computeAndGetGasteigerCharges(rdkMol, charges)) {
        return "Error: Gasteiger";
    }

    int maxRAtomIdx = -1;
    double maxR = -1.0;
    int maxPosQAtomIdx = -1;
    double maxPosQ = -1.0; // Initialize below possible positive charge

     for (const RDKit::Atom* atom : getHeavyAtoms(rdkMol)) {
          int idx = atom->getIdx();
          double r = getVdwRadius(atom);
          double q = charges[idx];

          if (r > maxR) {
              maxR = r;
              maxRAtomIdx = idx;
          }
          if (q > 0 && q > maxPosQ) {
              maxPosQ = q;
              maxPosQAtomIdx = idx;
          }
     }

     if (maxRAtomIdx == -1 || maxPosQAtomIdx == -1) return 0; // If no max R or no pos charge atom found
     if (maxRAtomIdx == maxPosQAtomIdx) return 0; // Same atom

     return getShortestPath(rdkMol, maxRAtomIdx, maxPosQAtomIdx); // Returns -1 if not connected
}

TopoDistMaxRMinNegChargeDescriptor::TopoDistMaxRMinNegChargeDescriptor()
    : Descriptor("TopoDistMaxRMinNegCharge", "Topological distance between max R atom and min negative charge atom") {}

std::variant<double, int, std::string> TopoDistMaxRMinNegChargeDescriptor::calculate(const Molecule& mol) const {
     if (!mol.isValid() || !mol.getMolecule()) return -1; // Return distance as int or -1
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
     if (!rdkMol || rdkMol->getNumHeavyAtoms() <= 1) return 0; // Return 0 for int descriptor

    std::vector<double> charges;
    if (!computeAndGetGasteigerCharges(rdkMol, charges)) {
        return "Error: Gasteiger";
    }

    int maxRAtomIdx = -1;
    double maxR = -1.0;
    int minNegQAtomIdx = -1;
    double minNegQ = 1.0; // Initialize above possible negative charge

     for (const RDKit::Atom* atom : getHeavyAtoms(rdkMol)) {
          int idx = atom->getIdx();
          double r = getVdwRadius(atom);
          double q = charges[idx];

          if (r > maxR) {
              maxR = r;
              maxRAtomIdx = idx;
          }
          if (q < 0 && q < minNegQ) {
              minNegQ = q;
              minNegQAtomIdx = idx;
          }
     }

     if (maxRAtomIdx == -1 || minNegQAtomIdx == -1) return 0; // If no max R or no neg charge atom found
     if (maxRAtomIdx == minNegQAtomIdx) return 0; // Same atom


    return getShortestPath(rdkMol, maxRAtomIdx, minNegQAtomIdx); // Returns -1 if not connected
}

AvgNeighborWeightRqDescriptor::AvgNeighborWeightRqDescriptor()
    : Descriptor("AvgNeighborWeightRq", "Avg neighbor weight (VdW Radius * abs(Charge)) per heavy atom") {}

std::variant<double, int, std::string> AvgNeighborWeightRqDescriptor::calculate(const Molecule& mol) const {
     if (!mol.isValid() || !mol.getMolecule()) return 0.0; // Changed NaN to 0.0
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
     if (!rdkMol) return 0.0; // Changed NaN to 0.0
     auto heavyAtoms = getHeavyAtoms(rdkMol);
     if (heavyAtoms.empty()) return 0.0;

    std::vector<double> charges;
    if (!computeAndGetGasteigerCharges(rdkMol, charges)) {
        return "Error: Gasteiger";
    }

    double totalAvgNeighborWeight = 0.0;

     for (const RDKit::Atom* atom : heavyAtoms) {
         double sumNeighborWeight = 0.0;
         int neighborCount = 0;
          for (const auto& nbr : rdkMol->atomNeighbors(atom)) {
              if (nbr->getAtomicNum() > 1) { // Consider only heavy neighbors
                  sumNeighborWeight += getVdwRadius(nbr) * std::abs(charges[nbr->getIdx()]);
                  neighborCount++;
              }
          }
          if (neighborCount > 0) {
              totalAvgNeighborWeight += sumNeighborWeight / neighborCount;
          }
     }

    return std::isfinite(totalAvgNeighborWeight / heavyAtoms.size()) ? totalAvgNeighborWeight / heavyAtoms.size() : 0.0;
}

TopoAutocorrRqDist2Descriptor::TopoAutocorrRqDist2Descriptor()
    : Descriptor("TopoAutocorrRqDist2", "Topological Autocorrelation (R*|q|) at distance 2") {}

std::variant<double, int, std::string> TopoAutocorrRqDist2Descriptor::calculate(const Molecule& mol) const {
     if (!mol.isValid() || !mol.getMolecule()) return 0.0; // Changed NaN to 0.0
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
     if (!rdkMol) return 0.0; // Changed NaN to 0.0
     auto heavyAtoms = getHeavyAtoms(rdkMol);
     if (heavyAtoms.size() < 2) return 0.0;


    std::vector<double> charges;
    if (!computeAndGetGasteigerCharges(rdkMol, charges)) {
        return "Error: Gasteiger";
    }

    double autocorrSum = 0.0;
    std::vector<double> RqValues(rdkMol->getNumAtoms());
     for (const RDKit::Atom* atom : heavyAtoms) {
         RqValues[atom->getIdx()] = getVdwRadius(atom) * std::abs(charges[atom->getIdx()]);
     }


     for (size_t i = 0; i < heavyAtoms.size(); ++i) {
         for (size_t j = i + 1; j < heavyAtoms.size(); ++j) {
              int dist = getShortestPath(rdkMol, heavyAtoms[i]->getIdx(), heavyAtoms[j]->getIdx());
              if (dist == 2) {
                  autocorrSum += RqValues[heavyAtoms[i]->getIdx()] * RqValues[heavyAtoms[j]->getIdx()];
              }
         }
     }
     return std::isfinite(autocorrSum) ? autocorrSum : 0.0;
}

TopoAutocorrRqDist3Descriptor::TopoAutocorrRqDist3Descriptor()
    : Descriptor("TopoAutocorrRqDist3", "Topological Autocorrelation (R*|q|) at distance 3") {}

std::variant<double, int, std::string> TopoAutocorrRqDist3Descriptor::calculate(const Molecule& mol) const {
     if (!mol.isValid() || !mol.getMolecule()) return 0.0; // Changed NaN to 0.0
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
     if (!rdkMol) return 0.0; // Changed NaN to 0.0
     auto heavyAtoms = getHeavyAtoms(rdkMol);
     if (heavyAtoms.size() < 2) return 0.0;


    std::vector<double> charges;
    if (!computeAndGetGasteigerCharges(rdkMol, charges)) {
        return "Error: Gasteiger";
    }

    double autocorrSum = 0.0;
    std::vector<double> RqValues(rdkMol->getNumAtoms());
     for (const RDKit::Atom* atom : heavyAtoms) {
         RqValues[atom->getIdx()] = getVdwRadius(atom) * std::abs(charges[atom->getIdx()]);
     }


     for (size_t i = 0; i < heavyAtoms.size(); ++i) {
         for (size_t j = i + 1; j < heavyAtoms.size(); ++j) {
              int dist = getShortestPath(rdkMol, heavyAtoms[i]->getIdx(), heavyAtoms[j]->getIdx());
              if (dist == 3) {
                  autocorrSum += RqValues[heavyAtoms[i]->getIdx()] * RqValues[heavyAtoms[j]->getIdx()];
              }
         }
     }
     return std::isfinite(autocorrSum) ? autocorrSum : 0.0;
}

TopoAutocorrRqDist4Descriptor::TopoAutocorrRqDist4Descriptor()
    : Descriptor("TopoAutocorrRqDist4", "Topological Autocorrelation (R*|q|) at distance 4") {}

std::variant<double, int, std::string> TopoAutocorrRqDist4Descriptor::calculate(const Molecule& mol) const {
     if (!mol.isValid() || !mol.getMolecule()) return 0.0; // Changed NaN to 0.0
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
     if (!rdkMol) return 0.0; // Changed NaN to 0.0
     auto heavyAtoms = getHeavyAtoms(rdkMol);
      if (heavyAtoms.size() < 2) return 0.0;

    std::vector<double> charges;
    if (!computeAndGetGasteigerCharges(rdkMol, charges)) {
        return "Error: Gasteiger";
    }

    double autocorrSum = 0.0;
    std::vector<double> RqValues(rdkMol->getNumAtoms());
     for (const RDKit::Atom* atom : heavyAtoms) {
         RqValues[atom->getIdx()] = getVdwRadius(atom) * std::abs(charges[atom->getIdx()]);
     }

     for (size_t i = 0; i < heavyAtoms.size(); ++i) {
         for (size_t j = i + 1; j < heavyAtoms.size(); ++j) {
              int dist = getShortestPath(rdkMol, heavyAtoms[i]->getIdx(), heavyAtoms[j]->getIdx());
              if (dist == 4) {
                  autocorrSum += RqValues[heavyAtoms[i]->getIdx()] * RqValues[heavyAtoms[j]->getIdx()];
              }
         }
     }
     return std::isfinite(autocorrSum) ? autocorrSum : 0.0;
}

AvgTopoDistPiWeightedRqDescriptor::AvgTopoDistPiWeightedRqDescriptor()
    : Descriptor("AvgTopoDistPiWeightedRq", "Avg Topo Dist between Pi atoms, weighted by sum(R*|q|)") {}

std::variant<double, int, std::string> AvgTopoDistPiWeightedRqDescriptor::calculate(const Molecule& mol) const {
     if (!mol.isValid() || !mol.getMolecule()) return 0.0; // Changed NaN to 0.0
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
     if (!rdkMol) return 0.0; // Changed NaN to 0.0

    std::vector<double> charges;
    if (!computeAndGetGasteigerCharges(rdkMol, charges)) {
        return "Error: Gasteiger";
    }

     std::vector<const RDKit::Atom*> piAtoms;
     std::vector<double> RqValues(rdkMol->getNumAtoms());
      for (const RDKit::Atom* atom : getHeavyAtoms(rdkMol)) {
          if (isInPiSystem(atom)) {
              piAtoms.push_back(atom);
          }
          RqValues[atom->getIdx()] = getVdwRadius(atom) * std::abs(charges[atom->getIdx()]);
      }

     if (piAtoms.size() < 2) return 0.0;

     double totalWeightedDist = 0.0;
     double totalWeight = 0.0;
     int pairCount = 0;

     for (size_t i = 0; i < piAtoms.size(); ++i) {
         for (size_t j = i + 1; j < piAtoms.size(); ++j) {
              int dist = getShortestPath(rdkMol, piAtoms[i]->getIdx(), piAtoms[j]->getIdx());
              if (dist > 0) {
                   double weight = RqValues[piAtoms[i]->getIdx()] + RqValues[piAtoms[j]->getIdx()];
                   totalWeightedDist += static_cast<double>(dist) * weight;
                   totalWeight += weight;
                   pairCount++;
              }
         }
     }

     if (totalWeight < 1e-9 || pairCount == 0) return 0.0;

     double result = totalWeightedDist / totalWeight;
     return std::isfinite(result) ? result : 0.0;
}

AvgTopoDistLPWeightedRqDescriptor::AvgTopoDistLPWeightedRqDescriptor()
    : Descriptor("AvgTopoDistLPWeightedRq", "Avg Topo Dist between Lone Pair atoms, weighted by sum(R*|q|)") {}

std::variant<double, int, std::string> AvgTopoDistLPWeightedRqDescriptor::calculate(const Molecule& mol) const {
     if (!mol.isValid() || !mol.getMolecule()) return 0.0; // Changed NaN to 0.0
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
     if (!rdkMol) return 0.0; // Changed NaN to 0.0

    std::vector<double> charges;
    if (!computeAndGetGasteigerCharges(rdkMol, charges)) {
        return "Error: Gasteiger";
    }

     std::vector<const RDKit::Atom*> lpAtoms;
     std::vector<double> RqValues(rdkMol->getNumAtoms());
      for (const RDKit::Atom* atom : getHeavyAtoms(rdkMol)) {
          if (hasLonePair(atom)) {
              lpAtoms.push_back(atom);
          }
           RqValues[atom->getIdx()] = getVdwRadius(atom) * std::abs(charges[atom->getIdx()]);
      }

     if (lpAtoms.size() < 2) return 0.0;

     double totalWeightedDist = 0.0;
     double totalWeight = 0.0;
     int pairCount = 0;

     for (size_t i = 0; i < lpAtoms.size(); ++i) {
         for (size_t j = i + 1; j < lpAtoms.size(); ++j) {
              int dist = getShortestPath(rdkMol, lpAtoms[i]->getIdx(), lpAtoms[j]->getIdx());
              if (dist > 0) {
                   double weight = RqValues[lpAtoms[i]->getIdx()] + RqValues[lpAtoms[j]->getIdx()];
                   totalWeightedDist += static_cast<double>(dist) * weight;
                   totalWeight += weight;
                   pairCount++;
              }
         }
     }

     if (totalWeight < 1e-9 || pairCount == 0) return 0.0;

     double result = totalWeightedDist / totalWeight;
     return std::isfinite(result) ? result : 0.0;
}

AvgTopoDistPiLpWeightedRqDescriptor::AvgTopoDistPiLpWeightedRqDescriptor()
    : Descriptor("AvgTopoDistPiLpWeightedRq", "Avg Topo Dist between Pi and Lone Pair atoms, weighted by sum(R*|q|)") {}

std::variant<double, int, std::string> AvgTopoDistPiLpWeightedRqDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0; // Changed NaN to 0.0
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return 0.0; // Changed NaN to 0.0

    std::vector<double> charges;
    if (!computeAndGetGasteigerCharges(rdkMol, charges)) {
        return "Error: Gasteiger";
    }

     std::vector<const RDKit::Atom*> piAtoms;
     std::vector<const RDKit::Atom*> lpAtoms;
     std::vector<double> RqValues(rdkMol->getNumAtoms());

     for (const RDKit::Atom* atom : getHeavyAtoms(rdkMol)) {
          if (isInPiSystem(atom)) piAtoms.push_back(atom);
          if (hasLonePair(atom)) lpAtoms.push_back(atom);
          RqValues[atom->getIdx()] = getVdwRadius(atom) * std::abs(charges[atom->getIdx()]);
     }


     if (piAtoms.empty() || lpAtoms.empty()) return 0.0;

     double totalWeightedDist = 0.0;
     double totalWeight = 0.0;
     int pairCount = 0;

     for (const auto* piAtom : piAtoms) {
         for (const auto* lpAtom : lpAtoms) {
              if (piAtom->getIdx() == lpAtom->getIdx()) continue; // Don't compare same atom
              int dist = getShortestPath(rdkMol, piAtom->getIdx(), lpAtom->getIdx());
              if (dist > 0) {
                   double weight = RqValues[piAtom->getIdx()] + RqValues[lpAtom->getIdx()];
                   totalWeightedDist += static_cast<double>(dist) * weight;
                   totalWeight += weight;
                   pairCount++;
              }
         }
     }

     if (totalWeight < 1e-9 || pairCount == 0) return 0.0;

     double result = totalWeightedDist / totalWeight;
     return std::isfinite(result) ? result : 0.0;
}

PathCountAlternatingRChargeDescriptor::PathCountAlternatingRChargeDescriptor()
    : Descriptor("PathCountAlternatingRCharge", "Count of paths with alternating R and charge patterns") {}

std::variant<double, int, std::string> PathCountAlternatingRChargeDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol || rdkMol->getNumHeavyAtoms() < 3) return 0.0;
    
    std::vector<double> charges;
    if (!computeAndGetGasteigerCharges(rdkMol, charges)) return 0.0;
    
    // Build a graph where alternating R and charge patterns can be analyzed
    std::vector<std::vector<int>> graph(rdkMol->getNumAtoms());
    for (const RDKit::Bond* bond : rdkMol->bonds()) {
        unsigned int begin = bond->getBeginAtomIdx();
        unsigned int end = bond->getEndAtomIdx();
        graph[begin].push_back(end);
        graph[end].push_back(begin);
    }
    
    // Count paths where R and charge alternate (one high, one low)
    int pathCount = 0;
    const double radiusThreshold = 1.5; // Example threshold for VdW radius
    const double chargeThreshold = 0.05; // Example threshold for charge magnitude
    
    for (const RDKit::Atom* atom1 : getHeavyAtoms(rdkMol)) {
        int idx1 = atom1->getIdx();
        bool highR1 = getVdwRadius(atom1) > radiusThreshold;
        bool highQ1 = std::abs(charges[idx1]) > chargeThreshold;
        
        // Use BFS to find paths with alternating properties
        std::vector<bool> visited(rdkMol->getNumAtoms(), false);
        std::queue<std::pair<int, std::pair<bool, bool>>> q; // (index, (high radius, high charge))
        q.push({idx1, {highR1, highQ1}});
        visited[idx1] = true;
        
        while (!q.empty()) {
            auto [currentIdx, properties] = q.front();
            auto [currentHighR, currentHighQ] = properties;
            q.pop();
            
            for (int neighborIdx : graph[currentIdx]) {
                if (visited[neighborIdx]) continue;
                
                const RDKit::Atom* neighbor = rdkMol->getAtomWithIdx(neighborIdx);
                if (neighbor->getAtomicNum() <= 1) continue; // Skip hydrogens
                
                bool neighborHighR = getVdwRadius(neighbor) > radiusThreshold;
                bool neighborHighQ = std::abs(charges[neighborIdx]) > chargeThreshold;
                
                // Check if we have alternating pattern (different from current)
                if (neighborHighR != currentHighR || neighborHighQ != currentHighQ) {
                    pathCount++;
                    q.push({neighborIdx, {neighborHighR, neighborHighQ}});
                    visited[neighborIdx] = true;
                }
            }
        }
    }
    
    return std::isfinite(pathCount) ? pathCount : 0.0;
}

RatioLongestWeightedPathRingChainDescriptor::RatioLongestWeightedPathRingChainDescriptor()
    : Descriptor("RatioLongestWeightedPathRingChain", "Ratio of longest weighted path in rings vs chains") {}

std::variant<double, int, std::string> RatioLongestWeightedPathRingChainDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return 0.0;
    
    const RDKit::RingInfo* ringInfo = rdkMol->getRingInfo();
    if (!ringInfo || !ringInfo->isInitialized()) {
        RDKit::MolOps::findSSSR(*const_cast<RDKit::ROMol*>(rdkMol));
        ringInfo = rdkMol->getRingInfo();
    }
    
    if (!ringInfo) return 0.0;
    
    std::vector<double> charges;
    if (!computeAndGetGasteigerCharges(rdkMol, charges)) return 0.0;
    
    // Find longest weighted path in rings and chains
    double longestRingPath = 0.0;
    double longestChainPath = 0.0;
    
    auto heavyAtoms = getHeavyAtoms(rdkMol);
    for (size_t i = 0; i < heavyAtoms.size(); ++i) {
        for (size_t j = i+1; j < heavyAtoms.size(); ++j) {
            int dist = getShortestPath(rdkMol, heavyAtoms[i]->getIdx(), heavyAtoms[j]->getIdx());
            if (dist <= 0) continue;
            
            double weight1 = getVdwRadius(heavyAtoms[i]) * std::abs(charges[heavyAtoms[i]->getIdx()]);
            double weight2 = getVdwRadius(heavyAtoms[j]) * std::abs(charges[heavyAtoms[j]->getIdx()]);
            double weightedPath = dist * weight1 * weight2;
            
            bool atom1InRing = ringInfo->numAtomRings(heavyAtoms[i]->getIdx()) > 0;
            bool atom2InRing = ringInfo->numAtomRings(heavyAtoms[j]->getIdx()) > 0;
            
            if (atom1InRing && atom2InRing) {
                longestRingPath = std::max(longestRingPath, weightedPath);
            } else if (!atom1InRing && !atom2InRing) {
                longestChainPath = std::max(longestChainPath, weightedPath);
            }
        }
    }
    
    if (longestChainPath < 1e-9) {
        return (longestRingPath > 1e-9) ? 9999.9 : 1.0;
    }
    
    double result = longestRingPath / longestChainPath;
    return std::isfinite(result) ? result : 0.0;
}

SumWeightRChargeDegreeGt3Descriptor::SumWeightRChargeDegreeGt3Descriptor()
    : Descriptor("SumWeightRChargeDegreeGt3", "Sum of (R+|q|) for atoms with degree > 3") {}

std::variant<double, int, std::string> SumWeightRChargeDegreeGt3Descriptor::calculate(const Molecule& mol) const {
     if (!mol.isValid() || !mol.getMolecule()) return 0.0; // Changed NaN to 0.0
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
     if (!rdkMol) return 0.0; // Changed NaN to 0.0

    std::vector<double> charges;
    if (!computeAndGetGasteigerCharges(rdkMol, charges)) {
        return "Error: Gasteiger";
    }

    double sumVal = 0.0;
    for (const RDKit::Atom* atom : getHeavyAtoms(rdkMol)) {
        if (atom->getDegree() > 3) {
            sumVal += getVdwRadius(atom) + std::abs(charges[atom->getIdx()]);
        }
    }
    return std::isfinite(sumVal) ? sumVal : 0.0;
}

RatioAvgRWeightedChargeTerminalDescriptor::RatioAvgRWeightedChargeTerminalDescriptor()
    : Descriptor("RatioAvgRWeightedChargeTerminal", "Ratio of Avg(R*|q|) for terminal vs non-terminal heavy atoms") {}

std::variant<double, int, std::string> RatioAvgRWeightedChargeTerminalDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0; // Changed NaN to 0.0
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return 0.0; // Changed NaN to 0.0

    std::vector<double> charges;
    if (!computeAndGetGasteigerCharges(rdkMol, charges)) {
        return "Error: Gasteiger";
    }

    double sumRqTerminal = 0.0;
    int countTerminal = 0;
    double sumRqInternal = 0.0;
    int countInternal = 0;

    for (const RDKit::Atom* atom : getHeavyAtoms(rdkMol)) {
         double val = getVdwRadius(atom) * std::abs(charges[atom->getIdx()]);
         if (atom->getDegree() <= 1) {
             sumRqTerminal += val;
             countTerminal++;
         } else {
             sumRqInternal += val;
             countInternal++;
         }
    }

     if(countTerminal == 0 && countInternal == 0) return 1.0; // No heavy atoms
     if(countInternal == 0) return (countTerminal > 0) ? 9999.9 : 1.0; // All terminal or none
     if(countTerminal == 0) return 0.0; // All internal

     double avgTerminal = sumRqTerminal / countTerminal;
     double avgInternal = sumRqInternal / countInternal;

      if (std::abs(avgInternal) < 1e-9) {
          return (std::abs(avgTerminal) > 1e-9) ? 9999.9 : 1.0; // Keep large num/1.0
      }

    double result = avgTerminal / avgInternal;
    return std::isfinite(result) ? result : 0.0;
}

EigenvalueWeightedConnectivityDescriptor::EigenvalueWeightedConnectivityDescriptor()
    : Descriptor("EigenvalueWeightedConnectivity", "Eigenvalue of connectivity matrix weighted by R*|q|") {}

std::variant<double, int, std::string> EigenvalueWeightedConnectivityDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol || rdkMol->getNumHeavyAtoms() < 2) return 0.0;
    
    std::vector<double> charges;
    if (!computeAndGetGasteigerCharges(rdkMol, charges)) return 0.0;
    
    auto heavyAtoms = getHeavyAtoms(rdkMol);
    size_t n = heavyAtoms.size();
    
    // Calculate R*|q| values for each atom
    std::vector<double> RqValues(rdkMol->getNumAtoms(), 0.0);
    for (size_t i = 0; i < n; ++i) {
        const RDKit::Atom* atom = heavyAtoms[i];
        RqValues[atom->getIdx()] = getVdwRadius(atom) * std::abs(charges[atom->getIdx()]);
    }
    
    // Create weighted adjacency matrix
    std::vector<std::vector<double>> adjMatrix(n, std::vector<double>(n, 0.0));
    for (size_t i = 0; i < n; ++i) {
        const RDKit::Atom* atom1 = heavyAtoms[i];
        for (size_t j = 0; j < n; ++j) {
            if (i == j) continue;
            const RDKit::Atom* atom2 = heavyAtoms[j];
            
            const RDKit::Bond* bond = rdkMol->getBondBetweenAtoms(atom1->getIdx(), atom2->getIdx());
            if (bond) {
                adjMatrix[i][j] = std::sqrt(RqValues[atom1->getIdx()] * RqValues[atom2->getIdx()]);
            }
        }
    }
    
    // Calculate eigenvalue using power iteration method
    std::vector<double> x(n, 1.0); // Initial guess
    std::vector<double> x_new(n, 0.0);
    double lambda = 0.0;
    
    const int MAX_ITER = 100;
    const double TOLERANCE = 1e-6;
    
    for (int iter = 0; iter < MAX_ITER; ++iter) {
        // Matrix-vector multiplication
        for (size_t i = 0; i < n; ++i) {
            x_new[i] = 0.0;
            for (size_t j = 0; j < n; ++j) {
                x_new[i] += adjMatrix[i][j] * x[j];
            }
        }
        
        // Calculate L2 norm
        double norm = 0.0;
        for (double val : x_new) {
            norm += val * val;
        }
        norm = std::sqrt(norm);
        
        if (norm < 1e-10) return 0.0; // Zero vector
        
        // Normalize
        for (size_t i = 0; i < n; ++i) {
            x_new[i] /= norm;
        }
        
        // Estimate eigenvalue (Rayleigh quotient)
        double lambda_new = 0.0;
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < n; ++j) {
                lambda_new += x_new[i] * adjMatrix[i][j] * x_new[j];
            }
        }
        
        // Check convergence
        if (std::abs(lambda_new - lambda) < TOLERANCE) {
            lambda = lambda_new;
            break;
        }
        
        lambda = lambda_new;
        x = x_new;
    }
    
    return std::isfinite(lambda) ? lambda : 0.0;
}

WeightedBranchPointComplexityDescriptor::WeightedBranchPointComplexityDescriptor()
    : Descriptor("WeightedBranchPointComplexity", "Complexity measure of branch points weighted by R*|q|") {}

std::variant<double, int, std::string> WeightedBranchPointComplexityDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return 0.0;
    
    std::vector<double> charges;
    if (!computeAndGetGasteigerCharges(rdkMol, charges)) return 0.0;
    
    std::vector<double> RqValues(rdkMol->getNumAtoms());
    for (const RDKit::Atom* atom : getHeavyAtoms(rdkMol)) {
        RqValues[atom->getIdx()] = getVdwRadius(atom) * std::abs(charges[atom->getIdx()]);
    }
    
    double totalComplexity = 0.0;
    
    for (const RDKit::Atom* atom : getHeavyAtoms(rdkMol)) {
        unsigned int degree = atom->getDegree();
        if (degree >= 3) {  // Branch point
            // Calculate complexity as product of degree and R*|q| value
            double branchComplexity = degree * RqValues[atom->getIdx()];
            
            // Add neighbor information (sum of neighbor R*|q| values)
            double neighborSum = 0.0;
            for (const auto& nbr : rdkMol->atomNeighbors(atom)) {
                neighborSum += RqValues[nbr->getIdx()];
            }
            
            totalComplexity += branchComplexity * (1.0 + 0.1 * neighborSum);
        }
    }
    
    return std::isfinite(totalComplexity) ? totalComplexity : 0.0;
}

ShannonEntropyRChargeDescriptor::ShannonEntropyRChargeDescriptor()
    : Descriptor("ShannonEntropyRCharge", "Shannon entropy of R*|q| distribution") {}

std::variant<double, int, std::string> ShannonEntropyRChargeDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return 0.0;
    
    std::vector<const RDKit::Atom*> heavyAtoms = getHeavyAtoms(rdkMol);
    if (heavyAtoms.empty()) return 0.0;
    
    std::vector<double> charges;
    if (!computeAndGetGasteigerCharges(rdkMol, charges)) return 0.0;
    
    std::vector<double> RqValues;
    RqValues.reserve(heavyAtoms.size());
    double totalRqSum = 0.0;
    
    for (const RDKit::Atom* atom : heavyAtoms) {
        double val = getVdwRadius(atom) * std::abs(charges[atom->getIdx()]);
        RqValues.push_back(val);
        totalRqSum += val;
    }
    
    if (totalRqSum < 1e-9) return 0.0;
    
    // Calculate probabilities
    std::vector<double> probabilities;
    probabilities.reserve(RqValues.size());
    for (double val : RqValues) {
        probabilities.push_back(val / totalRqSum);
    }
    
    // Calculate Shannon entropy
    double entropy = 0.0;
    for (double p : probabilities) {
        if (p > 1e-9) {
            entropy -= p * std::log2(p);
        }
    }
    
    return std::isfinite(entropy) ? entropy : 0.0;
}

AvgRadiusPlusChargePerDegreeDescriptor::AvgRadiusPlusChargePerDegreeDescriptor()
    : Descriptor("AvgRadiusPlusChargePerDegree", "Average (VdW Radius + abs(Charge)) / Degree for heavy atoms") {}

std::variant<double, int, std::string> AvgRadiusPlusChargePerDegreeDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0; // Changed NaN to 0.0
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return 0.0; // Changed NaN to 0.0
     auto heavyAtoms = getHeavyAtoms(rdkMol);
     if (heavyAtoms.empty()) return 0.0;


    std::vector<double> charges;
    if (!computeAndGetGasteigerCharges(rdkMol, charges)) {
        return "Error: Gasteiger";
    }

    double sumVal = 0.0;
    double epsilon = 1e-6; // Prevent division by zero degree for isolated atoms

    for (const RDKit::Atom* atom : heavyAtoms) {
         double degree = static_cast<double>(atom->getDegree());
         sumVal += (getVdwRadius(atom) + std::abs(charges[atom->getIdx()])) / (degree + epsilon) ;
    }

    return std::isfinite(sumVal / heavyAtoms.size()) ? sumVal / heavyAtoms.size() : 0.0;
}

SumDeltaRadiusChargeWeightedDescriptor::SumDeltaRadiusChargeWeightedDescriptor()
    : Descriptor("SumDeltaRadiusChargeWeighted", "Sum of (VdW Radius - Covalent Radius) * abs(Charge)") {}

std::variant<double, int, std::string> SumDeltaRadiusChargeWeightedDescriptor::calculate(const Molecule& mol) const {
     if (!mol.isValid() || !mol.getMolecule()) return 0.0; // Changed NaN to 0.0
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
     if (!rdkMol) return 0.0; // Changed NaN to 0.0

    std::vector<double> charges;
    if (!computeAndGetGasteigerCharges(rdkMol, charges)) {
        return "Error: Gasteiger";
    }

    double sumVal = 0.0;
     for (const RDKit::Atom* atom : getHeavyAtoms(rdkMol)) {
          double deltaR = getVdwRadius(atom) - getCovalentRadius(atom);
          sumVal += deltaR * std::abs(charges[atom->getIdx()]);
     }
    return std::isfinite(sumVal) ? sumVal : 0.0;
}

WeightedBuriedAtomCountDescriptor::WeightedBuriedAtomCountDescriptor()
    : Descriptor("WeightedBuriedAtomCount", "Placeholder: Count of 'buried' atoms weighted by R*|q|") {}

std::variant<double, int, std::string> WeightedBuriedAtomCountDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return 0.0;
    
    std::vector<double> charges;
    if (!computeAndGetGasteigerCharges(rdkMol, charges)) return 0.0;
    
    // Define a buried atom as one with at least 3 heavy atom neighbors
    double buriedAtomSum = 0.0;
    
    for (const RDKit::Atom* atom : getHeavyAtoms(rdkMol)) {
        int heavyNeighborCount = 0;
        for (const auto& nbr : rdkMol->atomNeighbors(atom)) {
            if (nbr->getAtomicNum() > 1) heavyNeighborCount++;
        }
        
        if (heavyNeighborCount >= 3) {
            double weight = getVdwRadius(atom) * std::abs(charges[atom->getIdx()]);
            buriedAtomSum += weight;
        }
    }
    
    return std::isfinite(buriedAtomSum) ? buriedAtomSum : 0.0;
}

RatioSumRqFormalChargeDescriptor::RatioSumRqFormalChargeDescriptor()
    : Descriptor("RatioSumRqFormalCharge", "Ratio of Sum(R*|q|) for formally charged vs uncharged heavy atoms") {}

std::variant<double, int, std::string> RatioSumRqFormalChargeDescriptor::calculate(const Molecule& mol) const {
     if (!mol.isValid() || !mol.getMolecule()) return 0.0; // Changed NaN to 0.0
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
     if (!rdkMol) return 0.0; // Changed NaN to 0.0

    std::vector<double> charges;
    if (!computeAndGetGasteigerCharges(rdkMol, charges)) {
        return "Error: Gasteiger";
    }

     double sumRqCharged = 0.0;
     double sumRqUncharged = 0.0;

      for (const RDKit::Atom* atom : getHeavyAtoms(rdkMol)) {
           double val = getVdwRadius(atom) * std::abs(charges[atom->getIdx()]);
           if (atom->getFormalCharge() != 0) {
               sumRqCharged += val;
           } else {
               sumRqUncharged += val;
           }
      }

      if (sumRqUncharged < 1e-9) {
          return (sumRqCharged > 1e-9) ? 9999.9 : 1.0; // Keep large num / 1.0
      }

    double result = sumRqCharged / sumRqUncharged;
    return std::isfinite(result) ? result : 0.0;
}

AvgRadiusHBDonorWeightedChargeDescriptor::AvgRadiusHBDonorWeightedChargeDescriptor()
    : Descriptor("AvgRadiusHBDonorWeightedCharge", "Avg VdW Radius of HBD heavy atoms, weighted by their charge") {}

std::variant<double, int, std::string> AvgRadiusHBDonorWeightedChargeDescriptor::calculate(const Molecule& mol) const {
     if (!mol.isValid() || !mol.getMolecule()) return 0.0; // Changed NaN to 0.0
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
     if (!rdkMol) return 0.0; // Changed NaN to 0.0

    std::vector<double> charges;
    if (!computeAndGetGasteigerCharges(rdkMol, charges)) {
        return "Error: Gasteiger";
    }

    double weightedSumRadius = 0.0;
    double sumWeights = 0.0; // Sum of |charges|
    int count = 0;

     for (const RDKit::Atom* atom : getHeavyAtoms(rdkMol)) {
         if (isHBondDonor(atom)) { // Check if the heavy atom is a donor // Corrected: use function from anon namespace
              double charge = charges[atom->getIdx()];
              double absCharge = std::abs(charge);
              weightedSumRadius += getVdwRadius(atom) * absCharge;
              sumWeights += absCharge;
              count++;
         }
     }

      if (count == 0 || sumWeights < 1e-9) return 0.0;

    double result = weightedSumRadius / sumWeights;
    return std::isfinite(result) ? result : 0.0;
}

AvgRadiusHBAcceptorWeightedChargeDescriptor::AvgRadiusHBAcceptorWeightedChargeDescriptor()
    : Descriptor("AvgRadiusHBAcceptorWeightedCharge", "Avg VdW Radius of HBA heavy atoms, weighted by their charge") {}

std::variant<double, int, std::string> AvgRadiusHBAcceptorWeightedChargeDescriptor::calculate(const Molecule& mol) const {
     if (!mol.isValid() || !mol.getMolecule()) return 0.0; // Changed NaN to 0.0
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
     if (!rdkMol) return 0.0; // Changed NaN to 0.0

    std::vector<double> charges;
    if (!computeAndGetGasteigerCharges(rdkMol, charges)) {
        return "Error: Gasteiger";
    }

    double weightedSumRadius = 0.0;
    double sumWeights = 0.0; // Sum of |charges|
    int count = 0;

     for (const RDKit::Atom* atom : getHeavyAtoms(rdkMol)) {
         if (isHBondAcceptor(atom)) { // Check if the heavy atom is an acceptor // Corrected: use function from anon namespace
             double charge = charges[atom->getIdx()];
              double absCharge = std::abs(charge);
              weightedSumRadius += getVdwRadius(atom) * absCharge;
              sumWeights += absCharge;
              count++;
         }
     }

     if (count == 0 || sumWeights < 1e-9) return 0.0;

    double result = weightedSumRadius / sumWeights;
    return std::isfinite(result) ? result : 0.0;
}

RatioSumRPolarNonpolarFragDescriptor::RatioSumRPolarNonpolarFragDescriptor()
    : Descriptor("RatioSumRPolarNonpolarFrag", "Placeholder: Ratio Sum(R) polar frags / Sum(R) nonpolar frags") {}

std::variant<double, int, std::string> RatioSumRPolarNonpolarFragDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) return 0.0;
    
    double sumRPolar = 0.0;
    double sumRNonpolar = 0.0;
    
    for (const RDKit::Atom* atom : getHeavyAtoms(rdkMol)) {
        int atomicNum = atom->getAtomicNum();
        double radius = getVdwRadius(atom);
        
        // Define polar atoms (N, O, S, P, halogens)
        if (atomicNum == 7 || atomicNum == 8 || atomicNum == 16 || 
            atomicNum == 15 || atomicNum == 9 || atomicNum == 17 || 
            atomicNum == 35 || atomicNum == 53) {
            sumRPolar += radius;
        } else {
            sumRNonpolar += radius;
        }
    }
    
    if (sumRNonpolar < 1e-9) {
        return (sumRPolar > 1e-9) ? 9999.9 : 1.0;
    }
    
    double result = sumRPolar / sumRNonpolar;
    return std::isfinite(result) ? result : 0.0;
}

SumRadiusSp2CWeightedChargeDescriptor::SumRadiusSp2CWeightedChargeDescriptor()
    : Descriptor("SumRadiusSp2CWeightedCharge", "Sum of VdW Radius for sp2 carbons, weighted by Gasteiger charge") {}

std::variant<double, int, std::string> SumRadiusSp2CWeightedChargeDescriptor::calculate(const Molecule& mol) const {
     if (!mol.isValid() || !mol.getMolecule()) return 0.0; // Changed NaN to 0.0
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
     if (!rdkMol) return 0.0; // Changed NaN to 0.0

    std::vector<double> charges;
    if (!computeAndGetGasteigerCharges(rdkMol, charges)) {
        return "Error: Gasteiger";
    }

    double sumVal = 0.0;
    for (const RDKit::Atom* atom : rdkMol->atoms()) {
        if (atom->getAtomicNum() == 6 && atom->getHybridization() == RDKit::Atom::HybridizationType::SP2) {
            sumVal += getVdwRadius(atom) * charges[atom->getIdx()]; // Use actual charge, not abs
        }
    }
    return std::isfinite(sumVal) ? sumVal : 0.0;
}

SumRadiusSp3CWeightedChargeDescriptor::SumRadiusSp3CWeightedChargeDescriptor()
    : Descriptor("SumRadiusSp3CWeightedCharge", "Sum of VdW Radius for sp3 carbons, weighted by Gasteiger charge") {}

std::variant<double, int, std::string> SumRadiusSp3CWeightedChargeDescriptor::calculate(const Molecule& mol) const {
     if (!mol.isValid() || !mol.getMolecule()) return 0.0; // Changed NaN to 0.0
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
     if (!rdkMol) return 0.0; // Changed NaN to 0.0

    std::vector<double> charges;
    if (!computeAndGetGasteigerCharges(rdkMol, charges)) {
        return "Error: Gasteiger";
    }

    double sumVal = 0.0;
    for (const RDKit::Atom* atom : rdkMol->atoms()) {
        if (atom->getAtomicNum() == 6 && atom->getHybridization() == RDKit::Atom::HybridizationType::SP3) {
             sumVal += getVdwRadius(atom) * charges[atom->getIdx()]; // Use actual charge
        }
    }
    return std::isfinite(sumVal) ? sumVal : 0.0;
}

AvgNeighborChargeRadiusPerDegreeDescriptor::AvgNeighborChargeRadiusPerDegreeDescriptor()
    : Descriptor("AvgNeighborChargeRadiusPerDegree", "Avg Ratio: (R_heavy * Avg |q|_neighbors) / Degree") {}

std::variant<double, int, std::string> AvgNeighborChargeRadiusPerDegreeDescriptor::calculate(const Molecule& mol) const {
     if (!mol.isValid() || !mol.getMolecule()) return 0.0; // Changed NaN to 0.0
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
     if (!rdkMol) return 0.0; // Changed NaN to 0.0
      auto heavyAtoms = getHeavyAtoms(rdkMol);
      if (heavyAtoms.empty()) return 0.0;


    std::vector<double> charges;
    if (!computeAndGetGasteigerCharges(rdkMol, charges)) {
        return "Error: Gasteiger";
    }

    double totalRatioSum = 0.0;
    double epsilon = 1e-6;

     for (const RDKit::Atom* atom : heavyAtoms) {
         double sumNeighborAbsQ = 0.0;
         int neighborCount = 0;
          for (const auto& nbr : rdkMol->atomNeighbors(atom)) {
               // Consider all neighbors for charge average? Or just heavy? Let's use all for now.
               sumNeighborAbsQ += std::abs(charges[nbr->getIdx()]);
               neighborCount++;
          }

          double avgNeighborAbsQ = (neighborCount > 0) ? sumNeighborAbsQ / neighborCount : 0.0;
          double degree = static_cast<double>(atom->getDegree());
          double radius = getVdwRadius(atom);

          totalRatioSum += (radius * avgNeighborAbsQ) / (degree + epsilon);
     }


    double result = totalRatioSum / heavyAtoms.size();
    return std::isfinite(result) ? result : 0.0;
}

KierShapeIndexVariant3Descriptor::KierShapeIndexVariant3Descriptor()
    : Descriptor("KierShapeIndexVariant3", "Placeholder: Kier Kappa3 index variant using R*|q|") {}

std::variant<double, int, std::string> KierShapeIndexVariant3Descriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) return 0.0;
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol || rdkMol->getNumHeavyAtoms() < 3) return 0.0;
    
    std::vector<double> charges;
    if (!computeAndGetGasteigerCharges(rdkMol, charges)) return 0.0;
    
    // Calculate alpha values (weighted by R*|q|)
    std::vector<double> alphaValues;
    double alphaSum = 0.0;
    
    for (const RDKit::Atom* atom : getHeavyAtoms(rdkMol)) {
        double weight = getVdwRadius(atom) * std::abs(charges[atom->getIdx()]);
        double alpha = weight / static_cast<double>(atom->getDegree() + 1);
        alphaValues.push_back(alpha);
        alphaSum += alpha;
    }
    
    // Calculate P3 (number of paths of length 3)
    int p3Count = 0;
    for (size_t i = 0; i < rdkMol->getNumAtoms(); ++i) {
        for (size_t j = 0; j < rdkMol->getNumAtoms(); ++j) {
            if (i == j) continue;
            int dist = getShortestPath(rdkMol, i, j);
            if (dist == 3) p3Count++;
        }
    }
    p3Count /= 2; // Each path counted twice
    
    // Calculate A3 (modified alpha sum for paths of length 3)
    double a3Sum = 0.0;
    auto heavyAtoms = getHeavyAtoms(rdkMol);
    for (size_t i = 0; i < heavyAtoms.size(); ++i) {
        for (size_t j = i+1; j < heavyAtoms.size(); ++j) {
            int dist = getShortestPath(rdkMol, heavyAtoms[i]->getIdx(), heavyAtoms[j]->getIdx());
            if (dist == 3) {
                a3Sum += alphaValues[i] * alphaValues[j];
            }
        }
    }
    
    if (p3Count == 0) return 0.0;
    double result = a3Sum / static_cast<double>(p3Count);
    return std::isfinite(result) ? result : 0.0;
}

} // namespace descriptors
} // namespace desfact