#pragma once

#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace desfact {
namespace descriptors {
namespace ElementProperties {

    // Pauling electronegativity values
    static const std::unordered_map<int, double> electronegativity = {
        {1, 2.20}, {6, 2.55}, {7, 3.04}, {8, 3.44}, {9, 3.98},
        {15, 2.19}, {16, 2.58}, {17, 3.16}, {35, 2.96}, {53, 2.66}
    };

    // Covalent radii in Angstroms
    static const std::unordered_map<int, double> covalentRadius = {
        {1, 0.31}, {6, 0.76}, {7, 0.71}, {8, 0.66}, {9, 0.57},
        {15, 1.07}, {16, 1.05}, {17, 1.02}, {35, 1.20}, {53, 1.39}
    };

    // Atomic polarizability in cubic Angstroms
    static const std::unordered_map<int, double> atomicPolarizability = {
        {1, 0.667}, {6, 1.76}, {7, 1.10}, {8, 0.802}, {9, 0.557},
        {15, 3.63}, {16, 2.90}, {17, 2.18}, {35, 3.05}, {53, 5.35}
    };

    // Ionization energy in eV
    static const std::unordered_map<int, double> ionizationEnergy = {
        {1, 13.6}, {6, 11.3}, {7, 14.5}, {8, 13.6}, {9, 17.4},
        {15, 10.5}, {16, 10.4}, {17, 13.0}, {35, 11.8}, {53, 10.5}
    };

    // Electron affinity in eV
    static const std::unordered_map<int, double> electronAffinity = {
        {1, 0.75}, {6, 1.26}, {7, -0.07}, {8, 1.46}, {9, 3.40},
        {15, 0.75}, {16, 2.08}, {17, 3.62}, {35, 3.36}, {53, 3.06}
    };

    // Van der Waals volume in cubic Angstroms
    static const std::unordered_map<int, double> vanDerWaalsVolume = {
        {1, 7.24}, {6, 20.58}, {7, 15.60}, {8, 14.71}, {9, 13.31},
        {15, 24.43}, {16, 24.43}, {17, 22.45}, {35, 27.07}, {53, 35.19}
    };

    // Van der Waals radius in Angstroms (Bondi values)
    static const std::unordered_map<int, double> vanDerWaalsRadius = {
        {1, 1.20}, {6, 1.70}, {7, 1.55}, {8, 1.52}, {9, 1.47}, {14, 2.10},
        {15, 1.80}, {16, 1.80}, {17, 1.75}, {35, 1.85}, {53, 1.98}
    };

    // Common oxidation states
    static const std::unordered_map<int, std::vector<int>> oxidationStates = {
        {1, {1, -1}}, {6, {4, 2, -4}}, {7, {5, 3, -3}}, {8, {-2}}, {9, {-1}},
        {15, {5, 3, -3}}, {16, {6, 4, 2, -2}}, {17, {7, 5, 3, 1, -1}},
        {35, {7, 5, 3, 1, -1}}, {53, {7, 5, 3, 1, -1}}
    };

    // Metal elements (simplified list)
    static const std::unordered_set<int> metals = {
        3, 4, 11, 12, 13, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31,
        37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50,
        55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71,
        72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84,
        87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103
    };

    // Metalloid elements
    static const std::unordered_set<int> metalloids = {
        5, 14, 32, 33, 51, 52, 84
    };

    // Get element property with default value
    template<typename T>
    T getProperty(const std::unordered_map<int, T>& propertyMap, int atomicNum, T defaultValue) {
        auto it = propertyMap.find(atomicNum);
        return (it != propertyMap.end()) ? it->second : defaultValue;
    }

} // namespace ElementProperties
} // namespace descriptors
} // namespace desfact 