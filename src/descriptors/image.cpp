#include "descriptors/image.hpp"
#include <GraphMol/GraphMol.h>
#include <GraphMol/Depictor/RDDepictor.h>
#include <GraphMol/AtomIterators.h>
#include <GraphMol/BondIterators.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/MolTransforms/MolTransforms.h>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <set>
#include <functional>
#include "utils.hpp"

namespace desfact {
namespace descriptors {

// Helper functions
namespace {

double getVdwRadius(int atomicNum) {
    switch (atomicNum) {
        case 1:  return 1.20; // H
        case 6:  return 1.70; // C
        case 7:  return 1.55; // N
        case 8:  return 1.52; // O
        case 9:  return 1.47; // F
        case 15: return 1.80; // P
        case 16: return 1.80; // S
        case 17: return 1.75; // Cl
        case 35: return 1.85; // Br
        case 53: return 1.98; // I
        default: return 1.70; // Default to C
    }
}

// Get RGB color for atom (using common convention)
std::tuple<int, int, int> getAtomColor(int atomicNum) {
    if (atomicNum == 1) return {255, 255, 255}; // H
    if (atomicNum == 6) return {51, 51, 51};   // C
    if (atomicNum == 7) return {48, 80, 248};  // N
    if (atomicNum == 8) return {255, 13, 13};  // O
    if (atomicNum == 9) return {144, 224, 80}; // F
    if (atomicNum == 15) return {255, 128, 0}; // P
    if (atomicNum == 16) return {255, 255, 48}; // S
    if (atomicNum == 17) return {31, 240, 31}; // Cl
    if (atomicNum == 35) return {166, 41, 41}; // Br
    if (atomicNum == 53) return {148, 0, 148}; // I
    return {128, 128, 128}; // Default
}

double getLuminance(const std::tuple<int, int, int>& rgb) {
    double r = std::get<0>(rgb) / 255.0;
    double g = std::get<1>(rgb) / 255.0;
    double b = std::get<2>(rgb) / 255.0;
    return 0.2126 * r + 0.7152 * g + 0.0722 * b;
}

// Helper: Find all simple cycles of size 5 or 6 (for aromatic-like rings)
std::vector<std::vector<int>> findSmallRings(
    const std::vector<std::vector<int>>& adjList, int maxRingSize = 6) {
    std::vector<std::vector<int>> rings;
    int n = adjList.size();
    std::set<std::set<int>> uniqueRings;
    std::function<void(int,int,std::vector<int>&,std::set<int>&)> dfs =
        [&](int start, int curr, std::vector<int>& path, std::set<int>& visited) {
            if (path.size() > maxRingSize) return;
            for (int nb : adjList[curr]) {
                if (nb == start && path.size() >= 5) {
                    std::set<int> ringSet(path.begin(), path.end());
                    if (uniqueRings.insert(ringSet).second)
                        rings.push_back(path);
                } else if (!visited.count(nb)) {
                    visited.insert(nb);
                    path.push_back(nb);
                    dfs(start, nb, path, visited);
                    path.pop_back();
                    visited.erase(nb);
                }
            }
        };
    for (int i = 0; i < n; ++i) {
        std::vector<int> path = {i};
        std::set<int> visited = {i};
        dfs(i, i, path, visited);
    }
    return rings;
}

// Calculate the convex hull using Graham scan
std::vector<std::pair<double, double>> calculateConvexHull(const std::vector<double>& x, const std::vector<double>& y) {
    struct Point {
        double x, y;
        bool operator<(const Point& p) const {
            return x < p.x || (x == p.x && y < p.y);
        }
    };
    
    int n = x.size();
    std::vector<Point> points(n);
    for (int i = 0; i < n; i++) {
        points[i] = {x[i], y[i]};
    }
    
    // Sort points lexicographically
    std::sort(points.begin(), points.end());
    
    std::vector<Point> hull;
    // Build lower hull
    for (int i = 0; i < n; i++) {
        while (hull.size() >= 2) {
            Point p1 = hull[hull.size() - 2];
            Point p2 = hull[hull.size() - 1];
            Point p3 = points[i];
            double cross = (p2.x - p1.x) * (p3.y - p1.y) - (p2.y - p1.y) * (p3.x - p1.x);
            if (cross <= 0) hull.pop_back();
            else break;
        }
        hull.push_back(points[i]);
    }
    
    // Build upper hull
    size_t lower_hull_size = hull.size();
    for (int i = n - 2; i >= 0; i--) {
        while (hull.size() > lower_hull_size) {
            Point p1 = hull[hull.size() - 2];
            Point p2 = hull[hull.size() - 1];
            Point p3 = points[i];
            double cross = (p2.x - p1.x) * (p3.y - p1.y) - (p2.y - p1.y) * (p3.x - p1.x);
            if (cross <= 0) hull.pop_back();
            else break;
        }
        hull.push_back(points[i]);
    }
    
    // Remove duplicate points
    if (hull.size() > 1) hull.pop_back();
    
    std::vector<std::pair<double, double>> result;
    for (const auto& p : hull) {
        result.push_back({p.x, p.y});
    }
    return result;
}

double calculateConvexHullArea(const std::vector<std::pair<double, double>>& hull) {
    if (hull.size() < 3) return 0.0;
    
    double area = 0.0;
    for (size_t i = 0; i < hull.size(); i++) {
        size_t j = (i + 1) % hull.size();
        area += hull[i].first * hull[j].second;
        area -= hull[j].first * hull[i].second;
    }
    return std::abs(area) / 2.0;
}

double calculateConvexHullPerimeter(const std::vector<std::pair<double, double>>& hull) {
    if (hull.size() < 2) return 0.0;
    
    double perimeter = 0.0;
    for (size_t i = 0; i < hull.size(); i++) {
        size_t j = (i + 1) % hull.size();
        double dx = hull[i].first - hull[j].first;
        double dy = hull[i].second - hull[j].second;
        perimeter += std::sqrt(dx*dx + dy*dy);
    }
    return perimeter;
}

bool isAcidicAtom(const RDKit::Atom* atom) {
    int atomicNum = atom->getAtomicNum();
    int charge = atom->getFormalCharge();
    int nH = atom->getTotalNumHs();
    
    if (atomicNum == 8 && nH >= 1) return true; // -OH
    if (atomicNum == 16 && nH >= 1) return true; // -SH
    if (atomicNum == 7 && atom->getImplicitValence() > 2) return true; // -NH3+, etc
    if (atomicNum == 15 && atom->getExplicitValence() >= 4) return true; // phosphate
    
    return false;
}

bool isBasicAtom(const RDKit::Atom* atom) {
    int atomicNum = atom->getAtomicNum();
    int charge = atom->getFormalCharge();
    
    if (atomicNum == 7 && charge <= 0 && !isAcidicAtom(atom)) return true; // N not in acidic group
    if (atomicNum == 8 && charge < 0) return true; // O-
    
    return false;
}

} // anonymous namespace

// Base class implementation
ImageDescriptorBase::ImageDescriptorBase(const std::string& name, const std::string& description)
    : Descriptor(name, description) {}

std::unordered_map<std::string, double> ImageDescriptorBase::computeImageDescriptors(const Molecule& mol) {
    std::unordered_map<std::string, double> results;
    
    if (!mol.isValid() || !mol.getMolecule()) {
        globalLogger.debug("ImageDescriptorBase: Invalid molecule input.");
        return results;
    }
    
    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) {
        globalLogger.error("ImageDescriptorBase: RDKit molecule is null despite Molecule being valid?");
        return results;
    }
    
    // Make a copy for coords computation
    RDKit::ROMol molCopy(*rdkMol);
    
    try {
        // Compute 2D coordinates if needed
        if (!molCopy.getNumConformers()) {
            RDDepict::compute2DCoords(molCopy);
        }
        
        if (molCopy.getNumConformers() == 0) {
            globalLogger.error("ImageDescriptorBase: Failed to compute 2D coordinates");
            return results;
        }
        
        const auto& conf = molCopy.getConformer();
        
        // Compute min/max for 2D coordinates
        double minX = std::numeric_limits<double>::max();
        double minY = std::numeric_limits<double>::max();
        double maxX = std::numeric_limits<double>::lowest();
        double maxY = std::numeric_limits<double>::lowest();
        
        for (unsigned int i = 0; i < molCopy.getNumAtoms(); ++i) {
            const auto& pos = conf.getAtomPos(i);
            if (pos.x < minX) minX = pos.x;
            if (pos.y < minY) minY = pos.y;
            if (pos.x > maxX) maxX = pos.x;
            if (pos.y > maxY) maxY = pos.y;
        }
        
        // Calculate scaling and offset
        double width = 200;
        double height = 150;
        double padding = 20;
        
        double xRange = maxX - minX;
        double yRange = maxY - minY;
        double scale = std::min((width - 2*padding) / xRange, (height - 2*padding) / yRange);
        
        double xOffset = (width - xRange * scale) / 2;
        double yOffset = (height - yRange * scale) / 2;
        
        // Declare parameters and vectors for data collection
        double bondWidth = 3;
        double minSvgRadius = 4.0;
        double maxSvgRadius = 8.0;
        
        // Find min/max vdw for scaling
        double minVdw = std::numeric_limits<double>::max();
        double maxVdw = std::numeric_limits<double>::lowest();
        for (const auto atom : molCopy.atoms()) {
            double vdw = getVdwRadius(atom->getAtomicNum());
            if (vdw < minVdw) minVdw = vdw;
            if (vdw > maxVdw) maxVdw = vdw;
        }
        
        // Collect atom and bond geometry data
        std::vector<double> atomXs, atomYs, atomRadii, atomLuminances;
        struct Bond { double x1, y1, x2, y2; int atom1Idx, atom2Idx; };
        std::vector<Bond> bonds;
        
        // Store atom coordinates and properties
        for (unsigned int i = 0; i < molCopy.getNumAtoms(); ++i) {
            const auto atom = molCopy.getAtomWithIdx(i);
            const RDGeom::Point3D& pos = conf.getAtomPos(i);
            double x = (pos.x - minX) * scale + xOffset;
            double y = (pos.y - minY) * scale + yOffset;
            atomXs.push_back(x);
            atomYs.push_back(y);
            
            int atomicNum = atom->getAtomicNum();
            double vdw = getVdwRadius(atomicNum);
            double atomRadius = minSvgRadius;
            if (maxVdw > minVdw) {
                atomRadius = minSvgRadius + (maxSvgRadius - minSvgRadius) * (vdw - minVdw) / (maxVdw - minVdw);
            }
            atomRadii.push_back(atomRadius);
            
            auto [r, g, b] = getAtomColor(atomicNum);
            atomLuminances.push_back(getLuminance({r, g, b}));
        }
        
        // Store bond geometry
        for (const auto bond : molCopy.bonds()) {
            int beginIdx = bond->getBeginAtomIdx();
            int endIdx = bond->getEndAtomIdx();
            const RDGeom::Point3D& beginPos = conf.getAtomPos(beginIdx);
            const RDGeom::Point3D& endPos = conf.getAtomPos(endIdx);
            double x1 = (beginPos.x - minX) * scale + xOffset;
            double y1 = (beginPos.y - minY) * scale + yOffset;
            double x2 = (endPos.x - minX) * scale + xOffset;
            double y2 = (endPos.y - minY) * scale + yOffset;
            bonds.push_back({x1, y1, x2, y2, beginIdx, endIdx});
        }
        
        // Build adjacency list
        std::vector<std::vector<int>> adjList(atomXs.size());
        for (const auto& b : bonds) {
            adjList[b.atom1Idx].push_back(b.atom2Idx);
            adjList[b.atom2Idx].push_back(b.atom1Idx);
        }
        
        // Calculate bond lengths
        std::vector<double> bondLens;
        for (const auto& b : bonds) {
            double dx = b.x2 - b.x1, dy = b.y2 - b.y1;
            bondLens.push_back(std::sqrt(dx*dx + dy*dy));
        }
        
        // Calculate atom-atom distances
        std::vector<double> atomDistances;
        double sumDist = 0.0;
        int nPairs = 0;
        for (size_t i = 0; i < atomXs.size(); ++i) {
            for (size_t j = i+1; j < atomXs.size(); ++j) {
                double dx = atomXs[i] - atomXs[j], dy = atomYs[i] - atomYs[j];
                double dist = std::sqrt(dx*dx + dy*dy);
                atomDistances.push_back(dist);
                sumDist += dist;
                nPairs++;
            }
        }
        
        // Calculate atom areas
        std::vector<double> atomAreas;
        double totalAtomArea = 0.0;
        for (double r : atomRadii) {
            double area = M_PI * r * r;
            atomAreas.push_back(area);
            totalAtomArea += area;
        }
        
        // Calculate convex hull
        auto hull = calculateConvexHull(atomXs, atomYs);
        double convexHullArea = calculateConvexHullArea(hull);
        double convexHullPerimeter = calculateConvexHullPerimeter(hull);
        
        // Compute centroid
        double cxCentroid = 0.0, cyCentroid = 0.0;
        for (size_t i = 0; i < atomXs.size(); ++i) {
            cxCentroid += atomXs[i];
            cyCentroid += atomYs[i];
        }
        cxCentroid /= atomXs.size();
        cyCentroid /= atomYs.size();
        
        // Calculate distances to centroid
        std::vector<double> distsToCentroid;
        for (size_t i = 0; i < atomXs.size(); ++i) {
            double dx = atomXs[i] - cxCentroid, dy = atomYs[i] - cyCentroid;
            distsToCentroid.push_back(std::sqrt(dx*dx + dy*dy));
        }
        
        // Count atoms on hull
        int atomsOnHull = 0;
        for (const auto& pt : hull) {
            for (size_t i = 0; i < atomXs.size(); ++i) {
                if (std::abs(pt.first - atomXs[i]) < 1e-6 && std::abs(pt.second - atomYs[i]) < 1e-6) {
                    atomsOnHull++;
                    break;
                }
            }
        }
        
        // Calculate bond angles
        std::vector<double> bondAngles;
        for (size_t i = 0; i < atomXs.size(); ++i) {
            std::vector<std::pair<double, double>> neighbors;
            for (const auto& b : bonds) {
                if (std::abs(b.x1 - atomXs[i]) < 1e-6 && std::abs(b.y1 - atomYs[i]) < 1e-6)
                    neighbors.push_back({b.x2 - atomXs[i], b.y2 - atomYs[i]});
                else if (std::abs(b.x2 - atomXs[i]) < 1e-6 && std::abs(b.y2 - atomYs[i]) < 1e-6)
                    neighbors.push_back({b.x1 - atomXs[i], b.y1 - atomYs[i]});
            }
            for (size_t j = 0; j+1 < neighbors.size(); ++j) {
                for (size_t k = j+1; k < neighbors.size(); ++k) {
                    double dx1 = neighbors[j].first, dy1 = neighbors[j].second;
                    double dx2 = neighbors[k].first, dy2 = neighbors[k].second;
                    double dot = dx1*dx2 + dy1*dy2;
                    double norm1 = std::sqrt(dx1*dx1 + dy1*dy1);
                    double norm2 = std::sqrt(dx2*dx2 + dy2*dy2);
                    // Avoid division by zero or NaN errors
                    if (norm1 > 0 && norm2 > 0) {
                        double cosAngle = dot/(norm1*norm2);
                        // Clamp to avoid floating point errors
                        cosAngle = std::max(-1.0, std::min(1.0, cosAngle));
                        double angle = std::acos(cosAngle) * 180.0 / M_PI;
                        bondAngles.push_back(angle);
                    }
                }
            }
        }
        
        // Identify acidic and basic centers
        std::vector<int> acidicAtoms, basicAtoms;
        for (unsigned int i = 0; i < molCopy.getNumAtoms(); ++i) {
            const auto atom = molCopy.getAtomWithIdx(i);
            if (isAcidicAtom(atom)) acidicAtoms.push_back(i);
            if (isBasicAtom(atom)) basicAtoms.push_back(i);
        }
        
        // Calculate degrees of atoms
        std::vector<int> degrees;
        for (unsigned int i = 0; i < molCopy.getNumAtoms(); ++i) {
            degrees.push_back(molCopy.getAtomWithIdx(i)->getDegree());
        }
        
        // ===== DESCRIPTOR CALCULATIONS =====
        
        // 1. BoundingBoxArea
        double bboxArea = (maxX - minX) * (maxY - minY);
        results["BoundingBoxArea"] = bboxArea;
        
        // 2. MoleculeWidth
        double molWidth = maxX - minX;
        results["MoleculeWidth"] = molWidth;
        
        // 3. MoleculeHeight
        double molHeight = maxY - minY;
        results["MoleculeHeight"] = molHeight;
        
        // 4. AspectRatio
        double aspect = molWidth / (molHeight > 0 ? molHeight : 1.0);
        results["AspectRatio"] = aspect;
        
        // 5. AtomDensity
        double atomDensity = atomXs.size() / (bboxArea > 0 ? bboxArea : 1.0);
        results["AtomDensity"] = atomDensity;
        
        // 6. MeanAtomRadius
        double meanAtomRadius = !atomRadii.empty() ?
            std::accumulate(atomRadii.begin(), atomRadii.end(), 0.0) / atomRadii.size() : 0.0;
        results["MeanAtomRadius"] = meanAtomRadius;
        
        // 7. MaxAtomRadius
        double maxAtomRadius = !atomRadii.empty() ?
            *std::max_element(atomRadii.begin(), atomRadii.end()) : 0.0;
        results["MaxAtomRadius"] = maxAtomRadius;
        
        // 8. MinAtomRadius
        double minAtomRadius = !atomRadii.empty() ?
            *std::min_element(atomRadii.begin(), atomRadii.end()) : 0.0;
        results["MinAtomRadius"] = minAtomRadius;
        
        // 9. MeanBondLength
        double meanBondLen = !bondLens.empty() ?
            std::accumulate(bondLens.begin(), bondLens.end(), 0.0) / bondLens.size() : 0.0;
        results["MeanBondLength"] = meanBondLen;
        
        // 10. MaxBondLength
        double maxBondLen = !bondLens.empty() ?
            *std::max_element(bondLens.begin(), bondLens.end()) : 0.0;
        results["MaxBondLength"] = maxBondLen;
        
        // 11. MinBondLength
        double minBondLen = !bondLens.empty() ?
            *std::min_element(bondLens.begin(), bondLens.end()) : 0.0;
        results["MinBondLength"] = minBondLen;
        
        // 12. BondLengthStd
        double bondLenStd = 0.0;
        if (!bondLens.empty()) {
            double sum = 0.0;
            for (double len : bondLens) {
                sum += (len - meanBondLen) * (len - meanBondLen);
            }
            bondLenStd = std::sqrt(sum / bondLens.size());
        }
        results["BondLengthStd"] = bondLenStd;
        
        // 13. MeanAtomAtomDist
        double meanAtomDist = nPairs > 0 ? sumDist / nPairs : 0.0;
        results["MeanAtomAtomDist"] = meanAtomDist;
        
        // 14. MeanAtomLuminance
        double meanLuminance = !atomLuminances.empty() ?
            std::accumulate(atomLuminances.begin(), atomLuminances.end(), 0.0) / atomLuminances.size() : 0.0;
        results["MeanAtomLuminance"] = meanLuminance;
        
        // 15. MaxAtomAtomDist
        double maxAtomDist = !atomDistances.empty() ?
            *std::max_element(atomDistances.begin(), atomDistances.end()) : 0.0;
        results["MaxAtomAtomDist"] = maxAtomDist;
        
        // 16. MinAtomAtomDist
        double minAtomDist = !atomDistances.empty() ?
            *std::min_element(atomDistances.begin(), atomDistances.end()) : 0.0;
        results["MinAtomAtomDist"] = minAtomDist;
        
        // 17. MeanAtomArea
        double meanAtomArea = !atomAreas.empty() ?
            std::accumulate(atomAreas.begin(), atomAreas.end(), 0.0) / atomAreas.size() : 0.0;
        results["MeanAtomArea"] = meanAtomArea;
        
        // 18. MedianAtomRadius
        double medianAtomRadius = 0.0;
        if (!atomRadii.empty()) {
            std::vector<double> sortedRadii = atomRadii;
            std::sort(sortedRadii.begin(), sortedRadii.end());
            medianAtomRadius = sortedRadii[sortedRadii.size()/2];
        }
        results["MedianAtomRadius"] = medianAtomRadius;
        
        // 19. MedianBondLength
        double medianBondLen = 0.0;
        if (!bondLens.empty()) {
            std::vector<double> sortedBondLens = bondLens;
            std::sort(sortedBondLens.begin(), sortedBondLens.end());
            medianBondLen = sortedBondLens[sortedBondLens.size()/2];
        }
        results["MedianBondLength"] = medianBondLen;
        
        // 20. AtomRadiusRange
        double atomRadiusRange = maxAtomRadius - minAtomRadius;
        results["AtomRadiusRange"] = atomRadiusRange;
        
        // 21. AtomAreaFraction
        double atomAreaFraction = bboxArea > 0 ? totalAtomArea / bboxArea : 0.0;
        results["AtomAreaFraction"] = atomAreaFraction;
        
        // 22. BondCoverageFraction
        double bboxPerimeter = 2 * (molWidth + molHeight);
        double totalBondLen = std::accumulate(bondLens.begin(), bondLens.end(), 0.0);
        double bondCoverageFraction = bboxPerimeter > 0 ? totalBondLen / bboxPerimeter : 0.0;
        results["BondCoverageFraction"] = bondCoverageFraction;
        
        // 23. AtomPackingDensity
        double atomPackingDensity = convexHullArea > 0 ? totalAtomArea / convexHullArea : 0.0;
        results["AtomPackingDensity"] = atomPackingDensity;
        
        // 24. AtomMassX
        double atomMassX = cxCentroid;
        results["AtomMassX"] = atomMassX;
        
        // 25. AtomMassY
        double atomMassY = cyCentroid;
        results["AtomMassY"] = atomMassY;
        
        // 26. AtomMassDist
        double imgCenterX = width/2, imgCenterY = height/2;
        double atomMassDist = std::sqrt((atomMassX-imgCenterX)*(atomMassX-imgCenterX) + 
                                         (atomMassY-imgCenterY)*(atomMassY-imgCenterY));
        results["AtomMassDist"] = atomMassDist;
        
        // 27. AtomRadiusStd
        double atomRadiusStd = 0.0;
        if (!atomRadii.empty()) {
            double sum = 0.0;
            for (double r : atomRadii) {
                sum += (r - meanAtomRadius) * (r - meanAtomRadius);
            }
            atomRadiusStd = std::sqrt(sum / atomRadii.size());
        }
        results["AtomRadiusStd"] = atomRadiusStd;
        
        // 28. MeanBondAngle
        double meanBondAngle = !bondAngles.empty() ?
            std::accumulate(bondAngles.begin(), bondAngles.end(), 0.0) / bondAngles.size() : 0.0;
        results["MeanBondAngle"] = meanBondAngle;
        
        // 29. BondAngleStd
        double bondAngleStd = 0.0;
        if (!bondAngles.empty()) {
            double sum = 0.0;
            for (double a : bondAngles) {
                sum += (a - meanBondAngle) * (a - meanBondAngle);
            }
            bondAngleStd = std::sqrt(sum / bondAngles.size());
        }
        results["BondAngleStd"] = bondAngleStd;
        
        // 30. AtomXStd
        double atomXStd = 0.0;
        if (!atomXs.empty()) {
            double sum = 0.0;
            for (double x : atomXs) {                sum += (x - atomMassX) * (x - atomMassX);
            }
            atomXStd = std::sqrt(sum / atomXs.size());
        }
        results["AtomXStd"] = atomXStd;
        
        // 31. AtomYStd
        double atomYStd = 0.0;
        if (!atomYs.empty()) {
            double sum = 0.0;
            for (double y : atomYs) {
                sum += (y - atomMassY) * (y - atomMassY);
            }
            atomYStd = std::sqrt(sum / atomYs.size());
        }
        results["AtomYStd"] = atomYStd;
        
        // 32. AtomLuminanceStd
        double atomLuminanceStd = 0.0;
        if (!atomLuminances.empty()) {
            double sum = 0.0;
            for (double l : atomLuminances) {
                sum += (l - meanLuminance) * (l - meanLuminance);
            }
            atomLuminanceStd = std::sqrt(sum / atomLuminances.size());
        }
        results["AtomLuminanceStd"] = atomLuminanceStd;
        
        // 33. BondLenMAD
        double bondLenMAD = 0.0;
        if (!bondLens.empty()) {
            std::vector<double> absdevs;
            for (double l : bondLens) {
                absdevs.push_back(std::abs(l - medianBondLen));
            }
            std::sort(absdevs.begin(), absdevs.end());
            bondLenMAD = absdevs[absdevs.size()/2];
        }
        results["BondLenMAD"] = bondLenMAD;
        
        // 34. AtomRadiusMAD
        double atomRadiusMAD = 0.0;
        if (!atomRadii.empty()) {
            std::vector<double> absdevs;
            for (double r : atomRadii) {
                absdevs.push_back(std::abs(r - medianAtomRadius));
            }
            std::sort(absdevs.begin(), absdevs.end());
            atomRadiusMAD = absdevs[absdevs.size()/2];
        }
        results["AtomRadiusMAD"] = atomRadiusMAD;
        
        // 35. AtomAreaStd
        double atomAreaStd = 0.0;
        if (!atomAreas.empty()) {
            double sum = 0.0;
            for (double a : atomAreas) {
                sum += (a - meanAtomArea) * (a - meanAtomArea);
            }
            atomAreaStd = std::sqrt(sum / atomAreas.size());
        }
        results["AtomAreaStd"] = atomAreaStd;
        
        // 36. AtomRadiusCV
        double atomRadiusCV = meanAtomRadius > 0 ? atomRadiusStd / meanAtomRadius : 0.0;
        results["AtomRadiusCV"] = atomRadiusCV;
        
        // 37. AtomAreaCV
        double atomAreaCV = meanAtomArea > 0 ? atomAreaStd / meanAtomArea : 0.0;
        results["AtomAreaCV"] = atomAreaCV;
        
        // 38. AtomLuminanceCV
        double atomLuminanceCV = meanLuminance > 0 ? atomLuminanceStd / meanLuminance : 0.0;
        results["AtomLuminanceCV"] = atomLuminanceCV;
        
        // 39. AtomAreaRange
        double minArea = M_PI * minAtomRadius * minAtomRadius;
        double maxArea = M_PI * maxAtomRadius * maxAtomRadius;
        double atomAreaRange = maxArea - minArea;
        results["AtomAreaRange"] = atomAreaRange;
        
        // 40. AtomAreaMedian
        double atomAreaMedian = 0.0;
        if (!atomAreas.empty()) {
            std::vector<double> sortedAreas = atomAreas;
            std::sort(sortedAreas.begin(), sortedAreas.end());
            atomAreaMedian = sortedAreas[sortedAreas.size()/2];
        }
        results["AtomAreaMedian"] = atomAreaMedian;
        
        // 41. AtomAreaMAD
        double atomAreaMAD = 0.0;
        if (!atomAreas.empty()) {
            std::vector<double> absdevs;
            for (double a : atomAreas) {
                absdevs.push_back(std::abs(a - atomAreaMedian));
            }
            std::sort(absdevs.begin(), absdevs.end());
            atomAreaMAD = absdevs[absdevs.size()/2];
        }
        results["AtomAreaMAD"] = atomAreaMAD;
        
        // 42. AtomAreaHullFrac
        double hullAtomArea = 0.0;
        for (const auto& pt : hull) {
            for (size_t i = 0; i < atomXs.size(); ++i) {
                if (std::abs(pt.first - atomXs[i]) < 1e-6 && std::abs(pt.second - atomYs[i]) < 1e-6) {
                    hullAtomArea += atomAreas[i];
                    break;
                }
            }
        }
        double atomAreaHullFrac = totalAtomArea > 0 ? hullAtomArea / totalAtomArea : 0.0;
        results["AtomAreaHullFrac"] = atomAreaHullFrac;
        
        // 43. AtomAreaCenterFrac
        double cx = (minX + maxX)/2, cy = (minY + maxY)/2;
        double centerArea = 0.0;
        for (size_t i = 0; i < atomXs.size(); ++i) {
            if (std::abs(atomXs[i] - cxCentroid) < molWidth/4 && std::abs(atomYs[i] - cyCentroid) < molHeight/4) {
                centerArea += atomAreas[i];
            }
        }
        double atomAreaCenterFrac = totalAtomArea > 0 ? centerArea / totalAtomArea : 0.0;
        results["AtomAreaCenterFrac"] = atomAreaCenterFrac;
        
        // 44. StdAllAtomDist
        double stdAllDist = 0.0;
        if (!atomDistances.empty()) {
            double sum = 0.0;
            for (double d : atomDistances) {
                sum += (d - meanAtomDist) * (d - meanAtomDist);
            }
            stdAllDist = std::sqrt(sum / atomDistances.size());
        }
        results["StdAllAtomDist"] = stdAllDist;
        
        // 45. MinBondAngle
        double minBondAngle = !bondAngles.empty() ? 
            *std::min_element(bondAngles.begin(), bondAngles.end()) : 0.0;
        results["MinBondAngle"] = minBondAngle;
        
        // 46. MaxBondAngle
        double maxBondAngle = !bondAngles.empty() ? 
            *std::max_element(bondAngles.begin(), bondAngles.end()) : 0.0;
        results["MaxBondAngle"] = maxBondAngle;
        
        // 47. MedianBondAngle
        double medianBondAngle = 0.0;
        if (!bondAngles.empty()) {
            std::vector<double> sortedAngles = bondAngles;
            std::sort(sortedAngles.begin(), sortedAngles.end());
            medianBondAngle = sortedAngles[sortedAngles.size()/2];
        }
        results["MedianBondAngle"] = medianBondAngle;
        
        // 48. MeanBondColorDiff
        double meanBondColorDiff = 0.0;
        int colorDiffCount = 0;
        for (const auto& b : bonds) {
            int idx1 = b.atom1Idx;
            int idx2 = b.atom2Idx;
            auto rgb1 = getAtomColor(molCopy.getAtomWithIdx(idx1)->getAtomicNum());
            auto rgb2 = getAtomColor(molCopy.getAtomWithIdx(idx2)->getAtomicNum());
            double r1 = std::get<0>(rgb1), g1 = std::get<1>(rgb1), b1 = std::get<2>(rgb1);
            double r2 = std::get<0>(rgb2), g2 = std::get<1>(rgb2), b2 = std::get<2>(rgb2);
            double diff = std::sqrt((r1-r2)*(r1-r2) + (g1-g2)*(g1-g2) + (b1-b2)*(b1-b2));
            meanBondColorDiff += diff;
            colorDiffCount++;
        }
        meanBondColorDiff = colorDiffCount > 0 ? meanBondColorDiff / colorDiffCount : 0.0;
        results["MeanBondColorDiff"] = meanBondColorDiff;
        
        // 49. MeanBondLuminanceDiff
        double meanBondLuminanceDiff = 0.0;
        int luminanceDiffCount = 0;
        for (const auto& b : bonds) {
            int idx1 = b.atom1Idx;
            int idx2 = b.atom2Idx;
            double diff = std::abs(atomLuminances[idx1] - atomLuminances[idx2]);
            meanBondLuminanceDiff += diff;
            luminanceDiffCount++;
        }
        meanBondLuminanceDiff = luminanceDiffCount > 0 ? meanBondLuminanceDiff / luminanceDiffCount : 0.0;
        results["MeanBondLuminanceDiff"] = meanBondLuminanceDiff;
        
        // 50. MeanBondLenDiffColor
        double meanBondLenDiffColor = 0.0;
        int bondLenDiffColorCount = 0;
        for (const auto& b : bonds) {
            int idx1 = b.atom1Idx;
            int idx2 = b.atom2Idx;
            if (atomLuminances[idx1] != atomLuminances[idx2]) {
                double dx = b.x2 - b.x1, dy = b.y2 - b.y1;
                meanBondLenDiffColor += std::sqrt(dx*dx + dy*dy);
                bondLenDiffColorCount++;
            }
        }
        meanBondLenDiffColor = bondLenDiffColorCount > 0 ? meanBondLenDiffColor / bondLenDiffColorCount : 0.0;
        results["MeanBondLenDiffColor"] = meanBondLenDiffColor;
        
        // 51. MeanBondAngleDiffColor
        double meanBondAngleDiffColor = 0.0;
        int bondAngleDiffColorCount = 0;
        
        for (size_t i = 0; i < atomXs.size(); ++i) {
            std::vector<std::pair<double, double>> neighbors;
            std::vector<int> neighborIdxs;
            
            for (const auto& b : bonds) {
                if (b.atom1Idx == i) {
                    neighbors.push_back({atomXs[b.atom2Idx] - atomXs[i], atomYs[b.atom2Idx] - atomYs[i]});
                    neighborIdxs.push_back(b.atom2Idx);
                } else if (b.atom2Idx == i) {
                    neighbors.push_back({atomXs[b.atom1Idx] - atomXs[i], atomYs[b.atom1Idx] - atomYs[i]});
                    neighborIdxs.push_back(b.atom1Idx);
                }
            }
            
            for (size_t j = 0; j + 1 < neighbors.size(); ++j) {
                for (size_t k = j + 1; k < neighbors.size(); ++k) {
                    if (atomLuminances[neighborIdxs[j]] != atomLuminances[neighborIdxs[k]]) {
                        double dx1 = neighbors[j].first, dy1 = neighbors[j].second;
                        double dx2 = neighbors[k].first, dy2 = neighbors[k].second;
                        double dot = dx1*dx2 + dy1*dy2;
                        double norm1 = std::sqrt(dx1*dx1 + dy1*dy1);
                        double norm2 = std::sqrt(dx2*dx2 + dy2*dy2);
                        
                        if (norm1 > 0 && norm2 > 0) {
                            double cosAngle = dot/(norm1*norm2);
                            cosAngle = std::max(-1.0, std::min(1.0, cosAngle));
                            double angle = std::acos(cosAngle) * 180.0 / M_PI;
                            meanBondAngleDiffColor += angle;
                            bondAngleDiffColorCount++;
                        }
                    }
                }
            }
        }
        
        meanBondAngleDiffColor = bondAngleDiffColorCount > 0 ? meanBondAngleDiffColor / bondAngleDiffColorCount : 0.0;
        results["MeanBondAngleDiffColor"] = meanBondAngleDiffColor;
        
        // 52. MeanDistSameColor
        double meanDistSameColor = 0.0;
        int sameColorCount = 0;
        for (size_t i = 0; i < atomXs.size(); ++i) {
            for (size_t j = i+1; j < atomXs.size(); ++j) {
                if (std::abs(atomLuminances[i] - atomLuminances[j]) < 1e-6) {
                    double dx = atomXs[i] - atomXs[j], dy = atomYs[i] - atomYs[j];
                    meanDistSameColor += std::sqrt(dx*dx + dy*dy);
                    sameColorCount++;
                }
            }
        }
        meanDistSameColor = sameColorCount > 0 ? meanDistSameColor / sameColorCount : 0.0;
        results["MeanDistSameColor"] = meanDistSameColor;
        
        // 53. MeanDistDiffColor
        double meanDistDiffColor = 0.0;
        int diffColorCount = 0;
        for (size_t i = 0; i < atomXs.size(); ++i) {
            for (size_t j = i+1; j < atomXs.size(); ++j) {
                if (std::abs(atomLuminances[i] - atomLuminances[j]) > 1e-6) {
                    double dx = atomXs[i] - atomXs[j], dy = atomYs[i] - atomYs[j];
                    meanDistDiffColor += std::sqrt(dx*dx + dy*dy);
                    diffColorCount++;
                }
            }
        }
        meanDistDiffColor = diffColorCount > 0 ? meanDistDiffColor / diffColorCount : 0.0;
        results["MeanDistDiffColor"] = meanDistDiffColor;
        
        // 54. MeanBondOrder
        double sumBondOrder = 0.0;
        int nDouble = 0, nTriple = 0, nAromatic = 0;
        for (const auto& bond : molCopy.bonds()) {
            double order = 1.0;
            if (bond->getBondType() == RDKit::Bond::BondType::DOUBLE) { order = 2.0; nDouble++; }
            else if (bond->getBondType() == RDKit::Bond::BondType::TRIPLE) { order = 3.0; nTriple++; }
            else if (bond->getIsAromatic()) { order = 1.5; nAromatic++; }
            sumBondOrder += order;
        }
        double meanBondOrder = bonds.size() > 0 ? sumBondOrder / bonds.size() : 0.0;
        results["MeanBondOrder"] = meanBondOrder;
        
        // 55. FracDoubleBonds
        double fracDouble = bonds.size() > 0 ? (double)nDouble / bonds.size() : 0.0;
        results["FracDoubleBonds"] = fracDouble;
        
        // 56. FracTripleBonds
        double fracTriple = bonds.size() > 0 ? (double)nTriple / bonds.size() : 0.0;
        results["FracTripleBonds"] = fracTriple;
        
        // 57. FracAromaticBonds
        double fracAromatic = bonds.size() > 0 ? (double)nAromatic / bonds.size() : 0.0;
        results["FracAromaticBonds"] = fracAromatic;
        
        // 58. MeanAtomDegree
        double meanDegree = !degrees.empty() ? 
            std::accumulate(degrees.begin(), degrees.end(), 0.0) / degrees.size() : 0.0;
        results["MeanAtomDegree"] = meanDegree;
        
        // 59. MaxAtomDegree
        int maxDegree = !degrees.empty() ? *std::max_element(degrees.begin(), degrees.end()) : 0;
        results["MaxAtomDegree"] = maxDegree;
        
        // 60. MinAtomDegree
        int minDegree = !degrees.empty() ? *std::min_element(degrees.begin(), degrees.end()) : 0;
        results["MinAtomDegree"] = minDegree;
        
        // 61. MeanDistToCentroid
        double meanDistToCentroid = !distsToCentroid.empty() ?
            std::accumulate(distsToCentroid.begin(), distsToCentroid.end(), 0.0) / distsToCentroid.size() : 0.0;
        results["MeanDistToCentroid"] = meanDistToCentroid;
        
        // 62. StdDistToCentroid
        double stdDistToCentroid = 0.0;
        if (!distsToCentroid.empty()) {
            double sum = 0.0;
            for (double d : distsToCentroid) {
                sum += (d - meanDistToCentroid) * (d - meanDistToCentroid);
            }
            stdDistToCentroid = std::sqrt(sum / distsToCentroid.size());
        }
        results["StdDistToCentroid"] = stdDistToCentroid;
        
        // Calculate mean bond lengths by type
        std::vector<double> singleLens, doubleLens, tripleLens, aromLens;
        for (size_t i = 0; i < molCopy.getNumBonds(); ++i) {
            const auto bond = molCopy.getBondWithIdx(i);
            int idx1 = bond->getBeginAtomIdx(), idx2 = bond->getEndAtomIdx();
            double dx = atomXs[idx1] - atomXs[idx2], dy = atomYs[idx1] - atomYs[idx2];
            double len = std::sqrt(dx*dx + dy*dy);
            
            if (bond->getBondType() == RDKit::Bond::BondType::SINGLE) singleLens.push_back(len);
            else if (bond->getBondType() == RDKit::Bond::BondType::DOUBLE) doubleLens.push_back(len);
            else if (bond->getBondType() == RDKit::Bond::BondType::TRIPLE) tripleLens.push_back(len);
            else if (bond->getIsAromatic()) aromLens.push_back(len);
        }
        
        // 63. MeanLenSingle
        double meanLenSingle = !singleLens.empty() ?
            std::accumulate(singleLens.begin(), singleLens.end(), 0.0) / singleLens.size() : 0.0;
        results["MeanLenSingle"] = meanLenSingle;
        
        // 64. MeanLenDouble
        double meanLenDouble = !doubleLens.empty() ?
            std::accumulate(doubleLens.begin(), doubleLens.end(), 0.0) / doubleLens.size() : 0.0;
        results["MeanLenDouble"] = meanLenDouble;
        
        // 65. MeanLenTriple
        double meanLenTriple = !tripleLens.empty() ?
            std::accumulate(tripleLens.begin(), tripleLens.end(), 0.0) / tripleLens.size() : 0.0;
        results["MeanLenTriple"] = meanLenTriple;
        
        // 66. MeanLenAromatic
        double meanLenArom = !aromLens.empty() ?
            std::accumulate(aromLens.begin(), aromLens.end(), 0.0) / aromLens.size() : 0.0;
        results["MeanLenAromatic"] = meanLenArom;
        
        // Calculate fractions of atoms by degree
        int nDeg1 = 0, nDeg2 = 0, nDeg3 = 0, nDeg4 = 0;
        for (int d : degrees) {
            if (d == 1) nDeg1++;
            else if (d == 2) nDeg2++;
            else if (d == 3) nDeg3++;
            else if (d >= 4) nDeg4++; // Group all degrees 4 and higher
        }
        
        // 67. FracDeg1
        double fracDeg1 = atomXs.size() > 0 ? (double)nDeg1 / atomXs.size() : 0.0;
        results["FracDeg1"] = fracDeg1;
        
        // 68. FracDeg2
        double fracDeg2 = atomXs.size() > 0 ? (double)nDeg2 / atomXs.size() : 0.0;
        results["FracDeg2"] = fracDeg2;
        
        // 69. FracDeg3
        double fracDeg3 = atomXs.size() > 0 ? (double)nDeg3 / atomXs.size() : 0.0;
        results["FracDeg3"] = fracDeg3;
        
        // 70. FracDeg4
        double fracDeg4 = atomXs.size() > 0 ? (double)nDeg4 / atomXs.size() : 0.0;
        results["FracDeg4"] = fracDeg4;
        
        // 71. MeanBondRadiiDiff
        double sumBondRadiiDiff = 0.0;
        for (const auto& bond : molCopy.bonds()) {
            int idx1 = bond->getBeginAtomIdx(), idx2 = bond->getEndAtomIdx();
            sumBondRadiiDiff += std::abs(atomRadii[idx1] - atomRadii[idx2]);
        }
        double meanBondRadiiDiff = bonds.size() > 0 ? sumBondRadiiDiff / bonds.size() : 0.0;
        results["MeanBondRadiiDiff"] = meanBondRadiiDiff;
        
        // 72. MeanBondLuminanceDiff2 (already calculated as MeanBondLuminanceDiff)
        results["MeanBondLuminanceDiff2"] = meanBondLuminanceDiff;
        
        // 73. MeanAngleHighDegree
        std::vector<double> anglesAtHighDegree;
        for (size_t i = 0; i < atomXs.size(); ++i) {
            if (degrees[i] < 3) continue;
            
            std::vector<std::pair<double, double>> neighbors;
            for (const auto& bond : molCopy.bonds()) {
                int idx1 = bond->getBeginAtomIdx(), idx2 = bond->getEndAtomIdx();
                if (idx1 == i) neighbors.push_back({atomXs[idx2] - atomXs[i], atomYs[idx2] - atomYs[i]});
                else if (idx2 == i) neighbors.push_back({atomXs[idx1] - atomXs[i], atomYs[idx1] - atomYs[i]});
            }
            
            for (size_t j = 0; j + 1 < neighbors.size(); ++j) {
                for (size_t k = j + 1; k < neighbors.size(); ++k) {
                    double dx1 = neighbors[j].first, dy1 = neighbors[j].second;
                    double dx2 = neighbors[k].first, dy2 = neighbors[k].second;
                    double dot = dx1*dx2 + dy1*dy2;
                    double norm1 = std::sqrt(dx1*dx1 + dy1*dy1);
                    double norm2 = std::sqrt(dx2*dx2 + dy2*dy2);
                    
                    if (norm1 > 0 && norm2 > 0) {
                        double cosAngle = dot/(norm1*norm2);
                        cosAngle = std::max(-1.0, std::min(1.0, cosAngle));
                        double angle = std::acos(cosAngle) * 180.0 / M_PI;
                        anglesAtHighDegree.push_back(angle);
                    }
                }
            }
        }
        
        double meanAngleHighDegree = !anglesAtHighDegree.empty() ?
            std::accumulate(anglesAtHighDegree.begin(), anglesAtHighDegree.end(), 0.0) / anglesAtHighDegree.size() : 0.0;
        results["MeanAngleHighDegree"] = meanAngleHighDegree;
        
        // 74. BoundaryAtomRatio (already calculated as fraction of atoms on hull)
        double boundaryAtomRatio = atomXs.size() > 0 ? (double)atomsOnHull / atomXs.size() : 0.0;
        results["BoundaryAtomRatio"] = boundaryAtomRatio;
        
        // 75. MeanAtomEccentricity
        double maxDistToCentroid = !distsToCentroid.empty() ?
            *std::max_element(distsToCentroid.begin(), distsToCentroid.end()) : 0.0;
        
        std::vector<double> atomEccentricities;
        for (double dist : distsToCentroid) {
            atomEccentricities.push_back(maxDistToCentroid > 0 ? dist / maxDistToCentroid : 0.0);
        }
        
        double meanAtomEccentricity = !atomEccentricities.empty() ?
            std::accumulate(atomEccentricities.begin(), atomEccentricities.end(), 0.0) / atomEccentricities.size() : 0.0;
        results["MeanAtomEccentricity"] = meanAtomEccentricity;
        
        // 76. MolecularDiameter (maxAtomDist)
        double molecularDiameter = maxAtomDist;
        results["MolecularDiameter"] = molecularDiameter;
        
        // 77. MolecularRadius (maxDistToCentroid)
        double molecularRadius = maxDistToCentroid;
        results["MolecularRadius"] = molecularRadius;
        
        // 78. RadiusOfGyration
        double radiusOfGyration = 0.0;
        if (!distsToCentroid.empty()) {
            double sum = 0.0;
            for (double d : distsToCentroid) {
                sum += d * d;
            }
            radiusOfGyration = std::sqrt(sum / distsToCentroid.size());
        }
        results["RadiusOfGyration"] = radiusOfGyration;
        
        // 79. MolecularSphericity
        double inscribedRadius = std::numeric_limits<double>::max();
        for (const auto& pt : hull) {
            double dx = pt.first - cxCentroid, dy = pt.second - cyCentroid;
            double distToCenter = std::sqrt(dx*dx + dy*dy);
            if (distToCenter < inscribedRadius) inscribedRadius = distToCenter;
        }
        double circumscribedRadius = molecularRadius;
        double molecularSphericity = circumscribedRadius > 0 ? inscribedRadius / circumscribedRadius : 0.0;
        results["MolecularSphericity"] = molecularSphericity;
        
        // 80. BondLenToAtomRadiusRatio
        double bondLenToAtomRadiusRatio = meanAtomRadius > 0 ? meanBondLen / meanAtomRadius : 0.0;
        results["BondLenToAtomRadiusRatio"] = bondLenToAtomRadiusRatio;
        
        // 81. EdgeDensity
        double edgeDensity = convexHullArea > 0 ? bonds.size() / convexHullArea : 0.0;
        results["EdgeDensity"] = edgeDensity;
        
        // 82. PlanarityMeasure (trivially zero in 2D)
        double planarityMeasure = 0.0;
        results["PlanarityMeasure"] = planarityMeasure;
        
        // 83. MeanNearestNeighborDist
        std::vector<double> nearestNeighborDists;
        for (size_t i = 0; i < atomXs.size(); ++i) {
            double minDist = std::numeric_limits<double>::max();
            for (size_t j = 0; j < atomXs.size(); ++j) {
                if (i == j) continue;
                double dx = atomXs[i] - atomXs[j], dy = atomYs[i] - atomYs[j];
                double dist = std::sqrt(dx*dx + dy*dy);
                if (dist < minDist) minDist = dist;
            }
            if (minDist < std::numeric_limits<double>::max()) {
                nearestNeighborDists.push_back(minDist);
            }
        }
        double meanNearestNeighborDist = !nearestNeighborDists.empty() ?
            std::accumulate(nearestNeighborDists.begin(), nearestNeighborDists.end(), 0.0) / nearestNeighborDists.size() : 0.0;
        results["MeanNearestNeighborDist"] = meanNearestNeighborDist;
        
        // 84. MeanAtomsInRadius
        double neighborhoodRadius = meanBondLen * 1.5; // 1.5 times mean bond length
        std::vector<int> atomsInRadius;
        for (size_t i = 0; i < atomXs.size(); ++i) {
            int count = 0;
            for (size_t j = 0; j < atomXs.size(); ++j) {
                if (i == j) continue;
                double dx = atomXs[i] - atomXs[j], dy = atomYs[i] - atomYs[j];
                double dist = std::sqrt(dx*dx + dy*dy);
                if (dist <= neighborhoodRadius) count++;
            }
            atomsInRadius.push_back(count);
        }
        double meanAtomsInRadius = !atomsInRadius.empty() ?
            std::accumulate(atomsInRadius.begin(), atomsInRadius.end(), 0.0) / atomsInRadius.size() : 0.0;
        results["MeanAtomsInRadius"] = meanAtomsInRadius;
        
        // 85. VarNearestNeighborDist
        double varNearestNeighborDist = 0.0;
        if (!nearestNeighborDists.empty()) {
            double sum = 0.0;
            for (double d : nearestNeighborDists) {
                sum += (d - meanNearestNeighborDist) * (d - meanNearestNeighborDist);
            }
            varNearestNeighborDist = sum / nearestNeighborDists.size();
        }
        results["VarNearestNeighborDist"] = varNearestNeighborDist;
        
        // 86. MeanBondToBondAngle
        std::vector<double> bondToBondAngles;
        for (size_t i = 0; i < atomXs.size(); ++i) {
            if (adjList[i].size() < 2) continue;
            
            for (size_t j = 0; j < adjList[i].size(); ++j) {
                for (size_t k = j+1; k < adjList[i].size(); ++k) {
                    int n1 = adjList[i][j], n2 = adjList[i][k];
                    
                    double dx1 = atomXs[n1] - atomXs[i], dy1 = atomYs[n1] - atomYs[i];
                    double dx2 = atomXs[n2] - atomXs[i], dy2 = atomYs[n2] - atomYs[i];
                    
                    double len1 = std::sqrt(dx1*dx1 + dy1*dy1);
                    double len2 = std::sqrt(dx2*dx2 + dy2*dy2);
                    
                    if (len1 > 0 && len2 > 0) {
                        double cosAngle = (dx1*dx2 + dy1*dy2) / (len1 * len2);
                        cosAngle = std::max(-1.0, std::min(1.0, cosAngle));
                        double angle = std::acos(cosAngle) * 180.0 / M_PI;
                        bondToBondAngles.push_back(angle);
                    }
                }
            }
        }
        double meanBondToBondAngle = !bondToBondAngles.empty() ?
            std::accumulate(bondToBondAngles.begin(), bondToBondAngles.end(), 0.0) / bondToBondAngles.size() : 0.0;
        results["MeanBondToBondAngle"] = meanBondToBondAngle;
        
        // 87. MomentOfInertia
        double momentOfInertia = 0.0;
        for (size_t i = 0; i < atomXs.size(); ++i) {
            double dx = atomXs[i] - cxCentroid, dy = atomYs[i] - cyCentroid;
            momentOfInertia += dx*dx + dy*dy;
        }
        results["MomentOfInertia"] = momentOfInertia;
        
        // 88. MeanAtomCentrality
        std::vector<double> atomCentrality;
        for (size_t i = 0; i < atomXs.size(); ++i) {
            double totalDist = 0.0;
            for (size_t j = 0; j < atomXs.size(); ++j) {
                if (i != j) {
                    double dx = atomXs[i] - atomXs[j], dy = atomYs[i] - atomYs[j];
                    totalDist += std::sqrt(dx*dx + dy*dy);
                }
            }
            atomCentrality.push_back(atomXs.size() > 1 ? totalDist / (atomXs.size() - 1) : 0.0);
        }
        double meanAtomCentrality = !atomCentrality.empty() ?
            std::accumulate(atomCentrality.begin(), atomCentrality.end(), 0.0) / atomCentrality.size() : 0.0;
        results["MeanAtomCentrality"] = meanAtomCentrality;
        
        // 89. TerminalAtomCount (= nDeg1)
        int terminalAtomCount = nDeg1;
        results["TerminalAtomCount"] = terminalAtomCount;
        
        // 90. JunctionAtomCount
        int junctionAtomCount = 0;
        for (int d : degrees) {
            if (d >= 3) junctionAtomCount++;
        }
        results["JunctionAtomCount"] = junctionAtomCount;
        
        // 91. TerminalToJunctionRatio
        double terminalToJunctionRatio = junctionAtomCount > 0 ? (double)terminalAtomCount / junctionAtomCount : 0.0;
        results["TerminalToJunctionRatio"] = terminalToJunctionRatio;
        
        // 92. WienerIndex
        double wienerIndex = 0.0;
        std::vector<std::vector<double>> dist(atomXs.size(), std::vector<double>(atomXs.size(), std::numeric_limits<double>::max()));
        
        // Initialize with direct connections
        for (size_t i = 0; i < atomXs.size(); ++i) {
            dist[i][i] = 0.0;
            for (int j : adjList[i]) {
                dist[i][j] = 1.0; // Unit distance for bonds
            }
        }
        
        // Floyd-Warshall algorithm
        for (size_t k = 0; k < atomXs.size(); ++k) {
            for (size_t i = 0; i < atomXs.size(); ++i) {
                for (size_t j = 0; j < atomXs.size(); ++j) {
                    if (dist[i][k] != std::numeric_limits<double>::max() && 
                        dist[k][j] != std::numeric_limits<double>::max()) {
                        dist[i][j] = std::min(dist[i][j], dist[i][k] + dist[k][j]);
                    }
                }
            }
        }
        
        // Sum up the shortest paths
        for (size_t i = 0; i < atomXs.size(); ++i) {
            for (size_t j = i+1; j < atomXs.size(); ++j) {
                if (dist[i][j] != std::numeric_limits<double>::max()) {
                    wienerIndex += dist[i][j];
                }
            }
        }
        results["WienerIndex"] = wienerIndex;
        
        // 93. RingCount
        int ringCount = std::max(0, (int)bonds.size() - (int)atomXs.size() + 1);
        results["RingCount"] = ringCount;
        
        // 94. AcidicCenterCount
        int acidicCenterCount = acidicAtoms.size();
        results["AcidicCenterCount"] = acidicCenterCount;
        
        // 95. BasicCenterCount
        int basicCenterCount = basicAtoms.size();
        results["BasicCenterCount"] = basicCenterCount;
        
        // 96. AcidicToBasicRatio
        double acidicToBasicRatio = basicCenterCount > 0 ? (double)acidicCenterCount / basicCenterCount : 0.0;
        results["AcidicToBasicRatio"] = acidicToBasicRatio;
        
        // 97. MeanAcidicAcidicDist
        double meanAcidicAcidicDist = 0.0;
        int acidicPairCount = 0;
        for (size_t i = 0; i < acidicAtoms.size(); ++i) {
            for (size_t j = i+1; j < acidicAtoms.size(); ++j) {
                int idx1 = acidicAtoms[i], idx2 = acidicAtoms[j];
                double dx = atomXs[idx1] - atomXs[idx2];
                double dy = atomYs[idx1] - atomYs[idx2];
                meanAcidicAcidicDist += std::sqrt(dx*dx + dy*dy);
                acidicPairCount++;
            }
        }
        meanAcidicAcidicDist = acidicPairCount > 0 ? meanAcidicAcidicDist / acidicPairCount : 0.0;
        results["MeanAcidicAcidicDist"] = meanAcidicAcidicDist;
        
        // 98. MeanBasicBasicDist
        double meanBasicBasicDist = 0.0;
        int basicPairCount = 0;
        for (size_t i = 0; i < basicAtoms.size(); ++i) {
            for (size_t j = i+1; j < basicAtoms.size(); ++j) {
                int idx1 = basicAtoms[i], idx2 = basicAtoms[j];
                double dx = atomXs[idx1] - atomXs[idx2];
                double dy = atomYs[idx1] - atomYs[idx2];
                meanBasicBasicDist += std::sqrt(dx*dx + dy*dy);
                basicPairCount++;
            }
        }
        meanBasicBasicDist = basicPairCount > 0 ? meanBasicBasicDist / basicPairCount : 0.0;
        results["MeanBasicBasicDist"] = meanBasicBasicDist;
        
        // 99. MeanAcidicBasicDist
        double meanAcidicBasicDist = 0.0;
        int acidicBasicPairCount = 0;
        for (int idxA : acidicAtoms) {
            for (int idxB : basicAtoms) {
                double dx = atomXs[idxA] - atomXs[idxB];
                double dy = atomYs[idxA] - atomYs[idxB];
                meanAcidicBasicDist += std::sqrt(dx*dx + dy*dy);
                acidicBasicPairCount++;
            }
        }
        meanAcidicBasicDist = acidicBasicPairCount > 0 ? meanAcidicBasicDist / acidicBasicPairCount : 0.0;
        results["MeanAcidicBasicDist"] = meanAcidicBasicDist;
        
        // 100. MinAcidicBasicDist
        double minAcidicBasicDist = std::numeric_limits<double>::max();
        for (int idxA : acidicAtoms) {
            for (int idxB : basicAtoms) {
                double dx = atomXs[idxA] - atomXs[idxB];
                double dy = atomYs[idxA] - atomYs[idxB];
                double dist = std::sqrt(dx*dx + dy*dy);
                if (dist < minAcidicBasicDist) minAcidicBasicDist = dist;
            }
        }
        if (minAcidicBasicDist == std::numeric_limits<double>::max()) minAcidicBasicDist = 0.0;
        results["MinAcidicBasicDist"] = minAcidicBasicDist;
        
        // 101. MeanAtomsNearAcidic
        double bondLengthRadius = meanBondLen * 2.0;
        double meanAtomsNearAcidic = 0.0;
        for (int idxA : acidicAtoms) {
            int count = 0;
            for (size_t i = 0; i < atomXs.size(); ++i) {
                if (idxA == i) continue;
                double dx = atomXs[idxA] - atomXs[i];
                double dy = atomYs[idxA] - atomYs[i];
                double dist = std::sqrt(dx*dx + dy*dy);
                if (dist <= bondLengthRadius) count++;
            }
            meanAtomsNearAcidic += count;
        }
        meanAtomsNearAcidic = acidicAtoms.empty() ? 0.0 : meanAtomsNearAcidic / acidicAtoms.size();
        results["MeanAtomsNearAcidic"] = meanAtomsNearAcidic;
        
        // 102. MeanAtomsNearBasic
        double meanAtomsNearBasic = 0.0;
        for (int idxB : basicAtoms) {
            int count = 0;
            for (size_t i = 0; i < atomXs.size(); ++i) {
                if (idxB == i) continue;
                double dx = atomXs[idxB] - atomXs[i];
                double dy = atomYs[idxB] - atomYs[i];
                double dist = std::sqrt(dx*dx + dy*dy);
                if (dist <= bondLengthRadius) count++;
            }
            meanAtomsNearBasic += count;
        }
        meanAtomsNearBasic = basicAtoms.empty() ? 0.0 : meanAtomsNearBasic / basicAtoms.size();
        results["MeanAtomsNearBasic"] = meanAtomsNearBasic;
        
        // 103. MeanAcidicRadius
        double meanAcidicRadius = 0.0;
        for (int idx : acidicAtoms) meanAcidicRadius += atomRadii[idx];
        meanAcidicRadius = acidicAtoms.empty() ? 0.0 : meanAcidicRadius / acidicAtoms.size();
        results["MeanAcidicRadius"] = meanAcidicRadius;
        
        // 104. MeanBasicRadius
        double meanBasicRadius = 0.0;
        for (int idx : basicAtoms) meanBasicRadius += atomRadii[idx];
        meanBasicRadius = basicAtoms.empty() ? 0.0 : meanBasicRadius / basicAtoms.size();
        results["MeanBasicRadius"] = meanBasicRadius;
        
        // 105. MeanAcidicCentroidDist
        double meanAcidicCentroidDist = 0.0;
        for (int idx : acidicAtoms) {
            double dx = atomXs[idx] - cxCentroid;
            double dy = atomYs[idx] - cyCentroid;
            meanAcidicCentroidDist += std::sqrt(dx*dx + dy*dy);
        }
        meanAcidicCentroidDist = acidicAtoms.empty() ? 0.0 : meanAcidicCentroidDist / acidicAtoms.size();
        results["MeanAcidicCentroidDist"] = meanAcidicCentroidDist;
        
        // 106. MeanBasicCentroidDist
        double meanBasicCentroidDist = 0.0;
        for (int idx : basicAtoms) {
            double dx = atomXs[idx] - cxCentroid;
            double dy = atomYs[idx] - cyCentroid;
            meanBasicCentroidDist += std::sqrt(dx*dx + dy*dy);
        }
        meanBasicCentroidDist = basicAtoms.empty() ? 0.0 : meanBasicCentroidDist / basicAtoms.size();
        results["MeanBasicCentroidDist"] = meanBasicCentroidDist;
        
        // 107. FracAcidicOnHull
        int acidicOnHull = 0;
        for (int idx : acidicAtoms) {
            for (const auto& pt : hull) {
                if (std::abs(pt.first - atomXs[idx]) < 1e-6 && std::abs(pt.second - atomYs[idx]) < 1e-6) {
                    acidicOnHull++;
                    break;
                }
            }
        }
        double fracAcidicOnHull = acidicAtoms.empty() ? 0.0 : (double)acidicOnHull / acidicAtoms.size();
        results["FracAcidicOnHull"] = fracAcidicOnHull;
        
        // 108. FracBasicOnHull
        int basicOnHull = 0;
        for (int idx : basicAtoms) {
            for (const auto& pt : hull) {
                if (std::abs(pt.first - atomXs[idx]) < 1e-6 && std::abs(pt.second - atomYs[idx]) < 1e-6) {
                    basicOnHull++;
                    break;
                }
            }
        }
        double fracBasicOnHull = basicAtoms.empty() ? 0.0 : (double)basicOnHull / basicAtoms.size();
        results["FracBasicOnHull"] = fracBasicOnHull;
        
        // 109. MeanAcidicAngle
        std::vector<double> acidicAngles;
        for (int idx : acidicAtoms) {
            std::vector<std::pair<double, double>> neighbors;
            for (const auto& b : bonds) {
                if (b.atom1Idx == idx)
                    neighbors.push_back({atomXs[b.atom2Idx] - atomXs[idx], atomYs[b.atom2Idx] - atomYs[idx]});
                else if (b.atom2Idx == idx)
                    neighbors.push_back({atomXs[b.atom1Idx] - atomXs[idx], atomYs[b.atom1Idx] - atomYs[idx]});
            }
            
            for (size_t j = 0; j + 1 < neighbors.size(); ++j) {
                for (size_t k = j + 1; k < neighbors.size(); ++k) {
                    double dx1 = neighbors[j].first, dy1 = neighbors[j].second;
                    double dx2 = neighbors[k].first, dy2 = neighbors[k].second;
                    double dot = dx1*dx2 + dy1*dy2;
                    double norm1 = std::sqrt(dx1*dx1 + dy1*dy1);
                    double norm2 = std::sqrt(dx2*dx2 + dy2*dy2);
                    
                    if (norm1 > 0 && norm2 > 0) {
                        double cosAngle = dot/(norm1*norm2);
                        cosAngle = std::max(-1.0, std::min(1.0, cosAngle));
                        double angle = std::acos(cosAngle) * 180.0 / M_PI;
                        acidicAngles.push_back(angle);
                    }
                }
            }
        }
        double meanAcidicAngle = !acidicAngles.empty() ?
            std::accumulate(acidicAngles.begin(), acidicAngles.end(), 0.0) / acidicAngles.size() : 0.0;
        results["MeanAcidicAngle"] = meanAcidicAngle;
        
        // 110. MeanBasicAngle
        std::vector<double> basicAngles;
        for (int idx : basicAtoms) {
            std::vector<std::pair<double, double>> neighbors;
            for (const auto& b : bonds) {
                if (b.atom1Idx == idx)
                    neighbors.push_back({atomXs[b.atom2Idx] - atomXs[idx], atomYs[b.atom2Idx] - atomYs[idx]});
                else if (b.atom2Idx == idx)
                    neighbors.push_back({atomXs[b.atom1Idx] - atomXs[idx], atomYs[b.atom1Idx] - atomYs[idx]});
            }
            
            for (size_t j = 0; j + 1 < neighbors.size(); ++j) {
                for (size_t k = j + 1; k < neighbors.size(); ++k) {
                    double dx1 = neighbors[j].first, dy1 = neighbors[j].second;
                    double dx2 = neighbors[k].first, dy2 = neighbors[k].second;
                    double dot = dx1*dx2 + dy1*dy2;
                    double norm1 = std::sqrt(dx1*dx1 + dy1*dy1);
                    double norm2 = std::sqrt(dx2*dx2 + dy2*dy2);
                    
                    if (norm1 > 0 && norm2 > 0) {
                        double cosAngle = dot/(norm1*norm2);
                        cosAngle = std::max(-1.0, std::min(1.0, cosAngle));
                        double angle = std::acos(cosAngle) * 180.0 / M_PI;
                        basicAngles.push_back(angle);
                    }
                }
            }
        }
        double meanBasicAngle = !basicAngles.empty() ?
            std::accumulate(basicAngles.begin(), basicAngles.end(), 0.0) / basicAngles.size() : 0.0;
        results["MeanBasicAngle"] = meanBasicAngle;
        
        // 111. MeanAcidicLuminance
        double meanAcidicLuminance = 0.0;
        for (int idx : acidicAtoms) meanAcidicLuminance += atomLuminances[idx];
        meanAcidicLuminance = acidicAtoms.empty() ? 0.0 : meanAcidicLuminance / acidicAtoms.size();
        results["MeanAcidicLuminance"] = meanAcidicLuminance;
        
        // 112. MeanBasicLuminance
        double meanBasicLuminance = 0.0;
        for (int idx : basicAtoms) meanBasicLuminance += atomLuminances[idx];
        meanBasicLuminance = basicAtoms.empty() ? 0.0 : meanBasicLuminance / basicAtoms.size();
        results["MeanBasicLuminance"] = meanBasicLuminance;
        
        // 113. AcidicBasicLuminanceDiff
        double acidicBasicLuminanceDiff = std::abs(meanAcidicLuminance - meanBasicLuminance);
        results["AcidicBasicLuminanceDiff"] = acidicBasicLuminanceDiff;
        
        // 114. FractalDimension
        double fractalDimension = 0.0;
        {
            std::vector<int> boxCounts;
            std::vector<double> boxSizes = {0.1, 0.2, 0.4, 0.8, 1.6};
            
            for (double boxSize : boxSizes) {
                std::set<std::pair<int, int>> boxes;
                
                for (size_t i = 0; i < atomXs.size(); i++) {
                    int boxX = (int)(atomXs[i] / boxSize);
                    int boxY = (int)(atomYs[i] / boxSize);
                    boxes.insert({boxX, boxY});
                }
                boxCounts.push_back(boxes.size());
            }
            
            // Linear regression for fractal dimension
            double sumX = 0.0, sumY = 0.0, sumXY = 0.0, sumX2 = 0.0;
            int n = boxSizes.size();
            for (int i = 0; i < n; i++) {
                double x = std::log(1.0/boxSizes[i]);
                double y = std::log((double)boxCounts[i]);
                sumX += x;
                sumY += y;
                sumXY += x * y;
                sumX2 += x * x;
            }
            
            if (n > 0 && (n * sumX2 - sumX * sumX) != 0) {
                fractalDimension = (n * sumXY - sumX * sumY) / (n * sumX2 - sumX * sumX);
                if (std::isnan(fractalDimension)) fractalDimension = 0.0;
            }
        }
        results["FractalDimension"] = fractalDimension;
        
        // 115. ForegroundRatio
        double foregroundRatio = bboxArea > 0 ? totalAtomArea / bboxArea : 0.0;
        results["ForegroundRatio"] = foregroundRatio;
        
        // 116. AverageColor
        double averageColor = meanLuminance;
        results["AverageColor"] = averageColor;
        
        // 117. ColorVariance
        double colorVariance = atomLuminanceStd * atomLuminanceStd;
        results["ColorVariance"] = colorVariance;
        
        // 118. ImageCenterX
        double imageCenterX = cxCentroid;
        results["ImageCenterX"] = imageCenterX;
        
        // 119. ImageCenterY
        double imageCenterY = cyCentroid;
        results["ImageCenterY"] = imageCenterY;
        
        // 120. ImageOrientation
        double xxSum = 0.0, xySum = 0.0, yySum = 0.0;
        for (size_t i = 0; i < atomXs.size(); i++) {
            double x = atomXs[i] - imageCenterX;
            double y = atomYs[i] - imageCenterY;
            xxSum += x * x;
            xySum += x * y;
            yySum += y * y;
        }
        double orientation = 0.0;
        if (xxSum != yySum || xySum != 0) {
            orientation = 0.5 * std::atan2(2 * xySum, xxSum - yySum) * 180.0 / M_PI;
        }
        results["ImageOrientation"] = orientation;
        
        // Return all calculated descriptors
        return results;
    }
    catch (const std::exception& e) {
        globalLogger.error("ImageDescriptorBase: Exception during calculation: " + std::string(e.what()));
        return results;
    }
}

// Define the individual descriptor classes

// 1. BoundingBoxArea
BoundingBoxAreaDescriptor::BoundingBoxAreaDescriptor()
    : ImageDescriptorBase("BoundingBoxArea", "The area of the molecule's bounding box") {}

std::variant<double, int, std::string> BoundingBoxAreaDescriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("BoundingBoxArea") ? descriptors["BoundingBoxArea"] : 0.0;
}

// 2. MoleculeWidth
MoleculeWidthDescriptor::MoleculeWidthDescriptor()
    : ImageDescriptorBase("MoleculeWidth", "The width of the molecule's bounding box") {}

std::variant<double, int, std::string> MoleculeWidthDescriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("MoleculeWidth") ? descriptors["MoleculeWidth"] : 0.0;
}

// 3. MoleculeHeight
MoleculeHeightDescriptor::MoleculeHeightDescriptor()
    : ImageDescriptorBase("MoleculeHeight", "The height of the molecule's bounding box") {}

std::variant<double, int, std::string> MoleculeHeightDescriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("MoleculeHeight") ? descriptors["MoleculeHeight"] : 0.0;
}

// 4. AspectRatio
AspectRatioDescriptor::AspectRatioDescriptor()
    : ImageDescriptorBase("AspectRatio", "The width-to-height ratio of the molecule's bounding box") {}

std::variant<double, int, std::string> AspectRatioDescriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("AspectRatio") ? descriptors["AspectRatio"] : 0.0;
}

// 5. AtomDensity
AtomDensityDescriptor::AtomDensityDescriptor()
    : ImageDescriptorBase("AtomDensity", "The number of atoms per unit area in the molecule's bounding box") {}

std::variant<double, int, std::string> AtomDensityDescriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("AtomDensity") ? descriptors["AtomDensity"] : 0.0;
}

// 6. MeanAtomRadius
MeanAtomRadiusDescriptor::MeanAtomRadiusDescriptor()
    : ImageDescriptorBase("MeanAtomRadius", "The average radius of atoms in the molecule") {}

std::variant<double, int, std::string> MeanAtomRadiusDescriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("MeanAtomRadius") ? descriptors["MeanAtomRadius"] : 0.0;
}

// 7. MaxAtomRadius
MaxAtomRadiusDescriptor::MaxAtomRadiusDescriptor()
    : ImageDescriptorBase("MaxAtomRadius", "The maximum radius of atoms in the molecule") {}

std::variant<double, int, std::string> MaxAtomRadiusDescriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("MaxAtomRadius") ? descriptors["MaxAtomRadius"] : 0.0;
}

// 8. MinAtomRadius
MinAtomRadiusDescriptor::MinAtomRadiusDescriptor()
    : ImageDescriptorBase("MinAtomRadius", "The minimum radius of atoms in the molecule") {}

std::variant<double, int, std::string> MinAtomRadiusDescriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("MinAtomRadius") ? descriptors["MinAtomRadius"] : 0.0;
}

// 9. MeanBondLength
MeanBondLengthDescriptor::MeanBondLengthDescriptor()
    : ImageDescriptorBase("MeanBondLength", "The average length of bonds in the molecule") {}

std::variant<double, int, std::string> MeanBondLengthDescriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("MeanBondLength") ? descriptors["MeanBondLength"] : 0.0;
}

// 10. MaxBondLength
MaxBondLengthDescriptor::MaxBondLengthDescriptor()
    : ImageDescriptorBase("MaxBondLength", "The maximum length of bonds in the molecule") {}

std::variant<double, int, std::string> MaxBondLengthDescriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("MaxBondLength") ? descriptors["MaxBondLength"] : 0.0;
}

// 11. MinBondLength
MinBondLengthDescriptor::MinBondLengthDescriptor()
    : ImageDescriptorBase("MinBondLength", "The minimum length of bonds in the molecule") {}

std::variant<double, int, std::string> MinBondLengthDescriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("MinBondLength") ? descriptors["MinBondLength"] : 0.0;
}

// 12. BondLengthStd
BondLengthStdDescriptor::BondLengthStdDescriptor()
    : ImageDescriptorBase("BondLengthStd", "The standard deviation of bond lengths in the molecule") {}

std::variant<double, int, std::string> BondLengthStdDescriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("BondLengthStd") ? descriptors["BondLengthStd"] : 0.0;
}

// 13. MeanAtomAtomDist
MeanAtomAtomDistDescriptor::MeanAtomAtomDistDescriptor()
    : ImageDescriptorBase("MeanAtomAtomDist", "The average distance between all pairs of atoms in the molecule") {}

std::variant<double, int, std::string> MeanAtomAtomDistDescriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("MeanAtomAtomDist") ? descriptors["MeanAtomAtomDist"] : 0.0;
}

// 14. MeanAtomLuminance
MeanAtomLuminanceDescriptor::MeanAtomLuminanceDescriptor()
    : ImageDescriptorBase("MeanAtomLuminance", "The average luminance of atoms in the molecule") {}

std::variant<double, int, std::string> MeanAtomLuminanceDescriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("MeanAtomLuminance") ? descriptors["MeanAtomLuminance"] : 0.0;
}

// 15. MaxAtomAtomDist
MaxAtomAtomDistDescriptor::MaxAtomAtomDistDescriptor()
    : ImageDescriptorBase("MaxAtomAtomDist", "The maximum distance between any pair of atoms in the molecule") {}

std::variant<double, int, std::string> MaxAtomAtomDistDescriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("MaxAtomAtomDist") ? descriptors["MaxAtomAtomDist"] : 0.0;
}

// 16. MinAtomAtomDist
MinAtomAtomDistDescriptor::MinAtomAtomDistDescriptor()
    : ImageDescriptorBase("MinAtomAtomDist", "The minimum distance between any pair of atoms in the molecule") {}

std::variant<double, int, std::string> MinAtomAtomDistDescriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("MinAtomAtomDist") ? descriptors["MinAtomAtomDist"] : 0.0;
}

// 17. MeanAtomArea
MeanAtomAreaDescriptor::MeanAtomAreaDescriptor()
    : ImageDescriptorBase("MeanAtomArea", "The average area of atoms in the molecule") {}

std::variant<double, int, std::string> MeanAtomAreaDescriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("MeanAtomArea") ? descriptors["MeanAtomArea"] : 0.0;
}

// 18. MedianAtomRadius
MedianAtomRadiusDescriptor::MedianAtomRadiusDescriptor()
    : ImageDescriptorBase("MedianAtomRadius", "The median radius of atoms in the molecule") {}

std::variant<double, int, std::string> MedianAtomRadiusDescriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("MedianAtomRadius") ? descriptors["MedianAtomRadius"] : 0.0;
}

// 19. MedianBondLength
MedianBondLengthDescriptor::MedianBondLengthDescriptor()
    : ImageDescriptorBase("MedianBondLength", "The median length of bonds in the molecule") {}

std::variant<double, int, std::string> MedianBondLengthDescriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("MedianBondLength") ? descriptors["MedianBondLength"] : 0.0;
}

// 20. AtomRadiusRange
AtomRadiusRangeDescriptor::AtomRadiusRangeDescriptor()
    : ImageDescriptorBase("AtomRadiusRange", "The range of atom radii in the molecule") {}

std::variant<double, int, std::string> AtomRadiusRangeDescriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("AtomRadiusRange") ? descriptors["AtomRadiusRange"] : 0.0;
}

// 21. AtomAreaFraction
AtomAreaFractionDescriptor::AtomAreaFractionDescriptor()
    : ImageDescriptorBase("AtomAreaFraction", "The fraction of the bounding box area occupied by atoms") {}

std::variant<double, int, std::string> AtomAreaFractionDescriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("AtomAreaFraction") ? descriptors["AtomAreaFraction"] : 0.0;
}

// 22. BondCoverageFraction
BondCoverageFractionDescriptor::BondCoverageFractionDescriptor()
    : ImageDescriptorBase("BondCoverageFraction", "The fraction of the bounding box perimeter covered by bonds") {}

std::variant<double, int, std::string> BondCoverageFractionDescriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("BondCoverageFraction") ? descriptors["BondCoverageFraction"] : 0.0;
}

// 23. AtomPackingDensity
AtomPackingDensityDescriptor::AtomPackingDensityDescriptor()
    : ImageDescriptorBase("AtomPackingDensity", "The density of atom packing within the convex hull") {}

std::variant<double, int, std::string> AtomPackingDensityDescriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("AtomPackingDensity") ? descriptors["AtomPackingDensity"] : 0.0;
}

// 24. AtomMassX
AtomMassXDescriptor::AtomMassXDescriptor()
    : ImageDescriptorBase("AtomMassX", "The x-coordinate of the molecule's center of mass") {}

std::variant<double, int, std::string> AtomMassXDescriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("AtomMassX") ? descriptors["AtomMassX"] : 0.0;
}

// 25. AtomMassY
AtomMassYDescriptor::AtomMassYDescriptor()
    : ImageDescriptorBase("AtomMassY", "The y-coordinate of the molecule's center of mass") {}

std::variant<double, int, std::string> AtomMassYDescriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("AtomMassY") ? descriptors["AtomMassY"] : 0.0;
}

// 26. AtomMassDist
AtomMassDistDescriptor::AtomMassDistDescriptor()
    : ImageDescriptorBase("AtomMassDist", "The distance of the molecule's center of mass from the image center") {}

std::variant<double, int, std::string> AtomMassDistDescriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("AtomMassDist") ? descriptors["AtomMassDist"] : 0.0;
}

// 27. AtomRadiusStd
AtomRadiusStdDescriptor::AtomRadiusStdDescriptor()
    : ImageDescriptorBase("AtomRadiusStd", "The standard deviation of atom radii in the molecule") {}

std::variant<double, int, std::string> AtomRadiusStdDescriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("AtomRadiusStd") ? descriptors["AtomRadiusStd"] : 0.0;
}

// 28. MeanBondAngle
MeanBondAngleDescriptor::MeanBondAngleDescriptor()
    : ImageDescriptorBase("MeanBondAngle", "The average angle between bonds in the molecule") {}

std::variant<double, int, std::string> MeanBondAngleDescriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("MeanBondAngle") ? descriptors["MeanBondAngle"] : 0.0;
}

// 29. BondAngleStd
BondAngleStdDescriptor::BondAngleStdDescriptor()
    : ImageDescriptorBase("BondAngleStd", "The standard deviation of bond angles in the molecule") {}

std::variant<double, int, std::string> BondAngleStdDescriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("BondAngleStd") ? descriptors["BondAngleStd"] : 0.0;
}

// 30. AtomXStd
AtomXStdDescriptor::AtomXStdDescriptor()
    : ImageDescriptorBase("AtomXStd", "The standard deviation of atom x-coordinates in the molecule") {}

std::variant<double, int, std::string> AtomXStdDescriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("AtomXStd") ? descriptors["AtomXStd"] : 0.0;
}

// 31. AtomYStd
AtomYStdDescriptor::AtomYStdDescriptor()
    : ImageDescriptorBase("AtomYStd", "The standard deviation of atom y-coordinates in the molecule") {}

std::variant<double, int, std::string> AtomYStdDescriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("AtomYStd") ? descriptors["AtomYStd"] : 0.0;
}

// 32. AtomLuminanceStd
AtomLuminanceStdDescriptor::AtomLuminanceStdDescriptor()
    : ImageDescriptorBase("AtomLuminanceStd", "The standard deviation of atom luminances in the molecule") {}

std::variant<double, int, std::string> AtomLuminanceStdDescriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("AtomLuminanceStd") ? descriptors["AtomLuminanceStd"] : 0.0;
}

// 33. BondLenMAD
BondLenMADDescriptor::BondLenMADDescriptor()
    : ImageDescriptorBase("BondLenMAD", "The median absolute deviation of bond lengths in the molecule") {}

std::variant<double, int, std::string> BondLenMADDescriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("BondLenMAD") ? descriptors["BondLenMAD"] : 0.0;
}

// 34. AtomRadiusMAD
AtomRadiusMADDescriptor::AtomRadiusMADDescriptor()
    : ImageDescriptorBase("AtomRadiusMAD", "The median absolute deviation of atom radii in the molecule") {}

std::variant<double, int, std::string> AtomRadiusMADDescriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("AtomRadiusMAD") ? descriptors["AtomRadiusMAD"] : 0.0;
}

// 35. AtomAreaStd
AtomAreaStdDescriptor::AtomAreaStdDescriptor()
    : ImageDescriptorBase("AtomAreaStd", "The standard deviation of atom areas in the molecule") {}

std::variant<double, int, std::string> AtomAreaStdDescriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("AtomAreaStd") ? descriptors["AtomAreaStd"] : 0.0;
}

// 36. AtomRadiusCV
AtomRadiusCVDescriptor::AtomRadiusCVDescriptor()
    : ImageDescriptorBase("AtomRadiusCV", "The coefficient of variation of atom radii in the molecule") {}

std::variant<double, int, std::string> AtomRadiusCVDescriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("AtomRadiusCV") ? descriptors["AtomRadiusCV"] : 0.0;
}

// 37. AtomAreaCV
AtomAreaCVDescriptor::AtomAreaCVDescriptor()
    : ImageDescriptorBase("AtomAreaCV", "The coefficient of variation of atom areas in the molecule") {}

std::variant<double, int, std::string> AtomAreaCVDescriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("AtomAreaCV") ? descriptors["AtomAreaCV"] : 0.0;
}

// 38. AtomLuminanceCV
AtomLuminanceCVDescriptor::AtomLuminanceCVDescriptor()
    : ImageDescriptorBase("AtomLuminanceCV", "The coefficient of variation of atom luminances in the molecule") {}

std::variant<double, int, std::string> AtomLuminanceCVDescriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("AtomLuminanceCV") ? descriptors["AtomLuminanceCV"] : 0.0;
}

// 39. AtomAreaRange
AtomAreaRangeDescriptor::AtomAreaRangeDescriptor()
    : ImageDescriptorBase("AtomAreaRange", "The range of atom areas in the molecule") {}

std::variant<double, int, std::string> AtomAreaRangeDescriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("AtomAreaRange") ? descriptors["AtomAreaRange"] : 0.0;
}

// 40. AtomAreaMedian
AtomAreaMedianDescriptor::AtomAreaMedianDescriptor()
    : ImageDescriptorBase("AtomAreaMedian", "The median area of atoms in the molecule") {}

std::variant<double, int, std::string> AtomAreaMedianDescriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("AtomAreaMedian") ? descriptors["AtomAreaMedian"] : 0.0;
}

// 41. AtomAreaMAD
AtomAreaMADDescriptor::AtomAreaMADDescriptor()
    : ImageDescriptorBase("AtomAreaMAD", "The median absolute deviation of atom areas in the molecule") {}

std::variant<double, int, std::string> AtomAreaMADDescriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("AtomAreaMAD") ? descriptors["AtomAreaMAD"] : 0.0;
}

// 42. AtomAreaHullFrac
AtomAreaHullFracDescriptor::AtomAreaHullFracDescriptor()
    : ImageDescriptorBase("AtomAreaHullFrac", "The fraction of total atom area contributed by atoms on the convex hull") {}

std::variant<double, int, std::string> AtomAreaHullFracDescriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("AtomAreaHullFrac") ? descriptors["AtomAreaHullFrac"] : 0.0;
}

// 43. AtomAreaCenterFrac
AtomAreaCenterFracDescriptor::AtomAreaCenterFracDescriptor()
    : ImageDescriptorBase("AtomAreaCenterFrac", "The fraction of total atom area contributed by atoms in the central region") {}

std::variant<double, int, std::string> AtomAreaCenterFracDescriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("AtomAreaCenterFrac") ? descriptors["AtomAreaCenterFrac"] : 0.0;
}

// 44. StdAllAtomDist
StdAllAtomDistDescriptor::StdAllAtomDistDescriptor()
    : ImageDescriptorBase("StdAllAtomDist", "The standard deviation of all atom-atom distances in the molecule") {}

std::variant<double, int, std::string> StdAllAtomDistDescriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("StdAllAtomDist") ? descriptors["StdAllAtomDist"] : 0.0;
}

// 45. MinBondAngle
MinBondAngleDescriptor::MinBondAngleDescriptor()
    : ImageDescriptorBase("MinBondAngle", "The minimum angle between bonds in the molecule") {}

std::variant<double, int, std::string> MinBondAngleDescriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("MinBondAngle") ? descriptors["MinBondAngle"] : 0.0;
}

// 46. MaxBondAngle
MaxBondAngleDescriptor::MaxBondAngleDescriptor()
    : ImageDescriptorBase("MaxBondAngle", "The maximum angle between bonds in the molecule") {}

std::variant<double, int, std::string> MaxBondAngleDescriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("MaxBondAngle") ? descriptors["MaxBondAngle"] : 0.0;
}

// 47. MedianBondAngle
MedianBondAngleDescriptor::MedianBondAngleDescriptor()
    : ImageDescriptorBase("MedianBondAngle", "The median angle between bonds in the molecule") {}

std::variant<double, int, std::string> MedianBondAngleDescriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("MedianBondAngle") ? descriptors["MedianBondAngle"] : 0.0;
}

// 48. MeanBondColorDiff
MeanBondColorDiffDescriptor::MeanBondColorDiffDescriptor()
    : ImageDescriptorBase("MeanBondColorDiff", "The average color difference between bonded atoms") {}

std::variant<double, int, std::string> MeanBondColorDiffDescriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("MeanBondColorDiff") ? descriptors["MeanBondColorDiff"] : 0.0;
}

// 49. MeanBondLuminanceDiff
MeanBondLuminanceDiffDescriptor::MeanBondLuminanceDiffDescriptor()
    : ImageDescriptorBase("MeanBondLuminanceDiff", "The average luminance difference between bonded atoms") {}

std::variant<double, int, std::string> MeanBondLuminanceDiffDescriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("MeanBondLuminanceDiff") ? descriptors["MeanBondLuminanceDiff"] : 0.0;
}

// 50. MeanBondLenDiffColor
MeanBondLenDiffColorDescriptor::MeanBondLenDiffColorDescriptor()
    : ImageDescriptorBase("MeanBondLenDiffColor", "The average length of bonds between atoms of different colors") {}

std::variant<double, int, std::string> MeanBondLenDiffColorDescriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("MeanBondLenDiffColor") ? descriptors["MeanBondLenDiffColor"] : 0.0;
}

// 51. MeanBondAngleDiffColor
MeanBondAngleDiffColorDescriptor::MeanBondAngleDiffColorDescriptor()
    : ImageDescriptorBase("MeanBondAngleDiffColor", "The average angle between bonds of different colors") {}

std::variant<double, int, std::string> MeanBondAngleDiffColorDescriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("MeanBondAngleDiffColor") ? descriptors["MeanBondAngleDiffColor"] : 0.0;
}

// 52. MeanDistSameColor
MeanDistSameColorDescriptor::MeanDistSameColorDescriptor()
    : ImageDescriptorBase("MeanDistSameColor", "The average distance between atoms of the same color") {}

std::variant<double, int, std::string> MeanDistSameColorDescriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("MeanDistSameColor") ? descriptors["MeanDistSameColor"] : 0.0;
}

// 53. MeanDistDiffColor
MeanDistDiffColorDescriptor::MeanDistDiffColorDescriptor()
    : ImageDescriptorBase("MeanDistDiffColor", "The average distance between atoms of different colors") {}

std::variant<double, int, std::string> MeanDistDiffColorDescriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("MeanDistDiffColor") ? descriptors["MeanDistDiffColor"] : 0.0;
}

// 54. MeanBondOrder
MeanBondOrderDescriptor::MeanBondOrderDescriptor()
    : ImageDescriptorBase("MeanBondOrder", "The average bond order in the molecule") {}

std::variant<double, int, std::string> MeanBondOrderDescriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("MeanBondOrder") ? descriptors["MeanBondOrder"] : 0.0;
}

// 55. FracDoubleBonds
FracDoubleBondsDescriptor::FracDoubleBondsDescriptor()
    : ImageDescriptorBase("FracDoubleBonds", "The fraction of double bonds in the molecule") {}

std::variant<double, int, std::string> FracDoubleBondsDescriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("FracDoubleBonds") ? descriptors["FracDoubleBonds"] : 0.0;
}

// 56. FracTripleBonds
FracTripleBondsDescriptor::FracTripleBondsDescriptor()
    : ImageDescriptorBase("FracTripleBonds", "The fraction of triple bonds in the molecule") {}

std::variant<double, int, std::string> FracTripleBondsDescriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("FracTripleBonds") ? descriptors["FracTripleBonds"] : 0.0;
}

// 57. FracAromaticBonds
FracAromaticBondsDescriptor::FracAromaticBondsDescriptor()
    : ImageDescriptorBase("FracAromaticBonds", "The fraction of aromatic bonds in the molecule") {}

std::variant<double, int, std::string> FracAromaticBondsDescriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("FracAromaticBonds") ? descriptors["FracAromaticBonds"] : 0.0;
}

// 58. MeanAtomDegree
MeanAtomDegreeDescriptor::MeanAtomDegreeDescriptor()
    : ImageDescriptorBase("MeanAtomDegree", "The average degree of atoms in the molecule") {}

std::variant<double, int, std::string> MeanAtomDegreeDescriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("MeanAtomDegree") ? descriptors["MeanAtomDegree"] : 0.0;
}

// 59. MaxAtomDegree
MaxAtomDegreeDescriptor::MaxAtomDegreeDescriptor()
    : ImageDescriptorBase("MaxAtomDegree", "The maximum degree of atoms in the molecule") {}

std::variant<double, int, std::string> MaxAtomDegreeDescriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("MaxAtomDegree") ? descriptors["MaxAtomDegree"] : 0;
}

// 60. MinAtomDegree
MinAtomDegreeDescriptor::MinAtomDegreeDescriptor()
    : ImageDescriptorBase("MinAtomDegree", "The minimum degree of atoms in the molecule") {}

std::variant<double, int, std::string> MinAtomDegreeDescriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("MinAtomDegree") ? descriptors["MinAtomDegree"] : 0;
}

// 61. MeanDistToCentroid
MeanDistToCentroidDescriptor::MeanDistToCentroidDescriptor()
    : ImageDescriptorBase("MeanDistToCentroid", "The average distance of atoms to the centroid of the molecule") {}

std::variant<double, int, std::string> MeanDistToCentroidDescriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("MeanDistToCentroid") ? descriptors["MeanDistToCentroid"] : 0.0;
}

// 62. StdDistToCentroid
StdDistToCentroidDescriptor::StdDistToCentroidDescriptor()
    : ImageDescriptorBase("StdDistToCentroid", "The standard deviation of distances of atoms to the centroid of the molecule") {}

std::variant<double, int, std::string> StdDistToCentroidDescriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("StdDistToCentroid") ? descriptors["StdDistToCentroid"] : 0.0;
}

// 63. MeanLenSingle
MeanLenSingleDescriptor::MeanLenSingleDescriptor()
    : ImageDescriptorBase("MeanLenSingle", "The average length of single bonds in the molecule") {}

std::variant<double, int, std::string> MeanLenSingleDescriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("MeanLenSingle") ? descriptors["MeanLenSingle"] : 0.0;
}

// 64. MeanLenDouble
MeanLenDoubleDescriptor::MeanLenDoubleDescriptor()
    : ImageDescriptorBase("MeanLenDouble", "The average length of double bonds in the molecule") {}

std::variant<double, int, std::string> MeanLenDoubleDescriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("MeanLenDouble") ? descriptors["MeanLenDouble"] : 0.0;
}

// 65. MeanLenTriple
MeanLenTripleDescriptor::MeanLenTripleDescriptor()
    : ImageDescriptorBase("MeanLenTriple", "The average length of triple bonds in the molecule") {}

std::variant<double, int, std::string> MeanLenTripleDescriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("MeanLenTriple") ? descriptors["MeanLenTriple"] : 0.0;
}

// 66. MeanLenAromatic
MeanLenAromaticDescriptor::MeanLenAromaticDescriptor()
    : ImageDescriptorBase("MeanLenAromatic", "The average length of aromatic bonds in the molecule") {}

std::variant<double, int, std::string> MeanLenAromaticDescriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("MeanLenAromatic") ? descriptors["MeanLenAromatic"] : 0.0;
}

// 67. FracDeg1
FracDeg1Descriptor::FracDeg1Descriptor()
    : ImageDescriptorBase("FracDeg1", "The fraction of atoms with degree 1 in the molecule") {}

std::variant<double, int, std::string> FracDeg1Descriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("FracDeg1") ? descriptors["FracDeg1"] : 0.0;
}

// 68. FracDeg2
FracDeg2Descriptor::FracDeg2Descriptor()
    : ImageDescriptorBase("FracDeg2", "The fraction of atoms with degree 2 in the molecule") {}

std::variant<double, int, std::string> FracDeg2Descriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("FracDeg2") ? descriptors["FracDeg2"] : 0.0;
}

// 69. FracDeg3
FracDeg3Descriptor::FracDeg3Descriptor()
    : ImageDescriptorBase("FracDeg3", "The fraction of atoms with degree 3 in the molecule") {}

std::variant<double, int, std::string> FracDeg3Descriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("FracDeg3") ? descriptors["FracDeg3"] : 0.0;
}

// 70. FracDeg4
FracDeg4Descriptor::FracDeg4Descriptor()
    : ImageDescriptorBase("FracDeg4", "The fraction of atoms with degree 4 or higher in the molecule") {}

std::variant<double, int, std::string> FracDeg4Descriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("FracDeg4") ? descriptors["FracDeg4"] : 0.0;
}

// 71. MeanBondRadiiDiff
MeanBondRadiiDiffDescriptor::MeanBondRadiiDiffDescriptor()
    : ImageDescriptorBase("MeanBondRadiiDiff", "The average absolute difference in radii between bonded atoms") {}

std::variant<double, int, std::string> MeanBondRadiiDiffDescriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("MeanBondRadiiDiff") ? descriptors["MeanBondRadiiDiff"] : 0.0;
}

// 72. MeanBondLuminanceDiff2 (already calculated as MeanBondLuminanceDiff)
MeanBondLuminanceDiff2Descriptor::MeanBondLuminanceDiff2Descriptor()
    : ImageDescriptorBase("MeanBondLuminanceDiff2", "The average luminance difference between bonded atoms") {}

std::variant<double, int, std::string> MeanBondLuminanceDiff2Descriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("MeanBondLuminanceDiff2") ? descriptors["MeanBondLuminanceDiff2"] : 0.0;
}

// 73. MeanAngleHighDegree
MeanAngleHighDegreeDescriptor::MeanAngleHighDegreeDescriptor()
    : ImageDescriptorBase("MeanAngleHighDegree", "The average angle between bonds at high degree atoms") {}

std::variant<double, int, std::string> MeanAngleHighDegreeDescriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("MeanAngleHighDegree") ? descriptors["MeanAngleHighDegree"] : 0.0;
}

// 74. BoundaryAtomRatio (already calculated as fraction of atoms on hull)
BoundaryAtomRatioDescriptor::BoundaryAtomRatioDescriptor()
    : ImageDescriptorBase("BoundaryAtomRatio", "The fraction of atoms on the convex hull") {}

std::variant<double, int, std::string> BoundaryAtomRatioDescriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("BoundaryAtomRatio") ? descriptors["BoundaryAtomRatio"] : 0.0;
}

// 75. MeanAtomEccentricity
MeanAtomEccentricityDescriptor::MeanAtomEccentricityDescriptor()
    : ImageDescriptorBase("MeanAtomEccentricity", "The average eccentricity of atoms") {}

std::variant<double, int, std::string> MeanAtomEccentricityDescriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("MeanAtomEccentricity") ? descriptors["MeanAtomEccentricity"] : 0.0;
}

// 76. MolecularDiameter (maxAtomDist)
MolecularDiameterDescriptor::MolecularDiameterDescriptor()
    : ImageDescriptorBase("MolecularDiameter", "The diameter of the molecule") {}

std::variant<double, int, std::string> MolecularDiameterDescriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("MolecularDiameter") ? descriptors["MolecularDiameter"] : 0.0;
}

// 77. MolecularRadius (maxDistToCentroid)
MolecularRadiusDescriptor::MolecularRadiusDescriptor()
    : ImageDescriptorBase("MolecularRadius", "The radius of the molecule") {}

std::variant<double, int, std::string> MolecularRadiusDescriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("MolecularRadius") ? descriptors["MolecularRadius"] : 0.0;
}

// 78. RadiusOfGyration
RadiusOfGyrationDescriptor::RadiusOfGyrationDescriptor()
    : ImageDescriptorBase("RadiusOfGyration", "The radius of gyration of the molecule") {}

std::variant<double, int, std::string> RadiusOfGyrationDescriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("RadiusOfGyration") ? descriptors["RadiusOfGyration"] : 0.0;
}

// 79. MolecularSphericity
MolecularSphericityDescriptor::MolecularSphericityDescriptor()
    : ImageDescriptorBase("MolecularSphericity", "The sphericity of the molecule") {}

std::variant<double, int, std::string> MolecularSphericityDescriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("MolecularSphericity") ? descriptors["MolecularSphericity"] : 0.0;
}

// 80. BondLenToAtomRadiusRatio
BondLenToAtomRadiusRatioDescriptor::BondLenToAtomRadiusRatioDescriptor()
    : ImageDescriptorBase("BondLenToAtomRadiusRatio", "The ratio of mean bond length to mean atom radius") {}

std::variant<double, int, std::string> BondLenToAtomRadiusRatioDescriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("BondLenToAtomRadiusRatio") ? descriptors["BondLenToAtomRadiusRatio"] : 0.0;
}

// 81. EdgeDensity
EdgeDensityDescriptor::EdgeDensityDescriptor()
    : ImageDescriptorBase("EdgeDensity", "The edge density of the molecule") {}

std::variant<double, int, std::string> EdgeDensityDescriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("EdgeDensity") ? descriptors["EdgeDensity"] : 0.0;
}

// 82. PlanarityMeasure (trivially zero in 2D)
PlanarityMeasureDescriptor::PlanarityMeasureDescriptor()
    : ImageDescriptorBase("PlanarityMeasure", "The planarity measure of the molecule") {}

std::variant<double, int, std::string> PlanarityMeasureDescriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("PlanarityMeasure") ? descriptors["PlanarityMeasure"] : 0.0;
}

// 83. MeanNearestNeighborDist
MeanNearestNeighborDistDescriptor::MeanNearestNeighborDistDescriptor()
    : ImageDescriptorBase("MeanNearestNeighborDist", "The average distance between nearest neighbors") {}

std::variant<double, int, std::string> MeanNearestNeighborDistDescriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("MeanNearestNeighborDist") ? descriptors["MeanNearestNeighborDist"] : 0.0;
}

// 84. MeanAtomsInRadius
MeanAtomsInRadiusDescriptor::MeanAtomsInRadiusDescriptor()
    : ImageDescriptorBase("MeanAtomsInRadius", "The average number of atoms within a radius of mean bond length") {}

std::variant<double, int, std::string> MeanAtomsInRadiusDescriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("MeanAtomsInRadius") ? descriptors["MeanAtomsInRadius"] : 0.0;
}

// 85. VarNearestNeighborDist
VarNearestNeighborDistDescriptor::VarNearestNeighborDistDescriptor()
    : ImageDescriptorBase("VarNearestNeighborDist", "The variance of distances between nearest neighbors") {}

std::variant<double, int, std::string> VarNearestNeighborDistDescriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("VarNearestNeighborDist") ? descriptors["VarNearestNeighborDist"] : 0.0;
}

// 86. MeanBondToBondAngle
MeanBondToBondAngleDescriptor::MeanBondToBondAngleDescriptor()
    : ImageDescriptorBase("MeanBondToBondAngle", "The average angle between bonds") {}

std::variant<double, int, std::string> MeanBondToBondAngleDescriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("MeanBondToBondAngle") ? descriptors["MeanBondToBondAngle"] : 0.0;
}

// 87. MomentOfInertia
MomentOfInertiaDescriptor::MomentOfInertiaDescriptor()
    : ImageDescriptorBase("MomentOfInertia", "The moment of inertia of the molecule") {}

std::variant<double, int, std::string> MomentOfInertiaDescriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("MomentOfInertia") ? descriptors["MomentOfInertia"] : 0.0;
}

// 88. MeanAtomCentrality
MeanAtomCentralityDescriptor::MeanAtomCentralityDescriptor()
    : ImageDescriptorBase("MeanAtomCentrality", "The average centrality of atoms") {}

std::variant<double, int, std::string> MeanAtomCentralityDescriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("MeanAtomCentrality") ? descriptors["MeanAtomCentrality"] : 0.0;
}

// 89. TerminalAtomCount (= nDeg1)
TerminalAtomCountDescriptor::TerminalAtomCountDescriptor()
    : ImageDescriptorBase("TerminalAtomCount", "The number of terminal atoms (degree 1)") {}

std::variant<double, int, std::string> TerminalAtomCountDescriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("TerminalAtomCount") ? descriptors["TerminalAtomCount"] : 0;
}

// 90. JunctionAtomCount
JunctionAtomCountDescriptor::JunctionAtomCountDescriptor()
    : ImageDescriptorBase("JunctionAtomCount", "The number of junction atoms (degree >= 3)") {}

std::variant<double, int, std::string> JunctionAtomCountDescriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("JunctionAtomCount") ? descriptors["JunctionAtomCount"] : 0;
}

// 91. TerminalToJunctionRatio
TerminalToJunctionRatioDescriptor::TerminalToJunctionRatioDescriptor()
    : ImageDescriptorBase("TerminalToJunctionRatio", "The ratio of terminal atoms to junction atoms") {}

std::variant<double, int, std::string> TerminalToJunctionRatioDescriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("TerminalToJunctionRatio") ? descriptors["TerminalToJunctionRatio"] : 0.0;
}

// 92. WienerIndex
WienerIndexDescriptor::WienerIndexDescriptor()
    : ImageDescriptorBase("WienerIndex", "The Wiener index of the molecule") {}

std::variant<double, int, std::string> WienerIndexDescriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("WienerIndex") ? descriptors["WienerIndex"] : 0.0;
}

// 93. RingCount
RingCountDescriptor::RingCountDescriptor()
    : ImageDescriptorBase("RingCount", "The number of rings in the molecule") {}

std::variant<double, int, std::string> RingCountDescriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("RingCount") ? descriptors["RingCount"] : 0;
}

// 94. AcidicCenterCount
AcidicCenterCountDescriptor::AcidicCenterCountDescriptor()
    : ImageDescriptorBase("AcidicCenterCount", "The number of acidic centers") {}

std::variant<double, int, std::string> AcidicCenterCountDescriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("AcidicCenterCount") ? descriptors["AcidicCenterCount"] : 0;
}

// 95. BasicCenterCount
BasicCenterCountDescriptor::BasicCenterCountDescriptor()
    : ImageDescriptorBase("BasicCenterCount", "The number of basic centers") {}

std::variant<double, int, std::string> BasicCenterCountDescriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("BasicCenterCount") ? descriptors["BasicCenterCount"] : 0;
}

// 96. AcidicToBasicRatio
AcidicToBasicRatioDescriptor::AcidicToBasicRatioDescriptor()
    : ImageDescriptorBase("AcidicToBasicRatio", "The ratio of acidic centers to basic centers") {}

std::variant<double, int, std::string> AcidicToBasicRatioDescriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("AcidicToBasicRatio") ? descriptors["AcidicToBasicRatio"] : 0.0;
}

// 97. MeanAcidicAcidicDist
MeanAcidicAcidicDistDescriptor::MeanAcidicAcidicDistDescriptor()
    : ImageDescriptorBase("MeanAcidicAcidicDist", "The average distance between acidic centers") {}

std::variant<double, int, std::string> MeanAcidicAcidicDistDescriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("MeanAcidicAcidicDist") ? descriptors["MeanAcidicAcidicDist"] : 0.0;
}

// 98. MeanBasicBasicDist
MeanBasicBasicDistDescriptor::MeanBasicBasicDistDescriptor()
    : ImageDescriptorBase("MeanBasicBasicDist", "The average distance between basic centers") {}

std::variant<double, int, std::string> MeanBasicBasicDistDescriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("MeanBasicBasicDist") ? descriptors["MeanBasicBasicDist"] : 0.0;
}

// 99. MeanAcidicBasicDist
MeanAcidicBasicDistDescriptor::MeanAcidicBasicDistDescriptor()
    : ImageDescriptorBase("MeanAcidicBasicDist", "The average distance between acidic and basic centers") {}

std::variant<double, int, std::string> MeanAcidicBasicDistDescriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("MeanAcidicBasicDist") ? descriptors["MeanAcidicBasicDist"] : 0.0;
}

// 100. MinAcidicBasicDist
MinAcidicBasicDistDescriptor::MinAcidicBasicDistDescriptor()
    : ImageDescriptorBase("MinAcidicBasicDist", "The minimum distance between acidic and basic centers") {}

std::variant<double, int, std::string> MinAcidicBasicDistDescriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("MinAcidicBasicDist") ? descriptors["MinAcidicBasicDist"] : 0.0;
}

// 101. MeanAtomsNearAcidic
MeanAtomsNearAcidicDescriptor::MeanAtomsNearAcidicDescriptor()
    : ImageDescriptorBase("MeanAtomsNearAcidic", "The average number of atoms near acidic centers") {}

std::variant<double, int, std::string> MeanAtomsNearAcidicDescriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("MeanAtomsNearAcidic") ? descriptors["MeanAtomsNearAcidic"] : 0.0;
}

// 102. MeanAtomsNearBasic
MeanAtomsNearBasicDescriptor::MeanAtomsNearBasicDescriptor()
    : ImageDescriptorBase("MeanAtomsNearBasic", "The average number of atoms near basic centers") {}

std::variant<double, int, std::string> MeanAtomsNearBasicDescriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("MeanAtomsNearBasic") ? descriptors["MeanAtomsNearBasic"] : 0.0;
}

// 103. MeanAcidicRadius
MeanAcidicRadiusDescriptor::MeanAcidicRadiusDescriptor()
    : ImageDescriptorBase("MeanAcidicRadius", "The average radius of acidic centers") {}

std::variant<double, int, std::string> MeanAcidicRadiusDescriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("MeanAcidicRadius") ? descriptors["MeanAcidicRadius"] : 0.0;
}

// 104. MeanBasicRadius
MeanBasicRadiusDescriptor::MeanBasicRadiusDescriptor()
    : ImageDescriptorBase("MeanBasicRadius", "Mean radius of atoms with basic properties") {
}

std::variant<double, int, std::string> MeanBasicRadiusDescriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors["MeanBasicRadius"];
}

// 105. MeanAcidicCentroidDist
MeanAcidicCentroidDistDescriptor::MeanAcidicCentroidDistDescriptor()
    : ImageDescriptorBase("MeanAcidicCentroidDist", "The average distance of acidic centers to the centroid of the molecule") {}

std::variant<double, int, std::string> MeanAcidicCentroidDistDescriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("MeanAcidicCentroidDist") ? descriptors["MeanAcidicCentroidDist"] : 0.0;
}

// 106. MeanBasicCentroidDist
MeanBasicCentroidDistDescriptor::MeanBasicCentroidDistDescriptor()
    : ImageDescriptorBase("MeanBasicCentroidDist", "The average distance of basic centers to the centroid of the molecule") {}

std::variant<double, int, std::string> MeanBasicCentroidDistDescriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("MeanBasicCentroidDist") ? descriptors["MeanBasicCentroidDist"] : 0.0;
}

// 107. FracAcidicOnHull
FracAcidicOnHullDescriptor::FracAcidicOnHullDescriptor()
    : ImageDescriptorBase("FracAcidicOnHull", "The fraction of acidic centers on the convex hull") {}

std::variant<double, int, std::string> FracAcidicOnHullDescriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("FracAcidicOnHull") ? descriptors["FracAcidicOnHull"] : 0.0;
}

// 108. FracBasicOnHull
FracBasicOnHullDescriptor::FracBasicOnHullDescriptor()
    : ImageDescriptorBase("FracBasicOnHull", "The fraction of basic centers on the convex hull") {}

std::variant<double, int, std::string> FracBasicOnHullDescriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("FracBasicOnHull") ? descriptors["FracBasicOnHull"] : 0.0;
}

// 109. MeanAcidicAngle
MeanAcidicAngleDescriptor::MeanAcidicAngleDescriptor()
    : ImageDescriptorBase("MeanAcidicAngle", "The average angle of bonds to acidic centers") {}

std::variant<double, int, std::string> MeanAcidicAngleDescriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("MeanAcidicAngle") ? descriptors["MeanAcidicAngle"] : 0.0;
}

// 110. MeanBasicAngle
MeanBasicAngleDescriptor::MeanBasicAngleDescriptor()
    : ImageDescriptorBase("MeanBasicAngle", "The average angle of bonds to basic centers") {}

std::variant<double, int, std::string> MeanBasicAngleDescriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("MeanBasicAngle") ? descriptors["MeanBasicAngle"] : 0.0;
}

// 111. MeanAcidicLuminance
MeanAcidicLuminanceDescriptor::MeanAcidicLuminanceDescriptor()
    : ImageDescriptorBase("MeanAcidicLuminance", "The average luminance of acidic atoms") {}

std::variant<double, int, std::string> MeanAcidicLuminanceDescriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("MeanAcidicLuminance") ? descriptors["MeanAcidicLuminance"] : 0.0;
}

// 112. MeanBasicLuminance
MeanBasicLuminanceDescriptor::MeanBasicLuminanceDescriptor()
    : ImageDescriptorBase("MeanBasicLuminance", "The average luminance of basic atoms") {}

std::variant<double, int, std::string> MeanBasicLuminanceDescriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("MeanBasicLuminance") ? descriptors["MeanBasicLuminance"] : 0.0;
}

// 113. AcidicBasicLuminanceDiff
AcidicBasicLuminanceDiffDescriptor::AcidicBasicLuminanceDiffDescriptor()
    : ImageDescriptorBase("AcidicBasicLuminanceDiff", "The difference in average luminance between acidic and basic atoms") {}

std::variant<double, int, std::string> AcidicBasicLuminanceDiffDescriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("AcidicBasicLuminanceDiff") ? descriptors["AcidicBasicLuminanceDiff"] : 0.0;
}

// 114. FractalDimension
FractalDimensionDescriptor::FractalDimensionDescriptor()
    : ImageDescriptorBase("FractalDimension", "The fractal dimension of the molecule") {}

std::variant<double, int, std::string> FractalDimensionDescriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("FractalDimension") ? descriptors["FractalDimension"] : 0.0;
}

// 115. ForegroundRatio
ForegroundRatioDescriptor::ForegroundRatioDescriptor()
    : ImageDescriptorBase("ForegroundRatio", "The ratio of the area occupied by atoms to the total area of the molecule") {}

std::variant<double, int, std::string> ForegroundRatioDescriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("ForegroundRatio") ? descriptors["ForegroundRatio"] : 0.0;
}

// 116. AverageColor
AverageColorDescriptor::AverageColorDescriptor()
    : ImageDescriptorBase("AverageColor", "The average color of the molecule") {}

std::variant<double, int, std::string> AverageColorDescriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("AverageColor") ? descriptors["AverageColor"] : 0.0;
}

// 117. ColorVariance
ColorVarianceDescriptor::ColorVarianceDescriptor()
    : ImageDescriptorBase("ColorVariance", "The variance of colors in the molecule") {}

std::variant<double, int, std::string> ColorVarianceDescriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("ColorVariance") ? descriptors["ColorVariance"] : 0.0;
}

// 118. ImageCenterX
ImageCenterXDescriptor::ImageCenterXDescriptor()
    : ImageDescriptorBase("ImageCenterX", "The x-coordinate of the image center") {}

std::variant<double, int, std::string> ImageCenterXDescriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("ImageCenterX") ? descriptors["ImageCenterX"] : 0.0;
}

// 119. ImageCenterY
ImageCenterYDescriptor::ImageCenterYDescriptor()
    : ImageDescriptorBase("ImageCenterY", "The y-coordinate of the image center") {}

std::variant<double, int, std::string> ImageCenterYDescriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("ImageCenterY") ? descriptors["ImageCenterY"] : 0.0;
}

// 120. ImageOrientation
ImageOrientationDescriptor::ImageOrientationDescriptor()
    : ImageDescriptorBase("ImageOrientation", "The orientation of the molecule") {}

std::variant<double, int, std::string> ImageOrientationDescriptor::calculate(const Molecule& mol) const {
    auto descriptors = computeImageDescriptors(mol);
    return descriptors.count("ImageOrientation") ? descriptors["ImageOrientation"] : 0.0;
}

} // namespace descriptors
} // namespace desfact
