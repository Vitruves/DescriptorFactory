#pragma once

#include "descriptors.hpp"

namespace desfact {
namespace descriptors {

class ImageDescriptorBase : public Descriptor {
public:
    ImageDescriptorBase(const std::string& name, const std::string& description);
    
protected:
    // Helper function to generate 2D coords and compute image descriptors
    static std::unordered_map<std::string, double> computeImageDescriptors(const Molecule& mol);
};

// Now defining all individual image descriptor classes

class BoundingBoxAreaDescriptor : public ImageDescriptorBase {
public:
    BoundingBoxAreaDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class MoleculeWidthDescriptor : public ImageDescriptorBase {
public:
    MoleculeWidthDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class MoleculeHeightDescriptor : public ImageDescriptorBase {
public:
    MoleculeHeightDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class AspectRatioDescriptor : public ImageDescriptorBase {
public:
    AspectRatioDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class AtomDensityDescriptor : public ImageDescriptorBase {
public:
    AtomDensityDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class MeanAtomRadiusDescriptor : public ImageDescriptorBase {
public:
    MeanAtomRadiusDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class MaxAtomRadiusDescriptor : public ImageDescriptorBase {
public:
    MaxAtomRadiusDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class MinAtomRadiusDescriptor : public ImageDescriptorBase {
public:
    MinAtomRadiusDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class MeanBondLengthDescriptor : public ImageDescriptorBase {
public:
    MeanBondLengthDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class MaxBondLengthDescriptor : public ImageDescriptorBase {
public:
    MaxBondLengthDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class MinBondLengthDescriptor : public ImageDescriptorBase {
public:
    MinBondLengthDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class BondLengthStdDescriptor : public ImageDescriptorBase {
public:
    BondLengthStdDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class MeanAtomAtomDistDescriptor : public ImageDescriptorBase {
public:
    MeanAtomAtomDistDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class MeanAtomLuminanceDescriptor : public ImageDescriptorBase {
public:
    MeanAtomLuminanceDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class MaxAtomAtomDistDescriptor : public ImageDescriptorBase {
public:
    MaxAtomAtomDistDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class MinAtomAtomDistDescriptor : public ImageDescriptorBase {
public:
    MinAtomAtomDistDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class MeanAtomAreaDescriptor : public ImageDescriptorBase {
public:
    MeanAtomAreaDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class MedianAtomRadiusDescriptor : public ImageDescriptorBase {
public:
    MedianAtomRadiusDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class MedianBondLengthDescriptor : public ImageDescriptorBase {
public:
    MedianBondLengthDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class AtomRadiusRangeDescriptor : public ImageDescriptorBase {
public:
    AtomRadiusRangeDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class AtomAreaFractionDescriptor : public ImageDescriptorBase {
public:
    AtomAreaFractionDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class BondCoverageFractionDescriptor : public ImageDescriptorBase {
public:
    BondCoverageFractionDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class AtomPackingDensityDescriptor : public ImageDescriptorBase {
public:
    AtomPackingDensityDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class AtomMassXDescriptor : public ImageDescriptorBase {
public:
    AtomMassXDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class AtomMassYDescriptor : public ImageDescriptorBase {
public:
    AtomMassYDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class AtomMassDistDescriptor : public ImageDescriptorBase {
public:
    AtomMassDistDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class AtomRadiusStdDescriptor : public ImageDescriptorBase {
public:
    AtomRadiusStdDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class MeanBondAngleDescriptor : public ImageDescriptorBase {
public:
    MeanBondAngleDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class BondAngleStdDescriptor : public ImageDescriptorBase {
public:
    BondAngleStdDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class AtomXStdDescriptor : public ImageDescriptorBase {
public:
    AtomXStdDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class AtomYStdDescriptor : public ImageDescriptorBase {
public:
    AtomYStdDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class AtomLuminanceStdDescriptor : public ImageDescriptorBase {
public:
    AtomLuminanceStdDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class BondLenMADDescriptor : public ImageDescriptorBase {
public:
    BondLenMADDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class AtomRadiusMADDescriptor : public ImageDescriptorBase {
public:
    AtomRadiusMADDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class AtomAreaStdDescriptor : public ImageDescriptorBase {
public:
    AtomAreaStdDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class AtomRadiusCVDescriptor : public ImageDescriptorBase {
public:
    AtomRadiusCVDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class AtomAreaCVDescriptor : public ImageDescriptorBase {
public:
    AtomAreaCVDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class AtomLuminanceCVDescriptor : public ImageDescriptorBase {
public:
    AtomLuminanceCVDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class AtomAreaRangeDescriptor : public ImageDescriptorBase {
public:
    AtomAreaRangeDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class AtomAreaMedianDescriptor : public ImageDescriptorBase {
public:
    AtomAreaMedianDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class AtomAreaMADDescriptor : public ImageDescriptorBase {
public:
    AtomAreaMADDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class AtomAreaHullFracDescriptor : public ImageDescriptorBase {
public:
    AtomAreaHullFracDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class AtomAreaCenterFracDescriptor : public ImageDescriptorBase {
public:
    AtomAreaCenterFracDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

// Additional descriptors

class StdAllAtomDistDescriptor : public ImageDescriptorBase {
public:
    StdAllAtomDistDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class MinBondAngleDescriptor : public ImageDescriptorBase {
public:
    MinBondAngleDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class MaxBondAngleDescriptor : public ImageDescriptorBase {
public:
    MaxBondAngleDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class MedianBondAngleDescriptor : public ImageDescriptorBase {
public:
    MedianBondAngleDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class MeanBondColorDiffDescriptor : public ImageDescriptorBase {
public:
    MeanBondColorDiffDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class MeanBondLuminanceDiffDescriptor : public ImageDescriptorBase {
public:
    MeanBondLuminanceDiffDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class MeanBondLenDiffColorDescriptor : public ImageDescriptorBase {
public:
    MeanBondLenDiffColorDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class MeanBondAngleDiffColorDescriptor : public ImageDescriptorBase {
public:
    MeanBondAngleDiffColorDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class MeanDistSameColorDescriptor : public ImageDescriptorBase {
public:
    MeanDistSameColorDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class MeanDistDiffColorDescriptor : public ImageDescriptorBase {
public:
    MeanDistDiffColorDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

// Bond order descriptors
class MeanBondOrderDescriptor : public ImageDescriptorBase {
public:
    MeanBondOrderDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class FracDoubleBondsDescriptor : public ImageDescriptorBase {
public:
    FracDoubleBondsDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class FracTripleBondsDescriptor : public ImageDescriptorBase {
public:
    FracTripleBondsDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class FracAromaticBondsDescriptor : public ImageDescriptorBase {
public:
    FracAromaticBondsDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

// Atom degree descriptors
class MeanAtomDegreeDescriptor : public ImageDescriptorBase {
public:
    MeanAtomDegreeDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class MaxAtomDegreeDescriptor : public ImageDescriptorBase {
public:
    MaxAtomDegreeDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class MinAtomDegreeDescriptor : public ImageDescriptorBase {
public:
    MinAtomDegreeDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

// Centroid descriptors
class MeanDistToCentroidDescriptor : public ImageDescriptorBase {
public:
    MeanDistToCentroidDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class StdDistToCentroidDescriptor : public ImageDescriptorBase {
public:
    StdDistToCentroidDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

// Bond length by type descriptors
class MeanLenSingleDescriptor : public ImageDescriptorBase {
public:
    MeanLenSingleDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class MeanLenDoubleDescriptor : public ImageDescriptorBase {
public:
    MeanLenDoubleDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class MeanLenTripleDescriptor : public ImageDescriptorBase {
public:
    MeanLenTripleDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class MeanLenAromaticDescriptor : public ImageDescriptorBase {
public:
    MeanLenAromaticDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

// Atom degree fraction descriptors
class FracDeg1Descriptor : public ImageDescriptorBase {
public:
    FracDeg1Descriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class FracDeg2Descriptor : public ImageDescriptorBase {
public:
    FracDeg2Descriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class FracDeg3Descriptor : public ImageDescriptorBase {
public:
    FracDeg3Descriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class FracDeg4Descriptor : public ImageDescriptorBase {
public:
    FracDeg4Descriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

// Bond property descriptors
class MeanBondRadiiDiffDescriptor : public ImageDescriptorBase {
public:
    MeanBondRadiiDiffDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

// Additional bond property descriptors
class MeanBondLuminanceDiff2Descriptor : public ImageDescriptorBase {
public:
    MeanBondLuminanceDiff2Descriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class MeanAngleHighDegreeDescriptor : public ImageDescriptorBase {
public:
    MeanAngleHighDegreeDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class BoundaryAtomRatioDescriptor : public ImageDescriptorBase {
public:
    BoundaryAtomRatioDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class MeanAtomEccentricityDescriptor : public ImageDescriptorBase {
public:
    MeanAtomEccentricityDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class MolecularDiameterDescriptor : public ImageDescriptorBase {
public:
    MolecularDiameterDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

// Molecular shape descriptors
class MolecularRadiusDescriptor : public ImageDescriptorBase {
public:
    MolecularRadiusDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class RadiusOfGyrationDescriptor : public ImageDescriptorBase {
public:
    RadiusOfGyrationDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class MolecularSphericityDescriptor : public ImageDescriptorBase {
public:
    MolecularSphericityDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class BondLenToAtomRadiusRatioDescriptor : public ImageDescriptorBase {
public:
    BondLenToAtomRadiusRatioDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class EdgeDensityDescriptor : public ImageDescriptorBase {
public:
    EdgeDensityDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

// Additional shape and topology descriptors
class PlanarityMeasureDescriptor : public ImageDescriptorBase {
public:
    PlanarityMeasureDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class MeanNearestNeighborDistDescriptor : public ImageDescriptorBase {
public:
    MeanNearestNeighborDistDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class MeanAtomsInRadiusDescriptor : public ImageDescriptorBase {
public:
    MeanAtomsInRadiusDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class VarNearestNeighborDistDescriptor : public ImageDescriptorBase {
public:
    VarNearestNeighborDistDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class MeanBondToBondAngleDescriptor : public ImageDescriptorBase {
public:
    MeanBondToBondAngleDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class MomentOfInertiaDescriptor : public ImageDescriptorBase {
public:
    MomentOfInertiaDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class MeanAtomCentralityDescriptor : public ImageDescriptorBase {
public:
    MeanAtomCentralityDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

// Graph topology descriptors
class TerminalAtomCountDescriptor : public ImageDescriptorBase {
public:
    TerminalAtomCountDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class JunctionAtomCountDescriptor : public ImageDescriptorBase {
public:
    JunctionAtomCountDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class TerminalToJunctionRatioDescriptor : public ImageDescriptorBase {
public:
    TerminalToJunctionRatioDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class WienerIndexDescriptor : public ImageDescriptorBase {
public:
    WienerIndexDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class RingCountDescriptor : public ImageDescriptorBase {
public:
    RingCountDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

// Chemical property descriptors
class AcidicCenterCountDescriptor : public ImageDescriptorBase {
public:
    AcidicCenterCountDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class BasicCenterCountDescriptor : public ImageDescriptorBase {
public:
    BasicCenterCountDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class AcidicToBasicRatioDescriptor : public ImageDescriptorBase {
public:
    AcidicToBasicRatioDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

// Chemical interaction descriptors
class MeanAcidicAcidicDistDescriptor : public ImageDescriptorBase {
public:
    MeanAcidicAcidicDistDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class MeanBasicBasicDistDescriptor : public ImageDescriptorBase {
public:
    MeanBasicBasicDistDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class MeanAcidicBasicDistDescriptor : public ImageDescriptorBase {
public:
    MeanAcidicBasicDistDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class MinAcidicBasicDistDescriptor : public ImageDescriptorBase {
public:
    MinAcidicBasicDistDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class MeanAtomsNearAcidicDescriptor : public ImageDescriptorBase {
public:
    MeanAtomsNearAcidicDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class MeanAtomsNearBasicDescriptor : public ImageDescriptorBase {
public:
    MeanAtomsNearBasicDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class MeanAcidicRadiusDescriptor : public ImageDescriptorBase {
public:
    MeanAcidicRadiusDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

// Additional chemical descriptors
class MeanBasicRadiusDescriptor : public ImageDescriptorBase {
public:
    MeanBasicRadiusDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class MeanAcidicCentroidDistDescriptor : public ImageDescriptorBase {
public:
    MeanAcidicCentroidDistDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class MeanBasicCentroidDistDescriptor : public ImageDescriptorBase {
public:
    MeanBasicCentroidDistDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class FracAcidicOnHullDescriptor : public ImageDescriptorBase {
public:
    FracAcidicOnHullDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class FracBasicOnHullDescriptor : public ImageDescriptorBase {
public:
    FracBasicOnHullDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class MeanAcidicAngleDescriptor : public ImageDescriptorBase {
public:
    MeanAcidicAngleDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

// Additional chemical descriptors
class MeanBasicAngleDescriptor : public ImageDescriptorBase {
public:
    MeanBasicAngleDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class MeanAcidicLuminanceDescriptor : public ImageDescriptorBase {
public:
    MeanAcidicLuminanceDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class MeanBasicLuminanceDescriptor : public ImageDescriptorBase {
public:
    MeanBasicLuminanceDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class AcidicBasicLuminanceDiffDescriptor : public ImageDescriptorBase {
public:
    AcidicBasicLuminanceDiffDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

// Statistical and geometric descriptors
class FractalDimensionDescriptor : public ImageDescriptorBase {
public:
    FractalDimensionDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class ForegroundRatioDescriptor : public ImageDescriptorBase {
public:
    ForegroundRatioDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

// Additional color and geometry descriptors
class AverageColorDescriptor : public ImageDescriptorBase {
public:
    AverageColorDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class ColorVarianceDescriptor : public ImageDescriptorBase {
public:
    ColorVarianceDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class ImageCenterXDescriptor : public ImageDescriptorBase {
public:
    ImageCenterXDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class ImageCenterYDescriptor : public ImageDescriptorBase {
public:
    ImageCenterYDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class ImageOrientationDescriptor : public ImageDescriptorBase {
public:
    ImageOrientationDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

} // namespace descriptors
} // namespace desfact
