#pragma once

#include "descriptors.hpp"
#include <string>
#include <variant>

namespace RDKit { // Forward declare
class ROMol;
class Atom;
} // namespace RDKit

namespace desfact {
namespace descriptors {

// Helper functions potentially needed (can be defined in .cpp)
// - calculate_gasteiger_charges (maybe cached?)
// - calculate_vdw_radius
// - calculate_statistics (mean, stddev, skew, kurtosis, range)

// --- Category 1: Global Descriptors based on Statistics of (Radius & Charge) ---

class SumRadiusChargeDescriptor : public Descriptor {
public:
    SumRadiusChargeDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class AvgRadiusChargeDescriptor : public Descriptor {
public:
    AvgRadiusChargeDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class StdDevRadiusChargeDescriptor : public Descriptor {
public:
    StdDevRadiusChargeDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class RangeRadiusChargeDescriptor : public Descriptor {
public:
    RangeRadiusChargeDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class SkewRadiusChargeDescriptor : public Descriptor {
public:
    SkewRadiusChargeDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class KurtosisRadiusChargeDescriptor : public Descriptor {
public:
    KurtosisRadiusChargeDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class SumRadiusPerAbsChargeDescriptor : public Descriptor {
public:
    SumRadiusPerAbsChargeDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class AvgRadiusPositiveChargeDescriptor : public Descriptor {
public:
    AvgRadiusPositiveChargeDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class AvgRadiusNegativeChargeDescriptor : public Descriptor {
public:
    AvgRadiusNegativeChargeDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class RatioSumRadiusChargedDescriptor : public Descriptor {
public:
    RatioSumRadiusChargedDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class WeightedAvgElectronegativityDescriptor : public Descriptor {
public:
    WeightedAvgElectronegativityDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class WeightedStdDevElectronegativityDescriptor : public Descriptor {
public:
    WeightedStdDevElectronegativityDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class SumRadiusSqAbsChargeDescriptor : public Descriptor {
public:
    SumRadiusSqAbsChargeDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class WeightedStdDevChargeDescriptor : public Descriptor {
public:
    WeightedStdDevChargeDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class TotalAbsChargePerRadiusDescriptor : public Descriptor {
public:
    TotalAbsChargePerRadiusDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class SumInvRadiusAbsChargeDescriptor : public Descriptor {
public:
    SumInvRadiusAbsChargeDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class SumEnRadiusSqDescriptor : public Descriptor {
public:
    SumEnRadiusSqDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class RatioAvgEnRadiusRingChainDescriptor : public Descriptor {
public:
    RatioAvgEnRadiusRingChainDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class AvgRatioChargeRadiusDescriptor : public Descriptor {
public:
    AvgRatioChargeRadiusDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class RatioAvgRChargeToAvgRPlusChargeDescriptor : public Descriptor {
public:
    RatioAvgRChargeToAvgRPlusChargeDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class AvgRadiusLonePairAtomsDescriptor : public Descriptor {
public:
    AvgRadiusLonePairAtomsDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class AvgRadiusPiSystemAtomsDescriptor : public Descriptor {
public:
    AvgRadiusPiSystemAtomsDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class StdDevRatioChargeRadiusDescriptor : public Descriptor {
public:
    StdDevRatioChargeRadiusDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class SumRadiusPositiveChargeThreshDescriptor : public Descriptor {
public:
    SumRadiusPositiveChargeThreshDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class SumRadiusNegativeChargeThreshDescriptor : public Descriptor {
public:
    SumRadiusNegativeChargeThreshDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

// Topological Descriptors (Many placeholders)
class WeightedPathCountDescriptor : public Descriptor {
public:
    WeightedPathCountDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override; // Placeholder
};

class AvgWeightedPathLengthDescriptor : public Descriptor {
public:
    AvgWeightedPathLengthDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override; // Placeholder
};

class BalabanLikeIndexRChargeDescriptor : public Descriptor {
public:
    BalabanLikeIndexRChargeDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override; // Placeholder
};

class WienerLikeIndexRChargeDescriptor : public Descriptor {
public:
    WienerLikeIndexRChargeDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override; // Placeholder
};

class EccentricityMaxRadiusAtomChargeWeightedDescriptor : public Descriptor {
public:
    EccentricityMaxRadiusAtomChargeWeightedDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class TopoDistMaxRMaxPosChargeDescriptor : public Descriptor {
public:
    TopoDistMaxRMaxPosChargeDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class TopoDistMaxRMinNegChargeDescriptor : public Descriptor {
public:
    TopoDistMaxRMinNegChargeDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class AvgNeighborWeightRqDescriptor : public Descriptor {
public:
    AvgNeighborWeightRqDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class TopoAutocorrRqDist2Descriptor : public Descriptor {
public:
    TopoAutocorrRqDist2Descriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class TopoAutocorrRqDist3Descriptor : public Descriptor {
public:
    TopoAutocorrRqDist3Descriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class TopoAutocorrRqDist4Descriptor : public Descriptor {
public:
    TopoAutocorrRqDist4Descriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class AvgTopoDistPiWeightedRqDescriptor : public Descriptor {
public:
    AvgTopoDistPiWeightedRqDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class AvgTopoDistLPWeightedRqDescriptor : public Descriptor {
public:
    AvgTopoDistLPWeightedRqDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class AvgTopoDistPiLpWeightedRqDescriptor : public Descriptor {
public:
    AvgTopoDistPiLpWeightedRqDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class PathCountAlternatingRChargeDescriptor : public Descriptor {
public:
    PathCountAlternatingRChargeDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override; // Placeholder
};

class RatioLongestWeightedPathRingChainDescriptor : public Descriptor {
public:
    RatioLongestWeightedPathRingChainDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override; // Placeholder
};

class SumWeightRChargeDegreeGt3Descriptor : public Descriptor {
public:
    SumWeightRChargeDegreeGt3Descriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class RatioAvgRWeightedChargeTerminalDescriptor : public Descriptor {
public:
    RatioAvgRWeightedChargeTerminalDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class EigenvalueWeightedConnectivityDescriptor : public Descriptor {
public:
    EigenvalueWeightedConnectivityDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override; // Placeholder (NaN)
};

class WeightedBranchPointComplexityDescriptor : public Descriptor {
public:
    WeightedBranchPointComplexityDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class ShannonEntropyRChargeDescriptor : public Descriptor {
public:
    ShannonEntropyRChargeDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class AvgRadiusPlusChargePerDegreeDescriptor : public Descriptor {
public:
    AvgRadiusPlusChargePerDegreeDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class SumDeltaRadiusChargeWeightedDescriptor : public Descriptor {
public:
    SumDeltaRadiusChargeWeightedDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class WeightedBuriedAtomCountDescriptor : public Descriptor {
public:
    WeightedBuriedAtomCountDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override; // Placeholder
};

class RatioSumRqFormalChargeDescriptor : public Descriptor {
public:
    RatioSumRqFormalChargeDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class AvgRadiusHBDonorWeightedChargeDescriptor : public Descriptor {
public:
    AvgRadiusHBDonorWeightedChargeDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class AvgRadiusHBAcceptorWeightedChargeDescriptor : public Descriptor {
public:
    AvgRadiusHBAcceptorWeightedChargeDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class RatioSumRPolarNonpolarFragDescriptor : public Descriptor {
public:
    RatioSumRPolarNonpolarFragDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override; // Placeholder
};

class SumRadiusSp2CWeightedChargeDescriptor : public Descriptor {
public:
    SumRadiusSp2CWeightedChargeDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class SumRadiusSp3CWeightedChargeDescriptor : public Descriptor {
public:
    SumRadiusSp3CWeightedChargeDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class AvgNeighborChargeRadiusPerDegreeDescriptor : public Descriptor {
public:
    AvgNeighborChargeRadiusPerDegreeDescriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override;
};

class KierShapeIndexVariant3Descriptor : public Descriptor {
public:
    KierShapeIndexVariant3Descriptor();
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override; // Placeholder
};

} // namespace descriptors
} // namespace desfact
