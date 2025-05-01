// ... existing code ...
#pragma once

#include "descriptors.hpp"
#include <string>
#include <variant>

namespace desfact {
namespace descriptors {

class HalogenCount : public Descriptor {
public:
    HalogenCount();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

class CarbonCount : public Descriptor {
public:
    CarbonCount();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

class HydrogenCount : public Descriptor {
public:
    HydrogenCount();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

class OxygenCount : public Descriptor {
public:
    OxygenCount();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

class NitrogenCount : public Descriptor {
public:
    NitrogenCount();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

class SulfurCount : public Descriptor {
public:
    SulfurCount();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

class AromaticAtomCount : public Descriptor {
public:
    AromaticAtomCount();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

class NonAromaticAtomCount : public Descriptor {
public:
    NonAromaticAtomCount();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

class DoubleBondCount : public Descriptor {
public:
    DoubleBondCount();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

class TripleBondCount : public Descriptor {
public:
    TripleBondCount();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

class SingleBondCount : public Descriptor {
public:
    SingleBondCount();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

class AcidicFunctionCount : public Descriptor {
public:
    AcidicFunctionCount();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

class BasicFunctionCount : public Descriptor {
public:
    BasicFunctionCount();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

class ENAtoms2BondsFromAcidic : public Descriptor {
public:
    ENAtoms2BondsFromAcidic();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

class ENAtoms2BondsFromBasic : public Descriptor {
public:
    ENAtoms2BondsFromBasic();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

class ENAtoms3BondsFromAcidic : public Descriptor {
public:
    ENAtoms3BondsFromAcidic();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

class ENAtoms3BondsFromBasic : public Descriptor {
public:
    ENAtoms3BondsFromBasic();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

class LongestCSequenceSmiles : public Descriptor {
public:
    LongestCSequenceSmiles();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

class UppercaseCountSmiles : public Descriptor {
public:
    UppercaseCountSmiles();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

class LowercaseCountSmiles : public Descriptor {
public:
    LowercaseCountSmiles();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

class AtomVolumeSum : public Descriptor {
public:
    AtomVolumeSum();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

class HeavyAtomCount : public Descriptor {
public:
    HeavyAtomCount();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

class HeteroatomCount : public Descriptor {
public:
    HeteroatomCount();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

class RingAtomCount : public Descriptor {
public:
    RingAtomCount();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

class RingCount : public Descriptor {
public:
    RingCount();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

class ChiralCenterCount : public Descriptor {
public:
    ChiralCenterCount();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

class FormalChargeCount : public Descriptor {
public:
    FormalChargeCount();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

class PositiveChargeCount : public Descriptor {
public:
    PositiveChargeCount();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

class NegativeChargeCount : public Descriptor {
public:
    NegativeChargeCount();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

class RotatableBondCount : public Descriptor {
public:
    RotatableBondCount();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

class BridgeheadAtomCount : public Descriptor {
public:
    BridgeheadAtomCount();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

} // namespace descriptors
} // namespace desfact