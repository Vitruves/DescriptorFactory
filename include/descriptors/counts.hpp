#pragma once

#include "descriptors.hpp"
#include <string>
#include <variant>
#include <GraphMol/ROMol.h>

namespace RDKit {
    class ROMol;
}

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

class PhosphorusCount : public Descriptor {
public:
    PhosphorusCount();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

class SiliconCount : public Descriptor {
public:
    SiliconCount();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

class BoronCount : public Descriptor {
public:
    BoronCount();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

class MetalAtomCount : public Descriptor {
public:
    MetalAtomCount();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

class NonMetalAtomCount : public Descriptor {
public:
    NonMetalAtomCount();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

class NumSpAtoms : public Descriptor {
public:
    NumSpAtoms();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

class NumSp2Atoms : public Descriptor {
public:
    NumSp2Atoms();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

class NumSp3Atoms : public Descriptor {
public:
    NumSp3Atoms();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

class RingCount3 : public Descriptor {
public:
    RingCount3();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

class RingCount4 : public Descriptor {
public:
    RingCount4();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

class RingCount5 : public Descriptor {
public:
    RingCount5();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

class RingCount6 : public Descriptor {
public:
    RingCount6();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

class RingCountLarge : public Descriptor {
public:
    RingCountLarge();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

class FusedBondCount : public Descriptor {
public:
    FusedBondCount();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

class SpiroAtomCount : public Descriptor {
public:
    SpiroAtomCount();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

class MacrocycleRingCount : public Descriptor {
public:
    MacrocycleRingCount();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

class CNBondCount : public Descriptor {
public:
    CNBondCount();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

class COBondCount : public Descriptor {
public:
    COBondCount();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

class CSBondCount : public Descriptor {
public:
    CSBondCount();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

class CHaloBondCount : public Descriptor {
public:
    CHaloBondCount();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
};

class AlcoholGroupCount : public Descriptor {
public:
    AlcoholGroupCount();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
private:
    static RDKit::ROMol* alcoholSmarts;
};

class EtherGroupCount : public Descriptor {
public:
    EtherGroupCount();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
private:
    static RDKit::ROMol* etherSmarts;
};

class AmineCount : public Descriptor {
public:
    AmineCount();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
private:
    static RDKit::ROMol* amineSmarts;
};

class EsterGroupCount : public Descriptor {
public:
    EsterGroupCount();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
private:
    static RDKit::ROMol* esterSmarts;
};

class KetoneGroupCount : public Descriptor {
public:
    KetoneGroupCount();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
private:
    static RDKit::ROMol* ketoneSmarts;
};

class AldehydeGroupCount : public Descriptor {
public:
    AldehydeGroupCount();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
private:
    static RDKit::ROMol* aldehydeSmarts;
};

class NitroGroupCount : public Descriptor {
public:
    NitroGroupCount();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
private:
    static RDKit::ROMol* nitroSmarts;
};

class SulfonylGroupCount : public Descriptor {
public:
    SulfonylGroupCount();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
private:
    static RDKit::ROMol* sulfonylSmarts;
};

class CyanoGroupCount : public Descriptor {
public:
    CyanoGroupCount();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
private:
    static RDKit::ROMol* cyanoSmarts;
};

class PhenylGroupCount : public Descriptor {
public:
    PhenylGroupCount();
    std::variant<double, int, std::string> calculate(const ::desfact::Molecule& mol) const override;
private:
    static RDKit::ROMol* phenylSmarts;
};

} // namespace descriptors
} // namespace desfact