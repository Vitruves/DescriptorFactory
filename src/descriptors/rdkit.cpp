#include "descriptors/rdkit.hpp"
#include "utils.hpp"
#include <GraphMol/GraphMol.h>
#include <GraphMol/Descriptors/MolDescriptors.h>
#include <GraphMol/Descriptors/Lipinski.h>

namespace desfact {
namespace descriptors {

RDKitDescriptor::RDKitDescriptor(const std::string& name, const std::string& description,
                                DescriptorFunction calcFunction)
    : Descriptor(name, description), calcFunction_(std::move(calcFunction)) {}

std::variant<double, int, std::string> RDKitDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) {
        globalLogger.debug(getName() + ": Invalid molecule input.");
        return 0.0;
    }

    const RDKit::ROMol* rdkMol = mol.getMolecule().get();
    if (!rdkMol) {
        globalLogger.error(getName() + ": RDKit molecule is null despite Molecule being valid?");
        return std::string("Error: Internal RDKit Null");
    }

    try {
        double result = calcFunction_(*rdkMol);
        return result;
    } catch (const std::exception& e) {
        globalLogger.error(getName() + ": RDKit calculation error: " + e.what());
        return std::string("Error: " + std::string(e.what()));
    } catch (...) {
        globalLogger.error(getName() + ": Unknown RDKit calculation error");
        return std::string("Error: Unknown RDKit error");
    }
}

// Molecular Weight
RDKitMWDescriptor::RDKitMWDescriptor()
    : RDKitDescriptor("rdkitMW", "Molecular weight using RDKit",
                     [](const RDKit::ROMol& mol) {
                         return RDKit::Descriptors::calcExactMW(mol);
                     }) {}

// TPSA
RDKitTPSADescriptor::RDKitTPSADescriptor()
    : RDKitDescriptor("rdkitTPSA", "Topological Polar Surface Area using RDKit",
                     [](const RDKit::ROMol& mol) {
                         return RDKit::Descriptors::calcTPSA(mol);
                     }) {}

// LogP
RDKitLogPDescriptor::RDKitLogPDescriptor()
    : RDKitDescriptor("rdkitLogP", "LogP using RDKit's Crippen algorithm",
                     [](const RDKit::ROMol& mol) {
                         double logp, mr;
                         RDKit::Descriptors::calcCrippenDescriptors(mol, logp, mr);
                         return logp;
                     }) {}

// MR
RDKitMRDescriptor::RDKitMRDescriptor()
    : RDKitDescriptor("rdkitMR", "Molar Refractivity using RDKit's Crippen algorithm",
                     [](const RDKit::ROMol& mol) {
                         double logp, mr;
                         RDKit::Descriptors::calcCrippenDescriptors(mol, logp, mr);
                         return mr;
                     }) {}

// Lipinski HBA
RDKitLipinskiHBADescriptor::RDKitLipinskiHBADescriptor()
    : RDKitDescriptor("rdkitLipinskiHBA", "Number of Lipinski hydrogen bond acceptors using RDKit",
                     [](const RDKit::ROMol& mol) {
                         return static_cast<double>(RDKit::Descriptors::calcLipinskiHBA(mol));
                     }) {}

// Lipinski HBD
RDKitLipinskiHBDDescriptor::RDKitLipinskiHBDDescriptor()
    : RDKitDescriptor("rdkitLipinskiHBD", "Number of Lipinski hydrogen bond donors using RDKit",
                     [](const RDKit::ROMol& mol) {
                         return static_cast<double>(RDKit::Descriptors::calcLipinskiHBD(mol));
                     }) {}

// Number of Rotatable Bonds
RDKitNumRotatableBondsDescriptor::RDKitNumRotatableBondsDescriptor()
    : RDKitDescriptor("rdkitNumRotatableBonds", "Number of rotatable bonds using RDKit",
                     [](const RDKit::ROMol& mol) {
                         return static_cast<double>(RDKit::Descriptors::calcNumRotatableBonds(mol));
                     }) {}

// Number of H-Bond Acceptors
RDKitNumHBondAcceptorsDescriptor::RDKitNumHBondAcceptorsDescriptor()
    : RDKitDescriptor("rdkitNumHBondAcceptors", "Number of H-bond acceptors using RDKit",
                     [](const RDKit::ROMol& mol) {
                         return static_cast<double>(RDKit::Descriptors::calcNumHBA(mol));
                     }) {}

// Number of H-Bond Donors
RDKitNumHBondDonorsDescriptor::RDKitNumHBondDonorsDescriptor()
    : RDKitDescriptor("rdkitNumHBondDonors", "Number of H-bond donors using RDKit",
                     [](const RDKit::ROMol& mol) {
                         return static_cast<double>(RDKit::Descriptors::calcNumHBD(mol));
                     }) {}

// Number of Heavy Atoms
RDKitNumHeavyAtomsDescriptor::RDKitNumHeavyAtomsDescriptor()
    : RDKitDescriptor("rdkitNumHeavyAtoms", "Number of heavy atoms using RDKit",
                     [](const RDKit::ROMol& mol) {
                         return static_cast<double>(mol.getNumHeavyAtoms());
                     }) {}

// Number of Rings
RDKitNumRingsDescriptor::RDKitNumRingsDescriptor()
    : RDKitDescriptor("rdkitNumRings", "Number of rings using RDKit",
                     [](const RDKit::ROMol& mol) {
                         return static_cast<double>(RDKit::Descriptors::calcNumRings(mol));
                     }) {}

// Number of Aromatic Rings
RDKitNumAromaticRingsDescriptor::RDKitNumAromaticRingsDescriptor()
    : RDKitDescriptor("rdkitNumAromaticRings", "Number of aromatic rings using RDKit",
                     [](const RDKit::ROMol& mol) {
                         return static_cast<double>(RDKit::Descriptors::calcNumAromaticRings(mol));
                     }) {}

// Number of Aliphatic Rings
RDKitNumAliphaticRingsDescriptor::RDKitNumAliphaticRingsDescriptor()
    : RDKitDescriptor("rdkitNumAliphaticRings", "Number of aliphatic rings using RDKit",
                     [](const RDKit::ROMol& mol) {
                         return static_cast<double>(RDKit::Descriptors::calcNumAliphaticRings(mol));
                     }) {}

// Fraction of sp3 carbons
RDKitFractionCsp3Descriptor::RDKitFractionCsp3Descriptor()
    : RDKitDescriptor("rdkitFractionCsp3", "Fraction of sp3 hybridized carbons using RDKit",
                     [](const RDKit::ROMol& mol) {
                         return RDKit::Descriptors::calcFractionCSP3(mol);
                     }) {}

// Number of chiral centers
RDKitChiralCentersDescriptor::RDKitChiralCentersDescriptor()
    : RDKitDescriptor("rdkitChiralCenters", "Number of chiral centers using RDKit",
                     [](const RDKit::ROMol& mol) {
                         return static_cast<double>(RDKit::Descriptors::numAtomStereoCenters(mol));
                     }) {}

// Labute ASA
RDKitLabuteASADescriptor::RDKitLabuteASADescriptor()
    : RDKitDescriptor("rdkitLabuteASA", "Labute ASA using RDKit",
                     [](const RDKit::ROMol& mol) {
                         return RDKit::Descriptors::calcLabuteASA(mol);
                     }) {}

// Hall-Kier Alpha
RDKitHallKierAlphaDescriptor::RDKitHallKierAlphaDescriptor()
    : RDKitDescriptor("rdkitHallKierAlpha", "Hall-Kier alpha using RDKit",
                     [](const RDKit::ROMol& mol) {
                         return RDKit::Descriptors::calcHallKierAlpha(mol);
                     }) {}

// Kappa1
RDKitKappaDescriptor1::RDKitKappaDescriptor1()
    : RDKitDescriptor("rdkitKappa1", "Kappa 1 index using RDKit",
                     [](const RDKit::ROMol& mol) {
                         return RDKit::Descriptors::calcKappa1(mol);
                     }) {}

// Kappa2
RDKitKappaDescriptor2::RDKitKappaDescriptor2()
    : RDKitDescriptor("rdkitKappa2", "Kappa 2 index using RDKit",
                     [](const RDKit::ROMol& mol) {
                         return RDKit::Descriptors::calcKappa2(mol);
                     }) {}

// Kappa3
RDKitKappaDescriptor3::RDKitKappaDescriptor3()
    : RDKitDescriptor("rdkitKappa3", "Kappa 3 index using RDKit",
                     [](const RDKit::ROMol& mol) {
                         return RDKit::Descriptors::calcKappa3(mol);
                     }) {}

// Chi0 descriptor
RDKitChi0vDescriptor::RDKitChi0vDescriptor()
    : RDKitDescriptor("rdkitChi0v", "Chi0v connectivity index using RDKit",
                     [](const RDKit::ROMol& mol) {
                         return RDKit::Descriptors::calcChi0v(mol);
                     }) {}

// Chi1 descriptor
RDKitChi1vDescriptor::RDKitChi1vDescriptor()
    : RDKitDescriptor("rdkitChi1v", "Chi1v connectivity index using RDKit",
                     [](const RDKit::ROMol& mol) {
                         return RDKit::Descriptors::calcChi1v(mol);
                     }) {}

// Chi2 descriptor
RDKitChi2vDescriptor::RDKitChi2vDescriptor()
    : RDKitDescriptor("rdkitChi2v", "Chi2v connectivity index using RDKit",
                     [](const RDKit::ROMol& mol) {
                         return RDKit::Descriptors::calcChi2v(mol);
                     }) {}

// Chi3 descriptor
RDKitChi3vDescriptor::RDKitChi3vDescriptor()
    : RDKitDescriptor("rdkitChi3v", "Chi3v connectivity index using RDKit",
                     [](const RDKit::ROMol& mol) {
                         return RDKit::Descriptors::calcChi3v(mol);
                     }) {}

// Chi4 descriptor
RDKitChi4vDescriptor::RDKitChi4vDescriptor()
    : RDKitDescriptor("rdkitChi4v", "Chi4v connectivity index using RDKit",
                     [](const RDKit::ROMol& mol) {
                         return RDKit::Descriptors::calcChi4v(mol);
                     }) {}

// MURCKO Framework Count
RDKitMurckoScaffoldCountDescriptor::RDKitMurckoScaffoldCountDescriptor()
    : RDKitDescriptor("rdkitMurckoCount", "Number of Murcko Scaffolds using RDKit",
                     [](const RDKit::ROMol& mol) {
                         // Placeholder: actual Murcko scaffold count may require more logic
                         return 0.0;
                     }) {}

// Ring Count by Size - 5-membered rings
RDKitRingCount5Descriptor::RDKitRingCount5Descriptor()
    : RDKitDescriptor("rdkitRingCount5", "Number of 5-membered rings using RDKit",
                     [](const RDKit::ROMol& mol) {
                         // Placeholder: actual 5-membered ring count logic
                         return 0.0;
                     }) {}

// Spiro atom count
RDKitSpiroAtomCountDescriptor::RDKitSpiroAtomCountDescriptor()
    : RDKitDescriptor("rdkitSpiroAtomCount", "Number of spiro atoms using RDKit",
                     [](const RDKit::ROMol& mol) {
                         return RDKit::Descriptors::calcNumSpiroAtoms(mol);
                     }) {}

// Bridgehead atom count
RDKitBridgeheadAtomCountDescriptor::RDKitBridgeheadAtomCountDescriptor()
    : RDKitDescriptor("rdkitBridgeheadAtomCount", "Number of bridgehead atoms using RDKit",
                     [](const RDKit::ROMol& mol) {
                         return RDKit::Descriptors::calcNumBridgeheadAtoms(mol);
                     }) {}

} // namespace descriptors
} // namespace desfact
