#include "descriptors.hpp"
#include "descriptors/fractional.hpp"
#include "descriptors/sum.hpp"
#include "descriptors/eigen.hpp"
#include "descriptors/vague3.hpp"
#include "descriptors/strings.hpp"
#include "descriptors/vague4.hpp"
#include "descriptors/electronic.hpp"
#include "descriptors/vague5.hpp"
#include "descriptors/vague6.hpp"
#include "descriptors/vague7.hpp"
#include "descriptors/morgan.hpp"
#include "descriptors/selfies.hpp"
#include "descriptors/pka.hpp"
#include "descriptors/vague8.hpp"
#include "descriptors/solubility.hpp"
#include "descriptors/counts.hpp"
#include "descriptors/image.hpp"
#include <GraphMol/Descriptors/MolDescriptors.h>
#include <GraphMol/GraphMol.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/Fingerprints/MorganFingerprints.h>
#include <GraphMol/Fingerprints/Fingerprints.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <algorithm>
#include "utils.hpp"
#ifdef WITH_TBB
#include <tbb/parallel_for.h>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/concurrent_vector.h>
#include <tbb/task_arena.h>
#include <tbb/global_control.h>
#include <tbb/blocked_range.h>
#endif

namespace desfact {

// Descriptor implementations
std::variant<double, int, std::string> MolecularWeightDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) {
        return 0.0;
    }
    
    try {
        double mw = RDKit::Descriptors::calcExactMW(*mol.getMolecule());
        return mw;
    } catch (const std::exception& e) {
        globalLogger.error("MW Calc Error: " + std::string(e.what()) + " for SMILES: " + mol.getOriginalSmiles());
        return "Error: MW Calc";
    } catch (...) {
        globalLogger.error("Unknown MW Calc Error for SMILES: " + mol.getOriginalSmiles());
        return "Error: MW Calc Unknown";
    }
}

std::variant<double, int, std::string> NumAtomsDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid() || !mol.getMolecule()) {
        return -1;
    }
    
    try {
        return static_cast<int>(mol.getMolecule()->getNumAtoms());
    } catch (const std::exception& e) {
        globalLogger.error("NumAtoms Calc Error: " + std::string(e.what()) + " for SMILES: " + mol.getOriginalSmiles());
        return "Error: NumAtoms Calc";
    } catch (...) {
        globalLogger.error("Unknown NumAtoms Calc Error for SMILES: " + mol.getOriginalSmiles());
        return "Error: NumAtoms Calc Unknown";
    }
}

// DescriptorFactory implementation
DescriptorFactory::DescriptorFactory() {
    globalLogger.info("Initializing descriptor factory...");
    // Register core descriptors
    registerDescriptor(std::make_unique<MolecularWeightDescriptor>());
    registerDescriptor(std::make_unique<NumAtomsDescriptor>());
    
    // Register fractional descriptors - Atomic Fraction Descriptors
    registerDescriptor(std::make_unique<descriptors::FcCDescriptor>());
    registerDescriptor(std::make_unique<descriptors::FcFDescriptor>());
    registerDescriptor(std::make_unique<descriptors::FcODescriptor>());
    registerDescriptor(std::make_unique<descriptors::FcClDescriptor>());
    registerDescriptor(std::make_unique<descriptors::FcBrDescriptor>());
    registerDescriptor(std::make_unique<descriptors::FcIDescriptor>());
    registerDescriptor(std::make_unique<descriptors::FcSDescriptor>());
    registerDescriptor(std::make_unique<descriptors::FcNDescriptor>());
    registerDescriptor(std::make_unique<descriptors::FcHDescriptor>());
    registerDescriptor(std::make_unique<descriptors::FcHaloDescriptor>());
    registerDescriptor(std::make_unique<descriptors::FcHeteroDescriptor>());
    registerDescriptor(std::make_unique<descriptors::FcPolarDescriptor>());
    registerDescriptor(std::make_unique<descriptors::FcApolarDescriptor>());
    registerDescriptor(std::make_unique<descriptors::FcMetalsDescriptor>());
    registerDescriptor(std::make_unique<descriptors::FcLowEADescriptor>());
    registerDescriptor(std::make_unique<descriptors::FcIEOddDescriptor>());
    registerDescriptor(std::make_unique<descriptors::FcHeteroLowPolzDescriptor>());
    registerDescriptor(std::make_unique<descriptors::FcEvenValenceAtomsDescriptor>());
    registerDescriptor(std::make_unique<descriptors::FcNonHeteroHighEADescriptor>());
    
    // Hybridization-Based Bond Fractions
    registerDescriptor(std::make_unique<descriptors::FcCSp3Descriptor>());
    registerDescriptor(std::make_unique<descriptors::FcCSp2Descriptor>());
    
    // Bond Polarity Descriptors
    registerDescriptor(std::make_unique<descriptors::FcUnpolDescriptor>());
    registerDescriptor(std::make_unique<descriptors::FcPolDescriptor>());
    
    // Electronegativity-Based Descriptors
    registerDescriptor(std::make_unique<descriptors::FcSumPolAtDescriptor>());
    registerDescriptor(std::make_unique<descriptors::FcSumPolMWDescriptor>());
    
    // Bond Valence Descriptors
    registerDescriptor(std::make_unique<descriptors::FcBondAtDescriptor>());
    registerDescriptor(std::make_unique<descriptors::FcBondNDescriptor>());
    registerDescriptor(std::make_unique<descriptors::FcBondODescriptor>());
    registerDescriptor(std::make_unique<descriptors::FcBondCDescriptor>());
    registerDescriptor(std::make_unique<descriptors::FcBondSDescriptor>());
    registerDescriptor(std::make_unique<descriptors::FcBondPDescriptor>());
    
    // Additional normalized descriptors
    registerDescriptor(std::make_unique<descriptors::FcHBDonorsDescriptor>());
    registerDescriptor(std::make_unique<descriptors::FcHBAcceptorsDescriptor>());
    registerDescriptor(std::make_unique<descriptors::FcAromaticAtomsDescriptor>());
    registerDescriptor(std::make_unique<descriptors::FcRingAtomsDescriptor>());
    registerDescriptor(std::make_unique<descriptors::FcBridgeAtomsDescriptor>());
    registerDescriptor(std::make_unique<descriptors::FcChargedAtomsDescriptor>());
    registerDescriptor(std::make_unique<descriptors::FcHeavyAtomsDescriptor>());
    
    // Electronegativity-Based Fractions
    registerDescriptor(std::make_unique<descriptors::FcENAboveAvgDescriptor>());
    registerDescriptor(std::make_unique<descriptors::FcENBelowAvgDescriptor>());
    registerDescriptor(std::make_unique<descriptors::FcENHighDescriptor>());
    registerDescriptor(std::make_unique<descriptors::FcENLowDescriptor>());
    
    // Atomic Radius
    registerDescriptor(std::make_unique<descriptors::FcSmallRDescriptor>());
    registerDescriptor(std::make_unique<descriptors::FcLargeRDescriptor>());
    
    // Polarizability
    registerDescriptor(std::make_unique<descriptors::FcLowPolzDescriptor>());
    registerDescriptor(std::make_unique<descriptors::FcHighPolzDescriptor>());
    
    // Ionization Energy
    // registerDescriptor(std::make_unique<descriptors::FcLowIEDescriptor>()); // Zero or near-zero variance
    
    // Electron Affinity
    registerDescriptor(std::make_unique<descriptors::FcHighEADescriptor>());
    
    // Van der Waals Volume
    registerDescriptor(std::make_unique<descriptors::FcSmallVdWDescriptor>());
    registerDescriptor(std::make_unique<descriptors::FcLargeVdWDescriptor>());
    
    // Electronegativity Contribution
    // registerDescriptor(std::make_unique<descriptors::FcENMWDescriptor>()); // Zero or near-zero variance
    registerDescriptor(std::make_unique<descriptors::FcENBondedDescriptor>());
    
    // Radius-Based
    registerDescriptor(std::make_unique<descriptors::FcVdWMWDescriptor>());
    registerDescriptor(std::make_unique<descriptors::FcRcovMWDescriptor>());
    
    // Valence / Oxidation State
    registerDescriptor(std::make_unique<descriptors::FcHighOxStateDescriptor>());
    
    // Metal / Metalloid / Nonmetal Classification
    registerDescriptor(std::make_unique<descriptors::FcMetalloidDescriptor>());
    
    // Combined descriptors
    registerDescriptor(std::make_unique<descriptors::FcHETpolDescriptor>());
    registerDescriptor(std::make_unique<descriptors::FcHALpolDescriptor>());
    registerDescriptor(std::make_unique<descriptors::FcHeavyPolDescriptor>());
    
    // Register Sum descriptors
    
    // Atomic property sums
    registerDescriptor(std::make_unique<descriptors::SumENDescriptor>());
    registerDescriptor(std::make_unique<descriptors::SumCovalentRadiiDescriptor>());
    registerDescriptor(std::make_unique<descriptors::SumVdWVolumeDescriptor>());
    registerDescriptor(std::make_unique<descriptors::SumPolarizabilityDescriptor>());
    registerDescriptor(std::make_unique<descriptors::SumIonizationEnergyDescriptor>());
    registerDescriptor(std::make_unique<descriptors::SumElectronAffinityDescriptor>());
    
    // Atom Class Sums
    registerDescriptor(std::make_unique<descriptors::SumAtomsDescriptor>());
    registerDescriptor(std::make_unique<descriptors::SumHeavyAtomsDescriptor>());
    registerDescriptor(std::make_unique<descriptors::SumHeteroatomsDescriptor>());
    registerDescriptor(std::make_unique<descriptors::SumHalogensDescriptor>());
    registerDescriptor(std::make_unique<descriptors::SumChargedAtomsDescriptor>());
    
    // Bond and Hybridization Sums
    registerDescriptor(std::make_unique<descriptors::SumBondsDescriptor>());
    registerDescriptor(std::make_unique<descriptors::SumDoubleBondsDescriptor>());
    registerDescriptor(std::make_unique<descriptors::SumTripleBondsDescriptor>());
    registerDescriptor(std::make_unique<descriptors::SumAromaticBondsDescriptor>());
    registerDescriptor(std::make_unique<descriptors::SumPolarBondsDescriptor>());
    registerDescriptor(std::make_unique<descriptors::SumUnpolarBondsDescriptor>());
    registerDescriptor(std::make_unique<descriptors::SumSp3BondsDescriptor>());
    registerDescriptor(std::make_unique<descriptors::SumSp2BondsDescriptor>());
    
    // Structural Sums
    registerDescriptor(std::make_unique<descriptors::SumRingsDescriptor>());
    registerDescriptor(std::make_unique<descriptors::SumAromaticRingsDescriptor>());
    registerDescriptor(std::make_unique<descriptors::SumRotatableBondsDescriptor>());
    registerDescriptor(std::make_unique<descriptors::SumBridgeAtomsDescriptor>());
    registerDescriptor(std::make_unique<descriptors::SumRingAtomsDescriptor>());
    
    // Physicochemical Contribution Sums
    registerDescriptor(std::make_unique<descriptors::SumENMWDescriptor>());
    registerDescriptor(std::make_unique<descriptors::SumPolMWDescriptor>());
    registerDescriptor(std::make_unique<descriptors::SumENRcovDescriptor>());
    registerDescriptor(std::make_unique<descriptors::SumENPolDescriptor>());
    
    // Additional Hybridization Sums
    registerDescriptor(std::make_unique<descriptors::SumSp3AtomsDescriptor>());
    registerDescriptor(std::make_unique<descriptors::SumSp2AtomsDescriptor>());
    
    // Additional Pharmacophore-related Sums
    registerDescriptor(std::make_unique<descriptors::SumHBDonorsDescriptor>());
    registerDescriptor(std::make_unique<descriptors::SumHBAcceptorsDescriptor>());
    
    // Additional Functional Group Sums
    registerDescriptor(std::make_unique<descriptors::SumAmideGroupsDescriptor>());
    registerDescriptor(std::make_unique<descriptors::SumCarboxylGroupsDescriptor>());
    
    // Additional Sum descriptors
    registerDescriptor(std::make_unique<descriptors::SumElectronegativityRingDescriptor>());
    registerDescriptor(std::make_unique<descriptors::SumElectronegativityBondedDescriptor>());
    registerDescriptor(std::make_unique<descriptors::SumAtomicRadiusHeavyDescriptor>());
    registerDescriptor(std::make_unique<descriptors::SumPolzENRatioDescriptor>());
    registerDescriptor(std::make_unique<descriptors::SumIEENRatioDescriptor>());
    registerDescriptor(std::make_unique<descriptors::SumEAMWENDescriptor>());
    registerDescriptor(std::make_unique<descriptors::SumSp2CAtomsDescriptor>());
    registerDescriptor(std::make_unique<descriptors::SumSpAtomsENDescriptor>());
    registerDescriptor(std::make_unique<descriptors::SumFormalChargeAbsDescriptor>());
    registerDescriptor(std::make_unique<descriptors::SumFormalChargeHeavyDescriptor>());
    registerDescriptor(std::make_unique<descriptors::SumFormalChargePositiveDescriptor>());
    registerDescriptor(std::make_unique<descriptors::SumFormalChargeNegativeDescriptor>());
    registerDescriptor(std::make_unique<descriptors::SumENMWRatioDescriptor>());
    registerDescriptor(std::make_unique<descriptors::SumPolzRingDescriptor>());
    registerDescriptor(std::make_unique<descriptors::SumIEBondedDescriptor>());
    registerDescriptor(std::make_unique<descriptors::SumEAHeavyBondedDescriptor>());
    registerDescriptor(std::make_unique<descriptors::SumSp3ENWeightedDescriptor>());
    registerDescriptor(std::make_unique<descriptors::SumHeteroMWENDescriptor>());
    registerDescriptor(std::make_unique<descriptors::SumENSp2HeavyDescriptor>());
    registerDescriptor(std::make_unique<descriptors::SumENRadicalAtomsDescriptor>());
    registerDescriptor(std::make_unique<descriptors::SumOxHeteroMWDescriptor>());
    registerDescriptor(std::make_unique<descriptors::SumENRingHeavyDescriptor>());
    registerDescriptor(std::make_unique<descriptors::SumPolENBondedDescriptor>());
    
    // New Fractional Descriptors
    // registerDescriptor(std::make_unique<descriptors::FcRadiusMWRatioAbove1Descriptor>()); // Zero or near-zero variance
    // registerDescriptor(std::make_unique<descriptors::FcOxStateOddDescriptor>()); // Zero or near-zero variance
    // registerDescriptor(std::make_unique<descriptors::FcLowRadicalENDescriptor>()); // Zero or near-zero variance
    // registerDescriptor(std::make_unique<descriptors::FcGroup17OxAbove1Descriptor>()); // Zero or near-zero variance
    registerDescriptor(std::make_unique<descriptors::FcSp3HeavyAtomsDescriptor>());
    registerDescriptor(std::make_unique<descriptors::FcSp2ENAboveAvgDescriptor>());
    registerDescriptor(std::make_unique<descriptors::FcOxENAboveThresholdDescriptor>());
    registerDescriptor(std::make_unique<descriptors::FcFormalChargeNonZeroDescriptor>());
    registerDescriptor(std::make_unique<descriptors::FcFormalChargePositiveDescriptor>());
    registerDescriptor(std::make_unique<descriptors::FcFormalChargeNegativeDescriptor>());
    registerDescriptor(std::make_unique<descriptors::FcENMWAbove2Descriptor>());
    registerDescriptor(std::make_unique<descriptors::FcGroup16AtomsDescriptor>());
    registerDescriptor(std::make_unique<descriptors::FcSpAtomsHighENDescriptor>());
    registerDescriptor(std::make_unique<descriptors::FcRingOxHighDescriptor>());
    registerDescriptor(std::make_unique<descriptors::FcHeavyFormalChargeDescriptor>());
    registerDescriptor(std::make_unique<descriptors::FcSp3PolarizabilityAbove10Descriptor>());
    registerDescriptor(std::make_unique<descriptors::FcHeavyOxNegativeDescriptor>());
    registerDescriptor(std::make_unique<descriptors::FcENAboveMoleculeAvgDescriptor>());
    
    // Register String-based descriptors
    
    // ASCII and Bit-Level Descriptors
    registerDescriptor(std::make_unique<descriptors::AsciiSumDescriptor>());
    registerDescriptor(std::make_unique<descriptors::AsciiAverageDescriptor>());
    registerDescriptor(std::make_unique<descriptors::AsciiVarianceDescriptor>());
    registerDescriptor(std::make_unique<descriptors::AsciiRangeDescriptor>());
    registerDescriptor(std::make_unique<descriptors::BitCountDescriptor>());
    registerDescriptor(std::make_unique<descriptors::BitDensityDescriptor>());
    // registerDescriptor(std::make_unique<descriptors::BitEntropyDescriptor>()); // Zero or near-zero variance
    registerDescriptor(std::make_unique<descriptors::BitTransitionRateDescriptor>());
    registerDescriptor(std::make_unique<descriptors::OddBitRatioDescriptor>());
    
    // Compression & Entropy Descriptors
    registerDescriptor(std::make_unique<descriptors::RunLengthEncodingSizeDescriptor>());
    registerDescriptor(std::make_unique<descriptors::CharacterRepetitionScoreDescriptor>());
    registerDescriptor(std::make_unique<descriptors::NibblePatternCountDescriptor>());
    registerDescriptor(std::make_unique<descriptors::BytePatternCountDescriptor>());
    registerDescriptor(std::make_unique<descriptors::EntropyPerByteDescriptor>());
    registerDescriptor(std::make_unique<descriptors::AsciiHistogramUniformityDescriptor>());
    registerDescriptor(std::make_unique<descriptors::FunctionalGroupEntropyDescriptor>());
    registerDescriptor(std::make_unique<descriptors::LargestRingDescriptor>());
    registerDescriptor(std::make_unique<descriptors::CarbonSkeletonComplexityDescriptor>());
    registerDescriptor(std::make_unique<descriptors::HeterocycleDescriptor>());
    registerDescriptor(std::make_unique<descriptors::SubstructureBalanceDescriptor>());
    
    // Character Sequence Descriptors
    registerDescriptor(std::make_unique<descriptors::LongestCharacterRunDescriptor>());
    registerDescriptor(std::make_unique<descriptors::SpecialCharacterDensityDescriptor>());
    registerDescriptor(std::make_unique<descriptors::CaseTransitionRateDescriptor>());
    registerDescriptor(std::make_unique<descriptors::AlternatingCharacterScoreDescriptor>());
    registerDescriptor(std::make_unique<descriptors::CharacterTrigramCountDescriptor>());
    registerDescriptor(std::make_unique<descriptors::LetterDigitRatioDescriptor>());
    registerDescriptor(std::make_unique<descriptors::CharacterBigramEntropyDescriptor>());
    
    // Bit Patterns Descriptors
    // registerDescriptor(std::make_unique<descriptors::ByteMirrorSymmetryDescriptor>()); // Zero or near-zero variance
    registerDescriptor(std::make_unique<descriptors::BytePalindromeScoreDescriptor>());
    registerDescriptor(std::make_unique<descriptors::BitPalindromeScoreDescriptor>());
    registerDescriptor(std::make_unique<descriptors::BitAutocorrelationDescriptor>());
    registerDescriptor(std::make_unique<descriptors::HammingWeightDistributionDescriptor>());
    registerDescriptor(std::make_unique<descriptors::ByteParityDistributionDescriptor>());
    registerDescriptor(std::make_unique<descriptors::BitRunEntropyDescriptor>());
    
    // Positional Heuristics Descriptors
    registerDescriptor(std::make_unique<descriptors::PositionWeightDescriptor>());
    registerDescriptor(std::make_unique<descriptors::FrontBackRatioDescriptor>());
    registerDescriptor(std::make_unique<descriptors::BracketAsymmetryDescriptor>());
    registerDescriptor(std::make_unique<descriptors::CharDispersionDescriptor>());
    registerDescriptor(std::make_unique<descriptors::DigitPositionMeanDescriptor>());
    registerDescriptor(std::make_unique<descriptors::SymmetryScoreDescriptor>());
    
    // Pattern / Fragment Pseudodescriptors
    registerDescriptor(std::make_unique<descriptors::AromaticRatioDescriptor>());
    registerDescriptor(std::make_unique<descriptors::HeteroatomRatioDescriptor>());
    registerDescriptor(std::make_unique<descriptors::ElementDiversityDescriptor>());
    registerDescriptor(std::make_unique<descriptors::SaturationIndexDescriptor>());
    registerDescriptor(std::make_unique<descriptors::PatternEntropyDescriptor>());
    registerDescriptor(std::make_unique<descriptors::BranchComplexityDescriptor>());
    
    // Basic SMILES Structure-Based Descriptors
    // registerDescriptor(std::make_unique<descriptors::RingPositionVarianceDescriptor>()); // Zero or near-zero variance
    registerDescriptor(std::make_unique<descriptors::BondPositionEntropyDescriptor>());
    registerDescriptor(std::make_unique<descriptors::ElementPositionBiasDescriptor>());
    registerDescriptor(std::make_unique<descriptors::BondComplexityDescriptor>());
    registerDescriptor(std::make_unique<descriptors::RingComplexityDescriptor>());
    registerDescriptor(std::make_unique<descriptors::ElementRunLengthVariance>());
    registerDescriptor(std::make_unique<descriptors::HeteroAtomSequenceLengthMax>());
    
    // Register Vague3 descriptors
    registerDescriptor(std::make_unique<descriptors::AtomCentralityVarianceDescriptor>());
    registerDescriptor(std::make_unique<descriptors::PeripheryCoreRatioDescriptor>());
    registerDescriptor(std::make_unique<descriptors::MinimumSpanningTreeDepthDescriptor>());
    registerDescriptor(std::make_unique<descriptors::TerminalBranchRatioDescriptor>());
    registerDescriptor(std::make_unique<descriptors::AverageAtomPathRedundancyDescriptor>());
    registerDescriptor(std::make_unique<descriptors::CarbonClusterDensityDescriptor>());
    registerDescriptor(std::make_unique<descriptors::HeteroatomClusterDensityDescriptor>());
    registerDescriptor(std::make_unique<descriptors::TerminalHeavyAtomDiversityDescriptor>());
    registerDescriptor(std::make_unique<descriptors::CentralHeteroatomRatioDescriptor>());
    registerDescriptor(std::make_unique<descriptors::ChainEndElementDiversityDescriptor>());
    registerDescriptor(std::make_unique<descriptors::HighENNeighborDensityDescriptor>());
    registerDescriptor(std::make_unique<descriptors::ElectronDeficientCarbonRatioDescriptor>());
    registerDescriptor(std::make_unique<descriptors::PolarizableAtomDispersionDescriptor>());
    registerDescriptor(std::make_unique<descriptors::DipoleMomentProxyDescriptor>());
    registerDescriptor(std::make_unique<descriptors::LongChainFragmentDensityDescriptor>());
    registerDescriptor(std::make_unique<descriptors::ShortChainDensityDescriptor>());
    registerDescriptor(std::make_unique<descriptors::SubstitutionDensityPerRingDescriptor>());
    registerDescriptor(std::make_unique<descriptors::LinearVsBranchedCarbonRatioDescriptor>());
    registerDescriptor(std::make_unique<descriptors::RingToBranchConnectivityRatioDescriptor>());
    registerDescriptor(std::make_unique<descriptors::IsolatedHBondSiteDensityDescriptor>());
    registerDescriptor(std::make_unique<descriptors::FunctionalGroupSeparationIndexDescriptor>());
    registerDescriptor(std::make_unique<descriptors::PeripheralHBondDensityDescriptor>());
    registerDescriptor(std::make_unique<descriptors::AromaticCoreRatioDescriptor>());
    registerDescriptor(std::make_unique<descriptors::ConjugationGapRatioDescriptor>());
    registerDescriptor(std::make_unique<descriptors::NonAromaticRingSubstitutionDescriptor>());
    registerDescriptor(std::make_unique<descriptors::AromaticToNonAromaticRingRatioDescriptor>());
    registerDescriptor(std::make_unique<descriptors::ConjugationTerminalDensityDescriptor>());
    registerDescriptor(std::make_unique<descriptors::ChargePairDensityDescriptor>());
    registerDescriptor(std::make_unique<descriptors::FormalChargeAccessibilityDescriptor>());
    registerDescriptor(std::make_unique<descriptors::ChargeSeparationIndexDescriptor>());
    registerDescriptor(std::make_unique<descriptors::NeutralAtomRatioDescriptor>());
    registerDescriptor(std::make_unique<descriptors::FunctionalGroupChargeBalanceDescriptor>());
    registerDescriptor(std::make_unique<descriptors::StereocenterDensityDescriptor>());
    registerDescriptor(std::make_unique<descriptors::DoubleBondConfigurationDensityDescriptor>());
    registerDescriptor(std::make_unique<descriptors::RingStereogenicAtomDensityDescriptor>());
    registerDescriptor(std::make_unique<descriptors::TerminalStereocenterRatioDescriptor>());
    registerDescriptor(std::make_unique<descriptors::AdjacentStereocenterDensityDescriptor>());
    registerDescriptor(std::make_unique<descriptors::GroupPeriodicDiversityDescriptor>());
    registerDescriptor(std::make_unique<descriptors::PeriodDiversityIndexDescriptor>());
    registerDescriptor(std::make_unique<descriptors::AtomicMassVarianceDescriptor>());
    registerDescriptor(std::make_unique<descriptors::RareElementDensityDescriptor>());
    registerDescriptor(std::make_unique<descriptors::AlkaliAlkalineEarthRatioDescriptor>());
    registerDescriptor(std::make_unique<descriptors::UnsaturationClustersDescriptor>());
    registerDescriptor(std::make_unique<descriptors::RingSpanningBondDensityDescriptor>());
    registerDescriptor(std::make_unique<descriptors::TerminalUnsaturatedBondRatioDescriptor>());
    registerDescriptor(std::make_unique<descriptors::HeavyAtomBondOrderVarianceDescriptor>());
    registerDescriptor(std::make_unique<descriptors::HeteroatomBondingDiversityDescriptor>());

    // Register Vague 4 Descriptors
    // 1. Structural and Connectivity-Based
    registerDescriptor(std::make_unique<descriptors::AtomicConnectivityImbalanceDescriptor>());
    // registerDescriptor(std::make_unique<descriptors::RingBridgeRatioDescriptor>()); // Near-zero variance
    registerDescriptor(std::make_unique<descriptors::SubstitutionPatternComplexityDescriptor>());
    registerDescriptor(std::make_unique<descriptors::OpenChainSaturationRatioDescriptor>());
    registerDescriptor(std::make_unique<descriptors::AtomValenceSpread>());
    registerDescriptor(std::make_unique<descriptors::MaxConsecutiveSingleBonds>());
    registerDescriptor(std::make_unique<descriptors::ValenceMismatchBondCount>());

    // 2. Atomic Neighborhood Patterns
    registerDescriptor(std::make_unique<descriptors::HeteroatomNeighborhoodDiversityDescriptor>());
    registerDescriptor(std::make_unique<descriptors::CarbonNeighborhoodUniformityDescriptor>());
    registerDescriptor(std::make_unique<descriptors::PolarAtomNeighborhoodRatioDescriptor>());
    registerDescriptor(std::make_unique<descriptors::AverageAtomDegreeRangeDescriptor>());
    registerDescriptor(std::make_unique<descriptors::AtomTypeNeighborDiversity>());

    // 3. Ring System Specific
    registerDescriptor(std::make_unique<descriptors::RingJunctionComplexityDescriptor>());
    registerDescriptor(std::make_unique<descriptors::NonFusedRingDensityDescriptor>());
    // registerDescriptor(std::make_unique<descriptors::RingBridgeRatioDescriptor>()); // Near-zero variance
    registerDescriptor(std::make_unique<descriptors::RingChainAttachmentDensityDescriptor>());
    registerDescriptor(std::make_unique<descriptors::RingTerminalSubstituentRatioDescriptor>());
    registerDescriptor(std::make_unique<descriptors::RingSaturationBalanceDescriptor>());
    registerDescriptor(std::make_unique<descriptors::RingHydrogenRatio>());
    registerDescriptor(std::make_unique<descriptors::SSSRtoRingAtomRatio>());
    registerDescriptor(std::make_unique<descriptors::MeanBondStereoStates>());
    registerDescriptor(std::make_unique<descriptors::RingSizeMedian>());
    registerDescriptor(std::make_unique<descriptors::RingIndexSum>());
    registerDescriptor(std::make_unique<descriptors::PercentageIsolatedRingSystems>());
    registerDescriptor(std::make_unique<descriptors::RingSystemCount>());
    registerDescriptor(std::make_unique<descriptors::MeanRingPerimeterDegree>());

    // 4. Electronic Influence
    registerDescriptor(std::make_unique<descriptors::PolarizabilityGradientDescriptor>());
    registerDescriptor(std::make_unique<descriptors::ElectronWithdrawingAtomDensityDescriptor>());
    registerDescriptor(std::make_unique<descriptors::ElectronDonatingAtomDensityDescriptor>());
    registerDescriptor(std::make_unique<descriptors::ElectronegativityGradientDensityDescriptor>());
    registerDescriptor(std::make_unique<descriptors::PeripheralElectronRichAtomRatioDescriptor>());
    registerDescriptor(std::make_unique<descriptors::MaxBondENGap>());

    // 5. Functional Group and Substitution Patterns
    registerDescriptor(std::make_unique<descriptors::FunctionalGroupIsolationIndexDescriptor>());
    registerDescriptor(std::make_unique<descriptors::HydroxylGroupDispersionDescriptor>());
    registerDescriptor(std::make_unique<descriptors::AlkylChainDiversityDescriptor>());
    registerDescriptor(std::make_unique<descriptors::SubstituentPositionVarianceDescriptor>());
    registerDescriptor(std::make_unique<descriptors::TerminalFunctionalGroupClusteringDescriptor>());
    registerDescriptor(std::make_unique<descriptors::AromaticSubstituentDiversity>());

    // 6. Bond-type and Hybridization
    registerDescriptor(std::make_unique<descriptors::ChainSaturationVariabilityDescriptor>());
    registerDescriptor(std::make_unique<descriptors::TripleBondTerminalRatioDescriptor>());
    registerDescriptor(std::make_unique<descriptors::AdjacentHybridizationTransitionDescriptor>());
    registerDescriptor(std::make_unique<descriptors::CyclicHybridizationHomogeneityDescriptor>());
    registerDescriptor(std::make_unique<descriptors::InternalChainUnsaturationDensityDescriptor>());
    registerDescriptor(std::make_unique<descriptors::MeanResonanceBondOrder>());
    //registerDescriptor(std::make_unique<descriptors::BondOrderEntropy>());
    registerDescriptor(std::make_unique<descriptors::HybridizationSymmetryIndex>());
    registerDescriptor(std::make_unique<descriptors::TerminalDoubleBondCount>());
    registerDescriptor(std::make_unique<descriptors::NonCarbonAtomAdjacencyCount>());
    registerDescriptor(std::make_unique<descriptors::ConsecutiveHeteroBondFraction>());

    // 7. Hydrogen Bonding Patterns
    registerDescriptor(std::make_unique<descriptors::HydrogenBondDonorClusteringDescriptor>());
    registerDescriptor(std::make_unique<descriptors::AcceptorDonorRatioImbalanceDescriptor>());
    registerDescriptor(std::make_unique<descriptors::PeripheralDonorAcceptorBalanceDescriptor>());
    registerDescriptor(std::make_unique<descriptors::IntraringHBondPotentialDescriptor>());
    registerDescriptor(std::make_unique<descriptors::ChainEndHydrogenBondDensityDescriptor>());

    // 8. Formal Charge Distribution
    registerDescriptor(std::make_unique<descriptors::FormalChargeNeighborhoodVarianceDescriptor>());
    registerDescriptor(std::make_unique<descriptors::OppositeChargeNeighborRatioDescriptor>());
    registerDescriptor(std::make_unique<descriptors::ChargeGradientAlongChainsDescriptor>());
    registerDescriptor(std::make_unique<descriptors::LocalizedChargeClustersDescriptor>());
    registerDescriptor(std::make_unique<descriptors::PeripheralChargeNeutralityDescriptor>());
    registerDescriptor(std::make_unique<descriptors::ChargeBalanceSkew>());

    // 9. Aromaticity & Conjugation Patterns
    registerDescriptor(std::make_unique<descriptors::InterRingConjugationRatioDescriptor>());
    registerDescriptor(std::make_unique<descriptors::AromaticChainDensityDescriptor>());
    registerDescriptor(std::make_unique<descriptors::ConjugationLengthVarianceDescriptor>());
    registerDescriptor(std::make_unique<descriptors::TerminalAromaticSubstitutionDescriptor>());
    registerDescriptor(std::make_unique<descriptors::CyclicVsChainAromaticRatioDescriptor>());

    // 10. Structural Diversity and Complexity
    registerDescriptor(std::make_unique<descriptors::UniqueElementPairRatioDescriptor>());
    registerDescriptor(std::make_unique<descriptors::HeavyAtomDegreeDiversityDescriptor>());
    registerDescriptor(std::make_unique<descriptors::AtomPathDiversityDescriptor>());
    registerDescriptor(std::make_unique<descriptors::InternalAtomComplexityDescriptor>());
    registerDescriptor(std::make_unique<descriptors::HeteroatomBondOrderVariabilityDescriptor>());
    registerDescriptor(std::make_unique<descriptors::CarbonIsotopeCount>());
    registerDescriptor(std::make_unique<descriptors::HeteroatomDegreeVariance>());
    registerDescriptor(std::make_unique<descriptors::HeavyAtomHydrogenRatio>());
    registerDescriptor(std::make_unique<descriptors::MeanAtomicMassPerDegree>());

    // Vague 6 Specific (from vague6.cpp)
    registerDescriptor(std::make_unique<descriptors::DegreeLogWeightedSum>());
    registerDescriptor(std::make_unique<descriptors::BondOrderSkewness>());
    registerDescriptor(std::make_unique<descriptors::PeripheralAtomTypeDiversity>());
    registerDescriptor(std::make_unique<descriptors::HeteroRingEdgeRatio>());
    registerDescriptor(std::make_unique<descriptors::BranchPointDensity>());
    registerDescriptor(std::make_unique<descriptors::PathLengthGini>());
    registerDescriptor(std::make_unique<descriptors::LongestHomoelementPath>());
    registerDescriptor(std::make_unique<descriptors::MeanBranchSeparation>());
    registerDescriptor(std::make_unique<descriptors::CarbonChainBranchingIndex>());
    registerDescriptor(std::make_unique<descriptors::BondTypeAlternationRatio>());
    registerDescriptor(std::make_unique<descriptors::HeteroatomPathFraction>());
    registerDescriptor(std::make_unique<descriptors::MaxRingDistance>());
    registerDescriptor(std::make_unique<descriptors::HalogenNeighborHybridRatio>());
    registerDescriptor(std::make_unique<descriptors::MeanAtomBetweennessCentrality>());
    registerDescriptor(std::make_unique<descriptors::LongestPathBondOrderProduct>());
    registerDescriptor(std::make_unique<descriptors::BranchDepthAverage>());
    registerDescriptor(std::make_unique<descriptors::TopoDistanceSkewness>());

    // Electronic Specific (from electronic.cpp)
    registerDescriptor(std::make_unique<descriptors::RatioAvgRChargeToAvgRPlusChargeDescriptor>());

    // Register Electronic Descriptors
    registerDescriptor(std::make_unique<descriptors::SumRadiusChargeDescriptor>());
    registerDescriptor(std::make_unique<descriptors::AvgRadiusChargeDescriptor>());
    registerDescriptor(std::make_unique<descriptors::StdDevRadiusChargeDescriptor>());
    registerDescriptor(std::make_unique<descriptors::RangeRadiusChargeDescriptor>());
    registerDescriptor(std::make_unique<descriptors::SkewRadiusChargeDescriptor>());
    registerDescriptor(std::make_unique<descriptors::KurtosisRadiusChargeDescriptor>());
    registerDescriptor(std::make_unique<descriptors::SumRadiusPerAbsChargeDescriptor>());
    registerDescriptor(std::make_unique<descriptors::AvgRadiusPositiveChargeDescriptor>());
    registerDescriptor(std::make_unique<descriptors::AvgRadiusNegativeChargeDescriptor>());
    registerDescriptor(std::make_unique<descriptors::RatioSumRadiusChargedDescriptor>());
    registerDescriptor(std::make_unique<descriptors::WeightedAvgElectronegativityDescriptor>());
    registerDescriptor(std::make_unique<descriptors::WeightedStdDevElectronegativityDescriptor>());
    registerDescriptor(std::make_unique<descriptors::SumRadiusSqAbsChargeDescriptor>());
    registerDescriptor(std::make_unique<descriptors::WeightedStdDevChargeDescriptor>());
    registerDescriptor(std::make_unique<descriptors::TotalAbsChargePerRadiusDescriptor>());
    registerDescriptor(std::make_unique<descriptors::SumInvRadiusAbsChargeDescriptor>());
    registerDescriptor(std::make_unique<descriptors::SumEnRadiusSqDescriptor>());
    registerDescriptor(std::make_unique<descriptors::RatioAvgEnRadiusRingChainDescriptor>());
    registerDescriptor(std::make_unique<descriptors::AvgRatioChargeRadiusDescriptor>());
    registerDescriptor(std::make_unique<descriptors::AvgRadiusLonePairAtomsDescriptor>());
    registerDescriptor(std::make_unique<descriptors::AvgRadiusPiSystemAtomsDescriptor>());
    registerDescriptor(std::make_unique<descriptors::StdDevRatioChargeRadiusDescriptor>());
    registerDescriptor(std::make_unique<descriptors::SumRadiusPositiveChargeThreshDescriptor>());
    registerDescriptor(std::make_unique<descriptors::SumRadiusNegativeChargeThreshDescriptor>());
    registerDescriptor(std::make_unique<descriptors::WeightedPathCountDescriptor>());
    registerDescriptor(std::make_unique<descriptors::AvgWeightedPathLengthDescriptor>());
    registerDescriptor(std::make_unique<descriptors::BalabanLikeIndexRChargeDescriptor>());
    registerDescriptor(std::make_unique<descriptors::WienerLikeIndexRChargeDescriptor>());
    registerDescriptor(std::make_unique<descriptors::EccentricityMaxRadiusAtomChargeWeightedDescriptor>());
    registerDescriptor(std::make_unique<descriptors::TopoDistMaxRMaxPosChargeDescriptor>());
    registerDescriptor(std::make_unique<descriptors::TopoDistMaxRMinNegChargeDescriptor>());
    registerDescriptor(std::make_unique<descriptors::AvgNeighborWeightRqDescriptor>());
    registerDescriptor(std::make_unique<descriptors::TopoAutocorrRqDist2Descriptor>());
    registerDescriptor(std::make_unique<descriptors::TopoAutocorrRqDist3Descriptor>());
    registerDescriptor(std::make_unique<descriptors::TopoAutocorrRqDist4Descriptor>());
    registerDescriptor(std::make_unique<descriptors::AvgTopoDistPiWeightedRqDescriptor>());
    registerDescriptor(std::make_unique<descriptors::AvgTopoDistLPWeightedRqDescriptor>());
    registerDescriptor(std::make_unique<descriptors::AvgTopoDistPiLpWeightedRqDescriptor>());
    registerDescriptor(std::make_unique<descriptors::PathCountAlternatingRChargeDescriptor>());
    registerDescriptor(std::make_unique<descriptors::RatioLongestWeightedPathRingChainDescriptor>());
    registerDescriptor(std::make_unique<descriptors::SumWeightRChargeDegreeGt3Descriptor>());
    registerDescriptor(std::make_unique<descriptors::RatioAvgRWeightedChargeTerminalDescriptor>());
    registerDescriptor(std::make_unique<descriptors::EigenvalueWeightedConnectivityDescriptor>());
    registerDescriptor(std::make_unique<descriptors::WeightedBranchPointComplexityDescriptor>());
    registerDescriptor(std::make_unique<descriptors::ShannonEntropyRChargeDescriptor>());
    registerDescriptor(std::make_unique<descriptors::AvgRadiusPlusChargePerDegreeDescriptor>());
    registerDescriptor(std::make_unique<descriptors::SumDeltaRadiusChargeWeightedDescriptor>());
    registerDescriptor(std::make_unique<descriptors::WeightedBuriedAtomCountDescriptor>());
    registerDescriptor(std::make_unique<descriptors::RatioSumRqFormalChargeDescriptor>());
    registerDescriptor(std::make_unique<descriptors::AvgRadiusHBDonorWeightedChargeDescriptor>());
    registerDescriptor(std::make_unique<descriptors::AvgRadiusHBAcceptorWeightedChargeDescriptor>());
    registerDescriptor(std::make_unique<descriptors::RatioSumRPolarNonpolarFragDescriptor>());
    registerDescriptor(std::make_unique<descriptors::SumRadiusSp2CWeightedChargeDescriptor>());
    registerDescriptor(std::make_unique<descriptors::SumRadiusSp3CWeightedChargeDescriptor>());
    registerDescriptor(std::make_unique<descriptors::AvgNeighborChargeRadiusPerDegreeDescriptor>());
    registerDescriptor(std::make_unique<descriptors::KierShapeIndexVariant3Descriptor>());

    // Register Eigen-based descriptors
    registerDescriptor(std::make_unique<descriptors::AdjacencyNonZeroEntries>());
    registerDescriptor(std::make_unique<descriptors::AdjacencyMatrixTrace>());
    registerDescriptor(std::make_unique<descriptors::AdjacencyFrobeniusNorm>());
    registerDescriptor(std::make_unique<descriptors::AdjacencySpectralRadius>());
    registerDescriptor(std::make_unique<descriptors::AdjacencySmallestEigenvalue>());
    registerDescriptor(std::make_unique<descriptors::AdjacencySumEigenvalues>());
    registerDescriptor(std::make_unique<descriptors::AdjacencySumSquaresEigenvalues>());
    registerDescriptor(std::make_unique<descriptors::AdjacencyEigenvalueVariance>());
    registerDescriptor(std::make_unique<descriptors::AdjacencyEigenvalueSkewness>());
    registerDescriptor(std::make_unique<descriptors::AdjacencyEigenvalueKurtosis>());
    registerDescriptor(std::make_unique<descriptors::AdjacencyPositiveEigenvalues>());
    registerDescriptor(std::make_unique<descriptors::AdjacencyNegativeEigenvalues>());
    registerDescriptor(std::make_unique<descriptors::AdjacencyZeroEigenvalues>());
    registerDescriptor(std::make_unique<descriptors::AdjacencyMaxDegree>());
    registerDescriptor(std::make_unique<descriptors::AdjacencyMinDegree>());
    registerDescriptor(std::make_unique<descriptors::AdjacencyMeanDegree>());
    registerDescriptor(std::make_unique<descriptors::AdjacencyDegreeVariance>());
    registerDescriptor(std::make_unique<descriptors::AdjacencyDegreeStdDev>());
    registerDescriptor(std::make_unique<descriptors::AdjacencyMatrixRank>());
    registerDescriptor(std::make_unique<descriptors::GraphEnergy>());
    registerDescriptor(std::make_unique<descriptors::LaplacianSpectralRadius>());
    registerDescriptor(std::make_unique<descriptors::LaplacianAlgebraicConnectivity>());
    registerDescriptor(std::make_unique<descriptors::LaplacianZeroEigenvalues>());
    registerDescriptor(std::make_unique<descriptors::LaplacianEnergy>());
    registerDescriptor(std::make_unique<descriptors::LaplacianMatrixTrace>());
    registerDescriptor(std::make_unique<descriptors::LaplacianMatrixDeterminant>());
    registerDescriptor(std::make_unique<descriptors::LaplacianTotalEffectiveResistance>());
    registerDescriptor(std::make_unique<descriptors::LaplacianKirchhoffIndex>());
    registerDescriptor(std::make_unique<descriptors::LaplacianEigenvalueVariance>());
    registerDescriptor(std::make_unique<descriptors::LaplacianEigenvalueSkewness>());
    registerDescriptor(std::make_unique<descriptors::NormalizedLaplacianSpectralRadius>());
    registerDescriptor(std::make_unique<descriptors::NormalizedLaplacianSmallestNonzero>());
    registerDescriptor(std::make_unique<descriptors::NormalizedLaplacianLargestEigenvalue>());
    registerDescriptor(std::make_unique<descriptors::NormalizedLaplacianEnergy>());
    registerDescriptor(std::make_unique<descriptors::NormalizedLaplacianTrace>());
    registerDescriptor(std::make_unique<descriptors::DegreeMatrixMaxDegree>());
    registerDescriptor(std::make_unique<descriptors::DegreeMatrixMinDegree>());
    registerDescriptor(std::make_unique<descriptors::DegreeMatrixAvgDegree>());
    registerDescriptor(std::make_unique<descriptors::DegreeMatrixVariance>());
    registerDescriptor(std::make_unique<descriptors::DegreeMatrixSkewness>());
    registerDescriptor(std::make_unique<descriptors::DegreeMatrixEntropy>());
    registerDescriptor(std::make_unique<descriptors::NumberOf2Walks>());
    registerDescriptor(std::make_unique<descriptors::NumberOf3Walks>());
    registerDescriptor(std::make_unique<descriptors::NumberOf4Walks>());
    registerDescriptor(std::make_unique<descriptors::MeanClosed3WalksPerNode>());
    registerDescriptor(std::make_unique<descriptors::MeanClosed4WalksPerNode>());
    registerDescriptor(std::make_unique<descriptors::Walk2Energy>());
    registerDescriptor(std::make_unique<descriptors::Walk3Energy>());
    registerDescriptor(std::make_unique<descriptors::Walk4Energy>());
    registerDescriptor(std::make_unique<descriptors::GraphIrregularityWalkCount>());
    registerDescriptor(std::make_unique<descriptors::TrianglesToPathsRatio>());
    registerDescriptor(std::make_unique<descriptors::MaxSingularValue>());
    registerDescriptor(std::make_unique<descriptors::MinNonZeroSingularValue>());
    registerDescriptor(std::make_unique<descriptors::ConditionNumber>());
    registerDescriptor(std::make_unique<descriptors::SumSingularValues>());
    registerDescriptor(std::make_unique<descriptors::FrobeniusNormSingularValues>());
    registerDescriptor(std::make_unique<descriptors::SingularValueEntropy>());
    registerDescriptor(std::make_unique<descriptors::SingularValueVariance>());
    registerDescriptor(std::make_unique<descriptors::SingularValueSkewness>());
    registerDescriptor(std::make_unique<descriptors::SpectralEffectiveRank>());
    registerDescriptor(std::make_unique<descriptors::NuclearNorm>());
    registerDescriptor(std::make_unique<descriptors::NormalizedEigenvalueGap>());
    registerDescriptor(std::make_unique<descriptors::SumNormalizedEigenvalues>());
    registerDescriptor(std::make_unique<descriptors::VarianceNormalizedEigenvalues>());
    registerDescriptor(std::make_unique<descriptors::CountNormalizedEigenvaluesAboveHalf>());
    registerDescriptor(std::make_unique<descriptors::NormalizedEnergy>());
    registerDescriptor(std::make_unique<descriptors::LargestNormalizedEigenvectorCentrality>());
    registerDescriptor(std::make_unique<descriptors::AverageNormalizedEigenvectorCentrality>());
    registerDescriptor(std::make_unique<descriptors::NormalizedAdjacencyMatrixRank>());
    registerDescriptor(std::make_unique<descriptors::NormalizedAdjacencyEntropy>());
    registerDescriptor(std::make_unique<descriptors::SignlessLaplacianSpectralRadius>());
    registerDescriptor(std::make_unique<descriptors::SignlessLaplacianSmallestEigenvalue>());
    registerDescriptor(std::make_unique<descriptors::SignlessLaplacianEnergy>());
    registerDescriptor(std::make_unique<descriptors::SignlessLaplacianTrace>());
    registerDescriptor(std::make_unique<descriptors::SignlessLaplacianDeterminant>());
    registerDescriptor(std::make_unique<descriptors::SignlessLaplacianZeroEigenvalues>());
    registerDescriptor(std::make_unique<descriptors::SignlessLaplacianPositiveEigenvalues>());
    registerDescriptor(std::make_unique<descriptors::SignlessLaplacianNegativeEigenvalues>());
    registerDescriptor(std::make_unique<descriptors::SignlessLaplacianEigenvalueVariance>());
    registerDescriptor(std::make_unique<descriptors::SignlessLaplacianEigenvalueSkewness>());
    registerDescriptor(std::make_unique<descriptors::WeightedAdjacencySpectralRadius>());
    registerDescriptor(std::make_unique<descriptors::MeanFirstPassageTime>());
    registerDescriptor(std::make_unique<descriptors::CommuteTimeDistance>());
    registerDescriptor(std::make_unique<descriptors::KirchhoffIndexVariance>());
    registerDescriptor(std::make_unique<descriptors::EffectiveGraphResistanceDistribution>());
    registerDescriptor(std::make_unique<descriptors::LocalClusteringCoefficientDistribution>());
    registerDescriptor(std::make_unique<descriptors::GraphRobustnessIndex>());
    registerDescriptor(std::make_unique<descriptors::NormalizedEstradaIndex>());
    registerDescriptor(std::make_unique<descriptors::GraphBipartivityIndex>());
    registerDescriptor(std::make_unique<descriptors::SpanningTreeEntropy>());
    registerDescriptor(std::make_unique<descriptors::GraphIrregularity>());
    registerDescriptor(std::make_unique<descriptors::WienerIndex>());
    registerDescriptor(std::make_unique<descriptors::EstradaIndex>());
    registerDescriptor(std::make_unique<descriptors::NumberSpanningTrees>());
    registerDescriptor(std::make_unique<descriptors::GraphEccentricity>());
    registerDescriptor(std::make_unique<descriptors::SpectralGap>());
    registerDescriptor(std::make_unique<descriptors::TraceMatrixPower2>());
    registerDescriptor(std::make_unique<descriptors::NumberTriangles>());
    registerDescriptor(std::make_unique<descriptors::GraphDiameter>());

    // Enhanced Statistical Features
    registerDescriptor(std::make_unique<descriptors::ChainLengthDescriptor>());
    registerDescriptor(std::make_unique<descriptors::DeviationSymmetryDescriptor>());
    registerDescriptor(std::make_unique<descriptors::LocalComplexityDescriptor>());
    registerDescriptor(std::make_unique<descriptors::CharacterKurtosisDescriptor>());
    registerDescriptor(std::make_unique<descriptors::CharacterSkewnessDescriptor>());

    // Sequence Analysis Features
    registerDescriptor(std::make_unique<descriptors::SequenceNoisinessDescriptor>());
    registerDescriptor(std::make_unique<descriptors::CharacterTransitionMatrixEntropyDescriptor>());
    registerDescriptor(std::make_unique<descriptors::SequenceFractalDimensionDescriptor>());
    registerDescriptor(std::make_unique<descriptors::SubsequenceRepetitionDescriptor>());
    registerDescriptor(std::make_unique<descriptors::MaxRepeatedSubstringDescriptor>());

    // Fragment-Based Features
    // registerDescriptor(std::make_unique<descriptors::FragmentSizeVarianceDescriptor>()); // Zero or near-zero variance
    registerDescriptor(std::make_unique<descriptors::NonCarbonComplexityDescriptor>());
    registerDescriptor(std::make_unique<descriptors::SideChainIndexDescriptor>());
    registerDescriptor(std::make_unique<descriptors::RingBranchRatioDescriptor>());
    registerDescriptor(std::make_unique<descriptors::AromaticAliphaticRatioDescriptor>());
    // Specialized Chemical Pattern Features
    
    // Advanced Entropy/Information Features
    registerDescriptor(std::make_unique<descriptors::MarkovEntropyDescriptor>());
    registerDescriptor(std::make_unique<descriptors::ContextualEntropyDescriptor>());
    registerDescriptor(std::make_unique<descriptors::SyntacticComplexityDescriptor>());
    registerDescriptor(std::make_unique<descriptors::RelativeEntropyDescriptor>());
    registerDescriptor(std::make_unique<descriptors::InformationDensityDescriptor>());

    // Register Vague 5 Descriptors
    registerDescriptor(std::make_unique<descriptors::LongestContinuousSp2FragmentLength>());
    registerDescriptor(std::make_unique<descriptors::ChainBranchDepthMaximum>());
    registerDescriptor(std::make_unique<descriptors::ChainBranchingFrequency>());
    registerDescriptor(std::make_unique<descriptors::LongestAlternatingSingleDoubleBondPath>());
    registerDescriptor(std::make_unique<descriptors::AverageRingStrainProxy>());
    registerDescriptor(std::make_unique<descriptors::FractionBackboneCarbon>());
    registerDescriptor(std::make_unique<descriptors::TopologicalAsymmetryIndex>());
    registerDescriptor(std::make_unique<descriptors::FractionConjugatedBondsBackbone>());
    registerDescriptor(std::make_unique<descriptors::MaxChainNoHetero>());
    registerDescriptor(std::make_unique<descriptors::MeanRingSizeWeightedDegree>());
    registerDescriptor(std::make_unique<descriptors::PolarClusterCount>());
    registerDescriptor(std::make_unique<descriptors::MaxElectronegativityGradientChain>());
    registerDescriptor(std::make_unique<descriptors::FractionPolarizableHeavy>());
    registerDescriptor(std::make_unique<descriptors::CumulativeENDifferenceBackbone>());
    registerDescriptor(std::make_unique<descriptors::HeavyAtomLocalChargeSkewness>());
    registerDescriptor(std::make_unique<descriptors::ElectronegativeEndChainCount>());
    registerDescriptor(std::make_unique<descriptors::FractionBondsLargeENGap>());
    registerDescriptor(std::make_unique<descriptors::SymmetryElectronegativeDistribution>());
    registerDescriptor(std::make_unique<descriptors::FractionHeteroatomsTerminal>());
    registerDescriptor(std::make_unique<descriptors::ConjugatedSystemSizeMax>());
    registerDescriptor(std::make_unique<descriptors::FractionAromaticAtomsConnectedChains>());
    registerDescriptor(std::make_unique<descriptors::MeanDistanceAromaticSystems>());
    registerDescriptor(std::make_unique<descriptors::LengthLargestFullyConjugatedRingSystem>());
    registerDescriptor(std::make_unique<descriptors::RatioConjugatedCarbonsTotalCarbons>());
    registerDescriptor(std::make_unique<descriptors::NormalizedConjugationPathwayDensity>());
    registerDescriptor(std::make_unique<descriptors::EstimatedMolecularAspectRatio>());
    registerDescriptor(std::make_unique<descriptors::TerminalHeavyAtomDispersion>());
    registerDescriptor(std::make_unique<descriptors::ConnectivityCorePeripheryGradient>());
    registerDescriptor(std::make_unique<descriptors::HeavyAtomPlanarityHeuristic>());
    registerDescriptor(std::make_unique<descriptors::RingCoreChainPeripherySizeRatio>());
    registerDescriptor(std::make_unique<descriptors::RingChainAlternationFrequency>());
    registerDescriptor(std::make_unique<descriptors::PotentialDeprotonationSiteCount>());
    registerDescriptor(std::make_unique<descriptors::PotentialProtonationSiteCount>());
    registerDescriptor(std::make_unique<descriptors::FractionNonTerminalChargedSites>());
    registerDescriptor(std::make_unique<descriptors::MeanTopologicalDistanceIonizableSites>());
    registerDescriptor(std::make_unique<descriptors::RingIonizableSiteDensity>());
    registerDescriptor(std::make_unique<descriptors::AromaticProtonDonorAcceptorRatioHeuristic>());
    registerDescriptor(std::make_unique<descriptors::HydrophilicAtomClusterSizeMean>());
    registerDescriptor(std::make_unique<descriptors::HydrophobicIslandSizeMax>());
    registerDescriptor(std::make_unique<descriptors::HydrophobicPolarSurfaceInterfaceEstimation>());
    registerDescriptor(std::make_unique<descriptors::NormalizedHydrophobicPolarRatioEstimate>());
    registerDescriptor(std::make_unique<descriptors::HydrophobicAtomPathLengthMean>());
    registerDescriptor(std::make_unique<descriptors::HeteroatomConnectivityIndex>());
    registerDescriptor(std::make_unique<descriptors::RingFusionDensity>());
    registerDescriptor(std::make_unique<descriptors::AverageBondPolarity>());
    registerDescriptor(std::make_unique<descriptors::StericHindranceProxy>());
    registerDescriptor(std::make_unique<descriptors::MolecularFlexibilityProxy>());

    // SELFIES descriptors
    registerDescriptor(std::make_unique<descriptors::SelfiesTokenCountDescriptor>());
    registerDescriptor(std::make_unique<descriptors::SelfiesTokenDistributionDescriptor>());
    registerDescriptor(std::make_unique<descriptors::SelfiesAverageComplexityDescriptor>());
    registerDescriptor(std::make_unique<descriptors::SelfiesBranchingDepthDescriptor>());
    registerDescriptor(std::make_unique<descriptors::SelfiesMaxBranchingFanoutDescriptor>());
    registerDescriptor(std::make_unique<descriptors::SelfiesEmptyBranchesDescriptor>());
    registerDescriptor(std::make_unique<descriptors::SelfiesRingTokenFrequencyDescriptor>());
    registerDescriptor(std::make_unique<descriptors::SelfiesTokenEntropyDescriptor>());
    registerDescriptor(std::make_unique<descriptors::SelfiesStringLengthDescriptor>());
    registerDescriptor(std::make_unique<descriptors::SelfiesBranchTokenRatioDescriptor>());
    registerDescriptor(std::make_unique<descriptors::SelfiesRingTokenRatioDescriptor>());
    registerDescriptor(std::make_unique<descriptors::SelfiesBigramCountDescriptor>());
    registerDescriptor(std::make_unique<descriptors::SelfiesTrigramCountDescriptor>());
    registerDescriptor(std::make_unique<descriptors::SelfiesTokenRepetitionRateDescriptor>());
    registerDescriptor(std::make_unique<descriptors::SelfiesAromaticTokenRatioDescriptor>());
    registerDescriptor(std::make_unique<descriptors::SelfiesChargeTokenRatioDescriptor>());
    registerDescriptor(std::make_unique<descriptors::SelfiesAliphaticChainLengthDescriptor>());
    registerDescriptor(std::make_unique<descriptors::SelfiesHeteroatomRatioDescriptor>());
    registerDescriptor(std::make_unique<descriptors::SelfiesMaxRingSizeDescriptor>());
    registerDescriptor(std::make_unique<descriptors::SelfiesStereoTokenRatioDescriptor>());
    registerDescriptor(std::make_unique<descriptors::SelfiesGrammarDepthComplexityDescriptor>());
    registerDescriptor(std::make_unique<descriptors::SelfiesVocabularyDiversityDescriptor>());
    registerDescriptor(std::make_unique<descriptors::SelfiesBondComplexityDescriptor>());
    registerDescriptor(std::make_unique<descriptors::SelfiesSequenceSymmetryDescriptor>());
    registerDescriptor(std::make_unique<descriptors::SelfiesTokenLengthVarianceDescriptor>());
    registerDescriptor(std::make_unique<descriptors::SelfiesFormalChargeEntropyDescriptor>());
    registerDescriptor(std::make_unique<descriptors::SelfiesElementDiversityDescriptor>());
    registerDescriptor(std::make_unique<descriptors::SelfiesLevenshteinCompressionDescriptor>());
    registerDescriptor(std::make_unique<descriptors::SelfiesLoopComplexityDescriptor>());
    registerDescriptor(std::make_unique<descriptors::SelfiesUnpairedElectronRatioDescriptor>());

    // Vague6 descriptors (Combined with Vague 4 registrations above where applicable)
    registerDescriptor(std::make_unique<descriptors::HeteroatomFractionOfRings>());
    registerDescriptor(std::make_unique<descriptors::AromaticRingRatio>());
    registerDescriptor(std::make_unique<descriptors::RingComplexityIndex>());
    registerDescriptor(std::make_unique<descriptors::ElectronegativeAtomDensity>());
    registerDescriptor(std::make_unique<descriptors::ChargeAsymmetryIndex>());
    registerDescriptor(std::make_unique<descriptors::TopologicalPolarSurfaceEfficiency>());
    registerDescriptor(std::make_unique<descriptors::HeterocyclicRingCount>());
    registerDescriptor(std::make_unique<descriptors::CarbocyclicRingCount>());
    registerDescriptor(std::make_unique<descriptors::FusedRingSystems>());
    registerDescriptor(std::make_unique<descriptors::AverageFusedRingSize>());
    registerDescriptor(std::make_unique<descriptors::HeterocyclicFusedRingRatio>());
    registerDescriptor(std::make_unique<descriptors::BridgedRingCount>());
    registerDescriptor(std::make_unique<descriptors::SpiroRingCount>());
    registerDescriptor(std::make_unique<descriptors::MacrocyclicRingCount>());
    registerDescriptor(std::make_unique<descriptors::HydrogenBondAcceptorDensity>());
    registerDescriptor(std::make_unique<descriptors::HydrogenBondDonorDensity>());
    registerDescriptor(std::make_unique<descriptors::RotatableBondDensity>());
    registerDescriptor(std::make_unique<descriptors::FunctionalGroupDiversity>());
    registerDescriptor(std::make_unique<descriptors::AtomTypeEntropy>());
    registerDescriptor(std::make_unique<descriptors::BondTypeEntropy>());
    registerDescriptor(std::make_unique<descriptors::RingTypeEntropy>());
    registerDescriptor(std::make_unique<descriptors::ElectronegativityVariance>());
    // registerDescriptor(std::make_unique<descriptors::TopologicalComplexityIndex>());
    // registerDescriptor(std::make_unique<descriptors::MolecularSymmetryIndex>());
    // registerDescriptor(std::make_unique<descriptors::AtomConnectivityEquality>());
    // registerDescriptor(std::make_unique<descriptors::BranchingIndex>());
    // registerDescriptor(std::make_unique<descriptors::ChargeDistributionSkewness>());
    //registerDescriptor(std::make_unique<descriptors::ElectronegativitySkewness>());
    // registerDescriptor(std::make_unique<descriptors::AromaticAtomCountToRingRatio>());
    // registerDescriptor(std::make_unique<descriptors::ChargeToMassRatio>());
    // registerDescriptor(std::make_unique<descriptors::HeteroatomDistribution>());
    // registerDescriptor(std::make_unique<descriptors::PolarizabilityRatio>());
    // registerDescriptor(std::make_unique<descriptors::RingAtomToNonRingAtomRatio>());
    //registerDescriptor(std::make_unique<descriptors::RingSizeDistributionVariance>());
    // registerDescrip(std::make_unique<descriptors::RingBondFraction>());
    // registerDescriptor(std::make_unique<descriptors::RingFusionDensity>()); // Removed duplicate registration
    //registerDescriptor(std::make_unique<descriptors::AtomElectronegativityImbalance>());
    //registerDescriptor(std::make_unique<descriptors::ConnectivityAsymmetryIndex>());
    //registerDescriptor(std::make_unique<descriptors::BondDistancePathLength>());
    //registerDescriptor(std::make_unique<descriptors::HeteroatomPathFraction>());
    //registerDescriptor(std::make_unique<descriptors::ElectronegativePathLength>());
    // registerDescriptor(std::make_unique<descriptors::ConsecutiveHeteroatomCount>());
    registerDescriptor(std::make_unique<descriptors::ChiralAtomSkew>());
    registerDescriptor(std::make_unique<descriptors::RingAtomUnsaturationRatio>());
    registerDescriptor(std::make_unique<descriptors::TerminalHeteroBondOrderSum>());
    registerDescriptor(std::make_unique<descriptors::AromaticNonAromaticBondCount>());
    registerDescriptor(std::make_unique<descriptors::RingAtomChargeVariance>());


    // *** Register New pKa Descriptors ***
    registerDescriptor(std::make_unique<descriptors::AcidityTypeDescriptor>());
    registerDescriptor(std::make_unique<descriptors::PkaEstimate1Descriptor>());
    registerDescriptor(std::make_unique<descriptors::PkaEstimate2Descriptor>());
    registerDescriptor(std::make_unique<descriptors::IsAcidDescriptor>());
    registerDescriptor(std::make_unique<descriptors::IsBaseDescriptor>());
    registerDescriptor(std::make_unique<descriptors::IsNeutralDescriptor>());
    registerDescriptor(std::make_unique<descriptors::IsZwitterionDescriptor>());

    // Register Morgan Descriptors
    registerDescriptor(std::make_unique<descriptors::MorganBitDensityRatioDescriptor>());
    registerDescriptor(std::make_unique<descriptors::MorganFragmentUniquenessScoreDescriptor>());
    registerDescriptor(std::make_unique<descriptors::MorganRadiusInformationRatioDescriptor>());
    registerDescriptor(std::make_unique<descriptors::MorganBitClusteringCoefficientDescriptor>());
    registerDescriptor(std::make_unique<descriptors::MorganFrequentFragmentEntropyDescriptor>());
    registerDescriptor(std::make_unique<descriptors::MorganFingerprintAsymmetryIndexDescriptor>());
    registerDescriptor(std::make_unique<descriptors::MorganBitTransitionRateDescriptor>());
    registerDescriptor(std::make_unique<descriptors::MorganFragmentSizeDistributionSkewnessDescriptor>());
    registerDescriptor(std::make_unique<descriptors::MorganLongestCommonSubstructureScoreDescriptor>());
    registerDescriptor(std::make_unique<descriptors::MorganPharmacophorePatternDensityDescriptor>());
    registerDescriptor(std::make_unique<descriptors::MorganFragmentDiversityScoreDescriptor>());
    registerDescriptor(std::make_unique<descriptors::MorganBitCorrelationCoefficientDescriptor>());
    registerDescriptor(std::make_unique<descriptors::MorganPatternRecurrenceFrequencyDescriptor>());
    registerDescriptor(std::make_unique<descriptors::MorganBitSimilarityToReferenceSetDescriptor>());
    registerDescriptor(std::make_unique<descriptors::MorganFragmentComplexityScoreDescriptor>());
    registerDescriptor(std::make_unique<descriptors::MorganInformationContentDensityDescriptor>());
    registerDescriptor(std::make_unique<descriptors::MorganEnvironmentVariabilityIndexDescriptor>());
    registerDescriptor(std::make_unique<descriptors::MorganBitPositionImportanceScoreDescriptor>());
    registerDescriptor(std::make_unique<descriptors::MorganRingSystemRepresentationScoreDescriptor>());
    registerDescriptor(std::make_unique<descriptors::MorganBitProbabilityDistributionEntropyDescriptor>());
    registerDescriptor(std::make_unique<descriptors::MorganFragmentElectronegativitySpectrumDescriptor>());
    registerDescriptor(std::make_unique<descriptors::MorganBitPositionCorrelationMatrixDescriptor>());
    registerDescriptor(std::make_unique<descriptors::MorganFragmentConnectivityPatternDescriptor>());
    registerDescriptor(std::make_unique<descriptors::MorganBitOccurrenceFrequencySkewnessDescriptor>());
    registerDescriptor(std::make_unique<descriptors::MorganSubstructureHeterogeneityIndexDescriptor>());
    registerDescriptor(std::make_unique<descriptors::MorganBitPolarityDistributionDescriptor>());
    registerDescriptor(std::make_unique<descriptors::MorganFragmentSimilarityNetworkDescriptor>());
    registerDescriptor(std::make_unique<descriptors::MorganBitPositionInformationGainDescriptor>());
    registerDescriptor(std::make_unique<descriptors::MorganSubstructureDiversityGradientDescriptor>());
    registerDescriptor(std::make_unique<descriptors::MorganFingerprintSymmetryScoreDescriptor>());

    // *** Register Vague8 Descriptors ***
    registerDescriptor(std::make_unique<descriptors::TopologicalChargeDistributionSkewnessDescriptor>());
    registerDescriptor(std::make_unique<descriptors::ElementNeighborhoodDiversityDescriptor>());
    registerDescriptor(std::make_unique<descriptors::BondOrderAlternationPatternDescriptor>());
    registerDescriptor(std::make_unique<descriptors::AvgShortestPathBetweenDiffFuncGroupsDescriptor>());
    registerDescriptor(std::make_unique<descriptors::AvgHeteroClusterSizeDescriptor>());
    registerDescriptor(std::make_unique<descriptors::EntropyOfRingSubstCountDistDescriptor>());
    registerDescriptor(std::make_unique<descriptors::EntropyOfFuncGroupDistancesDescriptor>());
    registerDescriptor(std::make_unique<descriptors::ChainBranchingPatternDescriptor>());
    registerDescriptor(std::make_unique<descriptors::ElectronegativeAtomNeighborhoodScoreDescriptor>());
    registerDescriptor(std::make_unique<descriptors::TopologicalChargeSeparationIndexDescriptor>());
    registerDescriptor(std::make_unique<descriptors::RingSystemConnectivityPatternDescriptor>());
    registerDescriptor(std::make_unique<descriptors::HeteroatomPositionEntropyDescriptor>());
    registerDescriptor(std::make_unique<descriptors::BondOrderTransitionFrequencyDescriptor>());
    registerDescriptor(std::make_unique<descriptors::AtomNeighborhoodElectronegativityGradientDescriptor>());
    registerDescriptor(std::make_unique<descriptors::ChainLengthDistributionEntropyDescriptor>());
    // registerDescriptor(std::make_unique<descriptors::RingSubstitutionSymmetryDescriptor>()); // Temporarily comment out
    registerDescriptor(std::make_unique<descriptors::HeteroatomSequencePatternsDescriptor>());
    registerDescriptor(std::make_unique<descriptors::FunctionalGroupIsolationTopologyDescriptor>());
    registerDescriptor(std::make_unique<descriptors::HydrophobicPatchConnectivityDescriptor>());
    registerDescriptor(std::make_unique<descriptors::BondTopologicalEnvironmentFingerprintDescriptor>());
    registerDescriptor(std::make_unique<descriptors::ChiralCenterTopologicalDistributionDescriptor>());
    registerDescriptor(std::make_unique<descriptors::RingFusionPatternCodeDescriptor>());
    registerDescriptor(std::make_unique<descriptors::ElectronegativityTopologicalMomentDescriptor>());
    registerDescriptor(std::make_unique<descriptors::AtomicRadiiVarianceInNeighborhoodsDescriptor>());
    registerDescriptor(std::make_unique<descriptors::HeteroatomBondingPatternCodeDescriptor>());
    registerDescriptor(std::make_unique<descriptors::SubstructureFrequencySpectrumDescriptor>());
    registerDescriptor(std::make_unique<descriptors::FormalChargeNeighborhoodPatternDescriptor>());

    // Register Solubility Descriptor
    registerDescriptor(std::make_unique<desfact::descriptors::SolubilityDescriptor>());

    // Register Count descriptors
    registerDescriptor(std::make_unique<desfact::descriptors::HalogenCount>());
    registerDescriptor(std::make_unique<desfact::descriptors::CarbonCount>());
    registerDescriptor(std::make_unique<desfact::descriptors::HydrogenCount>());
    registerDescriptor(std::make_unique<desfact::descriptors::OxygenCount>());
    registerDescriptor(std::make_unique<desfact::descriptors::NitrogenCount>());
    registerDescriptor(std::make_unique<desfact::descriptors::SulfurCount>());
    registerDescriptor(std::make_unique<desfact::descriptors::AromaticAtomCount>());
    registerDescriptor(std::make_unique<desfact::descriptors::NonAromaticAtomCount>());
    registerDescriptor(std::make_unique<desfact::descriptors::DoubleBondCount>());
    registerDescriptor(std::make_unique<desfact::descriptors::TripleBondCount>());
    registerDescriptor(std::make_unique<desfact::descriptors::SingleBondCount>());
    registerDescriptor(std::make_unique<desfact::descriptors::AcidicFunctionCount>());
    registerDescriptor(std::make_unique<desfact::descriptors::BasicFunctionCount>());
    registerDescriptor(std::make_unique<desfact::descriptors::ENAtoms2BondsFromAcidic>());
    registerDescriptor(std::make_unique<desfact::descriptors::ENAtoms2BondsFromBasic>());
    registerDescriptor(std::make_unique<desfact::descriptors::ENAtoms3BondsFromAcidic>());
    registerDescriptor(std::make_unique<desfact::descriptors::ENAtoms3BondsFromBasic>());

    // Register new counts
    registerDescriptor(std::make_unique<desfact::descriptors::LongestCSequenceSmiles>());
    registerDescriptor(std::make_unique<desfact::descriptors::UppercaseCountSmiles>());
    registerDescriptor(std::make_unique<desfact::descriptors::LowercaseCountSmiles>());
    registerDescriptor(std::make_unique<desfact::descriptors::AtomVolumeSum>());

    // Register new count descriptors
    registerDescriptor(std::make_unique<desfact::descriptors::HeavyAtomCount>());
    registerDescriptor(std::make_unique<desfact::descriptors::HeteroatomCount>());
    registerDescriptor(std::make_unique<desfact::descriptors::RingAtomCount>());
    registerDescriptor(std::make_unique<desfact::descriptors::RingCount>());
    registerDescriptor(std::make_unique<desfact::descriptors::ChiralCenterCount>());
    registerDescriptor(std::make_unique<desfact::descriptors::FormalChargeCount>());
    registerDescriptor(std::make_unique<desfact::descriptors::PositiveChargeCount>());
    registerDescriptor(std::make_unique<desfact::descriptors::NegativeChargeCount>());
    registerDescriptor(std::make_unique<desfact::descriptors::RotatableBondCount>());
    registerDescriptor(std::make_unique<desfact::descriptors::BridgeheadAtomCount>());

    // Image-based descriptors
    registerDescriptor(std::make_unique<descriptors::BoundingBoxAreaDescriptor>());
    registerDescriptor(std::make_unique<descriptors::MoleculeWidthDescriptor>());
    registerDescriptor(std::make_unique<descriptors::MoleculeHeightDescriptor>());
    registerDescriptor(std::make_unique<descriptors::AspectRatioDescriptor>());
    registerDescriptor(std::make_unique<descriptors::AtomDensityDescriptor>());
    registerDescriptor(std::make_unique<descriptors::MeanAtomRadiusDescriptor>());
    registerDescriptor(std::make_unique<descriptors::MaxAtomRadiusDescriptor>());
    registerDescriptor(std::make_unique<descriptors::MinAtomRadiusDescriptor>());
    registerDescriptor(std::make_unique<descriptors::MeanBondLengthDescriptor>());
    registerDescriptor(std::make_unique<descriptors::MaxBondLengthDescriptor>());
    registerDescriptor(std::make_unique<descriptors::MinBondLengthDescriptor>());
    registerDescriptor(std::make_unique<descriptors::BondLengthStdDescriptor>());
    registerDescriptor(std::make_unique<descriptors::MeanAtomAtomDistDescriptor>());
    registerDescriptor(std::make_unique<descriptors::MeanAtomLuminanceDescriptor>());
    registerDescriptor(std::make_unique<descriptors::MaxAtomAtomDistDescriptor>());
    registerDescriptor(std::make_unique<descriptors::MinAtomAtomDistDescriptor>());
    registerDescriptor(std::make_unique<descriptors::MeanAtomAreaDescriptor>());
    registerDescriptor(std::make_unique<descriptors::MedianAtomRadiusDescriptor>());
    registerDescriptor(std::make_unique<descriptors::MedianBondLengthDescriptor>());
    registerDescriptor(std::make_unique<descriptors::AtomRadiusRangeDescriptor>());
    registerDescriptor(std::make_unique<descriptors::AtomAreaFractionDescriptor>());
    registerDescriptor(std::make_unique<descriptors::BondCoverageFractionDescriptor>());
    registerDescriptor(std::make_unique<descriptors::AtomPackingDensityDescriptor>());
    registerDescriptor(std::make_unique<descriptors::AtomMassXDescriptor>());
    registerDescriptor(std::make_unique<descriptors::AtomMassYDescriptor>());
    registerDescriptor(std::make_unique<descriptors::AtomMassDistDescriptor>());
    registerDescriptor(std::make_unique<descriptors::AtomRadiusStdDescriptor>());
    registerDescriptor(std::make_unique<descriptors::MeanBondAngleDescriptor>());
    registerDescriptor(std::make_unique<descriptors::BondAngleStdDescriptor>());
    registerDescriptor(std::make_unique<descriptors::AtomXStdDescriptor>());
    registerDescriptor(std::make_unique<descriptors::AtomYStdDescriptor>());
    registerDescriptor(std::make_unique<descriptors::AtomLuminanceStdDescriptor>());
    registerDescriptor(std::make_unique<descriptors::BondLenMADDescriptor>());
    registerDescriptor(std::make_unique<descriptors::AtomRadiusMADDescriptor>());
    registerDescriptor(std::make_unique<descriptors::AtomAreaStdDescriptor>());
    registerDescriptor(std::make_unique<descriptors::AtomRadiusCVDescriptor>());
    registerDescriptor(std::make_unique<descriptors::AtomAreaCVDescriptor>());
    registerDescriptor(std::make_unique<descriptors::AtomLuminanceCVDescriptor>());
    registerDescriptor(std::make_unique<descriptors::AtomAreaRangeDescriptor>());
    registerDescriptor(std::make_unique<descriptors::AtomAreaMedianDescriptor>());
    registerDescriptor(std::make_unique<descriptors::AtomAreaMADDescriptor>());
    registerDescriptor(std::make_unique<descriptors::AtomAreaHullFracDescriptor>());
    registerDescriptor(std::make_unique<descriptors::AtomAreaCenterFracDescriptor>());
    registerDescriptor(std::make_unique<descriptors::StdAllAtomDistDescriptor>());
    registerDescriptor(std::make_unique<descriptors::MinBondAngleDescriptor>());
    registerDescriptor(std::make_unique<descriptors::MaxBondAngleDescriptor>());
    registerDescriptor(std::make_unique<descriptors::MedianBondAngleDescriptor>());
    registerDescriptor(std::make_unique<descriptors::MeanBondColorDiffDescriptor>());
    registerDescriptor(std::make_unique<descriptors::MeanBondLuminanceDiffDescriptor>());
    registerDescriptor(std::make_unique<descriptors::MeanBondLenDiffColorDescriptor>());
    registerDescriptor(std::make_unique<descriptors::MeanDistSameColorDescriptor>());
    registerDescriptor(std::make_unique<descriptors::MeanDistDiffColorDescriptor>());

    // Image-based descriptors (continued)
    registerDescriptor(std::make_unique<descriptors::MeanBondAngleDiffColorDescriptor>());
    registerDescriptor(std::make_unique<descriptors::MeanBondRadiiDiffDescriptor>());
    registerDescriptor(std::make_unique<descriptors::MeanBondLuminanceDiff2Descriptor>()); // Note: Duplicate functionality?
    registerDescriptor(std::make_unique<descriptors::MeanAngleHighDegreeDescriptor>());
    registerDescriptor(std::make_unique<descriptors::BoundaryAtomRatioDescriptor>());
    registerDescriptor(std::make_unique<descriptors::MeanAtomEccentricityDescriptor>());
    registerDescriptor(std::make_unique<descriptors::MolecularDiameterDescriptor>());
    registerDescriptor(std::make_unique<descriptors::MolecularRadiusDescriptor>());
    registerDescriptor(std::make_unique<descriptors::RadiusOfGyrationDescriptor>());
    registerDescriptor(std::make_unique<descriptors::MolecularSphericityDescriptor>());
    registerDescriptor(std::make_unique<descriptors::BondLenToAtomRadiusRatioDescriptor>());
    registerDescriptor(std::make_unique<descriptors::EdgeDensityDescriptor>());
    registerDescriptor(std::make_unique<descriptors::PlanarityMeasureDescriptor>());
    registerDescriptor(std::make_unique<descriptors::MeanNearestNeighborDistDescriptor>());
    registerDescriptor(std::make_unique<descriptors::MeanAtomsInRadiusDescriptor>());
    registerDescriptor(std::make_unique<descriptors::VarNearestNeighborDistDescriptor>());
    registerDescriptor(std::make_unique<descriptors::MeanBondToBondAngleDescriptor>());
    registerDescriptor(std::make_unique<descriptors::MomentOfInertiaDescriptor>());
    registerDescriptor(std::make_unique<descriptors::MeanAtomCentralityDescriptor>());
    registerDescriptor(std::make_unique<descriptors::TerminalAtomCountDescriptor>());
    registerDescriptor(std::make_unique<descriptors::JunctionAtomCountDescriptor>());
    registerDescriptor(std::make_unique<descriptors::TerminalToJunctionRatioDescriptor>());
    registerDescriptor(std::make_unique<descriptors::WienerIndexDescriptor>()); // Image-based Wiener Index
    registerDescriptor(std::make_unique<descriptors::RingCountDescriptor>()); // Image-based Ring Count
    registerDescriptor(std::make_unique<descriptors::AcidicCenterCountDescriptor>());
    registerDescriptor(std::make_unique<descriptors::BasicCenterCountDescriptor>());
    registerDescriptor(std::make_unique<descriptors::AcidicToBasicRatioDescriptor>());
    registerDescriptor(std::make_unique<descriptors::MeanAcidicAcidicDistDescriptor>());
    registerDescriptor(std::make_unique<descriptors::MeanBasicBasicDistDescriptor>());
    registerDescriptor(std::make_unique<descriptors::MeanAcidicBasicDistDescriptor>());
    registerDescriptor(std::make_unique<descriptors::MinAcidicBasicDistDescriptor>());
    registerDescriptor(std::make_unique<descriptors::MeanAtomsNearAcidicDescriptor>());
    registerDescriptor(std::make_unique<descriptors::MeanAtomsNearBasicDescriptor>());
    registerDescriptor(std::make_unique<descriptors::MeanAcidicRadiusDescriptor>());
    registerDescriptor(std::make_unique<descriptors::MeanBasicRadiusDescriptor>());
    registerDescriptor(std::make_unique<descriptors::MeanAcidicCentroidDistDescriptor>());
    registerDescriptor(std::make_unique<descriptors::MeanBasicCentroidDistDescriptor>());
    registerDescriptor(std::make_unique<descriptors::FracAcidicOnHullDescriptor>());
    registerDescriptor(std::make_unique<descriptors::FracBasicOnHullDescriptor>());
    registerDescriptor(std::make_unique<descriptors::MeanAcidicAngleDescriptor>());
    registerDescriptor(std::make_unique<descriptors::MeanBasicAngleDescriptor>());
    registerDescriptor(std::make_unique<descriptors::MeanAcidicLuminanceDescriptor>());
    registerDescriptor(std::make_unique<descriptors::MeanBasicLuminanceDescriptor>());
    registerDescriptor(std::make_unique<descriptors::AcidicBasicLuminanceDiffDescriptor>());
    registerDescriptor(std::make_unique<descriptors::FractalDimensionDescriptor>());
    registerDescriptor(std::make_unique<descriptors::ForegroundRatioDescriptor>());
    registerDescriptor(std::make_unique<descriptors::AverageColorDescriptor>());
    registerDescriptor(std::make_unique<descriptors::ColorVarianceDescriptor>());
    registerDescriptor(std::make_unique<descriptors::ImageCenterXDescriptor>());
    registerDescriptor(std::make_unique<descriptors::ImageCenterYDescriptor>());
    registerDescriptor(std::make_unique<descriptors::ImageOrientationDescriptor>());
    registerDescriptor(std::make_unique<descriptors::MeanLenSingleDescriptor>());
    registerDescriptor(std::make_unique<descriptors::MeanLenDoubleDescriptor>());
    registerDescriptor(std::make_unique<descriptors::MeanLenTripleDescriptor>());
    registerDescriptor(std::make_unique<descriptors::MeanLenAromaticDescriptor>());
    registerDescriptor(std::make_unique<descriptors::FracDeg1Descriptor>());
    registerDescriptor(std::make_unique<descriptors::FracDeg2Descriptor>());
    registerDescriptor(std::make_unique<descriptors::FracDeg3Descriptor>());
    registerDescriptor(std::make_unique<descriptors::FracDeg4Descriptor>());
    registerDescriptor(std::make_unique<descriptors::MeanBondOrderDescriptor>());
    registerDescriptor(std::make_unique<descriptors::FracDoubleBondsDescriptor>());
    registerDescriptor(std::make_unique<descriptors::FracTripleBondsDescriptor>());
    registerDescriptor(std::make_unique<descriptors::FracAromaticBondsDescriptor>());
    registerDescriptor(std::make_unique<descriptors::MeanAtomDegreeDescriptor>());
    registerDescriptor(std::make_unique<descriptors::MaxAtomDegreeDescriptor>());
    registerDescriptor(std::make_unique<descriptors::MinAtomDegreeDescriptor>());
    registerDescriptor(std::make_unique<descriptors::MeanDistToCentroidDescriptor>());
    registerDescriptor(std::make_unique<descriptors::StdDistToCentroidDescriptor>());


    globalLogger.info("Initialized descriptor factory with " +
                     std::to_string(descriptors.size()) + " descriptors");
}

void DescriptorFactory::registerDescriptor(std::unique_ptr<Descriptor> descriptor) {
    if (descriptor) {
        descriptors[descriptor->getName()] = std::move(descriptor);
    }
}

const Descriptor* DescriptorFactory::getDescriptor(const std::string& name) const {
    auto it = descriptors.find(name);
    if (it != descriptors.end()) {
        return it->second.get();
    }
    return nullptr;
}

std::vector<std::string> DescriptorFactory::getAvailableDescriptors() const {
    std::vector<std::string> names;
    names.reserve(descriptors.size());
    
    for (const auto& [name, _] : descriptors) {
        names.push_back(name);
    }
    
    return names;
}

std::vector<std::string> DescriptorFactory::getCudaEnabledDescriptors() const {
    std::vector<std::string> names;
    
    for (const auto& [name, descriptor] : descriptors) {
        if (descriptor->isCudaSupported()) {
            names.push_back(name);
        }
    }
    
    return names;
}

std::variant<double, int, std::string> DescriptorFactory::calculate(
    const std::string& descriptorName, const Molecule& mol) const {
    
    const Descriptor* descriptor = getDescriptor(descriptorName);
    if (!descriptor) {
        throw DescriptorException(
            "Unknown descriptor: " + descriptorName,
            ErrorCode::NOT_IMPLEMENTED
        );
    }
    
    return descriptor->calculate(mol);
}

std::unordered_map<std::string, std::vector<std::variant<double, int, std::string>>>
DescriptorFactory::calculateBatch(const std::vector<std::string>& descriptorNames,
                                const MoleculeBatch& batch) const {

    std::unordered_map<std::string, std::vector<std::variant<double, int, std::string>>> finalResults;
    const auto& molecules = batch.getMolecules();
    size_t numMolecules = molecules.size();

    if (numMolecules == 0) {
        return finalResults;
    }

    // Determine result type for each descriptor *once*
    std::unordered_map<std::string, std::type_index> resultTypes;
    const Molecule* firstValidMol = nullptr;

    for(const auto& m : molecules){
        if (m.isValid()) {
            firstValidMol = &m;
            break;
        }
    }

    if (firstValidMol) {
         for(const auto& name : descriptorNames) {
             try {
                 // Ensure the descriptor exists before calculating
                 const Descriptor* desc = getDescriptor(name);
                 if (!desc) {
                      globalLogger.warning("Cannot determine type for unknown descriptor: " + name);
                      continue; // Skip type determination for unknown descriptors
                 }
                 auto sampleResult = desc->calculate(*firstValidMol);
                 resultTypes.emplace(name, sampleResult.index() == 0 ? typeid(double) :
                                          sampleResult.index() == 1 ? typeid(int) :
                                          typeid(std::string));
             } catch (const std::exception& e) {
                  globalLogger.warning("Error determining type for descriptor " + name + " using first valid molecule: " + e.what());
                  // Could potentially default to string here or leave it out of resultTypes
             }
         }
    } else {
         globalLogger.warning("Could not find a valid molecule in the batch to determine descriptor types.");
         // All molecules are invalid, results will be filled with errors anyway.
    }


    // Initialize result vectors
    for (const auto& name : descriptorNames) {
        finalResults[name].resize(numMolecules);
        auto it = resultTypes.find(name); // Use find instead of [] or count
        if (it != resultTypes.end()) {
            const std::type_index& type = it->second; // Get type safely
            if (type == typeid(double)) {
                std::fill(finalResults[name].begin(), finalResults[name].end(), std::numeric_limits<double>::quiet_NaN());
            } else if (type == typeid(int)) {
                std::fill(finalResults[name].begin(), finalResults[name].end(), -1);
            } else { // Default to string error if type is string or unknown/failed
                std::fill(finalResults[name].begin(), finalResults[name].end(), std::string("Error: Init/Str"));
            }
        } else {
            // Type determination failed or descriptor was unknown
             globalLogger.debug("Type unknown for descriptor '" + name + "', initializing with string error.");
            std::fill(finalResults[name].begin(), finalResults[name].end(), std::string("Error: Type Unknown"));
        }
    }

    // Verify all descriptors exist once before processing (already done partly in type determination)
    for (const auto& name : descriptorNames) {
        if (!getDescriptor(name)) {
            globalLogger.error("Unknown descriptor requested in batch: " + name);
            throw DescriptorException(
                "Unknown descriptor requested in batch: " + name,
                ErrorCode::NOT_IMPLEMENTED
            );
        }
    }

#ifdef WITH_TBB
    int numThreads = globalConfig.numThreads > 0 ? globalConfig.numThreads : tbb::this_task_arena::max_concurrency();
    if (numThreads > 1) {
        tbb::parallel_for(tbb::blocked_range<size_t>(0, numMolecules),
            [&](const tbb::blocked_range<size_t>& range) {
                for (size_t i = range.begin(); i != range.end(); ++i) {
                    const auto& mol = molecules[i];
                    if (!mol.isValid()) {
                        // Error value already set during initialization
                        continue;
                    }

                    for (const auto& name : descriptorNames) {
                        const Descriptor* descriptor = getDescriptor(name);
                        if (!descriptor) continue;

                        try {
                            finalResults[name][i] = descriptor->calculate(mol);
                        } catch (const std::exception& e) {
                            globalLogger.debug("Calculation Error for SMILES " + mol.getOriginalSmiles() + ", Desc: " + name + " - " + e.what());
                            // Error value was pre-filled, just log
                        } catch (...) {
                            globalLogger.debug("Unknown Calculation Error for SMILES " + mol.getOriginalSmiles() + ", Desc: " + name);
                            // Error value was pre-filled, just log
                        }
                    }
                }
            }
        );
    } else
#endif
    { // Single-threaded calculation loop
        for (size_t i = 0; i < numMolecules; ++i) {
            const auto& mol = molecules[i];
             if (!mol.isValid()) {
                  // Error value already set during initialization
                  continue;
             }
            for (const auto& name : descriptorNames) {
                 const Descriptor* descriptor = getDescriptor(name);
                 if (!descriptor) continue;
                  try {
                      finalResults[name][i] = descriptor->calculate(mol);
                  } catch (const std::exception& e) {
                       globalLogger.debug("Calculation Error (single-thread) for SMILES " + mol.getOriginalSmiles() + ", Desc: " + name + " - " + e.what());
                   } catch (...) {
                        globalLogger.debug("Unknown Calculation Error (single-thread) for SMILES " + mol.getOriginalSmiles() + ", Desc: " + name);
                   }
            }
        }
    }
    return finalResults;
}
}
