#pragma once

#include "descriptors.hpp"
#include <string>
#include <unordered_map>
#include <cmath>

namespace desfact {
namespace descriptors {

// Base class for string-based descriptors
class StringDescriptor : public Descriptor {
public:
    StringDescriptor(const std::string& name, const std::string& description)
        : Descriptor(name, description) {}

    // Helper methods for common string operations
    static double calcAsciiSum(const std::string& str);
    static double calcAsciiAverage(const std::string& str);
    static double calcAsciiVariance(const std::string& str);
    static double calcBitCount(const std::string& str);
    static double calcEntropy(const std::string& str);
    static double calcCharacterFrequency(const std::string& str, char c);
    static std::unordered_map<char, int> getCharacterHistogram(const std::string& str);
    static int getLongestRun(const std::string& str, char c = 0);
    static double getSymmetryScore(const std::string& str);
    
    // Override to handle string-based calculation
    std::variant<double, int, std::string> calculate(const Molecule& mol) const override {
        if (!mol.isValid()) return 0.0;
        return calculateFromSmiles(mol.getSmiles());
    }
    
    // Each derived class implements this
    virtual std::variant<double, int, std::string> calculateFromSmiles(const std::string& smiles) const = 0;
};

// ASCII and Bit-Level Descriptors
class AsciiSumDescriptor : public StringDescriptor {
public:
    AsciiSumDescriptor() : StringDescriptor("ascii_sum", "Sum of ASCII values of all characters") {}
    std::variant<double, int, std::string> calculateFromSmiles(const std::string& smiles) const override;
};

class AsciiAverageDescriptor : public StringDescriptor {
public:
    AsciiAverageDescriptor() : StringDescriptor("ascii_average", "Average of ASCII values") {}
    std::variant<double, int, std::string> calculateFromSmiles(const std::string& smiles) const override;
};

class AsciiVarianceDescriptor : public StringDescriptor {
public:
    AsciiVarianceDescriptor() : StringDescriptor("ascii_variance", "Variance of ASCII values") {}
    std::variant<double, int, std::string> calculateFromSmiles(const std::string& smiles) const override;
};

class AsciiRangeDescriptor : public StringDescriptor {
public:
    AsciiRangeDescriptor() : StringDescriptor("ascii_range", "Range between highest and lowest ASCII values") {}
    std::variant<double, int, std::string> calculateFromSmiles(const std::string& smiles) const override;
};

class BitCountDescriptor : public StringDescriptor {
public:
    BitCountDescriptor() : StringDescriptor("bit_count", "Count of set bits in byte representation") {}
    std::variant<double, int, std::string> calculateFromSmiles(const std::string& smiles) const override;
};

class BitDensityDescriptor : public StringDescriptor {
public:
    BitDensityDescriptor() : StringDescriptor("bit_density", "Density of set bits (bit count / total bits)") {}
    std::variant<double, int, std::string> calculateFromSmiles(const std::string& smiles) const override;
};

// class BitEntropyDescriptor : public StringDescriptor {
// public:
//     BitEntropyDescriptor() : StringDescriptor("bit_entropy", "Shannon entropy of bit representation") {}
//     std::variant<double, int, std::string> calculateFromSmiles(const std::string& smiles) const override;
// };

class BitTransitionRateDescriptor : public StringDescriptor {
public:
    BitTransitionRateDescriptor() : StringDescriptor("bit_transition_rate", "Rate of 0->1 and 1->0 transitions") {}
    std::variant<double, int, std::string> calculateFromSmiles(const std::string& smiles) const override;
};

class OddBitRatioDescriptor : public StringDescriptor {
public:
    OddBitRatioDescriptor() : StringDescriptor("odd_bit_ratio", "Ratio of bits in odd positions that are set") {}
    std::variant<double, int, std::string> calculateFromSmiles(const std::string& smiles) const override;
};

// Compression & Entropy Descriptors
class RunLengthEncodingSizeDescriptor : public StringDescriptor {
public:
    RunLengthEncodingSizeDescriptor() : StringDescriptor("run_length_encoding_size", "Size after simple run-length encoding") {}
    std::variant<double, int, std::string> calculateFromSmiles(const std::string& smiles) const override;
};

class CharacterRepetitionScoreDescriptor : public StringDescriptor {
public:
    CharacterRepetitionScoreDescriptor() : StringDescriptor("character_repetition_score", "Score based on character repetition patterns") {}
    std::variant<double, int, std::string> calculateFromSmiles(const std::string& smiles) const override;
};

class NibblePatternCountDescriptor : public StringDescriptor {
public:
    NibblePatternCountDescriptor() : StringDescriptor("nibble_pattern_count", "Count of unique 4-bit patterns") {}
    std::variant<double, int, std::string> calculateFromSmiles(const std::string& smiles) const override;
};

class BytePatternCountDescriptor : public StringDescriptor {
public:
    BytePatternCountDescriptor() : StringDescriptor("byte_pattern_count", "Count of unique byte patterns") {}
    std::variant<double, int, std::string> calculateFromSmiles(const std::string& smiles) const override;
};

class EntropyPerByteDescriptor : public StringDescriptor {
public:
    EntropyPerByteDescriptor() : StringDescriptor("entropy_per_byte", "Shannon entropy normalized by string length") {}
    std::variant<double, int, std::string> calculateFromSmiles(const std::string& smiles) const override;
};

class AsciiHistogramUniformityDescriptor : public StringDescriptor {
public:
    AsciiHistogramUniformityDescriptor() : StringDescriptor("ascii_histogram_uniformity", "Uniformity of ASCII character distribution") {}
    std::variant<double, int, std::string> calculateFromSmiles(const std::string& smiles) const override;
};

// Character Sequence Descriptors
class LongestCharacterRunDescriptor : public StringDescriptor {
public:
    LongestCharacterRunDescriptor() : StringDescriptor("longest_character_run", "Length of longest run of same character") {}
    std::variant<double, int, std::string> calculateFromSmiles(const std::string& smiles) const override;
};

class SpecialCharacterDensityDescriptor : public StringDescriptor {
public:
    SpecialCharacterDensityDescriptor() : StringDescriptor("special_character_density", "Density of non-alphanumeric characters") {}
    std::variant<double, int, std::string> calculateFromSmiles(const std::string& smiles) const override;
};

class CaseTransitionRateDescriptor : public StringDescriptor {
public:
    CaseTransitionRateDescriptor() : StringDescriptor("case_transition_rate", "Rate of transitions between uppercase and lowercase") {}
    std::variant<double, int, std::string> calculateFromSmiles(const std::string& smiles) const override;
};

class AlternatingCharacterScoreDescriptor : public StringDescriptor {
public:
    AlternatingCharacterScoreDescriptor() : StringDescriptor("alternating_character_score", "Score for alternating character patterns") {}
    std::variant<double, int, std::string> calculateFromSmiles(const std::string& smiles) const override;
};

class CharacterTrigramCountDescriptor : public StringDescriptor {
public:
    CharacterTrigramCountDescriptor() : StringDescriptor("character_trigram_count", "Count of unique 3-character sequences") {}
    std::variant<double, int, std::string> calculateFromSmiles(const std::string& smiles) const override;
};

class LetterDigitRatioDescriptor : public StringDescriptor {
public:
    LetterDigitRatioDescriptor() : StringDescriptor("letter_digit_ratio", "Ratio of letters to digits") {}
    std::variant<double, int, std::string> calculateFromSmiles(const std::string& smiles) const override;
};

class CharacterBigramEntropyDescriptor : public StringDescriptor {
public:
    CharacterBigramEntropyDescriptor() : StringDescriptor("character_bigram_entropy", "Entropy of character bigrams") {}
    std::variant<double, int, std::string> calculateFromSmiles(const std::string& smiles) const override;
};

// Bit Patterns Descriptors
class BytePalindromeScoreDescriptor : public StringDescriptor {
public:
    BytePalindromeScoreDescriptor() : StringDescriptor("byte_palindrome_score", "Score for palindromic byte patterns") {}
    std::variant<double, int, std::string> calculateFromSmiles(const std::string& smiles) const override;
};

class BitPalindromeScoreDescriptor : public StringDescriptor {
public:
    BitPalindromeScoreDescriptor() : StringDescriptor("bit_palindrome_score", "Score for palindromic bit patterns") {}
    std::variant<double, int, std::string> calculateFromSmiles(const std::string& smiles) const override;
};

// class ByteMirrorSymmetryDescriptor : public StringDescriptor {
// public:
//     ByteMirrorSymmetryDescriptor() : StringDescriptor("byte_mirror_symmetry", "Symmetry score for byte patterns") {}
//     std::variant<double, int, std::string> calculateFromSmiles(const std::string& smiles) const override;
// };

class BitPeriodicityDescriptor : public StringDescriptor {
public:
    BitPeriodicityDescriptor() : StringDescriptor("bit_periodicity", "Periodicity in bit patterns") {}
    std::variant<double, int, std::string> calculateFromSmiles(const std::string& smiles) const override;
};

class BitAutocorrelationDescriptor : public StringDescriptor {
public:
    BitAutocorrelationDescriptor() : StringDescriptor("bit_autocorrelation", "Autocorrelation of bit patterns") {}
    std::variant<double, int, std::string> calculateFromSmiles(const std::string& smiles) const override;
};

class HammingWeightDistributionDescriptor : public StringDescriptor {
public:
    HammingWeightDistributionDescriptor() : StringDescriptor("hamming_weight_distribution", "Distribution of Hamming weights") {}
    std::variant<double, int, std::string> calculateFromSmiles(const std::string& smiles) const override;
};

class ByteParityDistributionDescriptor : public StringDescriptor {
public:
    ByteParityDistributionDescriptor() : StringDescriptor("byte_parity_distribution", "Distribution of byte parities") {}
    std::variant<double, int, std::string> calculateFromSmiles(const std::string& smiles) const override;
};

class BitRunEntropyDescriptor : public StringDescriptor {
public:
    BitRunEntropyDescriptor() : StringDescriptor("bit_run_entropy", "Entropy of bit runs") {}
    std::variant<double, int, std::string> calculateFromSmiles(const std::string& smiles) const override;
};

// Positional Heuristics Descriptors
class PositionWeightDescriptor : public StringDescriptor {
public:
    PositionWeightDescriptor() : StringDescriptor("position_weight", "Weighted sum based on character positions") {}
    std::variant<double, int, std::string> calculateFromSmiles(const std::string& smiles) const override;
};

class FrontBackRatioDescriptor : public StringDescriptor {
public:
    FrontBackRatioDescriptor() : StringDescriptor("front_back_ratio", "Ratio of characteristics front vs back of string") {}
    std::variant<double, int, std::string> calculateFromSmiles(const std::string& smiles) const override;
};

class BracketAsymmetryDescriptor : public StringDescriptor {
public:
    BracketAsymmetryDescriptor() : StringDescriptor("bracket_asymmetry", "Asymmetry in bracket distribution") {}
    std::variant<double, int, std::string> calculateFromSmiles(const std::string& smiles) const override;
};

class CharDispersionDescriptor : public StringDescriptor {
public:
    CharDispersionDescriptor() : StringDescriptor("char_dispersion", "Dispersion of characters across string") {}
    std::variant<double, int, std::string> calculateFromSmiles(const std::string& smiles) const override;
};

class DigitPositionMeanDescriptor : public StringDescriptor {
public:
    DigitPositionMeanDescriptor() : StringDescriptor("digit_position_mean", "Mean position of digits in string") {}
    std::variant<double, int, std::string> calculateFromSmiles(const std::string& smiles) const override;
};

class SymmetryScoreDescriptor : public StringDescriptor {
public:
    SymmetryScoreDescriptor() : StringDescriptor("symmetry_score", "Overall string symmetry score") {}
    std::variant<double, int, std::string> calculateFromSmiles(const std::string& smiles) const override;
};

// Pattern / Fragment Pseudodescriptors
class AromaticRatioDescriptor : public StringDescriptor {
public:
    AromaticRatioDescriptor() : StringDescriptor("aromatic_ratio", "Ratio of lowercase characters (aromatic notation)") {}
    std::variant<double, int, std::string> calculateFromSmiles(const std::string& smiles) const override;
};

class HeteroatomRatioDescriptor : public StringDescriptor {
public:
    HeteroatomRatioDescriptor() : StringDescriptor("heteroatom_ratio", "Ratio of non-C,H characters") {}
    std::variant<double, int, std::string> calculateFromSmiles(const std::string& smiles) const override;
};

class ElementDiversityDescriptor : public StringDescriptor {
public:
    ElementDiversityDescriptor() : StringDescriptor("element_diversity", "Diversity of elements based on parsing") {}
    std::variant<double, int, std::string> calculateFromSmiles(const std::string& smiles) const override;
};

class SaturationIndexDescriptor : public StringDescriptor {
public:
    SaturationIndexDescriptor() : StringDescriptor("saturation_index", "Index based on bond symbols") {}
    std::variant<double, int, std::string> calculateFromSmiles(const std::string& smiles) const override;
};

class PatternEntropyDescriptor : public StringDescriptor {
public:
    PatternEntropyDescriptor() : StringDescriptor("pattern_entropy", "Shannon entropy of SMILES tokens") {}
    std::variant<double, int, std::string> calculateFromSmiles(const std::string& smiles) const override;
};

class BranchComplexityDescriptor : public StringDescriptor {
public:
    BranchComplexityDescriptor() : StringDescriptor("branch_complexity", "Complexity based on branching") {}
    std::variant<double, int, std::string> calculateFromSmiles(const std::string& smiles) const override;
};

// Basic SMILES Structure-Based
// class RingPositionVarianceDescriptor : public StringDescriptor {
// public:
//     RingPositionVarianceDescriptor() : StringDescriptor("ring_position_variance", "Variance in position of ring indices") {}
//     std::variant<double, int, std::string> calculateFromSmiles(const std::string& smiles) const override;
// };

class BondPositionEntropyDescriptor : public StringDescriptor {
public:
    BondPositionEntropyDescriptor() : StringDescriptor("bond_position_entropy", "Entropy of bond symbol positions") {}
    std::variant<double, int, std::string> calculateFromSmiles(const std::string& smiles) const override;
};

class ElementPositionBiasDescriptor : public StringDescriptor {
public:
    ElementPositionBiasDescriptor() : StringDescriptor("element_position_bias", "Bias in element position distribution") {}
    std::variant<double, int, std::string> calculateFromSmiles(const std::string& smiles) const override;
};

class BondComplexityDescriptor : public StringDescriptor {
public:
    BondComplexityDescriptor() : StringDescriptor("bond_complexity", "Complexity based on bond symbols") {}
    std::variant<double, int, std::string> calculateFromSmiles(const std::string& smiles) const override;
};

class RingComplexityDescriptor : public StringDescriptor {
public:
    RingComplexityDescriptor() : StringDescriptor("ring_complexity", "Complexity based on ring closure symbols") {}
    std::variant<double, int, std::string> calculateFromSmiles(const std::string& smiles) const override;
};

// Enhanced Statistical Features
class ChainLengthDescriptor : public StringDescriptor {
public:
    ChainLengthDescriptor() : StringDescriptor("chain_length", "Estimates the length of the main chain") {}
    std::variant<double, int, std::string> calculateFromSmiles(const std::string& smiles) const override;
};

class DeviationSymmetryDescriptor : public StringDescriptor {
public:
    DeviationSymmetryDescriptor() : StringDescriptor("deviation_symmetry", "Measures asymmetry in character distribution") {}
    std::variant<double, int, std::string> calculateFromSmiles(const std::string& smiles) const override;
};

class LocalComplexityDescriptor : public StringDescriptor {
public:
    LocalComplexityDescriptor() : StringDescriptor("local_complexity", "Measures local complexity using sliding windows") {}
    std::variant<double, int, std::string> calculateFromSmiles(const std::string& smiles) const override;
};

class CharacterKurtosisDescriptor : public StringDescriptor {
public:
    CharacterKurtosisDescriptor() : StringDescriptor("character_kurtosis", "Measures peakedness of character distribution") {}
    std::variant<double, int, std::string> calculateFromSmiles(const std::string& smiles) const override;
};

class CharacterSkewnessDescriptor : public StringDescriptor {
public:
    CharacterSkewnessDescriptor() : StringDescriptor("character_skewness", "Measures asymmetry of character distribution") {}
    std::variant<double, int, std::string> calculateFromSmiles(const std::string& smiles) const override;
};

// Sequence Analysis Features
class SequenceNoisinessDescriptor : public StringDescriptor {
public:
    SequenceNoisinessDescriptor() : StringDescriptor("sequence_noisiness", "Measures randomness in character sequence") {}
    std::variant<double, int, std::string> calculateFromSmiles(const std::string& smiles) const override;
};

class CharacterTransitionMatrixEntropyDescriptor : public StringDescriptor {
public:
    CharacterTransitionMatrixEntropyDescriptor() : StringDescriptor("char_transition_matrix_entropy", "Entropy of transition matrix between characters") {}
    std::variant<double, int, std::string> calculateFromSmiles(const std::string& smiles) const override;
};

class SequenceFractalDimensionDescriptor : public StringDescriptor {
public:
    SequenceFractalDimensionDescriptor() : StringDescriptor("sequence_fractal_dimension", "Estimates fractal dimension using box-counting") {}
    std::variant<double, int, std::string> calculateFromSmiles(const std::string& smiles) const override;
};

class SubsequenceRepetitionDescriptor : public StringDescriptor {
public:
    SubsequenceRepetitionDescriptor() : StringDescriptor("subsequence_repetition", "Measures repetition of subsequences") {}
    std::variant<double, int, std::string> calculateFromSmiles(const std::string& smiles) const override;
};

class MaxRepeatedSubstringDescriptor : public StringDescriptor {
public:
    MaxRepeatedSubstringDescriptor() : StringDescriptor("max_repeated_substring", "Length of the longest repeated substring") {}
    std::variant<double, int, std::string> calculateFromSmiles(const std::string& smiles) const override;
};

// Fragment-Based Features
// class FragmentSizeVarianceDescriptor : public StringDescriptor {
// public:
//     FragmentSizeVarianceDescriptor() : StringDescriptor("fragment_size_variance", "Variance in fragment sizes separated by '.'") {}
//     std::variant<double, int, std::string> calculateFromSmiles(const std::string& smiles) const override;
// };

class NonCarbonComplexityDescriptor : public StringDescriptor {
public:
    NonCarbonComplexityDescriptor() : StringDescriptor("non_carbon_complexity", "Complexity measure focusing on heteroatoms") {}
    std::variant<double, int, std::string> calculateFromSmiles(const std::string& smiles) const override;
};

class SideChainIndexDescriptor : public StringDescriptor {
public:
    SideChainIndexDescriptor() : StringDescriptor("side_chain_index", "Measures prevalence of side chains") {}
    std::variant<double, int, std::string> calculateFromSmiles(const std::string& smiles) const override;
};

class RingBranchRatioDescriptor : public StringDescriptor {
public:
    RingBranchRatioDescriptor() : StringDescriptor("ring_branch_ratio", "Ratio of ring systems to branch systems") {}
    std::variant<double, int, std::string> calculateFromSmiles(const std::string& smiles) const override;
};

class AromaticAliphaticRatioDescriptor : public StringDescriptor {
public:
    AromaticAliphaticRatioDescriptor() : StringDescriptor("aromatic_aliphatic_ratio", "Ratio of aromatic to aliphatic fragments") {}
    std::variant<double, int, std::string> calculateFromSmiles(const std::string& smiles) const override;
};

// SMILES Syntax-Specific Features
// class DisjointFragmentCountDescriptor : public StringDescriptor {
// public:
//     DisjointFragmentCountDescriptor() : StringDescriptor("disjoint_fragment_count", "Count of disjoint fragments (separated by '.')") {}
//     std::variant<double, int, std::string> calculateFromSmiles(const std::string& smiles) const override;
// };

// class StereoDescriptorRatioDescriptor : public StringDescriptor {
// public:
//     StereoDescriptorRatioDescriptor() : StringDescriptor("stereo_descriptor_ratio", "Ratio of stereo descriptors") {}
//     std::variant<double, int, std::string> calculateFromSmiles(const std::string& smiles) const override;
// };

// Specialized Chemical Pattern Features
class FunctionalGroupEntropyDescriptor : public StringDescriptor {
public:
    FunctionalGroupEntropyDescriptor() : StringDescriptor("functional_group_entropy", "Entropy of simple functional group patterns") {}
    std::variant<double, int, std::string> calculateFromSmiles(const std::string& smiles) const override;
};

class LargestRingDescriptor : public StringDescriptor {
public:
    LargestRingDescriptor() : StringDescriptor("largest_ring", "Estimate of the largest ring size") {}
    std::variant<double, int, std::string> calculateFromSmiles(const std::string& smiles) const override;
};

class CarbonSkeletonComplexityDescriptor : public StringDescriptor {
public:
    CarbonSkeletonComplexityDescriptor() : StringDescriptor("carbon_skeleton_complexity", "Complexity of carbon skeleton") {}
    std::variant<double, int, std::string> calculateFromSmiles(const std::string& smiles) const override;
};

class HeterocycleDescriptor : public StringDescriptor {
public:
    HeterocycleDescriptor() : StringDescriptor("heterocycle", "Estimate of heterocyclic structures") {}
    std::variant<double, int, std::string> calculateFromSmiles(const std::string& smiles) const override;
};

class SubstructureBalanceDescriptor : public StringDescriptor {
public:
    SubstructureBalanceDescriptor() : StringDescriptor("substructure_balance", "Balance between different substructures") {}
    std::variant<double, int, std::string> calculateFromSmiles(const std::string& smiles) const override;
};

// Advanced Entropy/Information Features
class MarkovEntropyDescriptor : public StringDescriptor {
public:
    MarkovEntropyDescriptor() : StringDescriptor("markov_entropy", "Entropy estimated using Markov models") {}
    std::variant<double, int, std::string> calculateFromSmiles(const std::string& smiles) const override;
};

class ContextualEntropyDescriptor : public StringDescriptor {
public:
    ContextualEntropyDescriptor() : StringDescriptor("contextual_entropy", "Entropy based on preceding context") {}
    std::variant<double, int, std::string> calculateFromSmiles(const std::string& smiles) const override;
};

class SyntacticComplexityDescriptor : public StringDescriptor {
public:
    SyntacticComplexityDescriptor() : StringDescriptor("syntactic_complexity", "Complexity based on SMILES syntax rules") {}
    std::variant<double, int, std::string> calculateFromSmiles(const std::string& smiles) const override;
};

class RelativeEntropyDescriptor : public StringDescriptor {
public:
    RelativeEntropyDescriptor() : StringDescriptor("relative_entropy", "KL divergence from typical SMILES distribution") {}
    std::variant<double, int, std::string> calculateFromSmiles(const std::string& smiles) const override;
};

class InformationDensityDescriptor : public StringDescriptor {
public:
    InformationDensityDescriptor() : StringDescriptor("information_density", "Information content per character") {}
    std::variant<double, int, std::string> calculateFromSmiles(const std::string& smiles) const override;
};

} // namespace descriptors
} // namespace desfact 