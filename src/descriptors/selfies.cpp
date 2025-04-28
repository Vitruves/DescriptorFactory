#include "descriptors/selfies.hpp"
#include "utils.hpp"
#include <algorithm>
#include <cmath>
#include <sstream>
#include <unordered_map>
#include <numeric>
#include <unordered_set>

namespace desfact {
namespace descriptors {

SelfiesDescriptor::SelfiesDescriptor(const std::string& name, const std::string& description)
    : Descriptor(name, description) {}

std::vector<std::string> SelfiesDescriptor::tokenizeSelfies(const std::string& selfiesStr) {
    std::vector<std::string> tokens;
    std::string token;
    bool inToken = false;
    
    for (char c : selfiesStr) {
        if (c == '[') {
            if (inToken) {
                token += c;
            } else {
                inToken = true;
                token = "[";
            }
        } 
        else if (c == ']') {
            if (inToken) {
                token += c;
                tokens.push_back(token);
                token.clear();
                inToken = false;
            }
        }
        else if (inToken) {
            token += c;
        }
    }
    
    return tokens;
}

bool SelfiesDescriptor::isBranchToken(const std::string& token) {
    return token.find("[Branch") == 0;
}

bool SelfiesDescriptor::isRingToken(const std::string& token) {
    return token.find("[Ring") == 0;
}

bool SelfiesDescriptor::isAtomToken(const std::string& token) {
    if (isBranchToken(token) || isRingToken(token)) return false;
    
    // Basic check for atom tokens - could be enhanced with more specific rules
    static const std::vector<std::string> atomPrefixes = {
        "[C", "[N", "[O", "[S", "[P", "[F", "[Cl", "[Br", "[I", "[B"
    };
    
    for (const auto& prefix : atomPrefixes) {
        if (token.find(prefix) == 0) return true;
    }
    
    return false;
}

double SelfiesDescriptor::calculateTokenEntropy(const std::vector<std::string>& tokens) {
    if (tokens.empty()) return 0.0;
    
    std::unordered_map<std::string, int> counts;
    for (const auto& token : tokens) counts[token]++;
    
    double entropy = 0.0;
    for (const auto& [token, count] : counts) {
        double p = static_cast<double>(count) / tokens.size();
        entropy -= p * std::log2(p);
    }
    
    return entropy;
}

int SelfiesDescriptor::calculateBranchingDepth(const std::vector<std::string>& tokens) {
    int maxDepth = 0;
    int currentDepth = 0;
    
    for (const auto& token : tokens) {
        if (isBranchToken(token)) {
            currentDepth++;
            maxDepth = std::max(maxDepth, currentDepth);
        } else if (token.find("[=Branch") == 0) { // Branch closing token
            currentDepth = std::max(0, currentDepth - 1);
        }
    }
    
    return maxDepth;
}

std::variant<double, int, std::string> SelfiesDescriptor::calculate(const Molecule& mol) const {
    if (!mol.isValid()) {
        globalLogger.debug(getName() + ": Invalid molecule input.");
        return 0;
    }
    
    // Assuming the original SMILES is actually SELFIES format
    // In a real implementation, you'd need to ensure SELFIES input or conversion
    std::string selfiesStr = mol.getOriginalSmiles();
    if (selfiesStr.empty()) {
        globalLogger.debug(getName() + ": Empty SELFIES string.");
        return 0;
    }
    
    return calculateFromSelfies(selfiesStr);
}

// 1. Raw SELFIES token count
SelfiesTokenCountDescriptor::SelfiesTokenCountDescriptor()
    : SelfiesDescriptor("SelfiesTokenCount", "Raw count of tokens in the SELFIES string") {}

std::variant<double, int, std::string> SelfiesTokenCountDescriptor::calculateFromSelfies(const std::string& selfiesStr) const {
    return static_cast<int>(tokenizeSelfies(selfiesStr).size());
}

// 2. SELFIES token type distribution
SelfiesTokenDistributionDescriptor::SelfiesTokenDistributionDescriptor()
    : SelfiesDescriptor("SelfiesTokenDistribution", "Distribution of different token types in SELFIES") {}

std::variant<double, int, std::string> SelfiesTokenDistributionDescriptor::calculateFromSelfies(const std::string& selfiesStr) const {
    auto tokens = tokenizeSelfies(selfiesStr);
    std::unordered_map<std::string, int> distribution;
    
    for (const auto& token : tokens) distribution[token]++;
    
    std::stringstream ss;
    for (const auto& [token, count] : distribution) {
        ss << token << ":" << count << ";";
    }
    return ss.str();
}

// 3. Average token complexity
SelfiesAverageComplexityDescriptor::SelfiesAverageComplexityDescriptor()
    : SelfiesDescriptor("SelfiesAvgComplexity", "Average complexity of SELFIES tokens") {}

std::variant<double, int, std::string> SelfiesAverageComplexityDescriptor::calculateFromSelfies(const std::string& selfiesStr) const {
    auto tokens = tokenizeSelfies(selfiesStr);
    if (tokens.empty()) return 0.0;
    
    double totalComplexity = 0.0;
    for (const auto& token : tokens) {
        totalComplexity += token.length();
    }
    
    return totalComplexity / tokens.size();
}

// 4. Branching depth
SelfiesBranchingDepthDescriptor::SelfiesBranchingDepthDescriptor()
    : SelfiesDescriptor("SelfiesBranchingDepth", "Maximum branching depth in SELFIES") {}

std::variant<double, int, std::string> SelfiesBranchingDepthDescriptor::calculateFromSelfies(const std::string& selfiesStr) const {
    return calculateBranchingDepth(tokenizeSelfies(selfiesStr));
}

// 5. Maximum branching fan-out
SelfiesMaxBranchingFanoutDescriptor::SelfiesMaxBranchingFanoutDescriptor()
    : SelfiesDescriptor("SelfiesMaxBranchFanout", "Maximum branches from a single point") {}

std::variant<double, int, std::string> SelfiesMaxBranchingFanoutDescriptor::calculateFromSelfies(const std::string& selfiesStr) const {
    auto tokens = tokenizeSelfies(selfiesStr);
    int maxFanout = 0;
    
    for (size_t i = 0; i < tokens.size(); i++) {
        if (isAtomToken(tokens[i])) {
            int currentFanout = 0;
            for (size_t j = i + 1; j < tokens.size() && isBranchToken(tokens[j]); j++) {
                currentFanout++;
            }
            maxFanout = std::max(maxFanout, currentFanout);
        }
    }
    
    return maxFanout;
}

// 6. Empty branches
SelfiesEmptyBranchesDescriptor::SelfiesEmptyBranchesDescriptor()
    : SelfiesDescriptor("SelfiesEmptyBranches", "Count of empty branches in SELFIES") {}

std::variant<double, int, std::string> SelfiesEmptyBranchesDescriptor::calculateFromSelfies(const std::string& selfiesStr) const {
    auto tokens = tokenizeSelfies(selfiesStr);
    int emptyBranches = 0;
    
    for (const auto& token : tokens) {
        if (token.find("[Branch") == 0 && token.find("_0]") != std::string::npos) {
            emptyBranches++;
        }
    }
    
    return emptyBranches;
}

// 7. Ring closure token frequency
SelfiesRingTokenFrequencyDescriptor::SelfiesRingTokenFrequencyDescriptor()
    : SelfiesDescriptor("SelfiesRingTokenFreq", "Frequency of ring tokens in SELFIES") {}

std::variant<double, int, std::string> SelfiesRingTokenFrequencyDescriptor::calculateFromSelfies(const std::string& selfiesStr) const {
    auto tokens = tokenizeSelfies(selfiesStr);
    if (tokens.empty()) return 0.0;
    
    int ringTokens = 0;
    for (const auto& token : tokens) {
        if (isRingToken(token)) ringTokens++;
    }
    
    return static_cast<double>(ringTokens) / tokens.size();
}

// 8. Token entropy
SelfiesTokenEntropyDescriptor::SelfiesTokenEntropyDescriptor()
    : SelfiesDescriptor("SelfiesTokenEntropy", "Shannon entropy of SELFIES tokens") {}

std::variant<double, int, std::string> SelfiesTokenEntropyDescriptor::calculateFromSelfies(const std::string& selfiesStr) const {
    return calculateTokenEntropy(tokenizeSelfies(selfiesStr));
}

// 9. SELFIES string length
SelfiesStringLengthDescriptor::SelfiesStringLengthDescriptor()
    : SelfiesDescriptor("SelfiesStringLength", "Raw character length of SELFIES string") {}

std::variant<double, int, std::string> SelfiesStringLengthDescriptor::calculateFromSelfies(const std::string& selfiesStr) const {
    return static_cast<int>(selfiesStr.length());
}

// 10. Branch token ratio
SelfiesBranchTokenRatioDescriptor::SelfiesBranchTokenRatioDescriptor()
    : SelfiesDescriptor("SelfiesBranchTokenRatio", "Ratio of branch tokens to all tokens") {}

std::variant<double, int, std::string> SelfiesBranchTokenRatioDescriptor::calculateFromSelfies(const std::string& selfiesStr) const {
    auto tokens = tokenizeSelfies(selfiesStr);
    if (tokens.empty()) return 0.0;
    
    int branchTokens = 0;
    for (const auto& token : tokens) {
        if (isBranchToken(token)) branchTokens++;
    }
    
    return static_cast<double>(branchTokens) / tokens.size();
}

bool SelfiesDescriptor::isAromaticToken(const std::string& token) {
    static const std::vector<std::string> aromaticPrefixes = {
        "[c", "[n", "[o", "[s", "[p", "[as", "[se"
    };
    
    for (const auto& prefix : aromaticPrefixes) {
        if (token.find(prefix) == 0) return true;
    }
    
    return false;
}

bool SelfiesDescriptor::isChargeToken(const std::string& token) {
    return token.find("+") != std::string::npos || token.find("-") != std::string::npos;
}

bool SelfiesDescriptor::isHeteroatomToken(const std::string& token) {
    if (!isAtomToken(token)) return false;
    
    static const std::vector<std::string> carbonPrefixes = {
        "[C", "[c"
    };
    
    for (const auto& prefix : carbonPrefixes) {
        if (token.find(prefix) == 0) return false;
    }
    
    return true;
}

bool SelfiesDescriptor::isStereoToken(const std::string& token) {
    return token.find("/") != std::string::npos || 
           token.find("\\") != std::string::npos ||
           token.find("@") != std::string::npos;
}

bool SelfiesDescriptor::isUnpairedElectronToken(const std::string& token) {
    // Tokens with unpaired electrons (radicals)
    return token.find(".") != std::string::npos;
}

int SelfiesDescriptor::getRingSizeFromToken(const std::string& token) {
    if (isRingToken(token)) {
        // Extract ring size from [RingX] token where X is the size
        size_t startPos = token.find("Ring") + 4;
        size_t endPos = token.find("]");
        
        if (startPos < token.size() && endPos != std::string::npos) {
            try {
                return std::stoi(token.substr(startPos, endPos - startPos));
            } catch (...) {
                return 0;
            }
        }
    }
    return 0;
}

double SelfiesDescriptor::calculateSequenceSymmetry(const std::vector<std::string>& tokens) {
    if (tokens.size() <= 1) return 1.0; // Perfect symmetry for empty or single token
    
    size_t midpoint = tokens.size() / 2;
    int matchingTokens = 0;
    
    for (size_t i = 0; i < midpoint; ++i) {
        if (tokens[i] == tokens[tokens.size() - 1 - i]) {
            matchingTokens++;
        }
    }
    
    return static_cast<double>(matchingTokens) / midpoint;
}

int SelfiesDescriptor::calculateLevenshteinDistance(const std::string& s1, const std::string& s2) {
    const size_t len1 = s1.size();
    const size_t len2 = s2.size();
    
    std::vector<std::vector<int>> dp(len1 + 1, std::vector<int>(len2 + 1));
    
    for (size_t i = 0; i <= len1; i++) dp[i][0] = i;
    for (size_t j = 0; j <= len2; j++) dp[0][j] = j;
    
    for (size_t i = 1; i <= len1; i++) {
        for (size_t j = 1; j <= len2; j++) {
            int cost = (s1[i - 1] == s2[j - 1]) ? 0 : 1;
            dp[i][j] = std::min({
                dp[i - 1][j] + 1,      // deletion
                dp[i][j - 1] + 1,      // insertion
                dp[i - 1][j - 1] + cost // substitution
            });
        }
    }
    
    return dp[len1][len2];
}

// 1. Ring token ratio
SelfiesRingTokenRatioDescriptor::SelfiesRingTokenRatioDescriptor()
    : SelfiesDescriptor("SelfiesRingTokenRatio", "Ratio of ring tokens to all tokens") {}

std::variant<double, int, std::string> SelfiesRingTokenRatioDescriptor::calculateFromSelfies(const std::string& selfiesStr) const {
    auto tokens = tokenizeSelfies(selfiesStr);
    if (tokens.empty()) return 0.0;
    
    int ringTokens = 0;
    for (const auto& token : tokens) {
        if (isRingToken(token)) ringTokens++;
    }
    
    return static_cast<double>(ringTokens) / tokens.size();
}

// 2. Bigram count
SelfiesBigramCountDescriptor::SelfiesBigramCountDescriptor()
    : SelfiesDescriptor("SelfiesBigramCount", "Count of unique token bigrams in SELFIES") {}

std::variant<double, int, std::string> SelfiesBigramCountDescriptor::calculateFromSelfies(const std::string& selfiesStr) const {
    auto tokens = tokenizeSelfies(selfiesStr);
    if (tokens.size() < 2) return 0;
    
    std::unordered_map<std::string, int> bigrams;
    for (size_t i = 0; i < tokens.size() - 1; ++i) {
        bigrams[tokens[i] + "|" + tokens[i+1]]++;
    }
    
    return static_cast<int>(bigrams.size());
}

// 3. Trigram count
SelfiesTrigramCountDescriptor::SelfiesTrigramCountDescriptor()
    : SelfiesDescriptor("SelfiesTrigramCount", "Count of unique token trigrams in SELFIES") {}

std::variant<double, int, std::string> SelfiesTrigramCountDescriptor::calculateFromSelfies(const std::string& selfiesStr) const {
    auto tokens = tokenizeSelfies(selfiesStr);
    if (tokens.size() < 3) return 0;
    
    std::unordered_map<std::string, int> trigrams;
    for (size_t i = 0; i < tokens.size() - 2; ++i) {
        trigrams[tokens[i] + "|" + tokens[i+1] + "|" + tokens[i+2]]++;
    }
    
    return static_cast<int>(trigrams.size());
}

// 4. Token repetition rate
SelfiesTokenRepetitionRateDescriptor::SelfiesTokenRepetitionRateDescriptor()
    : SelfiesDescriptor("SelfiesTokenRepetitionRate", "Rate at which tokens repeat in SELFIES") {}

std::variant<double, int, std::string> SelfiesTokenRepetitionRateDescriptor::calculateFromSelfies(const std::string& selfiesStr) const {
    auto tokens = tokenizeSelfies(selfiesStr);
    if (tokens.empty()) return 0.0;
    
    std::unordered_map<std::string, int> tokenCounts;
    for (const auto& token : tokens) {
        tokenCounts[token]++;
    }
    
    int repeatedTokens = 0;
    for (const auto& [token, count] : tokenCounts) {
        if (count > 1) repeatedTokens += count;
    }
    
    return static_cast<double>(repeatedTokens) / tokens.size();
}

// 5. Aromatic token ratio
SelfiesAromaticTokenRatioDescriptor::SelfiesAromaticTokenRatioDescriptor()
    : SelfiesDescriptor("SelfiesAromaticTokenRatio", "Ratio of aromatic tokens in SELFIES") {}

std::variant<double, int, std::string> SelfiesAromaticTokenRatioDescriptor::calculateFromSelfies(const std::string& selfiesStr) const {
    auto tokens = tokenizeSelfies(selfiesStr);
    if (tokens.empty()) return 0.0;
    
    int aromaticTokens = 0;
    int atomTokens = 0;
    
    for (const auto& token : tokens) {
        if (isAtomToken(token)) {
            atomTokens++;
            if (isAromaticToken(token)) aromaticTokens++;
        }
    }
    
    return atomTokens > 0 ? static_cast<double>(aromaticTokens) / atomTokens : 0.0;
}

// 6. Charge token ratio
SelfiesChargeTokenRatioDescriptor::SelfiesChargeTokenRatioDescriptor()
    : SelfiesDescriptor("SelfiesChargeTokenRatio", "Ratio of tokens with formal charges in SELFIES") {}

std::variant<double, int, std::string> SelfiesChargeTokenRatioDescriptor::calculateFromSelfies(const std::string& selfiesStr) const {
    auto tokens = tokenizeSelfies(selfiesStr);
    if (tokens.empty()) return 0.0;
    
    int chargeTokens = 0;
    for (const auto& token : tokens) {
        if (isChargeToken(token)) chargeTokens++;
    }
    
    return static_cast<double>(chargeTokens) / tokens.size();
}

// 7. Aliphatic chain length
SelfiesAliphaticChainLengthDescriptor::SelfiesAliphaticChainLengthDescriptor()
    : SelfiesDescriptor("SelfiesAliphaticChainLength", "Longest run of aliphatic atom tokens") {}

std::variant<double, int, std::string> SelfiesAliphaticChainLengthDescriptor::calculateFromSelfies(const std::string& selfiesStr) const {
    auto tokens = tokenizeSelfies(selfiesStr);
    int maxChainLength = 0;
    int currentChainLength = 0;
    
    for (const auto& token : tokens) {
        if (isAtomToken(token) && !isAromaticToken(token) && !isHeteroatomToken(token)) {
            currentChainLength++;
            maxChainLength = std::max(maxChainLength, currentChainLength);
        } else if (isBranchToken(token) || isRingToken(token)) {
            // Continue chain through branch/ring tokens
            continue;
        } else {
            currentChainLength = 0;
        }
    }
    
    return maxChainLength;
}

// 8. Heteroatom ratio
SelfiesHeteroatomRatioDescriptor::SelfiesHeteroatomRatioDescriptor()
    : SelfiesDescriptor("SelfiesHeteroatomRatio", "Ratio of non-carbon atom tokens") {}

std::variant<double, int, std::string> SelfiesHeteroatomRatioDescriptor::calculateFromSelfies(const std::string& selfiesStr) const {
    auto tokens = tokenizeSelfies(selfiesStr);
    int heteroatomTokens = 0;
    int atomTokens = 0;
    
    for (const auto& token : tokens) {
        if (isAtomToken(token)) {
            atomTokens++;
            if (isHeteroatomToken(token)) {
                heteroatomTokens++;
            }
        }
    }
    
    return atomTokens > 0 ? static_cast<double>(heteroatomTokens) / atomTokens : 0.0;
}

// 9. Maximum ring size
SelfiesMaxRingSizeDescriptor::SelfiesMaxRingSizeDescriptor()
    : SelfiesDescriptor("SelfiesMaxRingSize", "Maximum ring size from SELFIES ring tokens") {}

std::variant<double, int, std::string> SelfiesMaxRingSizeDescriptor::calculateFromSelfies(const std::string& selfiesStr) const {
    auto tokens = tokenizeSelfies(selfiesStr);
    int maxRingSize = 0;
    
    for (const auto& token : tokens) {
        if (isRingToken(token)) {
            int ringSize = getRingSizeFromToken(token);
            maxRingSize = std::max(maxRingSize, ringSize);
        }
    }
    
    return maxRingSize;
}

// 10. Stereo token ratio
SelfiesStereoTokenRatioDescriptor::SelfiesStereoTokenRatioDescriptor()
    : SelfiesDescriptor("SelfiesStereoTokenRatio", "Ratio of tokens with stereochemistry indicators") {}

std::variant<double, int, std::string> SelfiesStereoTokenRatioDescriptor::calculateFromSelfies(const std::string& selfiesStr) const {
    auto tokens = tokenizeSelfies(selfiesStr);
    if (tokens.empty()) return 0.0;
    
    int stereoTokens = 0;
    for (const auto& token : tokens) {
        if (isStereoToken(token)) {
            stereoTokens++;
        }
    }
    
    return static_cast<double>(stereoTokens) / tokens.size();
}

// 11. Grammar depth complexity
SelfiesGrammarDepthComplexityDescriptor::SelfiesGrammarDepthComplexityDescriptor()
    : SelfiesDescriptor("SelfiesGrammarDepthComplexity", "Grammar nesting depth score of SELFIES") {}

std::variant<double, int, std::string> SelfiesGrammarDepthComplexityDescriptor::calculateFromSelfies(const std::string& selfiesStr) const {
    auto tokens = tokenizeSelfies(selfiesStr);
    int branchDepth = calculateBranchingDepth(tokens);
    
    // Average depth score weighted by token count
    return static_cast<double>(branchDepth) * std::sqrt(tokens.size()) / 10.0;
}

// 12. Vocabulary diversity
SelfiesVocabularyDiversityDescriptor::SelfiesVocabularyDiversityDescriptor()
    : SelfiesDescriptor("SelfiesVocabularyDiversity", "Unique tokens divided by total tokens") {}

std::variant<double, int, std::string> SelfiesVocabularyDiversityDescriptor::calculateFromSelfies(const std::string& selfiesStr) const {
    auto tokens = tokenizeSelfies(selfiesStr);
    if (tokens.empty()) return 0.0;
    
    std::unordered_set<std::string> uniqueTokens(tokens.begin(), tokens.end());
    return static_cast<double>(uniqueTokens.size()) / tokens.size();
}

// 13. Bond complexity
SelfiesBondComplexityDescriptor::SelfiesBondComplexityDescriptor()
    : SelfiesDescriptor("SelfiesBondComplexity", "Score based on bond types in tokens") {}

std::variant<double, int, std::string> SelfiesBondComplexityDescriptor::calculateFromSelfies(const std::string& selfiesStr) const {
    auto tokens = tokenizeSelfies(selfiesStr);
    double bondComplexity = 0.0;
    
    for (const auto& token : tokens) {
        if (token.find("=") != std::string::npos) bondComplexity += 2.0;
        else if (token.find("#") != std::string::npos) bondComplexity += 3.0;
        else if (token.find(":") != std::string::npos) bondComplexity += 1.5;
        else if (isAtomToken(token)) bondComplexity += 1.0;
    }
    
    return bondComplexity;
}

// 14. Sequence symmetry
SelfiesSequenceSymmetryDescriptor::SelfiesSequenceSymmetryDescriptor()
    : SelfiesDescriptor("SelfiesSequenceSymmetry", "Symmetry score of token sequences") {}

std::variant<double, int, std::string> SelfiesSequenceSymmetryDescriptor::calculateFromSelfies(const std::string& selfiesStr) const {
    auto tokens = tokenizeSelfies(selfiesStr);
    return calculateSequenceSymmetry(tokens);
}

// 15. Token length variance
SelfiesTokenLengthVarianceDescriptor::SelfiesTokenLengthVarianceDescriptor()
    : SelfiesDescriptor("SelfiesTokenLengthVariance", "Variance in character length of tokens") {}

std::variant<double, int, std::string> SelfiesTokenLengthVarianceDescriptor::calculateFromSelfies(const std::string& selfiesStr) const {
    auto tokens = tokenizeSelfies(selfiesStr);
    if (tokens.empty()) return 0.0;
    
    std::vector<double> lengths;
    lengths.reserve(tokens.size());
    
    for (const auto& token : tokens) {
        lengths.push_back(static_cast<double>(token.length()));
    }
    
    double mean = std::accumulate(lengths.begin(), lengths.end(), 0.0) / lengths.size();
    double variance = 0.0;
    
    for (double length : lengths) {
        variance += (length - mean) * (length - mean);
    }
    
    return variance / lengths.size();
}

// 16. Formal charge entropy
SelfiesFormalChargeEntropyDescriptor::SelfiesFormalChargeEntropyDescriptor()
    : SelfiesDescriptor("SelfiesFormalChargeEntropy", "Entropy of formal charge distribution") {}

std::variant<double, int, std::string> SelfiesFormalChargeEntropyDescriptor::calculateFromSelfies(const std::string& selfiesStr) const {
    auto tokens = tokenizeSelfies(selfiesStr);
    
    std::unordered_map<int, int> chargeCounts;
    int totalChargedTokens = 0;
    
    for (const auto& token : tokens) {
        if (isChargeToken(token)) {
            int charge = 0;
            
            if (token.find("+") != std::string::npos) {
                size_t pos = token.find("+");
                size_t count = 1;
                while (pos + count < token.size() && token[pos + count] == '+') count++;
                charge = count;
            } else if (token.find("-") != std::string::npos) {
                size_t pos = token.find("-");
                size_t count = 1;
                while (pos + count < token.size() && token[pos + count] == '-') count++;
                charge = -count;
            }
            
            chargeCounts[charge]++;
            totalChargedTokens++;
        }
    }
    
    if (totalChargedTokens == 0) return 0.0;
    
    double entropy = 0.0;
    for (const auto& [charge, count] : chargeCounts) {
        double p = static_cast<double>(count) / totalChargedTokens;
        entropy -= p * std::log2(p);
    }
    
    return entropy;
}

// 17. Element diversity
SelfiesElementDiversityDescriptor::SelfiesElementDiversityDescriptor()
    : SelfiesDescriptor("SelfiesElementDiversity", "Count of different element types in SELFIES") {}

std::variant<double, int, std::string> SelfiesElementDiversityDescriptor::calculateFromSelfies(const std::string& selfiesStr) const {
    auto tokens = tokenizeSelfies(selfiesStr);
    std::unordered_set<std::string> elements;
    
    for (const auto& token : tokens) {
        if (isAtomToken(token)) {
            // Extract element symbol from token (basic approach)
            size_t start = token.find("[") + 1;
            size_t end = token.find_first_of("@+-=:#.]", start);
            if (end == std::string::npos) end = token.find("]", start);
            
            if (start < token.size() && end != std::string::npos) {
                std::string element = token.substr(start, end - start);
                elements.insert(element);
            }
        }
    }
    
    return static_cast<int>(elements.size());
}

// 18. Levenshtein compression
SelfiesLevenshteinCompressionDescriptor::SelfiesLevenshteinCompressionDescriptor()
    : SelfiesDescriptor("SelfiesLevenshteinCompression", "Measure of redundancy in SELFIES") {}

std::variant<double, int, std::string> SelfiesLevenshteinCompressionDescriptor::calculateFromSelfies(const std::string& selfiesStr) const {
    // Split string into halves and compute Levenshtein distance
    if (selfiesStr.length() < 4) return 0.0;
    
    size_t midpoint = selfiesStr.length() / 2;
    std::string firstHalf = selfiesStr.substr(0, midpoint);
    std::string secondHalf = selfiesStr.substr(midpoint);
    
    int distance = calculateLevenshteinDistance(firstHalf, secondHalf);
    
    // Normalize by maximum possible distance (length of longer half)
    double maxDistance = std::max(firstHalf.length(), secondHalf.length());
    if (maxDistance == 0) return 0.0;
    
    return 1.0 - (distance / maxDistance);
}

// 19. Loop complexity
SelfiesLoopComplexityDescriptor::SelfiesLoopComplexityDescriptor()
    : SelfiesDescriptor("SelfiesLoopComplexity", "Complexity score based on ring tokens") {}

std::variant<double, int, std::string> SelfiesLoopComplexityDescriptor::calculateFromSelfies(const std::string& selfiesStr) const {
    auto tokens = tokenizeSelfies(selfiesStr);
    
    double complexity = 0.0;
    std::unordered_map<int, int> ringSizes;
    
    for (const auto& token : tokens) {
        if (isRingToken(token)) {
            int size = getRingSizeFromToken(token);
            if (size > 0) {
                ringSizes[size]++;
                complexity += std::log(size);
            }
        }
    }
    
    // Add a factor for diversity of ring sizes
    complexity *= (1.0 + 0.1 * ringSizes.size());
    
    return complexity;
}

// 20. Unpaired electron ratio
SelfiesUnpairedElectronRatioDescriptor::SelfiesUnpairedElectronRatioDescriptor()
    : SelfiesDescriptor("SelfiesUnpairedElectronRatio", "Ratio of tokens with unpaired electrons") {}

std::variant<double, int, std::string> SelfiesUnpairedElectronRatioDescriptor::calculateFromSelfies(const std::string& selfiesStr) const {
    auto tokens = tokenizeSelfies(selfiesStr);
    if (tokens.empty()) return 0.0;
    
    int unpaired = 0;
    for (const auto& token : tokens) {
        if (isUnpairedElectronToken(token)) {
            unpaired++;
        }
    }
    
    return static_cast<double>(unpaired) / tokens.size();
}

} // namespace descriptors
} // namespace desfact
