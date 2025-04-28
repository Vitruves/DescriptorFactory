#include "descriptors/strings.hpp"
#include <bitset>
#include <numeric>
#include <unordered_set>
#include <cmath>
#include <sstream>
#include <map> // Needed for ordered map for run entropy

namespace desfact {
namespace descriptors {

// Helper struct to store results from a single pass
struct StringStats {
    double asciiSum = 0.0;
    unsigned char minAscii = 255;
    unsigned char maxAscii = 0;
    int bitCount = 0;
    int oddPositionSetBits = 0;
    int totalOddBitPositions = 0;
    std::unordered_map<char, int> charHistogram;
    std::vector<uint8_t> bitVector;
    int lowerCaseCount = 0;
    int upperCaseCount = 0;
    int digitCount = 0;
    int specialCharCount = 0;
    int caseTransitions = 0;
    int openBrackets = 0;
    int closeBrackets = 0;
    int bondSymbols = 0;
    int ringDigits = 0;
    std::vector<std::string> simpleTokens; // For basic tokenization
};

// Function to compute stats in a single pass (Optimized Bit Handling)
StringStats computeStringStats(const std::string& str) {
    StringStats stats;
    if (str.empty()) return stats;

    stats.bitVector.reserve(str.length() * 8);
    char prevChar = 0;

    for (size_t i = 0; i < str.length(); ++i) {
        char c = str[i];
        unsigned char uc = static_cast<unsigned char>(c);

        // ASCII stats
        stats.asciiSum += uc;
        stats.minAscii = std::min(stats.minAscii, uc);
        stats.maxAscii = std::max(stats.maxAscii, uc);
        stats.charHistogram[c]++;

        // Bit stats (Optimized - using direct bitwise operations)
        for (int j = 0; j < 8; ++j) {
            bool bit_is_set = (uc >> j) & 1;
            stats.bitVector.push_back(bit_is_set ? uint8_t(1) : uint8_t(0));
            if (bit_is_set) {
                stats.bitCount++;
            }
            if (j % 2 != 0) { // Odd position (1, 3, 5, 7)
                stats.totalOddBitPositions++;
                if (bit_is_set) {
                    stats.oddPositionSetBits++;
                }
            }
        }

        // Character type counts
        if (std::islower(c)) stats.lowerCaseCount++;
        else if (std::isupper(c)) stats.upperCaseCount++;
        else if (std::isdigit(c)) {
            stats.digitCount++;
            stats.ringDigits++;
        } else if (!std::isalnum(c)) {
            stats.specialCharCount++;
            if (c == '(') stats.openBrackets++;
            else if (c == ')') stats.closeBrackets++;
            else if (c == '-' || c == '=' || c == '#') stats.bondSymbols++;
        }

        // Case transitions
        if (i > 0) {
             if ((std::islower(prevChar) && std::isupper(c)) ||
                 (std::isupper(prevChar) && std::islower(c))) {
                 stats.caseTransitions++;
             }
        }
        prevChar = c;

        // Basic Tokenization (remains the same for now)
        if (std::isalpha(c) || c == '[' || c == ']' || c == '(' || c == ')' || std::isdigit(c) || c == '-' || c == '=' || c == '#') {
            stats.simpleTokens.push_back(std::string(1, c));
        }
    }
    return stats;
}

// StringDescriptor helper methods implementation
double StringDescriptor::calcAsciiSum(const std::string& str) {
    double sum = 0.0;
    for (char c : str) {
        sum += static_cast<unsigned char>(c);
    }
    return sum;
}

double StringDescriptor::calcAsciiAverage(const std::string& str) {
    if (str.empty()) return 0.0;
    return calcAsciiSum(str) / str.length();
}

double StringDescriptor::calcAsciiVariance(const std::string& str) {
    if (str.empty()) return 0.0;
    
    double avg = calcAsciiAverage(str);
    double variance = 0.0;
    
    for (char c : str) {
        double diff = static_cast<unsigned char>(c) - avg;
        variance += diff * diff;
    }
    
    return variance / str.length();
}

double StringDescriptor::calcBitCount(const std::string& str) {
    int count = 0;
    for (char c : str) {
        std::bitset<8> bits(static_cast<unsigned char>(c));
        count += bits.count();
    }
    return static_cast<double>(count);
}

double StringDescriptor::calcEntropy(const std::string& str) {
    if (str.empty()) return 0.0;
    
    std::unordered_map<char, int> freqs;
    for (char c : str) {
        freqs[c]++;
    }
    
    double entropy = 0.0;
    double len = static_cast<double>(str.length());
    
    for (const auto& pair : freqs) {
        double p = pair.second / len;
        entropy -= p * std::log2(p);
    }
    
    return entropy;
}

double StringDescriptor::calcCharacterFrequency(const std::string& str, char c) {
    if (str.empty()) return 0.0;
    
    int count = 0;
    for (char ch : str) {
        if (ch == c) count++;
    }
    
    return static_cast<double>(count) / str.length();
}

std::unordered_map<char, int> StringDescriptor::getCharacterHistogram(const std::string& str) {
    std::unordered_map<char, int> histogram;
    for (char c : str) {
        histogram[c]++;
    }
    return histogram;
}

int StringDescriptor::getLongestRun(const std::string& str, char c) {
    int longest = 0;
    int current = 0;
    
    if (c == 0) {  // Find longest run of any repeated character
        if (str.empty()) return 0;
        
        char prev = str[0];
        current = 1;
        
        for (size_t i = 1; i < str.length(); ++i) {
            if (str[i] == prev) {
                current++;
            } else {
                longest = std::max(longest, current);
                current = 1;
                prev = str[i];
            }
        }
        
        longest = std::max(longest, current);
    } else {  // Find longest run of specific character
        for (char ch : str) {
            if (ch == c) {
                current++;
            } else {
                longest = std::max(longest, current);
                current = 0;
            }
        }
        
        longest = std::max(longest, current);
    }
    
    return longest;
}

double StringDescriptor::getSymmetryScore(const std::string& str) {
    if (str.length() <= 1) return 1.0;
    
    int matches = 0;
    int total = str.length() / 2;
    
    for (size_t i = 0; i < static_cast<size_t>(total); ++i) {
        if (str[i] == str[str.length() - 1 - i]) {
            matches++;
        }
    }
    
    return static_cast<double>(matches) / total;
}

// ASCII and Bit-Level Descriptors implementation (Refactored)
std::variant<double, int, std::string> AsciiSumDescriptor::calculateFromSmiles(const std::string& smiles) const {
    StringStats stats = computeStringStats(smiles);
    return stats.asciiSum;
}

std::variant<double, int, std::string> AsciiAverageDescriptor::calculateFromSmiles(const std::string& smiles) const {
    if (smiles.empty()) return 0.0;
    StringStats stats = computeStringStats(smiles);
    return smiles.length() > 0 ? stats.asciiSum / smiles.length() : 0.0;
}

std::variant<double, int, std::string> AsciiVarianceDescriptor::calculateFromSmiles(const std::string& smiles) const {
    if (smiles.empty()) return 0.0;
    StringStats stats = computeStringStats(smiles);
    double avg = stats.asciiSum / smiles.length();
    double variance = 0.0;
    for (char c : smiles) {
        double diff = static_cast<unsigned char>(c) - avg;
        variance += diff * diff;
    }
    return variance / smiles.length();
}

std::variant<double, int, std::string> AsciiRangeDescriptor::calculateFromSmiles(const std::string& smiles) const {
    if (smiles.empty()) return 0.0;
    StringStats stats = computeStringStats(smiles);
    return static_cast<double>(stats.maxAscii - stats.minAscii);
}

std::variant<double, int, std::string> BitCountDescriptor::calculateFromSmiles(const std::string& smiles) const {
    StringStats stats = computeStringStats(smiles);
    return static_cast<double>(stats.bitCount);
}

std::variant<double, int, std::string> BitDensityDescriptor::calculateFromSmiles(const std::string& smiles) const {
    if (smiles.empty()) return 0.0;
    StringStats stats = computeStringStats(smiles);
    return static_cast<double>(stats.bitCount) / (smiles.length() * 8.0);
}


std::variant<double, int, std::string> BitTransitionRateDescriptor::calculateFromSmiles(const std::string& smiles) const {
     if (smiles.length() * 8 <= 1) return 0.0; 
     StringStats stats = computeStringStats(smiles);
     int transitions = 0;
     for (size_t i = 0; i < stats.bitVector.size() - 1; ++i) {
         // Compare uint8_t values
         if (stats.bitVector[i] != stats.bitVector[i + 1]) { 
             transitions++;
         }
     }
     return static_cast<double>(transitions) / (stats.bitVector.size() - 1);
}

std::variant<double, int, std::string> OddBitRatioDescriptor::calculateFromSmiles(const std::string& smiles) const {
    if (smiles.empty()) return 0.0;
    StringStats stats = computeStringStats(smiles);
    if (stats.totalOddBitPositions == 0) return 0.0;
    return static_cast<double>(stats.oddPositionSetBits) / stats.totalOddBitPositions;
}

// Compression & Entropy Descriptors implementation
std::variant<double, int, std::string> RunLengthEncodingSizeDescriptor::calculateFromSmiles(const std::string& smiles) const {
    if (smiles.empty()) return 0.0;
    
    std::stringstream rle;
    char current = smiles[0];
    int count = 1;
    
    for (size_t i = 1; i < smiles.length(); ++i) {
        if (smiles[i] == current) {
            count++;
        } else {
            rle << current;
            if (count > 1) {
                rle << count;
            }
            current = smiles[i];
            count = 1;
        }
    }
    
    rle << current;
    if (count > 1) {
        rle << count;
    }
    
    std::string encoded = rle.str();
    return static_cast<double>(encoded.length());
}

std::variant<double, int, std::string> CharacterRepetitionScoreDescriptor::calculateFromSmiles(const std::string& smiles) const {
    if (smiles.empty()) return 0.0;
    StringStats stats = computeStringStats(smiles);
    
    double score = 0.0;
    for (const auto& pair : stats.charHistogram) {
        score += pair.second * pair.second; // Square the counts
    }
    
    return score / (smiles.length() * smiles.length());
}

std::variant<double, int, std::string> NibblePatternCountDescriptor::calculateFromSmiles(const std::string& smiles) const {
    if (smiles.empty()) return 0.0;
    
    std::unordered_set<int> patterns;
    
    for (char c : smiles) {
        unsigned char byte = static_cast<unsigned char>(c);
        // Extract high and low nibbles (4 bits each)
        int highNibble = (byte >> 4) & 0xF;
        int lowNibble = byte & 0xF;
        
        patterns.insert(highNibble);
        patterns.insert(lowNibble);
    }
    
    return static_cast<double>(patterns.size());
}

std::variant<double, int, std::string> BytePatternCountDescriptor::calculateFromSmiles(const std::string& smiles) const {
    if (smiles.empty()) return 0.0;
    
    std::unordered_set<unsigned char> patterns;
    
    for (char c : smiles) {
        patterns.insert(static_cast<unsigned char>(c));
    }
    
    return static_cast<double>(patterns.size());
}

std::variant<double, int, std::string> EntropyPerByteDescriptor::calculateFromSmiles(const std::string& smiles) const {
    if (smiles.empty()) return 0.0;
    
    return calcEntropy(smiles) / smiles.length();
}

std::variant<double, int, std::string> AsciiHistogramUniformityDescriptor::calculateFromSmiles(const std::string& smiles) const {
    if (smiles.empty()) return 0.0;
    StringStats stats = computeStringStats(smiles);
    const auto& histogram = stats.charHistogram;
    if (histogram.empty()) return 1.0; // Perfectly uniform if empty or single character type

    double mean = static_cast<double>(smiles.length()) / histogram.size();
    if (mean <= 1e-9) return 1.0; // Avoid division by zero if mean is effectively zero

    double variance = 0.0;
    for (const auto& pair : histogram) {
        double diff = pair.second - mean;
        variance += diff * diff;
    }
    variance /= histogram.size(); // Variance of counts

    // Coefficient of variation
    double cv = std::sqrt(variance) / mean;
    // Transform to [0,1] range where 1 is perfectly uniform
    return 1.0 / (1.0 + cv);
}

// Character Sequence Descriptors implementation
std::variant<double, int, std::string> LongestCharacterRunDescriptor::calculateFromSmiles(const std::string& smiles) const {
    return static_cast<double>(getLongestRun(smiles));
}

std::variant<double, int, std::string> SpecialCharacterDensityDescriptor::calculateFromSmiles(const std::string& smiles) const {
    if (smiles.empty()) return 0.0;
    StringStats stats = computeStringStats(smiles);
    return static_cast<double>(stats.specialCharCount) / smiles.length();
}

std::variant<double, int, std::string> CaseTransitionRateDescriptor::calculateFromSmiles(const std::string& smiles) const {
    if (smiles.length() <= 1) return 0.0;
    StringStats stats = computeStringStats(smiles);
    return static_cast<double>(stats.caseTransitions) / (smiles.length() - 1);
}

std::variant<double, int, std::string> AlternatingCharacterScoreDescriptor::calculateFromSmiles(const std::string& smiles) const {
    if (smiles.length() <= 2) return 0.0;
    
    int alternations = 0;
    size_t comparisons = smiles.length() - 2;
    
    for (size_t i = 0; i < comparisons; ++i) {
        if (smiles[i] == smiles[i+2] && smiles[i] != smiles[i+1]) {
            alternations++;
        }
    }
    
    return static_cast<double>(alternations) / comparisons;
}

std::variant<double, int, std::string> CharacterTrigramCountDescriptor::calculateFromSmiles(const std::string& smiles) const {
    if (smiles.length() < 3) return 0.0;
    
    std::unordered_set<std::string> trigrams;
    
    for (size_t i = 0; i <= smiles.length() - 3; ++i) {
        trigrams.insert(smiles.substr(i, 3));
    }
    
    // Normalize by possible number of trigrams? Or just count? Let's just count for now.
    return static_cast<double>(trigrams.size()); 
}

std::variant<double, int, std::string> LetterDigitRatioDescriptor::calculateFromSmiles(const std::string& smiles) const {
    if (smiles.empty()) return 0.0;
    StringStats stats = computeStringStats(smiles);
    int letters = stats.lowerCaseCount + stats.upperCaseCount;
    int digits = stats.digitCount;
    if (digits == 0) return static_cast<double>(letters); // Avoid division by zero
    return static_cast<double>(letters) / digits;
}

std::variant<double, int, std::string> CharacterBigramEntropyDescriptor::calculateFromSmiles(const std::string& smiles) const {
    if (smiles.length() < 2) return 0.0;
    
    std::unordered_map<std::string, int> bigramCounts;
    int totalBigrams = smiles.length() - 1;

    for (size_t i = 0; i < static_cast<size_t>(totalBigrams); ++i) {
        std::string bigram = smiles.substr(i, 2);
        bigramCounts[bigram]++;
    }
    
    double entropy = 0.0;    
    for (const auto& pair : bigramCounts) {
        double p = static_cast<double>(pair.second) / totalBigrams;
         if (p > 0) { // Avoid log2(0)
            entropy -= p * std::log2(p);
        }
    }
    
    return entropy;
}

// Bit Patterns Descriptors implementation
std::variant<double, int, std::string> BytePalindromeScoreDescriptor::calculateFromSmiles(const std::string& smiles) const {
    return getSymmetryScore(smiles);
}

std::variant<double, int, std::string> BitPalindromeScoreDescriptor::calculateFromSmiles(const std::string& smiles) const {
    if (smiles.empty()) return 0.0;
    StringStats stats = computeStringStats(smiles);
    const auto& bits = stats.bitVector; // Use alias for clarity
    if (bits.empty()) return 0.0; // Add check for empty vector

    int matches = 0;
    int total = bits.size() / 2;
    if (total == 0) return 1.0; // Single bit is a palindrome

    for (int i = 0; i < total; ++i) {
        // Compare uint8_t values
        if (bits[i] == bits[bits.size() - 1 - i]) { 
            matches++;
        }
    }
    
    return static_cast<double>(matches) / total;
}

std::variant<double, int, std::string> BitPeriodicityDescriptor::calculateFromSmiles(const std::string& smiles) const {
    if (smiles.length() < 2) return 0.0;
    StringStats stats = computeStringStats(smiles);
    const auto& bits = stats.bitVector;
    size_t n = bits.size(); // Use size_t for size
    if (n < 4) return 0.0;

    // Optimization: Limit maximum period check
    size_t maxPeriodToCheck = std::min(n / 2, static_cast<size_t>(128)); // Check periods up to 128 bits or n/2, whichever is smaller

    double bestMatch = 0.0;

    for (size_t period = 1; period <= maxPeriodToCheck; ++period) { // Use size_t for period
        int matches = 0;
        size_t comparisons = n - period; // Correct number of comparisons

        if (comparisons == 0) continue; // Avoid division by zero if period is too large relative to n

        for (size_t i = 0; i < comparisons; ++i) {
             if (bits[i] == bits[i + period]) {
                matches++;
            }
        }

        double match = static_cast<double>(matches) / comparisons;
        if (match > bestMatch) {
            bestMatch = match;
            // We're only returning the match quality, not the period itself
        }
        // Optimization: Early exit if match rate is extremely high (optional)
        // if (bestMatch > 0.95) { 
        //     break;
        // }
    }

    return bestMatch;
}

std::variant<double, int, std::string> BitAutocorrelationDescriptor::calculateFromSmiles(const std::string& smiles) const {
    if (smiles.empty()) return 0.0;
    StringStats stats = computeStringStats(smiles);
    const auto& bits = stats.bitVector;
    if (bits.size() < 2) return 0.0;

    int matches = 0;
    for (size_t i = 0; i < bits.size() - 1; ++i) {
        // Compare uint8_t values
        if (bits[i] == bits[i + 1]) { 
            matches++;
        }
    }
    
    return static_cast<double>(matches) / (bits.size() - 1);
}

std::variant<double, int, std::string> HammingWeightDistributionDescriptor::calculateFromSmiles(const std::string& smiles) const {
    if (smiles.empty()) return 0.0;

    std::vector<int> weights;
    weights.reserve(smiles.length()); 
    for (char c : smiles) {
        unsigned char uc = static_cast<unsigned char>(c);
        int weight = 0;
        // Count bits directly
        for(int k=0; k<8; ++k) {
            if((uc >> k) & 1) {
                weight++;
            }
        }
        weights.push_back(weight);
    }

    if (weights.empty()) return 0.0; 

    // Calculate variance of Hamming weights
    double sum = std::accumulate(weights.begin(), weights.end(), 0.0);
    double mean = sum / weights.size();
    double variance = 0.0;

    for (int w : weights) {
        double diff = w - mean;
        variance += diff * diff;
    }
    variance /= weights.size();

    // Max possible variance is 16 (half 0s, half 8s -> mean 4)
    return std::min(1.0, variance / 16.0); // Normalize to [0,1]
}

std::variant<double, int, std::string> ByteParityDistributionDescriptor::calculateFromSmiles(const std::string& smiles) const {
    if (smiles.empty()) return 0.0;

    int evenParityCount = 0;

    for (char c : smiles) {
        unsigned char uc = static_cast<unsigned char>(c);
        int weight = 0;
        // Count bits directly
        for(int k=0; k<8; ++k) {
            if((uc >> k) & 1) {
                weight++;
            }
        }
        if (weight % 2 == 0) {
            evenParityCount++;
        }
    }

    return static_cast<double>(evenParityCount) / smiles.length();
}

std::variant<double, int, std::string> BitRunEntropyDescriptor::calculateFromSmiles(const std::string& smiles) const {
    if (smiles.empty()) return 0.0;
    StringStats stats = computeStringStats(smiles);
    const auto& bits = stats.bitVector;
    if (bits.empty()) return 0.0;

    std::map<int, int> runLengths; // Use map for potentially ordered runs if needed
    int currentRun = 1;
    
    for (size_t i = 1; i < bits.size(); ++i) {
        // Compare uint8_t values
        if (bits[i] == bits[i-1]) { 
            currentRun++;
        } else {
            runLengths[currentRun]++;
            currentRun = 1;
        }
    }
    runLengths[currentRun]++;
    
    // Calculate entropy of run lengths
    double entropy = 0.0;
    int totalRuns = 0;
    
    for (const auto& pair : runLengths) {
        totalRuns += pair.second;
    }
    
    if (totalRuns == 0) return 0.0; // Avoid division by zero

    for (const auto& pair : runLengths) {
        double p = static_cast<double>(pair.second) / totalRuns;
        if (p > 0) { // Avoid log2(0)
            entropy -= p * std::log2(p);
        }
    }
    
    return entropy;
}

// Positional Heuristics Descriptors
std::variant<double, int, std::string> PositionWeightDescriptor::calculateFromSmiles(const std::string& smiles) const {
    if (smiles.empty()) return 0.0;
    StringStats stats = computeStringStats(smiles); // Get asciiSum from stats
    if (stats.asciiSum <= 1e-9) return 0.0; // Avoid division by zero

    double weightedSum = 0.0;
    double len_minus_1 = smiles.length() > 1 ? smiles.length() - 1 : 1.0; // Avoid division by zero for length 1

    for (size_t i = 0; i < smiles.length(); ++i) {
        double weight = static_cast<double>(i) / len_minus_1; // Normalize weight to [0, 1]
        weightedSum += static_cast<unsigned char>(smiles[i]) * weight;
    }
    
    return weightedSum / stats.asciiSum; // Normalize by total ASCII sum
}

std::variant<double, int, std::string> FrontBackRatioDescriptor::calculateFromSmiles(const std::string& smiles) const {
    if (smiles.length() < 2) return 1.0;
    
    size_t midpoint = smiles.length() / 2;
    double frontSum = 0.0;
    double backSum = 0.0;

    for(size_t i=0; i<midpoint; ++i) {
        frontSum += static_cast<unsigned char>(smiles[i]);
    }
     for(size_t i=midpoint; i<smiles.length(); ++i) {
        backSum += static_cast<unsigned char>(smiles[i]);
    }
    
    return (backSum > 1e-9) ? frontSum / backSum : frontSum; // Avoid division by zero
}

std::variant<double, int, std::string> BracketAsymmetryDescriptor::calculateFromSmiles(const std::string& smiles) const {
    if (smiles.empty()) return 0.0;
    
    double openSumPos = 0.0;
    double closeSumPos = 0.0;
    int openCount = 0;
    int closeCount = 0;
    
    for (size_t i = 0; i < smiles.length(); ++i) {
        if (smiles[i] == '(') {
            openSumPos += i;
            openCount++;
        } else if (smiles[i] == ')') {
            closeSumPos += i;
            closeCount++;
        }
    }
    
    if (openCount == 0 || closeCount == 0) return 0.0;
    
    double openAvg = openSumPos / openCount;
    double closeAvg = closeSumPos / closeCount;
    
    // Normalize by string length
    return std::abs(openAvg - closeAvg) / smiles.length(); 
}

std::variant<double, int, std::string> CharDispersionDescriptor::calculateFromSmiles(const std::string& smiles) const {
    if (smiles.length() <= 1) return 0.0;
    
    std::unordered_map<char, std::vector<size_t>> charPositions;
    for (size_t i = 0; i < smiles.length(); ++i) {
        charPositions[smiles[i]].push_back(i);
    }
    
    double totalAvgDistance = 0.0;
    int charTypesWithMultipleOccurrences = 0;
    
    for (const auto& pair : charPositions) {
        if (pair.second.size() > 1) {
            double sumDistances = 0.0;
            for (size_t i = 1; i < pair.second.size(); ++i) {
                sumDistances += pair.second[i] - pair.second[i-1];
            }
            totalAvgDistance += sumDistances / (pair.second.size() - 1);
            charTypesWithMultipleOccurrences++;
        }
    }
    
    if (charTypesWithMultipleOccurrences == 0) return 0.0; 
    // Normalize by average distance across types and string length
    return (totalAvgDistance / charTypesWithMultipleOccurrences) / smiles.length(); 
}

std::variant<double, int, std::string> DigitPositionMeanDescriptor::calculateFromSmiles(const std::string& smiles) const {
    if (smiles.empty()) return 0.0;
    
    double digitSumPos = 0.0;
    int digitCount = 0;
    
    for (size_t i = 0; i < smiles.length(); ++i) {
        if (std::isdigit(smiles[i])) {
            digitSumPos += i;
            digitCount++;
        }
    }
    
    if (digitCount == 0) return 0.0;
    
    double mean = digitSumPos / digitCount;
    return mean / smiles.length(); // Normalize by string length
}

std::variant<double, int, std::string> SymmetryScoreDescriptor::calculateFromSmiles(const std::string& smiles) const {
    return getSymmetryScore(smiles);
}

// Pattern / Fragment Pseudodescriptors
std::variant<double, int, std::string> AromaticRatioDescriptor::calculateFromSmiles(const std::string& smiles) const {
    if (smiles.empty()) return 0.0;
    StringStats stats = computeStringStats(smiles);
    int alphaCount = stats.lowerCaseCount + stats.upperCaseCount;
    if (alphaCount == 0) return 0.0;
    return static_cast<double>(stats.lowerCaseCount) / alphaCount;
}

std::variant<double, int, std::string> HeteroatomRatioDescriptor::calculateFromSmiles(const std::string& smiles) const {
    if (smiles.empty()) return 0.0;

    int heteroCount = 0;
    int elementCount = 0;
    
    for (size_t i = 0; i < smiles.length(); ++i) {
        if (std::isalpha(smiles[i])) {
            elementCount++;
            // Basic check - ignores bracket atoms for simplicity
            if (smiles[i] != 'C' && smiles[i] != 'c' && smiles[i] != 'H') { 
                 heteroCount++;
            }
            // Skip potential second lowercase letter of element symbol
            if (i + 1 < smiles.length() && std::islower(smiles[i+1])) {
                i++; 
            }
        } else if (smiles[i] == '[') { // Handle bracket atoms
             size_t endBracket = smiles.find(']', i);
             if (endBracket != std::string::npos) {
                 elementCount++;
                 // Extract element symbol within brackets
                 std::string bracketContent = smiles.substr(i + 1, endBracket - i - 1);
                 // Very basic check, assumes first letter(s) is element
                 if (bracketContent.length() > 0 && std::isalpha(bracketContent[0])) {
                     std::string elementSymbol;
                     elementSymbol += bracketContent[0];
                      if (bracketContent.length() > 1 && std::islower(bracketContent[1])) {
                          elementSymbol += bracketContent[1];
                      }
                      if (elementSymbol != "C" && elementSymbol != "H") {
                          heteroCount++;
                      }
                 }
                 i = endBracket; // Move past the bracket atom
             }
        }
    }

    return (elementCount > 0) ? static_cast<double>(heteroCount) / elementCount : 0.0;
}

std::variant<double, int, std::string> ElementDiversityDescriptor::calculateFromSmiles(const std::string& smiles) const {
    if (smiles.empty()) return 0.0;
    
    std::unordered_set<std::string> elements;
    for (size_t i = 0; i < smiles.length(); ++i) {
        if (smiles[i] == '[') { // Bracket atom
            size_t endBracket = smiles.find(']', i);
            if (endBracket != std::string::npos) {
                // Extract potential element symbol (first 1 or 2 letters)
                if (endBracket > i + 1 && std::isalpha(smiles[i+1])) {
                    std::string elementSymbol;
                    elementSymbol += smiles[i+1];
                    if (endBracket > i + 2 && std::islower(smiles[i+2])) {
                         elementSymbol += smiles[i+2];
                    }
                     elements.insert(elementSymbol);
                }
                 i = endBracket; 
            }
        } else if (std::isupper(smiles[i])) { // Standard element
            std::string elementSymbol;
            elementSymbol += smiles[i];
            if (i + 1 < smiles.length() && std::islower(smiles[i+1])) {
                 elementSymbol += smiles[i+1];
                 i++; // Skip the lowercase letter
            }
             elements.insert(elementSymbol);
        } else if (std::islower(smiles[i])) { // Aromatic element
             elements.insert(std::string(1, smiles[i]));
        }
        // Ignore other characters like bonds, digits, etc. for diversity count
    }
    
    return static_cast<double>(elements.size());
}

std::variant<double, int, std::string> SaturationIndexDescriptor::calculateFromSmiles(const std::string& smiles) const {
    if (smiles.empty()) return 0.0;
    
    int singleBonds = 0;
    int doubleBonds = 0;
    int tripleBonds = 0;
    int aromaticBonds = 0; // Approx based on lowercase letters between elements
    
    for (size_t i = 0; i < smiles.length(); ++i) {
        if (smiles[i] == '-') singleBonds++;
        else if (smiles[i] == '=') doubleBonds++;
        else if (smiles[i] == '#') tripleBonds++;
         // Rudimentary check for aromatic bonds (adjacent lowercase letters)
        else if (i > 0 && std::islower(smiles[i]) && std::isalpha(smiles[i-1])) {
             aromaticBonds++;
         }
    }
    
    double totalExplicitBonds = singleBonds + doubleBonds + tripleBonds + aromaticBonds;
    if (totalExplicitBonds == 0) return 0.5; // Default for molecules with no explicit bonds?
    
    // Saturation index: higher value means more saturated (fewer double/triple/aromatic)
    // Normalize: 1 = fully single, 0 = fully triple/aromatic (approx)
    return (singleBonds + 0.5*aromaticBonds) / totalExplicitBonds; 
}

std::variant<double, int, std::string> PatternEntropyDescriptor::calculateFromSmiles(const std::string& smiles) const {
    if (smiles.empty()) return 0.0;

    std::vector<std::string> tokens;
    std::string currentToken;

    for (size_t i = 0; i < smiles.length(); ++i) {
        char c = smiles[i];

        if (c == '[') { // Start bracket atom
            if (!currentToken.empty()) tokens.push_back(currentToken);
            currentToken = "[";
        } else if (c == ']') { // End bracket atom
            currentToken += "]";
            tokens.push_back(currentToken);
            currentToken = "";
        } else if (!currentToken.empty() && currentToken[0] == '[') { // Inside bracket atom
            currentToken += c;
        } else if (std::isupper(c)) { // Start uppercase element
             if (!currentToken.empty()) tokens.push_back(currentToken);
             currentToken = c;
             if (i + 1 < smiles.length() && std::islower(smiles[i+1])) {
                 currentToken += smiles[i+1];
                 i++; // Consume lowercase char
             }
             tokens.push_back(currentToken);
             currentToken = "";
        } else if (std::islower(c)) { // Aromatic element
             if (!currentToken.empty()) tokens.push_back(currentToken);
             tokens.push_back(std::string(1, c));
             currentToken = "";
        } else if (std::isdigit(c)) { // Ring number or bond configuration
            if (!currentToken.empty()) tokens.push_back(currentToken);
            currentToken = c;
             // Look ahead for multi-digit ring numbers (less common)
             // while (i + 1 < smiles.length() && std::isdigit(smiles[i+1])) {
             //     currentToken += smiles[i+1];
             //     i++;
             // }
             tokens.push_back(currentToken);
             currentToken = "";
        } else if (c == '%' && i + 2 < smiles.length()) { // Multi-digit ring number
             if (!currentToken.empty()) tokens.push_back(currentToken);
             tokens.push_back(smiles.substr(i, 3));
             currentToken = "";
             i += 2;
        } else if (c == '-' || c == '=' || c == '#' || c == ':' || c == '(' || c == ')' || c == '.' || c == '+' || c == '*' || c == '$') { // Single char tokens
             if (!currentToken.empty()) tokens.push_back(currentToken);
             tokens.push_back(std::string(1, c));
             currentToken = "";
        } else {
             // Append to current if it's part of a weird bracket atom or ignore?
             // For now, just ignore unexpected characters outside brackets
             if (!currentToken.empty() && currentToken[0] != '[') {
                  tokens.push_back(currentToken);
                  currentToken = "";
             }
        }
    }
     if (!currentToken.empty()) tokens.push_back(currentToken); // Add last token


    if (tokens.empty()) return 0.0;
    
    std::unordered_map<std::string, int> tokenCounts;
    for (const std::string& token : tokens) {
        tokenCounts[token]++;
    }
    
    double entropy = 0.0;
    double totalTokens = tokens.size();
    for (const auto& pair : tokenCounts) {
        double p = static_cast<double>(pair.second) / totalTokens;
         if (p > 0) {
            entropy -= p * std::log2(p);
        }
    }
    
    return entropy;
}

std::variant<double, int, std::string> BondComplexityDescriptor::calculateFromSmiles(const std::string& smiles) const {
    if (smiles.empty()) return 0.0;
    
    int singleBonds = 0;
    int doubleBonds = 0;
    int tripleBonds = 0;
    int aromaticBonds = 0; // Use colon for aromatic
    
    for (char c : smiles) {
        if (c == '-') singleBonds++;
        else if (c == '=') doubleBonds++;
        else if (c == '#') tripleBonds++;
        else if (c == ':') aromaticBonds++;
    }
    
    int totalExplicitBonds = singleBonds + doubleBonds + tripleBonds + aromaticBonds;
    if (totalExplicitBonds == 0) return 0.0; // No explicit bonds found
    
    // Complexity score: higher orders and diversity increase complexity
    double complexityScore = (1.0 * singleBonds + 1.5 * aromaticBonds + 2.0 * doubleBonds + 3.0 * tripleBonds) / totalExplicitBonds;
    
    int bondTypes = (singleBonds > 0) + (doubleBonds > 0) + (tripleBonds > 0) + (aromaticBonds > 0);
    
    // Normalize roughly to [0, 1] - max score around 3-4? Divide by 4.
    return std::min(1.0, (complexityScore * bondTypes) / 4.0); 
}

std::variant<double, int, std::string> RingComplexityDescriptor::calculateFromSmiles(const std::string& smiles) const {
    if (smiles.empty()) return 0.0;
    
    std::map<int, int> ringCounts; // Map digit value to count
    int totalRingMarkers = 0;
    
    for (size_t i = 0; i < smiles.length(); ++i) {
        if (std::isdigit(smiles[i])) {
             int ringNum = smiles[i] - '0';
             ringCounts[ringNum]++;
             totalRingMarkers++;
        } else if (smiles[i] == '%' && i + 2 < smiles.length() && std::isdigit(smiles[i+1]) && std::isdigit(smiles[i+2])) {
             // Handle multi-digit rings like %10, %11...
             int ringNum = (smiles[i+1] - '0') * 10 + (smiles[i+2] - '0');
             ringCounts[ringNum]++;
             totalRingMarkers += 2; // Count both digits after %
             i += 2; // Skip the digits
        }
    }
    
    if (totalRingMarkers == 0) return 0.0;
    
    int uniqueRingIndices = ringCounts.size();
    
    // Ensure pairs: divide total markers by 2 for number of closures
    int ringClosures = totalRingMarkers / 2; 
    
    // Complexity increases with more unique rings and more total closures
    // Normalize by string length (heuristic)
    return std::min(1.0, (static_cast<double>(uniqueRingIndices * ringClosures)) / (smiles.length() * 5.0)); 
}

// Basic SMILES Structure-Based Descriptors

std::variant<double, int, std::string> BondPositionEntropyDescriptor::calculateFromSmiles(const std::string& smiles) const {
    if (smiles.empty()) return 0.0;
    
    std::vector<size_t> bondPositions;
    
    for (size_t i = 0; i < smiles.length(); ++i) {
        if (smiles[i] == '-' || smiles[i] == '=' || smiles[i] == '#' || smiles[i] == ':') { // Include aromatic bond
            bondPositions.push_back(i);
        }
    }
    
    if (bondPositions.size() <= 1) return 0.0; // Need at least two bonds for distances
    
    // Calculate entropy of distances between bond symbols
    std::vector<size_t> distances;
    for (size_t i = 1; i < bondPositions.size(); ++i) {
        distances.push_back(bondPositions[i] - bondPositions[i-1]);
    }
    
    std::unordered_map<size_t, int> distanceCounts;
    for (size_t dist : distances) {
        distanceCounts[dist]++;
    }
    
    double entropy = 0.0;
    double totalDistances = distances.size();
    if (totalDistances == 0) return 0.0;

    for (const auto& pair : distanceCounts) {
        double p = static_cast<double>(pair.second) / totalDistances;
        if (p > 0) { // Avoid log2(0)
             entropy -= p * std::log2(p);
        }
    }
    
    return entropy;
}

std::variant<double, int, std::string> ElementPositionBiasDescriptor::calculateFromSmiles(const std::string& smiles) const {
    if (smiles.empty()) return 0.0;
    
    std::unordered_map<std::string, std::vector<size_t>> elementPositions;
    
    // Manual parsing for elements
     for (size_t i = 0; i < smiles.length(); ++i) {
        if (smiles[i] == '[') { 
            size_t endBracket = smiles.find(']', i);
            if (endBracket != std::string::npos) {
                 if (endBracket > i + 1 && std::isalpha(smiles[i+1])) {
                    std::string elementSymbol;
                    elementSymbol += smiles[i+1];
                    if (endBracket > i + 2 && std::islower(smiles[i+2])) {
                         elementSymbol += smiles[i+2];
                    }
                    elementPositions[elementSymbol].push_back(i); // Position of '[' ? or center? Use start for now.
                }
                 i = endBracket; 
            }
        } else if (std::isupper(smiles[i])) { 
            std::string elementSymbol;
            elementSymbol += smiles[i];
            if (i + 1 < smiles.length() && std::islower(smiles[i+1])) {
                 elementSymbol += smiles[i+1];
                 elementPositions[elementSymbol].push_back(i);
                 i++; 
            } else {
                 elementPositions[elementSymbol].push_back(i);
            }
        } else if (std::islower(smiles[i])) { // Aromatic element
             elementPositions[std::string(1, smiles[i])].push_back(i);
        }
    }

    if (elementPositions.empty()) return 0.0;

    // Calculate positional bias as average standard deviation of normalized positions
    double totalStdDev = 0.0;
    int elementTypesCounted = 0;
    
    for (const auto& pair : elementPositions) {
        if (pair.second.size() > 1) {
            std::vector<double> normPositions;
            for (size_t pos : pair.second) {
                normPositions.push_back(static_cast<double>(pos) / smiles.length());
            }
            
            double sum = std::accumulate(normPositions.begin(), normPositions.end(), 0.0);
            double mean = sum / normPositions.size();
            
            double variance = 0.0;
            for (double pos : normPositions) {
                double diff = pos - mean;
                variance += diff * diff;
            }
            variance /= normPositions.size();
            
            totalStdDev += std::sqrt(variance);
            elementTypesCounted++;
        }
    }
    
    return (elementTypesCounted > 0) ? totalStdDev / elementTypesCounted : 0.0;
}


std::variant<double, int, std::string> BranchComplexityDescriptor::calculateFromSmiles(const std::string& smiles) const {
    if (smiles.empty()) return 0.0;
    
    int maxNesting = 0;
    int currentNesting = 0;
    int totalBranches = 0; // Count number of '('
    
    for (char c : smiles) {
        if (c == '(') {
            currentNesting++;
            maxNesting = std::max(maxNesting, currentNesting);
            totalBranches++;
        } else if (c == ')') {
            currentNesting--;
             if (currentNesting < 0) { /* Handle potentially malformed SMILES */ currentNesting = 0;}
        }
    }
    
    // Complexity based on max nesting depth and number of branches per length
    double complexity = static_cast<double>(maxNesting) + static_cast<double>(totalBranches);
    return complexity / smiles.length();
}

std::variant<double, int, std::string> ChainLengthDescriptor::calculateFromSmiles(const std::string& smiles) const {
    if (smiles.empty()) return 0.0;
    
    int maxChainLength = 0;
    int currentChainLength = 0;
    int nestedLevel = 0;
    
    for (char c : smiles) {
        if (c == '(') {
            nestedLevel++;
        } else if (c == ')') {
            nestedLevel--;
        } else if (c == '.') {
            maxChainLength = std::max(maxChainLength, currentChainLength);
            currentChainLength = 0;
        } else if (nestedLevel == 0 && (std::isalpha(c) || c == '[' || c == ']')) {
            currentChainLength++;
        }
    }
    
    maxChainLength = std::max(maxChainLength, currentChainLength);
    return static_cast<double>(maxChainLength);
}

std::variant<double, int, std::string> DeviationSymmetryDescriptor::calculateFromSmiles(const std::string& smiles) const {
    if (smiles.length() <= 1) return 1.0;
    
    int len = smiles.length();
    int halfLen = len / 2;
    double avgDiff = 0.0;
    
    for (int i = 0; i < halfLen; ++i) {
        unsigned char leftChar = static_cast<unsigned char>(smiles[i]);
        unsigned char rightChar = static_cast<unsigned char>(smiles[len - 1 - i]);
        double normalizedDiff = std::abs(static_cast<double>(leftChar - rightChar)) / 255.0;
        avgDiff += normalizedDiff;
    }
    
    return 1.0 - (avgDiff / halfLen);
}

std::variant<double, int, std::string> LocalComplexityDescriptor::calculateFromSmiles(const std::string& smiles) const {
    if (smiles.length() < 5) return 0.0;
    
    int windowSize = std::min(5, static_cast<int>(smiles.length() / 2));
    double totalComplexity = 0.0;
    int windows = 0;
    
    for (size_t i = 0; i <= smiles.length() - windowSize; ++i) {
        std::string window = smiles.substr(i, windowSize);
        std::unordered_set<char> uniqueChars(window.begin(), window.end());
        double localEntropy = static_cast<double>(uniqueChars.size()) / windowSize;
        totalComplexity += localEntropy;
        windows++;
    }
    
    return totalComplexity / windows;
}

std::variant<double, int, std::string> CharacterKurtosisDescriptor::calculateFromSmiles(const std::string& smiles) const {
    if (smiles.length() < 4) return 0.0;
    
    StringStats stats = computeStringStats(smiles);
    double mean = stats.asciiSum / smiles.length();
    double variance = 0.0;
    double fourthMoment = 0.0;
    
    for (char c : smiles) {
        double diff = static_cast<unsigned char>(c) - mean;
        double squaredDiff = diff * diff;
        variance += squaredDiff;
        fourthMoment += squaredDiff * squaredDiff;
    }
    
    variance /= smiles.length();
    fourthMoment /= smiles.length();
    
    if (variance < 1e-10) return 0.0;
    double kurtosis = (fourthMoment / (variance * variance)) - 3.0;
    
    return kurtosis;
}

std::variant<double, int, std::string> CharacterSkewnessDescriptor::calculateFromSmiles(const std::string& smiles) const {
    if (smiles.length() < 3) return 0.0;
    
    StringStats stats = computeStringStats(smiles);
    double mean = stats.asciiSum / smiles.length();
    double variance = 0.0;
    double thirdMoment = 0.0;
    
    for (char c : smiles) {
        double diff = static_cast<unsigned char>(c) - mean;
        variance += diff * diff;
        thirdMoment += diff * diff * diff;
    }
    
    variance /= smiles.length();
    thirdMoment /= smiles.length();
    
    if (variance < 1e-10) return 0.0;
    double standardDeviation = std::sqrt(variance);
    double skewness = thirdMoment / std::pow(standardDeviation, 3);
    
    return skewness;
}

// Sequence Analysis Features implementations
std::variant<double, int, std::string> SequenceNoisinessDescriptor::calculateFromSmiles(const std::string& smiles) const {
    if (smiles.length() < 3) return 0.0;
    
    int transitions = 0;
    for (size_t i = 1; i < smiles.length(); ++i) {
        if (smiles[i] != smiles[i-1]) {
            transitions++;
        }
    }
    
    return static_cast<double>(transitions) / (smiles.length() - 1);
}

std::variant<double, int, std::string> CharacterTransitionMatrixEntropyDescriptor::calculateFromSmiles(const std::string& smiles) const {
    if (smiles.length() < 2) return 0.0;
    
    std::unordered_map<char, std::unordered_map<char, int>> transitionCounts;
    std::unordered_map<char, int> charCounts;
    
    for (size_t i = 0; i < smiles.length() - 1; ++i) {
        char current = smiles[i];
        char next = smiles[i + 1];
        transitionCounts[current][next]++;
        charCounts[current]++;
    }
    charCounts[smiles.back()]++;
    
    double entropy = 0.0;
    
    for (const auto& [from, transitions] : transitionCounts) {
        double fromProb = static_cast<double>(charCounts[from]) / smiles.length();
        double transitionEntropy = 0.0;
        
        for (const auto& [to, count] : transitions) {
            double transProb = static_cast<double>(count) / charCounts[from];
            if (transProb > 0) {
                transitionEntropy -= transProb * std::log2(transProb);
            }
        }
        
        entropy += fromProb * transitionEntropy;
    }
    
    return entropy;
}

std::variant<double, int, std::string> SequenceFractalDimensionDescriptor::calculateFromSmiles(const std::string& smiles) const {
    if (smiles.length() < 8) return 0.0;
    
    std::vector<int> boxCounts;
    std::vector<int> boxSizes = {2, 4, 8, 16};
    
    for (int boxSize : boxSizes) {
        if (static_cast<int>(smiles.length()) < boxSize) break;
        
        std::unordered_set<std::string> uniquePatterns;
        for (size_t i = 0; i <= smiles.length() - boxSize; ++i) {
            std::string pattern = smiles.substr(i, boxSize);
            uniquePatterns.insert(pattern);
        }
        
        boxCounts.push_back(uniquePatterns.size());
    }
    
    if (boxCounts.size() < 2) return 0.0;
    
    double sumX = 0.0, sumY = 0.0, sumXY = 0.0, sumX2 = 0.0;
    int n = boxCounts.size();
    
    for (size_t i = 0; i < boxCounts.size(); ++i) {
        double x = std::log(boxSizes[i]);
        double y = std::log(boxCounts[i]);
        
        sumX += x;
        sumY += y;
        sumXY += x * y;
        sumX2 += x * x;
    }
    
    double slope = (n * sumXY - sumX * sumY) / (n * sumX2 - sumX * sumX);
    return slope;
}

std::variant<double, int, std::string> SubsequenceRepetitionDescriptor::calculateFromSmiles(const std::string& smiles) const {
    if (smiles.length() < 4) return 0.0;
    
    std::unordered_map<std::string, int> subseqCounts;
    int totalSubseqs = 0;
    
    for (size_t len = 2; len <= std::min(smiles.length() / 2, size_t(4)); ++len) {
        for (size_t i = 0; i <= smiles.length() - len; ++i) {
            std::string subseq = smiles.substr(i, len);
            subseqCounts[subseq]++;
            totalSubseqs++;
        }
    }
    
    int repeatedSubseqs = 0;
    for (const auto& [subseq, count] : subseqCounts) {
        if (count > 1) {
            repeatedSubseqs += count;
        }
    }
    
    return static_cast<double>(repeatedSubseqs) / totalSubseqs;
}

std::variant<double, int, std::string> MaxRepeatedSubstringDescriptor::calculateFromSmiles(const std::string& smiles) const {
    if (smiles.length() < 4) return 0.0;
    
    int maxLength = 0;
    
    for (size_t len = 2; len <= smiles.length() / 2; ++len) {
        std::unordered_set<std::string> seen;
        bool foundRepeat = false;
        
        for (size_t i = 0; i <= smiles.length() - len; ++i) {
            std::string substring = smiles.substr(i, len);
            if (seen.find(substring) != seen.end()) {
                maxLength = len;
                foundRepeat = true;
                break;
            }
            seen.insert(substring);
        }
        
        if (!foundRepeat) break;
    }
    
    return static_cast<double>(maxLength);
}

std::variant<double, int, std::string> NonCarbonComplexityDescriptor::calculateFromSmiles(const std::string& smiles) const {
    if (smiles.empty()) return 0.0;
    
    double complexity = 0.0;
    
    for (size_t i = 0; i < smiles.length(); ++i) {
        char c = smiles[i];
        if (std::isupper(c) && c != 'C') {
            double weight = 1.0;
            if (i + 1 < smiles.length() && std::islower(smiles[i+1])) {
                weight = 1.5;  // Two-letter element
                i++;
            }
            complexity += weight;
        } else if (c == '[') {
            size_t closeBracket = smiles.find(']', i);
            if (closeBracket != std::string::npos) {
                std::string bracketContent = smiles.substr(i + 1, closeBracket - i - 1);
                if (!bracketContent.empty() && bracketContent[0] != 'C') {
                    complexity += 2.0;  // Non-carbon atom in bracket
                }
                i = closeBracket;
            }
        }
    }
    
    return complexity / smiles.length();
}

std::variant<double, int, std::string> SideChainIndexDescriptor::calculateFromSmiles(const std::string& smiles) const {
    if (smiles.empty()) return 0.0;
    
    int sideChainCount = 0;
    int depth = 0;
    int maxDepth = 0;
    
    for (char c : smiles) {
        if (c == '(') {
            depth++;
            sideChainCount++;
            maxDepth = std::max(maxDepth, depth);
        } else if (c == ')') {
            depth--;
        }
    }
    
    if (sideChainCount == 0) return 0.0;
    
    return static_cast<double>(maxDepth) * sideChainCount / smiles.length();
}

std::variant<double, int, std::string> RingBranchRatioDescriptor::calculateFromSmiles(const std::string& smiles) const {
    if (smiles.empty()) return 0.0;
    
    int ringMarkers = 0;
    int branchMarkers = 0;
    
    for (size_t i = 0; i < smiles.length(); ++i) {
        if (std::isdigit(smiles[i])) {
            ringMarkers++;
        } else if (smiles[i] == '%' && i + 2 < smiles.length() && std::isdigit(smiles[i+1]) && std::isdigit(smiles[i+2])) {
            ringMarkers++;
            i += 2;
        } else if (smiles[i] == '(' || smiles[i] == ')') {
            branchMarkers++;
        }
    }
    
    if (branchMarkers == 0) return ringMarkers > 0 ? 100.0 : 0.0;
    
    return static_cast<double>(ringMarkers) / branchMarkers;
}

std::variant<double, int, std::string> AromaticAliphaticRatioDescriptor::calculateFromSmiles(const std::string& smiles) const {
    if (smiles.empty()) return 0.0;
    
    int aromaticCount = 0;
    int aliphaticCount = 0;
    
    for (size_t i = 0; i < smiles.length(); ++i) {
        char c = smiles[i];
        
        if (c == 'c' || c == 'n' || c == 'o' || c == 's' || c == 'p') {
            aromaticCount++;
        } else if (std::isupper(c) && (c == 'C' || c == 'N' || c == 'O' || c == 'S' || c == 'P')) {
            aliphaticCount++;
            if (i + 1 < smiles.length() && std::islower(smiles[i+1])) {
                i++;
            }
        }
    }
    
    if (aliphaticCount == 0) return aromaticCount > 0 ? 100.0 : 0.0;
    
    return static_cast<double>(aromaticCount) / aliphaticCount;
}



// Specialized Chemical Pattern Features implementations
std::variant<double, int, std::string> FunctionalGroupEntropyDescriptor::calculateFromSmiles(const std::string& smiles) const {
    if (smiles.length() < 3) return 0.0;
    
    std::unordered_map<std::string, int> functionalGroups;
    int totalGroups = 0;
    
    // Simple pattern matching for common functional groups
    for (size_t i = 0; i < smiles.length() - 2; ++i) {
        std::string pattern;
        
        if (i + 3 <= smiles.length()) {
            pattern = smiles.substr(i, 3);
            if (pattern == "C=O" || pattern == "COO" || pattern == "COH" || 
                pattern == "C#N" || pattern == "C-N" || pattern == "N-C") {
                functionalGroups[pattern]++;
                totalGroups++;
            }
        }
        
        if (i + 2 <= smiles.length()) {
            pattern = smiles.substr(i, 2);
            if (pattern == "OH" || pattern == "NH" || pattern == "SH" || 
                pattern == "CO" || pattern == "CN" || pattern == "CS") {
                functionalGroups[pattern]++;
                totalGroups++;
            }
        }
    }
    
    if (totalGroups == 0) return 0.0;
    
    double entropy = 0.0;
    for (const auto& [group, count] : functionalGroups) {
        double p = static_cast<double>(count) / totalGroups;
        entropy -= p * std::log2(p);
    }
    
    return entropy;
}

std::variant<double, int, std::string> LargestRingDescriptor::calculateFromSmiles(const std::string& smiles) const {
    if (smiles.empty()) return 0.0;
    
    std::map<int, std::vector<size_t>> ringPositions;
    
    for (size_t i = 0; i < smiles.length(); ++i) {
        if (std::isdigit(smiles[i])) {
            int ringNumber = smiles[i] - '0';
            ringPositions[ringNumber].push_back(i);
        } else if (smiles[i] == '%' && i + 2 < smiles.length() && 
                  std::isdigit(smiles[i+1]) && std::isdigit(smiles[i+2])) {
            int ringNumber = (smiles[i+1] - '0') * 10 + (smiles[i+2] - '0');
            ringPositions[ringNumber].push_back(i);
            i += 2;
        }
    }
    
    int maxRingSize = 0;
    
    for (const auto& [ringNumber, positions] : ringPositions) {
        if (positions.size() == 2) {
            int size = std::abs(static_cast<int>(positions[1] - positions[0]));
            maxRingSize = std::max(maxRingSize, size);
        }
    }
    
    return static_cast<double>(maxRingSize) / 2.0;  // Rough estimate, dividing by 2
}

std::variant<double, int, std::string> CarbonSkeletonComplexityDescriptor::calculateFromSmiles(const std::string& smiles) const {
    if (smiles.empty()) return 0.0;
    
    int carbonCount = 0;
    int carbonChainTransitions = 0;
    bool prevWasCarbon = false;
    
    for (size_t i = 0; i < smiles.length(); ++i) {
        bool isCarbon = false;
        
        if (smiles[i] == 'C' || smiles[i] == 'c') {
            isCarbon = true;
            carbonCount++;
        } else if (smiles[i] == '[') {
            size_t closeBracket = smiles.find(']', i);
            if (closeBracket != std::string::npos) {
                std::string bracketContent = smiles.substr(i + 1, closeBracket - i - 1);
                if (bracketContent.length() > 0 && (bracketContent[0] == 'C' || bracketContent[0] == 'c')) {
                    isCarbon = true;
                    carbonCount++;
                }
                i = closeBracket;
            }
        }
        
        if (isCarbon != prevWasCarbon && i > 0) {
            carbonChainTransitions++;
        }
        
        prevWasCarbon = isCarbon;
    }
    
    if (carbonCount == 0) return 0.0;
    
    return static_cast<double>(carbonChainTransitions) / carbonCount;
}

std::variant<double, int, std::string> HeterocycleDescriptor::calculateFromSmiles(const std::string& smiles) const {
    if (smiles.empty()) return 0.0;
    
    int ringCount = 0;
    int heteroCycleCount = 0;
    std::map<int, std::vector<size_t>> ringPositions;
    std::map<int, bool> ringHasHeteroatom;
    
    for (size_t i = 0; i < smiles.length(); ++i) {
        // Track ring markers
        if (std::isdigit(smiles[i])) {
            int ringNumber = smiles[i] - '0';
            ringPositions[ringNumber].push_back(i);
        } else if (smiles[i] == '%' && i + 2 < smiles.length()) {
            int ringNumber = (smiles[i+1] - '0') * 10 + (smiles[i+2] - '0');
            ringPositions[ringNumber].push_back(i);
            i += 2;
        }
        
        // Check for heteroatoms
        if (std::isupper(smiles[i]) && smiles[i] != 'C') {
            // Look backward to see if we're in a ring context
            for (auto& [ringNum, positions] : ringPositions) {
                if (!positions.empty() && positions.back() < i && 
                    i - positions.back() < 6) {  // Arbitrary proximity check
                    ringHasHeteroatom[ringNum] = true;
                }
            }
        }
    }
    
    // Count rings and heterocycles
    for (const auto& [ringNum, positions] : ringPositions) {
        if (positions.size() >= 2) {
            ringCount++;
            if (ringHasHeteroatom[ringNum]) {
                heteroCycleCount++;
            }
        }
    }
    
    if (ringCount == 0) return 0.0;
    
    return static_cast<double>(heteroCycleCount) / ringCount;
}

// Continuing implementions for the remaining descriptors

std::variant<double, int, std::string> SubstructureBalanceDescriptor::calculateFromSmiles(const std::string& smiles) const {
    if (smiles.empty()) return 0.0;
    
    // Count some basic substructure types
    int aromatics = 0;
    int aliphatics = 0;
    int rings = 0;
    int branches = 0;
    int heteroatoms = 0;
    
    for (size_t i = 0; i < smiles.length(); ++i) {
        if (std::islower(smiles[i]) && smiles[i] != 'c') {
            aromatics++;
            heteroatoms++;
        } else if (smiles[i] == 'c') {
            aromatics++;
        } else if (std::isupper(smiles[i]) && smiles[i] != 'C') {
            aliphatics++;
            heteroatoms++;
        } else if (smiles[i] == 'C') {
            aliphatics++;
        } else if (std::isdigit(smiles[i])) {
            rings++;
        } else if (smiles[i] == '(' || smiles[i] == ')') {
            branches++;
        }
    }
    
    double total = aromatics + aliphatics + rings + branches + heteroatoms;
    if (total == 0) return 0.0;
    
    // Calculate balance as Shannon entropy of the distribution
    double balance = 0.0;
    double p;
    
    if (aromatics > 0) {
        p = aromatics / total;
        balance -= p * std::log2(p);
    }
    if (aliphatics > 0) {
        p = aliphatics / total;
        balance -= p * std::log2(p);
    }
    if (rings > 0) {
        p = rings / total;
        balance -= p * std::log2(p);
    }
    if (branches > 0) {
        p = branches / total;
        balance -= p * std::log2(p);
    }
    if (heteroatoms > 0) {
        p = heteroatoms / total;
        balance -= p * std::log2(p);
    }
    
    // Normalize to [0,1] - max entropy for 5 equal categories is log2(5)
    return balance / std::log2(5);
}

std::variant<double, int, std::string> MarkovEntropyDescriptor::calculateFromSmiles(const std::string& smiles) const {
    if (smiles.length() < 3) return 0.0;
    
    // Build first-order Markov model
    std::unordered_map<char, std::unordered_map<char, int>> transitions;
    std::unordered_map<char, int> stateCounts;
    
    for (size_t i = 0; i < smiles.length() - 1; ++i) {
        char current = smiles[i];
        char next = smiles[i + 1];
        transitions[current][next]++;
        stateCounts[current]++;
    }
    stateCounts[smiles.back()]++;
    
    // Calculate conditional entropy
    double entropy = 0.0;
    
    for (const auto& [state, count] : stateCounts) {
        if (count > 0) {
            double stateProb = static_cast<double>(count) / smiles.length();
            double stateEntropy = 0.0;
            
            for (const auto& [nextState, transCount] : transitions[state]) {
                double transProb = static_cast<double>(transCount) / count;
                if (transProb > 0) {
                    stateEntropy -= transProb * std::log2(transProb);
                }
            }
            
            entropy += stateProb * stateEntropy;
        }
    }
    
    return entropy;
}

std::variant<double, int, std::string> ContextualEntropyDescriptor::calculateFromSmiles(const std::string& smiles) const {
    if (smiles.length() < 4) return 0.0;
    
    // Build context-based probabilities with context window of 2
    std::unordered_map<std::string, std::unordered_map<char, int>> contextTransitions;
    std::unordered_map<std::string, int> contextCounts;
    
    for (size_t i = 0; i < smiles.length() - 2; ++i) {
        std::string context = smiles.substr(i, 2);
        char next = smiles[i + 2];
        contextTransitions[context][next]++;
        contextCounts[context]++;
    }
    
    // Calculate entropy from context to next symbol
    double entropy = 0.0;
    int totalContexts = 0;
    
    for (const auto& [context, count] : contextCounts) {
        double contextEntropy = 0.0;
        totalContexts += count;
        
        for (const auto& [next, nextCount] : contextTransitions[context]) {
            double prob = static_cast<double>(nextCount) / count;
            if (prob > 0) {
                contextEntropy -= prob * std::log2(prob);
            }
        }
        
        entropy += contextEntropy * (static_cast<double>(count) / (smiles.length() - 2));
    }
    
    return entropy;
}

std::variant<double, int, std::string> SyntacticComplexityDescriptor::calculateFromSmiles(const std::string& smiles) const {
    if (smiles.empty()) return 0.0;
    
    // Track different syntactic elements
    int atoms = 0;
    int bonds = 0;
    int rings = 0;
    int branches = 0;
    int stereo = 0;
    int charge = 0;
    int bracketAtoms = 0;
    int fragments = 1; // Start with 1 for the first fragment
    
    for (size_t i = 0; i < smiles.length(); ++i) {
        char c = smiles[i];
        
        if (std::isalpha(c)) {
            atoms++;
        } else if (c == '-' || c == '=' || c == '#' || c == ':') {
            bonds++;
        } else if (std::isdigit(c) || (c == '%' && i + 2 < smiles.length())) {
            rings++;
            if (c == '%') i += 2; // Skip the digits after %
        } else if (c == '(' || c == ')') {
            branches++;
        } else if (c == '/' || c == '\\' || c == '@') {
            stereo++;
        } else if (c == '+' || c == '-') {
            charge++;
        } else if (c == '[') {
            bracketAtoms++;
            // Skip to matching closing bracket
            size_t closeBracket = smiles.find(']', i);
            if (closeBracket != std::string::npos) {
                i = closeBracket;
            }
        } else if (c == '.') {
            fragments++;
        }
    }
    
    // Calculate complexity as weighted sum of different syntactic elements
    double complexity = 0.0;
    complexity += atoms * 1.0;
    complexity += bonds * 1.2;
    complexity += rings * 1.5;
    complexity += branches * 1.3;
    complexity += stereo * 1.8;
    complexity += charge * 1.4;
    complexity += bracketAtoms * 1.6;
    complexity += fragments * 1.1;
    
    // Normalize by a typical maximum - rough estimate
    return std::min(1.0, complexity / (smiles.length() * 2.0));
}

std::variant<double, int, std::string> RelativeEntropyDescriptor::calculateFromSmiles(const std::string& smiles) const {
    if (smiles.empty()) return 0.0;
    
    // Reference distribution for typical SMILES (simplified approximation)
    // Based on frequency analysis of common SMILES strings
    std::unordered_map<char, double> refProbs = {
        {'C', 0.25}, {'c', 0.12}, {'(', 0.05}, {')', 0.05},
        {'O', 0.07}, {'N', 0.06}, {'=', 0.04}, {'1', 0.03},
        {'2', 0.03}, {'3', 0.02}, {'[', 0.02}, {']', 0.02},
        {'S', 0.01}, {'P', 0.01}, {'F', 0.01}, {'-', 0.01},
        {'#', 0.01}, {':', 0.01}, {'.', 0.01}, {'/', 0.01},
        {'\\', 0.01}, {'@', 0.01}, {'+', 0.01}, {'n', 0.01},
        {'o', 0.01}, {'s', 0.01}, {'H', 0.01}, {'B', 0.01}
    };
    
    // Calculate actual distribution
    std::unordered_map<char, int> charCounts;
    for (char c : smiles) {
        charCounts[c]++;
    }
    
    // Calculate KL divergence
    double kl = 0.0;
    for (const auto& [c, count] : charCounts) {
        double p = static_cast<double>(count) / smiles.length();
        double q = refProbs.count(c) ? refProbs[c] : 0.001; // Small value for unseen chars
        
        if (p > 0 && q > 0) {
            kl += p * std::log2(p / q);
        }
    }
    
    // Normalize to [0,1] - typical range is 0-5 for KL divergence of SMILES
    return std::min(1.0, kl / 5.0);
}

std::variant<double, int, std::string> InformationDensityDescriptor::calculateFromSmiles(const std::string& smiles) const {
    if (smiles.empty()) return 0.0;
    
    // Calculate character frequencies
    std::unordered_map<char, int> charCounts;
    for (char c : smiles) {
        charCounts[c]++;
    }
    
    // Calculate Shannon entropy
    double entropy = 0.0;
    for (const auto& [c, count] : charCounts) {
        double p = static_cast<double>(count) / smiles.length();
        if (p > 0) {
            entropy -= p * std::log2(p);
        }
    }
    
    // Calculate information density based on entropy and semantics
    double semanticWeight = 0.0;
    for (size_t i = 0; i < smiles.length(); ++i) {
        char c = smiles[i];
        
        // Assign weights to different character types based on their semantic importance
        if (std::isalpha(c)) {
            semanticWeight += 1.2; // Element symbols
        } else if (c == '[' || c == ']') {
            semanticWeight += 0.8; // Bracket atoms
        } else if (c == '-' || c == '=' || c == '#') {
            semanticWeight += 1.0; // Bond symbols
        } else if (c == '(' || c == ')') {
            semanticWeight += 0.9; // Branch symbols
        } else if (std::isdigit(c) || c == '%') {
            semanticWeight += 0.7; // Ring markers
        } else if (c == '/' || c == '\\' || c == '@') {
            semanticWeight += 1.1; // Stereo symbols
        } else if (c == '+' || c == '-') {
            semanticWeight += 0.9; // Charges
        } else {
            semanticWeight += 0.5; // Other characters
        }
    }
    
    double density = (entropy * semanticWeight) / (smiles.length() * 1.2);
    
    // Normalize to [0,1] - typical max value is around 8-10
    return std::min(1.0, density / 8.0);
}

}  // namespace descriptors
}  // namespace desfact