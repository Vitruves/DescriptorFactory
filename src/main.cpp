#include "descriptors.hpp"
#include <cxxopts.hpp>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <chrono>
#include <atomic>
#include <thread> 
#include <random>
#include <utility>
#include <unordered_set>
#include <sys/ioctl.h> // For terminal width detection
#include <unistd.h>    // For STDOUT_FILENO
#include "utils.hpp"
#include "io.hpp"
#include <regex> // Required for potential future coloring


#ifdef WITH_TBB
#include <tbb/tbb.h>
#include <tbb/task_group.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_pipeline.h>
#include <tbb/concurrent_queue.h>
#include <tbb/flow_graph.h>
#else
#warning "TBB not found, multiprocessing disabled. Falling back to single-threaded execution."
#ifndef WITH_TBB
#define WITH_TBB 0
#endif
#endif

using namespace desfact;

// Helper function to add color to option group titles
std::string colorizeHelp(const std::string& helpText) {
    // Basic coloring example: Make group titles bold cyan
    // This is a simple approach; more sophisticated coloring might require
    // parsing the help string more carefully.
    std::regex groupRegex("^[A-Za-z /]+:$", std::regex_constants::multiline);
    std::string coloredHelp = std::regex_replace(helpText, groupRegex, "\033[1;36m$&\033[0m");
    return coloredHelp;
}

void printVersion() {
    std::cout << "\033[1;36mDescriptorFactory\033[0m (\033[1mdesfact\033[0m) v0.1.0" << std::endl;
}

void printHelp(const cxxopts::Options& options) {
    std::cout << colorizeHelp(options.help()) << std::endl;
}

// Add a simple file format checker
void inspectCsvFile(const std::string& filePath, const std::string& expectedSmilesColumn, bool hasHeader) {
    std::ifstream file(filePath);
    if (!file.is_open()) {
        globalLogger.error("Failed to open file for inspection: " + filePath);
        return;
    }
    
    std::string line;
    int lineCount = 0;
    
    // Read the first 5 lines to inspect the format
    while (lineCount < 5 && std::getline(file, line)) {
        if (lineCount == 0 && hasHeader) {
            globalLogger.info("Header line: " + line);
            // Check for SMILES column in header
            std::istringstream iss(line);
            std::string token;
            bool found = false;
            int index = 0;
            while (std::getline(iss, token, ',')) {
                if (token == expectedSmilesColumn) {
                    globalLogger.info("Found SMILES column at index " + std::to_string(index));
                    found = true;
                    break;
                }
                index++;
            }
            if (!found) {
                globalLogger.error("SMILES column '" + expectedSmilesColumn + "' not found in header");
            }
        } else {
            globalLogger.info("Data line " + std::to_string(lineCount) + ": " + line);
        }
        lineCount++;
    }
    
    file.close();
}

std::vector<std::string> extractSmiles(const std::vector<std::string>& lines, const CsvIO& csvHandler) {
    std::vector<std::string> smilesStrings;
    smilesStrings.reserve(lines.size());
    
    for (const auto& line : lines) {
        std::vector<std::string> cells = CsvIO::parseCsvLine(line, csvHandler.getDelimiter());
        int smilesIndex = csvHandler.getSmilesIndex();
        
        if (smilesIndex >= 0 && cells.size() > smilesIndex) {
            std::string smiles = cells[smilesIndex];
            if (!smiles.empty()) {
                smilesStrings.push_back(smiles);
            }
        }
    }
    return smilesStrings;
}

// Replace the countLinesInFileMT function with this safer version
size_t countLinesInFileMT(const std::string& filePath, int numThreads) {
    uintmax_t fileSize = 0;
    try {
        // Explicitly construct path object
        fileSize = std::filesystem::file_size(std::filesystem::path(filePath));
    } catch (const std::filesystem::filesystem_error& e) {
         globalLogger.error("Failed to get file size for line counting: " + filePath + " - " + e.what());
         return 0; // Cannot count if size unknown
    }

    if (fileSize == 0) return 0;

    size_t chunkSize = 4 * 1024 * 1024;
    size_t numChunks = (fileSize + chunkSize - 1) / chunkSize;
    numChunks = std::min(numChunks, static_cast<size_t>(numThreads * 4));
    if (numChunks == 0) numChunks = 1; // Ensure at least one chunk

    std::vector<size_t> chunkBoundaries;
    chunkBoundaries.reserve(numChunks + 1);

    for (size_t i = 0; i <= numChunks; i++) {
        // Use uintmax_t for intermediate calculation to avoid overflow
        chunkBoundaries.push_back(std::min(static_cast<size_t>(i * static_cast<uintmax_t>(fileSize) / numChunks), static_cast<size_t>(fileSize)));
    }

    std::atomic<size_t> lineCount{0};

    #ifdef WITH_TBB
    tbb::parallel_for(size_t(0), numChunks, [&](size_t chunkIndex) {
        std::ifstream threadFile(filePath, std::ios::binary);
        if (!threadFile.is_open()) return; // Check if file opened successfully

        size_t start = chunkBoundaries[chunkIndex];
        size_t end = chunkBoundaries[chunkIndex + 1];
        std::streampos currentPos = start; // Use streampos

        if (chunkIndex > 0) {
            threadFile.seekg(start); // Position before reading
            char c;
            // Correctly use tellg() which returns streampos
            while (threadFile.get(c) && c != '\n') {
                 currentPos = threadFile.tellg(); // Update position
                 if (currentPos >= static_cast<std::streampos>(end)) break; // Stop if we go past the chunk end
            }
            // If we found newline and didn't go past end, update start
            if (threadFile && c == '\n' && currentPos < static_cast<std::streampos>(end)) {
                start = static_cast<size_t>(currentPos);
            } else if (currentPos >= static_cast<std::streampos>(end)) {
                 // If seeking the boundary went past the end, this chunk is empty or malformed
                 return;
            } else {
                 // If we hit EOF or error before finding newline, reset start effectively skipping partial line
                 start = chunkBoundaries[chunkIndex]; // Or decide to return if file seems corrupt
            }
        }


        if (start >= end) return; // Check again after potential boundary adjustment

        threadFile.seekg(start); // Position at adjusted start

        std::string line;
        size_t localCount = 0;
        std::streampos readPos = start;

        while (threadFile.good() && readPos < static_cast<std::streampos>(end)) {
            if (std::getline(threadFile, line)) {
                 readPos = threadFile.tellg(); // Update position after getline
                 // Only count if the line *ended* before or at the chunk boundary
                 // Note: tellg() might be -1 on failure or EOF
                 if (readPos == std::streampos(-1) || readPos <= static_cast<std::streampos>(end)) {
                     localCount++;
                 } else {
                     // The line crossed the boundary, don't count it in this chunk
                     break;
                 }
            } else {
                 break; // getline failed (EOF or error)
            }
        }

        lineCount.fetch_add(localCount);
        threadFile.close(); // Close file handle for this thread
    });
    #else
    std::ifstream file(filePath);
    if (!file.is_open()) {
         globalLogger.error("Failed to open file for single-threaded line counting: " + filePath);
         return 0;
    }
    std::string line;
    size_t count = 0;
    while (std::getline(file, line)) {
        count++;
    }
    lineCount = count;
    #endif

    return lineCount;
}

#ifdef WITH_TBB
// Define token types for pipeline stages
struct BatchLines {
    std::vector<std::string> lines;
};

struct SmilesData {
    std::vector<std::string> smiles;
    std::vector<std::vector<std::string>> originalData;
};

struct BatchWithResults {
    std::shared_ptr<MoleculeBatch> batch;
    std::vector<std::vector<std::string>> originalData;
    std::unordered_map<std::string, std::vector<std::variant<double, int, std::string>>> results;
};
#endif

int main(int argc, char* argv[]) {
    // Parse command line arguments
    cxxopts::Options options("\033[1;36mdesfact\033[0m", 
                           "\033[1;37mChemical descriptor calculator for molecular structures\033[0m");
    
    // Organize options into groups with most important first
    options.add_options("Input")
        ("i,input", "Input CSV file path", cxxopts::value<std::string>())
        ("o,output", "Output CSV file path", cxxopts::value<std::string>())
        ("s,smiles", "Direct SMILES string input (alternative to input file)", cxxopts::value<std::string>())
        ("d,descriptors", "Comma-separated list of descriptors to calculate, or 'all'/'fast'/'no-rdkit' ('fast' excludes eigen, image, and morgan descriptors, 'no-rdkit' excludes RDKit-specific descriptors)", cxxopts::value<std::string>()->default_value("all"));
    
    options.add_options("CSV")
        ("sc,smiles-col", "Name/index of the column containing SMILES", cxxopts::value<std::string>()->default_value("SMILES"))
        ("delimiter", "CSV delimiter character", cxxopts::value<std::string>()->default_value(","))
        ("no-header", "Input CSV file has no header")
        ("escapechar", "CSV escape character (disables standard \"\" quoting)", cxxopts::value<std::string>()->default_value("\\"));
    
    options.add_options("Performance")
        ("b,batch-size", "Number of molecules to process per batch", cxxopts::value<size_t>())
        ("j,jobs", "Number of parallel jobs (0=auto)", cxxopts::value<int>())
        ("no-stream", "Process molecules in batches instead of streaming (increases memory usage)")
        ("verbose", "Enable detailed logging output");
        
    options.add_options("General")
        ("h,help", "Display help information")
        ("v,version", "Display version information")
        ("l,list", "List available descriptors");
    
    // Set positional help
    options.positional_help("\033[1m<input file>\033[0m \033[1m<output file>\033[0m");
    
    // Adjust width for potentially better alignment. Increase if needed.
    options.set_width(120); 
    options.custom_help("[OPTION...]"); // Add standard usage indicator
    
    // Fast help/version display
    if (argc == 1 || (argc > 1 && (strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0))) {
        printHelp(options);
        return 0;
    }
    if (argc > 1 && (strcmp(argv[1], "-v") == 0 || strcmp(argv[1], "--version") == 0)) {
        printVersion();
        return 0;
    }

    try {
        auto result = options.parse(argc, argv);

        // Handle help/version again after parsing in case they are not the first arg
        if (result.count("help")) { printHelp(options); return 0; }
        if (result.count("version")) { printVersion(); return 0; }

        // Configure global settings
        globalConfig.numThreads = result.count("jobs") ? result["jobs"].as<int>() : 0;
        globalConfig.verbose = result.count("verbose") > 0;
        globalConfig.batchSize = result.count("batch-size") ? result["batch-size"].as<size_t>() : 64;
        
        // Set default minimal level based on verbose flag
        if (globalConfig.verbose) {
            globalLogger.setMinLevel(LogLevel::DEBUG);
            globalLogger.info("Verbose mode enabled.");
            globalLogger.info("Using " + std::to_string(globalConfig.numThreads) + " parallel jobs.");
            globalLogger.info("Processing in batches of " + std::to_string(globalConfig.batchSize));
        } else {
            globalLogger.setMinLevel(LogLevel::WARNING); // Only show warnings and errors by default
        }

        // Automatically determine optimal job count if not specified
        if (globalConfig.numThreads <= 0) {
            // Use available CPU cores minus 1 to leave a core for UI/system
            int availableCores = std::thread::hardware_concurrency();
            globalConfig.numThreads = availableCores > 1 ? availableCores - 1 : 1;
            globalLogger.info("Auto-configured to use " + std::to_string(globalConfig.numThreads) + " parallel jobs.");
        }

        // Ensure TBB thread settings are applied
        #ifdef WITH_TBB
        tbb::global_control global_limit(
            tbb::global_control::max_allowed_parallelism, 
            globalConfig.numThreads
        );
        globalLogger.debug("TBB configured with " + std::to_string(globalConfig.numThreads) + " max parallel jobs");
        #endif

        // Initialize descriptor factory
        DescriptorFactory factory;

        if (result.count("list")) {
            auto descriptors = factory.getAvailableDescriptors();
            std::cout << "\033[1;36mAvailable descriptors:\033[0m" << std::endl;
            
            // Get terminal width for better formatting
            size_t termWidth = 80; // Default width
            try {
                struct winsize w;
                if (ioctl(STDOUT_FILENO, TIOCGWINSZ, &w) != -1) {
                    termWidth = w.ws_col;
                }
            } catch (...) {
                // Fallback to default if we can't get terminal width
            }
            
            // Calculate the description column start position and max width
            const size_t nameColumn = 20; // Width allocated for descriptor name column
            const size_t indentWidth = 4;  // Spaces for continuation lines
            const size_t maxDescWidth = termWidth - nameColumn - 2;
            
            for (const auto& desc : descriptors) {
                const auto* descriptor = factory.getDescriptor(desc);
                std::string description = descriptor->getDescription();
                
                // Print the descriptor name with color
                std::cout << "\033[1;32m" << std::left << std::setw(nameColumn) << desc << "\033[0m";
                
                // Handle description wrapping based on terminal width
                if (description.empty()) {
                    std::cout << std::endl;
                    continue;
                }
                
                // Print first line
                size_t printedChars = 0;
                size_t descLength = description.length();
                
                if (descLength <= maxDescWidth) {
                    // Short description fits on one line
                    std::cout << description << std::endl;
                } else {
                    // First line with descriptor name
                    std::cout << description.substr(0, maxDescWidth) << std::endl;
                    printedChars += maxDescWidth;
                    
                    // Remaining lines indented
                    while (printedChars < descLength) {
                        std::cout << std::string(nameColumn + indentWidth, ' ');
                        size_t charsToTake = std::min(maxDescWidth - indentWidth, descLength - printedChars);
                        std::cout << description.substr(printedChars, charsToTake) << std::endl;
                        printedChars += charsToTake;
                    }
                }
            }
            return 0;
        }

        // Check if we have direct SMILES or input file
        bool directSmilesMode = result.count("smiles") > 0;
        bool hasInputFile = result.count("input") > 0;
        bool hasOutputFile = result.count("output") > 0;
        bool hasDescriptors = result.count("descriptors") > 0;

        // Validate required arguments
        if (!directSmilesMode && (!hasInputFile || !hasOutputFile || !hasDescriptors)) {
            std::cerr << "\033[1;31mError:\033[0m Either --smiles or (--input, --output, and --descriptors) are required." << std::endl;
            printHelp(options);
            return 1;
        }

        if (directSmilesMode && !hasDescriptors) {
            std::cerr << "\033[1;31mError:\033[0m --descriptors is required when using --smiles." << std::endl;
            printHelp(options);
            return 1;
        }

        // Process descriptor list
        std::string descriptorsStr = result["descriptors"].as<std::string>();
        std::vector<std::string> descriptorNames;

        // Convert descriptorsStr to lowercase for case-insensitive comparison
        std::string descriptorsStrLower = descriptorsStr;
        std::transform(descriptorsStrLower.begin(), descriptorsStrLower.end(), descriptorsStrLower.begin(),
                       [](unsigned char c){ return std::tolower(c); });

        if (descriptorsStrLower == "all") {
            descriptorNames = factory.getAvailableDescriptors();
            globalLogger.info("Using all " + std::to_string(descriptorNames.size()) + " available descriptors.");
        } else if (descriptorsStrLower == "fast") {
            // Get all available descriptors first
            std::vector<std::string> allDescriptors = factory.getAvailableDescriptors();
            
            // Create sets of descriptors to exclude from eigen.cpp, image.cpp, and morgan.cpp
            std::unordered_set<std::string> excludeDescriptors = {
                // Morgan descriptors
                "MorganBitDensityRatio", "MorganFragUniqueness", "MorganRadiusInfoRatio", "MorganBitClustering", "MorganFragEntropy", "MorganFPAsymmetry", "MorganBitTransitionRate", "MorganFragSizeSkew", "MorganLCSScore", "MorganPharmacophoreDensity", "MorganFragDiversity", "MorganBitCorrelation", "MorganPatternRecurrence", "MorganSimilarityToRef", "MorganFragComplexity", "MorganInfoDensity", "MorganEnvVariability", "MorganBitImportance", "MorganRingSysScore", "MorganBitProbEntropy", "MorganFragENSpectrum", "MorganBitCorrMatrix", "MorganFragConnectivity", "MorganBitFreqSkew", "MorganSubstructHeterogeneity", "MorganBitPolarityDist", "MorganFragSimNetwork", "MorganBitInfoGain", "MorganSubstructDivGradient", "MorganFPSymmetry",
                
                // Image descriptors
                "BoundingBoxArea", "MoleculeWidth", "MoleculeHeight", "AspectRatio", "AtomDensity", "MeanAtomRadius", "MaxAtomRadius", "MinAtomRadius", "MeanBondLength", "MaxBondLength", "MinBondLength", "BondLengthStd", "MeanAtomAtomDist", "MeanAtomLuminance", "MaxAtomAtomDist", "MinAtomAtomDist", "MeanAtomArea", "MedianAtomRadius", "MedianBondLength", "AtomRadiusRange", "AtomAreaFraction", "BondCoverageFraction", "AtomPackingDensity", "AtomMassX", "AtomMassY", "AtomMassDist", "AtomRadiusStd", "MeanBondAngle", "BondAngleStd", "AtomXStd", "AtomYStd", "AtomLuminanceStd", "BondLenMAD", "AtomRadiusMAD", "AtomAreaStd", "AtomRadiusCV", "AtomAreaCV", "AtomLuminanceCV", "AtomAreaRange", "AtomAreaMedian", "AtomAreaMAD", "AtomAreaHullFrac", "AtomAreaCenterFrac", "StdAllAtomDist", "MinBondAngle", "MaxBondAngle", "MedianBondAngle", "MeanBondColorDiff", "MeanBondLuminanceDiff", "MeanBondLenDiffColor", "MeanBondAngleDiffColor", "MeanDistSameColor", "MeanDistDiffColor", "MeanBondOrder", "FracDoubleBonds", "FracTripleBonds", "FracAromaticBonds", "MeanAtomDegree", "MaxAtomDegree", "MinAtomDegree", "MeanDistToCentroid", "StdDistToCentroid", "MeanLenSingle", "MeanLenDouble", "MeanLenTriple", "MeanLenAromatic", "FracDeg1", "FracDeg2", "FracDeg3", "FracDeg4", "MeanBondRadiiDiff", "MeanBondLuminanceDiff2", "MeanAngleHighDegree", "BoundaryAtomRatio", "MeanAtomEccentricity", "MolecularDiameter", "MolecularRadius", "RadiusOfGyration", "MolecularSphericity", "BondLenToAtomRadiusRatio", "EdgeDensity", "PlanarityMeasure", "MeanNearestNeighborDist", "MeanAtomsInRadius", "VarNearestNeighborDist", "MeanBondToBondAngle", "MomentOfInertia", "MeanAtomCentrality", "TerminalAtomCount", "JunctionAtomCount", "TerminalToJunctionRatio", "WienerIndex", "RingCount", "AcidicCenterCount", "BasicCenterCount", "AcidicToBasicRatio", "MeanAcidicAcidicDist", "MeanBasicBasicDist", "MeanAcidicBasicDist", "MinAcidicBasicDist", "MeanAtomsNearAcidic", "MeanAtomsNearBasic", "MeanAcidicRadius", "MeanBasicRadius", "MeanAcidicCentroidDist", "MeanBasicCentroidDist", "FracAcidicOnHull", "FracBasicOnHull", "MeanAcidicAngle", "MeanBasicAngle", "MeanAcidicLuminance", "MeanBasicLuminance", "AcidicBasicLuminanceDiff", "FractalDimension", "ForegroundRatio", "AverageColor", "ColorVariance", "ImageCenterX", "ImageCenterY", "ImageOrientation",
                
                // Eigen descriptors
                "AdjNonZero", "AdjTrace", "AdjFrobNorm", "AdjSpecRad", "AdjMinEig", "AdjEigSum", "AdjEigSumSq", "AdjEigVar", "AdjEigSkew", "AdjEigKurt", "AdjPosEig", "AdjNegEig", "AdjZeroEig", "AdjMaxDeg", "AdjMinDeg", "AdjMeanDeg", "AdjDegVar", "AdjDegStd", "AdjRank", "GraphEnergy", "LapSpecRad", "LapAlgConn", "LapZeroEig", "LapEnergy", "LapTrace", "LapDet", "LapTotEffRes", "LapKirchhoff", "LapEigVar", "LapEigSkew", "NormLapSpecRad", "NormLapMinNonZero", "NormLapMaxEig", "NormLapEnergy", "NormLapTrace", "DegMatMaxDeg", "DegMatMinDeg", "DegMatAvgDeg", "DegMatVar", "DegMatSkew", "DegMatEnt", "GraphIrreg", "WienerIdx", "EstradaIdx", "NumSpanTrees", "GraphEcc", "SpecGap", "TraceMatPow2", "NumTriangles", "GraphDiam", "Num2Walks", "Num3Walks", "Num4Walks", "MeanClosed3Walks", "MeanClosed4Walks", "Walk2Energy", "Walk3Energy", "Walk4Energy", "GraphIrregWalk", "TriPathRatio", "MaxSingVal", "MinNonZeroSingVal", "CondNumber", "SumSingVal", "FrobNormSingVal", "SingValEntropy", "SingValVar", "SingValSkew", "SpecEffRank", "NuclearNorm", "NormAdjSpecRad", "NormEigGap", "SumNormEig", "VarNormEig", "CountNormEigAboveHalf", "NormEnergy", "MaxNormEigCent", "AvgNormEigCent", "NormAdjRank", "NormAdjEntropy", "SignLapSpecRad", "SignLapMinEig", "SignLapEnergy", "SignLapTrace", "SignLapDet", "SignLapZeroEig", "SignLapPosEig", "SignLapNegEig", "SignLapEigVar", "SignLapEigSkew", "WtAdjSpecRad", "MeanFirstPassage", "CommuteTime", "KirchhoffVar", "EffGraphResDist", "LocalClustCoefDist", "GraphRobustIdx", "NormEstradaIdx", "GraphBipart", "SpanTreeEntropy", "DistMatSpecRad", "DistMatMinEig", "DistMatEigSum", "DistMatEigVar", "DistMatEigSkew", "DistMatEigKurt", "DistMatEnergy", "DistMatFrobNorm", "DistMatRank", "DistMatDet", "DistMatConditionNumber", "DistMatSingValEntropy", "DistMatSpectralEffRank", "WtNumAdjSpecRad", "WtNumAdjMinEig", "WtNumAdjEigSum", "WtNumAdjEigVar", "WtNumAdjEnergy", "WtNumLapSpecRad", "WtNumLapAlgConn", "WtNumLapEigSum", "WtNumLapEigVar", "WtNumLapEnergy", "AdjPow2SpecRad", "AdjPow2Energy", "AdjPow3SpecRad", "AdjPow3Energy", "AdjPow4SpecRad", "AdjPow4Energy", "LapPow2SpecRad", "LapPow2AlgConn", "LapPow2Energy", "LapPow2Trace", "DistMatEstradaIdx", "WtNumAdjEstradaIdx", "WtNumLapEstradaIdx", "AdjDomEigVecVar", "AdjDomEigVecSkew", "AdjDomEigVecKurt", "NormLapDomEigVecVar", "SignLapDomEigVecVar", "SumFormanRicci", "LogDetLaplacian", "AvgDegreeDistance", "VarianceDegreeDistance"
            };
            
            // Filter out excluded descriptors
            for (const auto& descriptor : allDescriptors) {
                if (excludeDescriptors.find(descriptor) == excludeDescriptors.end()) {
                    descriptorNames.push_back(descriptor);
                }
            }
            
            globalLogger.info("Using 'fast' mode with " + std::to_string(descriptorNames.size()) + 
                              " descriptors (excluding eigen, image, and morgan descriptors).");
        } else if (descriptorsStrLower == "no-rdkit") {
            // Get all available descriptors first
            std::vector<std::string> allDescriptors = factory.getAvailableDescriptors();
            
            // Create set of RDKit descriptors to exclude
            std::unordered_set<std::string> rdkitDescriptors = {
                "rdkitMW", "rdkitTPSA", "rdkitLogP", "rdkitMR", "rdkitLipinskiHBA", "rdkitLipinskiHBD",
                "rdkitRotatableBonds", "rdkitHBondAcceptors", "rdkitHBondDonors", "rdkitHeavyAtoms",
                "rdkitNumRings", "rdkitAromaticRings", "rdkitAliphaticRings", "rdkitFractionCsp3",
                "rdkitChiralCenters", "rdkitLabuteASA", "rdkitHallKierAlpha", "rdkitKappa1", 
                "rdkitKappa2", "rdkitKappa3", "rdkitChi0v", "rdkitChi1v", "rdkitChi2v", 
                "rdkitChi3v", "rdkitChi4v", "rdkitMurckoCount", "rdkitRingCount5", 
                "rdkitSpiroAtomCount", "rdkitBridgeheadAtomCount"
            };
            
            // Filter out RDKit descriptors
            for (const auto& descriptor : allDescriptors) {
                if (rdkitDescriptors.find(descriptor) == rdkitDescriptors.end()) {
                    descriptorNames.push_back(descriptor);
                }
            }
            
            globalLogger.info("Using 'no-rdkit' mode with " + std::to_string(descriptorNames.size()) + 
                             " descriptors (excluding RDKit-specific descriptors).");
        } else {
            std::stringstream ss(descriptorsStr);
            std::string item;
            while (std::getline(ss, item, ',')) {
                 // Trim whitespace
                item.erase(0, item.find_first_not_of(" \t\n\r\f\v"));
                item.erase(item.find_last_not_of(" \t\n\r\f\v") + 1);
                if (!item.empty()) {
                    if (factory.getDescriptor(item)) {
                        descriptorNames.push_back(item);
                    } else {
                        globalLogger.warning("Unknown descriptor specified and ignored: " + item);
                    }
                }
            }
            if (descriptorNames.empty()) {
                globalLogger.error("No valid descriptors specified.");
                return 1;
            }
            globalLogger.info("Using " + std::to_string(descriptorNames.size()) + " specified descriptors.");
        }

        if (directSmilesMode) {
            // Handle direct SMILES input
            std::string smilesStr = result["smiles"].as<std::string>();
            std::string outputPath = hasOutputFile ? result["output"].as<std::string>() : "";
            
            // Process single SMILES
            Molecule molecule(smilesStr);
            if (!molecule.isValid()) {
                globalLogger.error("Invalid SMILES string: " + smilesStr);
                return 1;
            }
            
            std::cout << "\033[1;36mProcessing SMILES:\033[0m " << smilesStr << std::endl;
            
            // Calculate descriptors
            std::unordered_map<std::string, std::vector<std::variant<double, int, std::string>>> results;
            for (const auto& descriptorName : descriptorNames) {
                auto result = factory.calculate(descriptorName, molecule);
                results[descriptorName].push_back(result);
            }
            
            // Output results
            if (hasOutputFile) {
                // Write to CSV file
                std::ofstream outFile(outputPath);
                if (!outFile.is_open()) {
                    globalLogger.error("Failed to open output file: " + outputPath);
                    return 1;
                }
                
                // Write header
                outFile << "SMILES";
                for (const auto& desc : descriptorNames) {
                    outFile << "," << desc;
                }
                outFile << std::endl;
                
                // Write data
                outFile << smilesStr;
                for (const auto& desc : descriptorNames) {
                    outFile << ",";
                    const auto& value = results[desc][0];
                    // Output value based on variant type
                    std::visit([&outFile](auto&& arg) {
                        using T = std::decay_t<decltype(arg)>;
                        if constexpr (std::is_same_v<T, double>)
                            outFile << std::fixed << std::setprecision(6) << arg;
                        else if constexpr (std::is_same_v<T, int>)
                            outFile << arg;
                        else if constexpr (std::is_same_v<T, std::string>)
                            outFile << arg;
                    }, value);
                }
                outFile << std::endl;
                outFile.close();
                
                std::cout << "\033[1;32m✓\033[0m Results written to \033[1m" << outputPath << "\033[0m" << std::endl;
            } else {
                // Print to console
                std::cout << "\033[1;36mResults:\033[0m" << std::endl;
                for (const auto& desc : descriptorNames) {
                    std::cout << "\033[1m" << desc << ":\033[0m ";
                    const auto& value = results[desc][0];
                    // Output value based on variant type
                    std::visit([](auto&& arg) {
                        using T = std::decay_t<decltype(arg)>;
                        if constexpr (std::is_same_v<T, double>)
                            std::cout << std::fixed << std::setprecision(6) << arg;
                        else if constexpr (std::is_same_v<T, int>)
                            std::cout << arg;
                        else if constexpr (std::is_same_v<T, std::string>)
                            std::cout << arg;
                    }, value);
                    std::cout << std::endl;
                }
            }
            
            return 0;
        } else {
            // Continue with regular file processing
            const std::string inputPath = result["input"].as<std::string>();
            const std::string outputPath = result["output"].as<std::string>();
            std::string smilesColumn = result.count("smiles-col") ? result["smiles-col"].as<std::string>() : "SMILES";
            std::string delimiter = result.count("delimiter") ? result["delimiter"].as<std::string>() : ",";
            bool hasHeader = !result.count("no-header");
            std::string escapeChar = result.count("escapechar") ? result["escapechar"].as<std::string>() : "\\";
            
            // Use explicit path constructor
            if (!std::filesystem::exists(std::filesystem::path(inputPath))) {
                globalLogger.error("Input file does not exist: " + inputPath);
                return 1;
            }

            globalLogger.info("Processing input file: " + inputPath);
            globalLogger.info("Output will be written to: " + outputPath);

            // Move essential information to standard output
            if (!globalConfig.verbose) {
                std::cout << "\033[1;36mProcessing:\033[0m " << inputPath << " → " << outputPath << std::endl;
            }

            // We need to handle smilesColumn as either a numeric index or a column name
            int smilesColumnIndex = -1;
            std::string smilesColumnName = "SMILES"; // Default column name

            if (smilesColumn.size() > 0) {
                try {
                    // If smilesColumn is numeric, convert to int directly
                    smilesColumnIndex = std::stoi(smilesColumn);
                } catch (const std::invalid_argument&) {
                    // If it's not a number, it's a column name - we'll keep it as -1
                    // and use the column name
                    smilesColumnIndex = -1;
                    smilesColumnName = smilesColumn; // Set the column name
                    globalLogger.info("Looking for SMILES column named '" + smilesColumnName + "' in the header");
                } catch (const std::exception& e) {
                    globalLogger.error("Error processing SMILES column specifier: " + std::string(e.what()));
                    return 1;
                }
            }

            // Now create CsvIO with the proper types
            CsvIO csvInputHandler(inputPath, outputPath, delimiter, smilesColumnIndex, hasHeader, escapeChar);

            // We also need to tell CsvIO the name of the column if it's not an index
            if (smilesColumnIndex == -1 && hasHeader) {
                // Set the column name BEFORE parsing the header
                csvInputHandler.setSmilesColumnName(smilesColumnName);
                
                // Re-parse the header with the new column name
                if (!csvInputHandler.reparseHeader()) {
                    globalLogger.fatal("Cannot proceed: SMILES column '" + smilesColumnName + "' not found in input file header.");
                    return 1;
                }
            }

            // Output handler uses its own output path
            CsvIO::ResultWriter resultWriter = csvInputHandler.createResultWriter(descriptorNames, csvInputHandler.getEscapeChar());

            CsvIO::LineReader lineReader = csvInputHandler.createLineReader();

            // Initialize progress bar using LineReader's estimation
            size_t estimatedLines = 0;
            try {
                if (globalConfig.verbose) {
                    globalLogger.info("Estimating lines in input file...");
                }
                estimatedLines = lineReader.estimateTotalLines(); // Use estimation
                 if (estimatedLines == 0) { // Handle case where estimation failed or file is empty/small
                     globalLogger.warning("Could not estimate lines accurately, defaulting to 1000 for progress bar.");
                     estimatedLines = 1000;
                 }

                if (globalConfig.verbose) {
                    globalLogger.info("Estimated " + std::to_string(estimatedLines) + " data lines in CSV file for progress bar.");
                }
            } catch (const std::exception& e) {
                globalLogger.warning("Error estimating lines, defaulting to 1000: " + std::string(e.what()));
                estimatedLines = 1000; // Default value if estimation fails
            }
            
            globalLogger.info("Estimated total lines (approx): " + std::to_string(estimatedLines));
            ProgressBar progressBar(estimatedLines, "Processing", 50);
            progressBar.start();
            auto startTime = std::chrono::steady_clock::now();
            std::atomic<size_t> totalLinesProcessed{0};
            std::atomic<size_t> totalMoleculesProcessed{0};
            std::atomic<size_t> totalValidMoleculesWritten{0};

            // Main processing loop - batch based approach
            size_t effectiveBatchSize = globalConfig.batchSize;
            globalLogger.debug("Using batch size: " + std::to_string(effectiveBatchSize));
            size_t processedCount = 0;

            // Check if stream mode is enabled (default) or disabled
            bool streamMode = !result.count("no-stream");
            if (streamMode) {
                #ifdef WITH_TBB
                globalLogger.info("Using row-by-row streaming mode with parallel descriptor calculation on " + 
                                 std::to_string(globalConfig.numThreads) + " jobs");
                #else
                globalLogger.info("Using row-by-row streaming mode (single-threaded)");
                #endif
                
                // Initialize MoleculeStream
                MoleculeStream moleculeStream(inputPath);
                
                // Skip header if needed
                if (hasHeader) {
                    moleculeStream.skipHeader();
                }
                
                // Use a molecule queue to buffer between processing stages
                #ifdef WITH_TBB
                const size_t bufferSize = std::max(size_t(1), size_t(globalConfig.numThreads * 2));
                tbb::concurrent_bounded_queue<std::pair<Molecule, std::vector<std::string>>> moleculeQueue;
                moleculeQueue.set_capacity(bufferSize);
                
                std::atomic<bool> producerDone{false};
                std::atomic<size_t> consumedCount{0};
                
                // Start consumer thread group
                tbb::task_group descriptorGroup;
                
                // Function to process molecules from the queue
                auto processQueuedMolecules = [&]() {
                    std::pair<Molecule, std::vector<std::string>> item;
                    
                    while (!producerDone || !moleculeQueue.empty()) {
                        if (moleculeQueue.try_pop(item)) {
                            Molecule& molecule = item.first;
                            std::vector<std::string>& rowData = item.second;
                            
                            // Calculate descriptors for this molecule
                            std::unordered_map<std::string, std::vector<std::variant<double, int, std::string>>> results;
                            
                            for (const auto& descriptorName : descriptorNames) {
                                std::variant<double, int, std::string> result = factory.calculate(descriptorName, molecule);
                                results[descriptorName].push_back(result);
                            }
                            
                            // Create a single-element batch for writing (needs mutex for thread safety)
                            static std::mutex writeMutex;
                            {
                                std::lock_guard<std::mutex> lock(writeMutex);
                                std::vector<Molecule> moleculeBatch = {molecule};
                                std::vector<std::vector<std::string>> originalDataBatch = {rowData};
                                
                                // Write to output
                                resultWriter.writeBatch(moleculeBatch, originalDataBatch, results, descriptorNames);
                                
                                // Update progress (thread-safe)
                                progressBar.update(1);
                                consumedCount.fetch_add(1, std::memory_order_relaxed);
                            }
                        } else {
                            // No molecule available yet, yield to other threads
                            std::this_thread::yield();
                        }
                    }
                };
                
                // Start consumer threads
                for (int i = 0; i < globalConfig.numThreads; i++) {
                    descriptorGroup.run(processQueuedMolecules);
                }
                
                // Producer: read molecules and push to queue
                Molecule molecule;
                std::vector<std::string> rowData;
                
                while (moleculeStream.nextWithOriginalData(molecule, rowData, delimiter, csvInputHandler.getSmilesIndex(), escapeChar)) {
                    // Make copies for thread safety
                    std::pair<Molecule, std::vector<std::string>> item = {molecule, rowData};
                    moleculeQueue.push(item);
                    processedCount++;
                    
                    // Clear for next row
                    rowData.clear();
                }
                
                // Signal that production is complete
                producerDone = true;
                
                // Wait for all processing to complete
                descriptorGroup.wait();
                
                // Ensure final flush
                resultWriter.flush();
                #else
                // Single-threaded fallback for streaming mode
                Molecule molecule;
                std::vector<std::string> rowData;
                
                // Process one molecule at a time
                while (moleculeStream.nextWithOriginalData(molecule, rowData, delimiter, csvInputHandler.getSmilesIndex(), escapeChar)) {
                    // Calculate descriptors for this molecule
                    std::unordered_map<std::string, std::vector<std::variant<double, int, std::string>>> results;
                    
                    for (const auto& descriptorName : descriptorNames) {
                        std::variant<double, int, std::string> result = factory.calculate(descriptorName, molecule);
                        results[descriptorName].push_back(result);
                    }
                    
                    // Create a single-element batch for writing
                    std::vector<Molecule> moleculeBatch = {molecule};
                    std::vector<std::vector<std::string>> originalDataBatch = {rowData};
                    
                    // Write to output
                    resultWriter.writeBatch(moleculeBatch, originalDataBatch, results, descriptorNames);
                    
                    // Update progress
                    progressBar.update(1);
                    processedCount++;
                    
                    // Clear for next row
                    rowData.clear();
                }
                
                // Ensure final flush
                resultWriter.flush();
                #endif
            } else {
                // Original batch processing mode
                #ifdef WITH_TBB
                std::random_device rd;
                std::mt19937 gen(rd()); // Mersenne Twister for random index

                try {
                    tbb::parallel_pipeline(
                        globalConfig.numThreads,
                        // Stage 1: Read batch (serial)
                        tbb::make_filter<void, BatchLines>(
                            tbb::filter_mode::serial_in_order,
                            [&](tbb::flow_control& fc) -> BatchLines {
                                BatchLines data;
                                // Read raw lines using the line reader
                                if (!lineReader.readBatch(data.lines, effectiveBatchSize) || data.lines.empty()) {
                                    fc.stop();
                                    return data;
                                }
                                totalLinesProcessed += data.lines.size();
                                return data;
                            }
                        ) &
                        // Stage 2: Parse Lines and Extract SMILES (parallel)
                        tbb::make_filter<BatchLines, SmilesData>(
                            tbb::filter_mode::parallel,
                            [&](const BatchLines& input) -> SmilesData {
                                SmilesData data;
                                data.smiles.reserve(input.lines.size());
                                data.originalData.reserve(input.lines.size());
                                const std::string& delim = csvInputHandler.getDelimiter(); // Get delimiter once
                                int smilesIdx = csvInputHandler.getSmilesIndex();

                                for (const auto& line : input.lines) {
                                    std::vector<std::string> parsedLine = CsvIO::parseCsvLine(line, delim);
                                    if (smilesIdx >= 0 && static_cast<size_t>(smilesIdx) < parsedLine.size()) {
                                        data.smiles.push_back(parsedLine[smilesIdx]);
                                        data.originalData.push_back(std::move(parsedLine)); // Store full parsed line
                                    } else {
                                        // Handle line with missing smiles column - add placeholder?
                                        globalLogger.debug("Skipping line due to missing SMILES at index " + std::to_string(smilesIdx));
                                        data.smiles.push_back(""); // Add empty smiles
                                        data.originalData.push_back(std::move(parsedLine)); // Still store original data
                                    }
                                }
                                return data;
                            }
                        ) &
                        // Stage 3: Create MoleculeBatch & Update Display (parallel)
                        tbb::make_filter<SmilesData, std::pair<std::shared_ptr<MoleculeBatch>, std::vector<std::vector<std::string>>>>(
                           tbb::filter_mode::parallel,
                           [&](SmilesData input) -> std::pair<std::shared_ptr<MoleculeBatch>, std::vector<std::vector<std::string>>> {
                               // Create batch directly from SMILES vector
                               auto batch = std::make_shared<MoleculeBatch>(input.smiles.size());
                               batch->addSmilesBatch(input.smiles); // This handles parsing internally

                               // Pass the originalData vector along, using move semantics
                               return {batch, std::move(input.originalData)};
                           }
                        ) &
                        // Stage 4: Calculate descriptors (parallel)
                        tbb::make_filter<std::pair<std::shared_ptr<MoleculeBatch>, std::vector<std::vector<std::string>>>, BatchWithResults>(
                            tbb::filter_mode::parallel,
                             [&](std::pair<std::shared_ptr<MoleculeBatch>, std::vector<std::vector<std::string>>> input) -> BatchWithResults {
                                 auto results = factory.calculateBatch(descriptorNames, *(input.first));
                                 // Use move semantics when creating BatchWithResults
                                 return {std::move(input.first), std::move(input.second), std::move(results)};
                             }
                        ) &
                         // Stage 5: Write results (serial)
                         tbb::make_filter<BatchWithResults, void>(
                             tbb::filter_mode::serial_in_order,
                             [&](const BatchWithResults& data) {
                                 // Count valid molecules written in this batch
                                 size_t validWrittenInBatch = 0;
                                 const auto& moleculesInBatch = data.batch->getMolecules();
                                 for(const auto& m : moleculesInBatch) {
                                     if (m.isValid()) {
                                         validWrittenInBatch++;
                                     }
                                 }

                                 resultWriter.writeBatch(
                                     moleculesInBatch,
                                     data.originalData, // Pass the original data rows
                                     data.results,
                                     descriptorNames
                                 );
                                 resultWriter.flush(); // Make sure to flush after writing
                                 progressBar.update(moleculesInBatch.size());
                                 processedCount += moleculesInBatch.size();
                             }
                         )
                     );
                } catch (const std::exception& e) {
                    progressBar.finish(); // Ensure progress bar is stopped on error
                    globalLogger.fatal("Processing pipeline error: " + std::string(e.what()));
                    return 1; // Exit after fatal error
                }
                #else
                 // --- Single-threaded fallback ---
                  globalLogger.warning("TBB not enabled. Running single-threaded.");
                  std::vector<std::string> lines;
                  std::vector<std::string> smilesBatch;
                  std::vector<std::vector<std::string>> originalDataBatch;

                  while (lineReader.readBatchWithOriginalData(smilesBatch, originalDataBatch, effectiveBatchSize, csvInputHandler.getSmilesIndex(), csvInputHandler.getPropertyIndicesToKeep(), csvInputHandler.getDelimiter())) {
                      totalLinesProcessed += originalDataBatch.size(); // Count raw lines read

                      auto moleculeBatch = std::make_shared<MoleculeBatch>(smilesBatch.size());
                      moleculeBatch->addSmilesBatch(smilesBatch); // Parse SMILES

                      auto results = factory.calculateBatch(descriptorNames, *moleculeBatch);

                      size_t validWrittenInBatch = 0;
                      const auto& moleculesInBatch = moleculeBatch->getMolecules();
                       for(const auto& m : moleculesInBatch) {
                           if (m.isValid()) {
                               validWrittenInBatch++;
                           }
                       }

                      resultWriter.writeBatch(moleculesInBatch, originalDataBatch, results, descriptorNames);
                      resultWriter.flush(); // Make sure to flush after writing
                      progressBar.update(moleculesInBatch.size());
                      processedCount += moleculesInBatch.size();

                      // Clear for next batch
                      smilesBatch.clear();
                      originalDataBatch.clear();
                  }
                 #endif // WITH_TBB
            }
            
            progressBar.finish();
            auto endTime = std::chrono::steady_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);
            
            // Always display final summary regardless of verbose mode
            std::cout << "\033[1;32m✓\033[0m Processed \033[1m" << processedCount << "\033[0m valid molecules" << std::endl;
            std::cout << "\033[1;32m✓\033[0m Results written to \033[1m" << outputPath << "\033[0m" << std::endl;
            std::cout << "\033[1;32m✓\033[0m Total processing time: \033[1m" << (duration.count() / 1000.0) << "\033[0m seconds" << std::endl;

            // Keep the detailed logging for verbose mode
            if (globalConfig.verbose) {
                globalLogger.info("Processing complete. Processed " + std::to_string(processedCount) + " valid molecules.");
                globalLogger.info("Results written successfully to " + outputPath);
                globalLogger.info("Total processing time: " + std::to_string(duration.count() / 1000.0) + " seconds");
            }
        }
    } catch (const cxxopts::exceptions::parsing& e) {
        std::cerr << "\033[1;31mError parsing options:\033[0m " << e.what() << std::endl;
        printHelp(options);
        return 1;
    } catch (const DescriptorException& e) {
        std::cerr << "\033[1;31mDescriptor error:\033[0m " << e.what() << std::endl;
        return 1;
    } catch (const std::exception& e) {
        std::cerr << "\033[1;31mError:\033[0m " << e.what() << std::endl;
        return 1;
    } catch (...) {
        std::cerr << "\033[1;31mUnknown error occurred.\033[0m" << std::endl;
        return 1;
    }

    return 0;
}