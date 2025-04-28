#pragma once

#include "utils.hpp" // Include for Molecule, Logger, etc.
#include <string>
#include <vector>
#include <unordered_map>
#include <variant>
#include <fstream>
#include <atomic>
#include <mutex>
#include <unordered_set>

namespace desfact {

class CsvIO {
private:
    std::string inputPath;
    std::string outputPath;
    std::string delimiter;
    int smilesColumnIndex;
    bool hasHeaders;
    std::string escapeChar;
    std::string outFilePath; // Store separate output path
    std::string smilesColumn;
    bool hasHeader;
    std::vector<std::string> propertyColumnsToKeep; // Renamed for clarity
    std::vector<std::string> headerColumns;
    int smilesIndex = -1;
    std::unordered_set<size_t> propertyIndicesToKeep; // Renamed for clarity

    bool parseHeader(std::ifstream& file);
    void determinePropertyIndices(); // Helper to populate propertyIndicesToKeep

public:
    CsvIO(const std::string& inputPath, const std::string& outputPath,
          const std::string& delimiter, int smilesColumnIndex, bool hasHeaders,
          const std::string& escapeChar = "");

    void setPropertyColumnsToKeep(const std::vector<std::string>& columns);
    const std::vector<std::string>& getHeaderColumns() const;
    const std::vector<std::string>& getPropertyColumnsToKeep() const;
    int getSmilesIndex() const;
    std::string getDelimiter() const;
    bool getHasHeader() const;
    const std::unordered_set<size_t>& getPropertyIndicesToKeep() const;
    std::string getEscapeChar() const;

    void setSmilesColumnName(const std::string& name) { smilesColumn = name; }

    class LineReader {
    private:
        std::string filePath;
        std::mutex readMutex;
        std::ifstream fileStream;
        std::vector<char> buffer;
        uintmax_t fileSize;
        std::atomic<size_t> bytesRead{0};
        size_t estimatedLines;
        bool headerSkipped;
        bool readerHasHeader;

    public:
        LineReader(const std::string& filePath, bool hasHeaderFlag);
        ~LineReader();

        bool readBatch(std::vector<std::string>& lines, size_t batchSize);
        bool readBatchWithOriginalData(
            std::vector<std::string>& smilesLines,
            std::vector<std::vector<std::string>>& originalDataLines,
            size_t batchSize,
            int smilesColIndex,
            const std::unordered_set<size_t>& originalIndicesToKeep,
            const std::string& lineDelimiter,
            const std::string& csvIOEscapeChar = "");

        size_t estimateTotalLines();
        double getProgress() const;
    };

    LineReader createLineReader();

    class ResultWriter {
    private:
        std::string outputPath;
        std::ofstream fileStream;
        std::string delimiter;
        std::vector<std::string> outputHeaders;
        std::mutex writeMutex;
        bool headerWritten;
        const std::unordered_set<size_t>& propertyIndicesToKeepRef;
        int smilesColumnIndex;
        std::string smilesColumnName;
        std::vector<size_t> orderedOriginalIndices;
        size_t smilesOutputIndex = 0;
        std::vector<size_t> descriptorOutputIndices;
        std::string escapeChar;

        void buildOutputHeaders(const std::vector<std::string>& inputHeaders,
                                const std::vector<std::string>& descriptorNames);

    public:
        ResultWriter(const std::string& outFilePath, const std::string& delimiter,
                   int smilesColumnIndex, const std::vector<std::string>& headerColumns,
                   const std::vector<std::string>& descriptorNames,
                   const std::unordered_set<size_t>& propertyIndices,
                   const std::string& escapeChar = "");
        ~ResultWriter();

        bool writeBatch(
            const std::vector<Molecule>& molecules, // Molecules contain canonical SMILES and maybe props
            const std::vector<std::vector<std::string>>& originalData, // Original CSV data rows
            const std::unordered_map<std::string, std::vector<std::variant<double, int, std::string>>>& results, // Calculated descriptors
            const std::vector<std::string>& descriptorNames // Order of descriptors in results
        );

        void flush();
    };

    ResultWriter createResultWriter(const std::vector<std::string>& descriptorNames, 
                                   const std::string& escapeChar);

    static std::vector<std::string> parseCsvLine(const std::string& line, 
                                                const std::string& delimiter,
                                                const std::string& escapeChar = "");
    static std::string trimQuotes(const std::string& str);

    static const std::unordered_set<size_t>& getDefaultPropertyIndices() {
        static std::unordered_set<size_t> emptySet;
        return emptySet;
    }

    bool reparseHeader();
};

} // namespace desfact 