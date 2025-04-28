#include "io.hpp"
#include "utils.hpp" // Access to globalLogger
#include <cmath>
#include <filesystem>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <algorithm> // For std::find, std::sort

namespace desfact {

std::string CsvIO::trimQuotes(const std::string& str) {
    size_t first = str.find_first_not_of(" \t\n\r\f\v\"");
    if (std::string::npos == first) return ""; // Return empty if all quotes/whitespace
    size_t last = str.find_last_not_of(" \t\n\r\f\v\"");
    return str.substr(first, (last - first + 1));
}

std::vector<std::string> CsvIO::parseCsvLine(const std::string& line, const std::string& delimiter, const std::string& escapeCharStr) {
    std::vector<std::string> cells;
    if (line.empty()) return cells;

    std::string currentCell;
    char delimChar = delimiter[0];
    char quoteChar = '"';
    char escCharVal = escapeCharStr.empty() ? '\0' : escapeCharStr[0];
    bool inQuotes = false;

    for (size_t i = 0; i < line.length(); ++i) {
        char c = line[i];
        if (inQuotes) {
            if (escCharVal != '\0' && c == escCharVal) {
                // Handle custom escape character if set
                if (i + 1 < line.length()) {
                    char nextC = line[i+1];
                    // Append the escaped character literally.
                    if (nextC == 'n') currentCell += '\n'; // Handle literal \n
                    else if (nextC == 'r') currentCell += '\r'; // Handle literal \r
                    else if (nextC == 't') currentCell += '\t';
                    else if (nextC == '\\') currentCell += '\\';
                    else if (nextC == '"') currentCell += '"';
                    else currentCell += nextC; // Handle escaped quotes, escape char itself, etc.
                    i++; // Consume next char
                } else {
                    // Escape char at end, treat as literal
                    currentCell += c;
                }
            } else if (c == quoteChar) {
                // Handle quote character while in quotes.
                if (escCharVal == '\0' && i + 1 < line.length() && line[i+1] == quoteChar) {
                    // Standard CSV: "" means a literal " (only if no custom escape is set)
                    currentCell += quoteChar;
                    i++; // Consume the second quote
                } else {
                    // Found a single quote while inQuotes. This is the closing quote.
                    inQuotes = false; // Standard CSV closing quote
                    // Do NOT append the closing quote
                }
            } else {
                // Any other character inside quotes (including newlines read by getline)
                currentCell += c;
            }
        } else { // Not in quotes
            if (c == quoteChar) {
                // Found a quote while NOT in quotes. This must be the opening quote
                // of a field *if* it's the very first character of the field content.
                if (currentCell.empty()) { // Check if current field is empty so far
                    inQuotes = true;
                    // Do NOT append the opening quote
                } else {
                    // Quote appeared mid-field outside quotes - treat as literal quote (non-standard but possible)
                    currentCell += c;
                }
            } else if (c == delimChar) { // Delimiter outside quotes ends the field
                cells.push_back(currentCell);
                currentCell.clear();
            } else { // Any other character outside quotes
                currentCell += c;
            }
        }
    }
    
    cells.push_back(currentCell); // Add the last cell
    return cells;
}


CsvIO::CsvIO(const std::string& inputPath, const std::string& outputPath,
            const std::string& delimiter, int smilesColumnIndex, bool hasHeaders,
            const std::string& escapeChar)
    : inputPath(inputPath), outputPath(outputPath), delimiter(delimiter),
      smilesColumnIndex(smilesColumnIndex), hasHeaders(hasHeaders),
      escapeChar(escapeChar), smilesColumn("SMILES") { // Initialize smilesColumn with default name

    std::ifstream file(inputPath);
    if (!file.is_open()) {
        throw DescriptorException("CsvIO: Failed to open input file: " + inputPath, ErrorCode::IO_ERROR);
    }

    if (hasHeaders) {
        if (!parseHeader(file)) {
             // parseHeader logs the error, we throw to stop construction
             throw DescriptorException("CsvIO: Failed to find SMILES column '" + 
                 (smilesColumnIndex >= 0 ? std::to_string(smilesColumnIndex) : smilesColumn) + 
                 "' in header.", ErrorCode::PARSE_ERROR);
        }
        // After parsing header, determine which property indices to keep
        determinePropertyIndices();
    } else {
        // Assume SMILES is the first column if no header
        smilesIndex = 0;
        globalLogger.info("CsvIO: No header specified, assuming SMILES in column 0.");
        // Cannot determine property indices by name without a header
        // User must specify indices or handle columns numerically if needed
    }
    file.close();
}

bool CsvIO::parseHeader(std::ifstream& file) {
    std::string headerLine;
    if (std::getline(file, headerLine)) {
        // Handle potential UTF-8 BOM
        if (headerLine.size() >= 3 && headerLine.rfind("\xEF\xBB\xBF", 0) == 0) {
            headerLine = headerLine.substr(3);
        }
        // Remove trailing CR/LF
        while (!headerLine.empty() && (headerLine.back() == '\r' || headerLine.back() == '\n')) {
            headerLine.pop_back();
        }

        headerColumns = parseCsvLine(headerLine, delimiter, escapeChar);
        if (headerColumns.empty()) {
            globalLogger.error("CsvIO: Header line parsed into zero columns. Check delimiter or file format.");
            return false;
        }

        // First check if we have a numeric index
        if (smilesColumnIndex >= 0 && smilesColumnIndex < static_cast<int>(headerColumns.size())) {
            smilesIndex = smilesColumnIndex;
            smilesColumn = headerColumns[smilesIndex];
            globalLogger.info("CsvIO: Using SMILES column at index " + std::to_string(smilesIndex) + 
                             " with name '" + smilesColumn + "'");
            return true;
        }
        
        // If not a valid index, look for the column by name
        for (size_t i = 0; i < headerColumns.size(); ++i) {
            if (headerColumns[i] == smilesColumn) {
                smilesIndex = static_cast<int>(i);
                globalLogger.info("CsvIO: Found SMILES column '" + smilesColumn + "' at index " + 
                                 std::to_string(smilesIndex));
                return true; // Found smiles column
            }
        }
        
        // If loop finishes without finding smiles column
        globalLogger.error("CsvIO: SMILES column '" + smilesColumn + "' not found in header.");
        smilesIndex = -1;
        return false;
    } else {
        globalLogger.error("CsvIO: Failed to read header line from file: " + inputPath);
        smilesIndex = -1;
        return false;
    }
}

// Populate propertyIndicesToKeep based on propertyColumnsToKeep and headerColumns
void CsvIO::determinePropertyIndices() {
    propertyIndicesToKeep.clear();
    if (!hasHeaders) {
        globalLogger.warning("CsvIO: Cannot determine property indices by name without a header.");
        return;
    }

    if (propertyColumnsToKeep.empty()) {
         // If no specific columns requested, keep *all* non-SMILES columns by default
         for (size_t i = 0; i < headerColumns.size(); ++i) {
             if (static_cast<int>(i) != smilesIndex) {
                 propertyIndicesToKeep.insert(i);
             }
         }
         globalLogger.info("CsvIO: No specific property columns requested. Keeping all original columns except SMILES.");

    } else {
        // Keep only the explicitly requested columns
        for (const auto& colName : propertyColumnsToKeep) {
            bool found = false;
            for (size_t i = 0; i < headerColumns.size(); ++i) {
                if (headerColumns[i] == colName) {
                    if (static_cast<int>(i) == smilesIndex) {
                         globalLogger.warning("CsvIO: Requested property column '" + colName + "' is the same as the SMILES column. It will be included via the main SMILES output, not as a separate property.");
                    } else {
                        propertyIndicesToKeep.insert(i);
                        found = true;
                    }
                    break; // Found the column name, move to next requested column
                }
            }
            if (!found) {
                globalLogger.warning("CsvIO: Requested property column '" + colName + "' not found in header.");
            }
        }
         globalLogger.info("CsvIO: Keeping " + std::to_string(propertyIndicesToKeep.size()) + " specified original property columns.");
    }
}


void CsvIO::setPropertyColumnsToKeep(const std::vector<std::string>& columns) {
    propertyColumnsToKeep = columns;
    // Re-determine indices based on the new list and existing header
    if (hasHeaders && !headerColumns.empty()) {
        determinePropertyIndices();
    } else if (!hasHeaders && !columns.empty()) {
         globalLogger.warning("CsvIO::setPropertyColumnsToKeep: Cannot use column names without a header. Indices must be used or inferred.");
         // Clear indices if names are set without header context
         propertyIndicesToKeep.clear();
    }
}

const std::vector<std::string>& CsvIO::getHeaderColumns() const { return headerColumns; }
const std::vector<std::string>& CsvIO::getPropertyColumnsToKeep() const { return propertyColumnsToKeep; }
int CsvIO::getSmilesIndex() const { return smilesIndex; }
std::string CsvIO::getDelimiter() const { return delimiter; }
bool CsvIO::getHasHeader() const { return hasHeaders; }
const std::unordered_set<size_t>& CsvIO::getPropertyIndicesToKeep() const { return propertyIndicesToKeep; }

std::string CsvIO::getEscapeChar() const {
    return escapeChar;
}


// --- LineReader Implementation ---
CsvIO::LineReader::LineReader(const std::string& filePath, bool hasHeaderFlag)
    : readerHasHeader(hasHeaderFlag),
      headerSkipped(false),
      buffer(65536), // 64KB buffer
      estimatedLines(0) {

    fileStream.open(filePath, std::ios::binary); // Open in binary for accurate size/position
    if (!fileStream.is_open()) {
        throw DescriptorException("LineReader: Failed to open file: " + filePath, ErrorCode::IO_ERROR);
    }
    fileStream.rdbuf()->pubsetbuf(buffer.data(), buffer.size());

    try {
        fileSize = std::filesystem::file_size(filePath);
    } catch (const std::filesystem::filesystem_error& e) {
        fileSize = 0;
        globalLogger.warning("LineReader: Could not determine file size for " + filePath + ": " + std::string(e.what()));
    }

    estimatedLines = estimateTotalLines(); // Estimate lines on construction
}

CsvIO::LineReader::~LineReader() {
    if (fileStream.is_open()) {
        fileStream.close();
    }
}

bool CsvIO::LineReader::readBatch(std::vector<std::string>& lines, size_t batchSize) {
    std::lock_guard<std::mutex> lock(readMutex);
    lines.clear();
    lines.reserve(batchSize);

    if (!fileStream.is_open() || fileStream.eof()) {
        return false;
    }

    // Skip header only once
    if (readerHasHeader && !headerSkipped) {
        std::string headerLine;
        std::streampos beforeHeader = fileStream.tellg();
        if (std::getline(fileStream, headerLine)) {
             std::streampos afterHeader = fileStream.tellg();
             bytesRead += (afterHeader - beforeHeader); // More accurate byte count
             headerSkipped = true;
        } else {
            globalLogger.warning("LineReader: Failed reading or skipping header.");
             fileStream.clear(); // Clear error flags if header read failed at EOF
             fileStream.seekg(beforeHeader); // Reset position
            return false; // Indicate no lines read
        }
    }

    std::string line;
    line.reserve(256); // Pre-allocate line buffer
    size_t linesRead = 0;
    std::streampos batchStartPos = fileStream.tellg();

    while (linesRead < batchSize && std::getline(fileStream, line)) {
        // Basic check for empty lines or lines just containing whitespace/CR/LF
        if (line.find_first_not_of(" \t\r\n") != std::string::npos) {
            lines.push_back(line);
            linesRead++;
        } else {
             globalLogger.debug("LineReader: Skipping empty or whitespace-only line.");
        }
    }
    std::streampos batchEndPos = fileStream.tellg();
    if(batchEndPos > batchStartPos) {
       bytesRead += (batchEndPos - batchStartPos);
    } else if (fileStream.eof() && batchStartPos != std::streampos(-1)) {
       // If we hit EOF exactly, count bytes to the end
       bytesRead = fileSize; // Assume we read the rest
    }


    return !lines.empty();
}

// Reads full lines, extracts SMILES, and keeps original data for specified columns
bool CsvIO::LineReader::readBatchWithOriginalData(
    std::vector<std::string>& smilesLines,
    std::vector<std::vector<std::string>>& originalDataLines,
    size_t batchSize,
    int smilesColIndex,
    const std::unordered_set<size_t>& _originalIndicesToKeep,
    const std::string& lineDelimiter,
    const std::string& csvIOEscapeChar)
{
    std::lock_guard<std::mutex> lock(readMutex);
    smilesLines.clear();
    originalDataLines.clear();
    smilesLines.reserve(batchSize);
    originalDataLines.reserve(batchSize);

    if (!fileStream.is_open() || fileStream.eof() || smilesColIndex < 0) {
        if (smilesColIndex < 0) globalLogger.error("LineReader: Invalid SMILES column index provided.");
        return false;
    }

    if (readerHasHeader && !headerSkipped) {
        std::string headerLine;
        std::streampos beforeHeader = fileStream.tellg();
         if (std::getline(fileStream, headerLine)) {
              std::streampos afterHeader = fileStream.tellg();
              bytesRead += (afterHeader - beforeHeader);
              headerSkipped = true;
         } else {
             globalLogger.warning("LineReader: Failed reading or skipping header.");
              fileStream.clear();
              fileStream.seekg(beforeHeader);
             return false;
         }
    }

    std::string line;
    line.reserve(512);
    size_t linesProcessed = 0;
    std::streampos batchStartPos = fileStream.tellg();

    while (linesProcessed < batchSize && std::getline(fileStream, line)) {
         if (line.find_first_not_of(" \t\r\n") == std::string::npos) {
             globalLogger.debug("LineReader: Skipping empty or whitespace-only line during data extraction.");
             continue;
         }

        std::vector<std::string> parsedLine = CsvIO::parseCsvLine(line, lineDelimiter, csvIOEscapeChar);

        if (static_cast<size_t>(smilesColIndex) < parsedLine.size()) {
            smilesLines.push_back(parsedLine[smilesColIndex]);
            originalDataLines.push_back(std::move(parsedLine));
            linesProcessed++;
        } else {
            globalLogger.warning("LineReader: SMILES column index " + std::to_string(smilesColIndex) +
                                 " out of bounds for line: " + line.substr(0, 50) + "...");
        }
    }
    std::streampos batchEndPos = fileStream.tellg();
     if(batchEndPos > batchStartPos) {
        bytesRead += (batchEndPos - batchStartPos);
     } else if (fileStream.eof() && batchStartPos != std::streampos(-1)) {
        bytesRead = fileSize;
     }


    return !smilesLines.empty();
}


size_t CsvIO::LineReader::estimateTotalLines() {
    if (fileSize == 0) return 0;

    std::lock_guard<std::mutex> lock(readMutex);

    // Save original file position
    std::streampos originalPos = fileStream.tellg();
    fileStream.clear();
    fileStream.seekg(0, std::ios::beg);

    // Handle header
    std::string line;
    if (readerHasHeader) {
        if (!std::getline(fileStream, line)) {
            fileStream.clear();
            fileStream.seekg(originalPos);
            return 0;
        }
    }

    // Count all lines
    size_t lineCount = 0;
    while (std::getline(fileStream, line)) {
        if (line.find_first_not_of(" \t\r\n") != std::string::npos) {
            lineCount++;
        }
    }

    // Restore original position
    fileStream.clear();
    fileStream.seekg(originalPos);

    globalLogger.debug("Exact line count: " + std::to_string(lineCount));
    return lineCount;
}


double CsvIO::LineReader::getProgress() const {
    if (fileSize == 0) return 1.0; // Or 0.0? If no file, maybe progress is complete?
    size_t currentBytes = bytesRead.load(std::memory_order_relaxed);
    return std::min(static_cast<double>(currentBytes) / fileSize, 1.0);
}

CsvIO::LineReader CsvIO::createLineReader() {
     if (smilesIndex < 0 && hasHeaders) {
         throw DescriptorException("Cannot create LineReader: SMILES column index is invalid. Check CsvIO setup.", ErrorCode::PARSE_ERROR);
     }
    // Pass file path and header flag. Delimiter will be passed to readBatchWithOriginalData
    return LineReader(inputPath, hasHeaders);
}


// --- ResultWriter Implementation ---

// Helper to build the header and determine column order
void CsvIO::ResultWriter::buildOutputHeaders(
    const std::vector<std::string>& inputHeaders,
    const std::vector<std::string>& descriptorNames)
{
    outputHeaders.clear();
    orderedOriginalIndices.clear();
    descriptorOutputIndices.clear();

    // Just use the input headers directly without changes
    outputHeaders = inputHeaders;
    
    // Keep track of original indices for reference
    for (size_t i = 0; i < inputHeaders.size(); i++) {
        orderedOriginalIndices.push_back(i);
        if (static_cast<int>(i) == smilesColumnIndex) {
            smilesOutputIndex = i;
            smilesColumnName = inputHeaders[i];
        }
    }

    // Track where descriptor columns will begin
    size_t descriptorStartIndex = outputHeaders.size();
    for(size_t i = 0; i < descriptorNames.size(); ++i) {
        descriptorOutputIndices.push_back(descriptorStartIndex + i);
    }

    globalLogger.debug("ResultWriter: Output header order determined (" + 
                     std::to_string(outputHeaders.size()) + " columns plus " +
                     std::to_string(descriptorNames.size()) + " descriptors).");
}


CsvIO::ResultWriter::ResultWriter(const std::string& outFilePath, const std::string& delimiter,
                                int smilesColumnIndex, const std::vector<std::string>& headerColumns,
                                const std::vector<std::string>& descriptorNames,
                                const std::unordered_set<size_t>& propertyIndices,
                                const std::string& escapeChar)
    : outputPath(outFilePath), fileStream(), delimiter(delimiter), outputHeaders(),
      writeMutex(), headerWritten(false), 
      propertyIndicesToKeepRef(propertyIndices),
      smilesColumnIndex(smilesColumnIndex), 
      orderedOriginalIndices(), smilesOutputIndex(0),
      descriptorOutputIndices(), escapeChar(escapeChar)
{
    // Set the SMILES column name from the header columns if it exists
    if (smilesColumnIndex >= 0 && smilesColumnIndex < static_cast<int>(headerColumns.size())) {
        smilesColumnName = headerColumns[smilesColumnIndex];
    } else {
        // Default name if index is invalid
        smilesColumnName = "SMILES";
    }
    
    buildOutputHeaders(headerColumns, descriptorNames);
    
    fileStream.open(outFilePath, std::ios::out | std::ios::trunc | std::ios::binary);
    if (!fileStream.is_open()) {
        throw DescriptorException("ResultWriter: Failed to open output file: " + outFilePath, ErrorCode::IO_ERROR);
    }

    static char buffer[65536];
    fileStream.rdbuf()->pubsetbuf(buffer, sizeof(buffer));

    globalLogger.info("ResultWriter: Initialized for output file: " + outFilePath);
}


CsvIO::ResultWriter::~ResultWriter() {
    if (fileStream.is_open()) {
        // Ensure buffer is flushed before closing
        fileStream.flush();
        fileStream.close();
         globalLogger.debug("ResultWriter: Closed output file stream.");
    }
}


bool CsvIO::ResultWriter::writeBatch(
    const std::vector<Molecule>& molecules,
    const std::vector<std::vector<std::string>>& originalData,
    const std::unordered_map<std::string, std::vector<std::variant<double, int, std::string>>>& results,
    const std::vector<std::string>& descriptorNames) {

    std::lock_guard<std::mutex> lock(writeMutex);

    if (!fileStream.is_open()) {
        globalLogger.error("ResultWriter::writeBatch: Output file stream is not open.");
        return false;
    }

    std::ostringstream batchBuffer;

    // Write header if not already written
    if (!headerWritten) {
        // Write the original headers first, exactly as they appeared
        for (size_t i = 0; i < outputHeaders.size(); ++i) {
            if (i > 0) batchBuffer << delimiter;
            
            std::string header = outputHeaders[i];
            bool needsQuotes = header.find(delimiter) != std::string::npos || 
                             header.find('"') != std::string::npos || 
                             header.find_first_of(" \t\n\r") != std::string::npos;
            
            if (needsQuotes) batchBuffer << '"';
            for (char c : header) { // Handle escaped quotes within header
                if (c == '"') batchBuffer << "\"\"";
                else batchBuffer << c;
            }
            if (needsQuotes) batchBuffer << '"';
        }
        
        // Now append the descriptor names headers
        for (const auto& descName : descriptorNames) {
            batchBuffer << delimiter;
            
            bool needsQuotes = descName.find(delimiter) != std::string::npos || 
                             descName.find('"') != std::string::npos || 
                             descName.find_first_of(" \t\n\r") != std::string::npos;
            
            if (needsQuotes) batchBuffer << '"';
            for (char c : descName) { // Handle escaped quotes within header
                if (c == '"') batchBuffer << "\"\"";
                else batchBuffer << c;
            }
            if (needsQuotes) batchBuffer << '"';
        }
        
        batchBuffer << "\n"; // Use '\n' consistently
        headerWritten = true;
    }

    if (molecules.size() != originalData.size()) {
        globalLogger.fatal("ResultWriter::writeBatch: Mismatch between molecule count (" + 
                         std::to_string(molecules.size()) +
                         ") and original data row count (" + std::to_string(originalData.size()) + 
                         "). Data will be misaligned or dropped!");
    }

    // Write data rows for ALL molecules (valid or invalid)
    for (size_t i = 0; i < molecules.size(); ++i) {
        const auto& mol = molecules[i];

        // Ensure we have corresponding original data
        if (i >= originalData.size()) {
            globalLogger.warning("ResultWriter::writeBatch: Missing original data for molecule index " + std::to_string(i) + ". Skipping row.");
            continue;
        }
        const auto& rowData = originalData[i];

        std::ostringstream rowStream;
        
        // Write all original columns, replacing SMILES with canonical version if valid
        for (size_t j = 0; j < rowData.size(); ++j) {
            if (j > 0) rowStream << delimiter;
            
            if (j == static_cast<size_t>(smilesColumnIndex) && mol.isValid()) {
                // Write canonical SMILES for the SMILES column if valid
                rowStream << mol.getSmiles();
            } else {
                // For all other columns or invalid SMILES, write original data
                std::string cellData = rowData[j];
                // Remove trailing CR if present
                if (!cellData.empty() && cellData.back() == '\r') {
                    cellData.pop_back();
                }
                
                bool needsQuotes = cellData.find(delimiter) != std::string::npos || 
                                  cellData.find('"') != std::string::npos || 
                                  cellData.find_first_of(" \t\n\r") != std::string::npos;
                
                if (needsQuotes) rowStream << '"';
                for(char c : cellData) {
                    if (c == '\n') {
                        if (escapeChar.empty()) {
                            rowStream << '\n'; // Write actual newline in quoted field
                        } else {
                            rowStream << escapeChar << 'n'; // Write as escape sequence
                        }
                    } else if (c == '"') {
                        if (escapeChar.empty()) {
                            rowStream << "\"\""; // Standard escaping
                        } else {
                            rowStream << escapeChar << c;
                        }
                    } else {
                        rowStream << c;
                    }
                }
                if (needsQuotes) rowStream << '"';
            }
        }

        // Write descriptor columns
        for (const auto& descName : descriptorNames) {
            rowStream << delimiter;
            auto it = results.find(descName);
            if (it != results.end() && i < it->second.size()) {
                const auto& value = it->second[i];
                std::visit([&rowStream, &delim = delimiter, this](auto&& arg) {
                    using T = std::decay_t<decltype(arg)>;
                    if constexpr (std::is_same_v<T, double>) {
                        if (std::isnan(arg)) rowStream << "NA";
                        else rowStream << std::fixed << std::setprecision(6) << arg;
                    } else if constexpr (std::is_same_v<T, int>) {
                        if (arg == -1) rowStream << "NA";
                        else rowStream << arg;
                    } else if constexpr (std::is_same_v<T, std::string>) {
                        bool needsQuotes = arg.find(delim) != std::string::npos || 
                                         arg.find('"') != std::string::npos || 
                                         arg.find_first_of(" \t\n\r") != std::string::npos;
                        
                        if (needsQuotes) rowStream << '"';
                        for(char c : arg) {
                            if (c == '\n') {
                                if (escapeChar.empty()) {
                                    rowStream << '\n'; // Write actual newline in quoted field
                                } else {
                                    rowStream << escapeChar << 'n'; // Write as escape sequence
                                }
                            } else if (c == '"') {
                                if (escapeChar.empty()) {
                                    rowStream << "\"\""; // Standard escaping
                                } else {
                                    rowStream << escapeChar << c;
                                }
                            } else {
                                rowStream << c;
                            }
                        }
                        if (needsQuotes) rowStream << '"';
                    }
                }, value);
            } else {
                // Handle missing descriptor result for this molecule/descriptor
                globalLogger.debug("ResultWriter: Missing result for descriptor '" + descName + "' at index " + std::to_string(i));
                rowStream << "NA"; // Or empty string, or specific error value
            }
        }

        batchBuffer << rowStream.str() << "\n"; // Append row to batch buffer
    }

    // Create a string with consistent line endings
    std::string output = batchBuffer.str();
    // Remove any existing CR characters to prevent CRLF issues
    output.erase(std::remove(output.begin(), output.end(), '\r'), output.end());
    fileStream << output;

    // Check stream state
    if (!fileStream.good()) {
        globalLogger.error("ResultWriter::writeBatch: File stream encountered an error after writing batch.");
        // Attempt to recover or just return failure
        fileStream.clear(); // Clear error flags
        return false;
    }

    return true;
}


CsvIO::ResultWriter CsvIO::createResultWriter(const std::vector<std::string>& descriptorNames,
                                            const std::string& escapeChar) {
     if (smilesIndex < 0 && hasHeaders) {
          throw DescriptorException("Cannot create ResultWriter: SMILES column index is invalid. Check CsvIO setup.", ErrorCode::PARSE_ERROR);
     }
     if (outputPath.empty()) {
          throw DescriptorException("Cannot create ResultWriter: Output file path is not set.", ErrorCode::IO_ERROR);
     }
     if (smilesColumn.empty()) { // Also check if smilesColumn name is valid
          throw DescriptorException("Cannot create ResultWriter: SMILES column name is not set.", ErrorCode::PARSE_ERROR);
     }
     if (hasHeaders && headerColumns.empty()) {
         throw DescriptorException("Cannot create ResultWriter: Header expected but not parsed correctly.", ErrorCode::PARSE_ERROR);
     }

    // Pass the propertyIndicesToKeep from this CsvIO instance
    return ResultWriter(outputPath, delimiter, smilesIndex, headerColumns, 
                      descriptorNames, propertyIndicesToKeep, escapeChar);
}

bool CsvIO::reparseHeader() {
    std::ifstream file(inputPath);
    if (!file.is_open()) {
        globalLogger.error("CsvIO: Failed to open input file for header re-parsing: " + inputPath);
        return false;
    }
    bool result = parseHeader(file);
    file.close();
    if (result) {
        determinePropertyIndices(); // Re-determine property indices
    }
    return result;
}

void CsvIO::ResultWriter::flush() {
    if (fileStream.is_open()) {
        fileStream.flush();
    }
}

} // namespace desfact 