#pragma once

#include <string>
#include <vector>
#include <unordered_map>
#include <memory>
#include <stdexcept>
#include <chrono>
#include <mutex>
#include <atomic>
#include <thread>
#include <functional>
#include <utility>
#include <ostream>
#include <iostream> // For default std::cout in Logger
#include <fstream>

// RDKit Forward Declarations (if needed, or include specific headers)
namespace RDKit {
    class ROMol;
    class RWMol;
}

namespace desfact {

enum class ErrorCode {
    SUCCESS = 0,
    PARSE_ERROR,
    IO_ERROR,
    CALCULATION_ERROR,
    MEMORY_ERROR,
    NOT_IMPLEMENTED,
    UNKNOWN_ERROR
};

class DescriptorException : public std::runtime_error {
private:
    ErrorCode code;

public:
    DescriptorException(const std::string& message, ErrorCode code = ErrorCode::UNKNOWN_ERROR);
    ErrorCode getCode() const;
};

struct Config {
    int numThreads = 1;
    bool verbose = false;
    std::string logLevel = "INFO";
    size_t batchSize = 1000;
    std::string tempDir = "/tmp";
};

extern Config globalConfig;

namespace util {
    std::string getTimeStamp();
    template<typename T>
    void safeDelete(T*& ptr);
}


enum class LogLevel {
    DEBUG,
    INFO,
    WARNING,
    ERROR,
    FATAL
};

class ProgressBar;

class Logger {
private:
    LogLevel minLevel;
    std::mutex logMutex;
    std::ostream& out;
    std::ostream& err_out;
    bool colorEnabled;
    ProgressBar* activeProgressBar;
    std::vector<std::pair<LogLevel, std::string>> bufferedMessages;

    static const char* levelToString(LogLevel level);
    const char* levelToColor(LogLevel level);
    std::string formatMessage(LogLevel level, const std::string& message);

public:
    Logger(LogLevel minLevel = LogLevel::WARNING, 
          std::ostream& out_stream = std::cout, 
          std::ostream& err_stream = std::cerr, 
          bool colorEnabled = true);

    void log(LogLevel level, const std::string& message);
    void debug(const std::string& message);
    void info(const std::string& message);
    void warning(const std::string& message);
    void error(const std::string& message);
    void fatal(const std::string& message);

    void setMinLevel(LogLevel level);
    void enableColor(bool enable);

    void setActiveProgressBar(ProgressBar* progressBar);
    void clearActiveProgressBar();
    void flushBufferedMessages();
};

extern Logger globalLogger;

class ProgressBar {
private:
    size_t total;
    std::atomic<size_t> current{0};
    std::string description;
    size_t barWidth;
    std::mutex renderMutex;
    bool active = false;
    std::chrono::time_point<std::chrono::steady_clock> startTime;
    std::chrono::time_point<std::chrono::steady_clock> lastUpdateTime;
    size_t lastRenderCount;
    std::thread infoThread;
    int updateFrequency;
    std::chrono::time_point<std::chrono::steady_clock> lastRenderTime;
    double currentRate;
    std::mutex mutex;
    std::string barStyle = "gradient";
    std::thread displayThread;
    
    // Dynamic sizing
    int getTerminalWidth() const;
    int calculateBarWidth(int termWidth) const;
    int getMinBarWidth() const { return 10; }
    int getMaxBarWidth() const { return 40; }

    void displayLoop();
    void render();
    void renderFinal();
    void clearLine();
    std::string getColorGradient(double progress) const;
    std::string getBlockBar(double progress, int width) const;
    std::string getETA() const;

public:
    ProgressBar(size_t total,
                const std::string& description = std::string("Processing"),
                int updateFrequencyMs = 100);
    ~ProgressBar();

    void start();
    void update(size_t increment = 1);
    void finish();
    void setDescription(const std::string& desc);
    void setTotal(size_t newTotal);
    double getProgress() const;
    std::chrono::milliseconds getElapsedTime() const;
    std::chrono::milliseconds getEstimatedTimeRemaining() const;
    void temporarilyPauseDisplay(const std::function<void()>& callback);
    void forceRender();
};


class Molecule {
private:
    std::shared_ptr<RDKit::ROMol> mol;
    std::string smiles;
    std::string originalSmiles;
    std::unordered_map<std::string, std::string> properties;
    bool valid;
    std::string errorMessage;

public:
    Molecule();
    explicit Molecule(const std::string& smiles);
    explicit Molecule(RDKit::ROMol* mol, bool takeOwnership = true);

    bool parse(const std::string& smiles);
    bool isValid() const;
    const std::string& getErrorMessage() const;

    std::shared_ptr<RDKit::ROMol> getMolecule() const;
    const std::string& getSmiles() const;
    const std::string& getOriginalSmiles() const;

    void sanitize();
    void canonicalize();

    bool hasProperty(const std::string& name) const;
    std::string getProperty(const std::string& name) const;
    void setProperty(const std::string& name, const std::string& value);
    const std::unordered_map<std::string, std::string>& getProperties() const;

    int getNumAtoms() const;
    int getNumBonds() const;
    std::vector<int> getAtomicNumbers() const;
    std::vector<int> getBondTypes() const;

    std::string toJSON() const;
    static Molecule fromJSON(const std::string& json);

    void setOriginalSmiles(const std::string& original) { originalSmiles = original; }
    void setErrorMessage(const std::string& message) { errorMessage = message; }
};

class MoleculeBatch {
private:
    std::vector<Molecule> molecules;
    size_t capacity;
    std::atomic<size_t> processed;
    mutable std::mutex mutex;

public:
    explicit MoleculeBatch(size_t initialSize = 0);

    bool add(const Molecule& mol);
    bool add(const std::string& smiles);
    void addSmilesBatch(const std::vector<std::string>& smiles);

    size_t size() const;
    size_t getCapacity() const;
    size_t getProcessed() const;
    void incrementProcessed(size_t count = 1);

    const std::vector<Molecule>& getMolecules() const;
    std::vector<Molecule>& getMolecules();

    void clear();
    bool isFull() const;
};

class MoleculeStream {
private:
    std::ifstream inputFile;
    std::string currentLine;
    size_t processedCount;
    bool isOpen;

public:
    explicit MoleculeStream(const std::string& filename);
    ~MoleculeStream();

    bool open(const std::string& filename);
    void close();
    bool hasNext();
    bool next(Molecule& molecule);
    bool skipHeader();
    size_t getProcessedCount() const;
    
    bool nextWithOriginalData(Molecule& molecule, std::vector<std::string>& rowData, 
                             const std::string& delimiter, int smilesIndex,
                             const std::string& escapeChar = "");
};

} // namespace desfact 