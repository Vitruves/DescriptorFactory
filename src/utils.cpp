#include "utils.hpp"

// RDKit includes for implementation
#include <GraphMol/GraphMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/Descriptors/Property.h>
#include <GraphMol/Atom.h>
#include <GraphMol/Bond.h>
#include <GraphMol/MolStandardize/MolStandardize.h>
#include <GraphMol/SmilesParse/SmartsWrite.h> // If needed, e.g., for SMARTS matching
#include <GraphMol/Fingerprints/MorganFingerprints.h> // If needed for fingerprints
#include <GraphMol/Descriptors/MolDescriptors.h> // If needed for descriptors
#include <GraphMol/Substruct/SubstructMatch.h> // If needed for substructure matching

// RapidJSON includes for Molecule JSON
#include <rapidjson/document.h>
#include <rapidjson/stringbuffer.h>
#include <rapidjson/writer.h>

// Include IO for CsvIO
#include "io.hpp"

// TBB includes (conditionally)
#ifdef WITH_TBB
#include <tbb/parallel_for.h>
#include <tbb/concurrent_vector.h>
#include <tbb/blocked_range.h>
#endif

// Includes for ProgressBar display
#include <unistd.h>
#include <sys/ioctl.h>
#include <cmath>
#include <iomanip> // Keep here for getTimeStamp
#include <sstream> // Keep here for getTimeStamp
#include <functional> // For std::function
#include <cstdio> // For popen, pclose
#include <cstdlib> // For getenv

namespace desfact {

// --- Global Variables ---
Config globalConfig;
Logger globalLogger(LogLevel::WARNING, std::cout, std::cerr, true);

// --- DescriptorException ---
DescriptorException::DescriptorException(const std::string& message, ErrorCode code)
    : std::runtime_error(message), code(code) {}

ErrorCode DescriptorException::getCode() const { return code; }

// --- Utility Functions ---
namespace util {
    std::string getTimeStamp() {
        auto now = std::chrono::system_clock::now();
        auto in_time_t = std::chrono::system_clock::to_time_t(now);
        std::stringstream ss;
        ss << std::put_time(std::localtime(&in_time_t), "%Y-%m-%d %X");
        return ss.str();
    }

    template<typename T>
    void safeDelete(T*& ptr) {
        if (ptr) {
            delete ptr;
            ptr = nullptr;
        }
    }
    // Explicit template instantiation if needed, or keep in header for templates
    // template void safeDelete<SomeType>(SomeType*& ptr);
}

// --- Logger Implementation ---
Logger::Logger(LogLevel minLevel, std::ostream& out_stream, std::ostream& err_stream, bool colorEnabled)
    : minLevel(minLevel), out(out_stream), err_out(err_stream), 
      colorEnabled(colorEnabled), activeProgressBar(nullptr) {}

const char* Logger::levelToString(LogLevel level) {
    switch (level) {
        case LogLevel::DEBUG:   return "DEBUG";
        case LogLevel::INFO:    return "INFO";
        case LogLevel::WARNING: return "WARNING";
        case LogLevel::ERROR:   return "ERROR";
        case LogLevel::FATAL:   return "FATAL";
        default:                return "UNKNOWN";
    }
}

const char* Logger::levelToColor(LogLevel level) {
    // Check if output is a TTY before enabling colors potentially
    bool useEffectiveColor = colorEnabled && isatty(fileno(stderr)); // Example check
    if (!useEffectiveColor) return "";

    switch (level) {
        case LogLevel::DEBUG:   return "\033[38;5;250m"; // Lighter gray
        case LogLevel::INFO:    return "\033[38;5;44m";  // Sea green
        case LogLevel::WARNING: return "\033[38;5;208m"; // Soft orange
        case LogLevel::ERROR:   return "\033[38;5;203m"; // Soft red
        case LogLevel::FATAL:   return "\033[38;5;199m"; // Soft magenta
        default:                return "\033[0m";        // Reset
    }
}

std::string Logger::formatMessage(LogLevel level, const std::string& message) {
    const char* levelStr = levelToString(level);
    std::stringstream ss;
    const char* color = levelToColor(level);
    const char* reset = (color[0] == '\0') ? "" : "\033[0m"; // Only reset if color was applied
    ss << color << "[" << levelStr << "]" << reset << " " << message;
    return ss.str();
}

void Logger::log(LogLevel level, const std::string& message) {
    if (level < minLevel) return;
    std::lock_guard<std::mutex> lock(logMutex);

    std::ostream& target_out = (level >= LogLevel::WARNING) ? err_out : out;

    std::string formattedMessage = formatMessage(level, message);

    if (activeProgressBar) {
        if (level < LogLevel::ERROR) {
            bufferedMessages.push_back({level, message});
        } else {
            activeProgressBar->temporarilyPauseDisplay([this, &formattedMessage, &target_out]() {
                flushBufferedMessages();
                target_out << formattedMessage << std::endl;
            });
        }
    } else {
        target_out << formattedMessage << std::endl;
    }
}

void Logger::debug(const std::string& message) { log(LogLevel::DEBUG, message); }
void Logger::info(const std::string& message) { log(LogLevel::INFO, message); }
void Logger::warning(const std::string& message) { log(LogLevel::WARNING, message); }
void Logger::error(const std::string& message) { log(LogLevel::ERROR, message); }
void Logger::fatal(const std::string& message) { log(LogLevel::FATAL, message); }

void Logger::setMinLevel(LogLevel level) { minLevel = level; }
void Logger::enableColor(bool enable) { colorEnabled = enable; }

void Logger::setActiveProgressBar(ProgressBar* progressBar) {
    std::lock_guard<std::mutex> lock(logMutex);
    activeProgressBar = progressBar;
    bufferedMessages.clear();
}

void Logger::clearActiveProgressBar() {
    std::lock_guard<std::mutex> lock(logMutex);
    if (activeProgressBar) {
        flushBufferedMessages(); // Flush remaining messages
        activeProgressBar = nullptr;
    }
}

void Logger::flushBufferedMessages() {
    for (const auto& [level, message] : bufferedMessages) {
        std::string formattedMessage = formatMessage(level, message);
        std::ostream& target_out = (level >= LogLevel::WARNING) ? err_out : out;
        target_out << formattedMessage << std::endl;
    }
    bufferedMessages.clear();
}

// --- ProgressBar Implementation ---
ProgressBar::ProgressBar(size_t total, const std::string& description, int updateFrequencyMs)
    : total(total), current(0), description(description), active(false),
      updateFrequency(updateFrequencyMs), lastRenderCount(0), currentRate(0.0) {}

ProgressBar::~ProgressBar() {
    finish();
}

void ProgressBar::start() {
    if (active) return;
    startTime = std::chrono::steady_clock::now();
    lastRenderTime = startTime;
    lastRenderCount = 0;
    currentRate = 0.0;
    active = true;
    std::thread t(&ProgressBar::displayLoop, this);
    displayThread.swap(t);
    // Delay the first render to allow initialization to complete
    std::this_thread::sleep_for(std::chrono::milliseconds(10));
}

void ProgressBar::update(size_t increment) {
    size_t oldValue = current.load(std::memory_order_relaxed);
    current.fetch_add(increment, std::memory_order_relaxed);
    
    // Force render on significant progress changes (every 0.1% or more)
    size_t newValue = current.load(std::memory_order_relaxed);
    if (total > 0) {
        size_t oldPercent = (oldValue * 1000) / total;
        size_t newPercent = (newValue * 1000) / total;
        if (newPercent > oldPercent) {
            forceRender();
        }
    }
}

void ProgressBar::finish() {
    if (!active) return;
    active = false;
    if (displayThread.joinable()) {
        displayThread.join();
    }
    std::lock_guard<std::mutex> lock(mutex);
    current.store(total, std::memory_order_relaxed);
    clearLine();
    renderFinal();
    std::cout << std::endl;
    globalLogger.clearActiveProgressBar();
}

void ProgressBar::setDescription(const std::string& desc) {
    std::lock_guard<std::mutex> lock(mutex);
    description = desc;
}

void ProgressBar::setTotal(size_t newTotal) {
    std::lock_guard<std::mutex> lock(mutex);
    total = newTotal;
}

double ProgressBar::getProgress() const {
    size_t current_val = current.load(std::memory_order_relaxed);
    if (total == 0) return 0.0;
    return std::min(static_cast<double>(current_val) / total, 1.0);
}

std::chrono::milliseconds ProgressBar::getElapsedTime() const {
    auto now = std::chrono::steady_clock::now();
    return std::chrono::duration_cast<std::chrono::milliseconds>(now - startTime);
}

std::chrono::milliseconds ProgressBar::getEstimatedTimeRemaining() const {
    double progress = getProgress();
    auto elapsed_ms = getElapsedTime().count();
    if (progress <= 1e-6 || elapsed_ms <= 0) {
        return std::chrono::milliseconds(0);
    }
    double estimated_total_ms = static_cast<double>(elapsed_ms) / progress;
    double remaining_ms = estimated_total_ms - elapsed_ms;
    return std::chrono::milliseconds(static_cast<long long>(std::max(0.0, remaining_ms)));
}

void ProgressBar::displayLoop() {
    while (active) {
        render();
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
    }
}

void ProgressBar::clearLine() {
    std::cout << "\r\033[K" << std::flush;
}

int ProgressBar::getTerminalWidth() const {
    int termWidth = 80; // Default fallback width
    
    // Try using ioctl first (most reliable on Unix/Linux)
    struct winsize w;
    if (ioctl(STDOUT_FILENO, TIOCGWINSZ, &w) == 0 && w.ws_col > 0) {
        termWidth = w.ws_col;
    } else {
        // Try alternate methods
        const char* colsEnv = getenv("COLUMNS");
        if (colsEnv) {
            try {
                termWidth = std::stoi(colsEnv);
            } catch (...) {
                // Ignore conversion errors
            }
        } else {
            // Last resort: try using tput command
            FILE* pipe = popen("tput cols 2>/dev/null", "r");
            if (pipe) {
                char buffer[16];
                if (fgets(buffer, sizeof(buffer), pipe) != nullptr) {
                    try {
                        termWidth = std::stoi(buffer);
                    } catch (...) {
                        // Ignore conversion errors
                    }
                }
                pclose(pipe);
            }
        }
    }
    
    // Ensure reasonable width (between 40 and 250 columns)
    return std::max(40, std::min(termWidth, 250));
}

int ProgressBar::calculateBarWidth(int termWidth) const {
    // Adaptive calculation based on terminal width
    // - For narrow terminals (< 80 cols): use 20% of width
    // - For medium terminals (80-120 cols): use 30% of width
    // - For wide terminals (> 120 cols): use 40% of width
    float percentage;
    
    if (termWidth < 80) {
        percentage = 0.2f; // 20% for narrow terminals
    } else if (termWidth <= 120) {
        percentage = 0.3f; // 30% for medium terminals
    } else {
        percentage = 0.4f; // 40% for wide terminals
    }
    
    int calculatedWidth = static_cast<int>(termWidth * percentage);
    return std::max(getMinBarWidth(), std::min(getMaxBarWidth(), calculatedWidth));
}

void ProgressBar::render() {
    std::lock_guard<std::mutex> lock(mutex);
    auto now = std::chrono::steady_clock::now();
    auto timeSinceStart = std::chrono::duration_cast<std::chrono::milliseconds>(now - startTime);
    auto timeSinceLastRender = std::chrono::duration_cast<std::chrono::milliseconds>(now - lastRenderTime);
    size_t current_val = current.load(std::memory_order_relaxed);
    size_t itemsSinceLastRender = current_val - lastRenderCount;

    // Calculate rates with smoothing
    double overallRate = (timeSinceStart.count() > 1000)
        ? static_cast<double>(current_val) * 1000.0 / timeSinceStart.count()
        : 0.0;
        
    if (timeSinceLastRender.count() > 100 && itemsSinceLastRender > 0) {
        // Smooth the rate calculation
        double instantRate = static_cast<double>(itemsSinceLastRender) * 1000.0 / timeSinceLastRender.count();
        currentRate = (currentRate * 0.7) + (instantRate * 0.3);  // 70% old rate, 30% new rate
        lastRenderTime = now;
        lastRenderCount = current_val;
    } else if (timeSinceLastRender.count() > 500) {
        // If no updates for 0.5 seconds, gradually decrease the rate
        currentRate *= 0.98;  // Use 2% decay instead of 5% for smoother degradation
        lastRenderTime = now;
    }

    // Use a minimum rate to prevent display from dropping to zero
    double displayRate = std::max(currentRate, overallRate * 0.5);
    
    // Handle unknown total (max value) case
    bool isUnknownTotal = (total == std::numeric_limits<size_t>::max());
    double progress = isUnknownTotal ? 0.0 : (total > 0 ? static_cast<double>(current_val) / total : 0.0);

    // Get terminal width and calculate bar width
    int termWidth = getTerminalWidth();
    int barWidth = calculateBarWidth(termWidth);

    // Determine display mode based on available width
    enum class DisplayMode { Minimal, Compact, Normal, Full } mode;
    
    if (termWidth < 40) {
        mode = DisplayMode::Minimal;  // Just show percentage or simple indicator
    } else if (termWidth < 60) {
        mode = DisplayMode::Compact;  // Show percentage and simplified bar
    } else if (termWidth < 100) {
        mode = DisplayMode::Normal;   // Regular display without some extras
    } else {
        mode = DisplayMode::Full;     // Show all information
    }

    // Start building output
    std::stringstream ss;
    ss << "\r\033[K"; // Clear line

    if (mode == DisplayMode::Minimal) {
        // Extremely minimal display for very narrow terminals
        if (!isUnknownTotal) {
            int percentage = static_cast<int>(progress * 100);
            ss << percentage << "%";
        } else {
            ss << current_val;
        }
    } else {
        // For all other modes, include percentage and some form of bar
        if (!isUnknownTotal) {
            int percentage = static_cast<int>(progress * 100);
            ss << "\033[38;5;45m" << percentage << "%\033[0m" << " ";
            
            // Adjust bar width for compact mode
            int effectiveBarWidth = (mode == DisplayMode::Compact) ? 
                std::min(barWidth, 10) : barWidth;
                
            ss << getBlockBar(progress, effectiveBarWidth) << " ";
        } else {
            ss << "\033[38;5;45mProc\033[0m ";
        }
        
        // Add count information for all but minimal mode
        std::stringstream progressSection;
        if (isUnknownTotal) {
            progressSection << "\033[1m" << current_val << "\033[0m";
        } else {
            // In compact mode, only show current value, not total
            if (mode == DisplayMode::Compact) {
                progressSection << "\033[1m" << current_val << "\033[0m";
            } else {
                progressSection << "\033[1m" << current_val << "/" << total << "\033[0m";
            }
        }
        ss << progressSection.str() << " ";
        
        // Add rate for normal and full modes
        if (mode >= DisplayMode::Normal) {
            ss << "\033[38;5;208m" << std::fixed << std::setprecision(1) 
               << displayRate << "/s\033[0m" << " ";
        }
        
        // Add ETA for normal and full modes
        if (mode >= DisplayMode::Normal && !isUnknownTotal) {
            std::string eta = getETA();
            ss << "\033[38;5;105mETA: " << eta << "\033[0m";
        } else if (mode >= DisplayMode::Normal && isUnknownTotal) {
            ss << "\033[38;5;105m" << std::fixed << std::setprecision(1)
               << (timeSinceStart.count() / 1000.0) << "s\033[0m";
        }
    }
    
    std::cout << ss.str() << std::flush;
}

void ProgressBar::renderFinal() {
    auto elapsed = getElapsedTime();
    double seconds = elapsed.count() / 1000.0;
    size_t final_count = current.load(std::memory_order_relaxed);
    double overallItemsPerSecond = (seconds > 0.01)
        ? static_cast<double>(final_count) / seconds
        : 0.0;

    int termWidth = getTerminalWidth();
    
    // Determine display mode based on available width
    enum class DisplayMode { Minimal, Normal, Full } mode;
    
    if (termWidth < 40) {
        mode = DisplayMode::Minimal;
    } else if (termWidth < 80) {
        mode = DisplayMode::Normal;
    } else {
        mode = DisplayMode::Full;
    }
    
    std::stringstream ss;
    ss << "\r\033[K";
    
    // Always show completion mark and count
    ss << "\033[38;5;40m✓\033[0m ";
    
    bool isUnknownTotal = (total == std::numeric_limits<size_t>::max());
    if (isUnknownTotal) {
        ss << "\033[1m" << final_count << (mode == DisplayMode::Minimal ? "\033[0m" : " items\033[0m");
    } else {
        if (mode == DisplayMode::Minimal) {
            ss << "\033[1m" << final_count << "\033[0m";
        } else {
            ss << "\033[1m" << final_count << "/" << total << "\033[0m";
        }
    }
    
    // Show rate for normal and full modes
    if (mode >= DisplayMode::Normal) {
        ss << " \033[38;5;208m" << std::fixed << std::setprecision(1) 
           << overallItemsPerSecond << " it/s\033[0m";
    }
    
    // Show time for normal and full modes
    if (mode >= DisplayMode::Normal) {
        ss << " \033[38;5;105m" << std::fixed << std::setprecision(2) 
           << seconds << "s\033[0m";
    }
    
    // Add "Complete" text only for full mode
    if (mode == DisplayMode::Full) {
        ss << " \033[38;5;220mComplete\033[0m";
    }

    std::cout << ss.str();
}

void ProgressBar::forceRender() {
    render();
}

void ProgressBar::temporarilyPauseDisplay(const std::function<void()>& callback) {
    std::lock_guard<std::mutex> lock(mutex);
    clearLine();
    callback();
    render();
}

std::string ProgressBar::getColorGradient(double progress) const {
    int r = static_cast<int>(46 + (progress * 0));
    int g = static_cast<int>(126 + (progress * 129));
    int b = static_cast<int>(236 - (progress * 140));
    
    std::stringstream color;
    color << "\033[38;2;" << r << ";" << g << ";" << b << "m";
    return color.str();
}

std::string ProgressBar::getBlockBar(double progress, int width) const {
    double scaledProgress = progress * width;
    int fullBlocks = static_cast<int>(std::floor(scaledProgress));
    int remainder = static_cast<int>((scaledProgress - fullBlocks) * 8);
    
    const char* blocks[] = {" ", "▏", "▎", "▍", "▌", "▋", "▊", "▉", "█"};
    std::stringstream bar;
    
    if (barStyle == "gradient") {
        for (int i = 0; i < fullBlocks; i++) {
            double segmentProgress = static_cast<double>(i) / width;
            bar << getColorGradient(segmentProgress) << "█";
        }
        if (fullBlocks < width) {
            double segmentProgress = static_cast<double>(fullBlocks) / width;
            bar << getColorGradient(segmentProgress) << blocks[remainder] << "\033[0m";
            for (int i = fullBlocks + 1; i < width; i++) {
                bar << " ";
            }
        }
    } else {
        std::string color = getColorGradient(progress);
        bar << color;
        for (int i = 0; i < fullBlocks; i++) {
            bar << "█";
        }
        if (fullBlocks < width) {
            bar << blocks[remainder] << "\033[0m";
            for (int i = fullBlocks + 1; i < width; i++) {
                bar << " ";
            }
        } else {
            bar << "\033[0m";
        }
    }
    return bar.str();
}

std::string ProgressBar::getETA() const {
    auto remainingMs = getEstimatedTimeRemaining();
    auto seconds = remainingMs.count() / 1000;
    
    if (seconds <= 0) return "0s";
    
    std::stringstream eta;
    if (seconds >= 3600) {
        eta << seconds / 3600 << "h " << (seconds % 3600) / 60 << "m";
    } else if (seconds >= 60) {
        eta << seconds / 60 << "m " << seconds % 60 << "s";
    } else {
        eta << seconds << "s";
    }
    return eta.str();
}


// --- Molecule Implementation ---
Molecule::Molecule() : valid(false) {}

Molecule::Molecule(const std::string& smilesStr) : originalSmiles(smilesStr), valid(false) {
    parse(smilesStr);
}

Molecule::Molecule(RDKit::ROMol* rdkitMol, bool takeOwnership) : valid(false) {
    if (!rdkitMol) {
        errorMessage = "Input RDKit molecule pointer was null.";
        return;
    }
    try {
        if (takeOwnership) {
            mol = std::shared_ptr<RDKit::ROMol>(rdkitMol);
        } else {
            mol = std::shared_ptr<RDKit::ROMol>(new RDKit::RWMol(*rdkitMol));
        }
        smiles = RDKit::MolToSmiles(*mol);
        originalSmiles = smiles; // If created from RDKit mol, original IS the canonical form initially
        valid = true;
    } catch (const RDKit::MolSanitizeException& e) {
        errorMessage = "RDKit Sanity Exception during molecule creation: " + std::string(e.what());
        mol = nullptr;
        valid = false;
    } catch (const std::exception& e) {
        errorMessage = "Error creating molecule from RDKit molecule: " + std::string(e.what());
        mol = nullptr;
        valid = false;
    } catch (...) {
        errorMessage = "Unknown error creating molecule from RDKit molecule.";
        mol = nullptr;
        valid = false;
    }
}


bool Molecule::parse(const std::string& smilesStr) {
    originalSmiles = smilesStr; // Store the input smiles regardless of validity
    mol = nullptr;
    valid = false;
    errorMessage = "";

    if (smilesStr.empty()) {
        errorMessage = "Input SMILES string is empty.";
        return false;
    }

    try {
        RDKit::RWMol* rawMol = RDKit::SmilesToMol(smilesStr);
        if (!rawMol) {
            errorMessage = "RDKit failed to parse SMILES (returned null).";
            return false;
        }
        mol.reset(rawMol);
        smiles = RDKit::MolToSmiles(*mol); // Generate canonical SMILES after successful parse
        valid = true;
        return true;
    } catch (const RDKit::MolSanitizeException& e) {
        errorMessage = "RDKit Sanity Exception during SMILES parse: " + std::string(e.what());
        mol = nullptr; // Ensure mol is null on failure
        valid = false;
        return false;
    } catch (const std::exception& e) {
        errorMessage = "Error parsing SMILES: " + std::string(e.what());
        mol = nullptr;
        valid = false;
        return false;
    } catch (...) {
        errorMessage = "Unknown error parsing SMILES.";
        mol = nullptr;
        valid = false;
        return false;
    }
}

bool Molecule::isValid() const { return valid; }
const std::string& Molecule::getErrorMessage() const { return errorMessage; }
std::shared_ptr<RDKit::ROMol> Molecule::getMolecule() const { return mol; }
const std::string& Molecule::getSmiles() const { return smiles; }
const std::string& Molecule::getOriginalSmiles() const { return originalSmiles; }


void Molecule::sanitize() {
    if (!valid || !mol) return;
    try {
        RDKit::MolOps::sanitizeMol(*(RDKit::RWMol*)mol.get());
        smiles = RDKit::MolToSmiles(*mol); // Update canonical smiles after sanitization
    } catch (const RDKit::MolSanitizeException& e) {
        errorMessage = "Error sanitizing molecule: " + std::string(e.what());
        valid = false; // Mark as invalid if sanitization fails
    } catch (const std::exception& e) {
        errorMessage = "Error sanitizing molecule: " + std::string(e.what());
        valid = false;
    }
}

void Molecule::canonicalize() {
    if (!valid || !mol) return;
    try {
        RDKit::MolStandardize::CleanupParameters params;
        // Run cleanup directly on the RWMol if possible, or create a new one
        RDKit::RWMol* rwMol = dynamic_cast<RDKit::RWMol*>(mol.get());
        if (rwMol) {
            std::unique_ptr<RDKit::RWMol> cleaned(RDKit::MolStandardize::cleanup(*rwMol, params));
            if (cleaned) {
                 mol.reset(new RDKit::ROMol(*cleaned)); // Store as ROMol after cleanup
                 smiles = RDKit::MolToSmiles(*mol); // Update canonical smiles
            } else {
                 globalLogger.warning("Canonicalization returned null for: " + originalSmiles);
            }
        } else {
            // If it's already a ROMol, we might need to copy to RWMol first if cleanup requires it
             std::unique_ptr<RDKit::RWMol> tempMol(new RDKit::RWMol(*mol));
             std::unique_ptr<RDKit::RWMol> cleaned(RDKit::MolStandardize::cleanup(*tempMol, params));
             if (cleaned) {
                 mol.reset(new RDKit::ROMol(*cleaned));
                 smiles = RDKit::MolToSmiles(*mol);
             } else {
                 globalLogger.warning("Canonicalization returned null for: " + originalSmiles);
             }
        }

    } catch (const std::exception& e) {
        errorMessage = "Error canonicalizing molecule: " + std::string(e.what());
        globalLogger.warning(errorMessage + " for SMILES: " + originalSmiles);
        // Optionally mark as invalid? Depends on desired behavior.
        // valid = false;
    }
}


bool Molecule::hasProperty(const std::string& name) const {
    return properties.count(name);
}

std::string Molecule::getProperty(const std::string& name) const {
    auto it = properties.find(name);
    return (it != properties.end()) ? it->second : "";
}

void Molecule::setProperty(const std::string& name, const std::string& value) {
    properties[name] = value;
}

const std::unordered_map<std::string, std::string>& Molecule::getProperties() const {
    return properties;
}


int Molecule::getNumAtoms() const { return (valid && mol) ? mol->getNumAtoms() : 0; }
int Molecule::getNumBonds() const { return (valid && mol) ? mol->getNumBonds() : 0; }

std::vector<int> Molecule::getAtomicNumbers() const {
    std::vector<int> atomicNumbers;
    if (!valid || !mol) return atomicNumbers;
    atomicNumbers.reserve(mol->getNumAtoms());
    for (const auto atom : mol->atoms()) {
        atomicNumbers.push_back(atom->getAtomicNum());
    }
    return atomicNumbers;
}

std::vector<int> Molecule::getBondTypes() const {
    std::vector<int> bondTypes;
    if (!valid || !mol) return bondTypes;
    bondTypes.reserve(mol->getNumBonds());
    for (const auto bond : mol->bonds()) {
        bondTypes.push_back(static_cast<int>(bond->getBondType())); // Cast BondType enum
    }
    return bondTypes;
}

std::string Molecule::toJSON() const {
    rapidjson::StringBuffer buffer;
    rapidjson::Writer<rapidjson::StringBuffer> writer(buffer);
    writer.StartObject();
    writer.Key("smiles"); writer.String(smiles.c_str());
    writer.Key("original_smiles"); writer.String(originalSmiles.c_str());
    writer.Key("valid"); writer.Bool(valid);
    if (!valid) {
        writer.Key("error_message"); writer.String(errorMessage.c_str());
    }
    writer.Key("properties");
    writer.StartObject();
    for (const auto& prop : properties) {
        writer.Key(prop.first.c_str()); writer.String(prop.second.c_str());
    }
    writer.EndObject();
    writer.EndObject();
    return buffer.GetString();
}

Molecule Molecule::fromJSON(const std::string& json) {
    rapidjson::Document document;
    document.Parse(json.c_str());
    Molecule molInstance; // Changed variable name from 'mol'

    if (document.HasMember("smiles") && document["smiles"].IsString()) {
        molInstance.smiles = document["smiles"].GetString();
    }
    if (document.HasMember("original_smiles") && document["original_smiles"].IsString()) {
        molInstance.originalSmiles = document["original_smiles"].GetString();
    }
    if (document.HasMember("valid") && document["valid"].IsBool()) {
        molInstance.valid = document["valid"].GetBool();
    }
    if (document.HasMember("error_message") && document["error_message"].IsString()) {
        molInstance.errorMessage = document["error_message"].GetString();
    }
    if (document.HasMember("properties") && document["properties"].IsObject()) {
        const auto& props = document["properties"];
        for (auto it = props.MemberBegin(); it != props.MemberEnd(); ++it) {
            if (it->name.IsString() && it->value.IsString()) {
                molInstance.properties[it->name.GetString()] = it->value.GetString();
            }
        }
    }

    // Attempt to re-parse from the canonical SMILES if marked valid
    // If parsing fails here, it might indicate an issue with the stored SMILES
    if (molInstance.valid && !molInstance.smiles.empty()) {
        if (!molInstance.parse(molInstance.smiles)) {
             globalLogger.warning("Failed to re-parse valid molecule from JSON SMILES: " + molInstance.smiles);
             // Keep valid=true based on JSON, but log inconsistency.
        }
    } else if (!molInstance.originalSmiles.empty()) {
        // If not valid or smiles empty, try parsing originalSmiles for completeness
        molInstance.parse(molInstance.originalSmiles);
    }


    return molInstance;
}

// --- MoleculeBatch Implementation ---
MoleculeBatch::MoleculeBatch(size_t initialSize) {
    molecules.reserve(initialSize);
}

bool MoleculeBatch::add(const Molecule& mol) {
    std::lock_guard<std::mutex> lock(mutex);
    if (molecules.size() >= capacity) {
         globalLogger.warning("MoleculeBatch::add: Exceeded planned capacity (" + std::to_string(capacity) + "). Performance might degrade.");
         // Allow adding beyond capacity but warn
    }
    molecules.push_back(mol);
    return true;
}

bool MoleculeBatch::add(const std::string& smilesStr) {
    if (smilesStr.empty() || smilesStr == "SMILES" || smilesStr == "smiles") { // Simple header check
        return false;
    }
    Molecule mol; // Create on stack
    if (mol.parse(smilesStr)) {
        return add(mol); // Move semantics handled by push_back if mol is rvalue
    } else {
        globalLogger.debug("Skipping invalid SMILES: " + smilesStr + " - Error: " + mol.getErrorMessage());
        // Optionally create an "invalid" molecule entry to preserve row count?
        // add(mol); // Add the invalid molecule object
        return false; // Indicate failure to parse/add valid molecule
    }
}

void MoleculeBatch::addSmilesBatch(const std::vector<std::string>& smilesVec) {
    // Clear existing molecules and resize to match input lines
    molecules.clear();
    molecules.resize(smilesVec.size());

    #ifdef WITH_TBB
    // Use parallel_for to parse directly into the pre-sized vector elements
    tbb::parallel_for(tbb::blocked_range<size_t>(0, smilesVec.size()),
        [&](const tbb::blocked_range<size_t>& range) {
            for (size_t i = range.begin(); i != range.end(); ++i) {
                const std::string& smilesStr = smilesVec[i];
                molecules[i].parse(smilesStr);
                // Keep original smiles for invalid ones
                if (!molecules[i].isValid()) {
                    molecules[i].setOriginalSmiles(smilesStr);
                    molecules[i].setErrorMessage("Failed to parse SMILES");
                }
            }
        });
    #else
    for (size_t i = 0; i < smilesVec.size(); ++i) {
        molecules[i].parse(smilesVec[i]);
        if (!molecules[i].isValid()) {
            molecules[i].setOriginalSmiles(smilesVec[i]);
            molecules[i].setErrorMessage("Failed to parse SMILES (single-threaded)");
        }
    }
    #endif
}


void MoleculeBatch::incrementProcessed(size_t count) {
    processed.fetch_add(count, std::memory_order_relaxed);
}

size_t MoleculeBatch::size() const {
    std::lock_guard<std::mutex> lock(mutex); // Need lock if accessed concurrently
    return molecules.size();
}
size_t MoleculeBatch::getCapacity() const { return capacity; }
size_t MoleculeBatch::getProcessed() const { return processed.load(std::memory_order_relaxed); }

const std::vector<Molecule>& MoleculeBatch::getMolecules() const {
    // Returning const ref is generally safe without lock if additions don't invalidate refs
    // but locking provides stronger guarantees if needed.
    // std::lock_guard<std::mutex> lock(mutex);
    return molecules;
}
std::vector<Molecule>& MoleculeBatch::getMolecules() {
    // Returning non-const ref requires careful handling by the caller.
    // std::lock_guard<std::mutex> lock(mutex);
    return molecules;
}

void MoleculeBatch::clear() {
    std::lock_guard<std::mutex> lock(mutex);
    molecules.clear();
    processed.store(0, std::memory_order_relaxed);
}

bool MoleculeBatch::isFull() const {
     std::lock_guard<std::mutex> lock(mutex); // Need lock if checked concurrently with add
    return molecules.size() >= capacity;
}

// Implement MoleculeStream for row-by-row processing

MoleculeStream::MoleculeStream(const std::string& filename) 
    : processedCount(0), isOpen(false) {
    open(filename);
}

MoleculeStream::~MoleculeStream() {
    close();
}

bool MoleculeStream::open(const std::string& filename) {
    // Close existing file if open
    if (isOpen) {
        close();
    }
    
    inputFile.open(filename);
    isOpen = inputFile.is_open();
    
    if (!isOpen) {
        globalLogger.error("Failed to open file: " + filename);
    }
    
    processedCount = 0;
    return isOpen;
}

void MoleculeStream::close() {
    if (isOpen) {
        inputFile.close();
        isOpen = false;
    }
}

bool MoleculeStream::hasNext() {
    if (!isOpen) return false;
    return !inputFile.eof();
}

bool MoleculeStream::next(Molecule& molecule) {
    if (!isOpen) return false;
    
    if (std::getline(inputFile, currentLine)) {
        // Skip empty lines
        if (currentLine.empty()) {
            return next(molecule);
        }
        
        // Parse the SMILES from the line
        // Assuming format: "SMILES [optional_data]"
        std::string smiles = currentLine;
        
        // If line contains whitespace, extract just the SMILES part
        size_t firstSpace = currentLine.find_first_of(" \t");
        if (firstSpace != std::string::npos) {
            smiles = currentLine.substr(0, firstSpace);
        }
        
        // Parse the molecule
        bool success = molecule.parse(smiles);
        
        // If parsing succeeds, update the processed count
        if (success) {
            processedCount++;
            
            // Optionally store the whole line as a property
            molecule.setProperty("original_line", currentLine);
        } else {
            globalLogger.warning("Failed to parse SMILES: " + smiles);
        }
        
        return success;
    }
    
    return false;
}

bool MoleculeStream::skipHeader() {
    if (!isOpen) return false;
    
    if (std::getline(inputFile, currentLine)) {
        return true;
    }
    
    return false;
}

bool MoleculeStream::nextWithOriginalData(Molecule& molecule, std::vector<std::string>& rowData, 
                                          const std::string& delimiter, int smilesIndex, const std::string& escapeChar) {
    if (!isOpen) return false;
    
    if (std::getline(inputFile, currentLine)) {
        // Skip empty or whitespace-only lines
        if (currentLine.find_first_not_of(" \t\r\n") == std::string::npos) {
            return nextWithOriginalData(molecule, rowData, delimiter, smilesIndex, escapeChar);
        }
        
        // Parse CSV line into columns
        rowData = CsvIO::parseCsvLine(currentLine, delimiter, escapeChar);
        
        // Check if SMILES column index is valid
        if (smilesIndex >= 0 && static_cast<size_t>(smilesIndex) < rowData.size()) {
            // Parse the molecule from the SMILES column
            bool success = molecule.parse(rowData[smilesIndex]);
            
            // Update processed count on success
            if (success) {
                processedCount++;
            } else {
                globalLogger.warning("Failed to parse SMILES: " + rowData[smilesIndex]);
                // Still return true because we have data, even if molecule is invalid
            }
            
            return true;
        } else {
            globalLogger.warning("Invalid SMILES column index: " + std::to_string(smilesIndex));
            // Return true but molecule will be invalid
            return true;
        }
    }
    
    return false;
}

size_t MoleculeStream::getProcessedCount() const {
    return processedCount;
}

} // namespace desfact 