#include "mainwindow.h"  // Include from local directory
#include "descriptors.hpp" // Include your descriptor factory header
#include <QFileDialog>
#include <QMessageBox>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QGroupBox>
#include <QFormLayout>
#include <QThread> // For potential background worker
#include <QCoreApplication> // For processEvents
#include <QScrollBar>
#include <QDesktopWidget>
#include <QScreen>
#include <QPushButton>
#include <QLabel>
#include <QHeaderView>
#include <QTime>
#include <QFontDatabase>

// Worker class for background processing
class CalculationWorker : public QObject {
    Q_OBJECT
public:
    CalculationWorker(
        const QString& inputPath, const QString& outputPath, 
        const QString& smilesColumn, const QString& delimiter,
        bool hasHeader, const QStringList& descriptors,
        int threads, int batchSize, bool verbose, bool streamMode)
        : inputPath(inputPath), outputPath(outputPath), 
          smilesColumn(smilesColumn), delimiter(delimiter),
          hasHeader(hasHeader), descriptors(descriptors),
          threads(threads), batchSize(batchSize), 
          verbose(verbose), streamMode(streamMode) {}

public slots:
    void process() {
        emit logMessage("Starting descriptor calculation...");
        emit logMessage("Input: " + inputPath);
        emit logMessage("Output: " + outputPath);
        emit logMessage("SMILES Column: " + smilesColumn);
        emit logMessage("Descriptors: " + descriptors.join(", "));
        emit logMessage("Threads: " + QString::number(threads));
        emit logMessage("Batch Size: " + QString::number(batchSize));
        
        try {
            // Convert QStringList to std::vector<std::string>
            std::vector<std::string> descriptorList;
            for (const QString& desc : descriptors) {
                descriptorList.push_back(desc.toStdString());
            }
            
            // Create factory instance
            desfact::DescriptorFactory factory;
            
            // TODO: Implement actual processing using your backend classes
            // This is a placeholder simulation with progress updates
            int totalSteps = 100;
            for (int i = 0; i <= totalSteps; ++i) {
                if (QThread::currentThread()->isInterruptionRequested()) {
                    emit logMessage("Calculation interrupted");
                    return;
                }
                
                if (i % 10 == 0) {
                    emit logMessage("Processing step " + QString::number(i) + "...");
                }
                emit progressUpdated(static_cast<double>(i) / totalSteps);
                QThread::msleep(50); // Simulate work
            }
            
            emit logMessage("Calculation completed successfully");
        } catch (const std::exception& e) {
            emit logMessage("Error: " + QString::fromStdString(e.what()));
        } catch (...) {
            emit logMessage("Unknown error occurred");
        }
        
        emit progressUpdated(1.0);
        emit finished();
    }
    
signals:
    void progressUpdated(double progress);
    void logMessage(const QString& message);
    void finished();
    
private:
    QString inputPath;
    QString outputPath;
    QString smilesColumn;
    QString delimiter;
    bool hasHeader;
    QStringList descriptors;
    int threads;
    int batchSize;
    bool verbose;
    bool streamMode;
};

// Placeholder for backend execution logic
// This would ideally be in a separate class or thread
void executeBackendCalculation(
    const QString& inputPath,
    const QString& outputPath,
    const QString& smilesColumn,
    const QString& delimiter,
    bool hasHeader,
    const QStringList& descriptors,
    int threads,
    int batchSize,
    bool verbose,
    bool streamMode,
    std::function<void(double)> progressCallback,
    std::function<void(const QString&)> logCallback)
{
    // --- Placeholder Logic ---
    logCallback("Starting backend calculation (placeholder)...");
    logCallback("Input: " + inputPath);
    logCallback("Output: " + outputPath);
    logCallback("SMILES Col: " + smilesColumn);
    logCallback("Descriptors: " + descriptors.join(", "));
    logCallback("Threads: " + QString::number(threads));
    logCallback("Batch Size: " + QString::number(batchSize));
    logCallback("Verbose: " + QString(verbose ? "true" : "false"));
    logCallback("Stream Mode: " + QString(streamMode ? "true" : "false"));


    try {
        // TODO: Instantiate CsvIO, DescriptorFactory etc.
        // TODO: Adapt the logic from your existing src/main.cpp
        //       - Read configuration from GUI elements
        //       - Create CsvIO instance
        //       - Create LineReader/MoleculeStream/Batch processor
        //       - Create DescriptorFactory instance
        //       - Calculate descriptors
        //       - Create ResultWriter instance
        //       - Write results
        //       - IMPORTANT: Emit progress updates using progressCallback
        //       - IMPORTANT: Emit log messages using logCallback
        // Example Progress Update Loop:
        int totalSteps = 100;
        for (int i = 0; i <= totalSteps; ++i) {
             if (i % 10 == 0) {
                 logCallback("Processing step " + QString::number(i) + "...");
             }
             progressCallback(static_cast<double>(i) / totalSteps);
             QThread::msleep(50); // Simulate work
              // Allow GUI to update (essential if running in main thread)
             // QCoreApplication::processEvents(); 
        }

        logCallback("Backend calculation finished (placeholder).");

    } catch (const std::exception& e) {
        logCallback("Error during backend calculation: " + QString::fromStdString(e.what()));
    } catch (...) {
        logCallback("Unknown error during backend calculation.");
    }
    // Ensure progress reaches 100%
    progressCallback(1.0);

}

// --- MainWindow Implementation ---

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent), workerThread(nullptr)
{
    // Initialize backend components (if needed globally)
    // descriptorFactory = std::make_unique<DescriptorFactory>(); // Need to include header

    setupUi();
    populateDescriptorList();
    setupConnections();
    setStyleSheet(getStyleSheet());
}

MainWindow::~MainWindow()
{
    if (workerThread && workerThread->isRunning()) {
        workerThread->requestInterruption();
        workerThread->quit();
        workerThread->wait();
    }
}

void MainWindow::setupUi()
{
    // Set window properties
    setWindowTitle("Descriptor Factory");
    QScreen *screen = QGuiApplication::primaryScreen();
    QRect screenGeometry = screen->availableGeometry();
    int width = qMin(1024, screenGeometry.width() - 200);
    int height = qMin(768, screenGeometry.height() - 100);
    resize(width, height);
    
    QWidget *centralWidget = new QWidget(this);
    setCentralWidget(centralWidget);
    
    QVBoxLayout *mainLayout = new QVBoxLayout(centralWidget);
    mainLayout->setSpacing(15);
    mainLayout->setContentsMargins(15, 15, 15, 15);
    
    // Top section with file selection
    QGroupBox *fileGroup = new QGroupBox("Input & Output Files");
    QGridLayout *fileLayout = new QGridLayout(fileGroup);
    fileLayout->setSpacing(10);
    
    QLabel *inputLabel = new QLabel("Input CSV:");
    inputPathEdit = new QLineEdit;
    inputPathEdit->setPlaceholderText("Select input CSV file");
    QPushButton *inputBrowseBtn = new QPushButton("Browse");
    inputBrowseBtn->setFixedWidth(100);
    
    QLabel *outputLabel = new QLabel("Output CSV:");
    outputPathEdit = new QLineEdit;
    outputPathEdit->setPlaceholderText("Select output CSV file");
    QPushButton *outputBrowseBtn = new QPushButton("Browse");
    outputBrowseBtn->setFixedWidth(100);
    
    fileLayout->addWidget(inputLabel, 0, 0);
    fileLayout->addWidget(inputPathEdit, 0, 1);
    fileLayout->addWidget(inputBrowseBtn, 0, 2);
    fileLayout->addWidget(outputLabel, 1, 0);
    fileLayout->addWidget(outputPathEdit, 1, 1);
    fileLayout->addWidget(outputBrowseBtn, 1, 2);
    
    mainLayout->addWidget(fileGroup);
    
    // Middle section with options and descriptors
    QHBoxLayout *middleLayout = new QHBoxLayout();
    
    // Left side - Options
    QGroupBox *optionsGroup = new QGroupBox("Options");
    QFormLayout *optionsLayout = new QFormLayout(optionsGroup);
    optionsLayout->setSpacing(10);
    
    smilesColumnEdit = new QLineEdit("SMILES");
    optionsLayout->addRow("SMILES Column:", smilesColumnEdit);
    
    delimiterCombo = new QComboBox();
    delimiterCombo->addItems({",", "\t", ";", "|"});
    delimiterCombo->setEditable(true);
    optionsLayout->addRow("Delimiter:", delimiterCombo);
    
    headerCheckBox = new QCheckBox("Has Header");
    headerCheckBox->setChecked(true);
    optionsLayout->addRow("", headerCheckBox);
    
    threadsSpinBox = new QSpinBox();
    threadsSpinBox->setRange(1, QThread::idealThreadCount() * 2);
    threadsSpinBox->setValue(QThread::idealThreadCount());
    optionsLayout->addRow("Threads:", threadsSpinBox);
    
    batchSizeSpinBox = new QSpinBox();
    batchSizeSpinBox->setRange(100, 10000);
    batchSizeSpinBox->setSingleStep(100);
    batchSizeSpinBox->setValue(1000);
    optionsLayout->addRow("Batch Size:", batchSizeSpinBox);
    
    verboseCheckBox = new QCheckBox();
    optionsLayout->addRow("Verbose Logging:", verboseCheckBox);
    
    streamModeCheckBox = new QCheckBox();
    streamModeCheckBox->setChecked(true);
    optionsLayout->addRow("Stream Mode:", streamModeCheckBox);
    
    middleLayout->addWidget(optionsGroup, 1);
    
    // Right side - Descriptors
    QGroupBox *descriptorsGroup = new QGroupBox("Descriptors");
    QVBoxLayout *descriptorsLayout = new QVBoxLayout(descriptorsGroup);
    
    // Add search box
    QHBoxLayout *searchLayout = new QHBoxLayout();
    QLabel *searchLabel = new QLabel("Search:");
    descriptorSearchEdit = new QLineEdit();
    descriptorSearchEdit->setPlaceholderText("Filter descriptors...");
    searchLayout->addWidget(searchLabel);
    searchLayout->addWidget(descriptorSearchEdit);
    descriptorsLayout->addLayout(searchLayout);
    
    // Add select all/none buttons
    QHBoxLayout *selectionLayout = new QHBoxLayout();
    QPushButton *selectAllBtn = new QPushButton("Select All");
    QPushButton *selectNoneBtn = new QPushButton("Select None");
    selectionLayout->addWidget(selectAllBtn);
    selectionLayout->addWidget(selectNoneBtn);
    descriptorsLayout->addLayout(selectionLayout);
    
    // Add list widget
    descriptorListWidget = new QTableWidget();
    descriptorListWidget->setColumnCount(2);
    descriptorListWidget->setHorizontalHeaderLabels({"Descriptor", "Selected"});
    descriptorListWidget->horizontalHeader()->setSectionResizeMode(0, QHeaderView::Stretch);
    descriptorListWidget->horizontalHeader()->setSectionResizeMode(1, QHeaderView::ResizeToContents);
    descriptorListWidget->verticalHeader()->setVisible(false);
    descriptorListWidget->setSelectionBehavior(QAbstractItemView::SelectRows);
    descriptorListWidget->setAlternatingRowColors(true);
    descriptorListWidget->setSortingEnabled(true);
    descriptorsLayout->addWidget(descriptorListWidget);
    
    middleLayout->addWidget(descriptorsGroup, 2);
    mainLayout->addLayout(middleLayout);
    
    // Bottom section with log and progress
    QGroupBox *logGroup = new QGroupBox("Log Output");
    QVBoxLayout *logLayout = new QVBoxLayout(logGroup);
    
    logTextEdit = new QPlainTextEdit();
    logTextEdit->setReadOnly(true);
    QFont monoFont = QFontDatabase::systemFont(QFontDatabase::FixedFont);
    logTextEdit->setFont(monoFont);
    logTextEdit->setMaximumBlockCount(1000);
    logLayout->addWidget(logTextEdit);
    
    mainLayout->addWidget(logGroup);
    
    // Progress and control buttons at the bottom
    QHBoxLayout *controlLayout = new QHBoxLayout();
    
    progressBar = new QProgressBar();
    progressBar->setRange(0, 100);
    progressBar->setValue(0);
    progressBar->setTextVisible(true);
    
    runButton = new QPushButton("Run Calculation");
    runButton->setMinimumWidth(150);
    stopButton = new QPushButton("Stop");
    stopButton->setMinimumWidth(100);
    stopButton->setEnabled(false);
    
    controlLayout->addWidget(progressBar, 3);
    controlLayout->addWidget(runButton, 1);
    controlLayout->addWidget(stopButton, 1);
    
    mainLayout->addLayout(controlLayout);
    
    // Status bar
    statusLabel = new QLabel("Ready");
    statusBar()->addWidget(statusLabel, 1);
    
    // Set layout proportions
    mainLayout->setStretchFactor(fileGroup, 1);
    mainLayout->setStretchFactor(middleLayout, 3);
    mainLayout->setStretchFactor(logGroup, 2);
    mainLayout->setStretchFactor(controlLayout, 1);
}

void MainWindow::populateDescriptorList()
{
    descriptorListWidget->clearContents();
    descriptorListWidget->setRowCount(0);
    
    try {
        desfact::DescriptorFactory factory;
        std::vector<std::string> descriptors = factory.getAvailableDescriptors();
        std::sort(descriptors.begin(), descriptors.end());
        
        descriptorListWidget->setRowCount(descriptors.size());
        
        for (size_t i = 0; i < descriptors.size(); i++) {
            QTableWidgetItem *nameItem = new QTableWidgetItem(QString::fromStdString(descriptors[i]));
            nameItem->setFlags(nameItem->flags() & ~Qt::ItemIsEditable);
            
            QWidget *checkboxWidget = new QWidget();
            QHBoxLayout *layout = new QHBoxLayout(checkboxWidget);
            QCheckBox *checkbox = new QCheckBox();
            layout->addWidget(checkbox);
            layout->setAlignment(Qt::AlignCenter);
            layout->setContentsMargins(0, 0, 0, 0);
            
            descriptorListWidget->setItem(i, 0, nameItem);
            descriptorListWidget->setCellWidget(i, 1, checkboxWidget);
        }
        
        appendLogMessage("Loaded " + QString::number(descriptors.size()) + " descriptors");
    } catch (const std::exception& e) {
        QMessageBox::critical(this, "Error", 
            "Failed to load descriptors: " + QString::fromStdString(e.what()));
    }
}

void MainWindow::setupConnections()
{
    // File browsing
    connect(findChild<QPushButton*>("Browse"), &QPushButton::clicked, 
            this, &MainWindow::browseInputFile);
    connect(findChild<QPushButton*>("Browse", Qt::FindChildrenRecursively), &QPushButton::clicked,
            this, &MainWindow::browseOutputFile);
    
    // Descriptor search filter
    connect(descriptorSearchEdit, &QLineEdit::textChanged,
            this, &MainWindow::filterDescriptors);
    
    // Select All/None buttons
    QPushButton *selectAllBtn = findChild<QPushButton*>("Select All");
    QPushButton *selectNoneBtn = findChild<QPushButton*>("Select None");
    
    connect(selectAllBtn, &QPushButton::clicked, this, [this]() {
        for (int i = 0; i < descriptorListWidget->rowCount(); i++) {
            if (!descriptorListWidget->isRowHidden(i)) {
                QWidget *widget = descriptorListWidget->cellWidget(i, 1);
                QCheckBox *checkbox = widget->findChild<QCheckBox*>();
                if (checkbox) checkbox->setChecked(true);
            }
        }
    });
    
    connect(selectNoneBtn, &QPushButton::clicked, this, [this]() {
        for (int i = 0; i < descriptorListWidget->rowCount(); i++) {
            if (!descriptorListWidget->isRowHidden(i)) {
                QWidget *widget = descriptorListWidget->cellWidget(i, 1);
                QCheckBox *checkbox = widget->findChild<QCheckBox*>();
                if (checkbox) checkbox->setChecked(false);
            }
        }
    });
    
    // Run and stop buttons
    connect(runButton, &QPushButton::clicked, this, &MainWindow::runCalculation);
    connect(stopButton, &QPushButton::clicked, this, &MainWindow::stopCalculation);
}

void MainWindow::browseInputFile()
{
    QString filePath = QFileDialog::getOpenFileName(
        this, "Select Input CSV File", "",
        "CSV Files (*.csv *.tsv *.smi);;All Files (*)");
    
    if (!filePath.isEmpty()) {
        inputPathEdit->setText(filePath);
    }
}

void MainWindow::browseOutputFile()
{
    QString filePath = QFileDialog::getSaveFileName(
        this, "Select Output CSV File", "",
        "CSV Files (*.csv);;All Files (*)");
    
    if (!filePath.isEmpty()) {
        if (!filePath.endsWith(".csv", Qt::CaseInsensitive)) {
            filePath += ".csv";
        }
        outputPathEdit->setText(filePath);
    }
}

void MainWindow::filterDescriptors(const QString &text)
{
    for (int i = 0; i < descriptorListWidget->rowCount(); i++) {
        QTableWidgetItem *item = descriptorListWidget->item(i, 0);
        bool match = text.isEmpty() || 
                    item->text().contains(text, Qt::CaseInsensitive);
        descriptorListWidget->setRowHidden(i, !match);
    }
}

void MainWindow::appendLogMessage(const QString &message)
{
    if (QThread::currentThread() != thread()) {
        QMetaObject::invokeMethod(this, "appendLogMessage", 
                                  Qt::QueuedConnection,
                                  Q_ARG(QString, message));
        return;
    }
    
    QString timestamp = QTime::currentTime().toString("[hh:mm:ss.zzz]");
    logTextEdit->appendPlainText(timestamp + " " + message);
    
    // Auto-scroll to bottom
    QScrollBar *scrollBar = logTextEdit->verticalScrollBar();
    scrollBar->setValue(scrollBar->maximum());
}

void MainWindow::updateProgressBar(double value)
{
    if (QThread::currentThread() != thread()) {
        QMetaObject::invokeMethod(this, "updateProgressBar", 
                                 Qt::QueuedConnection,
                                 Q_ARG(double, value));
        return;
    }
    
    int percent = static_cast<int>(value * 100);
    progressBar->setValue(percent);
    progressBar->setFormat(QString("%1%").arg(percent));
}

void MainWindow::runCalculation()
{
    // Get settings from UI
    QString inputPath = inputPathEdit->text();
    QString outputPath = outputPathEdit->text();
    QString smilesCol = smilesColumnEdit->text();
    QString delimiter = delimiterCombo->currentText();
    bool hasHeader = headerCheckBox->isChecked();
    bool verbose = verboseCheckBox->isChecked();
    bool streamMode = streamModeCheckBox->isChecked();
    int threads = threadsSpinBox->value();
    int batchSize = batchSizeSpinBox->value();
    
    // Collect selected descriptors
    QStringList selectedDescriptors;
    for (int i = 0; i < descriptorListWidget->rowCount(); i++) {
        if (descriptorListWidget->isRowHidden(i)) continue;
        
        QWidget *widget = descriptorListWidget->cellWidget(i, 1);
        QCheckBox *checkbox = widget->findChild<QCheckBox*>();
        if (checkbox && checkbox->isChecked()) {
            selectedDescriptors.append(descriptorListWidget->item(i, 0)->text());
        }
    }
    
    // Basic validation
    if (inputPath.isEmpty()) {
        QMessageBox::warning(this, "Input Error", "Please select an input file.");
        return;
    }
    if (outputPath.isEmpty()) {
        QMessageBox::warning(this, "Input Error", "Please select an output file.");
        return;
    }
    if (selectedDescriptors.isEmpty()) {
        QMessageBox::warning(this, "Input Error", "Please select at least one descriptor.");
        return;
    }
    if (delimiter.isEmpty()) {
        QMessageBox::warning(this, "Input Error", "Please specify a delimiter.");
        return;
    }
    
    // Disable UI and prepare for calculation
    runButton->setEnabled(false);
    stopButton->setEnabled(true);
    progressBar->setValue(0);
    logTextEdit->clear();
    statusLabel->setText("Processing...");
    appendLogMessage("Starting calculation...");
    
    // Create and start worker thread
    if (workerThread) {
        if (workerThread->isRunning()) {
            workerThread->requestInterruption();
            workerThread->quit();
            workerThread->wait();
        }
        delete workerThread;
    }
    
    workerThread = new QThread();
    CalculationWorker *worker = new CalculationWorker(
        inputPath, outputPath, smilesCol, delimiter, hasHeader,
        selectedDescriptors, threads, batchSize, verbose, streamMode
    );
    worker->moveToThread(workerThread);
    
    connect(workerThread, &QThread::started, worker, &CalculationWorker::process);
    connect(worker, &CalculationWorker::finished, this, &MainWindow::onCalculationFinished);
    connect(worker, &CalculationWorker::progressUpdated, this, &MainWindow::updateProgressBar);
    connect(worker, &CalculationWorker::logMessage, this, &MainWindow::appendLogMessage);
    connect(worker, &CalculationWorker::finished, worker, &QObject::deleteLater);
    connect(workerThread, &QThread::finished, workerThread, &QObject::deleteLater);
    
    workerThread->start();
}

void MainWindow::stopCalculation()
{
    if (workerThread && workerThread->isRunning()) {
        appendLogMessage("Stopping calculation...");
        workerThread->requestInterruption();
    }
}

void MainWindow::onCalculationFinished()
{
    runButton->setEnabled(true);
    stopButton->setEnabled(false);
    statusLabel->setText("Calculation finished");
    appendLogMessage("Calculation finished");
    
    QMessageBox::information(this, "Finished", 
                            "Descriptor calculation completed successfully.");
}

QString MainWindow::getStyleSheet() const
{
    return R"(
        QMainWindow, QDialog {
            background-color: #f5f5f5;
        }
        
        QGroupBox {
            font-weight: bold;
            border: 1px solid #cccccc;
            border-radius: 5px;
            margin-top: 1em;
            padding-top: 10px;
        }
        
        QGroupBox::title {
            subcontrol-origin: margin;
            left: 10px;
            padding: 0 5px;
        }
        
        QPushButton {
            background-color: #2980b9;
            color: white;
            border: none;
            border-radius: 4px;
            padding: 5px 15px;
        }
        
        QPushButton:hover {
            background-color: #3498db;
        }
        
        QPushButton:pressed {
            background-color: #1c5c84;
        }
        
        QPushButton:disabled {
            background-color: #cccccc;
            color: #666666;
        }
        
        QLineEdit, QComboBox, QSpinBox {
            border: 1px solid #cccccc;
            border-radius: 4px;
            padding: 5px;
            background-color: white;
        }
        
        QTableWidget {
            border: 1px solid #cccccc;
            background-color: white;
            alternate-background-color: #f0f0f0;
        }
        
        QHeaderView::section {
            background-color: #e0e0e0;
            padding: 5px;
            border: 1px solid #cccccc;
            font-weight: bold;
        }
        
        QProgressBar {
            border: 1px solid #cccccc;
            border-radius: 4px;
            text-align: center;
            background-color: white;
        }
        
        QProgressBar::chunk {
            background-color: #27ae60;
            width: 1px;
        }
        
        QPlainTextEdit {
            background-color: #2c3e50;
            color: #ecf0f1;
            border-radius: 4px;
            font-family: monospace;
        }
    )";
}

#include "mainwindow.moc"