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
    : QMainWindow(parent)
{
    // Initialize backend components (if needed globally)
    // descriptorFactory = std::make_unique<DescriptorFactory>(); // Need to include header

    setupUi();
    populateDescriptorList();
    setupConnections();
}

MainWindow::~MainWindow()
{
    // Clean up resources if necessary
}

void MainWindow::setupUi()
{
    QWidget *centralWidget = new QWidget(this);
    setCentralWidget(centralWidget);

    QVBoxLayout *mainLayout = new QVBoxLayout(centralWidget);

    // --- Input/Output Group ---
    QGroupBox *ioGroup = new QGroupBox("Input / Output");
    QFormLayout *ioLayout = new QFormLayout(ioGroup);

    inputPathEdit = new QLineEdit;
    QPushButton *inputBrowseButton = new QPushButton("Browse...");
    QHBoxLayout *inputLayout = new QHBoxLayout;
    inputLayout->addWidget(inputPathEdit);
    inputLayout->addWidget(inputBrowseButton);
    ioLayout->addRow("Input CSV:", inputLayout);

    outputPathEdit = new QLineEdit;
    QPushButton *outputBrowseButton = new QPushButton("Browse...");
    QHBoxLayout *outputLayout = new QHBoxLayout;
    outputLayout->addWidget(outputPathEdit);
    outputLayout->addWidget(outputBrowseButton);
    ioLayout->addRow("Output CSV:", outputLayout);

    smilesColumnEdit = new QLineEdit("SMILES"); // Default
    ioLayout->addRow("SMILES Column:", smilesColumnEdit);

    delimiterCombo = new QComboBox;
    delimiterCombo->addItems({",", "\t", ";", "|"});
    delimiterCombo->setEditable(true); // Allow custom delimiters
    ioLayout->addRow("Delimiter:", delimiterCombo);

    mainLayout->addWidget(ioGroup);

    // --- Descriptor Selection Group ---
    QGroupBox *descGroup = new QGroupBox("Descriptors");
    QVBoxLayout *descLayout = new QVBoxLayout(descGroup);
    descriptorListWidget = new QListWidget;
    descriptorListWidget->setSelectionMode(QAbstractItemView::MultiSelection);
    descLayout->addWidget(descriptorListWidget);
    // Add Select All/None buttons? (Optional)
    mainLayout->addWidget(descGroup);

    // --- Options Group ---
    QGroupBox *optionsGroup = new QGroupBox("Options");
    QGridLayout *optionsLayout = new QGridLayout(optionsGroup);

    optionsLayout->addWidget(new QLabel("Threads:"), 0, 0);
    threadsSpinBox = new QSpinBox;
    threadsSpinBox->setRange(0, QThread::idealThreadCount() * 2); // 0 for auto
    threadsSpinBox->setValue(QThread::idealThreadCount() > 1 ? QThread::idealThreadCount() -1 : 1); // Default to n-1 cores
    optionsLayout->addWidget(threadsSpinBox, 0, 1);

    optionsLayout->addWidget(new QLabel("Batch Size:"), 1, 0);
    batchSizeSpinBox = new QSpinBox;
    batchSizeSpinBox->setRange(1, 100000);
    batchSizeSpinBox->setValue(1000); // Default
    optionsLayout->addWidget(batchSizeSpinBox, 1, 1);

    verboseCheckBox = new QCheckBox("Verbose Logging");
    optionsLayout->addWidget(verboseCheckBox, 0, 2);

    streamModeCheckBox = new QCheckBox("Stream Mode (Lower Memory)");
    streamModeCheckBox->setChecked(true); // Default to streaming
    optionsLayout->addWidget(streamModeCheckBox, 1, 2);

    mainLayout->addWidget(optionsGroup);


    // --- Log Output ---
    QGroupBox *logGroup = new QGroupBox("Log Output");
    QVBoxLayout *logLayout = new QVBoxLayout(logGroup);
    logTextEdit = new QPlainTextEdit;
    logTextEdit->setReadOnly(true);
    logLayout->addWidget(logTextEdit);
    mainLayout->addWidget(logGroup);


    // --- Progress and Run ---
    QHBoxLayout *progressRunLayout = new QHBoxLayout;
    progressBar = new QProgressBar;
    progressBar->setRange(0, 100);
    progressBar->setValue(0);
    progressRunLayout->addWidget(progressBar, 1); // Give progress bar more space

    runButton = new QPushButton("Run Calculation");
    progressRunLayout->addWidget(runButton);
    mainLayout->addLayout(progressRunLayout);

    // --- Status Bar ---
    statusLabel = new QLabel("Ready");
    statusBar()->addWidget(statusLabel);


    // Adjust main layout margins
    mainLayout->setContentsMargins(10, 10, 10, 10);
    mainLayout->setSpacing(15);
}

void MainWindow::populateDescriptorList()
{
    descriptorListWidget->clear();
    desfact::DescriptorFactory factory; // Temporary factory to list descriptors
    std::vector<std::string> descriptors = factory.getAvailableDescriptors();
    std::sort(descriptors.begin(), descriptors.end()); // Sort alphabetically

    for (const auto& name : descriptors) {
        QListWidgetItem *item = new QListWidgetItem(QString::fromStdString(name), descriptorListWidget);
        item->setFlags(item->flags() | Qt::ItemIsUserCheckable);
        item->setCheckState(Qt::Unchecked);
    }
}

void MainWindow::setupConnections()
{
    // File Browsing
    QPushButton *inputBrowseButton = nullptr;
    QPushButton *outputBrowseButton = nullptr;
    
    // Find browse buttons by traversing child widgets
    for (QObject *child : findChildren<QObject*>()) {
        if (QPushButton *btn = qobject_cast<QPushButton*>(child)) {
            if (btn->text() == "Browse...") {
                if (inputBrowseButton == nullptr) {
                    inputBrowseButton = btn;
                } else {
                    outputBrowseButton = btn;
                    break;
                }
            }
        }
    }
    
    // Connect browse buttons
    if (inputBrowseButton) {
        connect(inputBrowseButton, &QPushButton::clicked, this, &MainWindow::browseInputFile);
    }
    if (outputBrowseButton) {
        connect(outputBrowseButton, &QPushButton::clicked, this, &MainWindow::browseOutputFile);
    }

    // Run Button
    connect(runButton, &QPushButton::clicked, this, &MainWindow::runCalculation);
}

void MainWindow::browseInputFile()
{
    QString filePath = QFileDialog::getOpenFileName(this, "Select Input CSV File", "", "CSV Files (*.csv *.tsv *.smi);;All Files (*)");
    if (!filePath.isEmpty()) {
        inputPathEdit->setText(filePath);
    }
}

void MainWindow::browseOutputFile()
{
    QString filePath = QFileDialog::getSaveFileName(this, "Select Output CSV File", "", "CSV Files (*.csv);;All Files (*)");
    if (!filePath.isEmpty()) {
        // Ensure .csv extension if not provided
        if (!filePath.endsWith(".csv", Qt::CaseInsensitive)) {
            filePath += ".csv";
        }
        outputPathEdit->setText(filePath);
    }
}

void MainWindow::appendLogMessage(const QString& message)
{
    // Ensure this runs on the GUI thread if called from another thread
     if (QThread::currentThread() != thread()) {
         QMetaObject::invokeMethod(this, "appendLogMessage", Qt::QueuedConnection,
                                   Q_ARG(QString, message));
         return;
     }
     logTextEdit->appendPlainText(QTime::currentTime().toString("hh:mm:ss.zzz") + " - " + message);
}

void MainWindow::updateProgressBar(double value)
{
     // Ensure this runs on the GUI thread if called from another thread
     if (QThread::currentThread() != thread()) {
          QMetaObject::invokeMethod(this, "updateProgressBar", Qt::QueuedConnection,
                                   Q_ARG(double, value));
          return;
     }
     progressBar->setValue(static_cast<int>(value * 100));
}


void MainWindow::runCalculation()
{
    // --- Get settings from UI ---
    QString inputPath = inputPathEdit->text();
    QString outputPath = outputPathEdit->text();
    QString smilesCol = smilesColumnEdit->text();
    QString delimiter = delimiterCombo->currentText();
    // bool hasHeader = !findChild<QCheckBox*>("noHeaderCheckBox")->isChecked(); // Assuming you add this checkbox
    bool verbose = verboseCheckBox->isChecked();
    bool streamMode = streamModeCheckBox->isChecked();
    int threads = threadsSpinBox->value();
    int batchSize = batchSizeSpinBox->value();

    selectedDescriptors.clear();
    for (int i = 0; i < descriptorListWidget->count(); ++i) {
        QListWidgetItem *item = descriptorListWidget->item(i);
        if (item->checkState() == Qt::Checked) {
            selectedDescriptors.append(item->text());
        }
    }

    // --- Basic Validation ---
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

    // --- Disable UI and Run ---
    runButton->setEnabled(false);
    progressBar->setValue(0);
    logTextEdit->clear();
    statusLabel->setText("Processing...");
    appendLogMessage("Starting calculation...");


    // --- Execute Backend (Placeholder) ---
    // IMPORTANT: For a real application, run this in a separate QThread
    // to avoid freezing the GUI. Connect signals from the worker thread
    // back to the updateProgressBar and appendLogMessage slots.

    // For now, call the placeholder directly (will freeze GUI during execution)
    bool hasHeader = true; // Assume header for now, add checkbox later
    executeBackendCalculation(
        inputPath, outputPath, smilesCol, delimiter, hasHeader,
        selectedDescriptors, threads, batchSize, verbose, streamMode,
        // Lambda functions to pass progress/log callbacks
        [this](double progress) { this->updateProgressBar(progress); },
        [this](const QString& msg) { this->appendLogMessage(msg); }
    );


    // --- Re-enable UI ---
    runButton->setEnabled(true);
    statusLabel->setText("Calculation finished.");
    appendLogMessage("Calculation finished.");
    QMessageBox::information(this, "Finished", "Descriptor calculation complete.");


}