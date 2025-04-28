#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QtWidgets> // Include necessary Qt Widgets headers

// Forward declarations if needed later
class DescriptorFactory;
class CsvIO;

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

private slots:
    void browseInputFile();
    void browseOutputFile();
    void runCalculation();
    void updateProgressBar(double value); // Slot to receive progress updates
    void appendLogMessage(const QString& message); // Slot to receive log messages

private:
    void setupUi();
    void populateDescriptorList();
    void setupConnections();

    // UI Elements
    QLineEdit *inputPathEdit;
    QLineEdit *outputPathEdit;
    QLineEdit *smilesColumnEdit;
    QComboBox *delimiterCombo;
    QListWidget *descriptorListWidget;
    QSpinBox *threadsSpinBox;
    QSpinBox *batchSizeSpinBox;
    QCheckBox *verboseCheckBox;
    QCheckBox *streamModeCheckBox; // Assuming default is streaming (no-stream flag)
    QPushButton *runButton;
    QPlainTextEdit *logTextEdit;
    QProgressBar *progressBar;
    QLabel *statusLabel;

    // Backend integration (placeholders)
    // std::unique_ptr<DescriptorFactory> descriptorFactory; // Need proper initialization
    // std::unique_ptr<CsvIO> csvHandler; // Need proper initialization
    QStringList selectedDescriptors; // Store selected descriptor names
};

#endif // MAINWINDOW_H 