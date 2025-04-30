#include "mainwindow.h"
#include <QApplication>
#include <QSplashScreen>
#include <QPixmap>
#include <QTimer>

int main(int argc, char *argv[])
{
    QApplication app(argc, argv);
    app.setApplicationName("Descriptor Factory");
    app.setOrganizationName("MedChemInfoLab");
    app.setApplicationVersion("1.0.0");
    
    // Create splash screen (optional)
    QPixmap pixmap(":/splash.png");
    QSplashScreen splash(pixmap.isNull() ? QPixmap(400, 300) : pixmap);
    
    if (pixmap.isNull()) {
        // Create default splash if no image available
        splash.setStyleSheet("background-color: #2c3e50; color: white; font-size: 24px;");
        splash.showMessage("Loading Descriptor Factory...", Qt::AlignCenter | Qt::AlignBottom, Qt::white);
    }
    
    splash.show();
    app.processEvents();
    
    // Set application style
    app.setStyle("Fusion");
    
    // Modern dark palette
    QPalette palette;
    palette.setColor(QPalette::Window, QColor(240, 240, 240));
    palette.setColor(QPalette::WindowText, QColor(50, 50, 50));
    palette.setColor(QPalette::Base, QColor(255, 255, 255));
    palette.setColor(QPalette::AlternateBase, QColor(245, 245, 245));
    palette.setColor(QPalette::Button, QColor(240, 240, 240));
    palette.setColor(QPalette::ButtonText, QColor(50, 50, 50));
    palette.setColor(QPalette::Link, QColor(42, 130, 218));
    palette.setColor(QPalette::Highlight, QColor(42, 130, 218));
    palette.setColor(QPalette::HighlightedText, QColor(255, 255, 255));
    app.setPalette(palette);
    
    // Load main window with short delay for splash effect
    MainWindow mainWindow;
    QTimer::singleShot(1500, &splash, [&]() {
        mainWindow.show();
        splash.finish(&mainWindow);
    });
    
    return app.exec();
}