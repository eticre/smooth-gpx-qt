#include <QApplication>
#include <QMainWindow>
#include <QFileDialog>
#include <QFileInfo>
#include <QSlider>
#include <QLabel>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QMessageBox>
#include <QPushButton>
#include <QDoubleSpinBox>
#include <QTimer>
#include <cmath>

#include "gpxtrack.h"
#include "graphwidget.h"

int main(int argc, char **argv)
{
    QApplication app(argc, argv);
    app.setApplicationName("GPX Gradient Smoother");

    app.setStyleSheet(R"QSS(
        QWidget {
            background-color: #181b21;
            color: #e8eaef;
        }
        QPushButton {
            background-color: #2d3442;
            border: 1px solid #3a4252;
            border-radius: 4px;
            padding: 5px 16px;
        }
        QPushButton:hover  { background-color: #374053; }
        QPushButton:pressed { background-color: #242a36; }
        QPushButton:checked { background-color: #1d4ed8; border: 1px solid #60a5fa; }
        QPushButton:checked { background-color: #1d4ed8; border: 1px solid #60a5fa; }
        QDoubleSpinBox {
            background-color: #2d3442;
            border: 1px solid #3a4252;
            border-radius: 4px;
            padding: 4px 8px;
            color: #e8eaef;
        }
        QDoubleSpinBox:focus  { border: 1px solid #60a5fa; }
        QSlider::groove:horizontal  { height: 4px; background: #3a4252; border-radius: 2px; }
        QSlider::sub-page:horizontal { background: #60a5fa; border-radius: 2px; }
        QSlider::add-page:horizontal { background: #3a4252; border-radius: 2px; }
        QSlider::handle:horizontal {
            width: 14px; margin: -6px 0; border-radius: 7px; background: #8ab4ff;
        }
        QSlider::add-line:horizontal, QSlider::sub-line:horizontal { width: 0px; }
        QSlider::tick { color: #3a4252; height: 6px; }
        QToolTip { background-color: #2d3442; color: #e8eaef; border: 1px solid #3a4252; }
    )QSS");

    QMainWindow win;
    win.setWindowTitle("GPX Gradient Smoother");
    win.resize(1000, 620);

    auto *central   = new QWidget;
    auto *vbox      = new QVBoxLayout(central);
    vbox->setContentsMargins(10, 10, 10, 10);
    vbox->setSpacing(8);

    // ---- toolbar row ----
    auto *topRow    = new QHBoxLayout;
    auto *openBtn   = new QPushButton("Open GPX…");
    auto *customBtn = new QPushButton("CustomReg");
    customBtn->setCheckable(true);
    customBtn->setToolTip("Attiva/disattiva la modifica manuale della curva\n(click & trascina sul grafico)");
    //customBtn->hide();   // temporaneamente nascosto; da ripulire/rimuovere dopo
    auto *saveBtn   = new QPushButton("Save");
    auto *fileLabel = new QLabel("no file loaded");
    fileLabel->setStyleSheet("color:#9aa3b2;");
    topRow->addWidget(openBtn);
    topRow->addWidget(customBtn);
    topRow->addWidget(fileLabel, 1);
    topRow->addWidget(saveBtn);   // top right, same height as Open
    vbox->addLayout(topRow);

    // ---- graph ----
    auto *graph = new GraphWidget;
    vbox->addWidget(graph, 1);

    // ---- smoothing row ----
    auto *smRow    = new QHBoxLayout;
    auto *smLabel  = new QLabel("FFT smoothing (cutoff period):");
    auto *slider   = new QSlider(Qt::Horizontal);
    slider->setRange(0, 5000); // cutoff period in meters (0 = off)
    slider->setSingleStep(10);
    slider->setPageStep(100);
    slider->setValue(250);
    auto *winLabel = new QLabel("250 m");
    winLabel->setMinimumWidth(100);
    winLabel->setAlignment(Qt::AlignRight | Qt::AlignVCenter);
    winLabel->setStyleSheet("color:#8ab4ff; font-weight:600;");

    smRow->addWidget(smLabel);
    smRow->addWidget(slider, 1);
    smRow->addWidget(winLabel);

    // ---- grade limits (bottom right) ----
    smRow->addSpacing(24);
    auto *limLabel = new QLabel("Grade limits:");
    limLabel->setStyleSheet("color:#9aa3b2;");
    auto *minSpin = new QDoubleSpinBox;
    minSpin->setRange(-99.0, 99.0);
    minSpin->setDecimals(1);
    minSpin->setSingleStep(0.1);
    minSpin->setValue(-20.0);
    minSpin->setSuffix(" %");
    minSpin->setAccessibleName("min grade");
    auto *maxSpin = new QDoubleSpinBox;
    maxSpin->setRange(-99.0, 99.0);
    maxSpin->setDecimals(1);
    maxSpin->setSingleStep(0.1);
    maxSpin->setValue(20.0);
    maxSpin->setSuffix(" %");
    maxSpin->setAccessibleName("max grade");
    smRow->addWidget(limLabel);
    smRow->addWidget(minSpin);
    smRow->addWidget(maxSpin);
    vbox->addLayout(smRow);

    // ---- stats row ----
    auto *statsLabel = new QLabel("—");
    statsLabel->setStyleSheet("color:#b8c0cf; font-weight:600; font-size:14px;");
    vbox->addWidget(statsLabel);

    win.setCentralWidget(central);

    Track trackStore;
    QString currentFile;

    auto refreshStats = [&]() {
        const Track &t = graph->currentTrack();
        if (t.empty() || t.size() < 2) {
            statsLabel->setText("—");
            return;
        }
        double dist = t.back().dist;
        double up = 0, down = 0, gmax = 0, gmin = 0;
        for (std::size_t i = 1; i < t.size(); ++i) {
            const double dE = t[i].ele - t[i - 1].ele;
            if (dE > 0) up += dE; else down -= dE;
            gmax = std::max(gmax, t[i].grad);
            gmin = std::min(gmin, t[i].grad);
        }
        const double gavg = dist > 0 ? (t.back().ele - t.front().ele) / dist * 100.0 : 0.0;
        statsLabel->setText(
            QString("Distance %1 km    climb %2 m    drop %3 m    |    %4 pts    gradient %5% … %6%    avg %7%")
                .arg(dist / 1000.0, 0, 'f', 1)
                .arg(up, 0, 'f', 0)
                .arg(down, 0, 'f', 0)
                .arg(int(t.size()))
                .arg(gmin, 0, 'f', 1)
                .arg(gmax, 0, 'f', 1)
                .arg(gavg, 0, 'f', 1));
    };

    auto loadFile = [&](const QString &fn) {
        std::string err;
        Track t;
        if (!loadGpx(fn.toStdString(), t, &err)) {
            QMessageBox::critical(&win, "GPX error",
                                  QString("Could not load:\n%1").arg(err.c_str()));
            return;
        }
        trackStore = std::move(t);
        currentFile = fn;
        fileLabel->setText(fn);

        graph->setTrack(trackStore);
        graph->setWindowMeters(slider->value());
        refreshStats();
    };

    auto doLoad = [&]() {
        QString fn = QFileDialog::getOpenFileName(
            &win, "Select GPX file", QString(),
            "GPX files (*.gpx);;All files (*)");
        if (fn.isEmpty()) return;
        loadFile(fn);
    };

    auto doSave = [&]() {
        if (!graph->hasTrack()) {
            QMessageBox::information(&win, "Save", "No track loaded to save.");
            return;
        }
        QString start = currentFile.isEmpty() ? QString()
                       : QFileInfo(currentFile).absolutePath() + "/" +
                         QFileInfo(currentFile).completeBaseName() + "-saved.gpx";
        QString fn = QFileDialog::getSaveFileName(
            &win, "Save GPX as", start, "GPX files (*.gpx)");
        if (fn.isEmpty()) return;
        std::string err;
        if (!saveGpx(fn.toStdString(), graph->trackForSave(), &err)) {
            QMessageBox::critical(&win, "GPX error",
                                  QString("Could not save:\n%1").arg(err.c_str()));
            return;
        }
        currentFile = fn;
        fileLabel->setText(fn + "  (saved)");
        refreshStats();
    };

    QObject::connect(openBtn, &QPushButton::clicked, doLoad);
    QObject::connect(customBtn, &QPushButton::toggled, graph, &GraphWidget::setEditEnabled);
    QObject::connect(saveBtn, &QPushButton::clicked, doSave);
    QObject::connect(graph, &GraphWidget::edited, &win, refreshStats);

    if (argc > 1) loadFile(QString::fromLocal8Bit(argv[1]));

    if (qEnvironmentVariableIsSet("GPXSMOOTH_SMOKE"))
        QTimer::singleShot(500, &app, &QCoreApplication::quit);

    // valueChanged fires continuously during a drag, and setWindowMeters
    // calls update() right away, so the graph re-renders instantly
    // under the cursor.
    QObject::connect(slider, &QSlider::valueChanged, [&](int v) {
        graph->setWindowMeters(v);
        winLabel->setText(v == 0 ? QString("off") : QString("%1 m").arg(v));
    });

    // Changing either limit immediately re-clamps the curve, redistributing
    // the excess across the neighbouring points for a smooth result.
    auto updateLimits = [=]() {
        if (minSpin->value() >= maxSpin->value()) return;
        graph->setGradeLimits(minSpin->value(), maxSpin->value());
        refreshStats();
    };
    QObject::connect(minSpin, qOverload<double>(&QDoubleSpinBox::valueChanged), updateLimits);
    QObject::connect(maxSpin, qOverload<double>(&QDoubleSpinBox::valueChanged), updateLimits);

    graph->setGradeLimits(minSpin->value(), maxSpin->value());

    win.show();
    return app.exec();
}
