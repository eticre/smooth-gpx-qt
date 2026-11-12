#pragma once
#include <QWidget>
#include "gpxtrack.h"

class GraphWidget : public QWidget {
    Q_OBJECT
public:
    explicit GraphWidget(QWidget *parent = nullptr);

    void setTrack(const Track &track);
    void setWindowMeters(double windowMeters);
    void setGradeLimits(double minGrad, double maxGrad);

    // manual grab-and-drag editing of the curve (off by default)
    void setEditEnabled(bool on);
    bool editEnabled() const { return m_editEnabled; }

    bool hasTrack() const { return !m_raw.empty(); }
    const Track &currentTrack() const;   // the track currently displayed (edits applied)

    // Track da esportare: versione densificata a passo 1 m (un punto per metro)
    Track trackForSave() const;

signals:
    void edited();   // the user grabbed the curve and raised/lowered a slope

protected:
    void paintEvent(QPaintEvent *) override;
    void mousePressEvent(QMouseEvent *) override;
    void mouseMoveEvent(QMouseEvent *) override;
    void mouseReleaseEvent(QMouseEvent *) override;
    void leaveEvent(QEvent *) override;

private:
    void recomputeSmoothing();
    void applyGradeLimits(Track &t) const;   // clamp grad + redistribute to neighbours
    void ensureEditBase();      // seed m_edit from the currently displayed curve
    void updateElevProfile();   // recompute cumulative elevation + range from m_edit
    void applyEdit(int idx, double targetGrad);
    int  pickIndex(const QPoint &pos) const;
    int  distToX(double distMeters) const;
    int  gradToY(double grad) const;
    double yToGrad(int y) const;
    int  elevToY(double elev) const;
    void drawCurve(QPainter &p, const Track &t, const QColor &color,
                   int penWidth, bool dashed, bool elevation = false) const;
    static int axisStep(double range, int targetLines);

    Track  m_raw;
    Track  m_smoothed;
    Track  m_edit;        // user-editable working copy
    bool   m_hasEdit = false;
    bool   m_dragging = false;
    int    m_grabIdx = -1;
    double m_baseEle = 0.0;
    double m_windowMeters = 200.0;
    double m_gradMin = -100.0, m_gradMax = 100.0;
    bool   m_limitsActive = false;
    bool   m_axisFromLimits = false;
    bool   m_editEnabled = false;
    double m_maxDist = 0.0;
    double m_yMin = 0.0, m_yMax = 10.0;
    double m_elevMin = 0.0, m_elevMax = 1.0;
    bool   m_hover = false;
    QPoint m_hoverPos;

    // movable legend
    QPoint m_legendPos;      // custom top-left; invalid -> default (top right)
    bool   m_legendDrag = false;
    QPoint m_legendGrabOff;
    QRect  m_legendRect;     // last painted rect, used for hit-testing
    int    m_plotLeft   = 0;
    int    m_plotRight  = 0;
    int    m_plotTop    = 0;
    int    m_plotBottom = 0;
};
