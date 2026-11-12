#include "graphwidget.h"
#include <QPainter>
#include <QMouseEvent>
#include <QPainterPath>
#include <cmath>
#include <vector>
#include <algorithm>

// ---- dark palette ----
namespace {
const QColor kBgMain    (24, 27, 33);
const QColor kGrid      (45, 51, 62);
const QColor kBorder    (58, 65, 78);
const QColor kTextDim   (150, 157, 170);
const QColor kText      (176, 183, 198);
const QColor kZero      (214, 96, 96);
const QColor kRaw       (110, 118, 132);
const QColor kSmooth    (96, 165, 250);
const QColor kElev      (74, 222, 128);
const QColor kElevDim   (122, 183, 150);
const QColor kLegendBg  (38, 43, 53);
const QColor kLegendPen (62, 69, 82);
}

GraphWidget::GraphWidget(QWidget *parent) : QWidget(parent)
{
    setMinimumSize(500, 300);
    setMouseTracking(true);
}

void GraphWidget::setTrack(const Track &track)
{
    m_raw = track;
    m_edit.clear();
    m_hasEdit = false;
    m_dragging = false;
    m_grabIdx = -1;
    m_baseEle = m_raw.empty() ? 0.0 : m_raw.front().ele;
    recomputeSmoothing();
    update();
}

void GraphWidget::setEditEnabled(bool on)
{
    if (on == m_editEnabled)
        return;
    m_editEnabled = on;
    if (!on && m_hasEdit) {
        // Disabling CustomReg discards the manual edits: the displayed
        // curve goes back to the FFT-smoothed profile, so the smoothing
        // slider (and the grade limits) take full effect again.
        m_hasEdit = false;
        m_edit.clear();
        recomputeSmoothing();
        emit edited();   // keep the stats row in sync in the main window
        update();
    }
}

void GraphWidget::setWindowMeters(double windowMeters)
{
    if (qFuzzyCompare(windowMeters, m_windowMeters))
        return;
    m_windowMeters = windowMeters;
    recomputeSmoothing();
    update(); // instant repaint: no timer, no queue
}

void GraphWidget::setGradeLimits(double minGrad, double maxGrad)
{
    if (minGrad >= maxGrad)
        return;
    if (qFuzzyCompare(minGrad, m_gradMin) && qFuzzyCompare(maxGrad, m_gradMax))
        return;
    m_gradMin = minGrad;
    m_gradMax = maxGrad;
    m_limitsActive = true;
    m_axisFromLimits = true;   // y-axis now follows the grade boxes
    if (m_hasEdit)
        applyGradeLimits(m_edit);
    recomputeSmoothing();
    update();
}

void GraphWidget::applyGradeLimits(Track &t) const
{
    if (t.size() < 3) return;
    const double lo = m_gradMin, hi = m_gradMax;

    // Multi-pass projection: clamp every offending point to the allowed
    // range and push the excess onto its neighbours with the same decay
    // used for the interactive edit, so the resulting curve stays smooth,
    // uniform and continuous (no abrupt kinks at the boundary).
    for (int pass = 0; pass < 12; ++pass) {
        bool any = false;
        for (int i = 1; i + 1 < (int)t.size(); ++i) {
            double target = t[i].grad;
            if (target > hi)      target = hi;
            else if (target < lo) target = lo;
            if (std::abs(target - t[i].grad) < 1e-6) continue;
            any = true;
            const double delta = target - t[i].grad;
            const double decay = 0.55;
            for (int s = 1; s < (int)t.size(); ++s) {
                const double d = delta * std::pow(decay, double(s));
                if (std::abs(d) < 0.005 && s > 1) break;
                if (i - s >= 0)            t[i - s].grad += d;
                if (i + s < (int)t.size()) t[i + s].grad += d;
            }
            t[i].grad = target;
        }
        if (!any) break;
    }

    // Keep the elevation profile continuous with the limited slopes.
    for (std::size_t k = 1; k < t.size(); ++k)
        t[k].ele = t[k - 1].ele
                 + t[k].grad / 100.0 * (t[k].dist - t[k - 1].dist);
}

void GraphWidget::recomputeSmoothing()
{
    if (m_raw.empty()) return;
    m_smoothed = (m_windowMeters > 0.0) ? smoothTrack(m_raw, m_windowMeters) : m_raw;
    if (m_limitsActive)
        applyGradeLimits(m_smoothed);
    m_maxDist = m_raw.back().dist;
    if (m_maxDist <= 0) m_maxDist = 1.0;

    double maxAbs = 0.0;
    for (const auto &p : m_raw) maxAbs = std::max(maxAbs, std::abs(p.grad));
    m_yMax = std::max(5.0, maxAbs * 1.15);
    if (m_axisFromLimits)   // axis tracks the spinbox limits: [-max, +max]
        m_yMax = std::max(1.0, std::max(std::abs(m_gradMin), std::abs(m_gradMax)));

    double elevLo = 1e300, elevHi = -1e300;
    const Track &eleSrc = m_hasEdit ? m_edit
                       : (m_limitsActive ? m_smoothed : m_raw);
    for (const auto &p : eleSrc) {
        elevLo = std::min(elevLo, p.ele);
        elevHi = std::max(elevHi, p.ele);
    }
    m_elevMin = elevLo;
    m_elevMax = elevHi;
}

int GraphWidget::distToX(double distMeters) const
{
    if (m_maxDist <= 0) return m_plotLeft;
    return m_plotLeft
         + int(std::lround(distMeters / m_maxDist * (m_plotRight - m_plotLeft)));
}

int GraphWidget::gradToY(double grad) const
{
    const double t = (grad + m_yMax) / (2.0 * m_yMax); // 0..1
    return m_plotBottom - int(std::lround(t * (m_plotBottom - m_plotTop)));
}

int GraphWidget::elevToY(double ele) const
{
    if (m_elevMax - m_elevMin < 1e-9) return m_plotBottom;
    const double t = (ele - m_elevMin) / (m_elevMax - m_elevMin); // 0..1
    return m_plotBottom - int(std::lround(t * (m_plotBottom - m_plotTop)));
}

int GraphWidget::axisStep(double range, int targetLines)
{
    const double raw = range / targetLines;
    double pow10 = std::pow(10.0, std::floor(std::log10(raw)));
    for (double c : {1.0, 2.0, 5.0, 10.0}) {
        if (pow10 * c >= raw)
            return int(std::round(pow10 * c));
    }
    return int(std::round(pow10 * 10.0));
}

void GraphWidget::paintEvent(QPaintEvent *)
{
    QPainter p(this);
    p.setRenderHint(QPainter::Antialiasing);

    const int w = width(), h = height();
    const bool showElev = (m_elevMax - m_elevMin) > 1e-6; // needs real variation
    const int leftMargin = 68, rightMargin = showElev ? 64 : 10, topMargin = 10, botMargin = 26;
    m_plotLeft   = leftMargin;
    m_plotRight  = w - rightMargin;
    m_plotTop    = topMargin;
    m_plotBottom = h - botMargin;
    const int plotW = m_plotRight - m_plotLeft;
    const int plotH = m_plotBottom - m_plotTop;
    if (plotW < 20 || plotH < 20) return;

    p.fillRect(rect(), kBgMain);

    if (m_raw.empty()) {
        p.setPen(kTextDim);
        p.drawText(rect(), Qt::AlignCenter,
                   "Open a GPX file to see its gradient profile");
        return;
    }

    const int y0 = gradToY(0.0);

    // ---- grid + axes ----
    QFont f = font();
    f.setPointSize(std::max(7, f.pointSize() - 1));
    p.setFont(f);

    p.setPen(kTextDim);
    p.drawText(rect().adjusted(0, 24, 0, 22), Qt::AlignCenter,
               "Distance (km)");

    // vertical grid: distance, step in meters
    int step = axisStep(m_maxDist, std::max(2, plotW / 110));
    step = std::max(int(step), 1);
    for (double d = 0; d <= m_maxDist; d += step) {
        const int x = distToX(d);
        p.setPen(kGrid);
        p.drawLine(x, m_plotTop, x, m_plotBottom);
        p.setPen(kText);
        p.drawText(QRect(x - 20, m_plotBottom + 3, 40, 16),
                   Qt::AlignCenter, QString::number(d / 1000.0, 'g', 3));
    }
    // horizontal grid: percent
    int gstep = axisStep(m_yMax * 2.0, std::max(2, plotH / 45));
    for (double g = -m_yMax; g <= m_yMax + 1e-9; g += gstep) {
        if (std::abs(g) > m_yMax) break;
        const int y = gradToY(g);
        p.setPen(kGrid);
        p.drawLine(m_plotLeft, y, m_plotRight, y);
        p.setPen(kText);
        p.drawText(QRect(0, y - 11, m_plotLeft - 10, 22),
                   Qt::AlignRight | Qt::AlignVCenter,
                   QString::number(g, 'g', 3) + "%");
    }

    // right axis: elevation ticks (green, matches the elevation curve)
    if (showElev) {
        const int estep = std::max(1, axisStep(m_elevMax - m_elevMin, std::max(2, plotH / 45)));
        const double e0 = std::ceil(m_elevMin / estep) * estep;
        const int dec = estep < 1 ? 1 : 0;
        for (double e = e0; e <= m_elevMax + 1e-9; e += estep) {
            const int y = elevToY(e);
            p.setPen(kElevDim);
            p.drawText(QRect(m_plotRight + 6, y - 11, 56, 22),
                       Qt::AlignLeft | Qt::AlignVCenter,
                       QString::number(e, 'f', dec));
        }
        p.setPen(kElevDim);
        p.drawText(QRect(m_plotRight + 6, m_plotTop - 1, 56, 16),
                   Qt::AlignLeft | Qt::AlignTop, "elev (m)");
    }

    // zero line
    p.setPen(QPen(kZero, 1));
    p.drawLine(m_plotLeft, y0, m_plotRight, y0);

    // plot border
    p.setPen(kBorder);
    p.drawRect(m_plotLeft, m_plotTop, plotW, plotH);

    // ---- curves ----
    const Track &disp = currentTrack();
    // elevation first (green) so the gradient curves render on top
    if (showElev)
        drawCurve(p, disp, kElev, 2, false, true);
    if (m_windowMeters > 0.0 || m_hasEdit || m_limitsActive) {
        drawCurve(p, m_raw, kRaw, 1, true);
        drawCurve(p, disp, kSmooth, 2, false);
    } else {
        drawCurve(p, disp, kSmooth, 2, false);
    }

    // legend (top right, larger & readable)
    QFont lf = font();
    lf.setPointSize(std::max(9, lf.pointSize()));
    p.setFont(lf);
    const int nRows = ((m_windowMeters > 0.0 || m_hasEdit || m_limitsActive) ? 2 : 1) + (showElev ? 1 : 0);
    const int lw = 250;
    const int rowH = 26;
    const int lh = nRows * rowH + 12;
    const QSize legendSize(lw - 16, lh);
    QRect legend;
    if (!m_legendPos.isNull())   // user moved the legend
        legend = QRect(m_legendPos, legendSize);
    else
        legend = QRect(QPoint(m_plotRight - lw, m_plotBottom - 8 - legendSize.height()), legendSize);  // default: bottom right
    // clamp inside the widget
    legend.setX(qBound(4, legend.x(), std::max(4, rect().width() - legend.width() - 4)));
    legend.setY(qBound(4, legend.y(), std::max(4, rect().height() - legend.height() - 4)));
    m_legendRect = legend;
    p.setBrush(QBrush(QColor(38, 43, 53, 235)));
    p.setPen(kLegendPen);
    p.drawRoundedRect(legend, 5, 5);
    auto drawRow = [&](int rowIdx, const QColor &c, int wpx, Qt::PenStyle style,
                       const QString &label) {
        const int cy = legend.top() + 8 + rowIdx * rowH + rowH / 2;
        p.setPen(QPen(c, wpx, style));
        p.drawLine(legend.left() + 10, cy, legend.left() + 34, cy);
        p.setPen(kText);
        p.drawText(QRect(legend.left() + 42, cy - rowH / 2,
                         legend.width() - 50, rowH),
                   Qt::AlignLeft | Qt::AlignVCenter, label);
    };
    int row = 0;
    if (m_windowMeters > 0.0 || m_hasEdit || m_limitsActive) {
        drawRow(row++, kRaw, 1, Qt::DashLine, "original");
        drawRow(row++, kSmooth, 2, Qt::SolidLine,
                m_hasEdit ? QString("edited (grab \u0026 drag)")
                          : (m_limitsActive
                             ? QString("limited %1% to %2%")
                                   .arg(m_gradMin, 0, 'f', 1).arg(m_gradMax, 0, 'f', 1)
                             : QString("smoothed (\u00b1%1 m)").arg(int(m_windowMeters / 2.0))));
    } else {
        drawRow(row++, kSmooth, 2, Qt::SolidLine, QString("gradient (grab \u0026 drag)"));
    }
    if (showElev)
        drawRow(row++, kElev, 2, Qt::SolidLine, "elevation (m)");

    // ---- grabbed point marker ----
    if (m_dragging && m_hasEdit && m_grabIdx >= 0 && m_grabIdx < (int)m_edit.size()) {
        const GpxPoint &pt = m_edit[static_cast<std::size_t>(m_grabIdx)];
        const int px = distToX(pt.dist);
        const int py = gradToY(pt.grad);
        p.setBrush(QColor(240, 244, 255));
        p.setPen(QPen(kSmooth, 2));
        p.drawEllipse(QPoint(px, py), 6, 6);
        QFont tfont = font();
        tfont.setBold(true);
        tfont.setPointSize(std::max(10, tfont.pointSize() + 1));
        p.setFont(tfont);
        const QString val = QString("%1%  ").arg(pt.grad, 0, 'f', 1);
        const int vw = p.fontMetrics().horizontalAdvance(val) + 16;
        int vb = px + 10;
        if (vb + vw > m_plotRight) vb = px - 10 - vw;
        p.setBrush(QBrush(QColor(52, 58, 72, 245)));
        p.setPen(kLegendPen);
        p.drawRoundedRect(vb, py - 26, vw, 20, 4, 4);
        p.setPen(QColor(240, 242, 248));
        p.drawText(QRect(vb, py - 26, vw, 20), Qt::AlignCenter, val);
    }

    // ---- hover crosshair ----
    if (m_hover && !m_dragging && !m_legendDrag
               && m_hoverPos.x() >= m_plotLeft && m_hoverPos.x() <= m_plotRight
               && m_hoverPos.y() >= m_plotTop && m_hoverPos.y() <= m_plotBottom
               && !disp.empty()) {
        const double d = (m_hoverPos.x() - m_plotLeft) / double(plotW) * m_maxDist;
        auto it = std::lower_bound(
            disp.begin(), disp.end(), d,
            [](const GpxPoint &pt, double v) { return pt.dist < v; });
        if (it != disp.end()) {
            const GpxPoint &pt = *it;
            const int px = distToX(pt.dist);

            QFont tfont = font();
            tfont.setBold(true);
            tfont.setPointSize(std::max(10, tfont.pointSize() + 1));
            p.setFont(tfont);
            p.setPen(QPen(kSmooth, 1, Qt::DashLine));            p.drawLine(px, m_plotTop, px, m_plotBottom);
            p.setPen(QColor(220, 225, 235, 160));
            p.drawLine(m_plotLeft, m_hoverPos.y(), m_plotRight, m_hoverPos.y());

            const QString text = QString("%1 km    %2%    %3 m elev")
                .arg(pt.dist / 1000.0, 0, 'f', 2)
                .arg(pt.grad, 0, 'f', 1)
                .arg(pt.ele, 0, 'f', 0);
            const int tw = p.fontMetrics().horizontalAdvance(text) + 24;
            const int th = 30;
            int bx = px + 12;
            if (bx + tw > m_plotRight) bx = px - 12 - tw;
            int by = m_hoverPos.y() - 2 - th;
            if (by < m_plotTop) by = m_hoverPos.y() + 10;

            p.setBrush(QBrush(QColor(52, 58, 72, 245)));
            p.setPen(kLegendPen);
            p.drawRoundedRect(bx, by, tw, th, 5, 5);
            p.setPen(QColor(240, 242, 248));
            p.drawText(QRect(bx, by, tw, th), Qt::AlignCenter, text);
        }
    }
}

void GraphWidget::drawCurve(QPainter &painter, const Track &t, const QColor &color,
                            int penWidth, bool dashed,
                            bool elevation) const
{
    if (t.empty()) return;
    const int plotW = m_plotRight - m_plotLeft;
    painter.setClipRect(QRect(m_plotLeft, m_plotTop, plotW, m_plotBottom - m_plotTop));
    if (plotW <= 1) return;

    // 1) Decimate to one value per pixel column; the median keeps the curve
    //    faithful while resisting noisy individual samples.
    std::vector<std::vector<double>> col(static_cast<std::size_t>(plotW));
    std::vector<std::size_t> count(static_cast<std::size_t>(plotW), 0);
    for (const auto &pt : t) {
        const int cx = distToX(pt.dist) - m_plotLeft;
        if (cx < 0 || cx >= plotW) continue;
        ++count[static_cast<std::size_t>(cx)];
    }
    for (int c = 0; c < plotW; ++c)
        col[static_cast<std::size_t>(c)].reserve(count[static_cast<std::size_t>(c)]);
    for (const auto &pt : t) {
        const int cx = distToX(pt.dist) - m_plotLeft;
        if (cx < 0 || cx >= plotW) continue;
        col[static_cast<std::size_t>(cx)].push_back(elevation ? pt.ele : pt.grad);
    }

    std::vector<QPointF> pts;
    pts.reserve(1024);
    for (int c = 0; c < plotW; ++c) {
        std::vector<double> &v = col[static_cast<std::size_t>(c)];
        if (v.empty()) continue;
        std::nth_element(v.begin(), v.begin() + v.size() / 2, v.end());
        const double mid = v[v.size() / 2];
        pts.emplace_back(m_plotLeft + c, elevation ? elevToY(mid) : gradToY(mid));
    }

    const std::size_t n = pts.size();
    if (n == 0) return;

    // 2) Monotone cubic (Fritsch-Carlson) tangents: smooth curve joining
    //    the points, guaranteed not to overshoot beyond the data.
    std::vector<double> d(n ? n - 1 : 0), m(n, 0.0);
    if (n >= 2) {
        for (std::size_t i = 0; i + 1 < n; ++i) {
            const double dx = pts[i + 1].x() - pts[i].x();
            d[i] = (dx > 1e-9) ? (pts[i + 1].y() - pts[i].y()) / dx : 0.0;
        }
        m[0] = d[0];
        m[n - 1] = d[n - 2];
        for (std::size_t i = 1; i + 1 < n; ++i) {
            if (d[i - 1] * d[i] <= 0.0) m[i] = 0.0;                       // local extremum
            else m[i] = (2.0 * d[i - 1] * d[i]) / (d[i - 1] + d[i]);      // harmonic mean
        }
        for (std::size_t i = 0; i + 1 < n; ++i) {
            if (d[i] == 0.0) { m[i] = 0.0; m[i + 1] = 0.0; continue; }
            const double a = m[i] / d[i];
            const double b = m[i + 1] / d[i];
            const double s = a * a + b * b;
            if (s > 9.0) {
                const double tau = 3.0 / std::sqrt(s);
                m[i] = tau * a * d[i];
                m[i + 1] = tau * b * d[i];
            }
        }
    }

    // 3) Hermite segments -> cubic Bezier and draw.
    QPainterPath path;
    path.moveTo(pts.front());
    for (std::size_t i = 0; i + 1 < n; ++i) {
        const double dx3 = (pts[i + 1].x() - pts[i].x()) / 3.0;
        path.cubicTo(pts[i].x() + dx3, pts[i].y() + m[i] * dx3,
                     pts[i + 1].x() - dx3, pts[i + 1].y() - m[i + 1] * dx3,
                     pts[i + 1].x(), pts[i + 1].y());
    }

    QPen pen(color, penWidth);
    if (dashed) pen.setStyle(Qt::DashLine);
    painter.setPen(pen);
    painter.drawPath(path);
}

const Track &GraphWidget::currentTrack() const
{
    if (m_hasEdit) return m_edit;
    if (m_windowMeters > 0.0 || m_limitsActive) return m_smoothed;
    return m_raw;
}

Track GraphWidget::trackForSave() const
{
    Track out;
    if (m_windowMeters > 0.0) {
        // Esporta il profilo già filtrato con FFT direttamente sulla griglia
        // uniforme a 1 m: lat/lon interpolati, ele = profilo pulito.
        // Questo è più liscio e realistico di un semplice risamplamento
        // dei punti originali.
        out = smoothTrackDensified(m_raw, m_windowMeters, 3, 1.0);
        if (m_limitsActive)
            applyGradeLimits(out);
    } else {
        // Nessun smoothing attivo: densificazione lineare del track corrente.
        out = resampleTrackUniform(currentTrack(), 1.0);
    }
    return out;
}

double GraphWidget::yToGrad(int y) const
{
    y = qBound(m_plotTop, y, m_plotBottom);
    const double t = (m_plotBottom - y) / double(std::max(1, m_plotBottom - m_plotTop));
    return t * 2.0 * m_yMax - m_yMax;
}

int GraphWidget::pickIndex(const QPoint &pos) const
{
    const Track &t = currentTrack();
    if (t.empty()) return -1;
    const int plotW = std::max(1, m_plotRight - m_plotLeft);
    const double d = (pos.x() - m_plotLeft) / double(plotW) * m_maxDist;
    auto it = std::lower_bound(
        t.begin(), t.end(), d,
        [](const GpxPoint &pt, double v) { return pt.dist < v; });
    if (it == t.end()) return int(t.size()) - 1;
    if (it != t.begin() &&
        std::abs((it - 1)->dist - d) < std::abs(it->dist - d))
        return int((it - 1) - t.begin());
    return int(it - t.begin());
}

void GraphWidget::ensureEditBase()
{
    if (m_hasEdit || m_raw.empty()) return;
    // Edits start from whatever the blue curve currently shows
    // (the smoothed profile, or the raw one).
    m_edit = (m_windowMeters > 0.0) ? smoothTrack(m_raw, m_windowMeters) : m_raw;
    if (m_limitsActive)
        applyGradeLimits(m_edit);
    m_hasEdit = true;
    updateElevProfile();
}

void GraphWidget::updateElevProfile()
{
    if (m_edit.empty()) return;
    // Recompute the elevation profile so it stays continuous with the
    // edited slopes: elev[k] = elev[k-1] + grad[k]/100 * horizontal step.
    m_edit[0].ele = m_baseEle;
    for (std::size_t k = 1; k < m_edit.size(); ++k)
        m_edit[k].ele = m_edit[k - 1].ele
            + m_edit[k].grad / 100.0 * (m_edit[k].dist - m_edit[k - 1].dist);
    double lo = m_edit.front().ele, hi = lo;
    for (const auto &p : m_edit) {
        lo = std::min(lo, p.ele);
        hi = std::max(hi, p.ele);
    }
    m_elevMin = lo;
    m_elevMax = hi;
}

void GraphWidget::applyEdit(int idx, double targetGrad)
{
    if (idx <= 0 || idx >= (int)m_edit.size()) return;
    const double delta = targetGrad - m_edit[idx].grad;
    if (std::abs(delta) < 1e-9) return;
    m_grabIdx = idx;

    // Ripple the change to the neighbouring points, decaying at each step,
    // until the perturbation vanishes: the grabbed point moves fully to the
    // requested slope, adjacent points follow in decreasing measure, and
    // the elevation profile stays continuous everywhere.
    const double decay = 0.55;
    for (int j = 1; j < (int)m_edit.size(); ++j) {
        const double d = delta * std::pow(decay, static_cast<double>(j));
        if (std::abs(d) < 0.01) break;
        if (idx - j >= 0) m_edit[idx - j].grad += d;
        if (idx + j < (int)m_edit.size()) m_edit[idx + j].grad += d;
    }
    m_edit[idx].grad = targetGrad;

    if (m_limitsActive)
        applyGradeLimits(m_edit);   // clamp back into [min, max] + keep elev continuous
    else
        updateElevProfile();
    emit edited();
    update();
}

void GraphWidget::mousePressEvent(QMouseEvent *e)
{
    if (e->button() != Qt::LeftButton || m_raw.empty()) return;
    const auto pos = e->position().toPoint();

    // grabbing the legend: it can be dragged anywhere over the widget
    if (m_legendRect.contains(pos)) {
        m_legendDrag = true;
        m_legendPos = m_legendRect.topLeft();
        m_legendGrabOff = pos - m_legendRect.topLeft();
        setCursor(Qt::PointingHandCursor);
        return;
    }

    /** Curve editing only while CustomReg mode is on */
    if (!m_editEnabled)
        return;

    if (pos.x() < m_plotLeft || pos.x() > m_plotRight)
        return;
    const int idx = pickIndex(pos);
    ensureEditBase();
    m_dragging = true;
    setCursor(Qt::PointingHandCursor);
    applyEdit(idx, yToGrad(pos.y()));
}

void GraphWidget::mouseReleaseEvent(QMouseEvent *e)
{
    if (e->button() == Qt::LeftButton) {
        m_dragging = false;
        m_legendDrag = false;
        unsetCursor();
    }
}

void GraphWidget::mouseMoveEvent(QMouseEvent *e)
{
    m_hover = true;
    m_hoverPos = e->position().toPoint();
    if (m_legendDrag) {
        m_legendPos = m_hoverPos - m_legendGrabOff;
        update();
        return;
    }
    if (m_dragging && m_hoverPos.x() >= m_plotLeft && m_hoverPos.x() <= m_plotRight)
        applyEdit(pickIndex(m_hoverPos), yToGrad(m_hoverPos.y()));
    update();
}

void GraphWidget::leaveEvent(QEvent *)
{
    m_hover = false;
    update();
}
