#include "gpxtrack.h"
#include <QFile>
#include <QTextStream>
#include <QXmlStreamReader>
#include <QtGlobal>
#include <cmath>
#include <algorithm>
#include <complex>

#include "simplefft.h"

namespace {

constexpr double kEarthR = 6371000.0; // meters

double haversine(double lat1, double lon1, double lat2, double lon2)
{
    auto d = [](double a) { return a * M_PI / 180.0; };
    const double dLat = d(lat2 - lat1);
    const double dLon = d(lon2 - lon1);
    const double a = std::sin(dLat / 2) * std::sin(dLat / 2)
                   + std::cos(d(lat1)) * std::cos(d(lat2))
                     * std::sin(dLon / 2) * std::sin(dLon / 2);
    return 2.0 * kEarthR * std::asin(std::sqrt(a));
}

void finalize(Track &track)
{
    if (track.size() >= 2) {
        for (std::size_t i = 1; i < track.size(); ++i) {
            const double d = haversine(track[i - 1].lat, track[i - 1].lon,
                                       track[i].lat, track[i].lon);
            track[i].dist = track[i - 1].dist + d;
            if (d > 1e-6)
                track[i].grad = (track[i].ele - track[i - 1].ele) / d * 100.0;
        }
    }
}

} // namespace

bool loadGpx(const std::string &path, Track &track, std::string *error)
{
    QFile f(QString::fromStdString(path));
    if (!f.open(QIODevice::ReadOnly | QIODevice::Text)) {
        if (error) *error = "cannot open file";
        return false;
    }

    QXmlStreamReader in(&f);
    bool sawEleTag = false;
    bool inTrkpt = false;
    GpxPoint current;

    while (!in.atEnd()) {
        in.readNext();
        if (in.isStartElement()) {
            const QString name = in.name().toString();
            if (name == QLatin1String("trkpt") || name == QLatin1String("rtept")) {
                inTrkpt = true;
                current = GpxPoint();
                current.lat = in.attributes().value("lat").toString().toDouble();
                current.lon = in.attributes().value("lon").toString().toDouble();
            } else if (inTrkpt && name == QLatin1String("ele")) {
                current.ele = in.readElementText().toDouble();
                sawEleTag = true;
                inTrkpt = false;
                track.push_back(current);
                in.readNextStartElement(); // position back after <ele>
            }
        }
    }

    if (in.hasError()) {
        if (error) *error = in.errorString().toStdString();
        return false;
    }
    if (track.empty()) {
        if (error) *error = "no <trkpt> elements found (not a track file?)";
        return false;
    }
    if (!sawEleTag) {
        if (error) *error = "GPX contains no <ele> data";
        return false;
    }

    finalize(track);
    return true;
}

bool saveGpx(const std::string &path, const Track &track, std::string *error)
{
    QFile f(QString::fromStdString(path));
    if (!f.open(QIODevice::WriteOnly | QIODevice::Text)) {
        if (error) *error = "cannot open file for writing";
        return false;
    }

    QTextStream out(&f);
    out << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
        << "<gpx version=\"1.1\" creator=\"gpxsmooth\"\n"
        << "     xmlns=\"http://www.topografix.com/GPX/1/1\"\u003e\n"
        << "  <trk\u003e\n"
        << "    <trkseg\u003e\n";
    for (const auto &pt : track) {
        out << "      <trkpt lat=\"" << QString::number(pt.lat, 'g', 16)
            << "\" lon=\"" << QString::number(pt.lon, 'g', 16) << "\"\u003e\n"
            << "        <ele\u003e" << QString::number(pt.ele, 'f', 3) << "</ele\u003e\n"
            << "      </trkpt\u003e\n";
    }
    out << "    </trkseg\u003e\n"
        << "  </trk\u003e\n"
        << "</gpx\u003e\n";

    f.flush();
    f.close();
    return true;
}


namespace {

// ---------------------------------------------------------------------------
// C++ port of the filtering pipeline in smooth.py:
//   1) resample the elevation profile to a uniform step (np.arange + interp)
//   2) fft_lowpass_padded:
//        - reflect padding (avoids border artifacts)
//        - linear detrend (polyfit deg 1)
//        - FFT, low-pass gain g(f) = 1/(1+(|f|/fc)^4) raised to `passes`
//        - inverse FFT, re-add trend, trim padding
//   3) map the cleaned profile back onto the original sample distances
// ---------------------------------------------------------------------------

// np.interp (linear, increasing x; clamps outside the range)
std::vector<double> interp(const std::vector<double> &x,
                           const std::vector<double> &fp,
                           const std::vector<double> &xv)
{
    std::vector<double> out(xv.size());
    if (fp.size() == 1) {
        std::fill(out.begin(), out.end(), fp.front());
        return out;
    }
    std::size_t j = 1;
    for (std::size_t i = 0; i < xv.size(); ++i) {
        const double xi = xv[i];
        if (xi <= x.front()) { out[i] = fp.front(); j = 1; continue; }
        if (xi >= x.back())  { out[i] = fp.back();  j = x.size() - 1; continue; }
        while (j + 1 < x.size() && x[j] < xi) ++j;
        const double dx = x[j] - x[j - 1];
        out[i] = (dx < 1e-9) ? fp[j - 1]
               : fp[j - 1] + (fp[j] - fp[j - 1]) * (xi - x[j - 1]) / dx;
    }
    return out;
}

// np.pad(mode='reflect'): pads `pad` samples on each side, mirroring the
// array without repeating the end point (a b c d  ->  c b | a b c d | c b).
std::vector<double> reflectPad(const std::vector<double> &x, int pad)
{
    const long long n = static_cast<long long>(x.size());
    if (n < 2) return x;
    const long long period = 2 * n - 2;
    // Fold any virtual index onto [0, n-1] (triangular wave with turning
    // points at the ends: … x2 x1 x0 | x1 x2 x3 | x2 x1 …)
    auto at = [&](long long i) -> double {
        long long k = i % period;
        if (k < 0) k += period;
        if (k > n - 1) k = period - k;
        return x[static_cast<std::size_t>(k)];
    };
    std::vector<double> out;
    out.reserve(x.size() + 2 * static_cast<std::size_t>(pad));
    for (long long r = pad; r >= 1; --r) out.push_back(at(-r));
    for (double v : x) out.push_back(v);
    for (long long r = 1; r <= pad; ++r) out.push_back(at(n - 1 + r));
    return out;
}

// Port of fft_lowpass_padded() from smooth.py.
std::vector<double> fftLowpassPadded(const std::vector<double> &ele,
                                     double stepM,
                                     double cutoffPeriodM,
                                     int passes)
{
    const std::size_t n = ele.size();
    if (n < 2 || cutoffPeriodM <= 0.0 || passes <= 0)
        return ele;

    const int pad = static_cast<int>(std::max<std::size_t>(10, n / 5)); // int(n * 0.2)
    std::vector<double> padded = reflectPad(ele, pad);
    const std::size_t np = padded.size();

    // Linear trend (degree-1 least squares = polyfit deg 1)
    double mt = 0.0, me = 0.0;
    for (std::size_t i = 0; i < np; ++i) mt += static_cast<double>(i);
    mt /= np;
    for (double v : padded) me += v;
    me /= np;
    double cov = 0.0, var = 0.0;
    for (std::size_t i = 0; i < np; ++i) {
        const double dt = static_cast<double>(i) - mt;
        cov += dt * (padded[i] - me);
        var += dt * dt;
    }
    const double slope = cov / var;
    const double intercept = me - slope * mt;

    std::vector<std::complex<double>> spectrum(np);
    for (std::size_t i = 0; i < np; ++i)
        spectrum[i] = padded[i] - (intercept + slope * static_cast<double>(i));
    simplefft::fft(spectrum, false);

    // Butterworth-4 low-pass gain, repeated `passes` times
    const double fc = 1.0 / cutoffPeriodM;
    for (std::size_t i = 0; i < np; ++i) {
        const std::size_t k = (i <= np / 2) ? i : np - i;   // |freq| index
        const double f = static_cast<double>(k) / (np * stepM);
        double g = 1.0 / (1.0 + std::pow(f / fc, 4.0));
        spectrum[i] *= std::pow(g, static_cast<double>(passes));
    }

    simplefft::fft(spectrum, true);

    std::vector<double> out(n);
    for (std::size_t i = 0; i < n; ++i)
        out[i] = spectrum[pad + i].real()
               + (intercept + slope * static_cast<double>(pad + i));
    return out;
}

} // namespace

Track smoothTrack(const Track &src, double cutoffPeriodM,
                  int passes, double stepMeters)
{
    Track work = src;
    const std::size_t n = src.size();
    if (n < 2 || cutoffPeriodM <= 0.0 || src.back().dist <= 0.0)
        return work;

    // 1) resample to the uniform `stepMeters` grid (like resample_track)
    const double maxDist = src.back().dist;
    std::vector<double> distO(n), eleO(n), distU;
    for (std::size_t i = 0; i < n; ++i) {
        distO[i] = src[i].dist;
        eleO[i]  = src[i].ele;
    }
    for (double d = 0.0; d < maxDist; d += stepMeters)
        distU.push_back(d);
    if (distU.empty() || distU.back() < maxDist)
        distU.push_back(maxDist);

    // 2) FFT low-pass on the uniformly resampled elevations
    std::vector<double> eleU = interp(distO, eleO, distU);
    std::vector<double> cleaned =
        fftLowpassPadded(eleU, stepMeters, cutoffPeriodM, passes);

    // 3) map the cleaned profile back onto the original points
    std::vector<double> eleNew = interp(distU, cleaned, distO);
    work.front().ele = eleNew.front();
    for (std::size_t i = 1; i < n; ++i) {
        work[i].ele = eleNew[i];
        const double d = src[i].dist - src[i - 1].dist;
        if (d > 1e-9)
            work[i].grad = (work[i].ele - work[i - 1].ele) / d * 100.0;
    }
    return work;
}

// smoothTrackDensified: returns the FFT-smoothed track directly on the
// uniform `stepMeters` grid (lat/lon interpolated, ele = cleaned profile).
Track smoothTrackDensified(const Track &src, double cutoffPeriodM,
                           int passes, double stepMeters)
{
    Track out;
    const std::size_t n = src.size();
    if (n < 2 || cutoffPeriodM <= 0.0 || src.back().dist <= 0.0)
        return out;

    std::vector<double> distO(n), latO(n), lonO(n), eleO(n);
    for (std::size_t i = 0; i < n; ++i) {
        distO[i] = src[i].dist;
        latO[i]  = src[i].lat;
        lonO[i]  = src[i].lon;
        eleO[i]  = src[i].ele;
    }

    const double maxDist = src.back().dist;
    std::vector<double> distU;
    for (double d = 0.0; d < maxDist; d += stepMeters)
        distU.push_back(d);
    if (distU.empty() || distU.back() < maxDist)
        distU.push_back(maxDist);

    // Linear interpolation of lat/lon/ele onto the uniform grid (np.interp)
    std::vector<double> latU = interp(distO, latO, distU);
    std::vector<double> lonU = interp(distO, lonO, distU);
    std::vector<double> eleU = interp(distO, eleO, distU);

    // FFT low-pass directly on the uniform grid
    std::vector<double> cleaned =
        fftLowpassPadded(eleU, stepMeters, cutoffPeriodM, passes);

    out.reserve(distU.size());
    for (std::size_t i = 0; i < distU.size(); ++i) {
        GpxPoint p;
        p.lat = latU[i];
        p.lon = lonU[i];
        p.ele = cleaned[i];
        p.dist = distU[i];
        p.grad = 0.0;
        out.push_back(p);
    }

    // Recompute per-segment gradient on the densified samples.
    if (out.size() >= 2) {
        out.front().grad = 0.0;
        for (std::size_t i = 1; i < out.size(); ++i) {
            const double d = haversine(out[i - 1].lat, out[i - 1].lon,
                                       out[i].lat, out[i].lon);
            out[i].grad = (d > 1e-6)
                ? (out[i].ele - out[i - 1].ele) / d * 100.0
                : out[i - 1].grad;
        }
    }
    return out;
}

Track resampleTrackUniform(const Track &src, double stepMeters)
{
    Track out;
    const std::size_t n = src.size();
    if (n < 2 || stepMeters <= 0.0)
        return out;

    std::vector<double> distO(n), latO(n), lonO(n), eleO(n);
    for (std::size_t i = 0; i < n; ++i) {
        distO[i] = src[i].dist;
        latO[i]  = src[i].lat;
        lonO[i]  = src[i].lon;
        eleO[i]  = src[i].ele;
    }

    const double maxDist = src.back().dist;
    std::vector<double> distU;
    for (double d = 0.0; d < maxDist; d += stepMeters)
        distU.push_back(d);
    if (distU.empty() || distU.back() < maxDist)
        distU.push_back(maxDist);

    // Linear interpolation of lat/lon/ele onto the uniform grid (np.interp)
    std::vector<double> latU = interp(distO, latO, distU);
    std::vector<double> lonU = interp(distO, lonO, distU);
    std::vector<double> eleU = interp(distO, eleO, distU);

    out.reserve(distU.size());
    for (std::size_t i = 0; i < distU.size(); ++i) {
        GpxPoint p;
        p.lat = latU[i];
        p.lon = lonU[i];
        p.ele = eleU[i];
        p.dist = distU[i];
        p.grad = 0.0;
        out.push_back(p);
    }

    // Recompute per-segment gradient on the densified samples.
    if (out.size() >= 2) {
        out.front().grad = 0.0;
        for (std::size_t i = 1; i < out.size(); ++i) {
            const double d = haversine(out[i - 1].lat, out[i - 1].lon,
                                       out[i].lat, out[i].lon);
            out[i].grad = (d > 1e-6)
                ? (out[i].ele - out[i - 1].ele) / d * 100.0
                : out[i - 1].grad;
        }
    }
    return out;
}
