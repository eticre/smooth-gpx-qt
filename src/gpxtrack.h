#pragma once
#include <vector>
#include <string>

struct GpxPoint {
    double lat = 0.0;
    double lon = 0.0;
    double ele = 0.0;      // meters (default 0 if GPX has no ele)
    double dist = 0.0;     // cumulative distance from track start, meters
    double grad = 0.0;     // grade in percent (dElev / horizontal d * 100)
};

typedef std::vector<GpxPoint> Track;

// Parses a GPX file, extracts the first <trk> (or <rte>), computes
// cumulative distance and per-point gradient.
// Returns false and fills *error if the file cannot be parsed as a track.
bool loadGpx(const std::string &path, Track &track, std::string *error);

// Writes the track back to a GPX file (<gpx>/<trk>/<trkseg>/<trkpt>) with
// lat/lon/ele. Returns false and fills *error on I/O failure.
bool saveGpx(const std::string &path, const Track &track, std::string *error);

// Smoothing of the elevation profile (C++ port of the pipeline in smooth.py):
// the track is resampled to a uniform `stepMeters` spacing, reflect-padded
// to avoid border artifacts, detrended (linear), low-pass filtered in the
// frequency domain with a 4th-order Butterworth gain whose cutoff period is
// `cutoffPeriodM` meters, the filter being applied `passes` times in a row;
// the result is mapped back onto the original sample positions.
// Lat/lon are left unchanged; ele and grad get updated.
// cutoffPeriodM <= 0 returns a copy unchanged.
Track smoothTrack(const Track &src, double cutoffPeriodM,
                  int passes = 3, double stepMeters = 1.0);

// Resample a track to a uniform spacing of `stepMeters` meters along the
// cumulative distance, linearly interpolating lat/lon/ele. Produces one
// point every `stepMeters` metres, useful for exporting a densified GPX
// file. dist and grad are recomputed on the new samples.
Track resampleTrackUniform(const Track &src, double stepMeters = 1.0);

// Like smoothTrack but, instead of mapping the cleaned profile back onto the
// original sample positions, returns the track directly on the uniform
// `stepMeters` grid with the FFT-filtered elevations. This yields a smooth
// and realistic densified GPX (one point every metre) derived from the FFT
// pipeline rather than re-interpolating the original samples.
Track smoothTrackDensified(const Track &src, double cutoffPeriodM,
                           int passes = 3, double stepMeters = 1.0);
