# GPX Gradient Smoother

A Qt6 desktop application that loads a GPX track, smooths its elevation
profile with a frequency-domain FFT low-pass filter, and optionally exports a
densified GPX file with one sample every metre.

## What it does

* **Load** a GPX file (`<trkpt>` with `<ele>` elevation data).
* **Display** the raw and smoothed gradient profile (percent vs distance).
* **Smooth** the elevation profile using a 4th-order Butterworth low-pass
  filter applied in the FFT domain with reflect padding and linear detrending.
* **Export** a cleaned/densified GPX file. When smoothing is enabled, the
  output contains one point every **1 metre** along the track, with latitude
  and longitude linearly interpolated and elevation taken directly from the
  FFT-smoothed profile (not re-interpolated from the original sparse samples,
  which keeps the result smooth and realistic).
* **Interactively clamp** the gradient between user-defined min/max limits,
  redistributing the excess slope to neighbouring points for a continuous
  curve.
* **Manual editing** (CustomReg mode): click & drag on the curve to modify
  the gradient at a point; the change ripples to neighbours with exponential
  decay.

## Features

* FFT low-pass filtering with adjustable cutoff period (default 250 m).
* Reflect padding to avoid border artifacts.
* Linear trend removal before filtering.
* One-point-per-metre GPX export derived from the smoothed FFT profile.
* Grade limits with automatic redistribution.
* Real-time preview of the smoothed gradient and elevation profile.
* Statistics bar: distance, climb, drop, point count, gradient min/max/avg.

## Requirements

* Qt 6 (Widgets module)
* CMake 3.16+
* A C++20 compiler

## Build

```bash
mkdir build && cd build
cmake ..
cmake --build .
```

The executable is produced at `build/gpxsmooth`.

## Usage

### Graphical interface

1. Click **Open GPX…** and select a `.gpx` file containing `<trkpt><ele>` points.
2. The gradient profile is drawn automatically.
3. Use the **FFT smoothing (cutoff period)** slider to set the low-pass
   cutoff in metres (0 = off). Default is 250 m.
4. Adjust **Grade limits** (min/max %) to clamp the curve; excess slope is
   redistributed to neighbouring points.
5. Enable **CustomReg** to manually grab & drag points on the curve.
6. Click **Save** to write a GPX file. When smoothing is active, the output
   is densified to one point per metre with FFT-smoothed elevations.

### Command line

If a file path is passed as an argument, it is loaded automatically:

```bash
./gpxsmooth path/to/track.gpx
```

The `GPXSMOOTH_SMOKE` environment variable can be set to quit automatically
after 500 ms (useful for smoke testing).

## Technical details

The C++ pipeline (`src/gpxtrack.cpp`) is a port of the Python reference in
`smooth.py`:

1. **Resampling** – the track is resampled to a uniform grid with `stepMeters`
   spacing (1 m by default) using linear interpolation (`np.interp` equivalent)
   of lat/lon/ele.
2. **Reflect padding** – the array is mirrored at both ends to avoid FFT
   border artifacts.
3. **Detrending** – a degree-1 least-squares trend is removed before filtering.
4. **FFT low-pass** – a 4th-order Butterworth gain
   `g(f) = 1 / (1 + (|f|/fc)^4)` is applied `passes` times (default 3), where
   `fc = 1 / cutoffPeriodM`.
5. **Inverse FFT** – the trend is re-added and the padding is trimmed.
6. **Export** – when saving with smoothing enabled, the track is written
   directly on the uniform 1 m grid (`smoothTrackDensified`), so the output
   GPX contains one point per metre with smooth, realistic elevations.

## Project structure

```
gpxsmooth/
├── CMakeLists.txt          # CMake build configuration (Qt6, C++20)
├── README.md               # This file
├── smooth.py               # Reference Python implementation
├── test.gpx                # Sample GPX track (500 points, ~7 km)
├── src/
│   ├── main.cpp            # Qt application window and UI wiring
│   ├── gpxtrack.h          # Track structures and public API
│   ├── gpxtrack.cpp        # GPX load/save, FFT smoothing, densification
│   ├── graphwidget.h       # Graph widget header
│   ├── graphwidget.cpp     # Gradient plot, editing, grade limits
│   └── simplefft.h         # Simple FFT implementation
└── build/                  # Build directory (created by CMake)
```

## Notes

* Latitude and longitude are never modified by smoothing; only `ele` and `grad`
  are updated.
* The number of points in the exported GPX depends on the total track length
  (one point per metre). The first and last points coincide with the original
  track endpoints.
* If smoothing is disabled (cutoff = 0), the export uses plain linear
  densification of the current track without FFT filtering.
