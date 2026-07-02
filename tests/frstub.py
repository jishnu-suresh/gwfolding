"""Pure-Python stand-in for the ``Fr`` frame library used only in tests.

The real ``Fr`` module (frgetvect/frputvect) comes from the IGWN/conda
environment and is NOT a pip dependency. To exercise everything *except* the
frame I/O, this stub synthesises SID-frame channel data in memory so the folding
math in gwfolding runs end-to-end without any real .gwf files.

Nothing here is part of the shipped package.
"""

import numpy as np

# --- synthetic SID band definition (kept small for fast tests) -------------
FLOW = 20.0
FHIGH = 30.0
DELTAF = 0.25
NFREQ = int((FHIGH - FLOW) / DELTAF + 1e-9) + 1  # 41 bins

SEGMENT_DURATION = 192.0

# Scalar SID parameters. winFactor = w1w2bar^2 / w1w2squaredbar = 1.0,
# winRatio = 0.5 * w1w2ovlsquaredbar / w1w2squaredbar = 0.5.
_SCALARS = {
    "segmentDuration": SEGMENT_DURATION,
    "flow": FLOW,
    "fhigh": FHIGH,
    "deltaF": DELTAF,
    "w1w2bar": 1.0,
    "w1w2squaredbar": 1.0,
    "w1w2ovlsquaredbar": 1.0,
    "bias": 1.0,
}

# Records of what would have been written to frames, so tests can assert on it.
written_frames = []


def reset():
    written_frames.clear()


def _synthetic_vector(channel, gps_start):
    """Deterministic per-(channel, gps) synthetic frequency-series data."""
    seed = (abs(hash((channel, int(gps_start)))) % (2**32))
    rng = np.random.default_rng(seed)
    if channel.endswith(("AdjacentPSD", "LocalPSD")):
        # strictly positive PSD-like values
        return 1e-40 * (1.0 + 0.1 * rng.random(NFREQ))
    if channel.endswith(("ReCSD", "ImCSD")):
        return 1e-42 * (rng.random(NFREQ) - 0.5)
    # unknown vector channel -> zeros of the right length
    return np.zeros(NFREQ)


def frgetvect(filename, channel, start=-1, span=-1, verbose=False):
    """Mimic pylal-style frgetvect return: element [0] is the data array."""
    suffix = channel.split(":")[-1]
    if suffix in _SCALARS:
        data = np.array([_SCALARS[suffix]], dtype=float)
    else:
        data = np.asarray(_synthetic_vector(channel, start), dtype=float)
    # (data, gpsStart, x0, dx, xunit, yunit) — only [0] (and [0][0]) are used.
    return (data, float(start), 0.0, DELTAF, [], [])


def frputvect(filename, channellist):
    """Record the write instead of touching disk."""
    written_frames.append(
        {"filename": filename,
         "channels": [d.get("name") for d in channellist]}
    )
    return 0
