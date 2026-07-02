"""Unit tests for the dependency-free helpers in gwfolding.foldUtils."""

import numpy as np

from gwfolding import foldUtils


def test_gmst_in_range():
    # Greenwich Mean Sidereal Time should wrap into [0, 24) hours.
    for gps in [630720013, 1000000000, 1396809233, 1396809233 + 43082]:
        st = foldUtils.GPStoGreenwichMeanSiderealTime(gps)
        assert 0.0 <= st < 24.0


def test_gmst_unwrapped_monotonic():
    # Without the modulo, sidereal time increases with GPS time.
    a = foldUtils.GPStoGreenwichMeanSiderealTime(1396809233, usemod=False)
    b = foldUtils.GPStoGreenwichMeanSiderealTime(1396809233 + 3600, usemod=False)
    assert b > a


def test_smooth_length_preserved():
    a = np.arange(50, dtype=float)
    out = foldUtils.smooth(a, 5)
    assert out.size == a.size


def test_smooth_constant_signal():
    a = np.full(30, 7.0)
    out = foldUtils.smooth(a, 5)
    assert np.allclose(out, 7.0)
