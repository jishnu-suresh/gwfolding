"""End-to-end test of the folding pipeline using the synthetic ``Fr`` stub.

Exercises the real numpy folding math in foldSID -> loadSID -> writeFSID for both
folding branches (identicalNeighbors on/off). Frame I/O is faked, so no real
.gwf files or the IGWN ``Fr`` library are required.

The optional sigma-cut / weight-file code paths are intentionally NOT enabled
here (the sample parameter file does not use them).
"""

import numpy as np
import pytest

import frstub
from gwfolding import foldSID


N_SEGMENTS = 20            # number of synthetic SID segments
SEG_DURATION = 192         # must match frstub.SEGMENT_DURATION
FOLDED_SEG = SEG_DURATION / 2
GPS0 = 1396809233


def _write_inputs(tmp_path, identical_neighbors):
    # Metadata: one line per segment "GPSStart<TAB>SIDfile\n"
    meta = tmp_path / "SID_MetaData.txt"
    with open(meta, "w") as f:
        for i in range(N_SEGMENTS):
            gps = GPS0 + i * int(FOLDED_SEG)   # 96 s spacing -> adjacent segs
            f.write(f"{gps}\t{tmp_path}/H1L1-SID-{gps}.gwf\n")

    (tmp_path / "badGPS.txt").write_text("")   # empty bad-GPS list
    prefix = str(tmp_path) + "/FSID_"

    par = tmp_path / "foldSID.par"
    par.write_text(
        f"SIDMetaDataFile {meta}\n"
        f"segmentDuration {SEG_DURATION}\n"
        "ifo1 H1\n"
        "ifo2 L1\n"
        f"FSIDFramePrefix {prefix}\n"
        "segmentsPerFrame 1\n"
        f"FSIDGPSOffset {GPS0}\n"
        "version 1\n"
        "ovlWinCorrection 1\n"
        f"identicalNeighbors {identical_neighbors}\n"
        "backwardCompatible 1\n"
        f"FSIDJobFile {tmp_path}/jobfile.txt\n"
        f"badGPSTimesFile {tmp_path}/badGPS.txt\n"
    )
    return str(par)


@pytest.mark.parametrize("identical_neighbors", [1, 0])
def test_fold_end_to_end(tmp_path, identical_neighbors):
    frstub.reset()
    paramsFile = _write_inputs(tmp_path, identical_neighbors)

    foldedInvCov, foldedStatistic, params, foldParams = foldSID(
        paramsFile, writeFrames=1
    )

    nFreqBin = frstub.NFREQ
    nSegment = round(86164.1 / FOLDED_SEG)

    # Output shapes are consistent with the SID band + sidereal segmentation.
    assert foldedStatistic.shape == (nSegment, nFreqBin)
    assert foldedStatistic.dtype == complex
    for key in ("Diag", "Prev", "Next"):
        assert foldedInvCov[key].shape == (nSegment, nFreqBin)

    # Parameters read back from the (synthetic) first frame.
    assert params["flow"] == frstub.FLOW
    assert params["fhigh"] == frstub.FHIGH
    assert params["deltaF"] == frstub.DELTAF
    assert params["winFactor"] == pytest.approx(1.0)
    assert params["winRatio"] == pytest.approx(0.5)

    # Every good segment was folded somewhere: segDist sums to N_SEGMENTS.
    assert int(np.sum(foldParams["segDist"])) == N_SEGMENTS

    # Output is finite (no NaN/inf leaking from the folding arithmetic).
    assert np.all(np.isfinite(foldedInvCov["Diag"]))
    assert np.all(np.isfinite(foldedStatistic.view(float)))

    # Frames were actually written (frputvect called) with expected channels.
    assert len(frstub.written_frames) > 0
    chans = frstub.written_frames[0]["channels"]
    assert "H1L1:ReCSD" in chans and "H1L1:ImCSD" in chans


def test_diag_populated_where_segments_land(tmp_path):
    # Use the general (non-identicalNeighbors) branch: there the diagonal
    # accumulates isStationary * varSigmaSqInv (strictly positive) for every
    # segment that lands in a sidereal bin. (In the identicalNeighbors branch the
    # diagonal weight is 1 - prevExists - nextExists, which is 0 for mutually
    # adjacent segments, so that branch legitimately leaves Diag at zero here.)
    frstub.reset()
    paramsFile = _write_inputs(tmp_path, identical_neighbors=0)
    foldedInvCov, foldedStatistic, params, foldParams = foldSID(
        paramsFile, writeFrames=0
    )
    occupied = foldParams["segDist"].reshape(-1) > 0
    assert np.any(occupied)
    assert np.all(foldedInvCov["Diag"][occupied].sum(axis=1) > 0)
