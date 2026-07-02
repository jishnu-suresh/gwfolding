# Usage

## Folding SID frames

```python
from gwfolding import foldSID

foldedInvCov, foldedStatistic, params, foldParams = foldSID(
    "foldSID.par",   # path to the parameter file
    writeFrames=1,   # write folded SID frames
)
```

`foldSID` returns:

- `foldedInvCov` — dictionary of inverse-covariance arrays (`Diag`, `Prev`, `Next`).
- `foldedStatistic` — complex `(nSegment, nFreqBin)` array of the folded statistic.
- `params` — parameters read from the first SID frame.
- `foldParams` — parameters used for folding.

A worked example runner is provided in `run_foldSID.py`, and an example
parameter file in `foldSID.par.txt`.

## The parameter file

The parameter file lists one `name value` pair per line. Required entries are
`SIDMetaDataFile` and `segmentDuration` (plus `ifo1` / `ifo2`); the remaining
entries are optional. See `foldSID.par.txt` in the repository for a complete
example.
