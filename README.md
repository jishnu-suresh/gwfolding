# gwfolding

<p align="center">
  <img src="https://raw.githubusercontent.com/jishnu-suresh/gwfolding/main/assets/logo.png" alt="gwfolding logo" width="360">
</p>

<p align="center">
  <a href="https://pypi.org/project/gwfolding/"><img src="https://img.shields.io/pypi/v/gwfolding.svg" alt="PyPI version"></a>
  <a href="https://pypi.org/project/gwfolding/"><img src="https://img.shields.io/pypi/pyversions/gwfolding.svg" alt="Python versions"></a>
  <a href="https://gwfolding.readthedocs.io/en/latest/"><img src="https://readthedocs.org/projects/gwfolding/badge/?version=latest" alt="Documentation status"></a>
  <a href="https://github.com/jishnu-suresh/gwfolding/blob/main/LICENSE"><img src="https://img.shields.io/pypi/l/gwfolding.svg" alt="License: MIT"></a>
</p>

Folding of **Stochastic Intermediate Data (SID)** frames to one sidereal day for
stochastic gravitational-wave searches. The output can be written as folded SID
(`.gwf`) frames.

Mathematical details are described in the LIGO technical document **T0900093**.

## Prerequisite: the `Fr` frame library (not on PyPI)

This package reads and writes gravitational-wave frame files via the `Fr` module
(`frgetvect` / `frputvect`). `Fr` is **not** distributed on PyPI and is **not**
installed by `pip`. It must already be available in your Python environment,
typically from the **conda / IGWN** software stack used for LIGO/Virgo analysis.

Install or activate your IGWN/conda environment so that the following works
before using this package:

```python
from Fr import frgetvect, frputvect
```

If that import fails, install the frame library from your conda/IGWN channel
first.

## Installation

```bash
pip install gwfolding
```

This installs the pure-Python package and its PyPI dependencies (`numpy`,
`scipy`). The `Fr` frame library must be provided separately (see above).

## Usage

```python
from gwfolding import foldSID

foldedInvCov, foldedStatistic, params, foldParams = foldSID(
    "foldSID.par",   # path to the parameter file
    writeFrames=1,   # write folded SID frames
)
```

A worked example runner is provided in `run_foldSID.py`, and an example
parameter file in `foldSID.par.txt`.

## Authors

- **Authors:** Anirban Ain, Prathamesh Dalvi, Sanjit Mitra
- **Python translation:** Erik Floden
- **Data production & maintainer:** Jishnu Suresh (`jishnu.suresh@ligo.org`)

## License

MIT. See `LICENSE`.
