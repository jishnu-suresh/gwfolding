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
- **Data production and code maintenance:** Jishnu Suresh (`jishnu.suresh@ligo.org`)

## Citation

If you use this pipeline, please cite:

> A. Ain, P. Dalvi, and S. Mitra, "Fast gravitational wave radiometry using data
> folding," *Phys. Rev. D* **92** (July, 2015) 022003, arXiv:1504.01714 [gr-qc].

BibTeX:

```bibtex
@article{Ain:2015lea,
    author        = "Ain, Anirban and Dalvi, Prathamesh and Mitra, Sanjit",
    title         = "{Fast gravitational wave radiometry using data folding}",
    journal       = "Phys. Rev. D",
    volume        = "92",
    number        = "2",
    pages         = "022003",
    year          = "2015",
    doi           = "10.1103/PhysRevD.92.022003",
    eprint        = "1504.01714",
    archivePrefix = "arXiv",
    primaryClass  = "gr-qc"
}
```

## License

MIT. See `LICENSE`.
