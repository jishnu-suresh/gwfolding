# Installation

## Prerequisite: the `Fr` frame library (not on PyPI)

gwfolding reads and writes gravitational-wave frame files through the `Fr`
module (`frgetvect` / `frputvect`). `Fr` is **not** distributed on PyPI and is
**not** installed by `pip`; it must already be available in your Python
environment, typically from the **conda / IGWN** software stack used for
LIGO/Virgo analysis.

Before using gwfolding, confirm this import works:

```python
from Fr import frgetvect, frputvect
```

If it fails, install or activate your IGWN/conda environment first.

## Install from PyPI

```bash
pip install gwfolding
```

This installs the pure-Python package and its PyPI dependencies (`numpy`,
`scipy`). The `Fr` frame library must be provided separately (see above).

## Install from source

```bash
git clone https://github.com/jishnu-suresh/gwfolding.git
cd gwfolding
pip install .
```

## Running the tests

The test suite stubs out `Fr`, so it runs without the frame library:

```bash
pip install .[test]
pytest
```
