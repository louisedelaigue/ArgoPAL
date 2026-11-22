# Development

## Install for development

    git clone https://github.com/louisedelaigue/ArgoPAL.git
    cd ArgoPAL
    pip install -e ".[ais,docs]"
    pip install pytest

## Tests

    pytest                            # offline tests
    pytest -m network                 # tests that download from the GDAC
    sh .misc/update_coverage.sh       # all tests with coverage; updates the badge

The coverage badge and report (`.misc/coverage.svg`, `.misc/coverage.txt`)
are made locally with pytest-cov and genbadge, from all tests including
the GDAC ones; rerun the script before committing changes to the code.

The tests run on two floats:

- **7902280** (GDAC): small slices of its aux and core files in
  `tests/data/` (rebuild: `tests/make_fixtures.py`);
- **the 2025 sea-trial float** of Delaigue et al. (2026), which is not on the
  GDAC: its spectra, CTD, GPS fixes and real AIS around its track are in
  `tests/data/float2025/` (see the README there; rebuild:
  `tests/make_fixtures_2025.py`).  It covers full 23-band spectra
  (3.15 kHz), band levels not normalised on board, three park depths
  (500, 740, 1020 m) and real AIS (decoded CSV and raw NMEA).

The β port is also checked against the published MATLAB depth correction
(`tests/data/matlab_ligurian_*`).

## Documentation

This site is built with Sphinx from `docs/`, on Read the Docs at every
push.  To build it locally:

    pip install -e ".[docs]"
    sphinx-build -b html docs docs/_build/html

The tutorial pages are copies of the notebooks in `examples/`, shown with
their saved outputs.  Each notebook is built by a script,
`examples/make_*.py`: after changing it, run the script, then run the
notebook from `examples/` so that its outputs are saved.
