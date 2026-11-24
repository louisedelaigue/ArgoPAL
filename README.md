<p align="center">
  <img src="docs/img/argopal_logo.png" alt="ArgoPAL" width="480">
</p>

[![Tests](https://github.com/louisedelaigue/ArgoPAL/actions/workflows/run-tests.yml/badge.svg)](https://github.com/louisedelaigue/ArgoPAL/actions/workflows/run-tests.yml)
[![Coverage](.misc/coverage.svg)](.misc/coverage.txt)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

Wind from passive acoustics on BGC-Argo floats.

> [!NOTE]
> ArgoPAL is under active development.  Functions, options and results
> may change between versions, and the processing has so far been tested
> on five floats only.  Check the
> [status and limitations](#status-and-limitations) before relying on its
> output, and report problems on the
> [issue tracker](https://github.com/louisedelaigue/ArgoPAL/issues).

Documentation: https://argopal.readthedocs.io

ArgoPAL starts from a float's WMO number, downloads its passive acoustic
(PAL) spectra from the Argo GDAC, and processes them into surface wind
speed.

## Install

From GitHub:

    pip install git+https://github.com/louisedelaigue/ArgoPAL.git
    pip install "argopal[ais] @ git+https://github.com/louisedelaigue/ArgoPAL.git"  # with raw NMEA AIS

From a local copy, for development:

    pip install -e .

## Import and run

```python
import argopal

result = argopal.process_float(7902280)  # download from the GDAC and process
result["park"]      # one row per acoustic measurement: band levels, depth
                    # correction, spectral level, wind speed, CTD/BGC, position
result["profiles"]  # the float's CTD/BGC profiles

argopal.plot_all(result, save_dir="figures")  # one figure per step
argopal.export(result, "wmo7902280.nc", format="netcdf")
```

With reference winds (e.g. buoys), the wind models are refitted and each
wind gets an uncertainty:

```python
result = argopal.process_float(7902280, reference="reference_winds.csv")
```

## Documentation and tutorials

Everything else is in the documentation, [argopal.readthedocs.io](https://argopal.readthedocs.io/en/latest/):

- tutorials:
  [finding PAL floats and processing several](https://argopal.readthedocs.io/en/latest/tutorials/multi_float_tutorial.html),
  [getting the PAL data](https://argopal.readthedocs.io/en/latest/tutorials/pal_data_tutorial.html) (every
  band, no further processing),
  [from PAL data to wind speed](https://argopal.readthedocs.io/en/latest/tutorials/argopal_tutorial.html),
  [understanding the depth correction](https://argopal.readthedocs.io/en/latest/tutorials/depth_correction_tutorial.html);
- the user guide: [quick start](https://argopal.readthedocs.io/en/latest/quickstart.html),
  [options](https://argopal.readthedocs.io/en/latest/options.html),
  [outputs](https://argopal.readthedocs.io/en/latest/outputs.html),
  [how it works](https://argopal.readthedocs.io/en/latest/method.html),
  [figures](https://argopal.readthedocs.io/en/latest/figures.html);
- the [API reference](https://argopal.readthedocs.io/en/latest/api.html) and the
  [version history](https://argopal.readthedocs.io/en/latest/versions.html).

The tutorials are also in [`examples/`](examples/).

## Status and limitations

Version 0.1.0b1 is a beta: the processing chain is complete and tested,
but function names and options may still change before 1.0.  Changes are
listed in the [version history](https://argopal.readthedocs.io/en/latest/versions.html).

- **Calibration needs a reference.**  The published wind models were not
  made for floats and are biased on them (e.g. +1.4 m s⁻¹ for Nystuen on
  7902280, against three Météo-France buoys at 10 m).  Refitting needs reference winds (buoys) that the float
  passes close to; ArgoPAL does not download them.  A reference needs
  float data in 4 cycles to be calibrated (`calibration_min_cycles`); the
  optional bootstrap uncertainty of the coefficients needs about 10.  A
  refit is not to be trusted far outside the winds seen near the
  reference.
- **Low winds.**  Below about 2–3 m s⁻¹, few waves break and other noise
  sources set the 8 kHz level: the acoustic wind levels off while the
  true wind goes to zero.  Low acoustic winds are best read as "below
  about 3 m s⁻¹".
- **High winds.**  Above about 15 m s⁻¹, bubble clouds absorb part of
  the sound and floats usually see few such winds near a reference: the
  refitted wind tends to be low there (slope below 1 against the
  reference).
- **Delayed mode lags.**  The default mode uses delayed-mode profiles
  only; recent profiles are not in delayed mode yet, so the depth
  correction, and the wind, stop where delayed-mode data stop.
  `mode="raw"` covers the whole record with unadjusted data.
- **8 kHz only for wind.**  The three wind models are 8 kHz models.  Other
  bands are processed and depth-corrected but not turned into wind.
- **Band layouts.**  Only 23-band spectra and the 9-band subset are
  recognised; any other layout raises an error.
- **Ship noise needs AIS supplied by the user.**  Without AIS, ship noise
  is reduced by despiking only.
- **Rain is not flagged.**  Rain raises the 8 kHz level and the wind
  estimate with it.
- **Flat, loud spectra are not flagged.**  Some cycles hold spectra that
  are loud and flat across frequency (e.g. cycle 17 of 2904023), which is
  not wind noise; they give unrealistic winds (tens of m s⁻¹).
  Despiking removes short events, not whole cycles.
- **Tested floats.**  Tested on the four floats with PAL data on the GDAC
  (2904023, 2904024 and 7902280 in the Mediterranean, 7902388 in the
  eastern tropical Pacific; all NKE PROVOR with HTI hydrophones) and on
  the 2025 sea-trial float; other configurations may need adjustments.

## License

ArgoPAL is free software under the
[GNU General Public License v3.0](LICENSE).  It is developed and
maintained by Louise Delaigue.

## Citation

Work using ArgoPAL should cite the software:

> Delaigue, L. (2026). ArgoPAL: wind speed from passive acoustics on
> BGC-Argo floats (version 0.1.0b1).
> https://github.com/louisedelaigue/ArgoPAL

and, for the method:

> Delaigue, L., Cauchy, P., Cazau, D., Bonnel, J., Pensieri, S.,
> Bozzano, R., Gros-Martial, A., Schaeffer, C., David, A., Stil, P.,
> Poteau, A., Schmechtig, C., Leymarie, E. and Claustre, H. (2026).
> Passive acoustic monitoring from profiling floats as a pathway to
> scalable autonomous observations of global surface wind.  *Ocean
> Science*, 22, 101–117. https://doi.org/10.5194/os-22-101-2026
