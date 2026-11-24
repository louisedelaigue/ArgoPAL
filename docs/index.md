# ArgoPAL

**Wind from passive acoustics on BGC-Argo floats.**

```{note}
ArgoPAL is under active development: functions, options and results may
change between versions, and it has been tested on five floats so far.
See the [status and limitations](limitations.md), and report problems on
the [issue tracker](https://github.com/louisedelaigue/ArgoPAL/issues).
```

ArgoPAL takes a float's WMO number, downloads its passive acoustic (PAL)
spectra and its CTD/BGC data from the Argo GDAC, and turns them into
surface wind speed, following the method of Delaigue et al. (2026).

```python
import argopal

result = argopal.process_float(7902280)
result["park"]  # one row per acoustic measurement, with wind speed
```

## Installation

With Python 3.11 or later:

```bash
pip install git+https://github.com/louisedelaigue/ArgoPAL.git
```

For raw NMEA AIS files, install the `ais` extra:

```bash
pip install "argopal[ais] @ git+https://github.com/louisedelaigue/ArgoPAL.git"
```

`argopal.hello()` prints the installed version.

## How to use ArgoPAL

Start with the [quick-start guide](quickstart.md), then look up
[all options](options.md) and [all outputs](outputs.md).
[How it works](method.md) explains each processing step, the
[tutorials](tutorials/multi_float_tutorial.ipynb) run through real floats, and
the [API reference](api.md) lists every function.

## Citation

The method is described in:

> Delaigue, L., Cauchy, P., Cazau, D., Bonnel, J., Pensieri, S.,
> Bozzano, R., Gros-Martial, A., Schaeffer, C., David, A., Stil, P.,
> Poteau, A., Schmechtig, C., Leymarie, E. and Claustre, H. (2026).
> Passive acoustic monitoring from profiling floats as a pathway to
> scalable autonomous observations of global surface wind.  *Ocean
> Science*, 22, 101–117.
> [doi:10.5194/os-22-101-2026](https://doi.org/10.5194/os-22-101-2026)

Work using ArgoPAL should cite this paper and the software, with the
version used:

> Delaigue, L. (2026).  ArgoPAL: wind speed from passive acoustics on
> BGC-Argo floats (version 0.1.0b2).
> [github.com/louisedelaigue/ArgoPAL](https://github.com/louisedelaigue/ArgoPAL)

```{toctree}
:hidden:
:caption: User guide

quickstart
options
outputs
method
figures
```

```{toctree}
:hidden:
:caption: Tutorials

tutorials/multi_float_tutorial
tutorials/pal_data_tutorial
tutorials/argopal_tutorial
tutorials/depth_correction_tutorial
```

```{toctree}
:hidden:
:caption: Reference

api
limitations
versions
development
```
