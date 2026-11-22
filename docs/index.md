# ArgoPAL

```{image} img/argopal_logo.png
:alt: ArgoPAL
:width: 480px
:align: center
```

**Wind from passive acoustics on BGC-Argo floats.**

```{note}
ArgoPAL is under active development.  Functions, options and results may
change between versions, and the processing has so far been tested on two
floats only.  Check the [status and limitations](limitations.md) before
relying on its output, and report problems on the
[issue tracker](https://github.com/louisedelaigue/ArgoPAL/issues).
```

ArgoPAL starts from a float's WMO number, downloads its passive acoustic
(PAL) spectra and its CTD/BGC data from the Argo GDAC, and processes them
into surface wind speed: cleaning and despiking, depth correction from the
float's own profiles, published wind models, and calibration against
reference buoys with uncertainty.  Ship noise can be flagged with AIS.

```python
import argopal as ap

result = ap.process_float(7902280)                 # download + process
result = ap.process_float(
    7902280,
    reference="reference_winds.csv",               # buoy winds, optional
)
figs = ap.plot_all(result, save_dir="figures")     # one figure per step
ap.export(result, "wmo7902280.nc", format="netcdf")
```

```{toctree}
:maxdepth: 2
:caption: Getting started

installation
tutorial
limitations
```

```{toctree}
:maxdepth: 2
:caption: User guide

data
processing
wind
ais
figures
```

```{toctree}
:maxdepth: 2
:caption: Reference

api
versions
development
```
