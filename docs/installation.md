# Installation

```{include} ../README.md
---
start-after: "## Install"
end-before: "## Import"
relative-docs: "docs/"
relative-images: true
---
```

Then:

```python
import argopal
argopal.hello()  # prints the version
```

Requires Python 3.11 or later.  Dependencies: gsw, matplotlib, netCDF4,
numpy, pandas, scipy and xarray; pyais for raw NMEA AIS files (the `ais`
extra).
