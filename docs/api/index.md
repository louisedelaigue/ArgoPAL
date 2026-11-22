# Function reference

Every public function, grouped by processing step; each has its own page
with its parameters, what it returns, and the source code.  All are
available from the top-level package:

```python
import argopal as ap
```

```{toctree}
:maxdepth: 2

pipeline
gdac
reading
bands
cleaning
correction
wind
calibration
ais
figures
export
```

## Constants

| name | value |
|---|---|
| `ap.PAL_BANDS_HZ` | centre frequencies (Hz) of the 23 PAL bands |
| `ap.PAL_BANDS_9_HZ` | centre frequencies (Hz) of the 9-band subset |
| `ap.BAND_LAYOUTS` | number of values per spectrum -> band frequencies |
| `ap.DEFAULT_BANDS` | bands processed by default (8 kHz) |
| `ap.MODES` | data modes: `"delayed"`, `"adjusted"`, `"raw"` |
| `ap.DEFAULT_MODE` | `"delayed"` |
| `ap.COEFFICIENTS` | published coefficients of each wind model |
| `ap.LITERATURE_MODELS` | wind model name -> (function, band, reference) |

## Warnings

```{eval-rst}
.. autoclass:: argopal.BandLayoutWarning
.. autoclass:: argopal.BandOrderWarning
.. autoclass:: argopal.BandCoverageWarning
```
