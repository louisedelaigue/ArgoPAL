# Data from the Argo GDAC

## Where the acoustic data come from

NKE floats fitted with a PAL (Passive Acoustic Listener) hydrophone send
third-octave sound pressure level spectra. Coriolis decodes them into the
GDAC **auxiliary** folder, not the core one:

| File | Contents |
|---|---|
| `aux/<dac>/<wmo>/<wmo>_Rtraj_aux.nc` | `PASSIVE_ACOUSTICS_THIRD_OCTAVE_LEVEL_SPECTRUM` (23 or 9 bands, dB re 1 µPa) and `PASSIVE_ACOUSTICS_BAND_VARIANCE_CLASSIFICATION`, per measurement, with `JULD`, `PRES`, `CYCLE_NUMBER`, `MEASUREMENT_CODE` |
| `aux/<dac>/<wmo>/<wmo>_meta_aux.nc` | Hydrophone model and serial number, `CONFIG_Pal*` sampling and processing parameters |
| `aux/<dac>/<wmo>/<wmo>_tech_aux.nc` | Technical parameters |

The band centre frequencies are not stored in the files, so ArgoPAL
assigns them from the number of values per spectrum. The PAL sends one of
two layouts:

| Values | Bands (Hz) | Constant |
|---|---|---|
| 23 | 63, 100, 125, 160, 400, 500, 630, 800, 1000, 1250, 1600, 2000, 2500, 3150, 4000, 5000, 6300, 8000, 10000, 12500, 16000, 20000, 25000 | `ap.PAL_BANDS_HZ` |
| 9 | 63, 125, 400, 1000, 2000, 5000, 8000, 12500, 20000 | `ap.PAL_BANDS_9_HZ` |

Any other number of values is refused with an error rather than given
wrong frequencies. Each read also checks that the median spectrum falls
with frequency, as ambient noise does, and warns (`ap.BandOrderWarning`)
if it does not, which would mean the bands are in the wrong order.

## Band normalisation

`meta["band_normalization"]` reports `CONFIG_PalBandNormalization_LOGICAL`,
which says whether the float divides each band level by its bandwidth
before sending it. `spectral_level` reads it and subtracts the bandwidth
term only from floats that did not normalise, so the correction is never
applied twice; `band_normalization=` overrides it.

## Download and read

```{include} ../README.md
---
start-after: "## Get the acoustic data of a float"
end-before: "## Process a float: PAL + all other float data"
relative-docs: "docs/"
relative-images: true
---
```
