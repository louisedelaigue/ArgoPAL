# argopal

Wind from passive acoustics on BGC-Argo floats.

## Where the acoustic data come from

NKE floats fitted with a PAL (Passive Acoustic Listener) hydrophone send
third-octave sound pressure level spectra. Coriolis decodes them into the
GDAC **auxiliary** folder, not the core one:

| File | Contents |
|---|---|
| `aux/<dac>/<wmo>/<wmo>_Rtraj_aux.nc` | `PASSIVE_ACOUSTICS_THIRD_OCTAVE_LEVEL_SPECTRUM` (23 or 9 bands, dB re 1 µPa) and `PASSIVE_ACOUSTICS_BAND_VARIANCE_CLASSIFICATION`, per measurement, with `JULD`, `PRES`, `CYCLE_NUMBER`, `MEASUREMENT_CODE` |
| `aux/<dac>/<wmo>/<wmo>_meta_aux.nc` | Hydrophone model and serial number, `CONFIG_Pal*` sampling and processing parameters |
| `aux/<dac>/<wmo>/<wmo>_tech_aux.nc` | Technical parameters |

The band centre frequencies are not stored in the files, so argopal
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
before sending it. Check it before converting band levels to spectral
levels (dB re 1 µPa²/Hz), so the correction is not applied twice.

## Depth correction

β is computed from the float's own mean T/S profile (in the chosen data
mode): pressure → depth and sound speed with TEOS-10 (`gsw`), absorption
with Francois & Garrison (1982) in dB/m, rays traced from the surface to
each depth (Vagle et al. 1990; Cauchy et al. 2018) out to 8 km range.  The
power lost along a ray is 10^(−αL/10); the original MATLAB code used
exp(−αL), which treated α as Np/m and overestimated β by ≈ 4.3×.
