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
| 23 | 63, 100, 125, 160, 400, 500, 630, 800, 1000, 1250, 1600, 2000, 2500, 3150, 4000, 5000, 6300, 8000, 10000, 12500, 16000, 20000, 25000 | `argopal.PAL_BANDS_HZ` |
| 9 | 63, 125, 400, 1000, 2000, 5000, 8000, 12500, 20000 | `argopal.PAL_BANDS_9_HZ` |

Any other number of values is refused with an error rather than given
wrong frequencies. Each read also checks that the median spectrum falls
with frequency, as ambient noise does, and warns (`argopal.BandOrderWarning`)
if it does not, which would mean the bands are in the wrong order.

## Band normalisation

`meta["band_normalization"]` reports `CONFIG_PalBandNormalization_LOGICAL`,
which says whether the float divides each band level by its bandwidth
before sending it. `spectral_level` reads it and subtracts the bandwidth
term only from floats that did not normalise, so the correction is never
applied twice; `band_normalization=` overrides it.

## Which floats have a PAL

No GDAC index records which floats carry a PAL, so `find_pal_floats`
looks for the PAL spectra in every float folder of the GDAC's `aux/`
directory (a couple of minutes; the list is saved in
`argopal_data/pal_floats.csv` and reused for a week):

```python
floats = argopal.find_pal_floats()          # wmo, dac, dates, first/last positions
index = argopal.read_profile_index(floats["wmo"])  # every profile: date, position
```

## Download and read

```python
pal, meta = argopal.fetch_pal(7902280)
```

This downloads the float's auxiliary files from
`https://data-argo.ifremer.fr/aux/<dac>/<wmo>/` into
`argopal_data/<wmo>/` (re-downloading only files that changed on the
server) and returns:

- `pal`: a DataFrame with one row per acoustic measurement, with columns
  `cycle`, `measurement_code`, `datetime`, `pressure` and the
  third-octave band levels `f_63Hz` … `f_25000Hz` (dB re 1 µPa; 23 bands,
  or 9 for floats sending the reduced spectrum);
- `meta`: the hydrophone details and the float's PAL configuration
  (`CONFIG_Pal*` parameters).

The steps can also be run separately:

```python
argopal.find_dac(7902280)              # "coriolis"
paths = argopal.download_aux(7902280)  # aux files, incl. Rtraj_aux.nc
pal = argopal.read_pal(paths["7902280_Rtraj_aux.nc"])
meta = argopal.read_pal_meta(paths["7902280_meta_aux.nc"])
argopal.download_core(7902280)  # meta, Rtraj, Sprof (CTD, BGC, GPS)
```
