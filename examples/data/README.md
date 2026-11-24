# Example data

`reference_winds_2026.csv`: hourly wind speed (m s⁻¹) at three Météo-France
moored buoys in the north-western Mediterranean, January–June 2026, with
each buoy's position at every value.  Source: Météo-France marine
observations (open data).  Used in `argopal_tutorial.ipynb` as the reference
wind for float 7902280.

| reference_name | Météo-France name | WMO | area |
|---|---|---|---|
| DYFAMED | BOUEE_AZUR | 6100001 | Ligurian Sea |
| LION | BOUEE_LION | 6100002 | Gulf of Lion |
| PACA | BOUEE_PACA | 6101036 | off Provence |

Columns: `reference_name`, `datetime` (UTC), `latitude`, `longitude`,
`wind_speed` (as measured, m s⁻¹) and `height_m`, the anemometer height:
3.5 m for all three buoys (Météo-France gives 3–4 m).  ArgoPAL converts
the winds to 10 m when reading them (`argopal.wind_to_10m`, neutral
logarithmic profile; a factor of about 1.10 at 3.5 m).  The buoys give no
error per value, so the file has no `wind_speed_sigma` column; the
tutorial passes one error for all of them instead (`reference_sigma`).
