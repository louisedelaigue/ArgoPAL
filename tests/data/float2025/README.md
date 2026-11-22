# 2025 sea-trial float: test data

A second float for the tests, alongside 7902280. It is **not on the Argo
GDAC** (no WMO number), so these files are copied from the study's own data
and are used only by the tests, never by argopal itself.

Float: PROVOR_V, serial P53360-23FR001, Iridium login `lovuse024c`,
deployed near the DYFAMED buoy (Ligurian Sea) in February 2025 —
Delaigue et al. (2026), *Ocean Science* 22, 101–117,
https://doi.org/10.5194/os-22-101-2026.  Source: the study repository
(WIND-FROM-FLOAT-ACOUSTICS, `data/float/`, `data/AIS/`).

What makes it useful as a second float:

- full 23-band spectra, including 3.15 kHz (7902280 sends mostly 9 bands);
- band levels **not** normalised by bandwidth on board (7902280 normalises);
- three park depths: cycle 41 ≈ 500 m, cycle 42 ≈ 740 m, cycle 45 ≈ 1020 m;
- real AIS vessel positions along the track (none exist for 7902280).

| file | content |
|---|---|
| `pal.csv.gz` | PAL park-drift spectra of cycles 41, 42, 45 (`f_63Hz` … `f_25000Hz`, dB re 1 µPa, band levels) |
| `ctd.csv.gz` | SBE41 T/S of the same cycles: ascent profiles (`ASC`) and park samples (`PAR`), raw |
| `positions.csv` | surface GPS fixes (RData `technical$GPS`), cycles 41 (incl. pattern 0, before the first park) to 45 |
| `ais_positions.csv.gz` | decoded AIS positions (the study's `ais_positions_YYYYMMDD.csv` columns), 14 Feb 12:30 – 16 Feb 02:15 UTC, within 25 km of 43.36°N 7.77°E; further away, one position per minute only (keeps the record of when the AIS feed was receiving) |
| `ais_messages.nmea.gz` | the same, as raw NMEA sentences with tag blocks (`messages-YYYYMMDD.gz`), plus one base-station report (type 4) |
| `dyfamed_wind.csv` | DYFAMED buoy hourly wind (Météo-France station 6100001, `ff`, m s⁻¹), 9–17 Feb 2025 |

Rebuild with `python tests/make_fixtures_2025.py <WIND-FROM-FLOAT-ACOUSTICS>`
(needs R for the GPS fixes and pyais for the NMEA slice).
