# Ship noise (AIS)

Ship noise can dominate the 1-10 kHz band up to ~20 km away.  Given AIS
vessel positions, `process_float` flags each PAL measurement with a ship
within 20 km of the float and ±30 min (Delaigue et al. 2026, Sect. 2.5):

```python
result = argopal.process_float(7902280, ais="AIS/messages-*.gz")            # flag only
result = argopal.process_float(7902280, ais=my_table, ais_action="remove")  # blank flagged values
```

Adds `ships_nearby` (vessels counted once per MMSI; NaN when the AIS feed
received nothing at all in the window) and `ship_flag`.  `ais_radius_km`
and `ais_window` change the criteria.  Figure: `argopal.plot_ais(result["park"])`.

**What AIS data is needed.** Vessel positions covering the float's area
and dates, as any of:

| format | what ArgoPAL needs |
|---|---|
| pandas DataFrame | columns `datetime` (UTC), `mmsi`, `latitude`, `longitude` |
| CSV of decoded positions (`.csv`, `.csv.gz`) | a time column (`datetime`, `time`, `BaseDateTime`, `# Timestamp`, `tagblock_timestamp` = UNIX s, …), `mmsi`, latitude (`latitude`/`lat`/`y`), longitude (`longitude`/`lon`/`x`). Covers NOAA MarineCadastre, Danish Maritime Authority and the decoded files of the 2025 DYFAMED study. |
| raw NMEA 0183 (`!AIVDM`, plain or `.gz`) | NMEA 4.0 tag blocks with the receiver time, e.g. `\c:1739500274,s:2954*14\!AIVDM,...`; needs `pip install argopal[ais]` (pyais) |

Only vessel position reports (message types 1, 2, 3, 18, 19, 27) count:
base stations and aids to navigation are dropped when the CSV has a
message-type column.  Times must be the time each position was received:
do not date positions from file names, since daily AIS files often span
two calendar days.
