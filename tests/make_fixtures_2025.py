"""Rebuild tests/data/float2025 from the 2025 study repository.

    python tests/make_fixtures_2025.py /path/to/WIND-FROM-FLOAT-ACOUSTICS

Needs R (Rscript) to read the surface GPS fixes from the RData files, and
pyais to cut the raw NMEA slice.  See tests/data/float2025/README.md.
"""

import gzip
import os
import subprocess
import sys

import numpy as np
import pandas as pd


OUT = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "data", "float2025"
)
CYCLES_PAL = (41, 42, 45)  # park ≈ 500, 740, 1020 m
CYCLES_GPS = ((41, 0), (41, 1), (42, 1), (43, 1), (44, 1), (45, 1))
AIS_T0, AIS_T1 = (
    pd.Timestamp("2025-02-14 12:30"),
    pd.Timestamp("2025-02-16 02:15"),
)
AIS_CENTRE, AIS_RADIUS_KM = (43.36, 7.77), 25.0
AIS_DAYS = ("20250215", "20250216", "20250217")  # files span the day before
POSITION_TYPES = (1, 2, 3, 18, 19, 27)


def float_tables(repo):
    pal, ctd = [], []
    for c in CYCLES_PAL:
        d = pd.read_csv(
            os.path.join(
                repo, "data", "float", "csv", f"lovuse024c_{c:03d}_01.csv"
            ),
            low_memory=False,
        )
        d["datetime"] = pd.to_datetime(d["Date"], format="mixed")
        d["pressure"] = pd.to_numeric(d["Pressure_dbar"], errors="coerce")
        d["cycle"] = c
        bands = [k for k in d if k.startswith("f_")]
        p = d[(d["SensorType"] == "pal") & (d["PhaseName"] == "PAR")]
        pal.append(p[["cycle", "datetime", "pressure", *bands]])
        s = d[d["SensorType"] == "sbe41"].rename(
            columns={
                "PhaseName": "phase",
                "Temperature_degC": "temperature",
                "Salinity_PSU": "salinity",
            }
        )
        ctd.append(
            s[
                [
                    "cycle",
                    "datetime",
                    "phase",
                    "pressure",
                    "temperature",
                    "salinity",
                ]
            ]
        )
    pd.concat(pal).sort_values("datetime").to_csv(
        os.path.join(OUT, "pal.csv.gz"), index=False
    )
    pd.concat(ctd).sort_values("datetime").to_csv(
        os.path.join(OUT, "ctd.csv.gz"), index=False
    )


def gps_fixes(repo):
    rows = []
    for c, pattern in CYCLES_GPS:
        f = os.path.join(
            repo, "data", "float", f"lovuse024c_{c:03d}_{pattern:02d}.RData"
        )
        expr = (
            f'e<-new.env(); load("{f}", envir=e); g<-e$dataprofile$technical$GPS; '
            'cat(format(g$time, "%Y-%m-%d %H:%M:%S"), '
            'sprintf("%.6f %.6f", g[["lat (deg)"]], g[["lon (deg)"]]))'
        )
        out = subprocess.run(
            ["Rscript", "-e", expr], capture_output=True, text=True, check=True
        ).stdout.split()
        rows.append(
            (c, pattern, f"{out[0]} {out[1]}", float(out[2]), float(out[3]))
        )
    pd.DataFrame(
        rows, columns=["cycle", "pattern", "datetime", "latitude", "longitude"]
    ).to_csv(os.path.join(OUT, "positions.csv"), index=False)


def _box():
    lat0, lon0 = AIS_CENTRE
    dlat = AIS_RADIUS_KM / 111.2
    dlon = AIS_RADIUS_KM / (111.2 * np.cos(np.radians(lat0)))
    return lon0 - dlon, lon0 + dlon, lat0 - dlat, lat0 + dlat


def ais_tables(repo):
    from pyais.exceptions import AISBaseException
    from pyais.stream import IterMessages

    w, e, s, n = _box()
    cols = [
        "id", "mmsi", "x", "y", "timestamp", "tagblock_timestamp",
        "utc_hour", "utc_min", "file_year", "file_month", "file_day",
    ]  # fmt: skip
    parts = []
    for day in AIS_DAYS:
        df = pd.read_csv(
            os.path.join(
                repo, "data", "AIS", "parsed", f"ais_positions_{day}.csv"
            ),
            usecols=cols,
            low_memory=False,
        )
        t = pd.to_datetime(df["tagblock_timestamp"], unit="s")
        window = t.between(AIS_T0, AIS_T1)
        inside = window & df["x"].between(w, e) & df["y"].between(s, n)
        # Outside the box, one position per minute: keeps the record of when
        # the AIS feed was receiving (coverage) without the bulk.
        outside = window & ~inside & df["id"].isin(POSITION_TYPES)
        first = ~t[outside].dt.floor("min").duplicated()
        parts.append(df[inside])
        parts.append(df[outside][first.to_numpy()])
    pd.concat(parts).to_csv(
        os.path.join(OUT, "ais_positions.csv.gz"), index=False
    )

    keep_lines, minutes_seen, base_added = [], set(), False
    for day in AIS_DAYS:
        path = os.path.join(repo, "data", "AIS", f"messages-{day}.gz")
        with gzip.open(path, "rb") as f:
            for raw in (line.rstrip(b"\r\n") for line in f):
                if not raw.startswith(b"\\c:"):
                    continue  # multi-sentence continuation: not a position
                for msg in IterMessages([raw]):
                    try:
                        d = msg.decode()
                    except (AISBaseException, ValueError):
                        continue
                    msg.tag_block.init()
                    t = pd.Timestamp(
                        int(msg.tag_block.receiver_timestamp), unit="s"
                    )
                    if not AIS_T0 <= t <= AIS_T1:
                        continue
                    if d.msg_type not in (*POSITION_TYPES, 4):
                        continue  # no position (e.g. type 5, 24)
                    inside = w <= d.lon <= e and s <= d.lat <= n
                    if d.msg_type in POSITION_TYPES and inside:
                        keep_lines.append(raw)
                    elif d.msg_type in POSITION_TYPES:
                        if t.floor("min") not in minutes_seen:  # coverage
                            minutes_seen.add(t.floor("min"))
                            keep_lines.append(raw)
                    elif d.msg_type == 4 and not base_added:
                        keep_lines.append(raw)  # one base station
                        base_added = True
    with gzip.open(os.path.join(OUT, "ais_messages.nmea.gz"), "wb") as f:
        f.write(b"\n".join(keep_lines) + b"\n")


if __name__ == "__main__":
    repo = sys.argv[1]
    os.makedirs(OUT, exist_ok=True)
    float_tables(repo)
    gps_fixes(repo)
    ais_tables(repo)
