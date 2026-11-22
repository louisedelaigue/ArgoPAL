# %%
"""Finding PAL floats on the GDAC (requests faked) and the profile index."""

import os
import shutil

import pandas as pd
import pytest

import argopal
from argopal import search


DATA = os.path.join(os.path.dirname(__file__), "data")
INDEX = os.path.join(DATA, "synthetic_profile_index_sample.txt.gz")
GDAC = "https://gdac.test"

# a fake GDAC: two DACs; 2904023 has PAL, 1900722 has an Rtraj_aux file
# without PAL, 3900001 has no Rtraj_aux file at all
LISTINGS = {
    f"{GDAC}/aux/": '<a href="coriolis/">c</a> <a href="aoml/">a</a>',
    f"{GDAC}/aux/coriolis/": '<a href="2904023/">x</a> <a href="?C=N">s</a>',
    f"{GDAC}/aux/aoml/": '<a href="1900722/">x</a> <a href="3900001/">y</a>',
    f"{GDAC}/aux/coriolis/2904023/": '<a href="2904023_Rtraj_aux.nc">f</a>',
    f"{GDAC}/aux/aoml/1900722/": '<a href="1900722_Rtraj_aux.nc">f</a>',
    f"{GDAC}/aux/aoml/3900001/": '<a href="3900001_tech_aux.nc">f</a>',
}
HEADS = {
    f"{GDAC}/aux/coriolis/2904023/2904023_Rtraj_aux.nc": b"\x89HDF..."
    + search.PAL_VARIABLE
    + b"...",
    f"{GDAC}/aux/aoml/1900722/1900722_Rtraj_aux.nc": b"\x89HDF...DOXY...",
}


@pytest.fixture
def fake_gdac(monkeypatch):
    calls = []

    def get(url, nbytes=None):
        calls.append((url, nbytes))
        if url in LISTINGS:
            return LISTINGS[url].encode()
        assert nbytes == search.HEAD_BYTES  # files are only read in part
        return HEADS[url]

    def download_file(url, fname, update=True):
        assert url == f"{GDAC}/{search.PROFILE_INDEX}"
        os.makedirs(os.path.dirname(fname), exist_ok=True)
        shutil.copyfile(INDEX, fname)
        return True

    monkeypatch.setattr(search, "_get", get)
    monkeypatch.setattr(search, "download_file", download_file)
    return calls


def test_read_profile_index(fake_gdac, tmp_path):
    idx = argopal.read_profile_index(data_dir=tmp_path, gdac=GDAC)
    assert list(idx.columns) == [
        "wmo", "cycle", "direction", "dac", "datetime", "latitude",
        "longitude", "ocean", "profiler_type", "file",
    ]  # fmt: skip
    assert set(idx["wmo"]) == {2904023, 1900722}
    one = idx[(idx["wmo"] == 2904023) & (idx["cycle"] == 1)]
    assert set(one["direction"]) == {"A", "D"}
    first = one[one["direction"] == "A"].iloc[0]
    assert first["datetime"] == pd.Timestamp("2025-07-13 16:56:53")
    assert (first["latitude"], first["longitude"]) == (39.526, 12.926)
    only = argopal.read_profile_index(2904023, data_dir=tmp_path, gdac=GDAC)
    assert set(only["wmo"]) == {2904023}


def test_find_pal_floats(fake_gdac, tmp_path):
    found = argopal.find_pal_floats(
        data_dir=tmp_path, gdac=GDAC, workers=2, verbose=False
    )
    assert list(found["wmo"]) == [2904023]
    assert found["dac"].iloc[0] == "coriolis"
    assert found.attrs["n_checked"] == 3
    row = found.iloc[0]
    assert row["n_profiles"] == 9  # incl. the descending profile
    assert row["first_date"] == pd.Timestamp("2025-07-13 12:00:00")
    # saved, and reused without scanning again
    assert os.path.isfile(tmp_path / "pal_floats.csv")
    n_calls = len(fake_gdac)
    again = argopal.find_pal_floats(
        data_dir=tmp_path, gdac=GDAC, summary=False, verbose=False
    )
    assert len(fake_gdac) == n_calls
    assert list(again.columns) == ["wmo", "dac"]
    # unless asked to, or too old
    argopal.find_pal_floats(
        data_dir=tmp_path, gdac=GDAC, refresh=True, summary=False,
        verbose=False,
    )  # fmt: skip
    assert len(fake_gdac) > n_calls


def test_unreadable_folder_is_counted(fake_gdac, monkeypatch, tmp_path):
    import urllib.error

    real = search._get

    def flaky(url, nbytes=None):
        if url.endswith("/1900722/"):
            raise urllib.error.URLError("timeout")
        return real(url, nbytes)

    monkeypatch.setattr(search, "_get", flaky)
    found = argopal.find_pal_floats(
        data_dir=tmp_path, gdac=GDAC, summary=False, verbose=False
    )
    assert list(found["wmo"]) == [2904023]
    assert found.attrs["n_failed"] == 1


@pytest.mark.network
def test_profile_index_from_gdac(tmp_path):
    idx = argopal.read_profile_index(7902280, data_dir=tmp_path)
    assert len(idx) > 50
    assert idx["datetime"].min() >= pd.Timestamp("2026-01-01")
