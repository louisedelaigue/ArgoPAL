# %%
"""Tests that talk to the real GDAC.  Skipped by default; run with

pytest -m network
"""

import os

import pytest

import argopal as ap


WMO = 7902280
pytestmark = pytest.mark.network


def test_float_url():
    assert ap.gdac.float_url(WMO, "coriolis", aux=True) == (
        "https://data-argo.ifremer.fr/aux/coriolis/7902280/"
    )
    assert ap.gdac.float_url("7902280", "coriolis").endswith(
        "/dac/coriolis/7902280/"
    )


def test_find_dac():
    assert ap.find_dac(WMO) == "coriolis"


def test_find_dac_unknown():
    with pytest.raises(FileNotFoundError):
        ap.find_dac(1)


def test_list_aux_files():
    files = ap.list_files(WMO, dac="coriolis", aux=True)
    assert "7902280_Rtraj_aux.nc" in files
    assert "7902280_meta_aux.nc" in files


def test_fetch_pal(tmp_path):
    pal, meta = ap.fetch_pal(WMO, data_dir=str(tmp_path), dac="coriolis")
    assert len(pal) > 40000
    assert meta["sensor"] == "HYDROPHONE"
    traj = os.path.join(str(tmp_path), "7902280", "7902280_Rtraj_aux.nc")
    assert os.path.isfile(traj)
    # Second call: files are up to date, nothing is re-downloaded
    mtime = os.path.getmtime(traj)
    ap.fetch_pal(WMO, data_dir=str(tmp_path), dac="coriolis")
    assert os.path.getmtime(traj) == mtime


def test_process_float_and_export(tmp_path):
    result = ap.process_float(WMO, data_dir=str(tmp_path), dac="coriolis")
    park, prof = result["park"], result["profiles"]
    assert len(park) > 70000
    assert park.attrs["park"]["mode"] == "delayed"
    assert set(prof["TEMP_data_mode"]) <= {"D", ""}
    files = ap.export(result, str(tmp_path / "out.nc"), format="netcdf")
    back = ap.read_export(files[0])
    assert len(back["park"]) == len(park)


# test_fetch_pal(__import__("pathlib").Path("argopal_data_test"))
