# %%
import os

import numpy as np
import pandas as pd
import pytest

import argopal


DATA = os.path.join(os.path.dirname(__file__), "data")
TRAJ = os.path.join(DATA, "7902280_Rtraj_aux_first200.nc")
META = os.path.join(DATA, "7902280_meta_aux.nc")


def test_hello():
    argopal.hello()


def test_read_pal_shape_and_columns():
    pal = argopal.read_pal(TRAJ)
    assert len(pal) == 200
    bands = [f"f_{f}Hz" for f in argopal.PAL_BANDS_HZ]
    assert len(bands) == 23
    for col in ["cycle", "measurement_code", "datetime", "pressure"] + bands:
        assert col in pal.columns
    assert pd.api.types.is_datetime64_any_dtype(pal["datetime"])
    assert pal.attrs["wmo"] == 7902280


def test_read_pal_values_match_lov_decoding():
    # First park-drift spectrum of cycle 1, as in the LOV RData file
    # eriuse004b_013_01 (2026-01-12 04:49:00, 974.7 dbar).
    pal = argopal.read_pal(TRAJ)
    first = pal.iloc[0]
    assert first["cycle"] == 1
    assert first["measurement_code"] == 290
    assert abs(first["pressure"] - 974.7) < 0.01
    assert abs(first["datetime"] - pd.Timestamp("2026-01-12 04:49:00")) < (
        pd.Timedelta("10s")
    )
    assert first["f_63Hz"] == 86.5
    assert first["f_8000Hz"] == 48.0
    assert first["f_25000Hz"] == 36.5


def test_read_pal_sources():
    pal = argopal.read_pal(TRAJ)
    # This float has adjusted times but no adjusted pressure or cycles
    assert pal.attrs["time_source"] == "adjusted"
    assert pal.attrs["pressure_source"] == "raw"
    assert pal.attrs["cycle_source"] == "raw"


def test_read_pal_meta():
    meta = argopal.read_pal_meta(META)
    assert meta["wmo"] == 7902280
    assert meta["sensor"] == "HYDROPHONE"
    assert meta["sensor_model"] == "HTI_96_MIN"
    assert meta["band_normalization"] is True
    launch = meta["launch_config"]
    assert all(k.startswith("CONFIG_Pal") for k in launch)
    assert (
        launch["CONFIG_PalParkDriftPhaseDepthZone4SamplingPeriod_seconds"]
        == 300
    )
    assert isinstance(meta["mission_config"], pd.DataFrame)


def test_read_pal_with_meta_attrs():
    pal = argopal.read_pal(TRAJ, meta=argopal.read_pal_meta(META))
    assert pal.attrs["band_normalization"] is True
    assert pal.attrs["sensor_model"] == "HTI_96_MIN"


def test_read_pal_rejects_core_file(tmp_path):
    import xarray as xr

    fname = tmp_path / "no_pal.nc"
    xr.Dataset({"PRES": ("N_MEASUREMENT", np.zeros(3))}).to_netcdf(fname)
    with pytest.raises(KeyError):
        argopal.read_pal(str(fname))


# test_read_pal_values_match_lov_decoding()
