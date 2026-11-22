"""Shared test fixtures."""

import os

import matplotlib
import pytest


matplotlib.use("Agg")  # no windows during tests

DATA = os.path.join(os.path.dirname(__file__), "data")


@pytest.fixture(autouse=True)
def _close_figures():
    yield
    import matplotlib.pyplot as plt

    plt.close("all")


@pytest.fixture
def offline(monkeypatch):
    """Make process_float use the small files in tests/data instead of
    downloading from the GDAC."""
    from argopal import pipeline

    files_aux = {
        "7902280_Rtraj_aux.nc": os.path.join(
            DATA, "7902280_Rtraj_aux_first200.nc"
        ),
        "7902280_meta_aux.nc": os.path.join(DATA, "7902280_meta_aux.nc"),
    }
    files_core = {
        "7902280_Rtraj.nc": os.path.join(DATA, "7902280_Rtraj_cycles1-2.nc"),
        "7902280_Sprof.nc": os.path.join(DATA, "7902280_Sprof_4profiles.nc"),
    }
    monkeypatch.setattr(pipeline, "find_dac", lambda *a, **k: "coriolis")
    monkeypatch.setattr(pipeline, "download_aux", lambda *a, **k: files_aux)
    monkeypatch.setattr(pipeline, "download_core", lambda *a, **k: files_core)
