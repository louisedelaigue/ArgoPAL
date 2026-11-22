"""Underwater acoustics: sound speed, absorption and the depth correction β.

The depth correction follows Vagle et al. (1990) as implemented by Cauchy
et al. (2018) in MATLAB (beta_h.m), with one fix: the absorption
coefficient of Francois & Garrison (1982) is in dB/m, so the power lost
along a ray of length L is 10**(-alpha * L / 10).  The MATLAB code used
exp(-alpha * L), treating alpha as if it were in Np/m, which overestimates
absorption by a factor of 10 / ln(10) ≈ 4.34 (e.g. β at 1000 m, 8 kHz:
≈ 0.96 dB instead of ≈ 4.0 dB).
"""

import numpy as np
import pandas as pd


# --- Sound speed and absorption --------------------------------------------


def sound_speed_chen_millero(depth, temperature, salinity):
    """Sound speed (m/s), Chen & Millero (1977), as in celeriteChen.m.

    Used only inside the Francois-Garrison absorption formula, as in the
    original code.  `depth` in m (positive down), `temperature` in °C,
    `salinity` practical.
    """
    t = np.asarray(temperature, dtype=float)
    s = np.asarray(salinity, dtype=float)
    p = np.asarray(depth, dtype=float) / 10  # bar, as in celeriteChen.m
    t2, t3, t4, t5 = t**2, t**3, t**4, t**5
    p2, p3 = p**2, p**3
    c0 = (
        1402.388
        + 5.03711 * t
        - 5.80852e-2 * t2
        + 3.3420e-4 * t3
        - 1.47800e-6 * t4
        + 3.1464e-9 * t5
    )
    c1 = (
        0.153563
        + 6.8982e-4 * t
        - 8.1788e-6 * t2
        + 1.3621e-7 * t3
        - 6.1185e-10 * t4
    )
    c2 = (
        3.1260e-5
        - 1.7107e-6 * t
        + 2.5974e-8 * t2
        - 2.5335e-10 * t3
        + 1.0405e-12 * t4
    )
    c3 = -9.7729e-9 + 3.8504e-10 * t - 2.3643e-12 * t2
    a0 = 1.389 - 1.262e-2 * t + 7.164e-5 * t2 + 2.006e-6 * t3 - 3.21e-8 * t4
    a1 = (
        9.4742e-5
        - 1.2580e-5 * t
        - 6.4885e-8 * t2
        + 1.0507e-8 * t3
        - 2.0122e-10 * t4
    )
    a2 = -3.9064e-7 + 9.1041e-9 * t - 1.6002e-10 * t2 + 7.988e-12 * t3
    a3 = 1.100e-10 + 6.649e-12 * t - 3.389e-13 * t2
    b = -1.922e-2 - 4.42e-5 * t + (7.3637e-5 + 1.7945e-7 * t) * p
    c = -7.9836e-6 * p + 1.727e-3
    a = a0 + a1 * p + a2 * p2 + a3 * p3
    return (
        c0 + c1 * p + c2 * p2 + c3 * p3 + a * s + b * s * np.sqrt(s) + c * s**2
    )


def absorption_francois_garrison(frequency_khz, depth, temperature, salinity):
    """Sound absorption coefficient (dB/m), Francois & Garrison (1982).

    Port of equ_att_son_garrison.m (with FBOH3.m and FMgSO4.m): boric
    acid, magnesium sulphate and pure-water viscosity contributions.

    Parameters
    ----------
    frequency_khz : float
        Frequency in kHz.
    depth : array_like
        Depth in m (positive down).
    temperature, salinity : array_like
        In-situ temperature (°C) and practical salinity.

    Returns
    -------
    np.ndarray
        Absorption in dB/m.
    """
    f = float(frequency_khz)
    z = np.asarray(depth, dtype=float)
    t = np.asarray(temperature, dtype=float)
    s = np.asarray(salinity, dtype=float)
    c = sound_speed_chen_millero(z, t, s)
    f2 = f * f
    # Boric acid
    a1 = 154 / c
    f1 = 2.8 * np.sqrt(s / 35) * 10 ** (4 - 1245 / (t + 273))
    # Magnesium sulphate
    a2 = 21.44 * s / c * (1 + 0.025 * t)
    p2 = 1 - 1.37e-4 * z + 6.2e-9 * z**2
    f_mg = 8.17 * 10 ** (8 - 1990 / (273 + t)) / (1 + 1.8e-3 * (s - 35))
    # Pure water
    a3 = np.where(
        t <= 20,
        4.937e-4 - 2.59e-5 * t + 9.11e-7 * t**2 - 1.5e-8 * t**3,
        3.964e-4 - 1.146e-5 * t + 1.45e-7 * t**2 - 6.5e-10 * t**3,
    )
    p3 = 1 - 3.83e-5 * z + 4.9e-10 * z**2
    a_db_km = (
        a1 * f1 * f2 / (f1 * f1 + f2)
        + a2 * p2 * f_mg * f2 / (f2 + f_mg * f_mg)
        + a3 * p3 * f2
    )
    return a_db_km / 1000


# --- Ray geometry -----------------------------------------------------------


_X_MAX = 0.999999  # clamp on cos of the local grazing angle, as in beta_h.m


def _theta_grid(n=20000):
    """Surface grazing angles (rad): dense near 0 (long, shallow rays) and
    evenly spread up to vertical."""
    small = np.geomspace(1e-6, 0.5, n)
    large = np.linspace(0.01, np.pi / 2 - 1e-9, n)
    return np.unique(np.concatenate([small, large]))


def ray_tables(c, dz, receiver_index, theta=None):
    """Horizontal range and path length of rays from the surface.

    For a ray leaving the surface at grazing angle theta, Snell's law gives
    cos(theta_j) / c_j = cos(theta) / c_1 in layer j.  Summing over the
    layers above a receiver gives its horizontal range r(theta) and path
    length L(theta), as in beta_h.m.

    Parameters
    ----------
    c : array
        Sound speed (m/s) in layers of thickness `dz`, from the surface
        down (c[0] = surface layer).
    dz : float
        Layer thickness (m).
    receiver_index : array of int
        Index (0-based) of the deepest layer above each receiver.
    theta : array, optional
        Surface grazing angles; by default `_theta_grid()`.

    Returns
    -------
    theta : array (n_theta,)
    r : array (n_receivers, n_theta), horizontal range (m)
    length : array (n_receivers, n_theta), path length (m)
    """
    c = np.asarray(c, dtype=float)
    if theta is None:
        theta = _theta_grid()
    receiver_index = np.asarray(receiver_index, dtype=int)
    r = np.zeros((len(receiver_index), len(theta)))
    length = np.zeros_like(r)
    cos_theta = np.cos(theta)
    r_sum = np.zeros(len(theta))
    l_sum = np.zeros(len(theta))
    store = {int(j): k for k, j in enumerate(receiver_index)}
    for j in range(int(receiver_index.max()) + 1):
        x = np.minimum(c[j] / c[0] * cos_theta, _X_MAX)
        root = np.sqrt(1 - x * x)
        r_sum += dz * x / root
        l_sum += dz / root
        if j in store:
            r[store[j]] = r_sum
            length[store[j]] = l_sum
    return theta, r, length


def _rays_on_range_grid(theta, r_theta, l_theta, ranges):
    """Invert r(theta) to get theta and L at each horizontal range."""
    order = np.argsort(r_theta, kind="stable")
    r_sorted = r_theta[order]
    # r(theta) falls as theta grows; keep it strictly increasing for interp
    keep = np.concatenate([[True], np.diff(r_sorted) > 0])
    r_sorted = r_sorted[keep]
    th = np.interp(ranges, r_sorted, theta[order][keep])
    ll = np.interp(ranges, r_sorted, l_theta[order][keep])
    return th, ll


def beta_profile(
    depth,
    sound_speed,
    temperature,
    salinity,
    frequencies_hz,
    receiver_depths=None,
    r_max=8000.0,
    dr=1.0,
    absorption_units="dB",
):
    """Depth correction β (dB) for receivers at several depths.

    β(h, f) = -10 log10( 2 Σ_r r dr sin²θ(r) A(r) / L(r)² ), summing over
    horizontal ranges r from the receiver up to `r_max`, where θ is the
    surface grazing angle of the ray reaching the receiver from range r,
    L its path length and A the absorption loss along it (Vagle et al.
    1990, eq. as in beta_h.m).  Adding β to a level measured at depth h
    gives the equivalent level just below the surface.

    Parameters
    ----------
    depth : array
        Depth grid (m), evenly spaced and starting one step below the
        surface (e.g. 1, 2, ..., 1000).
    sound_speed, temperature, salinity : array
        Sound speed (m/s), in-situ temperature (°C) and practical salinity
        on `depth`.
    frequencies_hz : list of float
        Frequencies (Hz).
    receiver_depths : array, optional
        Depths (m) at which to compute β; must be on `depth`.  By default
        every 10 m down to the bottom of `depth`.
    r_max, dr : float, optional
        Largest horizontal range and range step (m), by default 8000 and 1,
        as in beta_h.m.
    absorption_units : {"dB", "Np"}, optional
        How the absorption coefficient is applied.  "dB" (default, correct):
        A = 10**(-alpha L / 10).  "Np" reproduces the original MATLAB code
        (A = exp(-alpha L)), which overestimates absorption; only for
        comparison with published values.

    Returns
    -------
    pd.DataFrame
        Columns "depth" and "beta_<f>Hz" for each frequency, plus
        "alpha_<f>Hz" (dB/m) at each receiver depth.
    """
    depth = np.asarray(depth, dtype=float)
    dz = depth[1] - depth[0]
    if not np.allclose(np.diff(depth), dz) or not np.isclose(depth[0], dz):
        raise ValueError(
            "depth must be evenly spaced, starting one step below the "
            "surface (e.g. 1, 2, ..., 1000)."
        )
    for name, arr in (
        ("sound_speed", sound_speed),
        ("temperature", temperature),
        ("salinity", salinity),
    ):
        if not np.isfinite(np.asarray(arr, dtype=float)).all():
            raise ValueError(f"{name} has missing values on the depth grid.")
    if absorption_units not in ("dB", "Np"):
        raise ValueError("absorption_units must be 'dB' or 'Np'.")
    if receiver_depths is None:
        receiver_depths = np.arange(10, depth[-1] + dz / 2, 10.0)
    receiver_depths = np.asarray(receiver_depths, dtype=float)
    index = np.rint(receiver_depths / dz).astype(int) - 1
    if (index < 0).any() or (index >= len(depth)).any():
        raise ValueError("receiver_depths must lie on the depth grid.")

    theta, r_theta, l_theta = ray_tables(sound_speed, dz, index)
    ranges = np.arange(dr, r_max + dr / 2, dr)
    # Ray geometry does not depend on frequency: invert r(θ) once per depth
    geometry = np.empty((len(index), len(ranges)))  # r dr sin²θ / L²
    path = np.empty((len(index), len(ranges)))  # L
    for k in range(len(index)):
        th, ll = _rays_on_range_grid(theta, r_theta[k], l_theta[k], ranges)
        geometry[k] = ranges * dr * np.sin(th) ** 2 / ll**2
        path[k] = ll
    out = pd.DataFrame({"depth": receiver_depths})
    for f in frequencies_hz:
        alpha = absorption_francois_garrison(
            f / 1000, depth, temperature, salinity
        )[index]
        if absorption_units == "dB":
            loss = 10 ** (-alpha[:, None] * path / 10)
        else:
            loss = np.exp(-alpha[:, None] * path)
        total = np.sum(geometry * loss, axis=1)
        out[f"beta_{int(f)}Hz"] = -10 * np.log10(2 * total)
        out[f"alpha_{int(f)}Hz"] = alpha
    out.attrs = {
        "r_max_m": float(r_max),
        "dr_m": float(dr),
        "absorption_units": absorption_units,
        "frequencies_hz": [float(f) for f in frequencies_hz],
    }
    return out
