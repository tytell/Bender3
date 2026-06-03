"""
Offline simulation of a cantilevered tube specimen for **Simulation mode** (no NI-DAQ).

Uses a solid circular section (25.4 mm OD) and linear elastic cantilever stiffness
``k = 3 E I / L^3`` with a small viscous term ``η \\dot{δ}`` for rate-dependent response.
Young's moduli span the soft-polymer ranges requested for polyurethane vs silicone.
"""
from __future__ import annotations

from typing import Any, Dict, Optional, Tuple

import numpy as np

# Specimen OD (m); user requirement: 25.4 mm tube.
TUBE_OUTER_DIAMETER_M = 25.4e-3

# Mid-range values within the stated bands (MPa → Pa).
MATERIALS: Dict[str, Dict[str, float]] = {
    'polyurethane': {
        'E_Pa': 35.0e6,
        'eta_frac_of_k_ms': 0.012,
    },
    'silicone': {
        'E_Pa': 3.0e6,
        'eta_frac_of_k_ms': 0.035,
    },
}


def _safe_float(value: Any, default: float, *, min_value: Optional[float] = None) -> float:
    """Return finite float with optional lower bound."""
    try:
        out = float(value)
    except (TypeError, ValueError):
        out = float(default)
    if not np.isfinite(out):
        out = float(default)
    if min_value is not None:
        out = max(out, float(min_value))
    return out


def _safe_int(value: Any, default: int, *, min_value: int = 1) -> int:
    """Return bounded int parsed from common numeric-like inputs."""
    try:
        out = int(float(value))
    except (TypeError, ValueError):
        out = int(default)
    if out < int(min_value):
        out = int(min_value)
    return out


def _safe_1d_float_array(values: Any) -> np.ndarray:
    try:
        arr = np.asarray(values, dtype=float).reshape(-1)
    except (TypeError, ValueError):
        return np.zeros(0, dtype=float)
    return arr[np.isfinite(arr)]


def second_moment_solid_circle(diameter_m: float) -> float:
    r = _safe_float(diameter_m, TUBE_OUTER_DIAMETER_M, min_value=1e-6) / 2.0
    return 0.25 * np.pi * r**4


def cantilever_linear_stiffness_N_per_m(E_Pa: float, I_m4: float, L_m: float) -> float:
    L = float(L_m)
    if L <= 0:
        raise ValueError('Cantilever length must be > 0.')
    return 3.0 * float(E_Pa) * float(I_m4) / (L**3)


def _length_m_from_mm(length_mm: Optional[float], fallback_mm: float = 30.0) -> float:
    Lmm = _safe_float(length_mm, fallback_mm, min_value=1.0)
    return Lmm / 1000.0


def specimen_effective_length_m(length_mm: Optional[float]) -> float:
    """Cantilever span (m) for simulation, from morphometrics ``dclamp`` / segment length."""
    return _length_m_from_mm(length_mm)


def tip_displacement_m(angle_deg: np.ndarray, L_m: float) -> np.ndarray:
    """Geometric proxy: δ ≈ L sin|θ| with sign(θ); θ is commanded bend angle (deg)."""
    th = np.deg2rad(_safe_1d_float_array(angle_deg))
    return L_m * np.sin(np.abs(th)) * np.sign(th)


def tip_velocity_m_s(angle_deg: np.ndarray, anglevel_deg_s: np.ndarray, L_m: float) -> np.ndarray:
    th = np.deg2rad(_safe_1d_float_array(angle_deg))
    wd = np.deg2rad(_safe_1d_float_array(anglevel_deg_s))
    n = min(th.size, wd.size)
    th = th[:n]
    wd = wd[:n]
    return L_m * np.cos(np.abs(th)) * wd * np.sign(th)


def resolve_material_key(name: Optional[str]) -> str:
    s = str(name or 'polyurethane').strip().lower()
    if s in MATERIALS:
        return s
    if 'silicone' in s:
        return 'silicone'
    return 'polyurethane'


def simulated_bending_force(
    angle_deg: np.ndarray,
    anglevel_deg_s: np.ndarray,
    *,
    length_mm: Optional[float],
    material_key: str,
    rng: Optional[np.random.Generator] = None,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Return ``(F_N, δ_m, δdot_m_s)`` aligned to the shortest common length of inputs.
    F = k δ + η \\dot{δ} (+ tiny noise).
    """
    key = resolve_material_key(material_key)
    mat = MATERIALS[key]
    L_m = _length_m_from_mm(length_mm)
    I = second_moment_solid_circle(TUBE_OUTER_DIAMETER_M)
    k = cantilever_linear_stiffness_N_per_m(mat['E_Pa'], I, L_m)
    eta = mat['eta_frac_of_k_ms'] * k * 0.001

    ang = _safe_1d_float_array(angle_deg)
    av = _safe_1d_float_array(anglevel_deg_s)
    n = int(min(ang.size, av.size))
    ang = ang[:n]
    av = av[:n]

    delta = tip_displacement_m(ang, L_m)
    deltad = tip_velocity_m_s(ang, av, L_m)
    F = k * delta + eta * deltad
    if rng is not None:
        scale = 0.003 * (np.nanmax(np.abs(F)) + 0.02)
        F = F + rng.normal(0.0, float(scale), size=n)
    return F, delta, deltad


def forcetorque_six_from_bending(F_N: np.ndarray, L_m: float) -> np.ndarray:
    """
    Map scalar bending reaction into a 6-axis wrench consistent with typical FT ordering:
    x,y,z forces (N) and x,y,z torques (N·mm).
    """
    F = np.asarray(F_N, dtype=float).reshape(-1)
    n = F.size
    T_Nm = F * L_m
    Mz_Nmm = T_Nm * 1000.0
    ft = np.zeros((6, n), dtype=float)
    ft[1, :] = F
    ft[2, :] = 0.15 * F
    ft[5, :] = Mz_Nmm
    return ft


def forcetorque_to_raw_voltages(ft_6xn: np.ndarray, calibration_6x6: np.ndarray) -> np.ndarray:
    """
    Invert ``ft = cal.T @ raw`` column-wise (see :meth:`bender_functions.Bender.apply_calibration_forcetorque`).
    """
    cal = np.asarray(calibration_6x6, dtype=float)
    ft = np.asarray(ft_6xn, dtype=float)
    if cal.shape != (6, 6):
        raise ValueError('Calibration must be 6x6 for simulation FT injection.')
    A = cal.T
    if ft.ndim != 2 or ft.shape[0] != 6:
        raise ValueError('ft must have shape (6, n).')
    try:
        return np.linalg.solve(A, ft)
    except np.linalg.LinAlgError:
        raw, _, _, _ = np.linalg.lstsq(A, ft, rcond=None)
        return raw


def force_displacement_series(
    angle_deg: np.ndarray,
    anglevel_deg_s: np.ndarray,
    *,
    length_mm: Optional[float],
    material_key: str,
    max_points: int = 8000,
) -> Tuple[np.ndarray, np.ndarray]:
    """Downsampled δ (mm) and F (N) for plotting."""
    F, delta, _ = simulated_bending_force(
        angle_deg,
        anglevel_deg_s,
        length_mm=length_mm,
        material_key=material_key,
        rng=None,
    )
    max_points = _safe_int(max_points, 8000, min_value=2)
    if F.size > max_points:
        idx = np.linspace(0, F.size - 1, num=max_points, dtype=int)
        return delta[idx] * 1000.0, F[idx]
    return delta * 1000.0, F


def static_stiffness_comparison_delta_grid(
    length_mm: Optional[float],
    delta_max_mm: float = 8.0,
    n_points: int = 120,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Quasi-static curves F = k δ (no rate term) for both materials — same δ grid, different E.
    Returns ``(delta_mm, F_poly, F_sil)``.
    """
    L_m = _length_m_from_mm(length_mm)
    I = second_moment_solid_circle(TUBE_OUTER_DIAMETER_M)
    dmax = _safe_float(delta_max_mm, 8.0, min_value=0.0)
    npts = _safe_int(n_points, 120, min_value=2)
    d_mm = np.linspace(0.0, dmax, npts)
    d_m = d_mm / 1000.0
    k_pu = cantilever_linear_stiffness_N_per_m(MATERIALS['polyurethane']['E_Pa'], I, L_m)
    k_si = cantilever_linear_stiffness_N_per_m(MATERIALS['silicone']['E_Pa'], I, L_m)
    return d_mm, k_pu * d_m, k_si * d_m


def cantilever_stiffness_N_per_m_for_geometry(
    *,
    length_mm: float,
    width_mm: float,
    material_key: str,
) -> Tuple[float, float, float]:
    """
    Linear bending stiffness ``k = 3 E I / L^3`` (N/m) for a **solid circular** section.

    ``width_mm`` is treated as outer diameter :math:`D`; :math:`I = \\pi D^4 / 64` scales as
    :math:`D^4` (stiffness grows sharply with diameter). Stiffness falls as :math:`1/L^3`.
    Returns ``(k_N_per_m, E_Pa, I_m4)``.
    """
    L_m = _safe_float(length_mm, 30.0, min_value=1.0) / 1000.0
    D_m = _safe_float(width_mm, 25.4, min_value=0.5) / 1000.0
    key = resolve_material_key(material_key)
    E = float(MATERIALS[key]['E_Pa'])
    I = float(second_moment_solid_circle(D_m))
    k = float(cantilever_linear_stiffness_N_per_m(E, I, L_m))
    return k, E, I


def force_displacement_comparison_curve(
    *,
    length_mm: float,
    width_mm: float,
    material_key: str,
    delta_max_mm: float = 12.0,
    n_points: int = 160,
    noise_std_N: float = 0.025,
    rng: Optional[np.random.Generator] = None,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Quasi-static ``F \\approx k \\, \\delta`` with Gaussian noise on force (N); ``\\delta`` in mm.

    Uses the same cantilever stiffness as :func:`cantilever_stiffness_N_per_m_for_geometry`.
    """
    k, _, _ = cantilever_stiffness_N_per_m_for_geometry(
        length_mm=length_mm,
        width_mm=width_mm,
        material_key=material_key,
    )
    dmax = _safe_float(delta_max_mm, 12.0, min_value=0.0)
    npts = _safe_int(n_points, 160, min_value=2)
    d_mm = np.linspace(0.0, dmax, npts)
    d_m = d_mm / 1000.0
    F = k * d_m
    gen = rng if rng is not None else np.random.default_rng(0)
    sig = _safe_float(noise_std_N, 0.025, min_value=0.0)
    F = F + gen.normal(0.0, sig, size=F.shape)
    return d_mm, F


def _trapz(y: np.ndarray, x: np.ndarray) -> float:
    if hasattr(np, 'trapezoid'):
        return float(np.trapezoid(y, x))
    return float(np.trapz(y, x))


def viscous_damping_n_s_per_m(material_key: str, k_N_per_m: float) -> float:
    """
    Kelvin–Voigt dashpot **c** (N·s/m) for tip motion.

    Polyurethane uses a **much larger** c than silicone so loops are wide (high damping);
    silicone stays **thin** (low damping). Scales weakly with **k** so geometry changes
    move elastic and viscous terms together.
    """
    key = resolve_material_key(material_key)
    k_r = max(float(k_N_per_m), 1.0)
    scale = np.sqrt(k_r / 8000.0)
    if key == 'polyurethane':
        return 950.0 * scale
    return 120.0 * scale


def oscillatory_viscoelastic_timeseries(
    *,
    length_mm: float,
    width_mm: float,
    material_key: str,
    freq_hz: float,
    amplitude_mm: float,
    n_cycles_shown: float = 4.0,
    points_per_cycle: int = 200,
    noise_std_N: float = 0.02,
    rng: Optional[np.random.Generator] = None,
) -> Dict[str, Any]:
    """
    Sinusoidal tip displacement δ(t)=A sin(ωt) with **Kelvin–Voigt** response::

        F(t) = k \\, \\delta(t) + c \\, \\dot{\\delta}(t)

    plus optional Gaussian noise on **F**. Phase lag between δ and F is
    φ = arctan(c ω / k). Dissipation over the **last full cycle** is
    ∫ F \\dot{δ} \\, dt (J).

    Derived series: root moment M ≈ F·L, small-angle θ ≈ δ/L, curvature κ = M/(EI).
    """
    k, E, I = cantilever_stiffness_N_per_m_for_geometry(
        length_mm=length_mm,
        width_mm=width_mm,
        material_key=material_key,
    )
    f = _safe_float(freq_hz, 1.0, min_value=0.05)
    omega = 2.0 * np.pi * f
    A_m = _safe_float(amplitude_mm, 2.0, min_value=0.01) / 1000.0
    L_m = _safe_float(length_mm, 30.0, min_value=1.0) / 1000.0
    c = viscous_damping_n_s_per_m(material_key, k)

    n_cyc = _safe_float(n_cycles_shown, 4.0, min_value=1.0)
    ppc = _safe_int(points_per_cycle, 200, min_value=20)
    n_tot = max(int(n_cyc * ppc), 80)
    t_end = n_cyc / f
    t = np.linspace(0.0, t_end, n_tot, endpoint=False)

    delta_m = A_m * np.sin(omega * t)
    delta_dot = A_m * omega * np.cos(omega * t)
    F = k * delta_m + c * delta_dot
    sig = _safe_float(noise_std_N, 0.02, min_value=0.0)
    if rng is not None and sig > 0.0:
        F = F + rng.normal(0.0, sig, size=F.shape)

    theta_rad = delta_m / max(L_m, 1e-12)
    M_root = F * L_m
    kappa = M_root / (E * I + 1e-30)

    tan_delta = (c * omega) / k if k > 1e-30 else float('nan')
    phase_lag_rad = float(np.arctan2(c * omega, k)) if k > 1e-30 else 0.0

    T = 2.0 * np.pi / omega
    t_start_last = float(t[-1]) - T
    mask = np.isfinite(t) & np.isfinite(F) & np.isfinite(delta_m) & (t >= t_start_last)
    if np.count_nonzero(mask) < 5:
        W_cycle = float('nan')
    else:
        tp = t[mask]
        Fp = F[mask]
        dm = delta_m[mask]
        if tp.size < 5 or np.nanmax(np.abs(np.diff(tp))) <= 0.0:
            W_cycle = float('nan')
        else:
            dp = np.gradient(dm, tp, edge_order=2)
            W_cycle = _trapz(Fp * dp, tp)

    return {
        't': t,
        'delta_m': delta_m,
        'delta_mm': delta_m * 1000.0,
        'delta_dot_m_s': delta_dot,
        'F': F,
        'theta_rad': theta_rad,
        'theta_deg': np.rad2deg(theta_rad),
        'M': M_root,
        'kappa': kappa,
        'k': float(k),
        'c': float(c),
        'E': float(E),
        'I': float(I),
        'omega': float(omega),
        'tan_delta': float(tan_delta) if np.isfinite(tan_delta) else float('nan'),
        'phase_lag_deg': float(np.rad2deg(phase_lag_rad)),
        'energy_dissipated_last_cycle_J': W_cycle,
    }
