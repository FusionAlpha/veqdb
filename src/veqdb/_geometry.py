"""Fusionprime diagnostic geometry kernels, lifted verbatim (88d0fa0).

Source: ``fusionprime.base.model.geometry`` and
``fusionprime.base.model.states.equilibrium``.  Every function body below is
byte-identical to the source module; only the module wrapper (docstring,
imports, ``MU0`` constant) is codec-side.  The compact writer's diagnostic
chain (``veqdb._diagnostics``) evaluates these kernels in the same
order as the fusionprime ``Equilibrium`` properties so that both
implementations emit bitwise-identical records for the same fit.

Provenance audit: see ``output/geqdsk_d1_audit_20260918/`` for the D1
verbatim methodology used by every lifted ``fitting`` module.
"""

from __future__ import annotations

import numpy as np
from numba import njit

MU0 = 4.0e-7 * np.pi

# ---- model/geometry.py (88d0fa0) ----

# :15
@njit(cache=True, nogil=True)
def _update_chebyshev_fields(out: np.ndarray, r: np.ndarray) -> None:
    """Fill Chebyshev values and two radial derivatives without temporaries."""

    orders = out.shape[1]
    for i in range(r.shape[0]):
        x = 2.0 * r[i] * r[i] - 1.0
        dx_dr = 4.0 * r[i]
        out[0, 0, i] = 1.0
        out[1, 0, i] = 0.0
        out[2, 0, i] = 0.0
        if orders == 1:
            continue
        value_previous = 1.0
        value = x
        dx_previous = 0.0
        dx = 1.0
        dxx_previous = 0.0
        dxx = 0.0
        out[0, 1, i] = value
        out[1, 1, i] = dx * dx_dr
        out[2, 1, i] = dxx * dx_dr * dx_dr + 4.0 * dx
        for order in range(1, orders - 1):
            value_next = 2.0 * x * value - value_previous
            dx_next = 2.0 * value + 2.0 * x * dx - dx_previous
            dxx_next = 4.0 * dx + 2.0 * x * dxx - dxx_previous
            out[0, order + 1, i] = value_next
            out[1, order + 1, i] = dx_next * dx_dr
            out[2, order + 1, i] = dxx_next * dx_dr * dx_dr + 4.0 * dx_next
            value_previous, value = value, value_next
            dx_previous, dx = dx, dx_next
            dxx_previous, dxx = dxx, dxx_next

# :49
@njit(cache=True, nogil=True)
def _evaluate_profile_at_axis(coefficients: np.ndarray) -> float:
    """Evaluate a Chebyshev coefficient series at ``x=2*0**2-1=-1``."""

    value = 0.0
    sign = 1.0
    for coefficient in coefficients:
        value += sign * coefficient
        sign = -sign
    return value

# :61
@njit(cache=True, nogil=True)
def _update_fourier_fields(out: np.ndarray, theta: np.ndarray) -> None:
    """Fill six Fourier value/derivative tables into one allocated buffer."""

    for order in range(out.shape[1]):
        harmonic = float(order)
        harmonic2 = harmonic * harmonic
        for j in range(theta.shape[0]):
            cosine = np.cos(harmonic * theta[j])
            sine = np.sin(harmonic * theta[j])
            out[0, order, j] = cosine
            out[1, order, j] = sine
            out[2, order, j] = harmonic * cosine
            out[3, order, j] = harmonic * sine
            out[4, order, j] = harmonic2 * cosine
            out[5, order, j] = harmonic2 * sine

# :79
@njit(cache=True, nogil=True, inline="always")
def _radial_power_terms(value: float, exponent: int) -> tuple[float, float, float]:
    """Return a radial power and its first two derivatives at one node."""

    if exponent == 0:
        return 1.0, 0.0, 0.0
    radial = value**exponent
    radial_r = exponent * value ** (exponent - 1)
    radial_rr = 0.0
    if exponent > 1:
        radial_rr = exponent * (exponent - 1) * value ** (exponent - 2)
    return radial, radial_r, radial_rr

# :93
@njit(cache=True, nogil=True)
def _update_basic_profile(
    out: np.ndarray,
    coefficients: np.ndarray,
    T: np.ndarray,
    T_r: np.ndarray,
    T_rr: np.ndarray,
    r: np.ndarray,
    lcfs: float,
) -> None:
    """Fill one enveloped profile through second order without work arrays."""

    for i in range(r.shape[0]):
        series = 0.0
        series_r = 0.0
        series_rr = 0.0
        for order in range(coefficients.shape[0]):
            coefficient = coefficients[order]
            series += coefficient * T[order, i]
            series_r += coefficient * T_r[order, i]
            series_rr += coefficient * T_rr[order, i]
        envelope = 1.0 - r[i] * r[i]
        envelope_r = -2.0 * r[i]
        out[0, i] = lcfs + envelope * series
        out[1, i] = envelope_r * series + envelope * series_r
        out[2, i] = -2.0 * series + 2.0 * envelope_r * series_r + envelope * series_rr

# :121
@njit(cache=True, nogil=True)
def _update_shaping_fields(
    out: np.ndarray,
    coefficients: np.ndarray,
    lcfs: np.ndarray,
    T: np.ndarray,
    T_r: np.ndarray,
    T_rr: np.ndarray,
    r: np.ndarray,
    powers: np.ndarray,
    first_harmonic: int,
) -> None:
    """Fill all cosine or sine shaping profiles without intermediate arrays."""

    out.fill(0.0)
    for coefficient_row in range(coefficients.shape[0]):
        harmonic = coefficient_row + first_harmonic
        exponent = int(powers[harmonic])
        for i in range(r.shape[0]):
            series = 0.0
            series_r = 0.0
            series_rr = 0.0
            for radial_order in range(coefficients.shape[1]):
                coefficient = coefficients[coefficient_row, radial_order]
                series += coefficient * T[radial_order, i]
                series_r += coefficient * T_r[radial_order, i]
                series_rr += coefficient * T_rr[radial_order, i]
            envelope = 1.0 - r[i] * r[i]
            amplitude = lcfs[coefficient_row] + envelope * series
            amplitude_r = -2.0 * r[i] * series + envelope * series_r
            amplitude_rr = -2.0 * series - 4.0 * r[i] * series_r + envelope * series_rr
            radial, radial_r, radial_rr = _radial_power_terms(r[i], exponent)
            out[harmonic, 0, i] = radial * amplitude
            out[harmonic, 1, i] = radial_r * amplitude + radial * amplitude_r
            out[harmonic, 2, i] = (
                radial_rr * amplitude + 2.0 * radial_r * amplitude_r + radial * amplitude_rr
            )

# :160
@njit(cache=True, nogil=True)
def _update_eta(
    out: np.ndarray,
    theta: np.ndarray,
    c: np.ndarray,
    s: np.ndarray,
    cos_mtheta: np.ndarray,
    sin_mtheta: np.ndarray,
) -> None:
    """Fill phase values directly from cached shaping and Fourier tables."""

    for i in range(out.shape[0]):
        for j in range(out.shape[1]):
            value = theta[j] + c[0, 0, i]
            for harmonic in range(1, c.shape[0]):
                value += c[harmonic, 0, i] * cos_mtheta[harmonic, j]
                value += s[harmonic, 0, i] * sin_mtheta[harmonic, j]
            out[i, j] = value

# :180
@njit(cache=True, nogil=True)
def _update_eta_gradient(
    out: np.ndarray,
    c: np.ndarray,
    s: np.ndarray,
    cos_mtheta: np.ndarray,
    sin_mtheta: np.ndarray,
    m_cos_mtheta: np.ndarray,
    m_sin_mtheta: np.ndarray,
) -> None:
    """Fill first radial and poloidal phase derivatives."""

    for i in range(out.shape[1]):
        for j in range(out.shape[2]):
            radial = c[0, 1, i]
            poloidal = 1.0
            for harmonic in range(1, c.shape[0]):
                radial += c[harmonic, 1, i] * cos_mtheta[harmonic, j]
                radial += s[harmonic, 1, i] * sin_mtheta[harmonic, j]
                poloidal -= c[harmonic, 0, i] * m_sin_mtheta[harmonic, j]
                poloidal += s[harmonic, 0, i] * m_cos_mtheta[harmonic, j]
            out[0, i, j] = radial
            out[1, i, j] = poloidal

# :205
@njit(cache=True, nogil=True)
def _update_eta_hessian(
    out: np.ndarray,
    c: np.ndarray,
    s: np.ndarray,
    cos_mtheta: np.ndarray,
    sin_mtheta: np.ndarray,
    m_cos_mtheta: np.ndarray,
    m_sin_mtheta: np.ndarray,
    m2_cos_mtheta: np.ndarray,
    m2_sin_mtheta: np.ndarray,
) -> None:
    """Fill the three independent second derivatives of phase."""

    for i in range(out.shape[1]):
        for j in range(out.shape[2]):
            radial2 = c[0, 2, i]
            mixed = 0.0
            poloidal2 = 0.0
            for harmonic in range(1, c.shape[0]):
                radial2 += c[harmonic, 2, i] * cos_mtheta[harmonic, j]
                radial2 += s[harmonic, 2, i] * sin_mtheta[harmonic, j]
                mixed -= c[harmonic, 1, i] * m_sin_mtheta[harmonic, j]
                mixed += s[harmonic, 1, i] * m_cos_mtheta[harmonic, j]
                poloidal2 -= c[harmonic, 0, i] * m2_cos_mtheta[harmonic, j]
                poloidal2 -= s[harmonic, 0, i] * m2_sin_mtheta[harmonic, j]
            out[0, i, j] = radial2
            out[1, i, j] = mixed
            out[2, i, j] = poloidal2

# :236
@njit(cache=True, nogil=True)
def _update_trigonometry(out: np.ndarray, phase: np.ndarray) -> None:
    """Fill cosine and sine rows for a one- or two-dimensional phase array."""

    flat_out = out.reshape(2, phase.size)
    flat_phase = phase.reshape(phase.size)
    for index in range(phase.size):
        flat_out[0, index] = np.cos(flat_phase[index])
        flat_out[1, index] = np.sin(flat_phase[index])

# :247
@njit(cache=True, nogil=True)
def _update_R(
    out: np.ndarray,
    R0: float,
    a: float,
    r: np.ndarray,
    h: np.ndarray,
    cosine_eta: np.ndarray,
) -> None:
    """Fill major-radius coordinates without broadcast temporaries."""

    for i in range(out.shape[0]):
        for j in range(out.shape[1]):
            out[i, j] = R0 + a * (h[i] + r[i] * cosine_eta[i, j])

# :263
@njit(cache=True, nogil=True)
def _update_R_gradient(
    out: np.ndarray,
    a: float,
    r: np.ndarray,
    h_r: np.ndarray,
    cosine_eta: np.ndarray,
    sine_eta: np.ndarray,
    eta_r: np.ndarray,
    eta_t: np.ndarray,
) -> None:
    """Fill the two major-radius first derivatives."""

    for i in range(out.shape[1]):
        for j in range(out.shape[2]):
            out[0, i, j] = a * (h_r[i] + cosine_eta[i, j] - r[i] * sine_eta[i, j] * eta_r[i, j])
            out[1, i, j] = -a * r[i] * sine_eta[i, j] * eta_t[i, j]

# :282
@njit(cache=True, nogil=True)
def _update_R_hessian(
    out: np.ndarray,
    a: float,
    r: np.ndarray,
    h_rr: np.ndarray,
    cosine_eta: np.ndarray,
    sine_eta: np.ndarray,
    eta_r: np.ndarray,
    eta_t: np.ndarray,
    eta_rr: np.ndarray,
    eta_rt: np.ndarray,
    eta_tt: np.ndarray,
) -> None:
    """Fill the three independent major-radius second derivatives."""

    for i in range(out.shape[1]):
        for j in range(out.shape[2]):
            cosine = cosine_eta[i, j]
            sine = sine_eta[i, j]
            radial_phase = eta_r[i, j]
            poloidal_phase = eta_t[i, j]
            out[0, i, j] = a * (
                h_rr[i]
                - 2.0 * sine * radial_phase
                - r[i] * (cosine * radial_phase * radial_phase + sine * eta_rr[i, j])
            )
            out[1, i, j] = -a * (
                sine * poloidal_phase
                + r[i] * (cosine * radial_phase * poloidal_phase + sine * eta_rt[i, j])
            )
            out[2, i, j] = (
                -a * r[i] * (cosine * poloidal_phase * poloidal_phase + sine * eta_tt[i, j])
            )

# :318
@njit(cache=True, nogil=True)
def _update_Z(
    out: np.ndarray,
    Z0: float,
    a: float,
    r: np.ndarray,
    v: np.ndarray,
    kappa: np.ndarray,
    sine_theta: np.ndarray,
) -> None:
    """Fill vertical coordinates on the tensor grid."""

    for i in range(out.shape[0]):
        for j in range(out.shape[1]):
            out[i, j] = Z0 + a * (v[i] - r[i] * kappa[i] * sine_theta[j])

# :335
@njit(cache=True, nogil=True)
def _update_Z_gradient(
    out: np.ndarray,
    a: float,
    r: np.ndarray,
    v_r: np.ndarray,
    kappa: np.ndarray,
    sine_theta: np.ndarray,
    cosine_theta: np.ndarray,
) -> None:
    """Fill the two vertical-coordinate first derivatives."""

    for i in range(out.shape[1]):
        radial_factor = kappa[0, i] + r[i] * kappa[1, i]
        for j in range(out.shape[2]):
            out[0, i, j] = a * (v_r[i] - radial_factor * sine_theta[j])
            out[1, i, j] = -a * r[i] * kappa[0, i] * cosine_theta[j]

# :354
@njit(cache=True, nogil=True)
def _update_Z_hessian(
    out: np.ndarray,
    a: float,
    r: np.ndarray,
    v_rr: np.ndarray,
    kappa: np.ndarray,
    sine_theta: np.ndarray,
    cosine_theta: np.ndarray,
) -> None:
    """Fill the three independent vertical-coordinate second derivatives."""

    for i in range(out.shape[1]):
        radial_factor = kappa[0, i] + r[i] * kappa[1, i]
        radial2_factor = 2.0 * kappa[1, i] + r[i] * kappa[2, i]
        for j in range(out.shape[2]):
            out[0, i, j] = a * (v_rr[i] - radial2_factor * sine_theta[j])
            out[1, i, j] = -a * radial_factor * cosine_theta[j]
            out[2, i, j] = a * r[i] * kappa[0, i] * sine_theta[j]

# :375
@njit(cache=True, nogil=True)
def _update_jacobian(
    out: np.ndarray,
    R_r: np.ndarray,
    R_t: np.ndarray,
    Z_r: np.ndarray,
    Z_t: np.ndarray,
) -> None:
    """Fill the oriented coordinate Jacobian."""

    for i in range(out.shape[0]):
        for j in range(out.shape[1]):
            out[i, j] = R_t[i, j] * Z_r[i, j] - R_r[i, j] * Z_t[i, j]

# :390
@njit(cache=True, nogil=True)
def _update_covariant_metric(
    out: np.ndarray,
    R_r: np.ndarray,
    R_t: np.ndarray,
    Z_r: np.ndarray,
    Z_t: np.ndarray,
) -> None:
    """Fill all three independent covariant metric components."""

    for i in range(out.shape[1]):
        for j in range(out.shape[2]):
            rr = R_r[i, j]
            rt = R_t[i, j]
            zr = Z_r[i, j]
            zt = Z_t[i, j]
            out[0, i, j] = rr * rr + zr * zr
            out[1, i, j] = rr * rt + zr * zt
            out[2, i, j] = rt * rt + zt * zt

# :427
@njit(cache=True, nogil=True, inline="always")
def _axis_linear_limit(value_1: float, value_2: float, r: np.ndarray) -> float:
    """Extrapolate a finite local field linearly to an axis grid point."""

    denominator = r[2] - r[1]
    if denominator == 0.0 or not np.isfinite(value_1) or not np.isfinite(value_2):
        return np.nan
    return value_1 + (value_2 - value_1) / denominator * (r[0] - r[1])

# :437
@njit(cache=True, nogil=True, inline="always")
def _axis_leading_radial_coefficient(
    value_1: float,
    value_2: float,
    r: np.ndarray,
) -> float:
    """Recover ``A`` from a regular field with leading behavior ``r*A``."""

    if r[1] == 0.0 or r[2] == 0.0:
        return np.nan
    return _axis_linear_limit(value_1 / r[1], value_2 / r[2], r)

# :450
@njit(cache=True, nogil=True)
def _update_gs_surface_fields(
    out: np.ndarray,
    r: np.ndarray,
    R: np.ndarray,
    R_r: np.ndarray,
    R_t: np.ndarray,
    R_rr: np.ndarray,
    R_rt: np.ndarray,
    R_tt: np.ndarray,
    Z_r: np.ndarray,
    Z_t: np.ndarray,
    Z_rr: np.ndarray,
    Z_rt: np.ndarray,
    Z_tt: np.ndarray,
) -> int:
    """Fill the four surface metric combinations used by the GS equation."""

    nr, nt = R.shape
    has_axis = nr >= 3 and abs(r[0]) < 1.0e-10
    for i in range(nr):
        for j in range(nt):
            radius = R[i, j]
            radius_r = R_r[i, j]
            radius_theta = R_t[i, j]
            radius_rr = R_rr[i, j]
            radius_rtheta = R_rt[i, j]
            radius_thetatheta = R_tt[i, j]
            vertical_r = Z_r[i, j]
            vertical_theta = Z_t[i, j]
            vertical_rr = Z_rr[i, j]
            vertical_rtheta = Z_rt[i, j]
            vertical_thetatheta = Z_tt[i, j]

            jacobian = radius_theta * vertical_r - radius_r * vertical_theta
            jacobian_r = (
                radius_rtheta * vertical_r
                + radius_theta * vertical_rr
                - radius_rr * vertical_theta
                - radius_r * vertical_rtheta
            )
            jacobian_theta = (
                radius_thetatheta * vertical_r
                + radius_theta * vertical_rtheta
                - radius_rtheta * vertical_theta
                - radius_r * vertical_thetatheta
            )
            jacobian_radius = jacobian * radius
            jacobian_radius_r = jacobian_r * radius + jacobian * radius_r
            jacobian_radius_theta = jacobian_theta * radius + jacobian * radius_theta
            grt = radius_r * radius_theta + vertical_r * vertical_theta
            grt_theta = (
                radius_rtheta * radius_theta
                + radius_r * radius_thetatheta
                + vertical_rtheta * vertical_theta
                + vertical_r * vertical_thetatheta
            )
            gtt = radius_theta * radius_theta + vertical_theta * vertical_theta
            gtt_r = 2.0 * (radius_theta * radius_rtheta + vertical_theta * vertical_rtheta)

            out[0, i, j] = jacobian / radius
            if has_axis and i == 0:
                # J remains exactly zero on an axis-owning coordinate row.  The
                # three ratios have removable limits reconstructed below from
                # resolved off-axis rows; the raw Jacobian is never clamped.
                out[1, i, j] = np.nan
                out[2, i, j] = np.nan
                out[3, i, j] = np.nan
                continue
            if jacobian_radius == 0.0 or not np.isfinite(jacobian_radius):
                return i + 1
            inverse = 1.0 / jacobian_radius
            out[1, i, j] = gtt * inverse
            out[2, i, j] = gtt_r * inverse - gtt * jacobian_radius_r * inverse * inverse
            out[3, i, j] = (grt_theta - grt * jacobian_radius_theta * inverse) * inverse

    if has_axis:
        for j in range(nt):
            leading = _axis_leading_radial_coefficient(out[1, 1, j], out[1, 2, j], r)
            out[1, 0, j] = r[0] * leading
            out[2, 0, j] = leading
            out[3, 0, j] = _axis_linear_limit(out[3, 1, j], out[3, 2, j], r)

    for field in range(out.shape[0]):
        for i in range(nr):
            for j in range(nt):
                if not np.isfinite(out[field, i, j]):
                    return i + 1
    return 0

# :541
@njit(cache=True, nogil=True)
def _update_gs_radial_fields(
    out: np.ndarray,
    gttdivJR: np.ndarray,
    gttdivJR_r: np.ndarray,
    JdivR: np.ndarray,
) -> None:
    """Take periodic means of the GS surface metric combinations."""

    scale = 1.0 / gttdivJR.shape[1]
    for i in range(gttdivJR.shape[0]):
        sum_gtt = 0.0
        sum_gtt_r = 0.0
        sum_jdivr = 0.0
        for j in range(gttdivJR.shape[1]):
            sum_gtt += gttdivJR[i, j]
            sum_gtt_r += gttdivJR_r[i, j]
            sum_jdivr += JdivR[i, j]
        out[0, i] = scale * sum_gtt
        out[1, i] = scale * sum_gtt_r
        out[2, i] = scale * sum_jdivr

# :596
@njit(cache=True, nogil=True)
def _update_poloidal_integral(
    out: np.ndarray,
    left: np.ndarray,
    right: np.ndarray,
    scale: float,
) -> None:
    """Fill scaled row reductions of one elementwise product."""

    for i in range(out.shape[0]):
        total = 0.0
        for j in range(left.shape[1]):
            total += left[i, j] * right[i, j]
        out[i] = scale * total

# :612
@njit(cache=True, nogil=True)
def _update_volume(
    out: np.ndarray,
    radius: np.ndarray,
    Z_t: np.ndarray,
    scale: float,
) -> None:
    """Fill toroidal volume from the surface line integral."""

    for i in range(out.shape[0]):
        total = 0.0
        for j in range(radius.shape[1]):
            total += radius[i, j] * radius[i, j] * Z_t[i, j]
        out[i] = scale * total

# :663
def _materialize_basic_profile(
    coefficients: np.ndarray,
    T: np.ndarray,
    T_r: np.ndarray,
    T_rr: np.ndarray,
    r: np.ndarray,
    *,
    lcfs: float = 0.0,
) -> np.ndarray:
    """Allocate one cache value and delegate its complete fill to Numba."""

    out = np.empty((3, r.size), dtype=np.float64)
    _update_basic_profile(out, coefficients, T, T_r, T_rr, r, lcfs)
    return out

# :707
def _uniform_theta(count: int) -> np.ndarray:
    """Construct the canonical uniform periodic poloidal axis."""

    return np.linspace(0.0, 2.0 * np.pi, count, endpoint=False, dtype=np.float64)

# ---- model/states/equilibrium.py (88d0fa0) ----

# :1887
@njit(cache=True, nogil=True)
def _update_matvec(out: np.ndarray, matrix: np.ndarray, values: np.ndarray) -> None:
    """Fill a dense matrix-vector product without a temporary result."""

    for i in range(out.shape[0]):
        total = 0.0
        for k in range(values.shape[0]):
            total += matrix[i, k] * values[k]
        out[i] = total

# :1898
@njit(cache=True, nogil=True)
def _update_scaled_profile(out: np.ndarray, values: np.ndarray, scale: float) -> None:
    for i in range(out.shape[0]):
        out[i] = scale * values[i]

# :1904
@njit(cache=True, nogil=True)
def _update_scaled_product(
    out: np.ndarray,
    left: np.ndarray,
    right: np.ndarray,
    scale: float,
) -> None:
    for i in range(out.shape[0]):
        out[i] = scale * left[i] * right[i]

# :1928
@njit(cache=True, nogil=True)
def _update_edge_anchored_integral(
    out: np.ndarray,
    derivative: np.ndarray,
    accumulator: np.ndarray,
    weights: np.ndarray,
    edge_value: float,
    derivative_scale: float,
) -> None:
    """Integrate a radial derivative while fixing the LCFS value."""

    edge_integral = 0.0
    for k in range(derivative.shape[0]):
        edge_integral += weights[k] * derivative[k]
    for i in range(out.shape[0]):
        prefix = 0.0
        for k in range(derivative.shape[0]):
            prefix += accumulator[i, k] * derivative[k]
        out[i] = edge_value + derivative_scale * (prefix - edge_integral)

# :1949
@njit(cache=True, nogil=True)
def _update_signed_sqrt(out: np.ndarray, values: np.ndarray, edge_value: float) -> None:
    sign = -1.0 if edge_value < 0.0 else 1.0
    for i in range(out.shape[0]):
        out[i] = sign * np.sqrt(values[i])

# :1956
@njit(cache=True, nogil=True)
def _update_G_fields(
    out: np.ndarray,
    R: np.ndarray,
    JdivR: np.ndarray,
    gttdivJR: np.ndarray,
    gttdivJR_r: np.ndarray,
    grtdivJR_t: np.ndarray,
    FF_psi: np.ndarray,
    P_psi: np.ndarray,
    psi_r: np.ndarray,
    psi_rr: np.ndarray,
) -> int:
    """Fill physical source, geometry, and total local GS defect fields."""

    for i in range(R.shape[0]):
        for j in range(R.shape[1]):
            radius = R[i, j]
            G1 = JdivR[i, j] * (FF_psi[i] + MU0 * radius * radius * P_psi[i])
            G2 = gttdivJR[i, j] * psi_rr[i] + (gttdivJR_r[i, j] - grtdivJR_t[i, j]) * psi_r[i]
            total = G1 + G2
            if not np.isfinite(G1) or not np.isfinite(G2) or not np.isfinite(total):
                return i + 1
            out[0, i, j] = G1
            out[1, i, j] = G2
            out[2, i, j] = total
    return 0

# :1985
@njit(cache=True, nogil=True)
def _integrate_G1(G1: np.ndarray, weights: np.ndarray) -> float:
    """Integrate the source contribution over normalized radius and theta."""

    total = 0.0
    theta_scale = 2.0 * np.pi / G1.shape[1]
    for i in range(G1.shape[0]):
        row = 0.0
        for j in range(G1.shape[1]):
            row += G1[i, j]
        total += weights[i] * row
    return theta_scale * total

# :1999
@njit(cache=True, nogil=True, inline="always")
def _axis_even_limit(values: np.ndarray, r: np.ndarray) -> None:
    """Repair one removable even-in-r axis sample from two off-axis rows."""

    if values.shape[0] < 3 or abs(r[0]) >= 1.0e-10:
        return
    x1 = r[1] * r[1]
    x2 = r[2] * r[2]
    if x2 == x1 or not np.isfinite(values[1]) or not np.isfinite(values[2]):
        return
    values[0] = values[1] + (values[2] - values[1]) / (x2 - x1) * (r[0] * r[0] - x1)

# :2012
@njit(cache=True, nogil=True)
def _update_q(
    out: np.ndarray,
    F: np.ndarray,
    Ln_r: np.ndarray,
    psi_r: np.ndarray,
    r: np.ndarray,
) -> int:
    for i in range(out.shape[0]):
        out[i] = np.nan if psi_r[i] == 0.0 else F[i] * Ln_r[i] / psi_r[i]
    _axis_even_limit(out, r)
    for i in range(out.shape[0]):
        if not np.isfinite(out[i]):
            return i + 1
    return 0

# :2029
@njit(cache=True, nogil=True)
def _update_betat(
    pressure: np.ndarray,
    V_r: np.ndarray,
    weights: np.ndarray,
    B0: float,
) -> float:
    numerator = 0.0
    denominator = 0.0
    for i in range(pressure.shape[0]):
        weighted_volume = weights[i] * V_r[i]
        numerator += weighted_volume * pressure[i]
        denominator += weighted_volume
    return 2.0 * MU0 * numerator / (denominator * B0 * B0)
