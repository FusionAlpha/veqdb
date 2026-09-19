"""Canonical negative-Z surface geometry and continuous LCFS projection."""

import numpy as np
from numba import njit


@njit(cache=True, nogil=True, inline="always")
def surface_point(p, theta, nc, ns):
    ct, st = np.cos(theta), np.sin(theta)
    before_c, value_c, before_s, value_s = 1.0, ct, 0.0, st
    eta, eta_t, eta_tt = theta + p[4], 1.0, 0.0
    for m in range(1, max(nc - 1, ns) + 1):
        if m < nc:
            c = p[4 + m]
            eta += c * value_c
            eta_t -= m * c * value_s
            eta_tt -= m * m * c * value_c
        if m <= ns:
            c = p[4 + nc + m - 1]
            eta += c * value_s
            eta_t += m * c * value_c
            eta_tt -= m * m * c * value_s
        before_c, value_c = value_c, 2 * ct * value_c - before_c
        before_s, value_s = value_s, 2 * ct * value_s - before_s
    ce, se = np.cos(eta), np.sin(eta)
    R = p[0] + p[2] * ce
    Z = p[1] - p[2] * p[3] * st
    Rt = -p[2] * se * eta_t
    Zt = -p[2] * p[3] * ct
    Rtt = -p[2] * (ce * eta_t * eta_t + se * eta_tt)
    Ztt = p[2] * p[3] * st
    return R, Z, Rt, Zt, Rtt, Ztt, ce, se


@njit(cache=True, nogil=True)
def initial_angles(points, p):
    n = len(points)
    bottom, top = np.argmin(points[:, 1]), np.argmax(points[:, 1])
    out = np.empty(n)
    for i in range(n):
        angle = np.arcsin(min(1.0, max(-1.0, -(points[i, 1] - p[1]) / (p[2] * p[3]))))
        out[i] = np.pi - angle if (i - bottom) % n <= (top - bottom) % n else 2 * np.pi + angle
    return out


@njit(cache=True, nogil=True)
def project(points, p, theta, nc, ns):
    angles = theta.copy()
    error = 0.0
    for i in range(len(points)):
        t = angles[i]
        for iteration in range(12):
            R, Z, Rt, Zt, Rtt, Ztt, _, _ = surface_point(p, t, nc, ns)
            dR, dZ = R - points[i, 0], Z - points[i, 1]
            current = dR * dR + dZ * dZ
            gradient = dR * Rt + dZ * Zt
            tangent2 = Rt * Rt + Zt * Zt
            if not np.isfinite(tangent2) or tangent2 == 0.0:
                angles[i] = np.nan
                return angles, np.inf
            hessian = tangent2 + dR * Rtt + dZ * Ztt
            if hessian <= 0.0:
                hessian = tangent2
            step = -gradient / hessian
            if abs(step) < 1e-10:
                break
            step = min(np.pi / 4, max(-np.pi / 4, step))
            accepted = False
            for backtrack in range(12):
                new = t + step
                (
                    Rn,
                    Zn,
                    _,
                    _,
                    _,
                    _,
                    _,
                    _,
                ) = surface_point(p, new, nc, ns)
                cost = (Rn - points[i, 0]) ** 2 + (Zn - points[i, 1]) ** 2
                if cost < current:
                    t = new
                    accepted = True
                    break
                step *= 0.5
            if not accepted:
                break
        angles[i] = t
        R, Z, _, _, _, _, _, _ = surface_point(p, t, nc, ns)
        error += (R - points[i, 0]) ** 2 + (Z - points[i, 1]) ** 2
    return angles, error


@njit(cache=True, nogil=True)
def projection_diagnostics(points, p, theta, nc, ns):
    """Qualify local closest-point stationarity independently of the outer LM stop."""
    maximum_step = 0.0
    failures = 0
    for i in range(len(points)):
        R, Z, Rt, Zt, Rtt, Ztt, _, _ = surface_point(p, theta[i], nc, ns)
        dR, dZ = R - points[i, 0], Z - points[i, 1]
        tangent2 = Rt * Rt + Zt * Zt
        gradient = dR * Rt + dZ * Zt
        hessian = tangent2 + dR * Rtt + dZ * Ztt
        step = abs(gradient) / tangent2 if tangent2 > 0.0 else np.inf
        maximum_step = max(maximum_step, step)
        if not np.isfinite(step) or step > np.sqrt(np.finfo(np.float64).eps) or hessian <= 0.0:
            failures += 1
    return maximum_step, failures


@njit(cache=True, nogil=True)
def normal_system(points, p, theta, nc, ns):
    matrix = np.empty((len(points), len(p)))
    target = np.empty(len(points))
    for i in range(len(points)):
        R, Z, Rt, Zt, _, _, ce, se = surface_point(p, theta[i], nc, ns)
        norm = np.sqrt(Rt * Rt + Zt * Zt)
        if not np.isfinite(norm) or norm == 0.0:
            matrix[i, :] = np.nan
            target[i] = np.nan
            continue
        nR, nZ = -Zt / norm, Rt / norm
        target[i] = nR * (R - points[i, 0]) + nZ * (Z - points[i, 1])
        matrix[i, 0] = nR
        matrix[i, 1] = nZ
        matrix[i, 2] = nR * ce - nZ * p[3] * np.sin(theta[i])
        matrix[i, 3] = -nZ * p[2] * np.sin(theta[i])
        ct, st = np.cos(theta[i]), np.sin(theta[i])
        matrix[i, 4] = -nR * p[2] * se
        before_c, value_c, before_s, value_s = 1.0, ct, 0.0, st
        for m in range(1, max(nc - 1, ns) + 1):
            if m < nc:
                matrix[i, 4 + m] = -nR * p[2] * se * value_c
            if m <= ns:
                matrix[i, 4 + nc + m - 1] = -nR * p[2] * se * value_s
            before_c, value_c = value_c, 2 * ct * value_c - before_c
            before_s, value_s = value_s, 2 * ct * value_s - before_s
    return matrix, target


@njit(cache=True, nogil=True)
def profile_basis(r, counts, powers):
    values = np.empty((len(r), len(counts), max(counts)))
    for i in range(len(r)):
        x = 2 * r[i] ** 2 - 1
        T = np.empty(max(counts))
        T[0] = 1.0
        if len(T) > 1:
            T[1] = x
        for j in range(2, len(T)):
            T[j] = 2 * x * T[j - 1] - T[j - 2]
        for k in range(len(counts)):
            for j in range(len(T)):
                values[i, k, j] = r[i] ** powers[k] * (1 - r[i] ** 2) * T[j]
    return values


@njit(cache=True, nogil=True)
def local_profiles(boundary, coefficients, r, values, powers, nc, ns):
    out = np.empty((len(r), 4 + nc + ns))
    for i in range(len(r)):
        p = np.empty(3 + nc + ns)
        for k in range(len(p)):
            edge = 0.0
            if k == 2:
                edge = boundary[3]
            elif k >= 3:
                edge = boundary[k + 1]
            p[k] = r[i] ** powers[k] * edge
            for j in range(coefficients.shape[1]):
                p[k] += values[i, k, j] * coefficients[k, j]
        out[i, 0] = boundary[0] + boundary[2] * p[0]
        out[i, 1] = boundary[1] + boundary[2] * p[1]
        out[i, 2] = boundary[2] * r[i]
        out[i, 3] = p[2]
        out[i, 4:] = p[3:]
    return out


@njit(cache=True, nogil=True)
def radial_derivatives(boundary, coefficients, r, powers, nc, ns):
    out = np.empty((len(r), 4 + nc + ns))
    for i in range(len(r)):
        x = 2 * r[i] ** 2 - 1
        T = np.empty(coefficients.shape[1])
        D = np.empty_like(T)
        T[0] = 1.0
        D[0] = 0.0
        if len(T) > 1:
            T[1] = x
            D[1] = 1.0
        for degree in range(2, len(T)):
            T[degree] = 2 * x * T[degree - 1] - T[degree - 2]
            D[degree] = 2 * T[degree - 1] + 2 * x * D[degree - 1] - D[degree - 2]
        p = np.empty(3 + nc + ns)
        for k in range(len(p)):
            rp = r[i] ** powers[k]
            drp = powers[k] * r[i] ** (powers[k] - 1) if powers[k] > 0 else 0.0
            edge = 0.0
            if k == 2:
                edge = boundary[3]
            elif k >= 3:
                edge = boundary[k + 1]
            p[k] = drp * edge
            for degree in range(len(T)):
                p[k] += coefficients[k, degree] * (
                    (drp * (1 - r[i] ** 2) - 2 * r[i] * rp) * T[degree]
                    + rp * (1 - r[i] ** 2) * D[degree] * 4 * r[i]
                )
        out[i, 0] = boundary[2] * p[0]
        out[i, 1] = boundary[2] * p[1]
        out[i, 2] = boundary[2]
        out[i, 3] = p[2]
        out[i, 4:] = p[3:]
    return out
