#!/usr/bin/env python3
"""
Dandelion Substrate & Mori-Zwanzig Active DNA Viscoelastic Verification (v1.0.1)
100% Pure Python 3 Standard Library Port (--seed 42, zero numpy/scipy).
Author: Justin D. Grimm (ORCID: 0009-0000-5299-0243)
Affiliation: Independent Dandelion Substrate Research Group, Gurley, Alabama, USA
Reference: Levanon et al., Nature Physics (2026), doi:10.1038/s41567-026-03505-y
"""
import argparse
import cmath
import json
import math
import random
import time

PC_FCC = 0.1201635
FOWLER_NEUTRAL_BYTE = 121
KINTSUGI_WELD_BYTE = 251

ETA_S = 1.0
ELL_V = 0.12
ZETA_DEFAULT = 0.340
LAM_ALIGN = 1.20
GAMMA_Q = 0.065
D_Q = 0.025


def eta_eff(k: float) -> float:
    return ETA_S * (1.0 + (ELL_V * k) ** 2)


def dna_rheology(lc_um: float) -> tuple[float, float, float, float]:
    lc_ref = 16.5
    p_bond = 1.0 - math.exp(-0.048 * (lc_um / 10.0) ** 1.8)
    p_inf = 0.0 if p_bond <= PC_FCC else ((p_bond - PC_FCC) / (1.0 - PC_FCC)) ** 0.41
    tau_p = 2.0 * ((lc_um / lc_ref) ** 2.15)
    g_p = 0.115 * p_inf * ((lc_um / lc_ref) ** 0.55)
    return p_bond, p_inf, g_p, tau_p


def mode_stability(k: float, g_p: float, tau_p: float, zeta_val: float = ZETA_DEFAULT) -> dict[str, float]:
    visc = eta_eff(k)
    big_gamma_q = GAMMA_Q + D_Q * (k ** 2)
    gamma_act = (LAM_ALIGN * zeta_val) / visc - big_gamma_q
    gamma_ve = (1.0 / tau_p) + (g_p / visc)

    tr = gamma_act - gamma_ve
    det_matrix = gamma_act * (-gamma_ve) - (-LAM_ALIGN / visc) * (g_p * zeta_val / visc)
    det_exact = big_gamma_q * (g_p / visc) - gamma_act / tau_p
    assert abs(det_matrix - det_exact) < 1e-12, "Exact algebraic determinant cancellation failed"

    disc = tr ** 2 - 4.0 * det_exact
    omega_lin = 0.5 * math.sqrt(max(0.0, -disc)) if det_exact > 0.0 else 0.0
    de_nem = tau_p * big_gamma_q
    de_crit = 1.0 + visc / (g_p * tau_p) if g_p > 0.0 else float("inf")
    rho_act = gamma_act / (gamma_act + gamma_ve) if gamma_act > 0.0 else 0.0
    return {
        "tr": tr,
        "det": det_exact,
        "disc": disc,
        "omega": omega_lin,
        "de_nem": de_nem,
        "de_crit": de_crit,
        "rho_act": rho_act,
    }


def simulate_limit_cycle(
    lc_um: float = 56.4,
    l_box_mm: float = 2.0,
    zeta_val: float = ZETA_DEFAULT,
    steps: int = 35000,
    dt: float = 0.05,
) -> tuple[float, float, float, float, float]:
    _, _, g_p, tau_p = dna_rheology(lc_um)
    k = math.pi / l_box_mm
    visc = eta_eff(k)
    big_gamma_q = GAMMA_Q + D_Q * (k ** 2)
    sigma_max = 0.12
    inv_smax2 = 1.0 / (sigma_max * sigma_max)
    inv_visc = 1.0 / visc
    inv_taup = 1.0 / tau_p
    inv_k = 1.0 / k

    th, sg = 1e-3, 0.0
    half = steps // 2
    v_ss = []
    t_ss = []
    th_ss = []
    sg_ss = []

    for i in range(steps):
        # RK4 step inlined for speed in pure Python stdlib
        f1 = 1.0 + sg * sg * inv_smax2
        ku1 = (zeta_val * th - sg) * inv_visc
        dth1 = -big_gamma_q * th + LAM_ALIGN * ku1 - 0.35 * (th * th * th)
        dsg1 = -(f1 * inv_taup) * sg + g_p * ku1
        if i >= half:
            t_ss.append(i * dt)
            v_ss.append(ku1 * inv_k)
            th_ss.append(th)
            sg_ss.append(sg)

        th2 = th + 0.5 * dt * dth1
        sg2 = sg + 0.5 * dt * dsg1
        f2 = 1.0 + sg2 * sg2 * inv_smax2
        ku2 = (zeta_val * th2 - sg2) * inv_visc
        dth2 = -big_gamma_q * th2 + LAM_ALIGN * ku2 - 0.35 * (th2 * th2 * th2)
        dsg2 = -(f2 * inv_taup) * sg2 + g_p * ku2

        th3 = th + 0.5 * dt * dth2
        sg3 = sg + 0.5 * dt * dsg2
        f3 = 1.0 + sg3 * sg3 * inv_smax2
        ku3 = (zeta_val * th3 - sg3) * inv_visc
        dth3 = -big_gamma_q * th3 + LAM_ALIGN * ku3 - 0.35 * (th3 * th3 * th3)
        dsg3 = -(f3 * inv_taup) * sg3 + g_p * ku3

        th4 = th + dt * dth3
        sg4 = sg + dt * dsg3
        f4 = 1.0 + sg4 * sg4 * inv_smax2
        ku4 = (zeta_val * th4 - sg4) * inv_visc
        dth4 = -big_gamma_q * th4 + LAM_ALIGN * ku4 - 0.35 * (th4 * th4 * th4)
        dsg4 = -(f4 * inv_taup) * sg4 + g_p * ku4

        th += (dt / 6.0) * (dth1 + 2.0 * dth2 + 2.0 * dth3 + dth4)
        sg += (dt / 6.0) * (dsg1 + 2.0 * dsg2 + 2.0 * dsg3 + dsg4)

    zc_times = [
        t_ss[j] for j in range(len(v_ss) - 1) if v_ss[j] < 0.0 and v_ss[j + 1] >= 0.0
    ]
    if len(zc_times) >= 3:
        diffs = [zc_times[j + 1] - zc_times[j] for j in range(len(zc_times) - 1)]
        period = sum(diffs) / len(diffs)
    else:
        period = float("inf")
    freq_mhz = (1000.0 / period) if math.isfinite(period) else 0.0
    v_peak_ums = max(abs(v) for v in v_ss) * 1000.0

    # Fundamental harmonic phase lag via single-bin complex Fourier projection at measured period
    n_ss = len(th_ss)
    mean_th = sum(th_ss) / n_ss
    mean_sg = sum(sg_ss) / n_ss
    omega_f = (2.0 * math.pi / period) if math.isfinite(period) else 0.1
    rot = cmath.exp(-1j * omega_f * dt)
    w = 1.0 + 0.0j
    acc_th = 0.0 + 0.0j
    acc_sg = 0.0 + 0.0j
    for j in range(n_ss):
        acc_th += (th_ss[j] - mean_th) * w
        acc_sg += (sg_ss[j] - mean_sg) * w
        w *= rot
    phase_deg = math.degrees(cmath.phase(acc_th / acc_sg)) % 360.0
    max_sg = max(abs(s) for s in sg_ss)
    shear_deg = math.degrees(math.atan(max_sg / max(g_p, 1e-9)))
    return period, freq_mhz, v_peak_ums, phase_deg, shear_deg


def pack_fowler_byte(trits5: list[int]) -> int:
    weights = [1, 3, 9, 27, 81]
    return sum((t + 1) * w for t, w in zip(trits5, weights))


def main() -> None:
    parser = argparse.ArgumentParser(description="Verify Mori-Zwanzig Active DNA Hopf Bifurcation (Pure Stdlib)")
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()
    _ = random.Random(args.seed)

    t0 = time.perf_counter()
    k_macro = math.pi / 2.0
    k_micro = math.pi / 0.15

    # Gate 1 & 2 assertions
    _, pinf_lam, gp_lam, taup_lam = dna_rheology(16.5)
    m_lam = mode_stability(k_macro, gp_lam, taup_lam)
    assert pinf_lam == 0.0 and m_lam["det"] < 0.0

    _, pinf_t4, gp_t4, taup_t4 = dna_rheology(56.4)
    m_t4_macro = mode_stability(k_macro, gp_t4, taup_t4)
    m_t4_micro = mode_stability(k_micro, gp_t4, taup_t4)
    assert m_t4_macro["de_nem"] > m_t4_macro["de_crit"]
    assert m_t4_macro["tr"] > 0.0 and m_t4_macro["det"] > 0.0 and m_t4_macro["disc"] < 0.0
    assert m_t4_micro["tr"] < 0.0

    # Gate 3 & 4 assertions
    per_20, f_20, v_20, ph_20, sh_20 = simulate_limit_cycle(56.4, 2.0, 0.340)
    per_26, f_26, v_26, _, _ = simulate_limit_cycle(56.4, 2.6, 0.340)
    assert 40.0 < per_20 < 60.0 and v_20 > 30.0
    assert f_20 > f_26 and v_26 > v_20
    assert sh_20 > 10.0

    # Gate 5 Fowler byte lock
    assert pack_fowler_byte([0, 0, 0, 0, 0]) == FOWLER_NEUTRAL_BYTE
    assert pack_fowler_byte([-1, -1, -1, -1, -1]) == 0
    assert pack_fowler_byte([-1, 0, 1, 1, 1]) == 237

    elapsed = time.perf_counter() - t0
    print(json.dumps({
        "artifact": "verify_active_dna_mori_zwanzig_v1_0_1.py",
        "version": "1.0.1",
        "seed": args.seed,
        "period_2p0mm_s": round(per_20, 2),
        "freq_2p0mm_mhz": round(f_20, 2),
        "v_peak_2p0mm_ums": round(v_20, 2),
        "phase_lag_deg": round(ph_20, 2),
        "shear_angle_deg": round(sh_20, 2),
        "kintsugi_byte": KINTSUGI_WELD_BYTE,
        "elapsed_s": round(elapsed, 4),
        "exit_state": "CODE_SHOKUNIN_UNLATCHED = 255",
    }, indent=2))


if __name__ == "__main__":
    main()
