Active-Flow DNA Remodeling & Mori-Zwanzig Viscoelastic Hopf Code Bundle (v1.0.0)
* Artifact ID: active_flow_dna_viscoelastic_code_bundle_v1.0.0.md
* Version: v1.0.0 (Strict 3-Digit SemVer — November 11, 2026 Thread-Closure Release)
* Author: Justin D. Grimm (ORCID: 0009-0000-5299-0243)
* Affiliation: Independent Dandelion Substrate Research Group, Gurley, Alabama 35748, USA
* Companion Pearl: [PEARL] Active_Flow_DNA_Remodeling_Mori_Zwanzig_Viscoelastic_Oscillations_v1.0.0
* Primary Empirical Reference: Levanon, Goldberg, Cohen, Bouchbinder, Adar & Tayar, Nature Physics (2026), doi:10.1038/s41567-026-03505-y · arXiv:2511.22589
* Verification Command: python3 verify_active_dna_mori_zwanzig_v1_0_0.py --seed 42 (< 2.0 s runtime, NumPy/stdlib only)
1. Python Verification & Nonlinear FENE Limit-Cycle Harness (verify_active_dna_mori_zwanzig_v1_0_0.py)
n
#!/usr/bin/env python3
"""
Dandelion Substrate & Mori-Zwanzig Active DNA Viscoelastic Verification (v1.0.0)
Author: Justin D. Grimm (ORCID: 0009-0000-5299-0243)
Affiliation: Independent Dandelion Substrate Research Group, Gurley, Alabama, USA
Reference: Levanon et al., Nature Physics (2026), doi:10.1038/s41567-026-03505-y
"""


import argparse
import time
import numpy as np


PC_FCC = 0.1201635          # 3D FCC bond percolation threshold
FOWLER_NEUTRAL_BYTE = 121   # {0, 0, 0, 0, 0} resting / Hopf keyhole byte
KINTSUGI_WELD_BYTE = 251    # >10.0 deg transient shear-knot weld byte


ETA_S = 1.0                 # Solvent + unentangled MT shear viscosity [normalized]
ELL_V = 0.12                # MT bundle vorticity screening length [mm]
ZETA_DEFAULT = 0.340        # Extensile MT-kinesin motor activity
LAM_ALIGN = 1.20            # Nematic flow-alignment coupling
GAMMA_Q = 0.065             # Homogeneous MT orientational relaxation [s^-1]
D_Q = 0.025                 # Frank / bundle orientational diffusivity [mm^2/s]


def eta_eff(k: float) -> float:
return ETA_S * (1.0 + (ELL_V * k) ** 2)


def dna_rheology(lc_um: float) -> tuple[float, float, float, float]:
lc_ref = 16.5  # Phage lambda DNA (48.5 kbp = 16.5 um)
p_bond = float(1.0 - np.exp(-0.048 * (lc_um / 10.0) ** 1.8))
p_inf = 0.0 if p_bond <= PC_FCC else float(((p_bond - PC_FCC) / (1.0 - PC_FCC)) ** 0.41)
tau_p = float(2.0 * ((lc_um / lc_ref) ** 2.15))
g_p = float(0.115 * p_inf * ((lc_um / lc_ref) ** 0.55))
return p_bond, p_inf, g_p, tau_p


def mode_stability(k: float, g_p: float, tau_p: float, zeta_val: float = ZETA_DEFAULT) -> dict[str, float]:
visc = eta_eff(k)
big_gamma_q = GAMMA_Q + D_Q * (k ** 2)
gamma_act = (LAM_ALIGN * zeta_val) / visc - big_gamma_q
gamma_ve = (1.0 / tau_p) + (g_p / visc)


tr = gamma_act - gamma_ve
det_matrix = gamma_act * (-gamma_ve) - (-LAM_ALIGN / visc) * (g_p * zeta_val / visc)
det_exact = big_gamma_q * (g_p / visc) - gamma_act / tau_p
assert abs(det_matrix - det_exact) < 1e-12, "Exact algebraic determinant cancellation failed"


disc = tr ** 2 - 4.0 * det_exact
omega_lin = 0.5 * float(np.sqrt(max(0.0, -disc))) if det_exact > 0.0 else 0.0
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

def simulate_limit_cycle(
lc_um: float = 56.4,
l_box_mm: float = 2.0,
zeta_val: float = ZETA_DEFAULT,
steps: int = 35000,
dt: float = 0.05,
) -> tuple[float, float, float, float, float]:
_, _, g_p, tau_p = dna_rheology(lc_um)
k = float(np.pi / l_box_mm)
visc = eta_eff(k)
big_gamma_q = GAMMA_Q + D_Q * (k ** 2)
sigma_max = 0.12


th, sg = 1e-3, 0.0
t_arr = np.zeros(steps)
v_arr = np.zeros(steps)
th_arr = np.zeros(steps)
sg_arr = np.zeros(steps)


def rhs(theta_v: float, sigma_v: float) -> tuple[float, float, float]:
    fene = 1.0 + (sigma_v / sigma_max) ** 2
    ku = (zeta_val / visc) * theta_v - (1.0 / visc) * sigma_v
    dth = -big_gamma_q * theta_v + LAM_ALIGN * ku - 0.35 * (theta_v ** 3)
    dsg = -(fene / tau_p) * sigma_v + g_p * ku
    return dth, dsg, (ku / k)


for i in range(steps):
    t_arr[i] = i * dt
    dth1, dsg1, v1 = rhs(th, sg)
    v_arr[i], th_arr[i], sg_arr[i] = v1, th, sg
    dth2, dsg2, _ = rhs(th + 0.5 * dt * dth1, sg + 0.5 * dt * dsg1)
    dth3, dsg3, _ = rhs(th + 0.5 * dt * dth2, sg + 0.5 * dt * dsg2)
    dth4, dsg4, _ = rhs(th + dt * dth3, sg + dt * dsg3)
    th += (dt / 6.0) * (dth1 + 2.0 * dth2 + 2.0 * dth3 + dth4)
    sg += (dt / 6.0) * (dsg1 + 2.0 * dsg2 + 2.0 * dsg3 + dsg4)


half = steps // 2
v_ss, t_ss = v_arr[half:], t_arr[half:]
th_ss, sg_ss = th_arr[half:], sg_arr[half:]


zc = np.where((v_ss[:-1] < 0.0) & (v_ss[1:] >= 0.0))[0]
period = float(np.mean(np.diff(t_ss[zc]))) if len(zc) >= 3 else float("inf")
freq_mhz = (1000.0 / period) if np.isfinite(period) else 0.0
v_peak_ums = float(np.max(np.abs(v_ss)) * 1000.0)


fft_th = np.fft.rfft(th_ss - np.mean(th_ss))
fft_sg = np.fft.rfft(sg_ss - np.mean(sg_ss))
pk = int(np.argmax(np.abs(fft_th[1:])) + 1)
phase_deg = float(np.degrees(np.angle(fft_th[pk] / fft_sg[pk]))) % 360.0
shear_deg = float(np.degrees(np.arctan(np.max(np.abs(sg_ss)) / max(g_p, 1e-9))))
return period, freq_mhz, v_peak_ums, phase_deg, shear_deg

def pack_fowler_byte(trits5: list[int]) -> int:
weights = [1, 3, 9, 27, 81]
return sum((t + 1) * w for t, w in zip(trits5, weights))


def main() -> None:
parser = argparse.ArgumentParser(description="Verify Mori-Zwanzig Active DNA Hopf Bifurcation")
parser.add_argument("--seed", type=int, default=42)
args = parser.parse_args()
_ = np.random.default_rng(args.seed)


t0 = time.per
f_counter()
k_macro = float(np.pi / 2.0)
k_micro = float(np.pi / 0.15)


# Gate 1 & 2 assertions
_, pinf_lam, gp_lam, taup_lam = dna_rheology(16.5)
m_lam = mode_stability(k_macro, gp_lam, taup_lam)
assert pinf_lam == 0.0 and m_lam["det"] < 0.0, "Lambda DNA must remain below percolation in chaotic regime"


_, pinf_t4, gp_t4, taup_t4 = dna_rheology(56.4)
m_t4_macro = mode_stability(k_macro, gp_t4, taup_t4)
m_t4_micro = mode_stability(k_micro, gp_t4, taup_t4)
assert m_t4_macro["de_nem"] > m_t4_macro["de_crit"], "T4 DNA must exceed Deborah-Hopf threshold"
assert m_t4_macro["tr"] > 0.0 and m_t4_macro["det"] > 0.0 and m_t4_macro["disc"] < 0.0
assert m_t4_micro["tr"] < 0.0, "Microscopic active turbulence must be damped by viscoelasticity"


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
print(f"[OK-121] T4-DNA (56.4 um) Period={per_20:.2f}s ({f_20:.2f} mHz), v_peak={v_20:.2f} um/s, "
      f"PhaseLag={ph_20:.2f} deg, Kintsugi={KINTSUGI_WELD_BYTE} | Verified in {elapsed:.4f}s")

if name == "main":
main()



## 2. Strict ISO C99 Zero-Allocation Fixed-Point Gate (active_dna_hopf_gate_v1_0_0.c)


```c
/*
 * active_dna_hopf_gate_v1_0_0.c
 * Strict ISO C99 Fixed-Point (Q16.16) Deborah-Hopf & Epigenetic Phase Classifier
 * Compile: gcc -std=c99 -O3 -Wall -Wextra -Wconversion -Werror active_dna_hopf_gate_v1_0_0.c -o active_dna_hopf_gate
 */
#include <stdint.h>
#include <stdio.h>


#define FOWLER_ARRESTED_BYTE  ((uint8_t)0U)
#define FOWLER_NEUTRAL_BYTE   ((uint8_t)121U)
#define FOWLER_TURBULENT_BYTE ((uint8_t)237U)
#define KINTSUGI_WELD_BYTE    ((uint8_t)251U)
#define Q16_ONE               ((int64_t)65536)
#define KRAMHAUS_10DEG_TAN_Q16 ((int64_t)11555) /* tan(10.0 deg) = 0.176327 * 65536 */


typedef struct {
    int64_t gamma_act_q16;  /* Bare active growth rate [s^-1 in Q16.16] */
    int64_t gamma_q_q16;    /* Passive Frank/nematic relaxation [s^-1 in Q16.16] */
    int64_t g_over_eta_q16; /* Viscoelastic relaxation rate G_p / eta_eff [Q16.16] */
    int64_t inv_taup_q16;   /* Polymer reptation decay 1 / tau_p [Q16.16] */


    int64_t peak_strain_q16;/* Dimensionless peak polymer shear |sigma_p| / G_p [Q16.16] */
|---|


static inline int64_t q16_mul(int64_t a, int64_t b) {
    return (a * b) >> 16;
}


uint8_t classify_viscoelastic_state(const ActiveViscoStateQ16 *s, int8_t *out_trit) {
    const int64_t gamma_ve = s->inv_taup_q16 + s->g_over_eta_q16;
    const int64_t tr_q16   = s->gamma_act_q16 - gamma_ve;
    const int64_t det_q16  = q16_mul(s->gamma_q_q16, s->g_over_eta_q16)
                           - q16_mul(s->gamma_act_q16, s->inv_taup_q16);


    if (det_q16 <= 0) {
        *out_trit = (int8_t)(+1);
        return FOWLER_TURBULENT_BYTE;
    }
    if (tr_q16 > 0) {
        const int64_t tr_sq = q16_mul(tr_q16, tr_q16);
        if (tr_sq < (4 * det_q16)) {
            *out_trit = (int8_t)0;
            if (s->peak_strain_q16 > KRAMHAUS_10DEG_TAN_Q16) {
                return KINTSUGI_WELD_BYTE;
            }
            return FOWLER_NEUTRAL_BYTE;
        }
    }
    *out_trit = (int8_t)(-1);
    return FOWLER_ARRESTED_BYTE;
}


int main(void) {
    /* T4-DNA (56.4 um) in L_sys = 2.0 mm chamber */
    const ActiveViscoStateQ16 t4_state = {
        .gamma_act_q16   = 17564, /* 0.2680 s^-1 */
        .gamma_q_q16     = 8303,  /* 0.1267 s^-1 */
        .g_over_eta_q16  = 11757, /* 0.1794 s^-1 */
        .inv_taup_q16    = 2332,  /* 0.0356 s^-1 */
        .peak_strain_q16 = 30146  /* 0.460 (> tan 10 deg) -> Kintsugi 251 */
    };
    int8_t trit = -1;
    uint8_t code = classify_viscoelastic_state(&t4_state, &trit);
    if (trit != 0 || code != KINTSUGI_WELD_BYTE) {
        return 1;
    }
    printf("[OK-C99] T4 Hopf Waveguide Locked: trit=%d byte=%u\n", (int)trit, (unsigned)code);
    return 0;
}
3. Julia 1.12 Relational Multiple-Dispatch Module (ActiveDNAWaveguide_v1_0_0.jl)
module ActiveDNAWaveguide


export ActiveMotorFrame, DNAViscoelasticMedium, RemodelingFlow, refract


struct ActiveMotorFrame
    l_sys_mm::Float64
    zeta::Float64
    lam_align::Float64
end


struct DNAViscoelasticMedium
    lc_um::Float64
    g_p::Float64
    tau_p::Float64
end


struct RemodelingFlow
    gamma_q::Float64
    d_q::Float64
    ell_v::Float64
end


"""
Simultaneous triadic relational dispatch: refract(Subject, Object, Action)
Maps (-1 Chamber/Motor Frame, 0 DNA Mori-Zwanzig Medium, +1 Flow Action) -> (trit, byte, de_nem, de_crit)
"""
function refract(sub::ActiveMotorFrame, obj::DNAViscoelasticMedium, act::RemodelingFlow)
    k = pi / sub.l_sys_mm
    visc = 1.0 + (act.ell_v * k)^2
    big_gamma_q = act.gamma_q + act.d_q * k^2
    gamma_act = (sub.lam_align * sub.zeta) / visc - big_gamma_q
    gamma_ve = (1.0 / obj.tau_p) + (obj.g_p / visc)


    tr = gamma_act - gamma_ve
    det_j = big_gamma_q * (obj.g_p / visc) - gamma_act / obj.tau_p
    de_nem = obj.tau_p * big_gamma_q
    de_crit = obj.g_p > 0.0 ? (1.0 + visc / (obj.g_p * obj.tau_p)) : Inf


    if det_j <= 0.0
        return (+1, UInt8(237), de_nem, de_crit)
    elseif tr > 0.0 && tr^2 < 4.0 * det_j
        return (0, UInt8(121), de_nem, de_crit)
    else
        return (-1, UInt8(0), de_nem, de_crit)
    end
end


end # module
4. Lazy Predator DataCite Metadata (.zenodo.json & CITATION.cff)
{
  "title": "Active-Flow DNA Remodeling as a Mori-Zwanzig Viscoelastic Hopf Waveguide: Exact Deborah Criterion and Balanced Ternary Epigenetic State Dispatch (v1.0.0)",
  "upload_type": "software",
  "version": "1.0.0",
  "publication_date": "2026-11-11",
  "creators": [
    {
      "name": "Grimm, Justin D.",
      "affiliation": "Independent Dandelion Substrate Research Group, Gurley, Alabama, USA",
      "orcid": "0009-0000-5299-0243"
    }
  ],
  "license": "MIT",
  "keywords": [
    "active matter",
    "DNA remodeling",
    "Mori-Zwanzig formalism",
    "viscoelastic Hopf bifurcation",
    "percolation",
    "balanced ternary",
    "epigenetic state dispatch"
  ],
  "related_identifiers": [
    {
      "identifier": "10.1038/s41567-026-03505-y",
      "relation": "cites",
      "resource_type": "publication-article"
    },
    {
      "identifier": "10.48550/arXiv.2511.22589",
      "relation": "cites",
      "resource_type": "publication-preprint"
    },
    {
      "identifier": "10.5281/zenodo.21255229",
      "relation": "continues",
      "resource_type": "publication-preprint"
    }
  ]
}

cff-version: 1.2.0
title: "Active-Flow DNA Remodeling & Mori-Zwanzig Viscoelastic Hopf Code Bundle"
version: "1.0.0"
date-released: "2026-11-11"
authors:
  - family-names: "Grimm"
    given-names: "Justin D."
    orcid: "https://orcid.org/0009-0000-5299-0243"
    affiliation: "Independent Dandelion Substrate Research Group, Gurley, Alabama, USA"
license: "MIT"

