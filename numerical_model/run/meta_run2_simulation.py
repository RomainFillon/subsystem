# -*- coding: utf-8 -*-

##############################################################################
#replication for "The need for regulation of climate subsystems" [simulation]
###############################################################################
###############################################################################
import sys
#sys.path.append("C:\\Users\\Fillon\\Desktop\\scientifique\\P2_Amazon\\replication_package\\numerical_model\\model\\")
#sys.path.append("C:\\Users\\Fillon\\Desktop\\scientifique\\P2_Amazon\\replication_package\\numerical_model\\run\\")
sys.path.append(r"/user/rbf2132/github/model")
sys.path.append(r"/user/rbf2132/github/run")

#Define runs to be simulated
#run_tobesimulated= 'final_amazon_tcre_run0014/'
#run_tobesimulated= 'final_amazon_tcre_run0005/'

comparative_run_tobesimulated = 'final_amazon_tcre_run0004/'
run_tobesimulated= 'final_amazon_tcre_run0006/'

#comparative_run_tobesimulated = 'final_amazon_tcre_run0014/'
#run_tobesimulated= 'final_amazon_tcre_run0016/'


import numbers
import csv
import numpy as np
import warnings
import pickle
import gzip
import random
import scipy.integrate as si
import scipy.special as sc
from scipy.stats import norm, truncnorm, beta 
from scipy.signal import savgol_filter
from scipy.optimize import minimize, fsolve, fmin, fmin_cg, fmin_bfgs
from operator import itemgetter
import math
import pyparsing as pypars
from functools import partial
import os
import functools
import copy
import sys
import pandas as pd
import multiprocess
from multiprocess import Process
from scipy.ndimage import uniform_filter1d

#defining the paths for folders
pathsep=os.sep
base_folder = r'/user/rbf2132/github'
preprod_folder   = os.path.join(base_folder, 'run') + pathsep
model_folder     = os.path.join(base_folder, 'model') + pathsep
param_folder     = os.path.join(base_folder, 'parameters') + pathsep
outputs_folder   = os.path.join(base_folder, 'outputs') + pathsep
run_folder = outputs_folder + run_tobesimulated
#comparative_folder = outputs_folder + comparative_run_tobesimulated

run_id      = os.path.basename(run_tobesimulated.rstrip('/\\'))[-4:]   # '0006'

if "comparative_run_tobesimulated" in globals():
    comparative_folder = outputs_folder + comparative_run_tobesimulated
    comparative_run_id = os.path.basename(comparative_run_tobesimulated.rstrip('/\\'))[-4:]

shared_folder = r"/shared/share_fillonwagner/"
temp_folder   = os.path.join(shared_folder, "temp_results" + run_id)   # .../temp_results0006


#reading parameters and calculating exogenous trends
exec(open(model_folder+'read_parameters.py').read())

print(params["nb_models"])
print(params["id_climate"])

#calling functions
#load functions for chebyshev interpolation
from all_functions import nodes_Chebyshev, Chebyshev
from all_functions import n_terms_cheb, ind_cheb, fill_cheb, fill_cheb2

#load functions for Bellman equation
from all_functions import Bell_max

#load functions for dynamics
from all_functions import utility, marginal_utility, production, abatement_cost, damage_factor, emissions, \
    growth_primary_tipping, maximand, partialc_maximand, law_motion_SEU

#load functions for stochastic risk
from all_functions import density_betalaw 

#load functions for simulation
from all_functions import EV_simulation_stochastic2, integrand1_simulation, integrand_terminal1, \
        quadrature_product_terminal_EV, terminal_value

#load functions for stochastic risk 
from all_functions import density_betalaw



import numpy as np

_V3_DIAG_DONE = False 

# =====================================================================
#  compute_scd_block_v3
#
#  FIX IN THIS VERSION � capacity normalisation.
#  The state impulse is normalised to d_{t0} = 1 in value(), i.e. every
#  channel is expressed PER UNIT OF COVER DESTROYED. The capacity source
#  was left in absolute units, so the permanence channel was understated
#  by a factor d_{t0} = p * A_{t0}, about 10 for p = 0.1. It is now
#  divided by d0. Reading: per unit of cover destroyed, the ground
#  sterilised is zeta / (p A_{t0}), which for phi = 1 equals 1 / A_{T_f}.
#
#  Everything else is unchanged: seven propagation conventions, both
#  carbon kernels (stock = correct under a reservoir law of motion, flow =
#  the alternative for forcing subsystems), and the structural freeze.
#  No parameter is introduced anywhere; each convention is a labelled
#  choice, each derivative comes from growth_primary_tipping.
# =====================================================================

def compute_scd_block_v3(u, x, params, monetary, marginal_utility_temp,
                         SCC_3, param_delta, climate_arg, shock1_vec,
                         FREEZE=38, FREEZE_STRUCT=38, verbose=True):
    """SCD under every convention and both carbon kernels.

    Returns absolute estimators in the SCCDS metric, the full Phi anatomy,
    and the diagnostics used by the sensitivity battery."""
    global _V3_DIAG_DONE
    do_print = bool(verbose) and (not _V3_DIAG_DONE)

    T   = params["T"]
    bet = params["bet"]
    LCS = params["carbonstock_amz"]
    mut = marginal_utility_temp
    R   = params["rTCRE"] / params["TCRE_mean"]
    nP  = T - 1

    # ---- conversion size and capacity withdrawal --------------------------
    A_ref = u[:, 2]
    A38   = A_ref[FREEZE] if FREEZE < T else A_ref[-1]
    EMP   = (1.0 - A38) / A38
    phi_z = params.get("scd_phi", 1.0)
    p_cv  = params.get("scd_p",   0.10)
    zeta  = p_cv * (1.0 + phi_z * EMP)          # >= p by the geometric assumption
    d0    = max(p_cv * float(u[0, 2]), 1e-12)   # d_{t0} = p * A_{t0}
    zeta_per_unit = zeta / d0                   # ground sterilised per unit cover lost

    # ---- per-period anatomy ------------------------------------------------
    g_net = np.zeros((nP, 1)); AnuA  = np.zeros((nP, 1))
    frag  = np.zeros((nP, 1)); space = np.zeros((nP, 1))
    Theta = np.zeros((nP, 1)); sce_v = np.zeros((nP, 1))
    src   = np.zeros((nP, 1))

    for t in range(1, nP):
        cov = float(u[t - 1, 2]); tmp = float(u[t - 1, 1])
        ca  = climate_arg[t] if hasattr(climate_arg, "__len__") else climate_arg

        g_here = float(growth_primary_tipping(t, cov, tmp, ca, params))
        gp = float(growth_primary_tipping(t, (1 + param_delta) * cov, tmp, ca, params))
        gm = float(growth_primary_tipping(t, (1 - param_delta) * cov, tmp, ca, params))
        dg_dA = (gp - gm) / (2 * param_delta * cov)

        # fragility / space split, exact from the regeneration term
        base = max((params["Upsilon"] * (1 - cov)) / params["beta0"], 1e-6)
        space_part = -params["growth0"] * (1.0 - base ** params["eta"])   # < 0, heals
        frag_part  = dg_dA - space_part                                   # > 0

        # carbon loop  Theta = psi * chi > 0
        xi_coef  = (params["climate_maxEU"] * ca) if params["stochastic"] == 2 \
                   else params["climate_meanEU"]
        xi_prime = xi_coef * R / ((tmp - params["T0"]) * R + 1.0)
        theta_t  = cov * xi_prime * (-float(SCC_3[t - 1]))

        try:
            sce_t = float(np.asarray(params["mean_scenarioEU"].iloc[t]).ravel()[0])
        except Exception:
            sce_t = 0.0

        g_net[t - 1] = g_here; AnuA[t - 1] = cov * dg_dA
        frag[t - 1]  = cov * frag_part; space[t - 1] = cov * space_part
        Theta[t - 1] = theta_t; sce_v[t - 1] = sce_t

        # CAPACITY SOURCE, per unit of cover destroyed (this is the fix)
        src[t - 1] = cov * params["growth0"] * (1 - base ** params["eta"]) * zeta_per_unit

    # ---- propagation conventions -------------------------------------------
    one = np.ones((nP, 1))
    PHI = {
        "M":            one + g_net + AnuA + Theta,
        "N":            one + AnuA + Theta,                              # = M - g
        "M_noTheta":    one + g_net + AnuA,
        "N_noTheta":    one + AnuA,
        "M_spaceSuppr": one + g_net + frag + Theta,
        "N_spaceSuppr": one + frag + Theta,
        "N_climOnly":   one + g_net + AnuA + Theta - (g_net + sce_v),     # removes (nu - xi)
    }
    # structural freeze: beyond FREEZE_STRUCT the law of motion is the identity
    for k_ in PHI:
        PHI[k_][FREEZE_STRUCT - 1:] = 1.0
    src[FREEZE_STRUCT - 1:] = 0.0

    # ---- valuation ----------------------------------------------------------
    def value(Phi, use_src, kernel, chi_vec=None):
        """Forward-recurse the deviation from a unit impulse, then value it.
           stock : one term per k (correct for a carbon reservoir)
           flow  : cumulative over j (correct for a forcing subsystem)"""
        chi = SCC_3 if chi_vec is None else chi_vec
        out = np.zeros((nP, 1))
        for t in range(nP):
            max_k = len(mut) - t - 1
            d = 1.0; run = 0.0; tot = 0.0
            for k in range(1, max_k):
                i = t + k - 1
                ph = float(Phi[i]) if i < nP else 1.0
                s_ = float(src[i]) if (use_src and i < nP) else 0.0
                d = ph * d + s_
                c_ = float(chi[t + k]) if (t + k) < len(chi) else 0.0
                if kernel == "stock":
                    tot += (bet ** k) * mut[t + k + 1] * c_ * d
                else:
                    run += c_ * d
                    tot += (bet ** k) * mut[t + k + 1] * run
            out[t] = tot
        return out

    # under the flow kernel the initial pulse is a separate channel; under the
    # stock kernel it is already inside value() (set d = 1, Phi = 1 to see it)
    standard_T = np.zeros((nP, 1))
    for t in range(nP):
        max_k = len(mut) - t - 1
        sk = 0.0
        c1 = float(SCC_3[t + 1]) if (t + 1) < len(SCC_3) else 0.0
        for k in range(0, max_k):
            sk += bet ** k * mut[t + k + 1] * c1
        standard_T[t] = sk

    sl = slice(0, T - 2)
    def scale(a):
        return (1.0 / LCS) * monetary[sl] * a[sl]
 
    # k = 0 : le carbone rel�ch� par la conversion (d_0 = 1) r�chauffe � t+1,
    # pond�r� par chi[t] = SCC_3[t]. value() d�marre � k=1 et l'omet ; on le
    # r�tablit ici pour le noyau STOCK uniquement (le FLOW le re�oit d�j� via
    # standard_T, donc l'ajouter l� double-compterait). Ind�pendant de la
    # convention et de la source : d_0 = 1 dans tous les cas.
    direct_T = np.zeros((nP, 1))
    for t in range(nP):
        if (t + 1) < len(mut) and t < len(SCC_3):
            direct_T[t] = float(mut[t + 1]) * float(SCC_3[t])      # bet**0 * mut[t+1] * chi[t] * 1
 
    E = {}
    for kern in ["stock", "flow"]:
        add_T = standard_T if kern == "flow" else direct_T          # STOCK : k=0 r�tabli
        for nm, Ph in PHI.items():
            E[f"rev_{nm}_{kern}"]  = scale(value(Ph, False, kern) + add_T)
            E[f"perm_{nm}_{kern}"] = scale(value(Ph, True,  kern) + add_T)
# -----------------------------------------------------------------------------


    # naive object: temperature weighted by dT/dS instead of dT/dA. Reported to
    # expose the sign trap; it is NOT a valid SCD.
    chi_wrong = np.asarray(shock1_vec[:len(SCC_3)]).reshape(-1, 1)
    E["naive_wrongT_stock"] = scale(value(PHI["M"], False, "stock", chi_vec=chi_wrong))

    # ---- diagnostics ---------------------------------------------------------
    if do_print:
        _V3_DIAG_DONE = True
        print("\n" + "=" * 78)
        print("  SCD v3 stock vs flow accounting, all conventions")
        print("=" * 78)
        print(f"  L = {LCS}   zeta = {zeta:.4f}   d0 = p*A_t0 = {d0:.4f}   "
              f"zeta per unit cover = {zeta_per_unit:.4f}")
        print(f"  structural freeze at t = {FREEZE_STRUCT}")
        print("\n  --- Phi anatomy (frozen to 1 after the freeze) ---")
        print("    t |     g     |  A*nuA (frag/space)   |  Theta |  src   | Phi_M   Phi_N")
        for tt in [0, 5, 10, 20, 30, 36, 37, 40, 60]:
            if tt < nP:
                print(f"  {tt:3d} | {float(g_net[tt]):+9.5f} | "
                      f"{float(AnuA[tt]):+7.4f}({float(frag[tt]):+.3f}/{float(space[tt]):+.3f}) | "
                      f"{float(Theta[tt]):+6.3f} | {float(src[tt]):+6.4f} | "
                      f"{float(PHI['M'][tt]):6.4f} {float(PHI['N'][tt]):6.4f}")
        print("\n  --- prod(Phi) : does the deviation heal? ---")
        for nm in ["M", "N", "M_noTheta", "M_spaceSuppr", "N_spaceSuppr", "N_climOnly"]:
            P = np.asarray(PHI[nm][:T - 2]).ravel()
            cp = np.cumprod(np.r_[1.0, P])
            pre = P[:max(FREEZE_STRUCT - 1, 1)]
            print(f"    Phi_{nm:13s} min={P.min():+.4f} max={P.max():+.4f} "
                  f"share>=1={np.mean(P >= 1.0):5.1%} (pre-freeze {np.mean(pre >= 1.0):5.1%})  "
                  f"prod@10/30/60={cp[min(10,len(cp)-1)]:.3f}/"
                  f"{cp[min(30,len(cp)-1)]:.3f}/{cp[min(60,len(cp)-1)]:.3f}")
        print("\n  --- ESTIMATORS at t0 : STOCK vs FLOW ---")
        print("    convention          |   reversible        |   permanent")
        print("                        |  stock      flow    |  stock      flow")
        for nm in PHI:
            print(f"    {nm:19s} | {float(E[f'rev_{nm}_stock'][0]):9.2f} "
                  f"{float(E[f'rev_{nm}_flow'][0]):9.2f} | "
                  f"{float(E[f'perm_{nm}_stock'][0]):9.2f} "
                  f"{float(E[f'perm_{nm}_flow'][0]):9.2f}")
        print(f"    naive_wrongT (stock) = {float(E['naive_wrongT_stock'][0]):.2f}  "
              "<- wrong temperature weight, exposes the sign trap")
        pm = float(E["perm_M_stock"][0]); rm = float(E["rev_M_stock"][0])
        print(f"\n  --- permanence gap (perm - rev)/perm at t0 = "
              f"{100*(pm-rm)/pm:+.3f}%   [was ~1% before the capacity fix]")
        print("\n  --- SIGNS over the horizon ---")
        for nm in ["perm_M_stock", "perm_N_stock", "perm_N_spaceSuppr_stock",
                   "perm_M_flow", "naive_wrongT_stock"]:
            v = np.asarray(E[nm]).ravel(); v = v[np.abs(v) > 1e-12]
            if v.size:
                print(f"    {nm:26s}: share<0 = {np.mean(v < 0):6.1%}  "
                      f"min={v.min():11.2f}  max={v.max():11.2f}")
        print("=" * 78 + "\n")

    return dict(
        SCD_M=E["perm_M_stock"], SCD_N=E["perm_N_stock"], dVdA=E["rev_M_stock"],
        estimators=E, PHI=PHI,
        Phi_M=PHI["M"], Phi_N=PHI["N"], Phi_bio=PHI["M_noTheta"],
        Phi_permM=PHI["M_spaceSuppr"], Phi_permN=PHI["N_spaceSuppr"],
        g_net=g_net, AnuA=AnuA, frag=frag, space=space, Theta=Theta,
        src_cap=src, zeta=zeta, d0=d0, zeta_per_unit=zeta_per_unit,
        SCD_T=standard_T[sl], SCD_A_M=value(PHI["M"], False, "stock")[sl],
        theta_ratio=float(np.max(np.abs(Theta)) / (np.max(np.abs(AnuA)) + 1e-12)))

# =====================================================================
#  scd_sensitivity  �  v3-compatible
#
#  FIX vs the previous version:
#    - v3 renamed the estimators: "perm_M" -> "perm_M_stock" / "perm_M_flow".
#      The headline key is now a parameter (HEAD), default "perm_M_stock".
#    - v3's signature takes shock1_vec, so the battery must carry it.
#    - guards added so a missing key degrades to a warning, not a KeyError
#      that kills the whole pool.
#
#  Replace the whole previous block (lines ~150-560) with this.
# =====================================================================

import numpy as np

HEAD = "perm_M_stock"          # headline estimator used for the elasticities


def _get(scd, key, fallback=None):
    E = scd.get("estimators", {})
    if key in E:
        return E[key]
    if fallback is not None and fallback in E:
        return E[fallback]
    return None


# ---------------------------------------------------------------------
def s1_tipping_active(u, params, scd, T):
    print("\n S1  IS THE TIPPING NON-LINEARITY ACTIVE ?")
    print(" " + "-" * 70)
    print("   b = Upsilon(1-A)/beta0 ; regeneration carries (1 - b^eta).")
    print("   b^eta ~ 0 => linear regime: no tipping, and fragility F ~ 0.")
    print("   t   |    A      |     b      |   b^eta    | 1-b^eta |   F_t    |   S_t")
    for tt in [0, 5, 10, 20, 30, 37, 40, 60, 80]:
        if tt < T - 2:
            A = float(u[tt, 2])
            b = (params["Upsilon"] * (1 - A)) / params["beta0"]
            be = b ** params["eta"]
            print(f"  {tt:4d} | {A:9.5f} | {b:10.4e} | {be:10.4e} | {1-be:7.4f} "
                  f"| {float(scd['frag'][tt]):+8.5f} | {float(scd['space'][tt]):+8.5f}")
    F = np.abs(np.asarray(scd["frag"][:T - 2])); S = np.abs(np.asarray(scd["space"][:T - 2]))
    print(f"   max|F| = {F.max():.3e}   max|S| = {S.max():.3e}   "
          f"F/S = {F.max()/max(S.max(), 1e-30):.4e}")
    print("   -> F/S < 1e-3 means the tipping mechanism is OFF along this path.")


# ---------------------------------------------------------------------
def s5_decay_budget(scd, T):
    print("\n S5  WHAT SETS THE HALF-LIFE OF THE DEVIATION ?")
    print(" " + "-" * 70)
    Phi = np.asarray(scd["Phi_M"][:T - 2]).ravel()
    g   = np.asarray(scd["g_net"][:T - 2]).ravel()
    F   = np.asarray(scd["frag"][:T - 2]).ravel()
    S   = np.asarray(scd["space"][:T - 2]).ravel()
    Th  = np.asarray(scd["Theta"][:T - 2]).ravel()
    mp = Phi.mean()
    hl = np.log(0.5) / np.log(mp) if 0 < mp < 1 else np.inf
    print(f"   mean Phi_M = {mp:.6f}   =>  half-life = {hl:.1f} periods")
    den = (1 - mp) if abs(1 - mp) > 1e-15 else np.nan
    print(f"   decomposition of (1 - Phi_M) = {1-mp:.6f} :")
    for nm, v in [("healing  -S", -S.mean()), ("fragility -F", -F.mean()),
                  ("growth   -g", -g.mean()), ("loop -Theta", -Th.mean())]:
        print(f"      {nm:14s} = {v:+.6f}   ({v/den:+7.1%})")
    print("   NOTE: Phi is forced to 1 after the structural freeze (t>=38),")
    print("         so the mean above mixes an active and a frozen regime.")
    pre = Phi[:min(37, len(Phi))]
    if pre.size:
        mpp = pre.mean()
        hlp = np.log(0.5) / np.log(mpp) if 0 < mpp < 1 else np.inf
        print(f"         pre-freeze only: mean Phi = {mpp:.6f}, half-life = {hlp:.1f}")


# ---------------------------------------------------------------------
def s2_param_elasticity(u, x, params, monetary, mut, SCC_3, param_delta,
                        climate_arg, shock1_vec, scd, SCC_fossil, T):
    print("\n S2  PROPAGATION-CHANNEL ELASTICITY  (baseline path held fixed)")
    print(" " + "-" * 70)
    base = _get(scd, HEAD)
    if base is None:
        print(f"   [skip] estimator '{HEAD}' absent. available: "
              f"{sorted(scd.get('estimators', {}).keys())[:6]} ...")
        return
    base_ratio = float(base[0]) / float(SCC_fossil[0])
    print(f"   reference {HEAD}: /SCC_foss = {base_ratio:.4f}")
    print("   param        factor |  mean Phi_M | prod@60 |  ratio  | elasticity")
    for key in ["growth0", "eta", "Upsilon", "beta0"]:
        if key not in params:
            continue
        for fac in [0.5, 1.5]:
            p2 = dict(params); p2[key] = params[key] * fac
            try:
                s2 = compute_scd_block_v3(u, x, p2, monetary, mut, SCC_3,
                                          param_delta, climate_arg, shock1_vec,
                                          verbose=False)
                v = _get(s2, HEAD)
                if v is None:
                    raise KeyError(HEAD)
            except Exception as e:
                print(f"   {key:10s}  x{fac:.1f}  -> failed: {type(e).__name__}: {e}")
                continue
            r  = float(v[0]) / float(SCC_fossil[0])
            P  = np.asarray(s2["Phi_M"][:T - 2]).ravel()
            pr = float(np.prod(P[:min(60, len(P))]))
            el = ((r - base_ratio) / base_ratio) / (fac - 1.0)
            print(f"   {key:10s}  x{fac:.1f} | {P.mean():11.6f} | {pr:7.4f} | "
                  f"{r:7.3f} | {el:+8.3f}")
    print("   -> elasticity = %change in ratio per %change in the parameter.")
    print("      This holds the baseline path fixed: it measures the PROPAGATION")
    print("      channel only. Re-solve the value function before quoting a range.")


# ---------------------------------------------------------------------
def s3_linearisation(u, params, scd, climate_arg, T, p_conv=None):
    print("\n S3  LINEARISATION CHECK  (exact non-linear gap vs prod(Phi))")
    print(" " + "-" * 70)
    p_conv = params.get("scd_p", 0.10) if p_conv is None else p_conv
    print(f"   perturbation p = {p_conv:.3f}; temperature held on the baseline path")
    print("   in BOTH worlds, so this isolates the non-linearity in A.")
    A = np.asarray(u[:, 2]).ravel()
    Ab = A[0]; Ac = (1.0 - p_conv) * A[0]
    d_nl = [Ab - Ac]
    for t in range(1, min(T - 2, len(A) - 1)):
        tmp = float(u[t - 1, 1])
        ca = climate_arg[t] if hasattr(climate_arg, "__len__") else climate_arg
        gb = float(growth_primary_tipping(t, Ab, tmp, ca, params))
        gc = float(growth_primary_tipping(t, max(Ac, 1e-9), tmp, ca, params))
        Ab = Ab * (1.0 + gb); Ac = Ac * (1.0 + gc)
        d_nl.append(Ab - Ac)
    d_nl = np.array(d_nl)
    Phi0 = np.asarray(scd["Phi_bio"]).ravel()
    d_lin = d_nl[0] * np.cumprod(np.r_[1.0, Phi0[:len(d_nl) - 1]])
    print("    t   |  d exact  |  d linear |  rel err")
    for tt in [0, 5, 10, 20, 30, 36, 40, 60]:
        if tt < len(d_nl):
            e = (d_lin[tt] - d_nl[tt]) / max(abs(d_nl[tt]), 1e-30)
            print(f"  {tt:5d} | {d_nl[tt]:9.6f} | {d_lin[tt]:9.6f} | {e:+8.3%}")
    n = min(len(d_nl), len(d_lin))
    rel = np.abs(d_lin[:n] - d_nl[:n]) / np.maximum(np.abs(d_nl[:n]), 1e-30)
    print(f"   max rel err = {np.nanmax(rel):.3%}")
    print("   NOTE: the non-linear sim above does NOT freeze A after t=38, while")
    print("         Phi does. Compare only over t < 38.")


# ---------------------------------------------------------------------
def s4_drought(u, x, params, monetary, mut, SCC_3, param_delta,
               climate_arg, shock1_vec, scd, SCC_fossil, T):
    print("\n S4  DROUGHT INTERMITTENCY  (Theta = psi*chi is 0 whenever psi = 0)")
    print(" " + "-" * 70)
    s2v = np.asarray(climate_arg).ravel()
    print(f"   shock2: n={s2v.size} min={s2v.min():.3e} max={s2v.max():.3e} "
          f"median={np.median(s2v):.3e}")
    print(f"   share<1e-20: {np.mean(s2v < 1e-20):.1%}   share>0.1: {np.mean(s2v > 0.1):.1%}")
    pos = s2v[s2v > 1e-20]
    if pos.size:
        q = np.percentile(np.log10(pos), [1, 25, 50, 75, 99])
        print(f"   log10 of non-negligible shocks q01/25/50/75/99: "
              f"{q[0]:.2f}/{q[1]:.2f}/{q[2]:.2f}/{q[3]:.2f}/{q[4]:.2f}")
    print("   -> a clean bimodal split = event-type shock (fine). A continuum from")
    print("      1e-300 to 1e0 = the shock law is underflowing.")
    print("\n   counterfactual intensities (LABELLED, not a recalibration):")
    print("   factor |  mean Phi_M | prod@60 |  ratio")
    for fac in [0.5, 1.0, 2.0]:
        try:
            s2b = compute_scd_block_v3(u, x, params, monetary, mut, SCC_3,
                                       param_delta,
                                       np.asarray(climate_arg, dtype=float) * fac,
                                       shock1_vec, verbose=False)
            v = _get(s2b, HEAD)
            if v is None:
                raise KeyError(HEAD)
        except Exception as e:
            print(f"   x{fac:.1f} -> failed: {type(e).__name__}: {e}")
            continue
        P = np.asarray(s2b["Phi_M"][:T - 2]).ravel()
        print(f"   x{fac:4.1f}  | {P.mean():11.6f} | "
              f"{float(np.prod(P[:min(60,len(P))])):7.4f} | "
              f"{float(v[0])/float(SCC_fossil[0]):7.3f}")


# ---------------------------------------------------------------------
def s6_kernel_gap(scd, SCC_fossil):
    """The single most important comparison: stock vs flow accounting."""
    print("\n S6  STOCK vs FLOW ACCOUNTING  (the kernel correction)")
    print(" " + "-" * 70)
    E = scd.get("estimators", {})
    f0 = float(SCC_fossil[0])
    print("   convention        | stock ratio | flow ratio | flow/stock")
    for nm in ["M", "N", "M_noTheta", "M_spaceSuppr", "N_spaceSuppr", "N_climOnly"]:
        a = E.get(f"perm_{nm}_stock"); b = E.get(f"perm_{nm}_flow")
        if a is None or b is None:
            continue
        ra, rb = float(a[0]) / f0, float(b[0]) / f0
        print(f"   perm_{nm:13s} | {ra:11.4f} | {rb:10.4f} | {rb/max(ra,1e-30):9.2f}")
    nv = E.get("naive_wrongT_stock")
    if nv is not None:
        print(f"   naive_wrongT (stock) ratio = {float(nv[0])/f0:+.4f}  "
              "<- sign trap, not a valid SCD")
    print("   -> the law of motion adds carbon via (TREE - TREEn), i.e. a STOCK.")
    print("      The stock column is the one consistent with that law.")


# ---------------------------------------------------------------------
def run_sensitivity_battery(u, x, params, monetary, marginal_utility_temp,
                            SCC_3, param_delta, climate_arg, scd, SCC_fossil,
                            shock1_vec=None):
    T = params["T"]
    if shock1_vec is None:
        shock1_vec = np.zeros(params["T"])
    print("\n" + "=" * 74)
    print("  SENSITIVITY BATTERY")
    print("=" * 74)
    for fn, args in [(s6_kernel_gap, (scd, SCC_fossil)),
                     (s1_tipping_active, (u, params, scd, T)),
                     (s5_decay_budget, (scd, T)),
                     (s2_param_elasticity, (u, x, params, monetary,
                                            marginal_utility_temp, SCC_3,
                                            param_delta, climate_arg,
                                            shock1_vec, scd, SCC_fossil, T)),
                     (s3_linearisation, (u, params, scd, climate_arg, T)),
                     (s4_drought, (u, x, params, monetary, marginal_utility_temp,
                                   SCC_3, param_delta, climate_arg, shock1_vec,
                                   scd, SCC_fossil, T))]:
        try:
            fn(*args)
        except Exception as e:
            print(f"\n [{fn.__name__} skipped: {type(e).__name__}: {e}]")
    print("\n NOT TESTABLE HERE (state explicitly in the paper):")
    print("   - spatial homogeneity: A is a scalar, competition is mean-field.")
    print("     Regrowth into a compact hole is slower than gamma0 implies,")
    print("     which biases the SCD DOWNWARD.")
    print("   - Pi^T = 1: released carbon never reabsorbs. Standard under TCRE,")
    print("     but it removes any sink by construction, biasing UP.")
    print("=" * 74 + "\n")
_SCDV2_DIAG_DONE = False   # print the heavy diagnostic once per process
_SCDV2_BR_DONE = False

def compute_scd_block_v2(u, x, params, monetary, marginal_utility_temp,
                         SCC_3, param_delta, climate_arg, FREEZE=38, verbose=True):
    """Return a dict of MANY SCD estimators (absolute, SCCDS metric) plus the
    full Phi anatomy, so we can see WHICH trap drives the numbers.
    SCC_3 carries -L*T_S (= -carbonstock_amz * shock1) exactly as in the branch.
    climate_arg is shock2 (drought) as passed to growth_primary_tipping."""
    global _SCDV2_DIAG_DONE
    do_print = bool(verbose) and (not _SCDV2_DIAG_DONE)

    T   = params["T"]
    bet = params["bet"]
    LCS = params["carbonstock_amz"]           # L : tC per unit cover (the "per tonne" base)
    mut = marginal_utility_temp
    R   = params["rTCRE"] / params["TCRE_mean"]

    # ---- zeta : capacity withdrawal, from YOUR existing scd_p / scd_phi -----
    A_ref = u[:, 2]
    A38   = A_ref[FREEZE] if FREEZE < T else A_ref[-1]
    EMP   = (1.0 - A38) / A38
    phi   = params.get("scd_phi", 1.0)
    p_cv  = params.get("scd_p",   0.10)
    zeta  = p_cv * (1.0 + phi * EMP)          # >= p by construction
    zeta_floor = p_cv                          # variant: zeta = p (no EMP amplification)

    # ---- per-period: growth, nu_A (frag/space split), Theta, all Phi --------
    g_net   = np.zeros((T - 1, 1))
    AnuA    = np.zeros((T - 1, 1))            # A * dg/dA  (Abar fixed: frag + space)
    frag    = np.zeros((T - 1, 1))            # A * fragility part  (positive)
    space   = np.zeros((T - 1, 1))            # A * space part      (negative, healing)
    Theta   = np.zeros((T - 1, 1))            # A * xi'(T) * T_S * L  (>=0)
    Phi_bio = np.zeros((T - 1, 1))            # 1 + g + A*nu_A                 (M, NO Theta)
    Phi_M   = np.zeros((T - 1, 1))            # 1 + g + A*nu_A + Theta         (M)
    Phi_N   = np.zeros((T - 1, 1))            # 1 +     A*nu_A + Theta         (N = M - g)
    Phi_Nb  = np.zeros((T - 1, 1))            # 1 +     A*nu_A                 (N, NO Theta)
    Phi_pM  = np.zeros((T - 1, 1))            # 1 + g + A*frag + Theta   (permanent: space suppressed)
    Phi_pN  = np.zeros((T - 1, 1))            # 1 +     A*frag + Theta   (permanent, N)
    src_cap = np.zeros((T - 1, 1))            # per-period cover lost from withdrawn ground
    src_cap0= np.zeros((T - 1, 1))            # same with zeta = p (floor variant)

    for t in range(1, T - 1):
        cov = float(u[t - 1, 2]); tmp = float(u[t - 1, 1])
        ca  = climate_arg[t] if hasattr(climate_arg, "__len__") else climate_arg

        g_here = float(growth_primary_tipping(t, cov, tmp, ca, params))
        gp = float(growth_primary_tipping(t, (1 + param_delta) * cov, tmp, ca, params))
        gm = float(growth_primary_tipping(t, (1 - param_delta) * cov, tmp, ca, params))
        dg_dA = (gp - gm) / (2 * param_delta * cov)            # Abar fixed: total nu_A

        # frag/space split, EXACT from your a-term (base stays on real cover):
        #   a = growth0 (1 - base^eta)(1 - cov),  base = Upsilon(1-cov)/beta0
        #   space part of da/dcov = -growth0 (1 - base^eta)   (= d a / d zeta)
        #   frag  part            = dg_dA - space part        (residual, positive)
        base = max((params["Upsilon"] * (1 - cov)) / params["beta0"], 1e-6)
        space_part = -params["growth0"] * (1.0 - base ** params["eta"])   # negative
        frag_part  = dg_dA - space_part                                    # positive

        # Theta : carbon -> temperature -> mortality feedback, folded into Phi.
        # xi = b = temp_regional(T) * xi_coef in growth_primary_tipping, so
        # xi'(T) = xi_coef * R / ((T - T0) R + 1). T_S * L = -SCC_3.
        xi_coef  = (params["climate_maxEU"] * ca) if params["stochastic"] == 2 \
                   else params["climate_meanEU"]
        xi_prime = xi_coef * R / ((tmp - params["T0"]) * R + 1.0)
        theta_t  = cov * xi_prime * (-float(SCC_3[t - 1]))

        g_net[t - 1]  = g_here
        AnuA[t - 1]   = cov * dg_dA
        frag[t - 1]   = cov * frag_part
        space[t - 1]  = cov * space_part
        Theta[t - 1]  = theta_t

        Phi_bio[t - 1] = 1.0 + g_here + cov * dg_dA
        Phi_M[t - 1]   = Phi_bio[t - 1] + theta_t
        Phi_N[t - 1]   = 1.0 + cov * dg_dA + theta_t
        Phi_Nb[t - 1]  = 1.0 + cov * dg_dA
        Phi_pM[t - 1]  = 1.0 + g_here + cov * frag_part + theta_t          # space healing removed
        Phi_pN[t - 1]  = 1.0 + cov * frag_part + theta_t

        src_cap[t - 1]  = cov * params["growth0"] * (1 - base ** params["eta"]) * zeta
        src_cap0[t - 1] = cov * params["growth0"] * (1 - base ** params["eta"]) * zeta_floor

    # ---- shared propagation kernel (state impulse or per-period capacity src)
    def propagate(Phi, src=None):
        out = np.zeros((T - 1, 1))
        for t in range(T - 1):
            max_k = len(mut) - t - 1
            sk = 0.0
            for k in range(1, max_k):
                sj = 0.0
                for j in range(1, k):
                    prod = np.prod(Phi[t + 1:t + j]) if j > 1 else 1.0
                    if src is None:
                        sj += SCC_3[t + j] * prod
                    elif (t + j - 1) < (T - 1):
                        sj += SCC_3[t + j] * float(src[t + j - 1]) * prod
                sk += (bet ** k) * sj * mut[t + k + 1]
            out[t] = sk
        return out

    # ---- temperature channel (direct pulse of the parcel on T) -------------
    standard_T = np.zeros((T - 1, 1))
    for t in range(T - 1):
        max_k = len(mut) - t - 1
        sk = 0.0
        for k in range(0, max_k):
            sk += bet ** k * mut[t + k + 1] * SCC_3[t + 1]
        standard_T[t] = sk

    # ---- state and capacity channels for every Phi variant -----------------
    st_M   = propagate(Phi_M);   cp_M   = propagate(Phi_M,  src_cap)
    st_N   = propagate(Phi_N);   cp_N   = propagate(Phi_N,  src_cap)
    st_bio = propagate(Phi_bio); cp_bio = propagate(Phi_bio, src_cap)
    st_Nb  = propagate(Phi_Nb);  cp_Nb  = propagate(Phi_Nb, src_cap)
    st_pM  = propagate(Phi_pM);  cp_pM  = propagate(Phi_pM, src_cap)
    st_pN  = propagate(Phi_pN);  cp_pN  = propagate(Phi_pN, src_cap)

    sl = slice(0, T - 2)
    def scale(x_):                                   # (1/L) * monetary * channel  -> $/tC metric
        return (1.0 / LCS) * monetary[sl] * x_[sl]

    # ---- THE ESTIMATOR MENU (each isolates a trap) -------------------------
    E = {}
    # T1/T12 : naive partial dU/dA = reversible on the multiplicative (M) law
    E["naive_dUdA"]      = scale(st_M + standard_T)
    # T2 : reversible (Abar fixed, space heals) � identical to naive by design
    E["reversible_M"]    = scale(st_M + standard_T)
    E["reversible_bio"]  = scale(st_bio + standard_T)          # T5: reversible WITHOUT Theta
    # T2/T3 : permanent (state + capacity), paper structure (shared Phi)
    E["perm_M"]          = scale(st_M  + cp_M  + standard_T)
    E["perm_N"]          = scale(st_N  + cp_N  + standard_T)
    # T5 : permanent WITHOUT Theta in the propagation
    E["perm_M_noTheta"]  = scale(st_bio + cp_bio + standard_T)
    E["perm_N_noTheta"]  = scale(st_Nb  + cp_Nb  + standard_T)
    # T4/T6 : permanent with space-suppressed Phi (capacity moves with state,
    #         so healing does NOT operate on the propagation) � the alternative
    #         that gives reversible and permanent genuinely DIFFERENT Phi.
    E["perm_M_spaceSuppr"] = scale(st_pM + cp_pM + standard_T)
    E["perm_N_spaceSuppr"] = scale(st_pN + cp_pN + standard_T)

    # channels kept separately for wiring / interpretation
    chan = dict(state_M=st_M, cap_M=cp_M, state_N=st_N, cap_N=cp_N, standard_T=standard_T)

    # ---- derived gaps (all as % of perm_M at t0) ---------------------------
    def pct(a, b):  # 100*(a-b)/b at t0, guarding zero
        b0 = float(b[0]); return np.nan if abs(b0) < 1e-12 else 100.0 * (float(a[0]) - b0) / b0
    gaps = {
        "capacity (perm_M - reversible)":       pct(E["perm_M"], E["reversible_M"]),
        "convention credit (perm_N - perm_M)":  pct(E["perm_N"], E["perm_M"]),
        "Theta effect (perm_M - perm_M_noTh)":  pct(E["perm_M"], E["perm_M_noTheta"]),
        "Phi structure (perm_M - spaceSuppr)":  pct(E["perm_M"], E["perm_M_spaceSuppr"]),
    }

    # ---------------------------- DIAGNOSTIC PRINTS -------------------------
    if do_print:
        _SCDV2_DIAG_DONE = True
        def prod_at(Phi, n):
            p = np.cumprod(np.r_[1.0, Phi.ravel()]); return float(p[min(n, len(p) - 1)])
        print("\n" + "=" * 78)
        print("  SCD FULL DIAGNOSTIC  (one draw; nothing ad hoc)")
        print("=" * 78)
        print(f"  L (carbonstock_amz) = {LCS:.4g}   [normaliser; note tC vs tCO2 = x44/12 if needed]")
        print(f"  zeta = {zeta:.4f}  (p={p_cv}, phi={phi}, EMP={EMP:.3f})   zeta_floor(=p) = {zeta_floor:.4f}")

        print("\n  --- Phi ANATOMY :  Phi_M = 1 + g + A*nu_A + Theta ;  Phi_N = Phi_M - g ---")
        print("    A*nu_A splits into fragility(+) and space(-, healing).")
        hdr = "    t |     g     | A*nuA (=frag+space)      | Theta  | Phi_bio Phi_M  Phi_N | Phi_permM"
        print(hdr)
        for tt in [0, 5, 10, 20, 30, 40, 50, 60]:
            if tt < T - 1:
                print(f"  {tt:3d} | {float(g_net[tt]):+8.4f} | "
                      f"{float(AnuA[tt]):+7.4f}({float(frag[tt]):+.3f}/{float(space[tt]):+.3f}) | "
                      f"{float(Theta[tt]):+6.3f} | {float(Phi_bio[tt]):6.3f} {float(Phi_M[tt]):6.3f} "
                      f"{float(Phi_N[tt]):6.3f} | {float(Phi_pM[tt]):6.3f}")

        print("\n  --- Phi SUMMARY (does the deviation heal <1 or accumulate >=1 ?) ---")
        for nm, P in [("Phi_bio", Phi_bio), ("Phi_M", Phi_M), ("Phi_N", Phi_N),
                      ("Phi_Nb", Phi_Nb), ("Phi_permM", Phi_pM), ("Phi_permN", Phi_pN)]:
            v = P[:T - 2]
            print(f"    {nm:9s}: min={float(v.min()):+.4f} max={float(v.max()):+.4f} "
                  f"mean={float(v.mean()):+.4f}  share>=1={float(np.mean(v >= 1.0)):5.1%}  "
                  f"prod@10/30/60 = {prod_at(v,10):.3g}/{prod_at(v,30):.3g}/{prod_at(v,60):.3g}")

        print("\n  --- theta_ratio (max|Theta| / max|A*nu_A|) = "
              f"{float(np.max(np.abs(Theta)) / (np.max(np.abs(AnuA)) + 1e-12)):.3f}")

        print("\n  --- ESTIMATORS at t0 (absolute, $/tC metric = (1/L)*monetary*channel) ---")
        for nm in ["naive_dUdA", "reversible_M", "reversible_bio",
                   "perm_M", "perm_N", "perm_M_noTheta", "perm_N_noTheta",
                   "perm_M_spaceSuppr", "perm_N_spaceSuppr"]:
            print(f"    {nm:20s} = {float(E[nm][0]):12.3f}")

        print("\n  --- DECOMPOSED GAPS (as % of perm_M at t0) ---")
        for k, v in gaps.items():
            print(f"    {k:38s}: {v:+8.3f}%")

        print("\n  --- ORDERING check (row-by-row, non-trunc rows) ---")
        m = (np.abs(E["perm_M"]) > 1e-9)
        o_rev = np.all(E["reversible_M"][m] <= E["perm_M"][m] + 1e-6)
        o_MN  = np.all(E["perm_M"][m]       <= E["perm_N"][m] + 1e-6)
        print(f"    reversible <= perm_M : {bool(o_rev)}    perm_M <= perm_N : {bool(o_MN)}")

        nz = np.sum(np.abs(E["perm_M"]) > 1e-9)
        print(f"\n  --- TRUNCATION : {int(T-2-nz)} tail rows are ~0 (horizon max_k->0), "
              f"ignore them in profiles. ---")
        print("=" * 78 + "\n")

    return dict(
        # headline objects used by the branch (paper convention M/N)
        dVdA=E["reversible_M"], SCD_M=E["perm_M"], SCD_N=E["perm_N"],
        SCD_T=standard_T[sl], SCD_A_M=chan["state_M"][sl], SCD_cap_M=chan["cap_M"][sl],
        # full estimator menu
        estimators=E, gaps=gaps, channels=chan,
        # Phi variants + anatomy (for further inspection / branch reuse)
        Phi_bio=Phi_bio, Phi_M=Phi_M, Phi_N=Phi_N, Phi_Nb=Phi_Nb,
        Phi_permM=Phi_pM, Phi_permN=Phi_pN,
        g_net=g_net, AnuA=AnuA, frag=frag, space=space, Theta=Theta,
        theta_ratio=float(np.max(np.abs(Theta)) / (np.max(np.abs(AnuA)) + 1e-12)),
        zeta=zeta, zeta_floor=zeta_floor, src_cap=src_cap)



_TESTS_DONE = False

def scd_sccds_tests(u, x, params, monetary, marginal_utility_temp,
                    SCC_2, SCC_3, SCC_4, standard_A, standard_T,
                    SCD_2, scd, shock_save1, shock_save2, param_delta,
                    SCCDS, SCC_fossil):
    """All per-draw robustness checks. Prints once per process."""
    global _TESTS_DONE
    if _TESTS_DONE:
        return
    _TESTS_DONE = True
    T = params["T"]; bet = params["bet"]; mut = marginal_utility_temp
    A = np.asarray(u[:, 2]).ravel()
    line = "-" * 74

    print("\n" + "=" * 74)
    print("  ROBUSTNESS SUITE  (one draw)")
    print("=" * 74)

    # ---- R1 : is the subsystem channel alive at all? -------------------
    print(f"\n R1  SUBSYSTEM CHANNEL ALIVE ?  (if SCC_4==0 then SCCDS==SCC exactly)\n{line}")
    s4 = np.asarray(SCC_4[:T - 2]).ravel()
    print(f"   SCC_4 (dg/dT): min={s4.min():+.4e} max={s4.max():+.4e} mean={s4.mean():+.4e}")
    print(f"   |SCC_4|<1e-12 in {int(np.sum(np.abs(s4) < 1e-12))}/{len(s4)} periods")
    _R = params["rTCRE"] / params["TCRE_mean"]
    print("   t  |   SCC_4 (fin.diff)   |   -xi'(T) (analytic)  |  ratio  | T_t-T0")
    for tt in [0, 1, 5, 10, 20, 30, 40, 60]:
        if tt < T - 2:
            tmp = float(u[tt, 1])
            ca = float(np.asarray(shock_save2).reshape(-1)[min(tt + 1, len(shock_save2) - 1)])
            xc = params["climate_maxEU"] * ca if params["stochastic"] == 2 else params["climate_meanEU"]
            xp = xc * _R / ((tmp - params["T0"]) * _R + 1.0)
            rat = float(s4[tt]) / (-xp) if abs(xp) > 1e-30 else np.nan
            print(f"  {tt:3d} | {float(s4[tt]):+18.6e}  | {-xp:+18.6e}  | {rat:+7.3f} | {tmp-params['T0']:+.5f}")
    print("   -> ratio should be ~1.0. If SCC_4 is 0 while -xi' is not, the finite")
    print("      difference is the problem (multiplicative (1+/-d)*T near T=T0).")

    # ---- R2 : does the path obey the law of motion? --------------------
    print(f"\n R2  LAW OF MOTION  A_(t+1) = A_t (1+g)\n{line}")
    g = np.asarray(scd["g_net"]).ravel()
    n = min(len(g), len(A) - 1)
    impl = A[1:n + 1] / np.maximum(A[0:n], 1e-12) - 1.0
    err = impl - g[:n]
    print(f"   max|implied g - model g| = {np.nanmax(np.abs(err)):.4e}   mean = {np.nanmean(np.abs(err)):.4e}")
    for tt in [0, 5, 10, 20, 30, 36, 37, 38, 45, 60]:
        if tt < n:
            print(f"   t={tt:3d}: implied={impl[tt]:+.6f}  model={g[tt]:+.6f}  err={err[tt]:+.3e}  A={A[tt]:.5f}")
    fz = np.where(np.abs(np.diff(A)) < 1e-14)[0]
    if len(fz):
        t_fz = int(fz[0]) + 1
        print(f"   STATE FREEZE at t={t_fz} (A={A[t_fz]:.5f}); {len(fz)}/{len(A)-1} periods frozen.")
        print("   -> Phi is built from the law, which the path does NOT follow there.")
        print("      Products of Phi crossing the freeze are not valid linearisations.")

    # ---- R3 : Phi behaviour and horizon mass ---------------------------
    print(f"\n R3  PROPAGATION Phi : heals (<1) or accumulates (>=1) ?\n{line}")
    for nm in ["Phi_bio", "Phi_M", "Phi_N", "Phi_permM", "Phi_permN"]:
        P = np.asarray(scd[nm][:T - 2]).ravel()
        cp = np.cumprod(np.r_[1.0, P])
        print(f"   {nm:10s} min={P.min():+.4f} max={P.max():+.4f} share>=1={np.mean(P>=1):5.1%} "
              f"prod@10/30/60={cp[min(10,len(cp)-1)]:.3f}/{cp[min(30,len(cp)-1)]:.3f}/{cp[min(60,len(cp)-1)]:.3f}")
    # where does the SCD mass come from? (fraction of the k-sum before the freeze)
    if len(fz):
        t_fz = int(fz[0]) + 1
        P = np.asarray(scd["Phi_M"]).ravel()
        cp = np.cumprod(np.r_[1.0, P])
        w = np.array([bet ** k * abs(float(mut[min(k + 1, len(mut) - 1)])) * cp[min(k, len(cp) - 1)]
                      for k in range(T - 2)])
        share_pre = w[:t_fz].sum() / max(w.sum(), 1e-30)
        print(f"   MASS: {share_pre:.1%} of the discounted propagation weight accrues BEFORE the")
        print(f"         freeze (t<{t_fz}); {1-share_pre:.1%} accrues in the frozen region.")

    # ---- R4 : sign of g flips -> ordering M<=N is conditional ----------
    print(f"\n R4  ORDERING  perm_M <= perm_N  requires g <= 0  (Phi_N = Phi_M - g)\n{line}")
    n_neg = int(np.sum(g[:T - 2] < 0)); n_pos = int(np.sum(g[:T - 2] > 0))
    print(f"   g<0 in {n_neg} periods, g>0 in {n_pos} periods "
          f"(first sign flip at t={int(np.argmax(g[:T-2] > 0)) if n_pos else -1})")
    E = scd["estimators"]
    m = np.abs(E["perm_M"]) > 1e-9
    viol = int(np.sum(E["perm_M"][m] > E["perm_N"][m] + 1e-6))
    print(f"   rows with perm_M > perm_N : {viol}/{int(m.sum())}")
    print("   -> the 'exoneration credit' becomes a PENALTY wherever g>0. Report the")
    print("      ordering proposition as conditional on the sign of g.")

    # ---- R5 : finite-difference step sensitivity -----------------------
    print(f"\n R5  FINITE-DIFFERENCE STEP (param_delta = {param_delta}) \n{line}")
    tt = 5
    cov = float(u[tt, 2]); tmp = float(u[tt, 1])
    ca = np.asarray(shock_save2).reshape(-1)[min(tt + 1, len(shock_save2) - 1)]
    for d in [param_delta * 0.1, param_delta, param_delta * 10]:
        gp = float(growth_primary_tipping(tt, (1 + d) * cov, tmp, ca, params))
        gm = float(growth_primary_tipping(tt, (1 - d) * cov, tmp, ca, params))
        dA = (gp - gm) / (2 * d * cov)
        gp2 = float(growth_primary_tipping(tt, cov, (1 + d) * tmp, ca, params))
        gm2 = float(growth_primary_tipping(tt, cov, (1 - d) * tmp, ca, params))
        dT_ = (gp2 - gm2) / (2 * d * tmp)
        print(f"   delta={d:.1e}:  dg/dA={dA:+.6e}   dg/dT={dT_:+.6e}")
    print("   -> both should be stable across a decade of delta. If dg/dT moves a lot,")
    print("      SCC_4 is a numerical artefact, not an economic zero.")

    # ---- R6 : L-invariance (units cannot flip the headline) ------------
    print(f"\n R6  L-INVARIANCE  (SCC_3 ~ L cancels the 1/L in the SCD scaling)\n{line}")
    print(f"   L = carbonstock_amz = {params['carbonstock_amz']} (TtC; = "
          f"{params['carbonstock_amz']*1000:.0f} GtC)")
    print("   -> ratio SCD/SCC is invariant to L up to the Theta channel (~1%).")
    print("      Rerun with L doubled to confirm the ratio moves by ~1%, not ~2x.")

    # ---- R7 : horizon truncation ---------------------------------------
    print(f"\n R7  HORIZON TRUNCATION\n{line}")
    sm = np.asarray(E["perm_M"]).ravel()
    nzc = int(np.sum(np.abs(sm) > 1e-9))
    print(f"   non-zero SCD rows: {nzc}/{len(sm)}   last non-zero at t={int(np.max(np.where(np.abs(sm)>1e-9))) if nzc else -1}")
    print(f"   SCD at t=T-10 relative to peak: {sm[max(0,len(sm)-10)]/max(np.abs(sm).max(),1e-30):.3f}")
    print("   -> the tail collapse is truncation, not economics. Quote SCD only over")
    print("      the interior of the horizon, or report the truncation-corrected profile.")

    # ---- R8 : estimator spread ------------------------------------------
    print(f"\n R8  ESTIMATOR SPREAD at t0 (which trap actually matters?)\n{line}")
    base = float(E["perm_M"][0])
    for nm in ["naive_dUdA", "reversible_M", "perm_M", "perm_N",
               "perm_M_noTheta", "perm_N_noTheta", "perm_M_spaceSuppr", "perm_N_spaceSuppr"]:
        v = float(E[nm][0])
        print(f"   {nm:20s} = {v:10.2f}   ({100*(v-base)/base:+7.3f}% vs perm_M)   "
              f"/SCC_foss = {v/float(SCC_fossil[0]):6.3f}")

    # ---- R9 : SCCDS reconstruction identity -----------------------------
    print(f"\n R9  SCCDS INTERNAL CONSISTENCY\n{line}")
    up = float(SCCDS[0] - SCC_fossil[0])
    print(f"   SCCDS t0 = {float(SCCDS[0]):.6f}   SCC_fossil t0 = {float(SCC_fossil[0]):.6f}")
    print(f"   subsystem uplift = {up:+.6e}  ({100*up/float(SCC_fossil[0]):+.5f}%)")
    print(f"   (II) standard_A t0 = {float(standard_A[0]):+.6e}")
    print(f"   loop  SCC_2*SCC_3*SCC_4 t0 = {float((SCC_2*SCC_3*SCC_4)[0]):+.6e}")
    print("   -> if uplift ~ 0 the paper's central claim is not reproduced by the run.")
    print("      Check R1 first: both terms carry SCC_4 as a common factor.")
    print("=" * 74 + "\n")


#load state space and chebyshev coefficients
u_max = np.genfromtxt(run_folder + '/u_max.csv', dtype=float, delimiter=';')
u_min = np.genfromtxt(run_folder + '/u_min.csv', dtype=float, delimiter=';')
params["u_max"]=u_max
params["u_min"]=u_min
coef_V = np.load(run_folder + 'coef_V.npy')

#initializing some value
#guess for minimizing Bell_max
#initialize (for  simulate.sce)
if params["dim"]==1:
    guess0=(1 - (params["bet"] * params["alpha"]))*np.ones((1,1))
    u0trend=params["K0"]/(((params["A0"])**(1/(1-params["alpha"])))*params["POP0"])*np.ones((T+2,1))
    u0=u0trend[[0],:]
if params["dim"]==3:
    guess0=np.array((0.1)).reshape(1,1)
    u0=np.concatenate((params["K0"]/((params["A0"]**(1/(1-params["alpha"])))*params["POP0"])*np.ones((1,1)),params["T0"]*np.ones((1,1)),params["TREE0"]*np.ones((1,1))),axis=1)
params["u0"]=u0
params["guess0"]=guess0

#degree interpolation depends on degree in each dimension
if params["dim"] == 1:
    degb = int(params["deg1"])
if params["dim"]==3:
    degb = int(np.max([params["deg1"], params["deg2"], params["deg3"]]))
degb = int(np.max([params["deg1"], params["deg2"], params["deg3"]]))
params['degb'] = degb

# matrix of weight (size: nb_terms_cheb *d)
# each line correspond to a multi-indice
weight_cheb = fill_cheb2(fill_cheb(params['degb'], params['dim']), params['dim'], params['deg1'], params['deg2'],params['deg3'])
params["weight_cheb"] = weight_cheb



#Define function for parallel simulation
def task(args):
    draw, params, coef_V = args
    final_doc= np.zeros((params["T"]-2, 37))
    x=np.zeros((params["T"],np.shape(params["guess0"])[1]))
    u=np.zeros((params["T"],params["dim"]))
    u[[0],:]=u0
    shock_save1=[0]*params["T"]
    shock_save2=[0]*params["T"]
    shock_save2_model=[0]*params["T"]

    #stochastic draws 	
    ss = seeds[draw]
    rng = np.random.default_rng(ss)
    seed_value = int(ss.generate_state(1)[0])
    random.seed(seed_value)

    # Seed Python's random module
    random.seed(seed_value)

    #from 1 to T
    for t_step in list(range(1,params["T"])):
        print(t_step)
        coefn=coef_V[[t_step],:]
        if t_step==1:
            guess=guess0
        else:
            guess=x[[t_step-1],:]
        if t_step <=20:#2100       
            (xopt,fopt, nb_iter, nb_funcalls, exitflag)=fmin(Bell_max,guess, args=(t_step,u, params, coefn),xtol=params["tol_simulation"], ftol= params["tol_simulation"],full_output=1,disp=0)#xopt is the minimum and fopt the value of the function at its minimum
            xopt=np.maximum(np.zeros(np.shape(xopt)),xopt)#transpose because x is vector
            xopt=np.minimum(np.ones(np.shape(xopt)),xopt)
            x[[t_step-1],:]=xopt
        else: #if >2100, check if we can just avoid minimization bellman.
            tolerance = 1e-2  # Tolerance de 2 ordres de grandeur
            if (all(np.isclose(x[t_step - 2, :], x[t_step - 3, :], rtol=tolerance)) and all(np.isclose(x[t_step - 4, :], x[t_step - 3, :], rtol=tolerance)) and all(np.isclose(x[t_step - 5, :], x[t_step - 4, :], rtol=tolerance))):
                x[[t_step-1],:]=x[[t_step-2],:]
            else:
                (xopt,fopt, nb_iter, nb_funcalls, exitflag)=fmin(Bell_max,guess, args=(t_step,u, params, coefn),xtol=params["tol_simulation"], ftol=params["tol_simulation"],full_output=1,disp=0)#xopt is the minimum and fopt the value of the function at its minimum
                xopt=np.maximum(np.zeros(np.shape(xopt)),xopt)#transpose because x is vector
                xopt=np.minimum(np.ones(np.shape(xopt)),xopt)
                x[[t_step-1],:]=xopt

    #random stochastic draw
    #law of motion state space
        if params["stochastic"]==0:
            u[[t_step], :] = law_motion_SEU(t_step, x[[t_step - 1]], u[[t_step - 1], :], params["TCRE_mean"],0, 0, params)
        if params["stochastic"]==1:
            shock1 = truncnorm.rvs(a=params["a_tcre"], b=params["b_tcre"],
                            loc=params["TCRE_mean"], scale=params["sdeviation"],
                            random_state=rng, size=1)
            shock_save1[t_step-1] = shock1[0]

            u[[t_step], :] = law_motion_SEU(t_step, x[[t_step-1]], u[[t_step-1], :], shock1, 0, 0, params)
        if params["stochastic"]==2:
            shock1 = truncnorm.rvs(params["a_tcre"], params["b_tcre"], loc=params["TCRE_mean"], scale=params["sdeviation"], random_state=rng, size=1)
            shock_save1[t_step-1] = shock1[0]
            shock2 = beta(params["beta_alpha"], params["beta_beta"]).rvs(random_state=rng)
            shock_save2[t_step-1] = shock2
            u[[t_step], :] = law_motion_SEU(t_step, x[[t_step - 1]], u[[t_step - 1], :], shock1,shock2, 0, params)
        if params["dim"]==1:
            if (u[[t_step], :] > params["u_max"][t_step]).any() or (u[[t_step], :] < params["u_min"][t_step]).any():
                break_file = os.path.join(run_folder, 'simulate_break.csv')
                if os.path.exists(break_file):
                    breaks = np.loadtxt(break_file, delimiter=";", ndmin=2)  # force 2D si 1 ligne
                    draws_to_skip = breaks[:, 0].astype(int)
                    t_steps = breaks[:, 1].astype(int)
                    print(f"Draws to skip (simulate_break.csv): {draws_to_skip}")
                else:
                    draws_to_skip = []
                return None
        else:
            if (u[[t_step],:]>params["u_max"][[t_step],:]).any() or (u[[t_step],:]<params["u_min"][[t_step],:]).any():
                print("state variable out of bounds in simulate, draw "+str(draw)+", time "+str(t_step))
                np.savetxt(run_folder+'simulate_break.csv', [draw,t_step], delimiter=";")

                return None

    final_doc[:, [0]] = x[0:params["T"] - 2]
    final_doc[:, [1]] = u[0:params["T"] - 2, [0]]
    final_doc[:, [2]] = u[0:params["T"] - 2, [1]]
    final_doc[:, [3]] = u[0:params["T"] - 2, [2]]
    
    #Now that we have control and states (stochastic)
    #We need SCC, SCCDS and SCD
    #We also need channels

    #marginal change for derivatives
    param_delta=0.001
    param_delta_A = 0.1
    dT = param_delta * u[:, [1]]
    dA = np.ones((params["T"], 1)) * param_delta_A
    Delta_u = np.hstack((np.zeros((params["T"], 1)), dT, np.zeros((params["T"], 1))))
    uplus = u + Delta_u
    uminus = u - Delta_u
    Delta_A = np.hstack((np.zeros((params["T"], 1)), np.zeros((params["T"], 1)), dA))
    vplus = u + Delta_A
    vminus = u - Delta_A

    #right dimensions
    sw_x = np.vstack((x[0:params["T"] - 1, :], x[[params["T"] - 2], :]))
    sw_time_horizon = params["time_horizon"][0:params["T"], :].reshape((params["T"]))

    #marginal derivative of utility wrt consumption, brought to present dollar terms
    monetary= params["bet"] *(1/partialc_maximand(sw_time_horizon, sw_x, u, params)[:params["T"] - 1, :])

    # marginal utility of temperature
    marginal_utility_temp = np.zeros((params["T"] - 1, 1))
    for t in range(0, params["T"] - 2):
        # du/dT
        marginal_utility_temp[t] = (maximand(t, x[[t], :], uplus[[t], :], params) - maximand(t, x[[t], :], uminus[[t], :], params)) / (2 * dT[[t], :])
        # dT/dS*du/dT
        if params["stochastic"] == 0:
            marginal_utility_temp[t] = marginal_utility_temp[t] * params["TCRE_mean"]
        else:
            marginal_utility_temp[t] = marginal_utility_temp[t] * shock_save1[t]

    if params["stochastic"] == 0:
        # ---- ENVELOPE PRICE, and nothing else --------------------------------
        # dU/dT from the value function. Note u_min is indexed with params["T"]
        # throughout; the original line used a bare `T` on one slice.
        SCC_1 = (approx_Chebyshev(coef_V[1:params["T"], :], params["weight_cheb"],
                                  uplus[1:params["T"], :],
                                  u_min[1:params["T"], :], u_max[1:params["T"], :])
                 - approx_Chebyshev(coef_V[1:params["T"], :], params["weight_cheb"],
                                    uminus[1:params["T"], :],
                                    u_min[1:params["T"], :], u_max[1:params["T"], :])
                 ) / (2 * dT[1:params["T"], :])
 
        # dU/dT . dT/dS, then to money. This is dV/dS: when the value function is
        # 3D the Amazon feedback is already inside it, so this IS the SCCDS; when
        # it is 2D it is the plain SCC. Either way it is the level to report and
        # nothing may be added on top without double counting.
        SCC_2      = SCC_1 * params["TCRE_mean"]
        SCC_fossil = - monetary * SCC_2
 
        _T2 = params["T"] - 2
        final_doc[:, [7]]  = SCC_fossil[0:_T2]     # SCCDS level
        final_doc[:, [18]] = SCC_fossil[0:_T2]     # same object, stochastic==2 convention
 
        # cross-draw arrays unused in this specification
        _sh3 = (params["T"] - 1, params["T"] - 1, params["T"] - 1)
        L1T_3D = np.zeros(_sh3); L2T_3D = np.zeros(_sh3)
        L1A_3D = np.zeros(_sh3); L2A_3D = np.zeros(_sh3)
 
 
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
    if params["stochastic"] == 1:
        # ---- ENVELOPE PRICE, and nothing else --------------------------------
        # dU/dT from the value function, integrated over the TCRE distribution.
        SCC_1 = np.zeros((params["T"] - 1, 1))
        for t in range(0, params["T"] - 1):
            cint, _ = si.quad(
                functools.partial(
                    integrand1_simulation,
                    t_step=t + 1, x_step=x[[t + 1], :], u_step=uplus[[t + 1], :],
                    coefn=coef_V[[t + 1], :], weight_cheb=params["weight_cheb"],
                    u_min=params["u_min"][[t + 1], :], u_max=params["u_max"][[t + 1], :],
                    params=params),
                params["TCRE_min"], params["TCRE_max"],
                epsabs=params["tolintabs"], epsrel=params["tolintrel"])
 
            cint2, _ = si.quad(
                functools.partial(
                    integrand1_simulation,
                    t_step=t + 1, x_step=x[[t + 1], :], u_step=uminus[[t + 1], :],
                    coefn=coef_V[[t + 1], :], weight_cheb=params["weight_cheb"],
                    u_min=params["u_min"][[t + 1], :], u_max=params["u_max"][[t + 1], :],
                    params=params),
                params["TCRE_min"], params["TCRE_max"],
                epsabs=params["tolintabs"], epsrel=params["tolintrel"])
 
            SCC_1[t] = (cint - cint2) / (2 * dT[[t + 1], :])
 
        # dU/dT . dT/dS with the drawn transient response, then to money.
        SCC_2 = np.zeros((params["T"] - 1, 1))
        for t in range(0, params["T"] - 2):
            SCC_2[t] = SCC_1[t] * shock_save1[t]
        SCC_fossil = - monetary * SCC_2
 
        _T2 = params["T"] - 2
        final_doc[:, [7]]  = SCC_fossil[0:_T2]     # SCCDS level, envelope, no double count
        final_doc[:, [18]] = SCC_fossil[0:_T2]     # same object, stochastic==2 convention
 
        # cross-draw arrays unused in this specification
        _sh3 = (params["T"] - 1, params["T"] - 1, params["T"] - 1)
        L1T_3D = np.zeros(_sh3); L2T_3D = np.zeros(_sh3)
        L1A_3D = np.zeros(_sh3); L2A_3D = np.zeros(_sh3)

# ----------------------------------------------------------------------------- 

    if params["stochastic"] == 2:

        # ---- (I) : dU/dT from the value function ----------------------------
        SCC_1 = np.zeros((params["T"] - 1, 1))
        for t in range(0, params["T"] - 1):
            cint = EV_simulation_stochastic2(
                t_step=t + 1, x_step=x[[t + 1], :], u_step=uplus[[t + 1], :],
                coefn=coef_V[[t + 1], :], weight_cheb=params["weight_cheb"],
                u_min=params["u_min"][[t + 1], :], u_max=params["u_max"][[t + 1], :],
                params=params, shock1=shock_save1[t + 1], shock2=shock_save2[t + 1])
            cint2 = EV_simulation_stochastic2(
                t_step=t + 1, x_step=x[[t + 1], :], u_step=uminus[[t + 1], :],
                coefn=coef_V[[t + 1], :], weight_cheb=params["weight_cheb"],
                u_min=params["u_min"][[t + 1], :], u_max=params["u_max"][[t + 1], :],
                params=params, shock1=shock_save1[t + 1], shock2=shock_save2[t + 1])
            SCC_1[t] = (cint - cint2) / (2 * dT[[t + 1], :])

        SCC_2 = np.zeros((params["T"] - 1, 1))                                # (I)
        for t in range(0, params["T"] - 2):
            SCC_2[t] = SCC_1[t] * shock_save1[t]
        SCC_3 = - np.ones((params["T"] - 1, 1)) * params["carbonstock_amz"]   # chi
        for t in range(0, params["T"] - 2):
            SCC_3[t] = SCC_3[t] * shock_save1[t]

        # ---- dA/dT : code stamp kept for comparison, PSI used everywhere -----
        SCC_4 = np.zeros((params["T"] - 1, 1))
        for t in range(1, params["T"] - 1):
            SCC_4[t - 1] = (growth_primary_tipping(t, u[[t], [2]], (1 + param_delta) * u[[t], [1]], shock_save2[t], params)
                            - growth_primary_tipping(t, u[[t], [2]], (1 - param_delta) * u[[t], [1]], shock_save2[t], params)) / (2 * param_delta * u[[t], [1]])
        SCC_4p = np.zeros((params["T"] - 1, 1))         # psi_t = A_t * dg/dT at state t
        for t in range(0, params["T"] - 1):
            cov_ = float(u[t, 2]); tmp_ = float(u[t, 1]); ca_ = shock_save2[t]
            gp = float(growth_primary_tipping(max(t, 1), cov_, (1 + param_delta) * tmp_, ca_, params))
            gm = float(growth_primary_tipping(max(t, 1), cov_, (1 - param_delta) * tmp_, ca_, params))
            SCC_4p[t] = cov_ * (gp - gm) / (2 * param_delta * tmp_)
        PSI = SCC_4p

        # ---- ENVELOPE derivatives, two directions ---------------------------
        #   joint     : A and the carbon it holds move together (the SCD object)
        #   stateonly : A moves at unchanged S, nothing released (sign flips)
        _A  = np.asarray(u[:, [2]], dtype=float)
        _Tm = np.asarray(u[:, [1]], dtype=float)
        _s1 = np.asarray(shock_save1, dtype=float).reshape(-1, 1)[:_A.shape[0]]
        dA_step = param_delta * np.abs(_A) + 1e-8
        _dS = _s1 * params["carbonstock_amz"] * dA_step

        uJp = np.array(u, dtype=float, copy=True); uJm = np.array(u, dtype=float, copy=True)
        uJp[:, [2]] = np.clip(_A + dA_step, 1e-6, 1.0); uJp[:, [1]] = _Tm - _dS
        uJm[:, [2]] = np.clip(_A - dA_step, 1e-6, 1.0); uJm[:, [1]] = _Tm + _dS
        uSp = np.array(u, dtype=float, copy=True); uSm = np.array(u, dtype=float, copy=True)
        uSp[:, [2]] = np.clip(_A + dA_step, 1e-6, 1.0)
        uSm[:, [2]] = np.clip(_A - dA_step, 1e-6, 1.0)
        dA_eff = uJp[:, [2]] - uJm[:, [2]]

        def _ev(u_step, t):
            return EV_simulation_stochastic2(
                t_step=t + 1, x_step=x[[t + 1], :], u_step=u_step[[t + 1], :],
                coefn=coef_V[[t + 1], :], weight_cheb=params["weight_cheb"],
                u_min=params["u_min"][[t + 1], :], u_max=params["u_max"][[t + 1], :],
                params=params, shock1=shock_save1[t + 1], shock2=shock_save2[t + 1])

        dUdA_joint = np.zeros((params["T"] - 1, 1))
        dUdA_state = np.zeros((params["T"] - 1, 1))
        for t in range(0, params["T"] - 1):
            _d = float(dA_eff[t + 1])
            if abs(_d) > 1e-14:
                dUdA_joint[t] = (float(np.asarray(_ev(uJp, t)).ravel()[0])
                                 - float(np.asarray(_ev(uJm, t)).ravel()[0])) / _d
                dUdA_state[t] = (float(np.asarray(_ev(uSp, t)).ravel()[0])
                                 - float(np.asarray(_ev(uSm, t)).ravel()[0])) / _d
        SCD_env_joint = (1.0 / params["carbonstock_amz"]) * monetary * dUdA_joint
        SCD_env_state = (1.0 / params["carbonstock_amz"]) * monetary * dUdA_state

        # ---- SCD -------------------------------------------------------------
        _climate_arg = np.asarray(shock_save2).reshape(-1)
        _shock1_vec  = np.asarray(shock_save1).reshape(-1)
        scd = compute_scd_block_v3(u, x, params, monetary, marginal_utility_temp,
                                   SCC_3, param_delta, _climate_arg, _shock1_vec,
                                   verbose=(draw == 0))
        SCD_2 = scd["Phi_bio"]     # 1+g+A*nu_A, no Theta: the carbon->T feedback is
                                   # carried explicitly by the loop channel here.
        EST = scd["estimators"]

        # ---- (I) as an explicit sum -------------------------------------------

        _nS1 = len(shock_save1)
        standard_T = np.zeros((params["T"] - 1, 1))
        for t in range(params["T"] - 1):
            max_k = len(marginal_utility_temp) - t - 1
            sk = 0.0
            for k in range(0, max_k):
                s1k = float(shock_save1[t + k]) if (t + k) < _nS1 else float(shock_save1[-1])
                sk += params["bet"] ** k * marginal_utility_temp[t + k + 1] * s1k
            standard_T[t] = sk

        # ---- (II), stock kernel ------------------------------------------------
        def _chanII(kernel):
            out = np.zeros((params["T"] - 1, 1))
            for t in range(params["T"] - 1):
                max_k = len(marginal_utility_temp) - t - 1
                d = 1.0; run = 0.0; tot = 0.0
                for k in range(1, max_k):
                    i = t + k - 1
                    d = d * (float(SCD_2[i]) if i < params["T"] - 1 else 1.0)
                    c_ = float(SCC_3[t + k]) if (t + k) < len(SCC_3) else 0.0
                    if kernel == "stock":
                        tot += (params["bet"] ** k) * marginal_utility_temp[t + k + 1] * c_ * d
                    else:
                        run += c_ * d
                        tot += (params["bet"] ** k) * marginal_utility_temp[t + k + 1] * run
                out[t] = tot * float(PSI[t]) * float(shock_save1[t])
            return out

        standard_A      = _chanII("stock")
        standard_A_flow = _chanII("flow")


        # ---- prices : ENVELOPPE (niveau) et EXPLICITE (decomposition) ---------
        SCC_fossil    = - monetary * (SCC_2)                    # dV/dS : niveau, GE+Amazonie inclus
        SCC_explicit  = - monetary * (standard_T)              # (I) SCC nu, equilibre partiel
        chanII_money  = - monetary * (standard_A)              # (II) sous-systeme
 
        # (I)xboucle, fold j=0 : EXPLICITE (standard_T), PAS SCC_2 (enveloppe),
        # sinon le canal loop est contamine par le GE.
        loop_j0_expl  = - monetary * (standard_T * SCC_3 * PSI)
 
        # SCCDS explicite : (I) + (II) + fold j=0 de la boucle.  La version pleine
        # ajoute (III_T + IV_T) au stage cross-draw (colonnes 10-13).
        SCCDS_expl_compact = SCC_explicit + loop_j0_expl + chanII_money
 
        # NIVEAU reporte = enveloppe pure (V est 3D -> contient deja la boucle).
        SCCDS = SCC_fossil
        SCCDS_compact = SCCDS_expl_compact

        # ---- uplift EXPLICITE (decomposable) ---------------------------------
        _e0 = np.where(np.abs(SCC_explicit) < 1e-9, np.nan, SCC_explicit)
        uplift_expl_compact = (SCCDS_expl_compact - SCC_explicit) / _e0
        #  uplift plein = (compact + III_T+IV_T)/SCC_explicit : a former a l'agreg,
        #  ou reutilise la colonne loop pleine deja calculee, divisee par SCC_explicit.
# -----------------------------------------------------------------------------
 
# ---- uplift ENVELOPPE (3D vs run0004 gele) : PAR DRAW, vecteur complet ----
        uplift_env_vec = np.full((params["T"] - 1, 1), np.nan)
        try:
            _rf2 = comparative_folder
            coef_V_2D = np.load(_rf2 + 'coef_V.npy')
            umin2 = np.genfromtxt(_rf2 + '/u_min.csv', delimiter=';')
            umax2 = np.genfromtxt(_rf2 + '/u_max.csv', delimiter=';')
            def _ev2(u_step, t):
                return EV_simulation_stochastic2(
                    t_step=t + 1, x_step=x[[t + 1], :], u_step=u_step[[t + 1], :],
                    coefn=coef_V_2D[[t + 1], :], weight_cheb=params["weight_cheb"],
                    u_min=umin2[[t + 1], :], u_max=umax2[[t + 1], :],
                    params=params, shock1=shock_save1[t + 1], shock2=shock_save2[t + 1])
            SCC_1_2D = np.zeros((params["T"] - 1, 1))
            for t in range(0, params["T"] - 1):
                SCC_1_2D[t] = (float(np.asarray(_ev2(uplus, t)).ravel()[0])
                               - float(np.asarray(_ev2(uminus, t)).ravel()[0])) / (2 * float(dT[t + 1]))
            SCC_2_2D      = SCC_1_2D * np.asarray(shock_save1).reshape(-1, 1)[:len(SCC_1_2D)]
            SCC_fossil_2D = - monetary * SCC_2_2D
            _f2v = np.where(np.abs(SCC_fossil_2D) < 1e-9, np.nan, SCC_fossil_2D)
            uplift_env_vec = 100.0 * (SCC_fossil - SCC_fossil_2D) / _f2v      # en %, par periode
        except Exception as _ex:
            if draw == 0:
                print("  [uplift_env] indisponible:", _ex)

        # ---- columns --------------------------------------------------------------
        _T2 = params["T"] - 2
        _s  = np.where(np.abs(SCCDS[0:_T2]) < 1e-9, np.nan, SCCDS[0:_T2])
        _f  = np.where(np.abs(SCC_fossil[0:_T2]) < 1e-9, np.nan, SCC_fossil[0:_T2])
        _e  = np.where(np.abs(SCC_explicit[0:_T2]) < 1e-9, np.nan, SCC_explicit[0:_T2])
        _rv = EST["rev_M_stock"]
        _rvd = np.where(np.abs(_rv) < 1e-12, np.nan, _rv)

        final_doc[:, [4]]  = EST["perm_M_stock"]
        final_doc[:, [5]]  = scd["SCD_T"]
        final_doc[:, [6]]  = scd["SCD_A_M"]
        final_doc[:, [7]]  = SCCDS[0:_T2]
        final_doc[:, [8]]  = standard_T[0:_T2]
        final_doc[:, [9]]  = standard_A[0:_T2]
        final_doc[:, [14]] = EST["rev_M_stock"]
        final_doc[:, [15]] = EST["perm_N_stock"]
        final_doc[:, [16]] = scd["zeta"] * np.ones((_T2, 1))
        # col 17 : credit d'exoneration, denominateur HOMOGENE (SCC explicite)
        final_doc[:, [18]] = SCC_fossil[0:_T2]              # NIVEAU : SCC enveloppe (GE inclus)
        # col 19 : ratio-titre SCD_M / SCC explicite (lu par X3/X4)
        final_doc[:, [19]] = EST["perm_M_stock"] / _e
        final_doc[:, [20]] = scd["g_net"][0:_T2]
        final_doc[:, [21]] = scd["Phi_M"][0:_T2]
        final_doc[:, [22]] = SCC_explicit[0:_T2]
        # col 23 : 2e convention centrale, SCD_N / SCC explicite
        final_doc[:, [23]] = EST["perm_N_stock"] / _e
        final_doc[:, [24]] = chanII_money[0:_T2]                        # (II)
        final_doc[:, [25]] = loop_j0_expl[0:_T2]                        # fold j=0, explicite
        final_doc[:, [26]] = SCCDS_expl_compact[0:_T2]                  # SCCDS explicite compact
        # uplift explicite compact, en % du SCC nu :
        final_doc[:, [17]] = 100.0 * uplift_expl_compact[0:_T2]         # <- uplift, pas credit

        final_doc[:, [27]] = np.asarray(monetary).reshape(-1, 1)[0:_T2]
        final_doc[:, [28]] = SCD_env_joint[0:_T2]
        final_doc[:, [29]] = SCD_env_state[0:_T2]
        final_doc[:, [30]] = SCD_env_joint[0:_T2] / _rvd
        final_doc[:, [31]] = SCD_env_joint[0:_T2] / _f      # ROBUSTESSE : ratio enveloppe/enveloppe
        # the conventions that bracket unity
        final_doc[:, [32]] = EST["perm_N_spaceSuppr_stock"]
        final_doc[:, [33]] = EST["rev_N_spaceSuppr_stock"]
        final_doc[:, [34]] = EST["perm_M_spaceSuppr_stock"]
        final_doc[:, [35]] = EST["perm_N_climOnly_stock"]
        final_doc[:, [36]] = uplift_env_vec[0:_T2]        # uplift enveloppe %, par draw -> pooled dans le CSV

        # ---- loop channel, MONETISED, for the cross-draw stage -----------------------
        # -monetary[t] is folded into L1T so III_T + IV_T comes out in SCCDS units.
        L1T_3D = np.zeros((params["T"] - 1, params["T"] - 1, params["T"] - 1))
        L2T_3D = np.zeros((params["T"] - 1, params["T"] - 1, params["T"] - 1))
        for t in range(params["T"] - 1):
            max_k = len(marginal_utility_temp) - t - 1
            m_t = - float(np.asarray(monetary).reshape(-1)[t])
            for k in range(1, max_k):
                for j in range(0, k):
                    if (t + k + 1) < len(SCC_3) and (t + j + 1) < len(PSI):
                        L1T_3D[t, k, j] = m_t * marginal_utility_temp[t + k + 1] * SCC_3[t + k + 1]
                        prod_AA = np.prod(SCD_2[t + j + 1:t + k]) if (k - 1) >= (j + 1) else 1.0
                        L2T_3D[t, k, j] = prod_AA * PSI[t + j + 1] * shock_save1[t]

        # ---- subsystem channel, MONETISED, like the T channel ------------------------
        # -monetary[t] is folded into L1A so that III_A + IV_A comes out in SCCDS
        # units, exactly as III_T + IV_T do. X1 therefore tests the identity against
        # column 24 (chanII_money), not against the raw column 9.
        # Stock kernel: ONE term per k at the j = 1 slot, running product identical
        # to _chanII("stock"), which is what makes the identity exact.
        L1A_3D = np.zeros((params["T"] - 1, params["T"] - 1, params["T"] - 1))
        L2A_3D = np.zeros((params["T"] - 1, params["T"] - 1, params["T"] - 1))
        for t in range(params["T"] - 1):
            max_k = len(marginal_utility_temp) - t - 1
            m_t = - float(np.asarray(monetary).reshape(-1)[t])
            d = 1.0
            for k in range(1, max_k):
                i = t + k - 1
                d = d * (float(SCD_2[i]) if i < params["T"] - 1 else 1.0)
                if k < params["T"] - 1 and (t + k) < len(SCC_3):
                    L1A_3D[t, k, 1] = m_t * marginal_utility_temp[t + k + 1] * SCC_3[t + k]
                    L2A_3D[t, k, 1] = d * PSI[t] * shock_save1[t]

    return final_doc, L1T_3D, L2T_3D, L1A_3D, L2A_3D


if __name__ == '__main__':
    BATCH_SIZE = params["CPU"]
#    temp_folder = os.path.join(run_folder, "temp_results")
    os.makedirs(temp_folder, exist_ok=True)

    existing_batches = sorted(
        f for f in os.listdir(temp_folder)
        if f.startswith("batch_") and f.endswith(".npz"))

    if existing_batches:
        batch_ranges = []
        for f in existing_batches:
            name = f.replace(".npz", "")
            _, start, end = name.split("_")
            batch_ranges.append((int(start), int(end)))
        last_done = max(end for (_, end) in batch_ranges)
        start_from = min(last_done + 1, params["draws"])
    else:
        start_from = 0

    zipped = [(draw, params, coef_V) for draw in range(start_from, params["draws"])]

    print(f"[RUN] cache = {temp_folder}")
    print(f"[RUN] existing batches = {len(existing_batches)}   start_from = {start_from}"
          f"   draws = {params['draws']}   to compute = {len(zipped)}")
    if len(zipped) == 0 and len(existing_batches) > 0:
        print("[RUN] WARNING: nothing to compute, the cache already covers every draw.")
        print("[RUN] Results below come entirely from cached .npz files, which may")
        print("[RUN] have been produced by an earlier version of the code.")
        print(f"[RUN] To force a recompute:  rm -f {temp_folder}/batch_*.npz")

    master_seed = 12345
    seeds = np.random.SeedSequence(master_seed).spawn(params["draws"])

    batch = []
    batch_indices = []

    with multiprocess.Pool(processes=min(64, params["draws"])) as pool:
        valid_runs = 0
        invalid_runs = 0

        for i, res in enumerate(pool.imap_unordered(task, zipped), start=start_from):
            if res is None:
                invalid_runs += 1
                continue
            valid_runs += 1
            batch.append(res)
            batch_indices.append(i)
            if len(batch) >= BATCH_SIZE or i == params["draws"] - 1:
                batch_file = os.path.join(
                    temp_folder, f"batch_{batch_indices[0]}_{batch_indices[-1]}.npz")
                np.savez_compressed(batch_file, *[np.array(r, dtype=object) for r in batch])
                batch.clear()
                batch_indices.clear()

    if len(batch) > 0:
        batch_file = os.path.join(temp_folder, f"batch_{batch_indices[0]}_{batch_indices[-1]}.npz")
        np.savez_compressed(batch_file, *[np.array(r, dtype=object) for r in batch])
        batch.clear()
        batch_indices.clear()

    # ---------------------------------------------------------------- pooling
    batch_files = sorted(
        os.path.join(temp_folder, f)
        for f in os.listdir(temp_folder)
        if f.endswith(".npz"))

    MAX_RESERVOIR = 500
    mean = None
    min_val = None
    max_val = None
    reservoir = []
    n = 0
    filtered_runs = 0
    mean_L1T = mean_L2T = M2_T = None
    mean_L1A = mean_L2A = M2_A = None
    n_cov = 0
    Tm1 = K = None

    for bf in batch_files:
        with np.load(bf, allow_pickle=True) as data:
            for key in data.files:
                res = data[key].tolist()
                X = np.asarray(res[0], dtype=float)          # (T-2, ncol)
                if np.any(X[0:20, 0] <= 0.05):
                    filtered_runs += 1
                    continue
                if mean is None:
                    mean = np.zeros_like(X)
                    min_val = np.full_like(X, np.inf)
                    max_val = np.full_like(X, -np.inf)
                n += 1
                mean += (X - mean) / n
                min_val = np.minimum(min_val, X)
                max_val = np.maximum(max_val, X)
                if len(reservoir) < MAX_RESERVOIR:
                    reservoir.append(X.copy())
                else:
                    j = np.random.randint(0, n)
                    if j < MAX_RESERVOIR:
                        reservoir[j] = X.copy()

                L1T = np.asarray(res[1], dtype=float)
                L2T = np.asarray(res[2], dtype=float)
                L1A = np.asarray(res[3], dtype=float)
                L2A = np.asarray(res[4], dtype=float)
                if mean_L1T is None:
                    Tm1, K, _ = L1T.shape
                    mean_L1T = np.zeros((Tm1, K, K)); mean_L2T = np.zeros((Tm1, K, K))
                    M2_T = np.zeros((Tm1, K, K))
                    mean_L1A = np.zeros((Tm1, K, K)); mean_L2A = np.zeros((Tm1, K, K))
                    M2_A = np.zeros((Tm1, K, K))
                n_cov += 1
                d1 = L1T - mean_L1T; mean_L1T += d1 / n_cov
                d2 = L2T - mean_L2T; mean_L2T += d2 / n_cov
                M2_T += d1 * (L2T - mean_L2T)
                d1 = L1A - mean_L1A; mean_L1A += d1 / n_cov
                d2 = L2A - mean_L2A; mean_L2A += d2 / n_cov
                M2_A += d1 * (L2A - mean_L2A)

    if n == 0:
        raise RuntimeError(f"no usable draw. batches={len(batch_files)}, "
                           f"filtered={filtered_runs}. Check {temp_folder}.")

    reservoir = np.stack(reservoir, axis=0)
    final_mean = mean
    final_min = min_val
    final_max = max_val
    final_05 = np.percentile(reservoir, 5, axis=0)
    final_95 = np.percentile(reservoir, 95, axis=0)

    beta = params["bet"]
    cov_T = M2_T / max(n_cov - 1, 1)
    cov_A = M2_A / max(n_cov - 1, 1)
    _T2 = params["T"] - 2

    IV_T = np.zeros(Tm1); IV_A = np.zeros(Tm1)
    III_T = np.zeros(Tm1); III_A = np.zeros(Tm1)
    for t in range(Tm1):
        for k in range(1, K):
            for j in range(0, k):
                IV_T[t]  += (beta ** k) * cov_T[t, k, j]
                III_T[t] += (beta ** k) * mean_L1T[t, k, j] * mean_L2T[t, k, j]
            for j in range(1, k + 1):
                IV_A[t]  += (beta ** k) * cov_A[t, k, j]
                III_A[t] += (beta ** k) * mean_L1A[t, k, j] * mean_L2A[t, k, j]

    final_mean[:, [10]] = III_T[:_T2].reshape(-1, 1)
    final_mean[:, [11]] = III_A[:_T2].reshape(-1, 1)
    final_mean[:, [12]] = IV_T[:_T2].reshape(-1, 1)
    final_mean[:, [13]] = IV_A[:_T2].reshape(-1, 1)

    _IIIA = III_A[:_T2].ravel(); _IVA = IV_A[:_T2].ravel()
    _IIIT = III_T[:_T2].ravel(); _IVT = IV_T[:_T2].ravel()

    # =================================================================== Y1
    # ===================================================================
    #  FULL SCCDS. L1T and L1A both carry -monetary (folded in task), so all of
    #  III_T, IV_T, III_A, IV_A are already in SCCDS units. The loop channel
    #  REPLACES the j = 0 fold, hence column 26 ((I) + (II)) rather than column 7.
# -----------------------------------------------------------------------------
 
    _loop_full  = _IIIT + _IVT
    _sccds_part = final_mean[:_T2, 26]
    SCCDS_full  = _sccds_part + _loop_full
    _foss       = final_mean[:_T2, 18]
    _den        = np.where(np.abs(_foss) < 1e-12, np.nan, _foss)

    print("\n" + "=" * 78)
    print(f"  RESULTS   (n = {n} draws pooled, {filtered_runs} filtered out, "
          f"reservoir = {reservoir.shape[0]})")
    print("=" * 78)

    print("\n Y1  SCCDS : compact (j=0 fold only) vs full loop channel")
    print("    t   |  SCC_foss  |  compact  |    full   | uplift compact | uplift full")
    for tt in [0, 5, 10, 20, 30, 40, 60]:
        if tt < _T2:
            uc = 100 * (final_mean[tt, 7] - _foss[tt]) / _den[tt]
            uf = 100 * (SCCDS_full[tt]    - _foss[tt]) / _den[tt]
            print(f"  {tt:5d} | {_foss[tt]:10.2f} | {final_mean[tt,7]:9.2f} | "
                  f"{SCCDS_full[tt]:9.2f} | {uc:+13.4f}% | {uf:+10.4f}%")
    print(f"    mean uplift compact = "
          f"{100*np.nanmean((final_mean[:_T2,7]-_foss)/_den):+.4f}%")
    print(f"    mean uplift full    = "
          f"{100*np.nanmean((SCCDS_full-_foss)/_den):+.4f}%")
    print("    -> the gap between the two columns is the contribution of the")
    print("       intermediate-date loop terms, i.e. (III)+(IV) of the paper.")

    # =================================================================== Y2
    print("\n Y2  SCD against the full SCCDS")
    _sf = np.where(np.abs(SCCDS_full) < 1e-12, np.nan, SCCDS_full)
    print("    t   |  dVdA/SCCDS | SCD_M/SCCDS | SCD_N/SCCDS | SCD_M/SCC_foss")
    for tt in [0, 5, 10, 20, 30, 40]:
        if tt < _T2:
            print(f"  {tt:5d} | {final_mean[tt,14]/_sf[tt]:11.4f} | "
                  f"{final_mean[tt,4]/_sf[tt]:11.4f} | "
                  f"{final_mean[tt,15]/_sf[tt]:11.4f} | {final_mean[tt,19]:14.4f}")
    print(f"    share of periods with SCD_M > SCCDS_full : "
          f"{np.nanmean(final_mean[:_T2,4] > _sf):.1%}")
    print(f"    share of periods with SCD_N > SCCDS_full : "
          f"{np.nanmean(final_mean[:_T2,15] > _sf):.1%}")

#    # =================================================================== X1
#    _II  = final_mean[:_T2, 9]
#    _rec = _IIIA + _IVA
#    _dII = np.where(np.abs(_II) < 1e-15, np.nan, _II)
#    print("\n X1  IDENTITY   (III)_A + (IV)_A  ==  (II) = standard_A  [column 9, raw]")
#    print(f"    max|III_A+IV_A - II| = {np.nanmax(np.abs(_rec - _II)):.4e}")
#    print(f"    max relative error   = {np.nanmax(np.abs((_rec - _II) / _dII)):.4e}")
#    for tt in [0, 5, 10, 20, 40]:
#        if tt < _T2:
#            print(f"    t={tt:3d}: II={_II[tt]:+.6e}  III_A={_IIIA[tt]:+.6e}  "
#                  f"IV_A={_IVA[tt]:+.6e}  sum={_rec[tt]:+.6e}")
#    print("    -> III_A/IV_A are the expectation/insurance SPLIT of (II) and must")
#    print("       NOT be added on top of it in the paper equation.")

    # Les quatre composantes 10-13 sont maintenant toutes en $/tC, donc
    # l'identit� se teste contre la colonne 24 (chanII mon�tis�), pas la 9.
    _II  = final_mean[:_T2, 24]
    _rec = _IIIA + _IVA
    _dII = np.where(np.abs(_II) < 1e-12, np.nan, _II)
    print("\n X1  IDENTITY   (II)^E + (II)^cov  ==  (II)   [column 24, $/tC]")
    print(f"    max|III_A+IV_A - II| = {np.nanmax(np.abs(_rec - _II)):.4e}")
    print(f"    max relative error   = {np.nanmax(np.abs((_rec - _II) / _dII)):.4e}")
    for tt in [0, 5, 10, 20, 40]:
        if tt < _T2:
            print(f"    t={tt:3d}: II={_II[tt]:+.6e}  III_A={_IIIA[tt]:+.6e}  "
                  f"IV_A={_IVA[tt]:+.6e}  sum={_rec[tt]:+.6e}")
    print("    -> III_A/IV_A are the expectation/covariance SPLIT of (II) and must")
    print("       NOT be added on top of it in the paper equation.")
    print("    -> all four of columns 10-13 are now in $/tC.")


    print("\n X2  INSURANCE SHARES")
    for tt in [0, 5, 10, 20, 40]:
        if tt < _T2:
            sh  = _IVA[tt] / _dII[tt] if np.isfinite(_dII[tt]) else np.nan
            dT_ = _IIIT[tt] + _IVT[tt]
            shT = _IVT[tt] / dT_ if abs(dT_) > 1e-15 else np.nan
            print(f"    t={tt:3d}: IV_A/(II) = {sh:+.4f}    IV_T/(III_T+IV_T) = {shT:+.4f}")
    print("    (only meaningful once X1 passes)")

#    print("\n X3  SCD DISTRIBUTION ACROSS DRAWS   (col19 = SCD_M / SCC_fossil)")
#    _r = reservoir[:, :, 19]
#    for tt in [0, 5, 10, 20, 40]:
#        if tt < _r.shape[1]:
#            v = _r[:, tt]; v = v[np.isfinite(v)]
#            if v.size:
#                print(f"    t={tt:3d}: mean={v.mean():8.3f} median={np.median(v):8.3f} "
#                      f"p5={np.percentile(v,5):8.3f} p95={np.percentile(v,95):8.3f} "
#                      f"max={v.max():9.3f}  share>1={np.mean(v>1):5.1%}  n={v.size}")
#    print("    -> quote the median and [p5,p95]; the mean is tail-driven.")


    print("\n X3  DISTRIBUTION DU SCD SUR LES DRAWS")
    print("     titre   col19 = SCD_M / SCC_explicite  (homogenous, Prop.4 exacte au gel)")
    print("     robuste col31 = SCD_env / SCC_enveloppe (GE inclus, plus bruit)")
    for tag, cc in [("expl", 19), ("env ", 31)]:
        _r = reservoir[:, :, cc]
        for tt in [0, 5, 10, 20, 40]:
            if tt < _r.shape[1]:
                v = _r[:, tt]; v = v[np.isfinite(v)]
                if v.size:
                    print(f"    [{tag}] t={tt:3d}: median={np.median(v):7.3f} "
                        f"p5={np.percentile(v,5):7.3f} p95={np.percentile(v,95):7.3f} "
                        f"share>1={np.mean(v>1):5.1%}")
    print("    -> cite la mdiane et [p5,p95] ; la moyenne est tire par la queue.")


    print("\n X4  DISPERSION at t0 (reservoir is a random subsample, not a path)")
    v = reservoir[:, 0, 19]; v = v[np.isfinite(v)]
    if v.size:
        print(f"    n={v.size}  mean={v.mean():.4f}  median={np.median(v):.4f}  "
              f"sd={v.std():.4f}  se={v.std()/np.sqrt(v.size):.4f}  "
              f"[p5,p95]=[{np.percentile(v,5):.3f},{np.percentile(v,95):.3f}]")

    print("\n X5  ORDERING ACROSS DRAWS   dVdA <= SCD_M <= SCD_N")
    d_, m_, n_ = reservoir[:, :, 14], reservoir[:, :, 4], reservoir[:, :, 15]
    msk = np.abs(m_) > 1e-9
    print(f"    share(dVdA<=SCD_M)  = {np.mean(d_[msk] <= m_[msk] + 1e-6):.4%}")
    print(f"    share(SCD_M<=SCD_N) = {np.mean(m_[msk] <= n_[msk] + 1e-6):.4%}")
    gg  = reservoir[:, :, 20]
    bad = (m_ > n_ + 1e-6) & msk
    if bad.any():
        print(f"    among {int(bad.sum())} violations, share with g>0 = "
              f"{np.mean(gg[bad] > 0):.4%}   (Phi_N = Phi_M - g)")

    print("\n X6  PHI, g, zeta across draws")
    ph = reservoir[:, :, 21]
    pre = ph[:, :37]
    print(f"    Phi_M: min={np.nanmin(ph):.4f} max={np.nanmax(ph):.4f} "
          f"share>=1 = {np.nanmean(ph >= 1.0):.2%}  (pre-freeze only: "
          f"{np.nanmean(pre >= 1.0):.2%})")
    print("      NB: Phi is forced to 1 after the structural freeze, so the raw")
    print("          share>=1 mostly measures the frozen fraction of the horizon.")
    print(f"    g    : min={np.nanmin(gg):+.5f} max={np.nanmax(gg):+.5f} "
          f"share>0 = {np.nanmean(gg > 0):.2%}")
    z = reservoir[:, :, 16]
    print(f"    zeta : min={np.nanmin(z):.5f} max={np.nanmax(z):.5f} "
          f"within-draw std = {np.nanmax(np.nanstd(z, axis=1)):.3e}  (must be ~0)")

    print("\n X7  SAMPLE")
    print(f"    valid_runs = {valid_runs}   returned None = {invalid_runs}   "
          f"filtered (X[0:20,0]<=0.05) = {filtered_runs}   pooled = {n}")
    if filtered_runs > 0:
        print("    -> check the filtered runs are not precisely the severe-dieback")
        print("       draws, otherwise the tail carrying the result is truncated.")

    print("\n X8  DENOMINATORS")
    _expl = final_mean[:_T2, 22]
    _de   = np.where(np.abs(_expl) < 1e-12, np.nan, _expl)
    print(f"    SCC envelope t0 = {_foss[0]:.3f}   SCC explicit sum t0 = {_expl[0]:.3f}"
          f"   ratio = {_foss[0]/_de[0]:.4f}")
    print(f"    mean ratio envelope/explicit = {np.nanmean(_foss/_de):.4f}")
    print(f"    SCD_M/SCC_explicit t0 = {final_mean[0,23]:.4f}   "
          f"vs /SCC_envelope = {final_mean[0,19]:.4f}")
    print("=" * 78 + "\n")

    # ---------------------------------------------------------------- outputs
    final95  = np.column_stack((final_95[:, 0], final_95[:, 2], final_95[:, 3]))
    final5   = np.column_stack((final_05[:, 0], final_05[:, 2], final_05[:, 3]))
    finalmin = np.column_stack((final_min[:, 0], final_min[:, 2], final_min[:, 3]))
    finalmax = np.column_stack((final_max[:, 0], final_max[:, 2], final_max[:, 3]))

    np.savetxt(run_folder + 'outputs_stochastic95.csv', final95, delimiter=';')
    np.savetxt(run_folder + 'outputs_stochastic5.csv', final5, delimiter=';')
    np.savetxt(run_folder + 'outputs_stochastic.csv', final_mean, delimiter=';')
    np.savetxt(run_folder + 'outputs_stochastic_max.csv', finalmax, delimiter=';')
    np.savetxt(run_folder + 'outputs_stochastic_min.csv', finalmin, delimiter=';')
    np.savetxt(run_folder + 'outputs_sccds_full.csv',
               np.column_stack((_foss, final_mean[:_T2, 7], SCCDS_full,
                                _IIIT, _IVT, _IIIA, _IVA)), delimiter=';')

    print(f"valid_runs {valid_runs}   pooled {n}")
