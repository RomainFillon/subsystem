# -*- coding: utf-8 -*-

#import some libraries
import numpy as np
from numpy import kron, tile
import functools
import scipy.integrate as si
from scipy.stats import truncnorm, beta, qmc
import math
import csv
from numpy.polynomial.hermite import hermgauss
from numpy.polynomial.legendre import leggauss


#Compute the expected continuation value at time t+1 under climate uncertainty using a two-dimensional Gaussian quadrature.
#dT/dS is integrated using Gauss�Hermite quadrature
#dA/dT is integrated using Gauss�Legendre quadrature on the unit interval. 
#The function evaluates the expected value of the Bellman equation by propagating the state variables
#through the stochastic law of motion and approximating the continuation value using Chebyshev polynomials.
def EV_quadrature_stochastic2(
    t_step, x_step, u_step,
    coefn, weight_cheb,
    u_min, u_max,
    params,
    N1=5,   # TCRE (normal)
    N2=5):
    # ---------- TCRE : Hermite ----------
    z1, w1 = hermgauss(N1)
    z1 = np.sqrt(2.0) * z1
    w1 = w1 / np.sqrt(np.pi)
    # ---------- S : [0,1] ----------
    z2, w2 = leggauss(N2)
    z2 = 0.5 * (z2 + 1.0)
    w2 = 0.5 * w2
    # ---------- bornes Chebyshev ----------
    u_min_adj = adjust_min(u_min, u_max, m=params["degb"], max_frac=0.1)
    u_max_adj = adjust_max(u_min, u_max, m=params["degb"], max_frac=0.1)
    EV = 0.0
    for i in range(N1):
        shock1 = params["TCRE_mean"] + params["sdeviation"] * z1[i]
        for j in range(N2):
            shock2 = z2[j]
            # law motion
            un = law_motion_SEU(
                t_step, x_step, u_step,
                shock1, shock2, 0, params)
            un_adj = np.clip(un, u_min_adj, u_max_adj)
            V_next = approx_Chebyshev(
                coefn, weight_cheb,
                un_adj, u_min_adj, u_max_adj)
            Vval = (
                maximand(np.array([t_step - 1]), x_step, u_step, params)
                + params["bet"] * V_next)
            if not np.isfinite(Vval).all():
                return np.nan
            EV += w1[i] * w2[j] * Vval
    return EV

#Compute the expected terminal continuation value under stochastic climate
#dynamics using two-dimensional Gaussian quadrature.
def quadrature_product_terminal_EV(
    t_step, x_step, u_step,
    params,
    N1=5,
    N2=5):
    z1, w1 = hermgauss(N1)
    z1 = np.sqrt(2) * z1
    w1 = w1 / np.sqrt(np.pi)
    z2, w2 = leggauss(N2)
    z2 = 0.5 * (z2 + 1.0)
    w2 = 0.5 * w2
    EV = 0.0
    for i in range(N1):
        shock1 = params["TCRE_mean"] + params["sdeviation"] * z1[i]
        for j in range(N2):
            shock2 = z2[j]
            un = law_motion_SEU(
                t_step, x_step, u_step,
                shock1, shock2, 0, params
            )
            Vval = (
                maximand(np.array([t_step - 1]), x_step, u_step, params)
                + params["bet"] * terminal_value(t_step + 1, un, params)
            )
            if not np.isfinite(Vval).all():
                return np.nan
            EV += w1[i] * w2[j] * Vval
    return EV

#terminal value function 
#for stochastic 1
def integrand_terminal1(shock1, t_step, x_step, u_step, params):
    # next-stage state variable should be in the domain of the approximation of next-stage value function
    un = law_motion_SEU(t_step, x_step, u_step, shock1,params["climate_meanEU"],0, params)
    v1 = maximand(np.array([t_step - 1]), x_step, u_step, params) + params["bet"] * terminal_value(t_step + 1,un, params)
    v2=truncnorm.pdf(shock1, params["a_tcre"], params["b_tcre"], loc=params["TCRE_mean"], scale=params["sdeviation"])
    v = v1 * v2
    return v;
    return v;

#compute interval for state variables
def compute_intervals(params):
    time_horizon = np.arange(params["T"]+2).reshape((params["T"]+2,1))

    if params["dim"] == 1:
        strip_width = 0.3
        ustat = (params["alpha"]*params["deltaT"]/(params["delta"]*np.ones((params["T"]+2,1)) + params["rho"]*np.ones((params["T"]+2,1)) + params["theta"]*params["deltaT"]*(params["GA"]/(1-params["alpha"]))))**(1/(1-params["alpha"]))
        u0trend = params["K0"] / (params["A0"]**(1/(1-params["alpha"])) * params["POP0"]) * np.ones((params["T"]+2,1))
        mix_param_min = 0.95
        mix_param_max = 0.95
        ucentered_min = ustat*(np.ones((params["T"]+2,1)) - (mix_param_min**(params["deltaT"]*time_horizon))) + u0trend*(mix_param_min**(params["deltaT"]*time_horizon))
        ucentered_max = ustat*(np.ones((params["T"]+2,1)) - (mix_param_max**(params["deltaT"]*time_horizon))) + u0trend*(mix_param_max**(params["deltaT"]*time_horizon))
        u_min = (1 - strip_width) * ucentered_min[0:params["T"],:]
        u_max = (1 + strip_width) * ucentered_max[0:params["T"],:]

    elif params["dim"] == 3:
        u_m = np.genfromtxt(params["prev_run_folder"]+'state_V_notstochastic.csv', dtype=float, delimiter=';')
        if u_m.shape == (params["T"],):
            k_ramsey = (params['A'][0:params["T"],:]**(1/(1-params["alpha"])))*params['L'][0:params["T"],:]*u_m.reshape((params["T"],1))
            k_min = 0.8*k_ramsey
            kpec_min = 0.8*u_m.reshape((params["T"],1))
            k_max = 1.2*k_ramsey
            kpec_max = 1.2*u_m.reshape((params["T"],1))
        else:
            k_ramsey = (params['A'][0:params["T"],:]**(1/(1-params["alpha"])))*params['L'][0:params["T"],:]*u_m[:,0].reshape((params["T"],1))
            k_min = 0.8*k_ramsey
            kpec_min = 0.8*u_m[:,0].reshape((params["T"],1))
            k_max = 1.2*k_ramsey
            kpec_max = 1.2*u_m[:,0].reshape((params["T"],1))
        emissions_max = np.zeros((params["T"],1))
        emissions_min = np.zeros((params["T"],1))
        ab_pess = 0.8*np.ones((params["T"],1))
        for j in range(1,params["T"]+1):
            emissions_max[j-1,:] = emissions(ab_pess[j-1,:], production(params['A'][j-1,:],k_max[j-1,:],params['L'][[j-1],:], params), params['sigm'][[j-1],:])
            emissions_min[j-1,:] = emissions(0.8, production(params['A'][j-1,:],k_min[j-1,:],params['L'][[j-1],:], params), params['sigm'][[j-1],:])

        s_max = 3*(params["T0"]/params["TCRE_mean"])*np.ones((params["T"],1))
        s_min = 0.9*(params["T0"]/params["TCRE_mean"])*np.ones((params["T"],1))
        for j in range(1,params["T"]):
            s_max[j,:] = s_max[j-1,:]+emissions_max[j-1,:]
        t_min = params["TCRE_mean"]*s_min
        t_max = params["TCRE_mean"]*s_max

        TREE_max = 1.00000001 * np.ones((params["T"], 1))
        TREE_max = TREE_max.reshape(-1, 1)  # forme (T,1)

        with open(params["prev_run_folder"] + 'param_model.csv', newline='') as f:
            reader = csv.reader(f, delimiter=';')
            param_model = {row[0]: row[1] for row in reader}

        dimprevious = int(param_model["dim"])
        if dimprevious==3:
            u_min = np.genfromtxt(params["prev_run_folder"] + 'u_min.csv', dtype=float, delimiter=';')
            TREE_min = u_min[:, [2]]
        else:
            n_mc = 10000  # nombre de tirages Monte-Carlo
            u_mc = np.zeros((n_mc, params["T"]))  # stockage des trajectoires
            u_mc[:, 0] = 1

            for i in range(n_mc):
                for t_step in range(1, params["T"]):
                    u_1 = (t_min[t_step - 1]+t_max[t_step-1])/2  # colonne 1 de u_init (temp�rature fixe)
                    input_val = u_mc[i, t_step - 1]  # valeur de la for�t pr�c�dente
                    # Tirage al�atoire de shock2 pour chaque p�riode
                    shock2 = beta(params["beta_alpha"], params["beta_beta"]).rvs()

                    u_mc[i, t_step] = law_motion_state_space(t_step, u_1, shock2, input_val, params)

            TREE_min = 0.9*np.min(u_mc, axis=0)
            TREE_min = TREE_min.reshape(-1, 1)

        u_max = np.concatenate((kpec_max, t_max, TREE_max), axis=1)
        u_min = np.concatenate((kpec_min, t_min, TREE_min), axis=1)
    return u_min, u_max;

#law motion for state space 
def law_motion_state_space(t_step, u_1, shock2, input_val, params):
    temp_regional = np.log(((u_1 - params["T0"]) * params["rTCRE"] / params["TCRE_mean"])+1)
    base = (params["Upsilon"] * (1 - input_val)) / params["beta0"]
    base = np.maximum(base, 1e-6)
    a = params["growth0"] * (1 - base**params["eta"]) * (1 - input_val)
    sce = params["mean_scenarioEU"].iloc[t_step][0]
    b = temp_regional * (params["climate_maxEU"] * shock2)
    growth = a - b - sce
    if t_step < 38:
        TREEn = input_val * (1 + growth)
    else:
        TREEn = input_val
    TREEn = np.clip(TREEn, 0, 1)
    return TREEn

# definition of utility function
def utility(consumption, params):
    if params["theta"] == 1:
        c = np.maximum(consumption, 10 ** (-8) * np.ones(np.shape(consumption)))
        u = np.log(c)
    else:
        c = np.maximum(consumption, 10 ** (-8) * np.ones(np.shape(consumption)))
        u = c ** (1 - params["theta"]) / (1 - params["theta"])
    return u;

#marginal utility of consumption
def marginal_utility(consumption, params):
    c = np.maximum(consumption, 10 ** (-8) * np.ones(np.shape(consumption)))
    du = c ** (-params["theta"])
    return du;

# total production without damages per period
def production(TFP, capital, labor, params):
    y = params["deltaT"] * TFP * (labor ** (1 - params["alpha"])) * capital ** params["alpha"]
    return y;

# abatement cost in percentage GDP
def abatement_cost(ab_level, total_cost, params):
    ab = np.maximum(ab_level, 0)
    ab = np.minimum(ab, 1)
    b = total_cost * ab ** params["theta2"]
    return b;

#emissions (endogenous)
def emissions(ab_level, production, carbon_intensity):  # emissions per period
    ab = np.maximum(ab_level, 0)
    ab = np.minimum(ab, 1)
    e = carbon_intensity * (1 - ab) * production
    return e;

#damage factor
def damage_factor(temp, params):  # damages factor, correcting TFP, from temperature and shock
    omega = np.maximum(10 ** (-8), 1 - (params["pi1"] * temp + params["pi2"] * temp ** params["pi3"]))
    return omega;

#growth tipping
def growth_primary_tipping(t, coverage, temp, shock, params):  # growth rate over period between t-1 and t
    temp_regional = np.log(((temp - params["T0"]) * params["rTCRE"] / params["TCRE_mean"])+1)  # feedback effect
    base = (params["Upsilon"] * (1 - coverage)) / params["beta0"]
    base = np.maximum(base, 1e-6)   
    a = params["growth0"] * (1 - base**params["eta"]) * (1 - coverage)
    sce = params["mean_scenarioEU"].iloc[t][0]
    if params["stochastic"] == 2:
        b = temp_regional * (params["climate_maxEU"] * shock)  # because rtcre
    else:
        b = temp_regional * (params["climate_meanEU"])
    e = a - b - sce  
    return (e);

def density_betalaw(x, a, b, params):
    if b < a:
        print("pb order parameters function density_beta")

    if a <= x <= b:
        # Use scipy's beta.pdf with loc and scale for [a, b]
        v = beta.pdf(
            x,
            params["beta_alpha"],
            params["beta_beta"],
            loc=a,
            scale=b - a
        )
    else:
        v = 0.0
    return v

# function to maximise
def maximand(t, x, u, params):
    c = 1 - (params["bet"] * params["alpha"])  # c is relative consumption
    if params["dim"]==1:
        c = 1 - (params["bet"] * params["alpha"])  # c is relative consumption
        k = (params["A"][t, :] ** (1 / (1 - params["alpha"]))) * params["L"][t, :] * u
        conso = production(params["A"][t, :], k, params["L"][t, :], params) / (params["L"][t, :]) * c
        v = 1 / params["L"][[0], :] * params["L"][t, :] * utility(conso, params)
    elif params["dim"]!=1:
        x = np.maximum(x, np.zeros(np.shape(x)))
        x = np.minimum(x, np.ones(np.shape(x)))
        a = x  # a is relative abatement
        k = (params["A"][t, :] ** (1 / (1 - params["alpha"]))) * params["L"][t, :] * u[:, [0]]  # u[:,[0]] is capital per efficient capita
        conso = damage_factor(u[:, [1]], params) * (1 - abatement_cost(a, params["theta1"][t, :], params)) * production(
            params["A"][t, :], k, params["L"][t, :], params) / params["L"][t, :] * c
        v = 1 / params["L"][[0], :] * params["L"][t, :] * utility(conso, params)
    return v;

# partial derivative with respect to consumption of maximand
def partialc_maximand(t, x, u, params):
    c = 1 - (params["bet"] * params["alpha"])  # c is relative consumption
    if params["dim"]==1:
        k = (params["A"][t, :] ** (1 / (1 - params["alpha"]))) * params["L"][t, :] * u
        conso = production(params["A"][t, :], k, params["L"][t, :], params) / (params["L"][t, :]) * c
        v = 1 / params["L"][[0], :] * marginal_utility(conso, params)
    elif params["dim"]!=1:
        x = np.maximum(x, np.zeros(np.shape(x)))
        x = np.minimum(x, np.ones(np.shape(x)))
        c = 1 - (params["bet"] * params["alpha"])  # c is relative consumption
        a = x  # a is relative abatement
        k = (params["A"][t, :] ** (1 / (1 - params["alpha"]))) * params["L"][t, :] * u[:, [0]]  # u[:,[0]] is capital per efficient capita
        conso = damage_factor(u[:, [1]], params) * (1 - abatement_cost(a, params["theta1"][t, :], params)) * production(
            params["A"][t, :], k, params["L"][t, :], params) / params["L"][t, :] * c
        v = 1 / params["L"][[0], :] * marginal_utility(conso, params)
    return v;

#law motion of state variables
def law_motion_SEU(t, x, u, shock1, shock2, model, params):
    if params["dim"] == 1:
        x = 1 - (params["bet"] * params["alpha"])  # c is relative consumption
        k = (params["A"][[t - 1], :] ** (1 / (1 - params["alpha"]))) * params["L"][[t - 1], :] * u
        kn = (1 - params["delta"]) * k + production(params["A"][[t - 1], :], k, params["L"][[t - 1], :], params) * (1 - x)
        un = kn / ((params["A"][[t], :] ** (1 / (1 - params["alpha"]))) * params["L"][[t], :])
    elif params["dim"] == 3:
        x = np.maximum(x, np.zeros(np.shape(x)))
        x = np.minimum(x, np.ones(np.shape(x)))
        c = 1 - (params["bet"] * params["alpha"])  # c is relative consumption
        a = x  # a is relative abatement
        TREE = u[:, [2]]  # TFP
        k = (params["A"][[t - 1], :] ** (1 / (1 - params["alpha"]))) * params["L"][[t - 1], :] * u[:, [0]]  # u[:,[0]] is capital per efficient capita
        y = production(params["A"][[t - 1], :], k, params["L"][[t - 1], :], params) * damage_factor(u[:, [1]], params)
        kn = ((1 - params["delta"])) * k * damage_factor(u[:, [1]], params) ** (params["damage_capital"] / params["alpha"]) + y * ( 1 - abatement_cost(a, params["theta1"][[t - 1], :], params)) * (1 - c)

        sn = u[:, [1]].copy()
        sn += shock1 * (emissions(a, y, params["sigm"][[t - 1], :]))

        print("alors")
        print(params["nb_models"])
        if params["nb_models"] == 1:
            TREEn = TREE
        else:
            if t < 38:  # year2200
                TREEn = TREE * (1 + growth_primary_tipping(t, TREE, sn, shock2, params))
            if t >= 38:  # year2200
                TREEn = TREE
        TREEn = np.clip(TREEn, 1e-6, 1)
        sn += shock1 * params["carbonstock_amz"] * (TREE - TREEn)
        un = np.concatenate((kn / ((params["A"][[t], :] ** (1 / (1 - params["alpha"]))) * params["L"][[t], :]), sn, TREEn), axis=1)
    return un;

# terminal value function
def terminal_value(t, u, params):
    if not (t == params["T"] + 1):
        print('Are you sure terminal_value function is correctly called?')
    if params["dim"]==1:
        k = (params["A"][[t - 1], :] ** (1 / (1 - params["alpha"]))) * params["L"][[t - 1], :] * u
        kn = (params["A"][[t], :] ** (1 / (1 - params["alpha"]))) * params["L"][[t], :] * u
        # consumption for constant capital per efficient capita
        conso = (production(params["A"][[t - 1], :], k, params["L"][[t - 1], :], params) + (1 - params["delta"]) * k - kn) / (
        params["L"][[t - 1], :])
        # adjustement in terminal constraint adapted from Barr Manne (GA is growth rate per year, not per period)
        v = 1 / params["L"][[0], :] * params["L"][[t - 1], :] * 1 / (1 - (params["bet"]) * ((1 + params["GA"][[t - 1], :]) ** (
                        params["deltaT"] * (1 - params["theta"]) / (1 - params["alpha"])))) * utility(conso, params)
    elif params["dim"]!=1:
        kpec = u[:, [0]]
        k = (params["A"][[t - 1], :] ** (1 / (1 - params["alpha"]))) * params["L"][[t - 1], :] * kpec
        kn = (params["A"][[t], :] ** (1 / (1 - params["alpha"]))) * params["L"][[t], :] * kpec  # constant capital per efficient capita
        a = 1
        # a=0
        conso = (damage_factor(u[:, [1]], params) * (1 - abatement_cost(a, params["theta1"][[t - 1], :], params)) * production(
            params["A"][[t - 1], :], k, params["L"][[t - 1], :], params) + ((1 - params["delta"])) * k * damage_factor(
            u[:, [1]], params) ** (params["damage_capital"] / params["alpha"]) - kn) / params["L"][[t - 1], :]
        # adjustement in terminal constraint adapted from Barr Manne (GA is growth rate per year, not per period)
        v = 1 / params["L"][[0], :] * params["L"][[t - 1], :] * 1 / (1 - params["bet"] * (
                    (1 + params["GA"][[t - 1], :]) ** (
                        params["deltaT"] * (1 - params["theta"]) / (1 - params["alpha"])))) * utility(conso, params)
#    return 1.1*v;
    return v;

# number of Chebyshev polynomials of total degree <= deg in dimension dim
def n_terms_cheb(deg, dim):
    z = np.prod(np.arange(1 + deg, 1 + deg + dim)) / np.prod(np.arange(1, 1 + dim))
    return z

# indicator of a weight (number of positive partial degree)
def ind_cheb(weight):
    d = np.minimum(weight, np.ones(np.shape(weight)))
    d = np.sum(d, axis=1).reshape((d.shape[0], 1))
    return d

# matrix of weight of Chebyshev polynomials of total degree <= deg in dimension dim
def fill_cheb(deg, dim):
    if dim == 1:
        weight_cheb = np.arange(deg + 1).reshape((deg + 1, 1))
    elif deg == 0:
        weight_cheb = np.zeros((1, dim))
    else:
        weight_cheb1 = fill_cheb(deg, dim - 1)
        s1 = weight_cheb1.shape[0]
        weight_cheb1 = np.concatenate((np.zeros((s1, 1)), weight_cheb1), axis=1)
        weight_cheb2 = fill_cheb(deg - 1, dim)
        s2 = weight_cheb2.shape[0]
        weight_cheb2 = np.concatenate((np.ones((s2, 1)), np.zeros((s2, dim - 1))), axis=1) + weight_cheb2
        weight_cheb = np.concatenate((weight_cheb1, weight_cheb2), axis=0)
    return weight_cheb
def fill_cheb3(deg, dim):
    return fill_cheb(deg, dim)  # identical to fill_cheb
def fill_cheb2(a, dim, deg1=None, deg2=None, deg3=None):
    if dim == 1:
        a = a[(a[:, 0] / deg1 <= 1)]
    elif dim == 2:
        a = a[(a[:, 0] / deg1 + a[:, 1] / deg2 <= 1)]
    elif dim == 3:
        a = a[(a[:, 0] / deg1 + a[:, 1] / deg2 + a[:, 2] / deg3 <= 1)]
    return a

# Chebyshev polynomials in one dimension
def Chebyshev_one_dim(weight, x):
    if weight == 0:
        z = 1
    elif weight == 1:
        z = x
    elif weight == 2:
        z = 2 * x**2 - 1
    elif weight == 3:
        z = 4 * x**3 - 3 * x
    elif weight == 4:
        z = 8 * x**4 - 8 * x**2 + 1
    elif weight == 5:
        z = 16 * x**5 - 20 * x**3 + 5 * x
    else:
        z = 2 * x * Chebyshev_one_dim(weight - 1, x) - Chebyshev_one_dim(weight - 2, x)
    return z

# Chebyshev polynomials in multiple dimensions
def Chebyshev(weight, x):
    z = 0
    d = weight.shape[1]
    if d == x.shape[1]:
        z = 1
        for i in range(d):
            z = z * Chebyshev_one_dim(weight[:, [i]], x[:, [i]])
    else:
        print('Incompatible dimensions in Chebyshev')
    return z

# Approximation by Chebyshev polynomials
def approx_Chebyshev(coef, weight, x, x_min, x_max):
    z = 0
    nb_terms = coef.shape[1]
    if weight.shape[0] == nb_terms:
        for i in range(nb_terms):
            z += coef[:, [i]] * Chebyshev(weight[[i], :], (2*x - x_min - x_max) / (x_max - x_min))
    else:
        print('Check Chebyshev approximation: dimension mismatch')
    return z

# Chebyshev nodes
def nodes_Chebyshev2(d):
    d = int(d)
    ind = np.arange(1, d + 2).reshape((d + 1, 1))
    z = -np.cos((2 * ind - 1) * np.pi / (2 * (d + 1)))
    return z

def nodes_Chebyshev(params):
    a = nodes_Chebyshev2(params["deg1"])
    if params["dim"] == 2:
        b = nodes_Chebyshev2(params["deg2"])
        b = np.kron(np.ones((int(params["deg1"]+1), 1)), b)
        a = np.kron(a, np.ones((int(params["deg2"]+1), 1)))
        a = np.concatenate((a, b), axis=1)
    if params["dim"] == 3:
        c = nodes_Chebyshev2(params["deg3"])
        c = np.kron(np.ones((int((params["deg2"] + 1)*(params["deg1"]+1)), 1)), c)
        b = nodes_Chebyshev2(params["deg2"])
        b = np.kron(b, np.ones((int(params["deg3"]+1), 1)))
        b = np.tile(b, (int(params["deg1"]+1), 1))
        a = np.kron(a, np.ones((int((params["deg2"] + 1)*(params["deg3"]+1)), 1)))
        a = np.concatenate((a, b), axis=1)
        a = np.concatenate((a, c), axis=1)
    return a

#Adjust the upper and lower bounds of the Chebyshev approximation domain to improve
#    numerical stability near the boundary.
#   The adjustment expands the lower bound based on the spacing of Chebyshev
#    nodes, with the expansion capped at a fixed fraction of the original domain
#    length.

def adjust_min(x_min, x_max, m, max_frac=0.1):
    x_min = np.array(x_min, dtype=float)
    x_max = np.array(x_max, dtype=float)
    z_adjust = -np.cos(np.pi / (2 * m))
    x_adjust = (z_adjust + 1) / (-2 * z_adjust) * (x_max - x_min)
    max_padding = max_frac * np.maximum(x_max - x_min, 1e-12)
    x_adjust = np.minimum(x_adjust, max_padding)
    return x_min - x_adjust

def adjust_max(x_min, x_max, m, max_frac=0.1):
    x_min = np.array(x_min, dtype=float)
    x_max = np.array(x_max, dtype=float)
    z_adjust = -np.cos(np.pi / (2 * m))
    x_adjust = (z_adjust + 1) / (-2 * z_adjust) * (x_max - x_min)
    max_padding = max_frac * np.maximum(x_max - x_min, 1e-12)
    x_adjust = np.minimum(x_adjust, max_padding)
    return x_max + x_adjust

#bellman function interpolated
def Bell_max(x, *args):
    x = np.asarray(x)
    if x.shape != (1,):
        raise ValueError(f"Bell_max received x of shape {x.shape}")
    
    x_step=x.reshape(np.shape(np.ones((1,1))))
    t_step=args[0]
    u=args[1]
    params=args[2]
    coefn=args[3]

    if t_step>params["T"]:
        print('Something wrong: t_step is beyond time horizon')
    else:
        #the state variable
        u_step=u[[t_step-1],:]

    #terminal value function 
    if t_step==params["T"]:
        #stochastic shocks
        if params["stochastic"]==2:
            v1 = quadrature_product_terminal_EV(t_step, x_step, u_step, params, N1=5, N2=5)
        if params["stochastic"]==1:
            cint,err=si.quad(functools.partial(integrand_terminal1, t_step=t_step, x_step=x_step,u_step=u_step, params=params),params["TCRE_min"],params["TCRE_max"],epsabs=params["tolintabs"], epsrel=params["tolintrel"])
            v1=cint
        #deterministic shocks
        if params["stochastic"]==0:
            un=law_motion_SEU(t_step,x_step,u_step,params["TCRE_mean"],params["climate_meanEU"],0, params)
            #value with penalities for x_step outside bounds
            v1=maximand(np.array([t_step-1]),x_step,u_step, params)+params["bet"]*terminal_value(t_step+1,un, params)
        #expected value
        v=float(-v1)#minus because we will minimize this function (so as to maximize the Bellman function)

    elif t_step<params["T"]:
        #coefficient approximating the next step value function
        #stochastic shocks
        if params["stochastic"]==2:
            v1 = EV_quadrature_stochastic2(
                t_step=t_step,
                x_step=x_step,
                u_step=u_step,
                coefn=coefn,
                weight_cheb=params["weight_cheb"],
                u_min=params["u_min"][[t_step], :],
                u_max=params["u_max"][[t_step], :],
                params=params,
                N1=5, N2=5)
        if params["stochastic"]==1:
            cint,err=si.quad(functools.partial(integrand1, t_step=t_step, x_step=x_step,u_step=u_step,coefn=coefn, params=params),params["TCRE_min"],params["TCRE_max"],epsabs=params["tolintabs"], epsrel=params["tolintrel"])
            v1=cint
        if params["stochastic"]==0:
            un=law_motion_SEU(t_step,x_step,u_step,params["TCRE_mean"],params["climate_meanEU"],0, params)
            if params["dim"] == 1:
                un_min = params["u_min"][t_step]
                un_max = params["u_max"][t_step]
            else:
                un_min = params["u_min"][[t_step], :]
                un_max = params["u_max"][[t_step], :]
            un_min_adj = adjust_min(un_min, un_max, m=params["degb"], max_frac=0.1)
            un_max_adj = adjust_max(un_min, un_max, m=params["degb"], max_frac=0.1)
            un_adj = np.clip(un, un_min_adj, un_max_adj)
            size=params["degb"]+1
            v1=maximand(np.array([t_step-1]),x_step,u_step, params)+params["bet"]*approx_Chebyshev(coefn,params["weight_cheb"],un_adj,un_min_adj,un_max_adj)
        v=float(-v1)#minus because we will minimize this function (so as to maximize the Bellman function)
    return v;

#stochastic for dT/dS only
def integrand1(shock1, t_step, x_step, u_step, coefn, params):
    un = law_motion_SEU(t_step, x_step, u_step, shock1,params["climate_meanEU"],0, params)
    un_min = params["u_min"][[t_step], :]
    un_max = params["u_max"][[t_step], :]
    un_min_adj = adjust_min(un_min, un_max, m=params["degb"], max_frac=0.1)
    un_max_adj = adjust_max(un_min, un_max, m=params["degb"], max_frac=0.1)
    un_adj = np.clip(un, un_min_adj, un_max_adj)
    size=params["degb"]+1
    v1 = maximand(np.array([t_step - 1]), x_step, u_step, params) + params["bet"] * approx_Chebyshev(coefn, params["weight_cheb"], un_adj,un_min_adj,un_max_adj)
    v2=truncnorm.pdf(shock1, params["a_tcre"], params["b_tcre"], loc=params["TCRE_mean"], scale=params["sdeviation"])
    v = v1 * v2
    return v;

#Evaluate the continuation value at time t+1 for given realizations of
#stochastic climate shocks.
#for stochastic=2
def EV_simulation_stochastic2(
    t_step, x_step, u_step,
    coefn, weight_cheb,
    u_min, u_max,
    params,
    shock1, shock2):
    un = law_motion_SEU(
        t_step, x_step, u_step,
        shock1, shock2, 1, params)
    u_min_adj = adjust_min(u_min, u_max, m=params["degb"], max_frac=0.1)
    u_max_adj = adjust_max(u_min, u_max, m=params["degb"], max_frac=0.1)
    un_adj = np.clip(un, u_min_adj, u_max_adj)
    V_next = approx_Chebyshev(
        coefn, weight_cheb,
        un_adj, u_min_adj, u_max_adj)
    Vval = (
        maximand(np.array([t_step - 1]), x_step, u_step, params)
        + params["bet"] * V_next)
    return Vval

#for stochastic=1
def integrand1_simulation(
    shock1, t_step, x_step, u_step,
    coefn, weight_cheb, u_min, u_max, params):
    un = law_motion_SEU(
        t_step, x_step, u_step,
        shock1, None, 1, params)
    u_min_adj = adjust_min(u_min, u_max, m=params["degb"], max_frac=0.1)
    u_max_adj = adjust_max(u_min, u_max, m=params["degb"], max_frac=0.1)
    un_adj = np.clip(un, u_min_adj, u_max_adj)
    V_next = approx_Chebyshev(
        coefn, weight_cheb,
        un_adj, u_min_adj, u_max_adj)
    Vval = (
        maximand(np.array([t_step - 1]), x_step, u_step, params)
        + params["bet"] * V_next)
    pdf = truncnorm.pdf(
        shock1,
        params["a_tcre"], params["b_tcre"],
        loc=params["TCRE_mean"],
        scale=params["sdeviation"])
    return Vval * pdf

#terminal value function 
#for stochastic 1
def integrand_terminal1(shock1, t_step, x_step, u_step, params):
    # next-stage state variable should be in the domain of the approximation of next-stage value function
    un = law_motion_SEU(t_step, x_step, u_step, shock1,params["climate_meanEU"],0, params)
    v1 = maximand(np.array([t_step - 1]), x_step, u_step, params) + params["bet"] * terminal_value(t_step + 1,un, params)
    v2=truncnorm.pdf(shock1, params["a_tcre"], params["b_tcre"], loc=params["TCRE_mean"], scale=params["sdeviation"])
    v = v1 * v2
    return v;
    return v;