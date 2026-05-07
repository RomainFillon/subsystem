# -*- coding: utf-8 -*-

##############################################################################
#replication for "The need for regulation of climate subsystems" [simulation]
###############################################################################
###############################################################################
import sys

#sys.path.append(r"")
#sys.path.append(r"")

#Define runs to be simulated
#run_tobesimulated= 'final_amazon_tcre_run0004/'
#run_tobesimulated= 'final_amazon_tcre_run0005/'
#run_tobesimulated= 'final_amazon_tcre_run0006/'
#run_tobesimulated= 'final_amazon_tcre_run0009/'
#run_tobesimulated= 'final_amazon_tcre_run0010/'
#run_tobesimulated= 'final_amazon_tcre_run0011/'
#run_tobesimulated= 'final_amazon_tcre_run0014/'
#run_tobesimulated= 'final_amazon_tcre_run0015/'
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
base_folder = r''
preprod_folder   = os.path.join(base_folder, 'run') + pathsep
model_folder     = os.path.join(base_folder, 'model') + pathsep
param_folder     = os.path.join(base_folder, 'parameters') + pathsep
outputs_folder   = os.path.join(base_folder, 'outputs') + pathsep
run_folder = outputs_folder + run_tobesimulated

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
    final_doc= np.zeros((params["T"]-2, 14))
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

    # compute SCC & SCD
    # deterministic vs stochastic (integrand + dT/dS for SCC)
    if params["stochastic"] == 0:
        # Temperature channel
        # Channel standard
        # dUt+1/dTt+1
        SCC_1 = (approx_Chebyshev(coef_V[1:params["T"], :], params["weight_cheb"], uplus[1:params["T"], :], u_min[1:T, :], u_max[1:params["T"], :]) - approx_Chebyshev(coef_V[1:params["T"], :],
                       params["weight_cheb"], uminus[1:params["T"], :], u_min[1:params["T"], :],u_max[1:params["T"],:])) / (2 * dT[1:params["T"], :])
        # dUt+1/dTt+1 . dTt/dSt
        SCC_2 = SCC_1 * params["TCRE_mean"]
        # Channel scaling
        # dT/dA
        SCC_3 = - params["TCRE_mean"] * np.ones((params["T"] - 1, 1)) * params["carbonstock_amz"]
        # dAt+1/dTt
        SCC_4 = np.zeros((params["T"] - 1, 1))
        # dAt+1/dTt
        for t in range(1, params["T"] - 1):
            SCC_4[t - 1] = (growth_primary_tipping(t, u[[t - 1], [2]], (1 + param_delta) * u[[t - 1], [1]],params["climate_meanEU"], params) - growth_primary_tipping(t, u[
                [t - 1], [2]], (1 - param_delta) * u[[t - 1], [1]], params["climate_meanEU"], params)) / (2 * param_delta * u[[t - 1], [1]])
        # dAt+1/dAt - assumed to be one because of short-term oscillation
        SCD_2 = np.ones((params["T"] - 1, 1))
        # log transformation of dAt+1/dAt against short-term oscillations
        SCD_2 = np.zeros((params["T"] - 1, 1))
        eps = 1e-12  # petit epsilon pour éviter log(0)
        for t in range(1, params["T"] - 1):
            num = (growth_primary_tipping(t, (1 + param_delta) * u[[t - 1], [2]], u[[t - 1], [1]],
                                          params["climate_meanEU"], params) - growth_primary_tipping(t, (1 - param_delta) *u[[t - 1], [2]],u[[t - 1], [1]],params["climate_meanEU"], params))
            den = 2 * param_delta * u[[t - 1], [2]]
            deriv = num / den  
            sign = np.sign(deriv) if np.abs(deriv) > eps else 1.0
            mag = np.log(np.abs(deriv) + eps)  # log-magnitude
            SCD_2[t - 1] = sign * np.exp(mag)

        # now SCD
        # SCD temperature channel
        standard_T = np.zeros((params["T"] - 1, 1))
        for t in range(params["T"] - 1):
            sum_k = 0.0
            max_k = len(marginal_utility_temp) - t - 1
            for k in range(0, max_k):
                sum_k += params["bet"] ** k * marginal_utility_temp[t + k + 1] * SCC_3[t + 1]
            standard_T[t] = sum_k

        # SCD subsystem channel
        standard_A = np.zeros((params["T"] - 1, 1))
        for t in range(params["T"] - 1):
            sum_k = 0.0
            max_k = len(marginal_utility_temp) - t - 1
            for k in range(1, max_k):
                sum_j = 0.0
                for j in range(1, k):
                    if j > 1:
                        prod_m = np.prod(SCD_2[t + 1:t + j])
                    else:
                        prod_m = 1.0
                    sum_j += SCC_3[t + j] * prod_m
                sum_j = sum_j * marginal_utility_temp[t + k + 1]
                sum_k += (params["bet"] ** k) * sum_j
            standard_A[t] = sum_k

        #this is SCD and SCD decomposition
        final_doc[:, [4]] = (1 / params["carbonstock_amz"]) * monetary[0:params["T"] - 2] * (standard_A[0:params["T"] - 2] + standard_T[0:params["T"] - 2])
        final_doc[:, [5]] = (standard_T[0:params["T"] - 2])
        final_doc[:, [6]] = (standard_A[0:params["T"] - 2])

        # cross propagation 
        cross_A = np.zeros(params["T"] - 1)
        
        for t in range(params["T"] - 1):
            max_k = len(marginal_utility_temp) - t - 1
            sum_k = 0.0
            for k in range(1, max_k):
                sum_j = 0.0
                for j in range(1, k):
                    prod_AA = np.prod(SCD_2[t + 1: t + j]) if j > 1 else 1.0
                    term = marginal_utility_temp[t + k + 1] * SCC_3[t + j] * prod_AA
                    sum_j += term
                sum_k += sum_j
            cross_A[t] = sum_k
        cross_A_col = cross_A.reshape(-1, 1)  # (99,1)
        SCC_5 = standard_A + cross_A_col

        #SCCDS with Amazon rainforest
        if params["id_climate"] == 1:
            SCCDS = - monetary * (SCC_2 * (1 + SCC_3 * SCC_4) + SCC_5 * SCC_4 * params["TCRE_mean"])
        #without
        else:
            SCCDS = - monetary * (SCC_2)
        final_doc[:, [7]] = SCCDS[0:params["T"] - 2]

        # now complete decomposition of SCCDS (for all periods)
        # standard propagation via T (and dT/dT=1)
        standard_T = np.zeros((params["T"] - 1, 1))
        for t in range(params["T"] - 1):
            sum_k = 0.0
            max_k = len(marginal_utility_temp) - t - 1
            for k in range(0, max_k):
                sum_k += params["bet"] ** k * marginal_utility_temp[t + k + 1] * params["TCRE_mean"]
            standard_T[t] = sum_k
        final_doc[:, [8]] = standard_T[0:params["T"] - 2]

        # standard propagation via A (and dT/dT=1)
        standard_A = np.zeros((params["T"] - 1, 1))
        for t in range(params["T"] - 1):
            sum_k = 0.0
            max_k = len(marginal_utility_temp) - t - 1
            for k in range(1, max_k):
                sum_j = 0.0
                for j in range(1, k):
                    if j > 1:
                        prod_m = np.prod(SCD_2[t + 1:t + j])
                    else:
                        prod_m = 1.0
                    sum_j += SCC_3[t + j] * prod_m
                sum_j = sum_j * marginal_utility_temp[t + k + 1] * SCC_4[t] * params["TCRE_mean"]
                sum_k += (params["bet"] ** k) * sum_j
            standard_A[t] = sum_k
        final_doc[:, [9]] = standard_A[0:params["T"] - 2]

        # expected cross via A (III_T)
        # Dimensions: T x T x T (t, k, j)
        L1T_3D = np.zeros((params["T"] - 1, params["T"] - 1, params["T"] - 1))
        L2T_3D = np.zeros((params["T"] - 1, params["T"] - 1, params["T"] - 1))

        for t in range(params["T"] - 1):
            max_k = len(marginal_utility_temp) - t - 1
            for k in range(1, max_k):
                for j in range(0, k):
                    L1T_3D[t, k, j] = marginal_utility_temp[t + k + 1] * SCC_3[t + k + 1]
                    if (k - 1) >= (j + 1):
                        prod_AA = np.prod(SCD_2[t + j + 1:t + k])
                    else:
                        prod_AA = 1.0
                    L2T_3D[t, k, j] = prod_AA * SCC_4[t + j + 1] * params["TCRE_mean"]

        # express crosss via T (III_A)
        L1A_3D = np.zeros((params["T"] - 1, params["T"] - 1, params["T"] - 1))
        L2A_3D = np.zeros((params["T"] - 1, params["T"] - 1, params["T"] - 1))

        for t in range(params["T"] - 1):
            max_k = len(marginal_utility_temp) - t - 1
            for k in range(1, max_k):
                for j in range(1, k):
                    L1A_3D[t, k, j] = marginal_utility_temp[t + k + 1] * SCC_3[t + j]
                    if j >= 1:
                        prod_AA = np.prod(SCD_2[t + 1:t + j])
                    else:
                        prod_AA = 1.0
                    L2A_3D[t, k, j] = prod_AA * SCC_4[t + 1] * params["TCRE_mean"]


    #now with dT/dS stochastic
    if params["stochastic"] == 1:
        SCC_1 = np.zeros((params["T"] - 1, 1))
        # x_t+1 and u_t+1 are already computed from stochastic draw
        for t in range(0, params["T"] - 1):
            cint, err = si.quad(
                functools.partial(
                    integrand1_simulation,
                    t_step=t + 1,
                    x_step=x[[t + 1], :],
                    u_step=uplus[[t + 1], :],
                    coefn=coef_V[[t + 1], :],
                    weight_cheb=params["weight_cheb"],
                    u_min=params["u_min"][[t + 1], :],
                    u_max=params["u_max"][[t + 1], :],
                    params=params
                ), params["TCRE_min"], params["TCRE_max"],epsabs=params["tolintabs"], epsrel=params["tolintrel"])

            cint2, err = si.quad(
                functools.partial(
                    integrand1_simulation,
                    t_step=t + 1,
                    x_step=x[[t + 1], :],
                    u_step=uminus[[t + 1], :],
                    coefn=coef_V[[t + 1], :],
                    weight_cheb=params["weight_cheb"],
                    u_min=params["u_min"][[t + 1], :],
                    u_max=params["u_max"][[t + 1], :],
                    params=params), params["TCRE_min"], params["TCRE_max"], epsabs=params["tolintabs"], epsrel=params["tolintrel"])
            SCC_1[t] = (cint - cint2) / (2 * dT[[t + 1], :])

        # dU/dT . dTt/dSt
        SCC_2 = np.zeros((params["T"] - 1, 1))

        for t in range(0, params["T"] - 2):
            SCC_2[t] = SCC_1[t] * shock_save1[t]
        # dS/dA (this shock is t+1) to have dT/dA
        SCC_3 = - np.ones((params["T"] - 1, 1)) * params["carbonstock_amz"]
        for t in range(0, params["T"] - 2):
            SCC_3[t] = SCC_3[t] * shock_save1[t]
        # dAt/dTt
        SCC_4 = np.zeros((params["T"] - 1, 1))
        for t in range(1, params["T"] - 1):
            SCC_4[t - 1] = (growth_primary_tipping(t, u[[t], [2]], (1 + param_delta) * u[[t], [1]], params["climate_meanEU"], params) - growth_primary_tipping(t, u[
                [t], [2]], (1 - param_delta) * u[[t], [1]], params["climate_meanEU"], params)) / (2 * param_delta * u[[t], [1]])
        # dA/dS #shock is t (vs t+1 above)
        SCC_4b = np.zeros((params["T"] - 1, 1))
        for t in range(0, params["T"] - 2):
            SCC_4b[t] = SCC_4[t] * shock_save1[t]

        SCD_2 = np.zeros((params["T"] - 1, 1))
        phi_vals = []
        for t in range(1, params["T"] - 1):
            num = (growth_primary_tipping(t,(1 + param_delta) * u[[t - 1], [2]], u[[t - 1], [1]], params["climate_meanEU"], params)
                    - growth_primary_tipping(t,(1 - param_delta) * u[[t - 1], [2]], u[[t - 1], [1]], params["climate_meanEU"],params))
            den = 2 * param_delta * u[[t - 1], [2]]
            deriv = num / den
            phi_vals.append(float(deriv))

        # expected transition operator
        phi_A = np.mean(phi_vals)
        # impose mild biological stability if needed
        phi_A = np.clip(phi_A, 0, 1)
        # fill SCD_2 with constant operator
        SCD_2[:] = phi_A
        # now SCD
        # standard propagation via T (and dT/dT=1) SCD1
        standard_T = np.zeros((params["T"] - 1, 1))
        for t in range(params["T"] - 1):
            sum_k = 0.0
            max_k = len(marginal_utility_temp) - t - 1
            for k in range(0, max_k):
                sum_k += params["bet"] ** k * marginal_utility_temp[t + k + 1] * SCC_3[t + 1]
            standard_T[t] = sum_k

        # updated propagation via A SCD2
        standard_A = np.zeros((params["T"] - 1, 1))
        for t in range(params["T"] - 1):
            sum_k = 0.0
            max_k = len(marginal_utility_temp) - t - 1
            for k in range(1, max_k):
                sum_j = 0.0
                for j in range(1, k):
                    if j > 1:
                        prod_m = np.prod(SCD_2[t + 1:t + j])
                    else:
                        prod_m = 1.0
                    sum_j += SCC_3[t + j] * prod_m
                #                            sum_j += SCC_3[t + 1 + j] * prod_m
                sum_j = sum_j * marginal_utility_temp[t + k + 1]
                sum_k += (params["bet"] ** k) * sum_j
            standard_A[t] = sum_k

        #save SCD and decomposition
        final_doc[:, [4]] = (1 / params["carbonstock_amz"]) * monetary[0:params["T"] - 2] * (standard_A[0:params["T"] - 2] + standard_T[0:params["T"] - 2])
        final_doc[:, [5]] = (standard_T[0:params["T"] - 2])
        final_doc[:, [6]] = (standard_A[0:params["T"] - 2])

        # cross propagation SCD3
        cross_A = np.zeros(params["T"] - 1)
        
        for t in range(params["T"] - 1):
            max_k = len(marginal_utility_temp) - t - 1
            sum_k = 0.0
            for k in range(1, max_k):
                sum_j = 0.0
                for j in range(1, k):
                    prod_AA = np.prod(SCD_2[t + 1: t + j]) if j > 1 else 1.0
                    term = marginal_utility_temp[t + k + 1] * SCC_3[t + j] * prod_AA
                    sum_j += term
                sum_k += (params["bet"] ** k) * sum_j
            cross_A[t] = sum_k
        
        cross_A_col = cross_A.reshape(-1, 1)  # (99,1)
        SCC_5 = standard_A + cross_A_col

        shock_save1 = np.array(shock_save1[0:params["T"] - 1]).reshape(-1, 1)  # shape (99,1)

        # SCCDS with/without amazon rainforest
        if params["id_climate"] == 1:
            SCCDS = - monetary * (SCC_2 * (1 + SCC_3 * SCC_4) + SCC_5 * SCC_4 * shock_save1)
        else:
            SCCDS = - monetary * (SCC_2)

        final_doc[:, [7]] = SCCDS[0:params["T"] - 2]

        # now complete decomposition of SCCDS (for all periods)
        # standard propagation via T (and dT/dT=1)
        standard_T = np.zeros((params["T"] - 1, 1))
        for t in range(params["T"] - 1):
            sum_k = 0.0
            max_k = len(marginal_utility_temp) - t - 1
            for k in range(0, max_k):
                sum_k += params["bet"] ** k * marginal_utility_temp[t + k + 1] * shock_save1[t]
            standard_T[t] = sum_k
        final_doc[:, [8]] = standard_T[0:params["T"] - 2]

        # standard propagation via A (and dT/dT=1)
        standard_A = np.zeros((params["T"] - 1, 1))
        for t in range(params["T"] - 1):
            sum_k = 0.0
            max_k = len(marginal_utility_temp) - t - 1
            for k in range(1, max_k):
                sum_j = 0.0
                for j in range(1, k):
                    if j > 1:
                        prod_m = np.prod(SCD_2[t + 1:t + j])
                    else:
                        prod_m = 1.0
                    sum_j += SCC_3[t + j] * prod_m
                sum_j = sum_j * marginal_utility_temp[t + k + 1] * SCC_4[t] * shock_save1[t]
                sum_k += (params["bet"] ** k) * sum_j
            standard_A[t] = sum_k
        final_doc[:, [9]] = standard_A[0:params["T"] - 2]

        # expected cross via A (III_T)
        # Dimensions: T x T x T (t, k, j)
        L1T_3D = np.zeros((params["T"] - 1, params["T"] - 1, params["T"] - 1))
        L2T_3D = np.zeros((params["T"] - 1, params["T"] - 1, params["T"] - 1))

        for t in range(params["T"] - 1):
            max_k = len(marginal_utility_temp) - t - 1
            for k in range(1, max_k):
                for j in range(0, k):
                    L1T_3D[t, k, j] = marginal_utility_temp[t + k + 1] * SCC_3[t + k + 1]
                    if (k - 1) >= (j + 1):
                        prod_AA = np.prod(SCD_2[t + j + 1:t + k])
                    else:
                        prod_AA = 1.0
                    L2T_3D[t, k, j] = prod_AA * SCC_4[t + j + 1] * shock_save1[t]

        # express crosss via T (III_A)
        L1A_3D = np.zeros((params["T"] - 1, params["T"] - 1, params["T"] - 1))
        L2A_3D = np.zeros((params["T"] - 1, params["T"] - 1, params["T"] - 1))

        for t in range(params["T"] - 1):
            max_k = len(marginal_utility_temp) - t - 1
            for k in range(1, max_k):
                for j in range(1, k):
                    L1A_3D[t, k, j] = marginal_utility_temp[t + k + 1] * SCC_3[t + j]
                    if j >= 1:
                        prod_AA = np.prod(SCD_2[t + 1:t + j])
                    else:
                        prod_AA = 1.0
                    L2A_3D[t, k, j] = prod_AA * SCC_4[t + 1] * shock_save1[t]

    #now both dT/dS and dA/dT stochastic
    if params["stochastic"] == 2:
        SCC_1 = np.zeros((params["T"] - 1, 1))
        for t in range(0, params["T"] - 1):
            cint = EV_simulation_stochastic2(
        t_step=t + 1,
        x_step=x[[t + 1], :],
        u_step=uplus[[t + 1], :],
        coefn=coef_V[[t + 1], :],
        weight_cheb=params["weight_cheb"],
        u_min=params["u_min"][[t + 1], :],
        u_max=params["u_max"][[t + 1], :],
        params=params,
        shock1=shock_save1[t + 1],
        shock2=shock_save2[t + 1])
            cint2 = EV_simulation_stochastic2(
        t_step=t + 1,
        x_step=x[[t + 1], :],
        u_step=uminus[[t + 1], :],
        coefn=coef_V[[t + 1], :],
        weight_cheb=params["weight_cheb"],
        u_min=params["u_min"][[t + 1], :],
        u_max=params["u_max"][[t + 1], :],
        params=params,
        shock1=shock_save1[t + 1],
        shock2=shock_save2[t + 1])
            SCC_1[t] = (cint - cint2) / (2 * dT[[t + 1], :])

        # dU/dT . dT/dS (and period t vs t+1)
        SCC_2 = np.zeros((params["T"] - 1, 1))
        for t in range(0, params["T"] - 2):
            SCC_2[t] = SCC_1[t] * shock_save1[t]
        # dS/dA (same period shock)
        SCC_3 = - np.ones((params["T"] - 1, 1)) * params["carbonstock_amz"]
        for t in range(0, params["T"] - 2):
            SCC_3[t] = SCC_3[t] * shock_save1[t]
        # dA/dT
        SCC_4 = np.zeros((params["T"] - 1, 1))
        for t in range(1, params["T"] - 1):
            SCC_4[t - 1] = (growth_primary_tipping(t, u[[t], [2]], (1 + param_delta) * u[[t], [1]], shock_save2[t],
                                                   params) - growth_primary_tipping(t, u[[t], [2]],(1 - param_delta) * u[[t], [1]],shock_save2[t], params)) / (2 * param_delta * u[[t], [1]])
        # then multiplied by dT/dS
        # but past period
        SCC_4b = np.zeros((params["T"] - 1, 1))
        for t in range(0, params["T"] - 2):
            SCC_4b[t] = SCC_4[t] * shock_save1[t]

        SCD_2 = np.zeros((params["T"] - 1, 1))
        phi_vals = []
        for t in range(1, params["T"] - 1):
            num = (growth_primary_tipping(t, (1 + param_delta) * u[[t - 1], [2]], u[[t - 1], [1]], params["climate_meanEU"],params)
                    - growth_primary_tipping(t, (1 - param_delta) * u[[t - 1], [2]], u[[t - 1], [1]], params["climate_meanEU"],params))
            den = 2 * param_delta * u[[t - 1], [2]]
            deriv = num / den
            phi_vals.append(float(deriv))

        # expected transition operator
        phi_A = np.mean(phi_vals)
        # impose mild biological stability if needed
        phi_A = np.clip(phi_A, 0, 1)
        # fill SCD_2 with constant operator
        SCD_2[:] = phi_A

        # now SCD
        # SCD temperature
        standard_T = np.zeros((params["T"] - 1, 1))
        for t in range(params["T"] - 1):
            sum_k = 0.0
            max_k = len(marginal_utility_temp) - t - 1
            for k in range(0, max_k):
                sum_k += params["bet"] ** k * marginal_utility_temp[t + k + 1] * SCC_3[t + 1]
            standard_T[t] = sum_k

        # SCD subsystem
        standard_A = np.zeros((params["T"] - 1, 1))
        for t in range(params["T"] - 1):
            sum_k = 0.0
            max_k = len(marginal_utility_temp) - t - 1
            for k in range(1, max_k):
                sum_j = 0.0
                for j in range(1, k):
                    if j > 1:
                        prod_m = np.prod(SCD_2[t + 1:t + j])
                    else:
                        prod_m = 1.0
                    sum_j += SCC_3[t + j] * prod_m
                sum_j = sum_j * marginal_utility_temp[t + k + 1]
                sum_k += (params["bet"] ** k) * sum_j
            standard_A[t] = sum_k

        #save SCD and decomposition
        final_doc[:, [4]] = (1 / params["carbonstock_amz"]) * monetary[0:params["T"] - 2] * (standard_A[0:params["T"] - 2] + standard_T[0:params["T"] - 2])
        final_doc[:, [5]] = (standard_T[0:params["T"] - 2])
        final_doc[:, [6]] = (standard_A[0:params["T"] - 2])

        # cross propagation 
        cross_A = np.zeros(params["T"] - 1)
  
        for t in range(params["T"] - 1):
            max_k = len(marginal_utility_temp) - t - 1
            sum_k = 0.0
            for k in range(1, max_k):
                sum_j = 0.0
                for j in range(1, k):
                    prod_AA = np.prod(SCD_2[t + 1: t + j]) if j > 1 else 1.0
                    term = marginal_utility_temp[t + k + 1] * SCC_3[t + j] * prod_AA
                    sum_j += term
                sum_k += (params["bet"] ** k) * sum_j
            cross_A[t] = sum_k
        
        cross_A_col = cross_A.reshape(-1, 1)  # (99,1)
        SCC_5 = standard_A + cross_A_col

        shock_save1 = np.array(shock_save1[0:params["T"] - 1]).reshape(-1, 1)  # shape (99,1)

        # SCCDS with/without amazon rainforest
        if params["id_climate"] == 1:
            SCCDS = - monetary * (SCC_2 * (1 + SCC_3 * SCC_4) + SCC_5 * SCC_4 * shock_save1)
        else:
            SCCDS = - monetary * (SCC_2)
        final_doc[:, [7]] = SCCDS[0:params["T"] - 2]

        # now complete decomposition for SCCDS (for all periods)
        # standard propagation via T (and dT/dT=1)
        standard_T = np.zeros((params["T"] - 1, 1))
        for t in range(params["T"] - 1):
            sum_k = 0.0
            max_k = len(marginal_utility_temp) - t - 1
            for k in range(0, max_k):
                sum_k += params["bet"] ** k * marginal_utility_temp[t + k + 1] * shock_save1[t]
            standard_T[t] = sum_k
        final_doc[:, [8]] = standard_T[0:params["T"] - 2]

        # standard propagation via A (and dT/dT=1)
        standard_A = np.zeros((params["T"] - 1, 1))
        for t in range(params["T"] - 1):
            sum_k = 0.0
            max_k = len(marginal_utility_temp) - t - 1
            for k in range(1, max_k):
                sum_j = 0.0
                for j in range(1, k):
                    if j > 1:
                        prod_m = np.prod(SCD_2[t + 1:t + j])
                    else:
                        prod_m = 1.0
                    sum_j += SCC_3[t + j] * prod_m
                sum_j = sum_j * marginal_utility_temp[t + k + 1] * SCC_4[t] * shock_save1[t]
                sum_k += (params["bet"] ** k) * sum_j
            standard_A[t] = sum_k
        final_doc[:, [9]] = standard_A[0:params["T"] - 2]

        # expected cross via A (III_T)
        # Dimensions: T x T x T (t, k, j)
        L1T_3D = np.zeros((params["T"] - 1, params["T"] - 1, params["T"] - 1))
        L2T_3D = np.zeros((params["T"] - 1, params["T"] - 1, params["T"] - 1))

        for t in range(params["T"] - 1):
            max_k = len(marginal_utility_temp) - t - 1
            for k in range(1, max_k):
                for j in range(0, k):
                    L1T_3D[t, k, j] = marginal_utility_temp[t + k + 1] * SCC_3[t + k + 1]
                    if (k - 1) >= (j + 1):
                        prod_AA = np.prod(SCD_2[t + j + 1:t + k])
                    else:
                        prod_AA = 1.0
                    L2T_3D[t, k, j] = prod_AA * SCC_4[t + j + 1] * shock_save1[t]

        # express crosss via T (III_A)
        L1A_3D = np.zeros((params["T"] - 1, params["T"] - 1, params["T"] - 1))
        L2A_3D = np.zeros((params["T"] - 1, params["T"] - 1, params["T"] - 1))

        for t in range(params["T"] - 1):
            max_k = len(marginal_utility_temp) - t - 1
            for k in range(1, max_k):
                for j in range(1, k):
                    L1A_3D[t, k, j] = marginal_utility_temp[t + k + 1] * SCC_3[t + j]
                    if j >= 1:
                        prod_AA = np.prod(SCD_2[t + 1:t + j])
                    else:
                        prod_AA = 1.0
                    L2A_3D[t, k, j] = prod_AA * SCC_4[t + 1] * shock_save1[t]

    return final_doc, L1T_3D, L2T_3D, L1A_3D, L2A_3D;

if __name__ == '__main__':
    BATCH_SIZE = params["CPU"]  
    temp_folder = os.path.join(run_folder, "temp_results")
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

    master_seed = 12345
    seeds = np.random.SeedSequence(master_seed).spawn(params["draws"])

    batch = []
    batch_indices = []

   #parallel computing on different simulation runs
   #save all runs to compute summary statistics, distributions, and covariances
   #handle memory with batch and for covariance computation
    with multiprocess.Pool(processes=min(128, params["draws"])) as pool:
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

    pool.close()
    pool.join()

    batch_files = sorted(os.path.join(temp_folder, f) for f in os.listdir(temp_folder) if f.endswith('.npz'))

    MAX_RESERVOIR = 500   # (RAM-safe)
    batch_files = sorted(
        os.path.join(temp_folder, f)
        for f in os.listdir(temp_folder)
        if f.endswith(".npz"))
    
    mean = None
    min_val = None
    max_val = None
    reservoir = []
    n = 0
    #covariance 
    mean_L1T = mean_L2T = M2_T = None
    mean_L1A = mean_L2A = M2_A = None
    n_cov = 0
    for bf in batch_files:
        with np.load(bf, allow_pickle=True) as data:
            for key in data.files:
                res = data[key].tolist()
                X = np.asarray(res[0], dtype=float)  # (T-2, K)
                if np.any(X[0:20, 0] <= 0.05):
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
                    mean_L1T = np.zeros((Tm1, K, K))
                    mean_L2T = np.zeros((Tm1, K, K))
                    M2_T = np.zeros((Tm1, K, K))
                    mean_L1A = np.zeros((Tm1, K, K))
                    mean_L2A = np.zeros((Tm1, K, K))
                    M2_A = np.zeros((Tm1, K, K))
                n_cov += 1
                #temperature
                d1 = L1T - mean_L1T
                mean_L1T += d1 / n_cov
                d2 = L2T - mean_L2T
                mean_L2T += d2 / n_cov
                M2_T += d1 * (L2T - mean_L2T)
    
                #subsystem
                d1 = L1A - mean_L1A
                mean_L1A += d1 / n_cov
                d2 = L2A - mean_L2A
                mean_L2A += d2 / n_cov
                M2_A += d1 * (L2A - mean_L2A)    
    reservoir = np.stack(reservoir, axis=0)
    
    final_mean = mean
    final_min = min_val
    final_max = max_val
    final_05 = np.percentile(reservoir, 5, axis=0)
    final_95 = np.percentile(reservoir, 95, axis=0)
    
    beta = params["bet"]
    cov_T = M2_T / (n_cov - 1)
    cov_A = M2_A / (n_cov - 1)

    IV_T = np.zeros(Tm1)
    IV_A = np.zeros(Tm1)
    for t in range(Tm1):
        for k in range(1, K):
            for j in range(0, k):
                IV_T[t] += (beta ** k) * cov_T[t, k, j]
            for j in range(1, k + 1):
                IV_A[t] += (beta ** k) * cov_A[t, k, j]

    final_mean[:, [12]] = IV_T[:params["T"] - 2].reshape(-1, 1)
    final_mean[:, [13]] = IV_A[:params["T"] - 2].reshape(-1, 1)
    
    III_T = np.zeros(Tm1)
    III_A = np.zeros(Tm1)
    
    for t in range(Tm1):
        for k in range(1, K):
            for j in range(0, k):
                III_T[t] += (beta ** k) * mean_L1T[t, k, j] * mean_L2T[t, k, j]
            for j in range(1, k + 1):
                III_A[t] += (beta ** k) * mean_L1A[t, k, j] * mean_L2A[t, k, j]

    final_mean[:, [10]] = III_T[:params["T"] - 2].reshape(-1, 1)
    final_mean[:, [11]] = III_A[:params["T"] - 2].reshape(-1, 1)
    
    final95 = np.column_stack((final_95[:, 0], final_95[:, 2], final_95[:, 3]))
    final5  = np.column_stack((final_05[:, 0], final_05[:, 2], final_05[:, 3]))
    finalmin = np.column_stack((final_min[:, 0], final_min[:, 2], final_min[:, 3]))
    finalmax = np.column_stack((final_max[:, 0], final_max[:, 2], final_max[:, 3]))
    
    np.savetxt(run_folder + 'outputs_stochastic95.csv', final95, delimiter=';')
    np.savetxt(run_folder + 'outputs_stochastic5.csv', final5, delimiter=';')
    np.savetxt(run_folder + 'outputs_stochastic.csv', final_mean, delimiter=';')
    np.savetxt(run_folder + 'outputs_stochastic_max.csv', finalmax, delimiter=';')
    np.savetxt(run_folder + 'outputs_stochastic_min.csv', finalmin, delimiter=';')

    print("valid_runs")
    print(valid_runs)
