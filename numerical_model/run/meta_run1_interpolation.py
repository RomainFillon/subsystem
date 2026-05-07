# -*- coding: utf-8 -*-
################################################################################
#replication for "The need for regulation of climate subsystems" [interpolation]
################################################################################
import sys
sys.path.append(r"")
sys.path.append(r"")

import faulthandler
faulthandler.enable()
import os
import numpy as np
import scipy.integrate as si
import scipy.special as sc
from scipy.stats import truncnorm, norm
import csv
from scipy.optimize import fsolve, fmin, fmin_cg, fmin_bfgs
from operator import itemgetter
import math
import pyparsing as pypars
import functools
from functools import partial
import copy
import pandas as pd
import multiprocess
from multiprocess import Process

#defining the paths for folders
pathsep=os.sep
base_folder = r''
preprod_folder   = os.path.join(base_folder, 'run') + pathsep
model_folder     = os.path.join(base_folder, 'model') + pathsep
param_folder     = os.path.join(base_folder, 'parameters') + pathsep
outputs_folder   = os.path.join(base_folder, 'outputs') + pathsep

#load functions for chebyshev interpolation
from all_functions import nodes_Chebyshev, Chebyshev

#load functions for Bellman equation
from all_functions import Bell_max

#load functions for dynamics
from all_functions import utility, marginal_utility, production, abatement_cost, damage_factor, emissions, \
    growth_primary_tipping, maximand, partialc_maximand, law_motion_SEU

#load functions for state space
from all_functions import law_motion_state_space, compute_intervals

#load functions for stochastic risk 
from all_functions import density_betalaw

#load functions for interpolation
from all_functions import quadrature_product_terminal_EV, EV_quadrature_stochastic2, integrand1, integrand_terminal1, terminal_value


#parallel task (on chebyshev nodes)
#interpolation from T to initial period
#solve the bellman problem at one chebyshev node
#and compute the corresponding chebyshev coefficients

def task(i, t_step, params, coefn, z, u_adjust):

    # Initial guess: fixed value at terminal period, warm start otherwise
    if t_step == params["T"]:
        guess = np.array((0.1)).reshape(1, 1)
    if t_step != params["T"]:
        guess = params["guess"][i]

    # Initialize state vector over the full horizon
    # map chebyshev nodes from [-1:1] to the admissible state space
    u = np.zeros((params["T"], params["dim"]))
    if params["dim"]==1:
        for j in range(0, params["dim"]):
            u[t_step - 1] = (z[[i], j] + 1) / 2 * (params['u_max'][t_step - 1] - params["u_min"][t_step - 1] + 2 * u_adjust[j]) + params["u_min"][t_step - 1] - u_adjust[j]
    else:
        for j in range(0,params["dim"]):
            u[[t_step - 1], j] = (z[[i], j] + 1) / 2 * (params["u_max"][[t_step - 1], j] - params["u_min"][[t_step - 1], j] + 2 * u_adjust[j]) + params["u_min"][[t_step - 1], j] - u_adjust[j]

    # Solve the Bellman maximization problem using Nelder�Mead (via fmin)
    (xopt, fopt, nb_iter, nb_funcalls, exitflag) = fmin(Bell_max, guess, args=(t_step,u, params, coefn), xtol=params["tolintabs"], ftol=params["tolintrel"], full_output=1, disp=0)  # xopt is the minimum and fopt the value of the function at its minimum

    # Check for numerical issues at the current Chebyshev node
    if np.any(np.isnan(xopt)) or np.any(np.isinf(xopt)) or np.isnan(fopt) or np.isinf(fopt):
        print(f"Warning: NaN or Inf detected at node {i}, t_step {t_step}")
        print("xopt:", xopt)
        print("fopt:", fopt)

    # Enforce admissible bounds on the control variables
    xopt = np.maximum(np.zeros(np.shape(xopt)), xopt)
    xopt = np.minimum(np.ones(np.shape(xopt)), xopt)

    # Evaluate the approximated value function at the optimum
    args = (t_step,u, params, coefn)
    approx_value = -Bell_max(xopt, *args)

    # Compute Chebyshev coefficients associated with the value function
    coef_V2 = np.zeros(params["nb_terms_cheb"])
    for j in range(params["nb_terms_cheb"]):
        w = np.atleast_2d(params["weight_cheb"][j])
        zi = np.atleast_2d(z[i])
        coef_V2[j] = (2 ** params["d_cheb"][j, 0] / params["nodes"]) * approx_value * Chebyshev(w, zi)

    # Store Chebyshev coefficients and optimal controls
    coef_V3 = np.zeros((params["nb_terms_cheb"], 2))
    coef_V3[:, 0] = coef_V2   
    coef_V3[:, 1] = xopt.copy()  
    return coef_V3


if __name__ == '__main__':

    #defining the study and the set of runs to compute or analyze #check is study well defined
    exec(open(preprod_folder+'study.py').read())
    table_study=np.loadtxt(preprod_folder+study+'.csv',dtype=str,delimiter=';')
    table_param_names=table_study[0,2:]
    run_id=np.loadtxt(preprod_folder+study+'.csv',dtype=int,delimiter=';',skiprows=1,usecols=(0,))
    if table_study[0,1]!='prev_run':
        print("You need to define the previous runs in the second column of the study file, please correct file "+study+".csv")
        exit()
    prev_run_id=np.loadtxt(preprod_folder+study+'.csv',dtype=int,delimiter=';',skiprows=1,usecols=(1,))
    table_param_values_temp=np.loadtxt(preprod_folder+study+'.csv',dtype=None,delimiter=';',skiprows=1)
    table_param_values=table_param_values_temp[:,2:]

    study_file = preprod_folder + study + '.csv'
    df_study = pd.read_csv(study_file, delimiter=';')
    table_param_values = df_study.iloc[:, 2:].select_dtypes(include=[np.number]).to_numpy()
    calibration_ids = df_study['scenario'].tolist()

    nb_digits_run_id=4
    to_run=int(to_run)

    print(to_run)

    exec(open(preprod_folder+'check_study.py').read())

    #finding the previous run results to use for approximation interval
    dim=int(table_param_values[run_id==to_run,table_param_names=='dim'].item())
    prev_run_id = np.loadtxt(preprod_folder + study + '.csv', dtype=int, delimiter=';', skiprows=1, usecols=(1,))

    if dim!=1:
        prev_run=int(prev_run_id[run_id==to_run].item())
        #checking if the output folder for previous run exists
        list_outputs=os.listdir(outputs_folder)
        already_done= [s for s in list_outputs if study+'_run'+str(prev_run).zfill(nb_digits_run_id) in s]
        if not already_done:
            print("Cannot compute run number "+str(to_run))
            print("Previous run outputs not available, please compute before the previous run, run number "+str(prev_run))
        elif np.shape(already_done)[0]!=1:
            print("There are more than one output folder for previous run, run number "+str(prev_run)+". Please clean the outputs folder.")
        else:
            prev_outputs= already_done[0]
            if os.path.isfile(outputs_folder+prev_outputs+pathsep+'simulate_break.csv'):
                print("Cannot compute run number "+str(to_run))
                print("Simulation of previous run (run number "+str(prev_run)+") has broken")
            else:
                prev_run_folder=outputs_folder+prev_outputs+pathsep

    #define run folder (store parameters, store results)
    run_folder = outputs_folder + study + '_run' + str(to_run).zfill(nb_digits_run_id) + pathsep
    os.makedirs(run_folder)

    #Create some files and define some values (especially approximation interval)
    #storing the parameters in the run folder
    exec(open(preprod_folder+'create_parameters_files.py').read())

    #define some values (functions and parameters)
    exec(open(model_folder + 'read_parameters.py').read())
    params["stochastic"]=int(table_param_values[run_id==to_run,table_param_names=='stochastic'].item())

    #add some values to iterable dictionary
    params["run_folder"]=run_folder
    try:
        params["prev_run_folder"] = prev_run_folder if prev_run_folder else "NA"
    except NameError:
        params["prev_run_folder"] = "NA"

    #run approx_interval
    from all_functions import compute_intervals
    params["u_min"], params["u_max"] = compute_intervals(params)

    #saving the boundary of interpolation
    np.savetxt(run_folder+'u_min.csv', params["u_min"], delimiter=";", header='inferior boundary of interpolation dim '+str(dim)+' stochastic '+str(params["stochastic"]))
    np.savetxt(run_folder+'u_max.csv', params["u_max"], delimiter=";", header='superior boundary of interpolation dim '+str(dim)+' stochastic '+str(params["stochastic"]))

    # Simplicial Chebyshev degrees
    if params["dim"] == 1:
        degb = int(params["deg1"])
    if params["dim"]==3:
        degb = int(np.max([params["deg1"], params["deg2"], params["deg3"]]))

    #Nodes for simplicial Chebyshev
    if params["dim"] == 1:
        nodes = int((params["deg1"] + 1))
    if params['dim']==3:
        nodes = int((params["deg1"] + 1) * (params["deg2"] + 1) * (params["deg3"] + 1))

    params['degb'] = degb
    params['nodes'] = nodes

    #import some functions
    from all_functions import n_terms_cheb, ind_cheb, fill_cheb, fill_cheb2

    # matrix of weight (size: nb_terms_cheb *d)
    # each line correspond to a multi-indice
    weight_cheb = fill_cheb2(fill_cheb(params['degb'], params['dim']), params['dim'], params['deg1'], params['deg2'], params['deg3'])
    nb_terms_cheb = len(weight_cheb)
    params['nb_terms_cheb'] = nb_terms_cheb
    params['weight_cheb'] = weight_cheb

    # indicatrice of the multi-indice (size nb_terms_cheb *1)
    # = number of positive indices
    d_cheb = ind_cheb(params['weight_cheb'])
    params['d_cheb'] = d_cheb
    coef_V2 = np.zeros((params['T'], params['nb_terms_cheb']))
    coef_V = np.zeros((params['T'], params['nb_terms_cheb']))

    # setting the boundary of approximation
    u = np.zeros((params["T"], params["dim"]))

    #Depends on the dimension for the simplicial algorithm
    z_adjust = np.zeros((1, params["dim"]))
    z_adjust[:, 0] = -np.cos(np.pi / (2 * (params["deg1"] + 1)))
    if dim>=2:
        z_adjust[:, 1] = -np.cos(np.pi / (2 * (params["deg2"] + 1)))
        z_adjust[:, 2] = -np.cos(np.pi / (2 * (params["deg3"] + 1)))
    params['z_adjust']=z_adjust
    z = nodes_Chebyshev(params)  # shape = (nodes, dim)

    #initialize guess for control variables
    params["guess"] = np.zeros((1, params["nodes"])) # one control var but we want to keep for all nodes

    print("nb_models1")
    print(params["nb_models"])
    print(params["id_climate"])


    for t_step in list(range(params["T"], 0, -1)):
        #define number of parallel CPU / standard = all CPU available
        u_adjust = np.zeros(params["dim"])
        if params["dim"] == 1:
            for j in range(params["dim"]):
                u_adjust[j] = (params['z_adjust'][:, j] + 1) / (-2 * params['z_adjust'][:, j]) * \
                              (params['u_max'][t_step - 1] - params['u_min'][t_step - 1])
        else:
            for j in range(params["dim"]):
                u_adjust[j] = (params['z_adjust'][:, j] + 1) / (-2 * params['z_adjust'][:, j]) * \
                              (params['u_max'][[t_step - 1], j] - params['u_min'][[t_step - 1], j])

        if t_step == params["T"]:
            coefn = None
        if t_step<params["T"]:
            coefn = coef_V[[t_step], :]

        task_with_params = partial(task, t_step=t_step, params=params, coefn=coefn, z=z, u_adjust=u_adjust)

        with multiprocess.Pool(processes=128) as pool:
            coef_V2 = pool.map(task_with_params, range(params["nodes"]))
            coef_V2 = np.array(coef_V2)  # shape (nodes, T, nb_terms_cheb, 2)
            coef = np.sum(coef_V2[:, :, 0], axis=0)  # shape = (nb_terms_cheb,)
            params["guess"] = coef_V2[:, 1, 1]  # keeps only the current period's guesses

            if t_step == T:
                coef_V = np.zeros((params["T"] - 1, params["nb_terms_cheb"]))
                coef_V = np.vstack([coef_V, coef])
            else:
                coef_V[t_step - 1, :] = coef
            np.save(params["run_folder"] + 'coef_V.npy', coef_V)

    #simulation is OK if deterministic
    #if stochastic, would need multiple draws
    if params["stochastic"]==0:
        exec(open(model_folder+'simulate.py').read()) #simulate for each run.
        np.savetxt(params["run_folder"] + 'state_V_notstochastic.csv', u, delimiter=";",header='state variables dim ' + str(params["dim"]) + ' stochastic ' + str(params["stochastic"]))
        np.savetxt(params["run_folder"] + 'control_V_notstochastic.csv', x, delimiter=";", header='control variables dim ' + str(params["dim"]) + ' stochastic ' + str(params["stochastic"]))

    print("nb_models")
    print(params["nb_models"])
    print(params["id_climate"])



