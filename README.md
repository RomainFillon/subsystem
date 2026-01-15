*Replication package for Fillon & Guivarch, "The need for regulation of climate subsystems"*

#Data


#Step 1 - calibration of Amazon dynamics
run analysis_figures/calibration_dynamics/analysis_calibration.R
this step use precipitation and temperature data from ISIMIP 2B
compute for each climate model the coefficient of climate shock
then, internal calibration to match Kriegler estimates

#Step 2 - solve optimal intertemporal program
We run 6 specificiations:
run 1: deterministic (stochastic=0), without temperature damages, exogenous paths of variables without control
run 2: deterministic (stochastic=0) with temperature damages, without Amazon
run 3: deterministic (stochastic=0) with temperature damages, with Amazon
run 4: stochastic=1 (one risk dT/dS only) with temperature damages, without Amazon
run 5: stochastic=1 (one risk dT/dS only) with temperature damages, with Amazon
run 6: stochastic=2 (two risks dT/dS and dA/dT) with temperature damages, with Amazon

#Step 2A - interpolation of value function
open study.csv and define the run number
run 1 then 2 should be run first because state space for all runs is created at run 2
run meta_run1_interpolation.py

#Step 2B - simulation of stochastic paths
only for runs 4 to 6
open meta_run2_simulation.py once chebyshev meta_run1_interpolation.py is done for this run
define the run number 
run meta_run2_simulation.py




