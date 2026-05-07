**Replication package for Fillon & Guivarch, "The need for regulation of climate subsystems"**

Data    
Most data are in analysis_figures/data/  
Some heavy datasets are available only from this Google Drive and should then be placed in analysis_figures/data with these names rainfall, temperature and earthdata_carbon2010  
These datasets can also be downloaded from [ISIMIP](https://data.isimip.org/) and [EarthData](https://www.earthdata.nasa.gov/data/catalog/ornl-cloud-global-maps-c-density-2010-1763-1)  
  
**Step 1 - calibration of Amazon dynamics**  
**Step 1A**  
run analysis_figures/calibration_dynamics/estimation.R  
this step uses temperature and precipitation data from ISIMIP 2B  
compute for each climate model the coefficient of climate shock  

**Step 1B**  
run analysis_figures/calibration_dynamics/calibration.R  
this step uses previous coefficients to match expert estimates  
three calibration: one benchmark, two counterfactuals  

**Step 2 - solve optimal intertemporal program**  
We run 6 specificiations as benchmark:  
run 1: deterministic (stochastic=0), without temperature damages, exogenous paths of variables without control  
run 2: deterministic (stochastic=0) with temperature damages, without Amazon  
run 3: deterministic (stochastic=0) with temperature damages, with Amazon  
run 4: stochastic=1 (one risk dT/dS only) with temperature damages, without Amazon  
run 5: stochastic=1 (one risk dT/dS only) with temperature damages, with Amazon  
run 6: stochastic=2 (two risks dT/dS and dA/dT) with temperature damages, with Amazon  
run 7 to 11 : first counterfactual - same specifications as above, but with a different calibration from expert elicitations  
run 12 to 16 : second counterfactual  - same specifications as above, but with a different calibration from expert elicitations  

**Step 2A - interpolation of value function**   
open study.csv and define the run number  
run 1 then 2 should be run first because state space for all runs is created at run 2  
run meta_run1_interpolation.py  

**Step 2B - simulation of stochastic paths**   
only for runs 4 to 6 (for counterfactuals - 9 to 11, 14 to 16)  
open meta_run2_simulation.py once chebyshev meta_run1_interpolation.py is done for this run  
define the run number   
run meta_run2_simulation.py  




