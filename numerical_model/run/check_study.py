# -*- coding: utf-8 -*-

#checks if the study is well-defined
well_defined=1

#checking is all runs have different run_id
if len(set(run_id)) != len(run_id):
    run_id_duplicates=list(set([x for x in run_id if list(run_id).count(x) > 1])) 
    print("There are several runs with the same run_id, "+str(run_id_duplicates)+", please correct first column in study")   
    well_defined=0
    exit()

#checking if the study specific parameters names are all defined in the parameters files
param_gen_names=np.loadtxt(param_folder+'param_gen.csv',dtype=str,delimiter=';',usecols=(0,))
param_model_names=np.loadtxt(param_folder+'param_model.csv',dtype=str,delimiter=';',usecols=(0,))

for p in table_param_names:
    if p not in param_gen_names:
        if p not in param_model_names:
            print(p + " from the study specific parameters names is not a valid parameter name (does not exist in the files in " + param_folder + ")")
            well_defined = 0

#checking if the study specific parameters names contain all the parameters to be checked to find the previous run (list_prev1 and list_perv2)
for p in list_prev1:
    if p not in table_param_names:
        print(p+" from list_prev1 is not in the study specific parameters names")
        well_defined=0

for p in list_prev2:
    if p not in table_param_names:
        print(p+" from list_prev2 is not in the study specific parameters names")
        well_defined=0

try:
    prev_run_id # does a exist in the current namespace
except NameError:
	exec(open(preprod_folder+'def_prev_run.py').read())
    
#same for idealine run (ie dim 1 run)   
ideal_run_id=np.zeros((1))
dim=int(table_param_values[run_id==to_run,table_param_names=='dim'])
if dim!=1:
    run_param_values=table_param_values[run_id==to_run,:]
    ideal_run_param_values=run_param_values
    is_ideal=np.ones(np.shape(run_id))
    ideal_run_param_values[:,table_param_names=='dim']=1
    ideal_run_param_values[:,table_param_names=='stochastic']=0
    for i in list(range(np.shape(is_ideal)[0])):
        for ptest in list_prev1:
            is_ideal[i]=is_ideal[i]*(table_param_values[i,table_param_names==ptest]==ideal_run_param_values[:,table_param_names==ptest])
    if np.sum(is_ideal)==0:
        print("No idealine run defined for run number "+str(to_run))
        well_defined=0
    elif np.sum(is_ideal)>1:
        print("Several idealine runs correpond to run number "+str(to_run)+": runs "+str(run_id[is_ideal==True]))
        #well_defined=0
    else:
        ideal_run_id[0]=run_id[is_ideal==True]
