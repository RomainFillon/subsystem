#import some libraries
import numpy as np
import pandas as pd
import csv
import itertools

#open parameter values
params = {}
with open(run_folder+'param_gen.csv', 'r') as f:
    reader = csv.reader(f, delimiter=';')
    for row in reader:
        key = row[0]
        val = row[1]
        if key in ['stochastic', 'draws', 'deg']:
            params[key] = int(val)
        elif key == 'model':
            params[key] = val
        else:
            params[key] = float(val)

with open(run_folder+'param_model.csv', 'r') as f:
    reader = csv.reader(f, delimiter=';')
    for row in reader:
        key = row[0]
        val = row[1]
        if key in ['dim', 'T', 'deltaT']:
            params[key] = int(val)
        else:
            params[key] = float(val)

param_file = os.path.join(param_folder, 'epsilons.csv')
with open(param_file, 'r') as f:
    reader = csv.reader(f, delimiter=';')  # ; au lieu de ,
    for row in reader:
        if not row or len(row) < 2:
            continue
        key = row[0].strip()
        val = row[1].strip()
        params[key] = float(val)


param_file = os.path.join(param_folder, 'calibration_kriegler.csv')
df = pd.read_csv(param_file)
df_calibration = df[df['calibration_id'] == params['scenario']]
params.update({
    'beta_beta': float(df_calibration['beta_beta']),
    'beta_alpha': float(df_calibration['beta_alpha']),
    'growth0': float(df_calibration['growth0']),
    'eta': float(df_calibration['eta']),
    'beta0': float(df_calibration['beta0']),
    'Upsilon': float(df_calibration['Upsilon']),
})

#adapt to 5y period

params['rho'] = (1 + params['rho1y'])**params['deltaT'] - 1
params['delta'] = 1 - (1 - params['delta1y'])**params['deltaT']
params['bet'] = 1 / (1 + params['rho'])

#exogenous population path
T = params['T']
params['time_horizon'] = np.arange(T + 2).reshape((T + 2, 1))
L = params['POP0'] * np.ones((T + 2, 1))
for t in range(T + 1):
    L[t + 1, 0] = L[t, 0] * (params['POPASYM'] / L[t, 0])**params['GPOP0']
params['L'] = L

#exogenous productivity path
#exogenous decarbonation cost
GA = params['GA0'] * np.exp(-params['DELA'] * params['deltaT'] * params['time_horizon'])
A = params['A0'] * np.ones((T + 2, 1))
for t in range(T + 1):
    A[t + 1, 0] = A[t, 0] / (1 - GA[t, 0])
params['GA'] = GA / params['deltaT']
params['A'] = A
GSIGMA = -params['GSIGMA'] * np.ones((T + 2, 1))
for t in range(T + 1):
    GSIGMA[t + 1, 0] = GSIGMA[t, 0] * (1 - params['DSIG'])**params['deltaT']
sigm = params['SIG0'] * np.ones((T + 2, 1))
for t in range(T + 1):
    sigm[t + 1, 0] = sigm[t, 0] * np.exp(GSIGMA[t, 0] * params['deltaT'])
params['GSIGMA'] = GSIGMA
params['sigm'] = sigm
PBACKTIME = params['PBACK'] * (1 - params['GBACK'])**(params['time_horizon'])
params['theta1'] = PBACKTIME * sigm / params['theta2']

#climate impact on the rainforest dA/dT
shock_mean = params['beta_alpha']/(params['beta_alpha']+ params['beta_beta'])
climate_mean = [shock_mean*params['epsilon_max1'], shock_mean*params['epsilon_max2'], shock_mean*params['epsilon_max3'], shock_mean*params['epsilon_max4']]
climate_max = [params['epsilon_max1'], params['epsilon_max2'], params['epsilon_max3'], params['epsilon_max4']]
params['climate_mean'] = climate_mean
params['climate_max'] = climate_max

params['climate_meanEU']= np.mean(climate_mean)
params['climate_maxEU'] = params['epsilon_max']

#deforestation and degradation scenarii
with open(param_folder+'param_scenariosEU.csv','r') as f:
    df_scenarioEU = pd.DataFrame(list(csv.reader(f, delimiter=';'))).astype(float)
params['df_scenarioEU'] = df_scenarioEU
mean_scenarioEU = params['ratio_deg'] * params['deltaT'] * df_scenarioEU / (params['total_area']*1000)
params['mean_scenarioEU'] = mean_scenarioEU

a = list(range(len(climate_mean))) if params['id_climate']==1 else [0]
b = [0]; c = [0]
params['models'] = list(itertools.product(a, b, c))
params['nb_models'] = len(params['models']) if params['id_climate'] else 1

params['a_tcre'] = (params['TCRE_min'] - params['TCRE_mean']) / params['sdeviation']
params['b_tcre'] = (params['TCRE_max'] - params['TCRE_mean']) / params['sdeviation']
