# -*- coding: utf-8 -*-
#this file is used to simulate control and state variables for deterministic program.

#initializing some value
#guess for minimizing Bell_max
#initialize 
if params["dim"]==1:
    guess0=(1 - (params["bet"] * params["alpha"]))*np.ones((1,1))
    u0trend=params["K0"]/(((params["A0"])**(1/(1-params["alpha"])))*params["POP0"])*np.ones((T+2,1))
    u0=u0trend[[0],:]

if params["dim"]==3:
    guess0=np.array((0.1)).reshape(1,1)
    u0=np.concatenate((params["K0"]/((params["A0"]**(1/(1-params["alpha"])))*params["POP0"])*np.ones((1,1)),params["T0"]*np.ones((1,1)),params["TREE0"]*np.ones((1,1))),axis=1)

#Declaration of matrixes to store results
x=np.zeros((T,1))
u=np.zeros((T,dim))
shock_save1=[0]*params["T"]
shock_save2=[0]*params["T"]
shock_save2_model=[0]*params["T"]

#Initialization of state variables
u[[0],:]=u0

#Iteration
guess=guess0
for t_step in list(range(1,params["T"])):
    
    #guess for control
    if t_step==1:
        guess=guess0
    else:
        guess=x[[t_step-2]]

    #Compute optimal control variable (consumption and abattement) at time t by the minimization of Bellman function
    (xopt,fopt, nb_iter, nb_funcalls, exitflag)=fmin(Bell_max,guess, args=(t_step,u, params, coef_V[[t_step],:]), xtol=params["tolintabs"], ftol=params["tolintrel"],full_output=1,disp=0)#xopt is the minimum and fopt the value of the function at its minimum
    guess=xopt.copy()
    xopt=np.maximum(np.zeros(np.shape(xopt)),xopt)#transpose because x is vector
    xopt=np.minimum(np.ones(np.shape(xopt)),xopt)
    x[[t_step-1]]=xopt

    # random stochastic draw
    if params["stochastic"]==0:
        u[[t_step], :] = law_motion_SEU(t_step, x[[t_step - 1]], u[[t_step - 1], :], params["TCRE_mean"],params["climate_meanEU"], 0, params)

    if params["stochastic"]==1:
        shock1 = truncnorm.rvs(params["a_tcre"], params["b_tcre"], loc=params["TCRE_mean"], scale=params["sdeviation"], size=1)
        shock_save1[t_step] = shock1
        u[[t_step], :] = law_motion_SEU(t_step, x[[t_step - 1]], u[[t_step - 1], :], shock1, 0, 0, params)

    if params["stochastic"]==2:
        shock1 = truncnorm.rvs(params["a_tcre"], params["b_tcre"], loc=params["TCRE_mean"], scale=params["sdeviation"], size=1)
        shock_save1[t_step] = shock1
        shock2 = random.beta(params["beta_alpha"], params["beta_beta"])
        shock_save2[t_step] = shock2
        u[[t_step], :] = law_motion_SEU(t_step, x[[t_step - 1]], u[[t_step - 1], :], shock1,params["climate_meanEU"], 0, params)
            
    if params["dim"]==1:
        if (u[[t_step], :] > params["u_max"][t_step]).any() or (u[[t_step], :] < params["u_min"][t_step]).any():
            out_of_bounds = 1
            break
    else:
        if (u[[t_step],:]>params["u_max"][[t_step],:]).any() or (u[[t_step],:]<params["u_min"][[t_step],:]).any():
            out_of_bounds=1
            break

