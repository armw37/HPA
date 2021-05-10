# -*- coding: utf-8 -*-
"""
Created on Wed Mar 31 08:28:48 2021

@author: mskwh
"""

import multiprocessing
import numpy as np
from copy import deepcopy
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve
import lmfit
import csv
import matplotlib.pyplot as plt


"""
Building up to a good fit

"""

##########################################################################
# Full HPA System                                                        #
##########################################################################


# 0
def crf(t, x, v):
    rhs = (
           v['k1']  *  
           # Positive hill feedback by cort_p
           (1 + v['zeta'] * (((x[2]/v['v_p']) ** v['alpha'])
             / ((x[2]/v['v_p']) ** v['alpha'] + v['c'] ** v['alpha']))
           # Negative hill feedback by cort_p
           - (v['psi'] * (((x[2]/v['v_p']) ** v['gamma'])
             / ((x[2]/v['v_p']) ** v['gamma'] + v['c3'] ** v['gamma']))))
           # Generic waste
           - v['w1'] * x[0]
           # following the bolus
          )
    return rhs


# 1
def acth(t, x, v):
    rhs = (
           v['k2'] *
           # Negative hill feedback by cort_p
           (1 - v['rho'] * ((x[2]/v['v_p']) ** v['alpha'])
             /((x[2]/v['v_p']) ** v['alpha'] + v['c'] ** v['alpha']))
           # Saturable CRH Stimulation
           #* (((x[0]/v['v_p']) ** v['omega'])
           #/ ((x[0]/v['v_p']) ** v['omega'] + v['cc'] ** v['omega']))
            * v['v_p']  * (x[0]/v['v_p'])
           # Generic waste
           - v['w2'] * x[1]
          )
    return rhs


##########################################################################
# Plasma Cortisol functions                                              #
##########################################################################

# 2
def cort_p(t, x, v):
    rhs = (
           # stimulation
           v['k3'] * (x[1]/v['v_p']) * v['v_p']
           - v['ucp'] * x[2]    
          )
    return rhs


##########################################################################
# Bolus                                                                  #
##########################################################################

# Bolus
Bolps = {'A' : 36.840467,'B': 19.50568645, 'alp': 0.1118, 'bet': 0.01359}

def CRF_bolus(t, v): # Units ng/ml
    """ CRF Bolus Function """
    crf = v['A']*np.exp(-v['alp'] *t) +\
            v['B']*np.exp(-v['bet'] *t)
    return crf


##########################################################################
# SYS                                                                    #
##########################################################################

def sys(t, x, v):
    CRH, ACTH, CORT = crf(t, x, v), acth(t, x, v), cort_p(t, x, v)
   
    return [CRH, ACTH, CORT]


##########################################################################
# Parameters                                                             #
##########################################################################


params = lmfit.Parameters()
# Constraints and extras
params.add('x0bar', value=24984.7752, min=22486, max=27483)            
params.add('x1bar', value=68496.12, min=55000,
           max=88718.784)
params.add('x2bar', value=9980.8632, min=6898.5378,
           max=12000)


params.add('v_p', value=3261.72, vary=False)

# CRH Equations
params.add('alpha', value=7, min=1, max=8)
params.add('c', value=3, min=0.5, max=7.5)
params.add('zeta', value=2, max=2, min=0.5)
params.add('c3', value=1.42, min=0.5, max=7.5)
params.add('gamma', value=3, min=1, max=8)
params.add('psi', value=0.5, min=0.1, max=0.9)
params.add('w1', value=0.02, min=1e-4, max=1) #value=0.1731, min=1e-4, max=1)
params.add('k1', value = 1000, min=1, max=5e3)

# ACTH
params.add('rho', value=0.8, min=0.1, max=0.9)
params.add('w2', value=0.0348, min=1e-4, max=1)
params.add('omega', value=3, min=1, max=8)
params.add('cc', value=7, min=5, max=15)
params.add('k2', value=5, min=1e-5, max= 20) #value=0.1271, min=1e-5, max= 10

# Cort Plasma
params.add('ucp', value=0.0016, min= 1e-6, max=0.1)
params.add('k3', value=0.01, min=1e-6, max=0.1)        


##########################################################################
# Generic functions                                                      #
##########################################################################

def toCortT_ugdl(x, v):
    totalCort = x / (v['v_p'].value * 0.04) * 0.1
    return totalCort

def get_csv_data(file_str, first_row_header=True):
    """ Read in the csv data """
    with open(file_str) as ac:
        reader = csv.reader(ac, delimiter=',')
        if first_row_header:
            next(reader)
       
        arr = np.array([row for row in reader]).astype('float64')
    return arr


def sample_with_copy(parms):
    ps = deepcopy(parms)
    for p in ps:
        if ps[p].vary != False:
            ps[p].value = np.random.uniform(ps[p].min, ps[p].max)
    return ps


def toConc(x, v):
    conc = x/v['v_p']
    return conc

##########################################################################
# Galluci                                                                #
##########################################################################


######################### model and objective ############################

def g(fun,t, ts, x0, ps):
    """ Wrapper for solving an IVP """
    sol = solve_ivp(fun, t, x0, method='LSODA',
            t_eval=ts, args=(ps,))
    return sol


def residual_galluci(ps, fun, t, ts, data, x0, indices):
    """
    Galluci objective function
    ps: Parameters object
    fun: system function
    t: interval
    ts: evaluation times
    data: experimental data points
    """
    z = g(fun, t, ts, x0, ps)
    if z.status != 0:
        print(z.message)
    # steady state mass
    v = ps.valuesdict()
    ssM = np.array([ps['x0bar'].value, ps['x1bar'].value,
                   ps['x2bar'].value])
    ss = [crf(0, ssM, v), acth(0, ssM, v), cort_p(0, ssM, v)]
    # dynamics 
    CRH = (z.y.T[indices[0],0]/ps['v_p'].value)/1000
    ACTH = z.y.T[indices[1],1]/ps['v_p'].value
    CORT = toCortT_ugdl(z.y.T[indices[2],2], ps)
    model = np.column_stack((CRH, ACTH, CORT))
    dataarr= (model - data).ravel()
   
    return np.concatenate([dataarr, ss])


# Curve fitting function
def curveFunction(cost_func, params, func_args, op_method):
    """ Wrapper for fitting a model """
    try:
        minner = lmfit.Minimizer(cost_func, params, fcn_args=(func_args))
        result = minner.minimize(method=op_method)
    except RuntimeError as e:
        print(e)
        result = None
    return result


##########################################################################
# Galluci                                                                #
##########################################################################


def plot_sols(t, sols, labels1, title1):
    """
    Plot mutliple sols of the hpa system
    """
    fig, ax = plt.subplots(1, 3, figsize=(9,7), sharex=True)
    lines = []
    for i in range(sols.shape[-1]):
        ax[0].plot(t, sols[:,0,i])
        ax[0].set_title('CRH (pg/ml)')
       
        ax[1].plot(t, sols[:,1,i])
        ax[1].set_title('ACTH (pg/ml)')        
       
        line, = ax[2].plot(t, sols[:,2,i])
        ax[2].set_title('CortP (ng/ml)')
       
        lines.append(line)

    leg1 = fig.legend(tuple(lines), labels1, title=title1 ,
            bbox_to_anchor=(1.005,0.85),loc='upper right')
    fig.tight_layout()
    fig.subplots_adjust(right=0.83)
    return [fig, leg1]


##########################################################################
# Main                                                                   #
##########################################################################

# Main
def main():
   
    ############################# file data ##############################
    acthFile = ('/home/anderson/Downloads/biphasic_ACTH.csv')
    cortFile = ('/home/anderson/Downloads/biphasic_Cort.csv')

    acthData = get_csv_data(acthFile); cortData = get_csv_data(cortFile)

    # Select the 3 ug/kg test
    acthD = acthData[:,4:]; cortD = cortData[:,4:]

    ############################ Master data #############################
    # set up initial number of runs
    num_runs = 4
    stardata = []
   
    # normalize data to the mean steady state of run m
    ssACTH = 21
    ssCortT = (3.06 / 0.04) * 0.1

    # Add the steady state
    acthD[:,1] = ssACTH + acthD[:,1]
    cortD[:,1] = ssCortT + cortD[:,1]

    # Change time to minutes
    acthD[:,0] = acthD[:,0] * 60
    cortD[:,0] = cortD[:,0] * 60


    # Add bolus data at equal increments over total length
    tspan = [acthD[0,0], acthD[-1,0]]
    ts = np.linspace(tspan[0], tspan[-1], acthD.shape[0])
    bolus = CRF_bolus(ts, Bolps)
    
    # unique time points to evaluate at
    uniquets, index = np.unique(np.concatenate([ts, acthD[:,0], 
        cortD[:,0]]), return_inverse=True)
    indices = [index[:acthD.shape[0]], 
        index[acthD.shape[0]:2*acthD.shape[0]],
        index[2*acthD.shape[0]:]]
    
    # Data 
    time_data = np.column_stack((bolus, acthD[:,1], cortD[:,1]))

    # Set up Initial conditions
    x0 = [7.66, 21*params['v_p'].value, 3.06* params['v_p'].value]
    x0[0] = CRF_bolus(0, Bolps) * 1000 * params['v_p'].value
   
    residual_galluci(params, sys, tspan, uniquets, time_data, x0, indices)
    # Create the data necessary for starmap
    for j in range(num_runs):
        p = sample_with_copy(params)
        op_method = 'leastsq'
        func_args = (sys, tspan, uniquets, time_data, x0, indices)
        CurveFitVals = (residual_galluci, p, func_args, op_method)
        stardata.append(CurveFitVals)
   

   
    ############################# Multiprocess ##############################
   
    cpus = multiprocessing.cpu_count()
    with multiprocessing.Pool(2*cpus-1) as pool:
        out = pool.starmap(curveFunction, stardata)
        pool.close()
        pool.join()
   
    # do top 10 fits
    res_arr = np.zeros((num_runs,2))
    for i, o in enumerate(out):
        res_arr[i,0], res_arr[i,1] = i, o.chisqr
   
    top5 = res_arr[res_arr[:,1].argsort()][:5,:]
    for i in top5[:,0]:
        lmfit.report_fit(out[int(i)])
       
        # Plot fit
        p1, = plt.plot(ts, time_data[:,0], 'o', label='CRH Data ng/ml')
        p2, = plt.plot(ts, time_data[:,1], 'o', label='ACTH Data pg/ml')
        p3, = plt.plot(ts, time_data[:,2], 'o', label='Cort Data ug/dl')
       
       
        final = out[int(i)].residual[:45].reshape(time_data.shape) + time_data
       
        p4, = plt.plot(ts, final[:,0], label='CRH model ng/ml')
        p5, = plt.plot(ts, final[:,1], label='ACTH model pg/ml')
        p6, = plt.plot(ts, final[:,2], label='Cort model ug/dl')
       

        plt.legend(handles=[p1,p2,p3,p4,p5,p6],
            bbox_to_anchor=(1.1, 1.05), fontsize='small')
        plt.xlabel('Time (min)'); plt.ylabel('Concentration')
        plt.show(block=True)
           
        # check steady state
        ps = out[int(i)].params.valuesdict()
       
        # Create the solver function
        y = lambda x : sys(0, x, ps)
   
        # ensure ICS
        # Set up Initial conditions
        x0 = [7.66, 21*params['v_p'].value, 3.06* params['v_p'].value]
        x0[0] = CRF_bolus(0, Bolps) * 1000 * params['v_p'].value
   
        # Solve
        steady, info, ier, mesg  = fsolve(y, x0,
                                                  full_output=True)
       
        if ier == 1:
            print('Steady State 1: {:<12.2f} {:<12.2f} {:<12.2f}'.format(
                *toConc(steady, ps)))
            print()
           
        else:
            print(mesg)
            print()

       
        # plot longer sol
        # initial CRH, last ps is fine for v_p since it doesnt change
        CRH_IC = np.array([8,100, 500, 1000]) * ps['v_p']
        # tspan and evals
        tspan = [0, 500]; num_evals = 1500
        # set title
        title = 'CRH: '
   

        sols = np.zeros((1500, 3, len(CRH_IC)))
        evals = np.linspace(tspan[0],tspan[1],num_evals)
        ic = steady
        for j, i in enumerate(CRH_IC):
            ic[0] = i
            sol = solve_ivp(sys, tspan, ic, method='LSODA', t_eval=evals,
            args=(ps,))
            sols[:,:,j] = sol.y.T

        labels = [str(int(i/ps['v_p']))  for i in CRH_IC]
        ConcSols = toConc(sols, ps)
        fig = plot_sols(evals, ConcSols, tuple(labels), title)
        fig[0].subplots_adjust(bottom=0.1)
        fig[0].text(0.5, 0.02, 'Time (min)', ha='center')
        plt.show(block=True)
   
   
    return 0

if __name__ == '__main__':
    main()


