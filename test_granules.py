import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

# worth a try

def acthSecretion(t, x, v, CRH):
    """
    CRHDelayM: x[0]
    M: x[1]
    CRHDelayD: x[2]
    D: x[3]
    CRHDelayF: x[4]
    F: x[5]
    ACTH: x[6]
    """



    #mhat
    mhat = (v['m'] * (1 - (x[1]/v['msat']))
            * ((x[0] - v['cm'])**v['nm']
                /((x[0] - v['cm'])**v['nm'] 
                    + (v['cm50'] - v['cm'])**v['nm'])) 
            * np.heaviside(x[0] - v['cm'], 1))

    #dhat
    dhat = (v['d'] 
            * ((x[2] - v['cd'])**v['nd']
                /((x[2] - v['cd'])**v['nd'] 
                    + (v['cd50'] - v['cd'])**v['nd'])) 
            * np.heaviside(x[2] - v['cd'], 1))

    #fhat
    fhat = (v['f']
            * ((x[4] - v['cf'])**v['nf']
                /((x[4] - v['cf'])**v['nf'] 
                    + (v['cf50'] - v['cf'])**v['nf'])) 
            * np.heaviside(x[4] - v['cf'], 1))

    CRHDelayM = v['tau1'] * (CRH(t,v) - x[0])
    CRHDelayD = v['tau2'] * (CRH(t,v) - x[2])
    CRHDelayF = v['tau3'] * (CRH(t,v) - x[4])

    Mdot = v['M0'] + mhat - (v['D0'] + dhat) * x[1]
    Ddot = (v['D0'] + dhat) * x[1] - (v['F0'] + fhat) * x[3]
    Fdot = (v['F0'] + fhat) * x[3] - v['s'] * x[5]
    ACTH = v['s'] * x[5] - v['w1'] * x[6]

    return [CRHDelayM, Mdot, CRHDelayD, Ddot, CRHDelayF, Fdot, ACTH]


# Bolus
def CRF_bolus(t, v): # Units ng/ml
    """ CRF Bolus Function """
    crf = v['A']*np.exp(-v['alp'] *t) +\
            v['B']*np.exp(-v['bet'] *t)
    return 100 * crf

def CRF_steady(t, v): # Units ng/ml
    return 7.66


# params
# for now I assume as much is the same as possible

ps = {
        'A'    : 36.840467,
        'B'    : 19.50568645, 
        'alp'  : 0.1118, 
        'bet'  : 0.01359,
        'm'    : 2,
        'msat' : 4,
        'cm'   : 5,
        'nm'   : 6,
        'cm50' : 100,
        'd'    : 0.1,
        'cd'   : 5,
        'nd'   : 3,
        'cd50' : 20,
        'f'    : 2,
        'cf'   : 5,
        'nf'   : 2, 
        'cf50' : 10,
        'tau1' : 0.5,
        'tau2' : 0.01,
        'tau3' : 0.05,
        'M0'   : 1,
        'D0'   : 0.1,
        'F0'   : 0.1,
        's'    : 2,
        'w1'   : 0.0348,
     }


# Find its steady state to start at
x0 = [7.66,1,7.66,1,7.66,1,20]
 
# Create the solver function
y = lambda x : acthSecretion(0, x, ps, CRF_steady)

y(x0)
# Solve

def steady_solve(y, x0):
    steady, info, ier, mesg  = fsolve(y, x0, full_output=True)
       
    if ier == 1:
        strFormat = '{:<6.2f} ' * len(x0)
        print('Steady State 1: ' + strFormat.format(*steady))
    else:
        print(mesg)
    return steady

# test different steady CRH

# crf steady state values
CRF = [7.66, 10, 15, 20, 25, 30, 50]
steadys = []
for c in CRF:
    crf = lambda t, v: c
    y = lambda x : acthSecretion(0, x, ps, crf)
    steadys.append(steady_solve(y, x0))

# Test a bolus



def plot_sols(t, sol):
    """
    Plot mutliple sols of the hpa system
    """
    fig, ax = plt.subplots(2, 3, figsize=(9,7), sharex=True)
    ax = ax.ravel()
    for i in range(len(ax)):
        ax[i].plot(t, sol[:,i+1])
       
    fig.tight_layout()
    return fig

sol = solve_ivp(acthSecretion, [0, 400], steadys[0], 
        method='LSODA', args=(ps,CRF_bolus))

plot_sols(sol.t, sol.y.T)
plt.show(block=True)
