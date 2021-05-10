import numpy as np
import lmfit
from lmfit.lineshapes import gaussian, lorentzian
import multiprocessing
import copy
import time

def residual(pars, x, data):
    model = (gaussian(x, pars['amp_g'], pars['cen_g'], pars['wid_g']) +
             lorentzian(x, pars['amp_l'], pars['cen_l'], pars['wid_l']))
    return model - data

# Generate Data
np.random.seed(0)
x = np.linspace(0, 20.0, 601)
data = (gaussian(x, 21, 6.1, 1.2) + lorentzian(x, 10, 9.6, 1.3) +
        np.random.normal(scale=0.1, size=x.size))

# Parameters
pfit = lmfit.Parameters()
pfit.add(name='amp_g', value=10, min=0.1, max=20)
pfit.add(name='amp_l', value=10, min=0.1, max=20)
pfit.add(name='cen_g', value=5, min=0.1, max=10)
pfit.add(name='peak_split', value=2.5, min=0, max=5, vary=True)
pfit.add(name='cen_l', expr='peak_split+cen_g')
pfit.add(name='wid_g', value=1, min=0.1, max=2)
pfit.add(name='wid_l', expr='wid_g')


def uniform_p_sample(parms):
    ps = copy.deepcopy(parms)
    for p in ps:
        if ps[p].vary != False:
            ps[p].value = np.random.uniform(ps[p].min, 
                    ps[p].max)
    return ps



# Creating curve fit function called
def curvefunction(residual, p, x, data):
    try:
        minner = lmfit.Minimizer(residual, p, fcn_args=(x, data))
        result = minner.minimize(method='leastsq')
    except RuntimeError as e:
        print(e)
        result = np.full(5, np.inf)
    return result

# Monte Carlo method to find best fit
d = []; n_runs = 1000

for j in range(n_runs):
    # Create list of tuples for starmap
    pf = uniform_p_sample(pfit)    
    CurvefitValues = (residual, pf, x, data)
    d.append(CurvefitValues)

if __name__ == '__main__':
    # Use Pool and curvefit to get all the parameters
    cpus = multiprocessing.cpu_count()
    with multiprocessing.Pool(2*cpus-1) as pool:
        start_time = time.time()
        out = pool.starmap(curvefunction, d)
        pool.close()
        pool.join()
        print("--- %s seconds ---" % (time.time() - start_time))




