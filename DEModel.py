import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve
from scipy.interpolate import CubicHermiteSpline as cs
import re
import copy
import pint
import matplotlib.pyplot as plt
import csv
import types
import pprint

class DeEquation:
    def __init__(self, eqStr, var_parm_dict = {}):
        self.eqStr = eqStr
        self.vp_dict = var_parm_dict


    def make_fun(self, funcName, funcHeader,
                 funcReturn, vars_dict):
        """
        deinst: The variable name of the DeEquation instance
        Returns
        -------
        Executable mthematical function given by eqStr,
        vp_dict determines what is a state variable and
        what is a parameter

        """
        
        replicateEq = copy.copy(self.eqStr)
        self.vp_dict[funcName] = vars_dict
        for key, var in vars_dict.items():
            replicateEq = re.sub(key, var, replicateEq)

        execStr = funcHeader + replicateEq + funcReturn
        exec(execStr)
        setattr(self, funcName, 
                types.MethodType(eval(funcName), self))
        #return inst
    

u = pint.UnitRegistry()

##########################################################################
# Parameter dictionary                                                   #
##########################################################################


Params = {# CRH
          'k1'       : 2801.81   * u.pg /  u.min,    # multiplied by volume.   
          'zeta'     : 2,
          'alpha'    : 3,
          'c'        : 3.06      * u.ng / u.ml,
          'psi'      : 0.5,
          'gamma'    : 3,
          'c3'       : 1.42      * u.ng / u.ml,
          'w1'       : 0.1731    * (1 / u.min),
          # ACTH
          'k2'       : 0.1271    * (1 / u.min),  
          'rho'      : 0.5,
          'w2'       : 0.0348    * 1 / u.min,
          'cc'       : 7.66      * u.pg / u.ml,
          'omega'    : 3,
          # Cortisol
          'k3'       : 0.4       * 1 / u.min, 
          'v_p'      : 3261.72   * u.ml,
          'v_i'      : 12153.9   * u.ml, 
          'Dcort_pi' : 10.8      * u.ml / u.min,
          'Dcort_pk' : 9        * u.ml / u.min,
          'ucp'      : 0.0016    * (1 / u.min), 
          'konA'     : 2079      * 1 / ((u.mol/u.l) * u.s),
          'koffA'    : 0.693     * (1 / u.s),
          'MWCort'   : 362.46    * u.g / u.mol,
          'MWAlb'    : 66437     * u.g / u.mol,
          'MWAC'     : 66799.46  * u.g / u.mol,
          'konC'     : 1650      * 1 / ((u.nmol/u.ml) * u.min),
          'koffC'    : 52.8      * 1 / u.min,
          'MWCBG'    : 52000     * u.g / u.mol,
          'MWCC'     : 52362.46  * u.g / u.mol,
          # Albumin
          'k4'       : 0.005     * u.g / u.min, 
          'Klp'      : 0.1,
          'sigA'     : 0.8,
          'Deg_a'    : 0.005     * u.g / u.min,
          'GFR'      : 103.66    * u.ml / u.min,
          'AlbPerm'  : 0, 
          'ModPerm'  : 0,
          'TL_LW_R'  : 7.908     * u.ml / u.min,
          'TM_LW_R'  : 14.41     * u.ml / u.min,
          'TU_LW_R'  : 6.24      * u.ml / u.min,
          'v_tl'     : 3806.42   * u.ml,
          'v_tm'     : 5790.4    * u.ml,
          'v_tu'     : 2554.52   * u.ml,
          'LAPerm'   : 5.4       * u.ml / u.min,
          'MAPerm'   : 6.4       * u.ml / u.min,
          'UAPerm'   : 5.4       * u.ml / u.min,
          'Temp_K'   : 309.772   * u.ml / u.mol, # Just for cancellation purposes
          'IgpCOP'   : 4.49,
          'HetapCOP' : 0,
          'AGPpCOP'  : 0.9426,
          # CBG
          'k5'       : 0.00423   * u.ug / u.min,
          'sigC'     : 0.6,                               
          'half'     : 9.62e-05  * 1 / u.min, 
          # Cort_k
          'v_k'      : 5.024     * u.ml, 
          'Vm2'      : 40        * u.nmol / (u.l * u.min), # no idea here
          'km2'      : 0.04      * u.umol / u.l,         
          'Kd_mr'    : 0.75      * u.nmol / u.l, 
          'Ald'      : 12.1      * u.ng / u.dl,            
          'Vm1'      : 4        * u.nmol / (u.l * u.min), # match / play around a bit
          'km1'      : 21        * u.umol / u.l,           
          'usp'      : 0.00503    * 1 / u.min,              # guesses
          'tonano'   : 1e9,
          'mictonano': 1e3,
          'kcp'      : 0.001     * u.ng/ u.min
          }


# Need to add Tissue albumin
stateVars = {'CRHc'         : 7.66     * u.pg/u.ml, 
             'ACTHc'        : 21       * u.pg/u.ml,
             'Cort_pc'      : 3.06     * u.ng/u.ml,     
             'Cort_ic'      : 3.06     * u.ng/u.ml,
             'Cort_kc'      : 0.121    * u.ng/u.ml,             
             'Cortisone_kc' : 2        * u.ng/u.ml,            
             'Cortisone_pc' : 3.06     * u.ng/u.ml,          
             'Alb_tlc'      : 0.015    * u.g/u.ml,
             'Alb-Cort_tlc' : 380      * u.ng/u.ml,
             'Alb_tmc'      : 0.01     * u.g/u.ml,
             'Alb-Cort_tmc' : 255      * u.ng/u.ml, 
             'Alb_tuc'      : 0.018    * u.g/u.ml,
             'Alb-Cort_tuc' : 450      * u.ng/u.ml,
             'COPEffect'    : 1,
             'Alb_pc'       : 0.045    * u.g/u.ml,
             'Alb-Cort_pc'  : 1150     * u.ng/u.ml,
             'CBG_tlc'      : 8        * u.ug/u.ml,
             'CBG-Cort_tlc' : 270      * u.ng/u.ml,
             'CBG_tmc'      : 5.5      * u.ug/u.ml,
             'CBG-Cort_tmc' : 175      * u.ng/u.ml, 
             'CBG_tuc'      : 10       * u.ug/u.ml,
             'CBG-Cort_tuc' : 335      * u.ng/u.ml,
             'CBG_pc'       : 34       * u.ug/u.ml,
             'CBG-Cort_pc'  : 1180     * u.ng/u.ml,            
             }


stateVarsMass = {'CRH'        : stateVars['CRHc']         * Params['v_p'],
                'ACTH'        : stateVars['ACTHc']        * Params['v_p'],
                'Cort_p'      : stateVars['Cort_pc']      * Params['v_p'],
                'Cort_i'      : stateVars['Cort_ic']      * Params['v_i'],
                'Cort_k'      : stateVars['Cort_kc']      * Params['v_k'],             
                'Cortisone_k' : stateVars['Cortisone_kc'] * Params['v_k'],
                'Cortisone_p' : stateVars['Cortisone_pc'] * Params['v_p'],          
                'Alb_tl'      : stateVars['Alb_tlc']      * Params['v_tl'], 
                'Alb-Cort_tl' : stateVars['Alb-Cort_tlc'] * Params['v_tl'],
                'Alb_tm'      : stateVars['Alb_tmc']      * Params['v_tm'], 
                'Alb-Cort_tm' : stateVars['Alb-Cort_tmc'] * Params['v_tm'],               
                'Alb_tu'      : stateVars['Alb_tuc']      * Params['v_tu'],
                'Alb-Cort_tu' : stateVars['Alb-Cort_tuc'] * Params['v_tu'],
                'COPEffect'   : stateVars['COPEffect'],
                'Alb_p'       : stateVars['Alb_pc']       * Params['v_p'],
                'Alb-Cort_p'  : stateVars['Alb-Cort_pc']  * Params['v_p'],
                'CBG_tl'      : stateVars['CBG_tlc']      * Params['v_tl'],
                'CBG-Cort_tl' : stateVars['CBG-Cort_tlc'] * Params['v_tl'],
                'CBG_tm'      : stateVars['CBG_tmc']      * Params['v_tm'],
                'CBG-Cort_tm' : stateVars['CBG-Cort_tmc'] * Params['v_tm'], 
                'CBG_tu'      : stateVars['CBG_tuc']      * Params['v_tu'],
                'CBG-Cort_tu' : stateVars['CBG-Cort_tuc'] * Params['v_tu'],
                'CBG_p'       : stateVars['CBG_pc']       * Params['v_p'],
                'CBG-Cort_p'  : stateVars['CBG-Cort_pc']  * Params['v_p'], 
                }


##########################################################################
# Female Parameter Dictionary                                            #
##########################################################################

Paramsf = {# CRH
          'k1'       : 2801.81   * u.pg /  u.min,
          'zeta'     : 2,
          'alpha'    : 3,
          'c'        : 3.06      * u.ng / u.ml,
          'psi'      : 0.5,
          'gamma'    : 3,
          'c3'       : 1.42      * u.ng / u.ml,
          'w1'       : 0.1731    * (1 / u.min),
          # ACTH
          'k2'       : 0.1271    * (1 / u.min),  
          'rho'      : 0.5,
          'w2'       : 0.0348    * 1 / u.min,
          'cc'       : 7.66      * u.pg / u.ml,
          'omega'    : 3,
          # Cortisol
          'k3'       : 0.4       * 1 / u.min, 
          'v_p'      : 3261.72   * u.ml,                # this hummod
          'v_i'      : 12153.9   * u.ml,                # this hummod
          'Dcort_pi' : 10.8      * u.ml / u.min,
          'Dcort_pk' : 9         * u.ml / u.min,
          'ucp'      : 0.0016    * (1 / u.min), 
          'konA'     : 2079      * 1 / ((u.mol/u.l) * u.s),
          'koffA'    : 0.693     * (1 / u.s),
          'MWCort'   : 362.46    * u.g / u.mol,
          'MWAlb'    : 66437     * u.g / u.mol,
          'MWAC'     : 66799.46  * u.g / u.mol,
          'konC'     : 1650      * 1 / ((u.nmol/u.ml) * u.min),
          'koffC'    : 52.8      * 1 / u.min,
          'MWCBG'    : 52000     * u.g / u.mol,
          'MWCC'     : 52362.46  * u.g / u.mol,
          # Albumin
          'k4'       : 0.005     * u.g / u.min, 
          'Klp'      : 0.1,
          'sigA'     : 0.8,
          'Deg_a'    : 0.005     * u.g / u.min,
          'GFR'      : 103.66    * u.ml / u.min,
          'AlbPerm'  : 0, 
          'ModPerm'  : 0,
          'TL_LW_R'  : 7.908     * u.ml / u.min,       # this hummod
          'TM_LW_R'  : 14.41     * u.ml / u.min,       # this hummod
          'TU_LW_R'  : 6.24      * u.ml / u.min,       # this hummod
          'v_tl'     : 3806.42   * u.ml,               # this hummod
          'v_tm'     : 5790.4    * u.ml,               # this hummod
          'v_tu'     : 2554.52   * u.ml,               # this hummod
          'LAPerm'   : 5.4       * u.ml / u.min,
          'MAPerm'   : 6.4       * u.ml / u.min,
          'UAPerm'   : 5.4       * u.ml / u.min,
          'Temp_K'   : 309.772   * u.ml / u.mol,
          'IgpCOP'   : 4.49,                           # this hummod
          'HetapCOP' : 0,                              # this hummod
          'AGPpCOP'  : 0.9426,                         # this hummod
          # CBG
          'k5'       : 0.00423   * u.ug / u.min,
          'sigC'     : 0.6,                               
          'half'     : 9.62e-05  * 1 / u.min, 
          # Cort_k
          'v_k'      : 5.024     * u.ml,               # this hummod
          'Vm2'      : 40        * u.nmol / (u.l * u.min),
          'km2'      : 0.04      * u.umol / u.l,         
          'Kd_mr'    : 0.75      * u.nmol / u.l, 
          'Ald'      : 12.1      * u.ng / u.dl,        # this hummod    
          'Vm1'      : 4        * u.nmol / (u.l * u.min), 
          'km1'      : 21        * u.umol / u.l,           
          'usp'      : 0.00503    * 1 / u.min,      
          'tonano'   : 1e3,
          'mictonano': 1e9,
          'kcp'      : 0.001     * u.ng/ u.min
          }


# Need to add Tissue albumin
stateVarsf = {'CRHc'         : 7.66     * u.pg/u.ml, 
             'ACTHc'        : 21       * u.pg/u.ml,
             'Cort_pc'      : 3.06     * u.ng/u.ml,     
             'Cort_ic'      : 3.06     * u.ng/u.ml,
             'Cort_kc'      : 0.121    * u.ng/u.ml,             
             'Cortisone_kc' : 2        * u.ng/u.ml,            
             'Cortisone_pc' : 3.06     * u.ng/u.ml,          
             'Alb_tlc'      : 0.015    * u.g/u.ml,
             'Alb-Cort_tlc' : 380      * u.ng/u.ml,
             'Alb_tmc'      : 0.01     * u.g/u.ml,
             'Alb-Cort_tmc' : 255      * u.ng/u.ml, 
             'Alb_tuc'      : 0.018    * u.g/u.ml,
             'Alb-Cort_tuc' : 450      * u.ng/u.ml,
             'COPEffect'    : 1,
             'Alb_pc'       : 0.045    * u.g/u.ml,
             'Alb-Cort_pc'  : 1150     * u.ng/u.ml,
             'CBG_tlc'      : 8        * u.ug/u.ml,
             'CBG-Cort_tlc' : 270      * u.ng/u.ml,
             'CBG_tmc'      : 5.5      * u.ug/u.ml,
             'CBG-Cort_tmc' : 175      * u.ng/u.ml, 
             'CBG_tuc'      : 10       * u.ug/u.ml,
             'CBG-Cort_tuc' : 335      * u.ng/u.ml,
             'CBG_pc'       : 34       * u.ug/u.ml,         # This might be the only one here
             'CBG-Cort_pc'  : 1180     * u.ng/u.ml,            
             }


stateVarsMassf = {'CRH'        : stateVarsf['CRHc']         * Paramsf['v_p'],
                 'ACTH'        : stateVarsf['ACTHc']        * Paramsf['v_p'],
                 'Cort_p'      : stateVarsf['Cort_pc']      * Paramsf['v_p'],
                 'Cort_i'      : stateVarsf['Cort_ic']      * Paramsf['v_i'],
                 'Cort_k'      : stateVarsf['Cort_kc']      * Paramsf['v_k'],             
                 'Cortisone_k' : stateVarsf['Cortisone_kc'] * Paramsf['v_k'],
                 'Cortisone_p' : stateVarsf['Cortisone_pc'] * Paramsf['v_p'],          
                 'Alb_tl'      : stateVarsf['Alb_tlc']      * Paramsf['v_tl'], 
                 'Alb-Cort_tl' : stateVarsf['Alb-Cort_tlc'] * Paramsf['v_tl'],
                 'Alb_tm'      : stateVarsf['Alb_tmc']      * Paramsf['v_tm'], 
                 'Alb-Cort_tm' : stateVarsf['Alb-Cort_tmc'] * Paramsf['v_tm'],               
                 'Alb_tu'      : stateVarsf['Alb_tuc']      * Paramsf['v_tu'],
                 'Alb-Cort_tu' : stateVarsf['Alb-Cort_tuc'] * Paramsf['v_tu'],
                 'COPEffect'   : stateVarsf['COPEffect'],
                 'Alb_p'       : stateVarsf['Alb_pc']       * Paramsf['v_p'],
                 'Alb-Cort_p'  : stateVarsf['Alb-Cort_pc']  * Paramsf['v_p'],
                 'CBG_tl'      : stateVarsf['CBG_tlc']      * Paramsf['v_tl'],
                 'CBG-Cort_tl' : stateVarsf['CBG-Cort_tlc'] * Paramsf['v_tl'],
                 'CBG_tm'      : stateVarsf['CBG_tmc']      * Paramsf['v_tm'],
                 'CBG-Cort_tm' : stateVarsf['CBG-Cort_tmc'] * Paramsf['v_tm'], 
                 'CBG_tu'      : stateVarsf['CBG_tuc']      * Paramsf['v_tu'],
                 'CBG-Cort_tu' : stateVarsf['CBG-Cort_tuc'] * Paramsf['v_tu'],
                 'CBG_p'       : stateVarsf['CBG_pc']       * Paramsf['v_p'],
                 'CBG-Cort_p'  : stateVarsf['CBG-Cort_pc']  * Paramsf['v_p'], 
                 }


##########################################################################
# Set up unitless dictionaries                                           #
##########################################################################

# for males
Paramsu = {key : (val.magnitude if isinstance(val, pint.Quantity)
        else val) for key, val in Params.items()}

stateVarsu = {key: val.magnitude if isinstance(val, pint.Quantity)
        else val for key, val in stateVars.items()}

stateVarsMassu = {key: val.magnitude if isinstance(val, pint.Quantity)
        else val for key, val in stateVarsMass.items()}

# for females
Paramsuf = {key : (val.magnitude if isinstance(val, pint.Quantity)
        else val) for key, val in Paramsf.items()}

stateVarsuf = {key: val.magnitude if isinstance(val, pint.Quantity)
        else val for key, val in stateVarsf.items()}

stateVarsMassuf = {key: val.magnitude if isinstance(val, pint.Quantity)
        else val for key, val in stateVarsMassf.items()}


# Convert units for male
Paramsu['konA'] = Params['konA'].to(1 / ((u.nmol/u.ml) * u.min)).magnitude
Paramsu['koffA'] = Params['koffA'].to(1 / u.min).magnitude
Paramsu['konC'] = Params['konC'].to(1 / ((u.nmol/u.ml) * u.min)).magnitude
Paramsu['koffA'] = Params['koffA'].to(1 / u.min).magnitude
Paramsu['Vm2'] = Params['Vm2'].to(u.nmol/(u.ml *u.min)).magnitude
Paramsu['km2'] = Params['km2'].to(u.nmol/u.ml).magnitude
Paramsu['Ald'] = Params['Ald'].to(u.ng/u.ml).magnitude
Paramsu['Vm1'] = Params['Vm1'].to(u.nmol/(u.ml *u.min)).magnitude
Paramsu['km1'] = Params['km1'].to(u.nmol/u.ml).magnitude
Paramsu['Kd_mr'] = Params['Kd_mr'].to(u.nmol/u.ml).magnitude


# Convert some units for feamle
Paramsuf['konA'] = Paramsf['konA'].to(1 / ((u.nmol/u.ml) * u.min)).magnitude
Paramsuf['koffA'] = Paramsf['koffA'].to(1 / u.min).magnitude
Paramsuf['konC'] = Paramsf['konC'].to(1 / ((u.nmol/u.ml) * u.min)).magnitude
Paramsuf['koffA'] = Paramsf['koffA'].to(1 / u.min).magnitude
Paramsuf['Vm2'] = Paramsf['Vm2'].to(u.nmol/(u.ml *u.min)).magnitude
Paramsuf['km2'] = Paramsf['km2'].to(u.nmol/u.ml).magnitude
Paramsuf['Ald'] = Paramsf['Ald'].to(u.ng/u.ml).magnitude
Paramsuf['Vm1'] = Paramsf['Vm1'].to(u.nmol/(u.ml *u.min)).magnitude
Paramsuf['km1'] = Paramsf['km1'].to(u.nmol/u.ml).magnitude
Paramsuf['Kd_mr'] = Paramsf['Kd_mr'].to(u.nmol/u.ml).magnitude



# Return string
ret_str = '\n\treturn rhs'


# Variable order 
var_order = ['CRH', 'ACTH', 'Cort_p', 'Cort_i', 'Cort_k',
             'Cortisone_p', 'Cortisone_k', 'Alb_tl', 'Alb-Cort_tl',
             'Alb_tm', 'Alb-Cort_tm', 'Alb_tu','Alb-Cort_tu',
             'COPEffect', 'Alb_p', 'Alb-Cort_p', 'CBG_tl',
             'CBG-Cort_tl', 'CBG_tm', 'CBG-Cort_tm', 'CBG_tu',
             'CBG-Cort_tu', 'CBG_p', 'CBG-Cort_p']


def createMaster(ps, svs, svms):
    
    # Create Full Data dictionary
    # Recreate the Conversions
    ps['tonano']    = ps['tonano'] * 1e9
    ps['mictonano'] = ps['mictonano'] * 1e3
    ps['konA']      = ps['konA'].to(1 / ((u.nmol/u.ml) * u.min))
    ps['koffA']     = ps['koffA'].to(1 / u.min)
    ps['konC']      = ps['konC'].to(1 / ((u.nmol/u.ml) * u.min))
    ps['koffA']     = ps['koffA'].to(1 / u.min)
    ps['Vm2']       = ps['Vm2'].to(u.nmol/(u.ml *u.min))
    ps['km2']       = ps['km2'].to(u.nmol/u.ml)
    ps['Ald']       = ps['Ald'].to(u.ng/u.ml)
    ps['Vm1']       = ps['Vm1'].to(u.nmol/(u.ml *u.min))
    ps['km1']       = ps['km1'].to(u.nmol/u.ml)
    ps['Kd_mr']     = ps['Kd_mr'].to(u.nmol/u.ml)
    
    
    Master = {}

    def add_sub(M, Sub, Type):
        for key, val in Sub.items():
            sub_dict = {}
            sub_dict['Type'] = Type
            if isinstance(val, pint.Quantity):
                sub_dict['Magnitude'] = val.magnitude
                sub_dict['Unit']      = str(val.units)
            else:
                sub_dict['Magnitude'] = val
                sub_dict['Unit']      = 'Dimensionless'
    
            M[key] = sub_dict
        return M

    add_sub(Master, ps, 'Parameter')
    add_sub(Master, svs, 'State Concentration')
    add_sub(Master, stateVarsMass, 'State Mass')

    return Master
