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

"""
HPA first order fitting.
Variables are to be solved. Not necessarily variables in the ODEs

"""

class DeEquation:
    def __init__(self, eqStr, var_parm_dict = {}):
        self.eqStr = eqStr
        self.vp_dict = var_parm_dict


    def make_fun(self, deinst, funcName, funcHeader,
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
        #print(execStr)
        exec(deinst + '.' + funcName +
             '= types.MethodType(' + funcName + ', ' + deinst
             +')')

u = pint.UnitRegistry()

##########################################################################
# Parameter dictionary                                                   #
##########################################################################





Params = {# CRH
          'k1'       : 0.859     * u.pg / (u.ml * u.min),  
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
          # Cortisol
          'k3'       : 0.001721  * 1 / u.min, 
          'f_lymp_w' : 28.65     * (u.ml / u.min),
          'v_p'      : 3261.72   * u.ml,
          'v_i'      : 12153.9   * u.ml, 
          'f_cap_w'  : 29.06     * (u.ml / u.min),
          'Dcort_pi' : 10.8         * u.ml / u.min,
          'ucp'      : 0.0016    * (1 / u.min), 
          'konA'     : 2310      * 1 / ((u.mol/u.l) * u.s),
          'koffA'    : 0.693     * (1 / u.s),
          'MWCort'   : 362.46    * u.g / u.mol,
          'MWAlb'    : 66437     * u.g / u.mol,
          'MWAC'     : 66862.46  * u.g / u.mol,
          'konC'     : 3.5e9     * 1 / ((u.mol/u.ml) * u.min),
          'koffC'    : 0.88      * (1 / u.min), 
          'MWCBG'    : 52000     * u.g / u.mol,
          'MWCC'     : 52362.46  * u.g / u.mol,
          # Albumin
          'k4'       : 0.005     * u.g / u.min, 
          'Klp'      : 0.1,
          'CapP'     : 101.47    * u.ml / u.min,
          'sigA'     : 0.8,
          'Deg_a'    : 0.005     * u.g / u.min,
          'GFR'      : 103.66    * u.ml / u.min,
          'AlbPerm'  : 0, 
          'ModPerm'  : 0,
          'TL_LW_R'  : 7.908     * u.ml / u.min,
          'TM_LW_R'  : 14.41     * u.ml / u.min,
          'TU_LW_R'  : 6.24      * u.ml / u.min,
          'TL_IW_V'  : 3806.42   * u.ml,
          'TM_IW_V'  : 5790.4    * u.ml,
          'TU_IW_V'  : 2554.52   * u.ml,
          'LAPerm'   : 5.4,
          'MAPerm'   : 6.4,
          'UAPerm'   : 5.4, 
          'Temp_K'   : 309.772,
          'IgpCOP'   : 4.49,
          'HetapCOP' : 0,
          'AGPpCOP'  : 0.9426,
          # CBG
          'k5'       : 0.00423   * u.ug / u.min,
          'sigC'     : 0.8,                               
          'half'     : 9.62e-05  * 1 / u.min, 
          'Dcort_ic' : 15        * u.ml / u.min, # need to add cell surface
          'v_k'      : 5.024     * u.ml, 
          'v_o'      : 25954.97  * u.ml, 
          'Vm2'      : 4        * u.nmol / (u.l * u.min), # no idea here
          'km2'      : 0.04      * u.umol / u.l,         
          'Kd_mr'    : 0.75      * u.nmol / u.l, 
          'Ald'      : 12.1      * u.ng / u.dl,            
          'Vm1'      : 48        * u.nmol / (u.l * u.min), # match / play around a bit
          'km1'      : 21        * u.umol / u.l,           
          'mcort'    : 0.5     * 1 / u.min,              # guesses
          'usp'      : 0.0016    * 1 / u.min,              # guesses
          'msone'    : 0.005     * 1 / u.min,              # guesses
          'tonano'   : 1,
          'mictonano': 1,
          'MWAld'    : 360.44    * u.g / u.mol,
          'MWsone'   : 360.45    * u.g / u.mol
          }


stateVars = {'CRH'         : 7.66 * u.pg / u.ml,               
             'ACTH'        : 21 * u.pg / u.ml,                     
             'Cort_p'      : 3.06 * u.ng / u.ml,         
             'Alb_p'       : 0.046 * u.g / u.ml,       
             'Alb-Cort_p'  : 1301.5 * u.ng / u.ml,        
             'CBG_p'       : 35 * u.ug / u.ml,         
             'CBG-Cort_p'  : 1180 * u.ng / u.ml,     
             'Cort_i'      : 3.06 * u.ng / u.ml,               
             'Alb_i'       : 0.0404 * u.g / u.ml,     
             'Alb-Cort_i'  : 1143 * u.ng / u.ml,          
             'CBG_i'       : 30 * u.ug / u.ml,                             
             'CBG-Cort_i'  : 1000 * u.ng / u.ml,             
             'Cort_k'      : 0.5 * u.ng / u.ml,             
             'Cort_o'      : 2 * u.ng / u.ml,               
             'Cortisone_k' : 0.5 * u.ng / u.ml,              
             'Cortisone_i' : 1 * u.ng / u.ml,          
             'Cortisone_p' : 1 * u.ng / u.ml,          
             'Cortisone_o' : 0.5 * u.ng / u.ml,         
            }

sVars = {"v\['CRH'\]": 'x[0]', "v\['ACTH'\]": 'x[1]',
        "v\['Cort_p'\]": 'x[2]', "v\['Alb_tl'\]": 'x[3]',
        "v\['Ald_tm'\]": 'x[4]', "v\['Alb_tu'\]": 'x[5]',
        "v\['Ald_p'\]": 'x[6]',
        "v\['Alb-Cort_p'\]": 'x[7]', "v\['CBG_p'\]": 'x[8]',
        "v\['CBG-Cort_p'\]": 'x[9]', "v\['Cort_i'\]": 'x[10]',
        "v\['Alb_i'\]": 'x[11]',
        "v\['Alb-Cort_i'\]": 'x[12]', "v\['CBG_i'\]": 'x[13]',
        "v\['CBG-Cort_i'\]": 'x[14]', "v\['Cort_k'\]": 'x[15]',
        "v\['Cort_o'\]": 'x[16]',
        "v\['Cortisone_k'\]": 'x[17]', "v\['Cortisone_i'\]": 'x[18]',
        "v\['Cortisone_p'\]": 'x[19]', "v\['Cortisone_o'\]": 'x[20]'}


testVars = {"v\['Alb_tl'\]": 'x[0]', "v\['Alb_tm'\]": 'x[1]', 
        "v\['Alb_tu'\]": 'x[2]', "v\['Alb_p'\]": 'x[3]',
        "v\['COPEffect'\]": 'x[4]'}


testVars2 = {"v\['CBG_p'\]": 'x[0]', "v\['CBG_i'\]": 'x[1]'}
        
# Create the unitless Parameters and state variables
Paramsu = {key : (val.magnitude if isinstance(val, pint.Quantity)
        else val) for key, val in Params.items()}

stateVarsu = {key: val.magnitude if isinstance(val, pint.Quantity)
        else val for key, val in stateVars.items()}


# Convert some units for numerical constants

Paramsu['tonano'] = Params['tonano'] * 1e9
Paramsu['mictonano'] = Params['mictonano'] * 1e3

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

##########################################################################
# CRF and ACTH functions                                                 #
##########################################################################

# Numerical is good
CRF_head = "def crf(self, t, x, v):\n\trhs = "
CRF_body = ("v['k1']  * "  
           # Positive hill feedback by cort_p
           "(1 + v['zeta'] * ((v['Cort_p'] ** v['alpha'])"
           " / (v['Cort_p'] ** v['alpha'] + v['c'] ** v['alpha']))" 
           # Negative hill feedback by cort_p
           " - (v['psi'] * ((v['Cort_p'] ** v['gamma'])"
           "  / (v['Cort_p'] ** v['gamma'] + v['c3'] ** v['gamma']))))"
           # Generic waste
           "- v['w1'] * v['CRH']"
           )
CRF_ret = '\n\treturn rhs'


CRH = DeEquation(CRF_body)
CRH.make_fun('CRH', 'crf', CRF_head, CRF_ret, sVars)

# Numerical is good
ACTH_head = "def acth(self, t, x, v):\n\trhs = "
ACTH_body = ("v['k2'] * "
           # Negative hill feedback by cort_p
           "(1 - v['rho'] * (v['Cort_p'] ** v['alpha'])"
           " /(v['Cort_p'] ** v['alpha'] + v['c'] ** v['alpha']))"
           # Saturable CRH Stimulation
           #"* ((v['CRH'] ** v['omega'])"
           #"/ (v['CRH'] ** v['omega'] + v['cc'] ** v['omega']))"
           "* v['CRH'] "
           # Generic waste
           "- v['w2'] * v['ACTH']"
           )

ACTH = DeEquation(ACTH_body)
ACTH.make_fun('ACTH', 'acth', ACTH_head, CRF_ret, sVars)


# Add the basic c to c3 binding relation    
Ctoc3_head = "def ctoc3(self, t, x, v):\n\trhs = "
Ctoc3_body = "((v['c'] ** v['alpha'])/10) - (v['c3'] ** v['gamma'])"

Ctoc3 = DeEquation(Ctoc3_body)
Ctoc3.make_fun('Ctoc3', 'ctoc3', Ctoc3_head, CRF_ret, sVars)

##########################################################################
# Plasma Cortisol functions                                              #
##########################################################################

# k3 looks weird in the units cause it doesnt absorb pg to ng of stim
# by acth

# Numerical is good
Cortp_head = "def cort_p(self, t, x, v):\n\trhs = "
Cortp_body = ("v['k3'] * v['ACTH']"   # stimulation
            "+ ((v['f_lymp_w'] / v['v_p']) * v['Cort_i'])"  # lymph
            "- (v['f_cap_w'] / v['v_p']) * v['Cort_p']"  # capillary
            # diffusion
            "- (v['Dcort_pi']/v['v_p']) * (v['Cort_p'] - v['Cort_i'])"
            "-  v['ucp'] * v['Cort_p']"  # Urine
            # Albumin Binding
            #"- v['konA'] *(1/v['MWAlb']) * v['Alb_p'] * v['Cort_p']"
            #" * v['tonano']"
            # albumin off binding
            #"+ v['koffA'] * v['Alb-Cort_p'] * (v['MWCort']/v['MWAC'])" 
            # CBG binding
            #"- v['konC'] * (1/v['MWCBG']) * v['CBG_p'] * v['Cort_p']"
            #" * v['mictonano']"
            # CBG off binding
            #"+ v['koffC'] * v['CBG-Cort_p'] * (v['MWCort']/v['MWCBG'])"  
            )

Cort_p = DeEquation(Cortp_body)
Cort_p.make_fun('Cort_p', 'cort_p', Cortp_head, CRF_ret, sVars)


##########################################################################
# Albumin                                                                #
##########################################################################

# Numerical is good
Albp_head = "def alb_p(self, t, x, v):\n\trhs = "
Albp_body = (#"(1/v['v_p'] * (v['k4'] "  # stimulation
            # lymph
            "+ v['f_lymp_w'] * (3* v['Klp'] * v['Alb_p']"  
            " + (1 - v['Klp']) * v['Alb_i'])"  
            # cap
            #"- v['CapP'] * v['sigA'] * (v['Alb_p'] - v['Alb_i'])"  
            #"- v['Deg_a']"  # Degredation maybe make it linear
            #"- v['GFR'] * v['gfrPa'] * v['Alb_p'])) "  # Urine
            # albumin binding
            #"- v['konA'] * (1/v['MWCort']) * v['Alb_p']"  
            #"* v['Cort_p'] " 
            # albumin off binding
            #"+ v['koffA'] * v['Alb-Cort_p'] * (v['MWAlb']/v['MWAC'])"
            #" * (1/v['tonano'])"
            )

Alb_p = DeEquation(Albp_body)
Alb_p.make_fun('Alb_p', 'alb_p', Albp_head, CRF_ret, sVars)

# numerical good, binding is a little weird
Albcortp_head = "def albcort_p(self, t, x, v):\n\trhs = "
Albcortp_body = ("(1/v['v_p'])"
                 # Lymph
                 " * (v['f_lymp_w'] * (3 *v['Klp'] * v['Alb-Cort_p']" 
                 " + (1 - v['Klp']) * v['Alb-Cort_i'])" 
                 # cap
                 "-v['CapP'] * v['sigA'] * (v['Alb-Cort_p'] " 
                 "-v['Alb-Cort_i'])" 
                 "- v['GFR'] *  v['gfrPa'] * v['Alb-Cort_p']) " # Urine
                 # albumin binding
                 #"+ v['konA'] * (1/v['MWAlb']) * v['Alb_p']"  
                 #"* v['Cort_p'] * (1/v['MWCort']) * (v['MWAC']) "
                 #"* v['tonano']"
                 # off binding
                 #"- v['koffA'] * v['Alb-Cort_p']"
                 )

AlbCort_p = DeEquation(Albcortp_body)
AlbCort_p.make_fun('AlbCort_p', 'albcort_p', Albcortp_head,CRF_ret, sVars)


Albtoconc_head = "def albtoconc(self, x, v):\n\trhs = "
Albtoconc_body = (
                 "[v['Alb_tl'] / v['TL_IW_V'], "
                 " v['Alb_tm'] / v['TM_IW_V'], "
                 " v['Alb_tu'] / v['TU_IW_V'], "
                 " v['Alb_p'] / v['v_p']]"
                 )

Albtoconc = DeEquation(Albtoconc_body)
Albtoconc.make_fun('Albtoconc', 'albtoconc', Albtoconc_head, 
        CRF_ret, testVars)



# Lower Lymph Alb
TL_LymphAlbRate_head = "def tl_lymphalbrate(self, t, x, v):\n\trhs = "
TL_LymphAlbRate_body = (
                        "((v['Klp'] * v['Alb_p'])"
                        "+ ( (1 - v['Klp']) * v['Alb_tl'] ))"
                        " * v['TL_LW_R']"
                       )

TL_LymphAlbRate = DeEquation(TL_LymphAlbRate_body)
TL_LymphAlbRate.make_fun('TL_LymphAlbRate', 'tl_lymphalbrate', 
        TL_LymphAlbRate_head, CRF_ret, testVars)


# Middle Lymph Alb
TM_LymphAlbRate_head = "def tm_lymphalbrate(self, t, x, v):\n\trhs = "
TM_LymphAlbRate_body = (
                        "((v['Klp'] * v['Alb_p'])"
                        "+ ( (1 - v['Klp']) * v['Alb_tm'] ))"
                        " * v['TM_LW_R']"
                       )

TM_LymphAlbRate = DeEquation(TM_LymphAlbRate_body)
TM_LymphAlbRate.make_fun('TM_LymphAlbRate', 'tm_lymphalbrate', 
        TM_LymphAlbRate_head, CRF_ret, testVars)


# Upper Lymph Alb
TU_LymphAlbRate_head = "def tu_lymphalbrate(self, t, x, v):\n\trhs = "
TU_LymphAlbRate_body = (
                        "((v['Klp'] * v['Alb_p'])"
                        "+ ( (1 - v['Klp']) * v['Alb_tu'] ))"
                        " * v['TU_LW_R']"
                       )

TU_LymphAlbRate = DeEquation(TU_LymphAlbRate_body)
TU_LymphAlbRate.make_fun('TU_LymphAlbRate', 'tu_lymphalbrate', 
        TU_LymphAlbRate_head, CRF_ret, testVars)


# Total Lymph Alb
Total_LymphAlbRate_head = "def total_lymphalbrate(self, t, x, v):\n\trhs = "
Total_LymphAlbRate_body = (
                           "TL_LymphAlbRate.tl_lymphalbrate(t, x, v)"
                           " + TM_LymphAlbRate.tm_lymphalbrate(t, x, v)"
                           " + TU_LymphAlbRate.tu_lymphalbrate(t, x, v)"
                          )


Total_LymphAlbRate = DeEquation(Total_LymphAlbRate_body)
Total_LymphAlbRate.make_fun('Total_LymphAlbRate', 'total_lymphalbrate', 
        Total_LymphAlbRate_head, CRF_ret, testVars)


# Lower Cap Alb
TL_CapAlbRate_head = "def tl_capalbrate(self, t, x, v):\n\trhs = "
TL_CapAlbRate_body = (
                      "v['LAPerm'] * v['sigA'] "
                      " * (v['Alb_p'] - v['Alb_tl'])"
                       )

TL_CapAlbRate = DeEquation(TL_CapAlbRate_body)
TL_CapAlbRate.make_fun('TL_CapAlbRate', 'tl_capalbrate', 
        TL_CapAlbRate_head, CRF_ret, testVars)


# Middle Cap Alb
TM_CapAlbRate_head = "def tm_capalbrate(self, t, x, v):\n\trhs = "
TM_CapAlbRate_body = (
                      "v['MAPerm'] * v['sigA'] " 
                      "* (v['Alb_p'] - v['Alb_tm'])"
                     )

TM_CapAlbRate = DeEquation(TM_CapAlbRate_body)
TM_CapAlbRate.make_fun('TM_CapAlbRate', 'tm_capalbrate', 
        TM_CapAlbRate_head, CRF_ret, testVars)


# Upper Cap Alb
TU_CapAlbRate_head = "def tu_capalbrate(self, t, x, v):\n\trhs = "
TU_CapAlbRate_body = (
                      "v['UAPerm'] * v['sigA']"
                      "* (v['Alb_p'] - v['Alb_tu'])"
                       )

TU_CapAlbRate = DeEquation(TU_CapAlbRate_body)
TU_CapAlbRate.make_fun('TU_CapAlbRate', 'tu_capalbrate', 
        TU_CapAlbRate_head, CRF_ret, testVars)


# Total Cap Alb
Total_CapAlbRate_head = "def total_capalbrate(self, t, x, v):\n\trhs = "
Total_CapAlbRate_body = (
                         "TL_CapAlbRate.tl_capalbrate(t, x, v)"
                         " + TM_CapAlbRate.tm_capalbrate(t, x, v)"
                         " + TU_CapAlbRate.tu_capalbrate(t, x, v)"
                        )


Total_CapAlbRate = DeEquation(Total_CapAlbRate_body)
Total_CapAlbRate.make_fun('Total_CapAlbRate', 'total_capalbrate', 
        Total_CapAlbRate_head, CRF_ret, testVars)


# Filtration
CD_AlbOutflow_head = "def cd_alboutflow(self, t, x, v):\n\trhs = "
CD_AlbOutflow_body = (
                      "v['Alb_p'] * (v['ModPerm'] +v['AlbPerm'])"
                      "* v['GFR'] / 1000"
                     )

CD_AlbOutflow = DeEquation(CD_AlbOutflow_body)
CD_AlbOutflow.make_fun('CD_AlbOutflow', 'cd_alboutflow', 
        CD_AlbOutflow_head, CRF_ret, testVars)

# Albumin partial colliod osmotic
Albpcop_head = "def albpcop(self, x, v):\n\trhs = "
Albpcop_body = (
                 "2.1 * (v['Alb_p'] * 1000 / v['MWAlb']) * 62.36367"
                 "* v['Temp_K'] * v['sigA']"
                 )

Albpcop = DeEquation(Albpcop_body)
Albpcop.make_fun('Albpcop', 'albpcop', Albpcop_head, 
        CRF_ret, testVars)


# colloid osmotic pressure
COP_head = "def cop(self, x, v):\n\trhs = "
COP_body = (
            "Albpcop.albpcop(x, v) + v['IgpCOP'] + v['HetapCOP']"
            " + v['AGPpCOP']"
            )

COP = DeEquation(COP_body)
COP.make_fun('COP', 'cop', COP_head, 
        CRF_ret, testVars)


# COPeffect 
COPE_head = "def cope(self, x, v):\n\trhs = "
COPE_body = (
             "0.015 * (copeCurve(COP.cop(x,v)) - v['COPEffect'])"
            )

COPE = DeEquation(COPE_body)
COPE.make_fun('COPE', 'cope', COPE_head, 
        CRF_ret, testVars)

# set up the cop effect curve
copefCurve = cs([20, 28, 40],[1.15, 1, 0],[0, -0.07, 0])

def copeCurve(cop):
    if cop <= 20.0:
        return 1.15
    elif cop >= 40.0:
        return 0
    else:
        return copefCurve.__call__(cop)



# Torso Lower Interstitial Alb
Albtl_head = "def alb_tl(self, t, x, v):\n\trhs = "
Albtl_body = (
              "TL_CapAlbRate.tl_capalbrate(t, x, v)"
              "- TL_LymphAlbRate.tl_lymphalbrate(t, x, v)"
             )

Albtl = DeEquation(Albtl_body)
Albtl.make_fun('Albtl', 'alb_tl', 
        Albtl_head, CRF_ret, testVars)


# Torso Middle Interstitial Protein Alb
Albtm_head = "def alb_tm(self, t, x, v):\n\trhs = "
Albtm_body = (
              "TM_CapAlbRate.tm_capalbrate(t, x, v)"
              "- TM_LymphAlbRate.tm_lymphalbrate(t, x, v)"
             )

Albtm = DeEquation(Albtm_body)
Albtm.make_fun('Albtm', 'alb_tm', 
        Albtm_head, CRF_ret, testVars)


# Torso_Upper_InterstitialProtein.[Alb]
Albtu_head = "def alb_tu(self, t, x, v):\n\trhs = "
Albtu_body = (
              "TU_CapAlbRate.tu_capalbrate(t, x, v)"
              "- TU_LymphAlbRate.tu_lymphalbrate(t, x, v)"
             )

Albtu = DeEquation(Albtu_body)
Albtu.make_fun('Albtu', 'alb_tu', 
        Albtu_head, CRF_ret, testVars)


# Plasma 
Albp_head = "def alb_p(self, t, x, v):\n\trhs = "
Albp_body = (
             "v['k4'] * v['COPEffect']" 
             " + Total_LymphAlbRate.total_lymphalbrate(t, x, v)"
             "- (v['Deg_a'] "
             "+ Total_CapAlbRate.total_capalbrate(t, x, v)" 
             " + CD_AlbOutflow.cd_alboutflow(t, x, v))"
            )

Albp = DeEquation(Albp_body)
Albp.make_fun('Albp', 'alb_p', 
        Albp_head, CRF_ret, testVars)





##########################################################################
# plasma CBG and CBG bound Cortisol functions                            #
##########################################################################

# Numerical done
CBGp_head = "def cbg_p(self, t, x, v):\n\trhs = "
CBGp_body = ("(1/v['v_p']) * (v['k5'] "  # stimulation
            # lymph
            "+ v['f_lymp_w'] * (3 * v['Klp'] * v['CBG_p']"  
            " + (1 - v['Klp']) * v['CBG_i'])"  
            # cap
            "- v['CapP'] * v['sigC'] * (v['CBG_p'] - v['CBG_i']))"  
            "- v['half'] * v['CBG_p'] "  # half life rate
            # CBG binding
            #"- v['konC'] * (1/v['MWCort']) * v['CBG_p']"  
            #"* v['Cort_p']"
            # Off
            #"+ v['koffC'] * v['CBG-Cort_p'] * (1/v['mictonano'])"
            )  # off binding

CBG_p = DeEquation(CBGp_body)
CBG_p.make_fun('CBG_p', 'cbg_p', CBGp_head, CRF_ret, testVars2)


# numerical done
CBGCortp_head = "def cbgcort_p(self, t, x, v):\n\trhs = "
CBGCortp_body = ("(1/v['v_p'])"
                 # Lymph
                 " * (v['f_lymp_w'] * (3 *v['Klp'] * v['CBG-Cort_p']" 
                 " + (1 - v['Klp']) * v['CBG-Cort_i'])" 
                 # cap
                 "-v['CapP'] * v['sigC'] * (v['CBG-Cort_p'] " 
                 "-v['CBG-Cort_i']))" 
                 # binding
                 #"+ v['konC'] * (1/v['MWCBG']) * v['CBG_p']"  
                 #"* v['Cort_p'] * (1/v['MWCort']) *v['MWCC'] "
                 #" * v['mictonano']"
                 #"- v['koffC'] * v['CBG-Cort_p'] " #  off binding
                 )

CBGCort_p = DeEquation(CBGCortp_body)
CBGCort_p.make_fun('CBGCort_p', 'cbgcort_p', CBGCortp_head, CRF_ret,
        sVars)


##########################################################################
# Interstitial Cortisol functions                                        #
##########################################################################

# numerical done but derivative is very low 
Corti_head = "def cort_i(self, t, x, v):\n\trhs = "
Corti_body = ("- (v['f_lymp_w'] / v['v_i']) * v['Cort_i'] " # lymph
            #"+ (v['f_cap_w'] / v['v_i']) * v['Cort_p'] " # capillary
            # diffusion
            #"+ (v['Dcort_pi']/v['v_i']) * (v['Cort_p'] - v['Cort_i']) "
            # diffusion
            #"- ( 100* v['Dcort_ic']/v['v_i']) * (v['Cort_i'] - v['Cort_k']) "
            # diffusion
            #"- (v['Dcort_ic']/v['v_i']) * (v['Cort_i'] - v['Cort_o']) "
            # Albumin Binding
            #"- v['konA'] *(1/v['MWAlb']) * v['Alb_i'] * v['Cort_i']"
            #" * v['tonano']"
            # albumin off binding
            #"+ v['koffA'] * v['Alb-Cort_i'] * (v['MWCort']/v['MWAC'])" 
            # CBG binding
            #"- v['konC'] * (1/v['MWCBG']) * v['CBG_i'] * v['Cort_i']"
            #" * v['mictonano']"
            # CBG off binding
            #"+ v['koffC'] * v['CBG-Cort_p'] * (v['MWCort']/v['MWCBG'])"  
            )


Cort_i = DeEquation(Corti_body)
Cort_i.make_fun('Cort_i', 'cort_i', Corti_head, CRF_ret,
        sVars)


##########################################################################
# Interstitial Albumin and Albumin bound Cortisol functions              #
##########################################################################

# Numerical done
Albi_head = "def alb_i(self, t, x, v):\n\trhs = "
Albi_body = ("1/v['v_i'] * "
            # lymph
            "(- v['f_lymp_w'] * (3 * v['Klp'] * v['Alb_p']"  
            " + (1 - v['Klp']) * v['Alb_i'])"  
            # cap
            "+ v['CapP'] * v['sigA'] * (v['Alb_p'] - v['Alb_i']))"  
            # albumin binding
            #"- v['konA'] * (1/v['MWCort']) * v['Alb_i']"  
            #"* v['Cort_i'] " 
            # albumin off binding
            #"+ v['koffA'] * v['Alb-Cort_i'] * (v['MWAlb']/v['MWAC'])"
            #" * (1/v['tonano'])"
            )


Alb_i = DeEquation(Albi_body)
Alb_i.make_fun('Alb_i', 'alb_i', Albi_head, CRF_ret, sVars)

# Numerical Done but derive is big
Albcorti_head = "def albcort_i(self, t, x, v):\n\trhs = "
Albcorti_body = ("(1/v['v_i'])"
                 # Lymph
                 " * (- v['f_lymp_w'] * (3 *v['Klp'] * v['Alb-Cort_p']" 
                 " + (1 - v['Klp']) * v['Alb-Cort_i'])" 
                 # cap
                 "+ v['CapP'] * v['sigA'] * (v['Alb-Cort_p'] " 
                 "-v['Alb-Cort_i']))" 
                 # albumin binding
                 #"+ v['konA'] * (1/v['MWAlb']) * v['Alb_i']"  
                 #"* v['Cort_i'] * (1/v['MWCort']) * (v['MWAC']) "
                 #"* v['tonano']"
                 # off binding
                 #"- v['koffA'] * v['Alb-Cort_i']"
                 )  


AlbCort_i = DeEquation(Albcorti_body)
AlbCort_i.make_fun('AlbCort_i', 'albcort_i', Albcorti_head, CRF_ret, sVars)





##########################################################################
# CBG and CBG bound Cortisol functions                                   #
##########################################################################

# Numerical Done
CBGi_head = "def cbg_i(self, t, x, v):\n\trhs = "
CBGi_body = ("(1/v['v_i']) *"  # stimulation
            # lymph
            "(- v['f_lymp_w'] * (3 * v['Klp'] * v['CBG_p']"  
            " + (1 - v['Klp']) * v['CBG_i'])"  
            # cap
            "+ v['CapP'] * v['sigC'] * (v['CBG_p'] - v['CBG_i']))"  
            # CBG binding
            #"- v['konC'] * (1/v['MWCort']) * v['CBG_i']"  
            #"* v['Cort_i']"
            # Off
            #"+ v['koffC'] * v['CBG-Cort_i'] * (1/v['mictonano'])"
            )  # off binding


CBG_i = DeEquation(CBGi_body)
CBG_i.make_fun('CBG_i', 'cbg_i', CBGi_head, CRF_ret, testVars2)



CBGCorti_head = "def cbgcort_i(self, t, x, v):\n\trhs = "
CBGCorti_body = ("(1/v['v_i'])"
                 # Lymph
                 " * (- v['f_lymp_w'] * (3 *v['Klp'] * v['CBG-Cort_p']" 
                 " + (1 - v['Klp']) * v['CBG-Cort_i'])" 
                 # cap
                 "+ v['CapP'] * v['sigC'] * (v['CBG-Cort_p'] " 
                 "- v['CBG-Cort_i']))" 
                 # albumin binding
                 #"+ v['konC'] * (1/v['MWCort']) * v['CBG_i']"  
                 #"* v['Cort_i'] * v['mictonano']"
                 #"- v['koffC'] * v['CBG-Cort_i'] " #  off binding
                 )

CBGCort_i = DeEquation(CBGCorti_body)
CBGCort_i.make_fun('CBGCort_i', 'cbgcort_i', CBGCorti_head, CRF_ret,
        sVars)


##########################################################################
# Cortisol in kidney and other cells                                     #
##########################################################################


MR_head = "def mr_b(self, t, x, v):\n\trhs = "
MR_body = ("v['Cort_k'] * (1/v['MWCort']) "
           "/ (v['Kd_mr'] + v['Cort_k'] * (1/v['MWCort']) "
           "+ v['Ald'] * (1/v['MWAld']))"
          )

MR_B = DeEquation(MR_body)
MR_B.make_fun('MR_B', 'mr_b', MR_head, CRF_ret,
        sVars)

MRA_head = "def mra_b(self, t, x, v):\n\trhs = "
MRA_body = ("v['Ald'] * (1/v['MWAld']) "
           "/ (v['Kd_mr'] + v['Cort_k'] * (1/v['MWCort']) "
           "+ v['Ald'] * (1/v['MWAld']))"
          )

MR_A = DeEquation(MRA_body)
MR_A.make_fun('MR_A', 'mra_b', MRA_head, CRF_ret,
        sVars)


# Numerical good
Cortk_head = "def cort_k(self, t, x, v):\n\trhs = "
Cortk_body = (# Diffusion
              #"v['Dcort_ic'] * (v['Cort_i'] - v['Cort_k'])/v['v_k'] "
              # 11B-HSD2
              "- v['Vm2'] * ((1 - MR_B.mr_b(0,x,v)) * v['Cort_k'])" 
              "/((1 - MR_B.mr_b(0, x,v)) * v['Cort_k'] * (1/v['MWCort'])"
              " + v['km2'])"
             )


Cort_k = DeEquation(Cortk_body)
Cort_k.make_fun('Cort_k', 'cort_k', Cortk_head, CRF_ret, sVars)

# Numerical Good
Corto_head = "def cort_o(self, t, x, v):\n\trhs = "
Corto_body = ( # Diffusion
             "1000* v['Dcort_ic'] * (v['Cort_i'] - v['Cort_o'])/v['v_o'] " 
             # 11B-HSD1
             "+ v['Vm1'] * v['MWCort'] "
             "* (v['Cortisone_o'] * (1/v['MWsone']))"
             " / (v['Cortisone_o'] * (1/v['MWsone']) + v['km1'])"
             # General Metabolism
             "- v['mcort'] * v['Cort_o']"
             ) 

Cort_o = DeEquation(Corto_body)
Cort_o.make_fun('Cort_o', 'cort_o', Corto_head, CRF_ret, sVars)



##########################################################################
# Cortisone  in kidney and other cells                                   #
##########################################################################

# Numerical Good
Cortisonek_head = "def cortisone_k(self, t, x, v):\n\trhs = "
Cortisonek_body = ( # Diffusion
                  "(v['Dcort_ic']/v['v_k']) *(v['Cortisone_i'] -"
                  " v['Cortisone_k']) "
                  # 11B-HSD2
                  "+ v['Vm2'] * ((1 - MR_B.mr_b(0,x,v)) * v['Cort_k'])" 
                  "/((1 - MR_B.mr_b(0, x,v)) * v['Cort_k'] " 
                  "* (1/v['MWCort']) + v['km2'])"
                  )

Cortisone_k = DeEquation(Cortisonek_body)
Cortisone_k.make_fun('Cortisone_k', 'cortisone_k',
        Cortisonek_head, CRF_ret, sVars)

# Numerical Good
Cortisonei_head = "def cortisone_i(self, t, x, v):\n\trhs = "
Cortisonei_body = ( # lymph
                   "- (v['f_lymp_w'] / v['v_i']) * v['Cortisone_i'] "
                   # capillary
                   "+ (v['f_cap_w'] / v['v_p']) * v['Cortisone_p']"
                   # Diffusion
                   "+ v['Dcort_pi'] * (v['Cortisone_p'] - "
                   "v['Cortisone_i'])/v['v_i'] "
                   # Diffusion
                   "- v['Dcort_ic'] * (v['Cortisone_i'] - "
                   "v['Cortisone_k'])/v['v_i'] "
                   # Diffusion
                   "- v['Dcort_ic'] * (v['Cortisone_i'] - "
                   "v['Cortisone_o'])/v['v_i'] "
                   )

Cortisone_i = DeEquation(Cortisonei_body)
Cortisone_i.make_fun('Cortisone_i', 'cortisone_i', Cortisonei_head,
                     CRF_ret, sVars)

# Numerical Good
Cortisonep_head = "def cortisone_p(self, t, x, v):\n\trhs = "
Cortisonep_body = (# lymph 
                   "(v['f_lymp_w'] / v['v_i']) * v['Cortisone_i'] "
                   # lymph
                   "- (v['f_cap_w'] / v['v_p']) * v['Cortisone_p'] "
                   # Diffusion
                   "- v['Dcort_pi'] * (v['Cortisone_p'] - "
                   "v['Cortisone_i'])/v['v_p'] "
                   # Urine
                   "- v['usp'] * v['Cortisone_p']"
                   )

Cortisone_p= DeEquation(Cortisonep_body)
Cortisone_p.make_fun('Cortisone_p', 'cortisone_p', Cortisonep_head,
        CRF_ret, sVars)


Cortisoneo_head = "def cortisone_o(self, t, x, v):\n\trhs = "
Cortisoneo_body = (# Diffusion
                   "v['Dcort_ic'] * (v['Cortisone_i'] - "
                   "v['Cortisone_o'])/v['v_o'] " 
                   # 11B-HSD1
                   "+ v['Vm1'] * v['MWCort'] "
                   "* (v['Cortisone_o'] * (1/v['MWsone']))"
                   " / (v['Cortisone_o'] * (1/v['MWsone']) + v['km1'])"
                   # General Metabolism
                   "- v['msone'] * v['Cortisone_o']"     
                   )

Cortisone_o= DeEquation(Cortisoneo_body)
Cortisone_o.make_fun('Cortisone_o', 'cortisone_o', Cortisoneo_head,
        CRF_ret, sVars)




##########################################################################
# Full function                                                          #
##########################################################################


def tester():
    testvec = [val for val in stateVars.values()]
    print('CRH:', CRH.crf(0, testvec, Params))
    # ACTH
    print('ACTH:',ACTH.acth(0, testvec, Params))
    # Cort_p
    print('Cort_p:', Cort_p.cort_p(0, testvec, 
        Params).to(u.ng/(u.ml * u.min)))
    # Alb_p
    print('Alb_p:', Alb_p.alb_p(0, testvec, 
        Params).to(u.g/(u.ml * u.min)))
    # Alb-Cort_p
    print('Alb-Cort_p:', AlbCort_p.albcort_p(0, testvec, 
        Params).to(u.ng/(u.ml * u.min)))
    # CBG_p
    print('CBG_p:', CBG_p.cbg_p(0, testvec, 
        Params).to(u.ug/(u.ml * u.min)))
    # CBG-Cort_p
    print('CBG-Cort_p:', CBGCort_p.cbgcort_p(0, testvec, 
        Params).to(u.ng/(u.ml * u.min)))
    # Cort_i
    print('Cort_i:', Cort_i.cort_i(0, testvec, 
        Params).to(u.ng/(u.ml * u.min)))
    # Alb_i
    print('Alb_i:', Alb_i.alb_i(0, testvec, 
        Params).to(u.g/(u.ml * u.min)))
    # Alb-Cort_i
    print('Alb-Cort_i:', AlbCort_i.albcort_i(0, testvec, 
        Params).to(u.ng/(u.ml * u.min)))
    # CBG_i
    print('CBG_i:', CBG_i.cbg_i(0, testvec, 
        Params).to(u.ug/(u.ml * u.min)))
    # CBG-Cort_i
    print('CBG-Cort_i:', CBGCort_i.cbgcort_i(0, testvec, 
        Params).to(u.ng/(u.ml * u.min)))
    # Cort_k
    print("Cort_k:", Cort_k.cort_k(0, testvec, 
        Params).to(u.ng/(u.ml * u.min)))
    # Cort_o
    print("Cort_o:", Cort_o.cort_o(0, testvec, 
        Params))
    # Cortisone_k
    print("Cortisone_k:", Cortisone_k.cortisone_k(0, testvec, 
        Params).to(u.ng/(u.ml * u.min)))
    # Cortisone_i
    print("Cortisone_i:", Cortisone_i.cortisone_i(0, testvec, Params))
    # Cortisone_p
    print("Cortisone_p:", Cortisone_p.cortisone_p(0, testvec, Params))
    # Cortisone_o
    print("Cortisone_o:", Cortisone_o.cortisone_o(0, testvec, Params))

    return 0


def testeru():
    testvec = [val for val in stateVarsu.values()]
    # CRH
    print('CRH Unitless:', CRH.crf(0, testvec, Paramsu))
    # ACTH
    print('ACTH:',ACTH.acth(0, testvec, Paramsu))
    # Cort_p
    print('Cort_p:', Cort_p.cort_p(0, testvec, Paramsu))
    # Alb_p
    print('Alb_p:', Alb_p.alb_p(0, testvec, Paramsu))
    # Alb-Cort_p
    print('Alb-Cort_p:', AlbCort_p.albcort_p(0, testvec, Paramsu))
    # CBG_p
    print('CBG_p:', CBG_p.cbg_p(0, testvec, Paramsu))
    # CBG-Cort_p
    print('CBG-Cort_p:', CBGCort_p.cbgcort_p(0, testvec, Paramsu))
    # Cort_i
    print('Cort_i:', Cort_i.cort_i(0, testvec, Paramsu))
    # Alb_i
    print('Alb_i:', Alb_i.alb_i(0, testvec, Paramsu))
    # Alb-Cort_i
    print('Alb-Cort_p:', AlbCort_i.albcort_i(0, testvec, Paramsu))
    # CBG_i
    print('CBG_i:', CBG_i.cbg_i(0, testvec, Paramsu))
    # CBG-Cort_i
    print('CBG-Cort_i:', CBGCort_i.cbgcort_i(0, testvec, Paramsu))
    # Cort_k
    print("Cort_k:", Cort_k.cort_k(0, testvec, Paramsu))
    # Cort_o
    print("Cort_o:", Cort_o.cort_o(0, testvec, Paramsu))
    # Cortisone_k
    print("Cortisone_k:", Cortisone_k.cortisone_k(0, testvec, Paramsu))
    # Cortisone_i
    print("Cortisone_i:", Cortisone_i.cortisone_i(0, testvec, Paramsu))
    # Cortisone_p
    print("Cortisone_p:", Cortisone_p.cortisone_p(0, testvec, Paramsu))
    # Cortisone_o
    print("Cortisone_o:", Cortisone_o.cortisone_o(0, testvec, Paramsu))
    
    return 0


def alb_sys(t, w, v):
    # Broader model I will need to unpack
    x  = Albtoconc.albtoconc(w, v)
    x += [w[4]]
    
    ddt = [
           Albtl.alb_tl(t, x, v), Albtm.alb_tm(t, x, v),
           Albtu.alb_tu(t, x, v), Albp.alb_p(t, x, v),
           COPE.cope(x, v)
          ]

    return ddt


def CGB_sys(t, x, v):
    ddt = [CBGp.cbg_p(t, x, v), CBGi.cbg_i(t, x, v)]
    return ddt



def DE(t, x, v):
   
    ddt = [CRH.crf(t, x, v), ACTH.acth(t, x, v), Cort_p.cort_p(t, x, v),
            Alb_p.alb_p(t, x, v), AlbCort_p.albcort_p(t, x, v),
            CBG_p.cbg_p(t, x, v), CBGCort_p.cbgcort_p(t, x, v), 
            Cort_i.cort_i(t, x, v), Alb_i.alb_i(t, x, v), 
            AlbCort_i.albcort_i(t, x, v), CBG_i.cbg_i(t, x, v), 
            CBGCort_i.cbgcort_i(t, x, v), Cort_k.cort_k(t, x, v), 
            Cort_o.cort_o(t, x, v), Cortisone_k.cortisone_k(t, x, v), 
            Cortisone_i.cortisone_i(t, x, v), 
            Cortisone_p.cortisone_p(t, x, v), 
            Cortisone_o.cortisone_o(t, x, v)]

    return ddt

color_cycle = plt.rcParams['axes.prop_cycle']()

def plot_sols(t, sols):
    """ 
    Plot mutliple sols of the hpa system
    """
    fig, ax = plt.subplots(3,3, figsize=(9,7), sharex=True)
    # CRH ACTH
    crh = ax[0,0].plot(t, sols[:, 0], label='CRH', **next(color_cycle)) 
    acth = ax[0,0].plot(t, sols[:, 1], label='ACTH', **next(color_cycle)) 
    ax[0,0].set_ylabel('pg/ml')
    ax[0,0].legend(bbox_to_anchor=(0,1.02,1,0.2), loc="lower left",
                mode="expand", borderaxespad=0, ncol=2)    
    
    # Cort_p cort_i
    cortp = ax[0,1].plot(t, sols[:, 2], label=r'$Cort_p$', **next(color_cycle)) 
    corti = ax[0,1].plot(t, sols[:, 7], label=r'$Cort_i$', **next(color_cycle)) 
    ax[0,1].set_ylabel('ng/ml')
    ax[0,1].legend(bbox_to_anchor=(0,1.02,1,0.2), loc="lower left",
                mode="expand", borderaxespad=0, ncol=2)    

    # Cort_k cort_o
    cortk = ax[0,2].plot(t, sols[:, 12], label=r'$Cort_k$', **next(color_cycle)) 
    corto = ax[0,2].plot(t, sols[:, 13], label=r'$Cort_o$', **next(color_cycle)) 
    ax[0,2].set_ylabel('ng/ml')
    ax[0,2].legend(bbox_to_anchor=(0,1.02,1,0.2), loc="lower left",
                mode="expand", borderaxespad=0, ncol=2)

    
    # alb_p alb_i
    albp = ax[1,0].plot(t, sols[:, 3], label=r'$Alb_p$', **next(color_cycle)) 
    albi = ax[1,0].plot(t, sols[:, 8], label=r'$Alb_i$', **next(color_cycle))
    ax[1,0].set_ylabel('g/ml')
    ax[1,0].legend(bbox_to_anchor=(0,1.02,1,0.2), loc="lower left",
                mode="expand", borderaxespad=0, ncol=2)       
 
    # CBG_p CBG_i
    CBGp = ax[1,1].plot(t, sols[:, 5], label=r'$CBG_p$', **next(color_cycle)) 
    CBGi = ax[1,1].plot(t, sols[:, 10], label=r'$CBG_i$', **next(color_cycle))
    ax[1,1].set_ylabel(r'$\mu$g/ml')
    ax[1,1].legend(bbox_to_anchor=(0,1.02,1,0.2), loc="lower left",
                mode="expand", borderaxespad=0, ncol=2)       

    # albcort_p albcort_i
    acp = ax[1,2].plot(t, sols[:, 4], label=r'$AlbCort_p$', **next(color_cycle)) 
    aci = ax[1,2].plot(t, sols[:, 9], label=r'$AlbCort_i$', **next(color_cycle))
    ax[1,2].set_ylabel('ng/ml')
    ax[1,2].legend(bbox_to_anchor=(0,1.02,1,0.2), loc="lower left",
                mode="expand", borderaxespad=0, ncol=2, prop={'size': 8})       


    # cbgcort_p cbgcort_i
    ccp = ax[2,0].plot(t, sols[:, 6], label=r'$CBGCort_p$', **next(color_cycle)) 
    cci = ax[2,0].plot(t, sols[:, 11], label=r'$CBGCort_i$', **next(color_cycle))
    ax[2,0].set_ylabel('ng/ml')
    ax[2,0].legend(bbox_to_anchor=(0,1.02,1,0.2), loc="lower left",
                mode="expand", borderaxespad=0, ncol=2, prop={'size': 8})       

    # sone k sone o
    sonek = ax[2,1].plot(t, sols[:, 14], label=r'$Csone_k$', **next(color_cycle)) 
    soneo = ax[2,1].plot(t, sols[:, -1], label=r'$Csone_o$', **next(color_cycle))
    ax[2,1].set_ylabel('ng/ml')
    ax[2,1].legend(bbox_to_anchor=(0,1.02,1,0.2), loc="lower left",
                mode="expand", borderaxespad=0, ncol=2, prop={'size': 8})       

    # sone p sone i
    sonep = ax[2,2].plot(t, sols[:, 16], label=r'$Csone_p$', **next(color_cycle)) 
    sonei = ax[2,2].plot(t, sols[:, 15], label=r'$Csone_i$', **next(color_cycle))
    ax[2,2].set_ylabel('ng/ml')
    ax[2,2].legend(bbox_to_anchor=(0,1.02,1,0.2), loc="lower left",
                mode="expand", borderaxespad=0, ncol=2, prop={'size': 8})       

    plt.subplots_adjust(wspace=0.4, hspace=0.4) 
    fig.text(0.5, 0.005, 'Time (min)', ha='center')


    return fig


##########################################################################
# main                                                                   #
##########################################################################

def main():
   
    i0 = [35, 30]

    sol = solve_ivp(CGB_sys, [0, 6000], i0, args=(Paramsu,),
             method='LSODA')

    z = sol.y.T
    
    
    """
    i0 = [53, 53, 42, 150.8, 1]
    #i0c = Albtoconc.albtoconc(i0[:-1], Paramsu)
    #i0c = i0c + [1]

    #print(Albpcop.albpcop(i0c, Paramsu))
    #print(COP.cop(i0c, Paramsu))
    #print(COPE.cope(i0c, Paramsu))
    #print(Albp.alb_p(0, i0c, Paramsu))

    print(alb_sys(0, i0, Paramsu))
    sol = solve_ivp(alb_sys, [0, 6000], i0, args=(Paramsu,),
             method='LSODA')

    z = sol.y.T  
    f = lambda x: Albtoconc.albtoconc(x, Paramsu)
    c = z
    #c = np.apply_along_axis(f, 1, z) 
 
    fig, axes = plt.subplots(2,2, sharex=True) 
 
    axes[0,0].plot(sol.t, c[:,3]) 
    axes[0,0].set_title('Plasma Concentration') 
    axes[0,1].plot(sol.t, c[:,0]) 
    axes[0,1].set_title('Lower Tissue') 
    axes[1,0].plot(sol.t, c[:,1]) 
    axes[1,0].set_title('Middle Tissue') 
    axes[1,1].plot(sol.t, c[:,2]) 
    axes[1,1].set_title('Upper Tissue') 
    plt.tight_layout() 
    plt.show(block=True) 
 
    y = lambda x : alb_sys(0, x, Paramsu) 
    steady = fsolve(y, i0) 
    #print(steady) 
    
    
    with open('foo.csv', 'r') as foo:
        ll = []; reader = csv.reader(foo)
        for row in reader:
            ll.append(row)

    test = np.array(ll).astype(float)
    
    # Cort_p
    ft = lambda x : Cort_k.cort_k(0, x, Paramsu)

    print(np.apply_along_axis(ft, axis=1, arr=test)) 
    

    i0 = [val for val in stateVarsu.values()]   
    
    sol = solve_ivp(DE, [0, 6000], i0, args=(Paramsu,),
             method='LSODA')
    
    np.savetxt("foo.csv", sol.y.T, delimiter=',')
    fig = plot_sols(sol.t, sol.y.T)
    plt.show(block=True)
    """

    return 0


if __name__ == '__main__':
    main() 





