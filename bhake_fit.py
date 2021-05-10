import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
import lmfit
import pdb
import copy
import csv
import DEModel as DEM
import json
import Albumin as ALB
import CBG2 as CBG
import pprint 

"""
Script to run the HPA fit and and analysis
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
           * (((x[0]/v['v_p']) ** v['omega'])
           / ((x[0]/v['v_p']) ** v['omega'] + v['cc'] ** v['omega']))
           * v['v_p'] #* (x[0]/v['v_p']) * 
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
           # diffusion to interstitium
           - v['Dcort_pi'] * ((x[2]/v['v_p']) - (x[3]/v['v_i']))
           # diffusion to the kidney
           - v['Dcort_pk'] * ((x[2]/v['v_p']) - (x[4]/v['v_k']))            
           # Source from Vm1
           +  (v['Vm1'] * v['MWCort'] * v['v_p'] 
             * (x[5]/(v['v_p'] * v['MWsone']))
             / ((x[5]/(v['v_p'] * v['MWsone'])) + v['km1']))
           # Urine
           - v['ucp'] * x[2]     
           # Albumin Binding
           + (-v['konA'] *(1/v['MWAlb']) * (x[14]/v['v_p']) * (x[2]/v['v_p'])
             * v['tonano']
           # albumin off binding
           + v['koffA'] * (x[15]/v['v_p']) * (v['MWCort']/v['MWAC'])
           # CBG binding
           - v['konC'] * (1/v['MWCBG']) * (x[22]/v['v_p']) * (x[2]/v['v_p'])
             * v['mictonano']
           # CBG off binding
           + v['koffC'] * (x[23]/v['v_p']) * (v['MWCort']/v['MWCC']))
             * v['v_p']
          )
    return rhs



##########################################################################
# Interstitial Cortisol functions                                        #
##########################################################################

# 3
def cort_i(t, x, v):
    rhs = (
           # diffusion
           + v['Dcort_pi'] * ((x[2]/v['v_p']) - (x[3]/v['v_i']))
                         
           ##################### Albumin #####################

           # ON TL
           + (- v['konA'] * (1/v['MWAlb']) * (x[7]/v['v_tl'])
           * (x[3]/v['v_i']) * v['tonano']
           # OFF TL
           + v['koffA'] * (x[8]/v['v_tl']) * (v['MWCort']/v['MWAC']))
           * v['v_tl']
              
           # ON TM
           + (- v['konA'] * (1/v['MWAlb']) * (x[9]/v['v_tm'])
             * (x[3]/v['v_i']) * v['tonano']
           # OFF TM
           + v['koffA'] * (x[10]/v['v_tm']) * (v['MWCort']/v['MWAC']))
             * v['v_tm']

           # ON TU
           + (- v['konA'] * (1/v['MWAlb']) * (x[11]/v['v_tu'])
            * (x[3]/v['v_i'])  * v['tonano']
           # OFF TU
           + v['koffA'] * (x[12]/v['v_tu']) * (v['MWCort']/v['MWAC']))
            * v['v_tu']

           ##################### CBG #####################
             
           # ON TL
           + (- v['konC'] * (1/v['MWCBG']) * (x[16]/v['v_tl'])
             * (x[3]/v['v_i']) * v['mictonano']
           # OFF TL
           + v['koffC'] * (x[17]/v['v_tl']) * (v['MWCort']/v['MWCC']))
             * v['v_tl']

           # ON TM
           + (- v['konC'] * (1/v['MWCBG']) * (x[18]/v['v_tm'])
             * (x[3]/v['v_i']) * v['mictonano']
           # ON TM
           + v['koffC'] * (x[19]/v['v_tm']) * (v['MWCort']/v['MWCC']))
             * v['v_tm']

           # ON TU
           + (- v['konC'] * (1/v['MWCBG']) * (x[20]/v['v_tu'])
             * (x[3]/v['v_i']) * v['mictonano'] 
           # OFF TU
           + v['koffC'] * (x[21]/v['v_tu']) * (v['MWCort']/v['MWCC']))
             * v['v_tu']
          )
    return rhs


##########################################################################
# Cortisol in kidney and other cells                                     #
##########################################################################

# 4
def cort_k(t, x, v):
    rhs = (
           # Diffusion from plasma
           + v['Dcort_pk'] * ((x[2]/v['v_p']) - (x[4]/v['v_k']))
           # 11B-HSD2
           - (v['Vm2'] * v['v_k'] * v['MWCort']
           * ((x[4]/(v['v_k'] * v['MWCort'])) 
           /((x[4]/(v['v_k'] * v['MWCort'])) + v['km2'])))
          )
    return rhs


##########################################################################
# Cortisone in kidney and other cells                                     #
##########################################################################

# 5
def Cortisonep(t, x, v):
    rhs = (
           # Constant source 
           v['kcp']
           # Diffusion from kidney
           - v['Dcort_pk'] *((x[5]/v['v_p']) 
             - (x[6]/v['v_k']))
           # Vm1 sink
           - (v['Vm1'] * v['MWCort'] * v['v_p'] 
             * (x[5]/(v['v_p'] * v['MWsone']))
             / ((x[5]/(v['v_p'] * v['MWsone'])) + v['km1']))
           # Urine
           - v['usp'] * x[5]
          )
    return rhs
        
# 6        
def Cortisonek(t, x, v):
    rhs= (
          # Diffusion to plasma
          + v['Dcort_pk'] *((x[5]/v['v_p']) - (x[6]/v['v_k']))
          # 11B-HSD2
          + (v['Vm2'] * v['v_k'] * v['MWsone']
            * ((x[4]/(v['v_k'] * v['MWCort'])) 
            /((x[4]/(v['v_k'] * v['MWCort'])) + v['km2'])))
         )
    return rhs

##########################################################################
# albumin eqs                                                            #
##########################################################################  

Vars = {
        "v\['Cort_pc'\]"      : "(x[2]/v['v_p'])", 
        "v\['Cort_ic'\]"      : "(x[3]/v['v_i'])",
        "v\['Alb_tlc'\]"      : "(x[7]/v['v_tl'])", 
        "v\['Alb-Cort_tlc'\]" : "(x[8]/v['v_tl'])",
        "v\['Alb_tmc'\]"      : "(x[9]/v['v_tm'])",
        "v\['Alb-Cort_tmc'\]" : "(x[10]/v['v_tm'])",
        "v\['Alb_tuc'\]"      : "(x[11]/v['v_tu'])",
        "v\['Alb-Cort_tuc'\]" : "(x[12]/v['v_tu'])",
        "v\['Alb_pc'\]"       : "(x[14]/v['v_p'])",
        "v\['Alb-Cort_pc'\]"  : "(x[15]/v['v_p'])",       
        "v\['CBG_tlc'\]"      : "(x[16]/v['v_tl'])",
        "v\['CBG-Cort_tlc'\]" : "(x[17]/v['v_tl'])",
        "v\['CBG_tmc'\]"      : "(x[18]/v['v_tm'])",
        "v\['CBG-Cort_tmc'\]" : "(x[19]/v['v_tm'])",
        "v\['CBG_tuc'\]"      : "(x[20]/v['v_tu'])",
        "v\['CBG-Cort_tuc'\]" : "(x[21]/v['v_tu'])",
        "v\['CBG_pc'\]"       : "(x[22]/v['v_p'])",
        "v\['CBG-Cort_pc'\]"  : "(x[23]/v['v_p'])",
        "v\['COPEffect'\]"    : 'x[13]',
        "v\['Alb_tl'\]": 'x[7]', "v\['Alb-Cort_tl'\]": 'x[8]',
        "v\['Alb_tm'\]": 'x[9]', "v\['Alb-Cort_tm'\]": 'x[10]',
        "v\['Alb_tu'\]": 'x[11]', "v\['Alb-Cort_tu'\]": 'x[12]',
        "v\['Alb_p'\]" : 'x[14]', "v\['Alb-Cort_p'\]" : 'x[15]',       
        "v\['CBG_tl'\]": 'x[16]', "v\['CBG-Cort_tl'\]": 'x[17]',
        "v\['CBG_tm'\]": 'x[18]', "v\['CBG-Cort_tm'\]": 'x[19]',
        "v\['CBG_tu'\]": 'x[20]', "v\['CBG-Cort_tu'\]": 'x[21]',
        "v\['CBG_p'\]" : 'x[22]', "v\['CBG-Cort_p'\]" : 'x[23]',
        } 


# Torso Lower Lymph 
TL_LymphAlbRate = DEM.DeEquation(ALB.TL_LymphAlbRate_body)
TL_LymphAlbRate.make_fun('tl_lymphalbrate', 
    ALB.TL_LymphAlbRate_head, DEM.ret_str, Vars)

TL_LymphACRate = DEM.DeEquation(ALB.TL_LymphACRate_body)
TL_LymphACRate.make_fun('tl_lymphacrate', 
    ALB.TL_LymphACRate_head, DEM.ret_str, Vars)

# Torso Middle Lymph
TM_LymphAlbRate = DEM.DeEquation(ALB.TM_LymphAlbRate_body)
TM_LymphAlbRate.make_fun('tm_lymphalbrate', 
    ALB.TM_LymphAlbRate_head, DEM.ret_str, Vars)

TM_LymphACRate = DEM.DeEquation(ALB.TM_LymphACRate_body)
TM_LymphACRate.make_fun('tm_lymphacrate', 
    ALB.TM_LymphACRate_head, DEM.ret_str, Vars)

# Torso Lower Lymph
TU_LymphAlbRate = DEM.DeEquation(ALB.TU_LymphAlbRate_body)
TU_LymphAlbRate.make_fun('tu_lymphalbrate', 
    ALB.TU_LymphAlbRate_head, DEM.ret_str, Vars)

TU_LymphACRate = DEM.DeEquation(ALB.TU_LymphACRate_body)
TU_LymphACRate.make_fun('tu_lymphacrate', 
    ALB.TU_LymphACRate_head, DEM.ret_str, Vars)


# Torso Lower Cap
TL_CapAlbRate = DEM.DeEquation(ALB.TL_CapAlbRate_body)
TL_CapAlbRate.make_fun('tl_capalbrate', 
    ALB.TL_CapAlbRate_head, DEM.ret_str, Vars)

TL_CapACRate = DEM.DeEquation(ALB.TL_CapACRate_body)
TL_CapACRate.make_fun('tl_capacrate', 
    ALB.TL_CapACRate_head, DEM.ret_str, Vars)

# Torso Middle Cap
TM_CapAlbRate = DEM.DeEquation(ALB.TM_CapAlbRate_body)
TM_CapAlbRate.make_fun('tm_capalbrate', 
    ALB.TM_CapAlbRate_head, DEM.ret_str, Vars)

TM_CapACRate = DEM.DeEquation(ALB.TM_CapACRate_body)
TM_CapACRate.make_fun('tm_capacrate', 
    ALB.TM_CapACRate_head, DEM.ret_str, Vars)

# Torso Upper Cap
TU_CapAlbRate = DEM.DeEquation(ALB.TU_CapAlbRate_body)
TU_CapAlbRate.make_fun('tu_capalbrate', 
    ALB.TU_CapAlbRate_head, DEM.ret_str, Vars)

TU_CapACRate = DEM.DeEquation(ALB.TU_CapACRate_body)
TU_CapACRate.make_fun('tu_capacrate', 
    ALB.TU_CapACRate_head, DEM.ret_str, Vars)

# Partial pressure
Albpcop = DEM.DeEquation(ALB.Albpcop_body)
Albpcop.make_fun('albpcop', ALB.Albpcop_head, 
    DEM.ret_str, Vars)

# Colloid Oncotic Pressure
COP = DEM.DeEquation(ALB.COP_body)
COP.make_fun('cop', ALB.COP_head, 
    DEM.ret_str, Vars)

# COP effect
COPE = DEM.DeEquation(ALB.COPE_body)
COPE.make_fun('cope', ALB.COPE_head, 
    DEM.ret_str, Vars)

# Torso Lower Balance
Albtl = DEM.DeEquation(ALB.Albtl_body)
Albtl.make_fun('alb_tl', ALB.Albtl_head, DEM.ret_str, Vars)

AlbCorttl = DEM.DeEquation(ALB.AlbCorttl_body)
AlbCorttl.make_fun('albcort_tl', 
    ALB.AlbCorttl_head, DEM.ret_str, Vars)

# Torso Middle
Albtm = DEM.DeEquation(ALB.Albtm_body)
Albtm.make_fun('alb_tm', 
    ALB.Albtm_head, DEM.ret_str, Vars)

AlbCorttm = DEM.DeEquation(ALB.AlbCorttm_body)
AlbCorttm.make_fun('albcort_tm', 
    ALB.AlbCorttm_head, DEM.ret_str, Vars)

# Torso Upper
Albtu = DEM.DeEquation(ALB.Albtu_body)
Albtu.make_fun('alb_tu', 
    ALB.Albtu_head, DEM.ret_str, Vars)

AlbCorttu = DEM.DeEquation(ALB.AlbCorttu_body)
AlbCorttu.make_fun('albcort_tu', 
        ALB.AlbCorttu_head, DEM.ret_str, Vars)

# Plasma
Albp = DEM.DeEquation(ALB.Albp_body)
Albp.make_fun('alb_p', 
    ALB.Albp_head, DEM.ret_str, Vars)

AlbCortp = DEM.DeEquation(ALB.AlbCortp_body)
AlbCortp.make_fun('albcort_p', 
    ALB.AlbCortp_head, DEM.ret_str, Vars)
    
    
Albrates = [TL_LymphAlbRate, TL_LymphACRate, TM_LymphAlbRate, 
            TM_LymphACRate, TU_LymphAlbRate, TU_LymphACRate, 
            TL_CapAlbRate, TL_CapACRate, TM_CapAlbRate, TM_CapACRate, 
            TU_CapAlbRate, TU_CapACRate, Albpcop, COP]

##########################################################################
# CBG eqs                                                                #
##########################################################################
    

# Torso Lower Lymph 
TL_LymphCBGRate = DEM.DeEquation(CBG.TL_LymphCBGRate_body)
TL_LymphCBGRate.make_fun('tl_lymphcbgrate', 
    CBG.TL_LymphCBGRate_head, DEM.ret_str, Vars)

TL_LymphCCRate = DEM.DeEquation(CBG.TL_LymphCCRate_body)
TL_LymphCCRate.make_fun('tl_lymphccrate', 
    CBG.TL_LymphCCRate_head, DEM.ret_str, Vars)

# Torso Middle Lymph
TM_LymphCBGRate = DEM.DeEquation(CBG.TM_LymphCBGRate_body)
TM_LymphCBGRate.make_fun('tm_lymphcbgrate', 
    CBG.TM_LymphCBGRate_head, DEM.ret_str, Vars)

TM_LymphCCRate = DEM.DeEquation(CBG.TM_LymphCCRate_body)
TM_LymphCCRate.make_fun('tm_lymphccrate', 
    CBG.TM_LymphCCRate_head, DEM.ret_str, Vars)

# Torso Lower Lymph
TU_LymphCBGRate = DEM.DeEquation(CBG.TU_LymphCBGRate_body)
TU_LymphCBGRate.make_fun('tu_lymphcbgrate', 
    CBG.TU_LymphCBGRate_head, DEM.ret_str, Vars)

TU_LymphCCRate = DEM.DeEquation(CBG.TU_LymphCCRate_body)
TU_LymphCCRate.make_fun('tu_lymphccrate', 
    CBG.TU_LymphCCRate_head, DEM.ret_str, Vars)


# Torso Lower Cap
TL_CapCBGRate = DEM.DeEquation(CBG.TL_CapCBGRate_body)
TL_CapCBGRate.make_fun('tl_capcbgrate', 
    CBG.TL_CapCBGRate_head, DEM.ret_str, Vars)

TL_CapCCRate = DEM.DeEquation(CBG.TL_CapCCRate_body)
TL_CapCCRate.make_fun('tl_capccrate', 
    CBG.TL_CapCCRate_head, DEM.ret_str, Vars)

# Torso Middle Cap
TM_CapCBGRate = DEM.DeEquation(CBG.TM_CapCBGRate_body)
TM_CapCBGRate.make_fun('tm_capcbgrate', 
    CBG.TM_CapCBGRate_head, DEM.ret_str, Vars)

TM_CapCCRate = DEM.DeEquation(CBG.TM_CapCCRate_body)
TM_CapCCRate.make_fun('tm_capccrate', 
    CBG.TM_CapCCRate_head, DEM.ret_str, Vars)

# Torso Upper Cap
TU_CapCBGRate = DEM.DeEquation(CBG.TU_CapCBGRate_body)
TU_CapCBGRate.make_fun('tu_capcbgrate', 
    CBG.TU_CapCBGRate_head, DEM.ret_str, Vars)

TU_CapCCRate = DEM.DeEquation(CBG.TU_CapCCRate_body)
TU_CapCCRate.make_fun('tu_capccrate', 
    CBG.TU_CapCCRate_head, DEM.ret_str, Vars)


# Torso Lower Balance
CBGtl = DEM.DeEquation(CBG.CBGtl_body)
CBGtl.make_fun('cbg_tl', CBG.CBGtl_head, DEM.ret_str, Vars)

CBGCorttl = DEM.DeEquation(CBG.CBGCorttl_body)
CBGCorttl.make_fun('cbgcort_tl', 
    CBG.CBGCorttl_head, DEM.ret_str, Vars)

# Torso Middle
CBGtm = DEM.DeEquation(CBG.CBGtm_body)
CBGtm.make_fun('cbg_tm', 
    CBG.CBGtm_head, DEM.ret_str, Vars)

CBGCorttm = DEM.DeEquation(CBG.CBGCorttm_body)
CBGCorttm.make_fun('cbgcort_tm', 
    CBG.CBGCorttm_head, DEM.ret_str, Vars)

# Torso Upper
CBGtu = DEM.DeEquation(CBG.CBGtu_body)
CBGtu.make_fun('cbg_tu', 
    CBG.CBGtu_head, DEM.ret_str, Vars)

CBGCorttu = DEM.DeEquation(CBG.CBGCorttu_body)
CBGCorttu.make_fun('cbgcort_tu', 
    CBG.CBGCorttu_head, DEM.ret_str, Vars)

# Plasma
CBGp = DEM.DeEquation(CBG.CBGp_body)
CBGp.make_fun('cbg_p', 
    CBG.CBGp_head, DEM.ret_str, Vars)

CBGCortp = DEM.DeEquation(CBG.CBGCortp_body)
CBGCortp.make_fun('cbgcort_p', 
    CBG.CBGCortp_head, DEM.ret_str, Vars)
    
    
CBGrates = [TL_LymphCBGRate, TL_LymphCCRate, TM_LymphCBGRate, 
            TM_LymphCCRate, TU_LymphCBGRate, TU_LymphCCRate, 
            TL_CapCBGRate, TL_CapCCRate, TM_CapCBGRate, TM_CapCCRate, 
            TU_CapCBGRate, TU_CapCCRate]


##########################################################################
# Make sys                                                               #
##########################################################################

def sys(t, x, params):
    
    v = params.valuesdict()   
    # Albumin Rates
    TL_LAR  = TL_LymphAlbRate.tl_lymphalbrate(0, x, v)
    TL_LACR = TL_LymphACRate.tl_lymphacrate(0, x, v)
    TM_LAR  = TM_LymphAlbRate.tm_lymphalbrate(0, x, v)
    TM_LACR = TM_LymphACRate.tm_lymphacrate(0, x, v)
    TU_LAR  = TU_LymphAlbRate.tu_lymphalbrate(0, x, v)
    TU_LACR = TU_LymphACRate.tu_lymphacrate(0, x, v)
    TL_CAR  = TL_CapAlbRate.tl_capalbrate(0, x, v)
    TL_CACR = TL_CapACRate.tl_capacrate(0, x, v)
    TM_CAR  = TM_CapAlbRate.tm_capalbrate(0, x, v)
    TM_CACR = TM_CapACRate.tm_capacrate(0, x, v)
    TU_CAR  = TU_CapAlbRate.tu_capalbrate(0, x, v)
    TU_CACR = TU_CapACRate.tu_capacrate(0, x, v)
    Apcop   = Albpcop.albpcop(x, v)
    ACOP    = COP.cop(x, v, Apcop)

    # Albumin Base equations
    albtl     = Albtl.alb_tl(0, x, v, TL_CAR, TL_LAR)
    albcorttl = AlbCorttl.albcort_tl(0, x, v, TL_CACR, TL_LACR)
    albtm     = Albtm.alb_tm(0, x, v, TM_CAR, TM_LAR)
    albcorttm = AlbCorttm.albcort_tm(0, x, v, TM_CACR, TM_LACR)
    albtu     = Albtu.alb_tu(0, x, v, TU_CAR, TU_LAR)
    albcorttu = AlbCorttu.albcort_tu(0, x, v, TU_CACR, TU_LACR)
    cope      = COPE.cope(0, x, v, ALB.copeCurve, ACOP)
    T_LAR     = TL_LAR + TM_LAR + TU_LAR
    T_CAR     = TL_CAR + TM_CAR + TU_CAR
    albp      = Albp.alb_p(0, x, v, T_LAR,T_CAR)
    T_LACR    = TL_LACR + TM_LACR + TU_LACR
    T_CACR    = TL_CACR + TM_CACR + TU_CACR
    albcortp  = AlbCortp.albcort_p(0, x, v, T_LACR, T_CACR)    

    # CBG Rates
    TL_LCR  = TL_LymphCBGRate.tl_lymphcbgrate(0, x, v)
    TL_LCCR = TL_LymphCCRate.tl_lymphccrate(0, x, v)
    TM_LCR  = TM_LymphCBGRate.tm_lymphcbgrate(0, x, v)
    TM_LCCR = TM_LymphCCRate.tm_lymphccrate(0, x, v)
    TU_LCR  = TU_LymphCBGRate.tu_lymphcbgrate(0, x, v)
    TU_LCCR = TU_LymphCCRate.tu_lymphccrate(0, x, v)
    TL_CCR  = TL_CapCBGRate.tl_capcbgrate(0, x, v)
    TL_CCCR = TL_CapCCRate.tl_capccrate(0, x, v)
    TM_CCR  = TM_CapCBGRate.tm_capcbgrate(0, x, v)
    TM_CCCR = TM_CapCCRate.tm_capccrate(0, x, v)
    TU_CCR  = TU_CapCBGRate.tu_capcbgrate(0, x, v)
    TU_CCCR = TU_CapCCRate.tu_capccrate(0, x, v)

    # CBG Base equations
    cbgtl     = CBGtl.cbg_tl(0, x, v, TL_CCR, TL_LCR)
    cbgcorttl = CBGCorttl.cbgcort_tl(0, x, v, TL_CCCR, TL_LCCR)
    cbgtm     = CBGtm.cbg_tm(0, x, v, TM_CCR, TM_LCR)
    cbgcorttm = CBGCorttm.cbgcort_tm(0, x, v, TM_CCCR, TM_LCCR)
    cbgtu     = CBGtu.cbg_tu(0, x, v, TU_CCR, TU_LCR)
    cbgcorttu = CBGCorttu.cbgcort_tu(0, x, v, TU_CCCR, TU_LCCR)
    T_LCR     = TL_LCR + TM_LCR + TU_LCR
    T_CCR     = TL_CCR + TM_CCR + TU_CCR
    cbgp      = CBGp.cbg_p(0, x, v, T_LCR,T_CCR)
    T_LCCR    = TL_LCCR + TM_LCCR + TU_LCCR
    T_CCCR    = TL_CCCR + TM_CCCR + TU_CCCR
    cbgcortp  = CBGCortp.cbgcort_p(0, x, v, T_LCCR, T_CCCR) 


    CRH   = crf(t, x, v)
    ACTH  = acth(t, x, v)
    Cortp = cort_p(t, x, v)
    Corti = cort_i(t, x, v)
    Cortk = cort_k(t, x, v)
    Sonep = Cortisonep(t, x, v)
    Sonek = Cortisonek(t, x, v)

    return [
            CRH, ACTH, Cortp, Corti, Cortk, Sonep, Sonek,
            albtl, albcorttl, albtm, albcorttm, albtu, albcorttu,
            cope, albp, albcortp, cbgtl, cbgcorttl, cbgtm, cbgcorttm,
            cbgtu, cbgcorttu, cbgp, cbgcortp
           ]


def solve_steady(func, guess, args):
    # assume DE starts with t argument
    f = lambda x, v: func(0, x, v)
    sol = sp.optimize.fsolve(f, guess, args)
    return sol

##########################################################################
# Parameters                                                             #
##########################################################################


k1expr =  (
           '(w1 * x0bar * v_p)' 
           '/ (1 + zeta * ((x2bar ** alpha)'
           '  / (x2bar ** alpha + c ** alpha))'
           # Negative hill feedback by cort_p
           '- (psi * ((x2bar ** gamma)'
           '/ (x2bar ** gamma + c3 ** gamma))))' 
          )


k2expr =  (
           '(w2 * x1bar * v_p) /((1 - ((rho * x2bar ** alpha) /'
           '(x2bar ** alpha + c ** alpha)))'
            '* ((x0bar ** omega)'
            '/ (x0bar ** omega + cc ** omega)) * v_p)'
          )


k3expr = (
           '-1 * (- Dcort_pi * (x2bar - x3bar)'
           # diffusion to the kidney
           '- Dcort_pk * (x2bar - x4bar)'
           # Source from Vm1
           '+ (Vm1 * MWCort * v_p'
           '  * (x5bar/(v_p * MWsone))'
           '  / ((x5bar/(v_p * MWsone)) + km1))'
           # Urine
           '- ucp * x2bar * v_p'
           # Albumin Binding
           '+ (-konA *(1/MWAlb) * x14bar * x2bar'
           '  * tonano'
           # albumin off binding
           '+ koffA * x15bar * (MWCort/MWAC)'
           # CBG binding
           '- konC * (1/MWCBG) * x22bar * x2bar'
           '  * mictonano'
           # CBG off binding
           '+ koffC * x23bar * (MWCort/MWCC))'
           '  * v_p) / (x1bar * v_p)'
         )



params = lmfit.Parameters()
# Constraints
params.add('x0bar', value=7.66, vary=False)             
params.add('x1bar', value=21, min=14.8, max=27.2)
params.add('x2bar', value=3.06, min=2.115, max=3.995)
params.add('x3bar', value=61.1, vary=False) # inter cort
params.add('x4bar', value=61.1, vary=False) # kidney cort
params.add('x5bar', value=61.1, vary=False ) # plasma cortisone MASS
params.add('x14bar', value=61.1, vary=False) # Alb
params.add('x15bar', value=61.1, vary=False) # AlbCort
params.add('x22bar', value=61.1, vary=False) # CBG
params.add('x23bar', value=61.1, vary=False) # CBGCort

# MWs and misc
params.add('MWAC', value=66799.46, vary=False)
params.add('MWAlb', value=66437.0, vary=False)
params.add('MWCBG', value=52000.0, vary=False)
params.add('MWCC', value=52362.46, vary=False)
params.add('MWCort', value=362.46, vary=False)
params.add('MWsone', value=360.45, vary=False)
params.add('mictonano', value=1000.0, vary=False)
params.add('tonano', value=1000000000.0, vary=False)
params.add('v_i', value=12153.9, vary=False)
params.add('v_k', value=5.024, vary=False)
params.add('v_p', value=3261.72, vary=False)
params.add('v_tl', value=3806.42, vary=False)
params.add('v_tm', value=5790.4, vary=False)
params.add('v_tu', value=2554.52, vary=False)
params.add('Ald', value=0.121, vary=False) # Pointless
params.add('Kd_mr', value=0.75, vary=False) # Pointless

# CRH Equations
params.add('alpha', value=7, min=1, max=8)
params.add('c', value=3, min=0.5, max=7.5)
params.add('zeta', value=2, vary=False)
params.add('c3', value=1.42, vary=False)
params.add('gamma', value=3, vary=False)
params.add('psi', value=0.5, vary=False)
params.add('w1', value=0.1731, vary=False)
params.add('k1', expr=k1expr)

# ACTH
params.add('rho', value=0.8, min=0, max=0.9)
params.add('w2', value=0.0348, vary=False)
params.add('omega', value=3, min=1, max=8)
params.add('cc', value=7, min=5, max=15)
params.add('k2', expr=k2expr)

# Cort Plasma
params.add('Dcort_pi', value=10.8, vary=False)
params.add('Dcort_pk', value=9.0, min=5, max=12)
params.add('Vm1', value=0.004, vary=False)
params.add('koffA', value=41.58, vary=False)
params.add('koffC', value=0.88, vary=False)
params.add('konA', value=0.12474, vary=False)
params.add('konC', value=3.5, vary=False)
params.add('km1', value=21, vary=False)
params.add('ucp', value=0.0016, vary=False)
params.add('k3', value=0.00065, min=1e-6, max=1)                #expr=k3expr)

# Cort Kidney
params.add('Vm2', value=0.04, vary=False)
params.add('km2', value=0.04, vary=False)

# Cortisone plasma
params.add('kcp', value=0.001, vary=False)
params.add('usp', value=0.00503, vary=False)

# Albumin
params.add('AGPpCOP', value=0.9426, vary=False)
params.add('AlbPerm', value=0, vary=False)
params.add('Deg_a', value=0.005, vary=False)
params.add('GFR', value=103.66, vary=False)
params.add('HetapCOP', value=0, vary=False)
params.add('IgpCOP', value=4.49, vary=False)
params.add('Klp', value=0.1, vary=False)
params.add('LAPerm', value=5.4, vary=False)
params.add('MAPerm', value=6.4, vary=False)
params.add('ModPerm', value=0, vary=False)
params.add('TL_LW_R', value=7.908, vary=False)
params.add('TM_LW_R', value=14.41, vary=False)
params.add('TU_LW_R', value=6.24, vary=False)
params.add('Temp_K', value=309.772, vary=False)
params.add('UAPerm', value=5.4, vary=False)
params.add('k4', value=0.005, vary=False)
params.add('sigA', value=0.8, vary=False)

# CBG
params.add('half', value=9.62e-05, vary=False)
params.add('k5', value=0.00423, vary=False)
params.add('sigC', value=0.6, vary=False)

# Bolus
params.add('A', value=6.467, vary=False)
params.add('B', value=3.455, vary=False)            
params.add('alp', value=0.081, vary=False)          
params.add('bet', value=0.01359, vary=False)           

# Ensure the correct values are in place
# Import Master Data
with open('MaleMaster.json') as json_file:
    Male_Master = json.load(json_file) 

pss = {key : val['Magnitude'] for key, val in Male_Master.items()
        if val['Type'] == 'Parameter'}

for key, val in pss.items():
    params[key].value = val

"""
print('{:<15}   {:<15} | {:<15}   {:<15}'.format('Params', 
    'Value', 'Master', 'Value'))
for key, val in params.valuesdict().items():
    try:
        pss[key]
        print('{:<15}   {:<15.4f} | {:<15}   {:<15.4f}'.format(key, 
            val, key, pss[key]))
    except KeyError:
        pass
"""
##########################################################################
# Generic functions                                                      #
##########################################################################

def toCortT_ugdl(x, v):
    cort = x[:,2] / (v['MWCort'].value * v['v_p'].value)# to moles
    albcort = x[:,15] / (v['MWAC'].value * v['v_p'].value)# to moles
    cbgcort = x[:,-1] / (v['MWCC'].value * v['v_p'].value)# to moles
    # Total cortisol converted to ug/dl
    totalCort = (cort + albcort + cbgcort) * v['MWCort'].value * 0.1

    return totalCort

def get_csv_data(file_str, first_row_header=True): 
    """ Read in the csv data """
    with open(file_str) as ac:
        reader = csv.reader(ac, delimiter=',')
        if first_row_header:
            next(reader)
        
        arr = np.array([row for row in reader]).astype('float64')
    return arr


def m_f_data(data_arr, male=[1]):
    """ Get male and female data based on columns """
    # assumes first colum is time
    female = [i for i in range(data_arr.shape[1]) 
            if i not in male and i != 0]
    male_arr = np.take(data_arr, male, axis=1) 
    female_arr = np.take(data_arr, female, axis=1)
    time = data_arr[:,0]
    return [time, male_arr, female_arr]



def uniform_p_sample(parms):
    for p in parms:
        if parms[p].vary != False:
            parms[p].value = np.random.uniform(parms[p].min, 
                    parms[p].max)
    return parms


##########################################################################
# Galluci                                                                #
##########################################################################

############################# file data ###################################
cort_galluci_d = ('/home/anderson/Documents/
        UW_AM_Research/Code/HPA/Cort_walluci.csv')


gall_data = get_csv_data(cort_galluci_d)

ts, gall_male, gall_female =  m_f_data(gall_data, male=[1, 3])

# normalize data to the mean steady state of 
ssACTH = 21
ssCortF = 3.06 
ssAlbCort = 1176.367 * (pss['MWCort']/pss['MWAC'])
ssCBGCort = 10008.7841 * (pss['MWCort']/pss['MWCC'])
ssCortT = (ssCortF + ssAlbCort + ssCBGCort) * 0.1
gall_male[:,0] = (ssACTH - gall_male[0,0]) + gall_male[:,0]
gall_male[:,1] = (ssCortT - gall_male[0,1]) + gall_male[:,1]
gall_female[:,0] = (21 - gall_female[0,0]) + gall_female[:,0]
gall_female[:,1] = (6.11 - gall_female[0,1]) + gall_female[:,1]


######################### model and objective ############################

def g(fun,t, ts, x0, ps):
    """ Wrapper for solving an IVP """
    sol = sp.integrate.solve_ivp(fun, t, x0, method='LSODA',
            t_eval=ts, args=(ps,))
    return sol


def fit(cost_func, params, func_args, op_method):
    """ Wrapper for fitting a model """
    minner = lmfit.Minimizer(cost_func, params, fcn_args=(func_args))
    result = minner.minimize(method=op_method)
    return result


def residual_galluci(ps, fun, t, ts, data, x0):
    """ 
    Galluci objective function 
    ps: Parameters object
    fun: system function
    t: interval
    ts: evaluation times
    data: experimental data points
    """
    z = g(fun, t, ts, x0, ps)
    acth = z.y.T[:,1]/ps['v_p'].value
    # Total cortisol converted to ug/dl
    totalCort = toCortT_ugdl(z.y.T, ps)
    # Model results
    model = np.column_stack((acth, totalCort))
    return (model - data).ravel()



def fit_galluci(gall_data, t, ts, x0, params, method='leastsq'):
    # fit the male model
    result = fit(residual_galluci, params, 
        (sys, t, ts, gall_data, x0), method)
    return result


##########################################################################
# Main                                                                   #
##########################################################################


def main(params, ts):
    # set up initial number of runs
    num_runs = 1; num_parm = len(params)
    
    # set up initial conditions
    with open('MaleSteady.json') as f:
        x0 = json.load(f)
    
    x0_m = [x0['Mass'][var] for var in DEM.var_order]
    # Set CRH IC
    x0_m[0] = CRF_bolus(0, params.valuesdict()) * 1000 * params['v_p'].value  # Units pg/ml
    
    # Set the mean params for the expresions
    params['x0bar'].value  = x0['Concentration']['CRH']      
    params['x1bar'].value  = x0['Concentration']['ACTH']
    params['x2bar'].value  = x0['Concentration']['Cort_p']
    params['x3bar'].value  = x0['Concentration']['Cort_i']
    params['x4bar'].value  = x0['Concentration']['Cort_k']
    params['x5bar'].value  = x0['Mass']['Cortisone_p'] 
    params['x14bar'].value = x0['Concentration']['Alb_p']
    params['x15bar'].value = x0['Concentration']['Alb-Cort_p']
    params['x22bar'].value = x0['Concentration']['CBG_p']
    params['x23bar'].value = x0['Concentration']['CBG-Cort_p']
    
    
    # Set 
    #x0_f = [gall_female[0,0], gall_female[0,1] * 10] # convert units

    # main loop
    for i in range(num_runs):
        # sample parms 
        params = uniform_p_sample(params)
        params['k3'].value = 0.06
        # Fit
        male_res = fit_galluci(gall_male, [ts[0], ts[-1]], 
                ts, x0_m, params)
        
        plt.plot(ts, gall_male[:,0], 'o', label='ACTH Data pg/ml')
        plt.plot(ts, gall_male[:,1], 'o', label='Cort Data ug/dl')
        final = male_res.residual.reshape(gall_male.shape) + gall_male
        
        plt.plot(ts, final[:,0], label='ACTH model pg/ml')
        plt.plot(ts, final[:,1], label='Cort model ug/dl')
        
        plt.xlabel('Time (min)'); plt.ylabel('Concentration')
        plt.legend(); plt.show(block=True)
    
        print('\nMale Fit')
        lmfit.report_fit(male_res)
        
        return [male_res.params, final]

def cort_p2(t, x, v):
    # stimulation
    t1 = v['k3'] * (x[1]/v['v_p']) * v['v_p']
    # diffusion to interstitium
    t2 = - v['Dcort_pi'] * ((x[2]/v['v_p']) - (x[3]/v['v_i']))
    # diffusion to the kidney
    t3 = - v['Dcort_pk'] * ((x[2]/v['v_p']) - (x[4]/v['v_k']))            
    # Source from Vm1
    t4 =  (v['Vm1'] * v['MWCort'] * v['v_p'] 
             * (x[5]/(v['v_p'] * v['MWsone']))
             / ((x[5]/(v['v_p'] * v['MWsone'])) + v['km1']))
    # Urine
    t5 = - v['ucp'] * x[2]     
    # Albumin Binding
    t6 = ((-v['konA'] *(1/v['MWAlb']) * (x[14]/v['v_p']) * (x[2]/v['v_p'])
             * v['tonano']
    # albumin off binding
    + v['koffA'] * (x[15]/v['v_p']) * (v['MWCort']/v['MWAC'])
    # CBG binding
    - v['konC'] * (1/v['MWCBG']) * (x[22]/v['v_p']) * (x[2]/v['v_p'])
    * v['mictonano']
    # CBG off binding
    + v['koffC'] * (x[23]/v['v_p']) * (v['MWCort']/v['MWCC']))
    * v['v_p'])
    s = t1 + t2 + t3 + t4 + t5 + t6    

    return [s, [t1, t2, t3, t4, t5, t6]]



# set up initial conditions
with open('MaleSteady.json') as f:
    x0 = json.load(f)

x0_m = [x0['Mass'][var] for var in DEM.var_order]
# Set CRH IC
x0_m[0] = CRF_bolus(0, params.valuesdict()) * 1000 * params['v_p'].value  # Units pg/ml

p, x = main(params, ts)

# Simulate the System
sols = DEM.np.zeros((3000, 24, 1))
tspan = [0,120]
evals = DEM.np.linspace(0,120,3000)
for i in range(1):
    sol = DEM.solve_ivp(sys, tspan, x0_m, method='LSODA', t_eval=evals,
            args=(p,))
    sols[:,:,i] = sol.y.T

# Calculate the concentrations

def calc_conc(x, v):
    c = DEM.np.zeros(x.shape)
    c[:,0,:]  = x[:,0,:]/v['v_p']  ; c[:,1,:]  = x[:,1,:]/v['v_p']
    c[:,2,:]  = x[:,2,:]/v['v_p']  ; c[:,3,:]  = x[:,3,:]/v['v_i']
    c[:,4,:]  = x[:,4,:]/v['v_k']  ; c[:,5,:]  = x[:,5,:]/v['v_p']
    c[:,6,:]  = x[:,6,:]/v['v_k']  ; c[:,7,:]  = x[:,7,:]/v['v_tl']
    c[:,8,:]  = x[:,8,:]/v['v_tl'] ; c[:,9,:]  = x[:,9,:]/v['v_tm']
    c[:,10,:] = x[:,10,:]/v['v_tm']; c[:,11,:] = x[:,11,:]/v['v_tu']
    c[:,12,:] = x[:,12,:]/v['v_tu']; c[:,13,:] = x[:,13,:]
    c[:,14,:] = x[:,14,:]/v['v_p'] ; c[:,15,:] = x[:,15,:]/v['v_p']
    c[:,16,:] = x[:,16,:]/v['v_tl']; c[:,17,:] = x[:,17,:]/v['v_tl']
    c[:,18,:] = x[:,18,:]/v['v_tm']; c[:,19,:] = x[:,19,:]/v['v_tm']
    c[:,20,:] = x[:,20,:]/v['v_tu']; c[:,21,:] = x[:,21,:]/v['v_tu']
    c[:,22,:] = x[:,22,:]/v['v_p'] ; c[:,23,:] = x[:,23,:]/v['v_p']

    return c


# Plot the system
def plot_sols(t, sols, labels1, title1):
    """ 
    Plot mutliple sols of the hpa system
    """
    fig, ax = DEM.plt.subplots(4,4, figsize=(9,7), sharex=True)
    lines = []
    for i in range(sols.shape[-1]):
        ax[0,0].plot(t, sols[:,0,i])
        ax[0,0].set_title('CRH (pg/ml)')
        
        ax[0,1].plot(t, sols[:,1,i])
        ax[0,1].set_title('ACTH (pg/ml)')        
        
        ax[0,2].plot(t, sols[:,2,i])
        ax[0,2].set_title('CortP (ng/ml)')
        
        ax[0,3].plot(t, sols[:,3,i])
        ax[0,3].set_title('CortI (ng/ml)')       
        
        ax[1,0].plot(t, sols[:,4,i])
        ax[1,0].set_title('CortK (ng/ml)')
        
        ax[1,1].plot(t, sols[:,5,i])
        ax[1,1].set_title('CortisoneP (ng/ml)')       
        
        ax[1,2].plot(t, sols[:,7,i] + sols[:,9,i] + sols[:,11,i])
        ax[1,2].set_ylim([0.02, 0.06])
        ax[1,2].set_title('AlbI (g/ml)')      
        
        ax[1,3].plot(t, sols[:,8,i] + sols[:,10,i] + sols[:,12,i])
        ax[1,3].set_title('AlbCortI (ng/ml)')
        
        ax[2,0].plot(t, sols[:,14,i])
        ax[2,0].set_ylim([0.02, 0.06])
        ax[2,0].set_title('AlbP (g/ml)')       
        
        ax[2,1].plot(t, sols[:,15,i])
        ax[2,1].set_title('AlbCortP (ng/ml)')
        
        ax[2,2].plot(t, sols[:,16,i] + sols[:,18,i] + sols[:,20,i])
        ax[2,2].set_title('CBGI (ug/ml)')       
        
        ax[2,3].plot(t, sols[:,17,i] + sols[:,19,i] + sols[:,21,i])
        ax[2,3].set_title('CBGCortI (ng/ml)')
        
        ax[3,0].plot(t, sols[:,22,i])
        ax[3,0].set_title('CBGP (pg/ml)')       
        
        ax[3,1].plot(t, sols[:,23,i])
        ax[3,1].set_title('CBGCortP (pg/ml)')
        
        ax[3,2].plot(t, sols[:,13,i])
        ax[3,2].set_ylim([-2, 2])
        ax[3,2].set_title('COPE')    
        
        line, = ax[3,3].plot(t, sols[:,6,i]) 
        ax[3,3].set_title('CortisoneK (ng/ml)')
        
        lines.append(line)

    leg1 = fig.legend(tuple(lines), labels1, title=title1 , 
            bbox_to_anchor=(1.005,0.85),loc='upper right')
    fig.tight_layout()
    fig.subplots_adjust(right=0.87)
    return [fig, leg1]


labels1 = ['1 ug/kg']
title1 = 'CRH IC'
ConcSols = calc_conc(sols, p.valuesdict())
male = plot_sols(evals, ConcSols, tuple(labels1), title1)
male[0].text(0.5, 0.04, 'Time (min)', ha='center')

male[0].savefig('test.png',
        bbox_extra_artists=tuple(male[1:]), bbox_inches='tight')


if __name__ == '__main__':
    pass
    #main(params, ts)

