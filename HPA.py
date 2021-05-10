import DEModel as DEM
import Albumin as ALB
import CBG2 as CBG
import json

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
            #"* ((v['CRH'] ** v['omega'])"
            #"/ (v['CRH'] ** v['omega'] + v['cc'] ** v['omega']))"
           * (x[0]/v['v_p']) * v['v_p']
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
#########################################################################  

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

def sys(t, x, v):
    
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

##########################################################################
# Test the steady state                                                  #
##########################################################################

def solve_steady(func, guess, args):
    # assume DE starts with t argument
    f = lambda x, v: func(0, x, v)
    sol = DEM.fsolve(f, guess, args)
    return sol


# Import Master Data
with open('MaleMaster.json') as json_file:
    Male_Master = json.load(json_file) 

# Calculate Steady State
var_order = ['CRH', 'ACTH', 'Cort_p', 'Cort_i', 'Cort_k',
             'Cortisone_p', 'Cortisone_k', 'Alb_tl', 'Alb-Cort_tl',
             'Alb_tm', 'Alb-Cort_tm', 'Alb_tu','Alb-Cort_tu',
             'COPEffect', 'Alb_p', 'Alb-Cort_p', 'CBG_tl',
             'CBG-Cort_tl', 'CBG_tm', 'CBG-Cort_tm', 'CBG_tu',
             'CBG-Cort_tu', 'CBG_p', 'CBG-Cort_p']

i0 = [Male_Master[i]['Magnitude'] for i in var_order 
        if Male_Master[i]['Type'] == 'State Mass']
ps = {key : val['Magnitude'] for key, val in Male_Master.items()
        if val['Type'] == 'Parameter'}
sys(0, i0, ps)
sol = solve_steady(sys, i0, ps)

# Calculate the Concentration
def calc_conc(x, v):
    c = [i for i in x]
    c[0]  = x[0]/v['v_p']  ; c[1]  = x[1]/v['v_p']
    c[2]  = x[2]/v['v_p']  ; c[3]  = x[3]/v['v_i']
    c[4]  = x[4]/v['v_k']  ; c[5]  = x[5]/v['v_p']
    c[6]  = x[6]/v['v_k']  ; c[7]  = x[7]/v['v_tl']
    c[8]  = x[8]/v['v_tl'] ; c[9]  = x[9]/v['v_tm']
    c[10] = x[10]/v['v_tm']; c[11] = x[11]/v['v_tu']
    c[12] = x[12]/v['v_tu']; c[13] = x[13]
    c[14] = x[14]/v['v_p'] ; c[15] = x[15]/v['v_p']
    c[16] = x[16]/v['v_tl']; c[17] = x[17]/v['v_tl']
    c[18] = x[18]/v['v_tm']; c[19] = x[19]/v['v_tm']
    c[20] = x[20]/v['v_tu']; c[21] = x[21]/v['v_tu']
    c[22] = x[22]/v['v_p'] ; c[23] = x[23]/v['v_p']

    return c


# Check Binding
def calcBinding(x, v):
    # Plasma Mass
    PF = x[2]
    AlbCortP = x[15] * (v['MWCort']/v['MWAC'])
    CBGCortP = x[23] * (v['MWCort']/v['MWCC'])
    total = PF + AlbCortP + CBGCortP
    PFfrac = PF/total; AlbCortPfrac = AlbCortP/total
    CBGCortPfrac = CBGCortP/total
    # Interstitial Mass
    IF = x[3]
    AlbCortI = (x[8] + x[10] + x[12]) * (v['MWCort']/v['MWAC'])
    CBGCortI = (x[17] + x[19] + x[21]) * (v['MWCort']/v['MWCC'])
    totalI = IF + AlbCortI + CBGCortI
    IFfrac = IF/totalI; AlbCortIfrac = AlbCortI/totalI
    CBGCortIfrac = CBGCortI/totalI

    return [[PF, AlbCortP, CBGCortP, PFfrac, AlbCortPfrac, CBGCortPfrac], 
            [IF, AlbCortI, CBGCortI, IFfrac, AlbCortIfrac, CBGCortIfrac]]


def print_ss_report(x, v, var_order, Master):
    strFormat = '{:<25} {:<25} {:<25}'
    # header
    print(strFormat.format(*('Species', 'Value', 'Unit')))
    # get conc
    c = calc_conc(x, v) 
    # print conc
    for i, j in enumerate(var_order):
        print(strFormat.format(*(j, c[i], Master[j]['Unit'])))

    # Binding 
    bind = calcBinding(c, v)
    strFormat2 = '{:<12} {:<12} {:<12} {:<12} {:<12} {:<12}'
    strFormat2N = '{:<12.2f} {:<12.2f} {:<12.2f} {:<12.2f} {:<12.2f} {:<12.2f}'
    Cstates = ['Free', 'Alb Bound', 'CBG Bound', 
            'Free %', 'Alb %', 'CBG %']
    print('\nPlasma\n======')
    print(strFormat2.format(*Cstates))
    print(strFormat2N.format(*bind[0]))
    print('\nInterstitial\n============')
    print(strFormat2.format(*Cstates))
    print(strFormat2N.format(*bind[1]))
   
#print_ss_report(sol, ps, var_order, Male_Master)


# Save the steady state
c = calc_conc(sol, ps)
Male_steady = {
               'Mass': {key: sol[i] for i, key in enumerate(var_order)},
               'Concentration' : {key: c[i] for i, key in 
                   enumerate(var_order)}
              }

with open('MaleSteady.json', 'w') as f:
    json.dump(Male_steady, f)


# Test CRH Bolus
"""
CRH_ICS = [7.66] #2, 50, 500, 1000] # pg/ml

# Set Initial Conditions
i0Sets = [DEM.copy.copy(sol) for i in range(len(CRH_ICS))]
for i, s in enumerate(i0Sets):
    i0Sets[i][0] = CRH_ICS[i] * ps['v_p']

# Simulate the System
sols = DEM.np.zeros((3000, 24, len(CRH_ICS)))
tspan = [0,300]
evals = DEM.np.linspace(0,300,3000)
for j, i in enumerate(i0Sets):
    sol = DEM.solve_ivp(sys, tspan, i, method='LSODA', t_eval=evals,
            args=(ps,))
    sols[:,:,j] = sol.y.T
"""
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

"""
labels1 = [str(float(i))  for i in CRH_ICS]
title1 = 'CRH IC'
ConcSols = calc_conc(sols, ps)
male = plot_sols(evals, ConcSols, tuple(labels1), title1)
male[0].text(0.5, 0.04, 'Time (min)', ha='center')

male[0].savefig('hpaMaleCRHStim.png',
        bbox_extra_artists=tuple(male[1:]), bbox_inches='tight')
"""

def toCortT_ugdl(x, v):
    cort = x[:,2,:] / (v['MWCort'] * v['v_p'])# to moles
    albcort = x[:,15,:] / (v['MWAC'] * v['v_p'])# to moles
    cbgcort = x[:,-1,:] / (v['MWCC'] * v['v_p'])# to moles
    # Total cortisol converted to ug/dl
    totalCort = (cort + albcort + cbgcort) * 1000

    return totalCort


# ACTH Bolus

ACTH_ICS = [141.4, 352.6] #2, 50, 500, 1000] # pg/ml

# Set Initial Conditions
i0Sets = [DEM.copy.copy(sol) for i in range(len(ACTH_ICS))]
for i, s in enumerate(i0Sets):
    i0Sets[i][1] = ACTH_ICS[i] * ps['v_p']

# Simulate the System
sols = DEM.np.zeros((600, 24, len(ACTH_ICS)))
tspan = [0,60]
evals = DEM.np.linspace(0,60,600)
for j, i in enumerate(i0Sets):
    sol = DEM.solve_ivp(sys, tspan, i, method='LSODA', t_eval=evals,
            args=(ps,))
    sols[:,:,j] = sol.y.T


labels1 = [str(float(i))  for i in ACTH_ICS]
title1 = 'ACTH IC'
ConcSols = calc_conc(sols, ps)
ConcSols[:,6,:] = toCortT_ugdl(sols, ps)
male = plot_sols(evals, ConcSols, tuple(labels1), title1)
male[0].text(0.5, 0.04, 'Time (min)', ha='center')

male[0].savefig('hpaMaleACTHStim.png',
        bbox_extra_artists=tuple(male[1:]), bbox_inches='tight')

