import DEModel as DEM
import Albumin as ALB
import CBG2 as CBG

##########################################################################
# Solve for unknowns  of Albumin and CBG                                 #
##########################################################################

"""

Solution Params:
'v_p'      : 
'MWCort'   : 
'MWAlb'    : 
'MWAC'     : 
'konC'     : 
'koffC'    :  
'MWCBG'    : 
'MWCC'     : 
'k4'       :  
'Klp'      : 
'sigA'     : 
'Deg_a'    : 
'GFR'      : 
'AlbPerm'  :  
'ModPerm'  : 
'TL_LW_R'  : 
'TM_LW_R'  : 
'TU_LW_R'  : 
'v_tl'     : 
'v_tm'     : 
'v_tu'     : 
'LAPerm'   : 
'MAPerm'   : 
'UAPerm'   : 
'Temp_K'   : 
'IgpCOP'   : 
'HetapCOP' : 
'AGPpCOP'  : 
'half'     : 
'tonano'   : 
'mictonano': 
'COPEffect':


Solution Vars
'Alb_tl'      : 
'Alb-Cort_tl' : 
'Alb_tm'      : 
'Alb-Cort_tm' :                
'Alb_tu'      : 
'Alb-Cort_tu' : 
'Alb_p'       : 
'Alb-Cort_p'  : 
'CBG_tlc'     : 
'CBG-Cort_tlc': 
'CBG_tmc'     : 
'CBG-Cort_tmc':  
'CBG_tuc'     : 
'CBG-Cort_tuc': 
'CBG_pc'      : 
'CBG-Cort_pc' :  
'konA'        : 
'koffA'       : 
'k5'          : 
'sigC'        : 

EQS     = 17
Unkowns = 20 


Extra Constraints:
Alb_p + Alb-Cort_p = 0.046 * v_p
CBG_p + CBG-Cort_p = 35 * v_p
koffa / kona = 3e-5


"""

def main():
    
    # Vars
    Vars = {"v\['Alb_tl'\]": 'x[0]' , "v\['Alb-Cort_tl'\]": 'x[1]' ,
            "v\['Alb_tm'\]": 'x[2]' , "v\['Alb-Cort_tm'\]": 'x[3]' ,
            "v\['Alb_tu'\]": 'x[4]' , "v\['Alb-Cort_tu'\]": 'x[5]' , 
            "v\['Alb_p'\]" : 'x[6]' , "v\['Alb-Cort_p'\]" : 'x[7]' ,
            "v\['konA'\]"  : 'x[8]' , "v\['koffA'\]"      : 'x[9]',
            "v\['k5'\]"    : 'x[10]', "v\['sigC'\]"       : 'x[11]',
            "v\['CBG_tl'\]": 'x[12]', "v\['CBG-Cort_tl'\]": 'x[13]',
            "v\['CBG_tm'\]": 'x[14]', "v\['CBG-Cort_tm'\]": 'x[15]',
            "v\['CBG_tu'\]": 'x[16]', "v\['CBG-Cort_tu'\]": 'x[17]',
            "v\['CBG-Cort_p'\]": 'x[18]'}


    ######################################################################
    # albumin eqs                                                        #
    ######################################################################  

    Albtoconc = DEM.DeEquation(ALB.Albtoconc_body)
    Albtoconc.make_fun('albtoconc', ALB.Albtoconc_head, 
        DEM.ret_str, Vars)

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

    ######################################################################
    # CBG eqs                                                            #
    ###################################################################### 
    
    CBGtoconc = DEM.DeEquation(CBG.CBGtoconc_body)
    CBGtoconc.make_fun('cbgtoconc', CBG.CBGtoconc_head, 
        DEM.ret_str, Vars)

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


    ######################################################################
    # Make sys                                                           #
    ######################################################################

    def sys(t, x, v):
    
        # Calculate concentration for albumin components
        v = Albtoconc.albtoconc(x, v)
    
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

        # Calculate concentration for CBG components
        v = CBGtoconc.cbgtoconc(x, v)

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

        # KDA constraint
        KDA       = x[8] / x[9] - 3e-3
        # Additional Constraint of interstitial CBG 75% of plasma
        InterCBG = (
                    # 66% of plasma    
                    0.66 * (v['CBG_pc'])
                    # is equal to the total 
                    - ((x[12]/v['v_tl'] + (x[14]/v['v_tm'])  
                        + (x[16]/v['v_tu'])))
                   )

        lhs = [albtl, albcorttl, albtm, albcorttm, albtu, albcorttu, 
               cope, albp, albcortp, KDA, cbgp, cbgtl, InterCBG, 
               cbgcorttl, cbgtm, cbgcorttm, cbgtu, cbgcorttu, cbgcortp]
        
        return lhs

    # Parameters
    full_dict = {**DEM.Paramsu, **DEM.stateVarsu, **DEM.stateVarsMassu}
    full_dict['CBG_p'] = 37.675 * full_dict['v_p']
    full_dict['CBG_pc'] = 37.675
    full_dict['Cort_ic'] = 3.06
    full_dict['Cort_pc'] = 3.06
    # Create the solver function
    y = lambda x : sys(0, x, full_dict)
    # create initial Guess 
    guess = DEM.np.zeros(19)
    guess[0]  = full_dict['Alb_tl']
    guess[1]  = full_dict['Alb-Cort_tl']
    guess[2]  = full_dict['Alb_tm']
    guess[3]  = full_dict['Alb-Cort_tm']
    guess[4]  = full_dict['Alb_tu']
    guess[5]  = full_dict['Alb-Cort_tu']
    guess[6]  = full_dict['Alb_p']
    guess[7]  = full_dict['Alb-Cort_p']
    guess[8]  = full_dict['konA']
    guess[9]  = full_dict['koffA']   
    guess[10] = full_dict['k5']
    guess[11] = 0.65
    guess[12] = 3000
    guess[13] = 1000
    guess[14] = 1000
    guess[15] = 1000
    guess[16] = 1000
    guess[17] = 1000
    guess[18] = 1000


    # Solve
    steady, info, ier, mesg  = DEM.fsolve(y, guess, full_output=True)
    if ier == 1:
        print('The Steady State solution is:')
        print(steady)
        print('Checking the Sys: ')
        print(sys(0, steady, full_dict))
    else:
        print('Solution was not found')
        print('Error:', mesg)


    return 0

if __name__ == '__main__':
    main()
