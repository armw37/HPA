import DEModel as DEM

"""
Working well for no binding interpreation of the sigma c is weird

"""

##########################################################################
# plasma CBG and CBG bound Cortisol functions                            #
##########################################################################

# Convert to concentration
CBGtoconc_head = "def cbgtoconc(self, x, v):\n\t"
CBGtoconc_body = (
                 "v['CBG_tlc']      = v['CBG_tl'] / v['v_tl']; "
                 "v['CBG_tmc']      = v['CBG_tm'] / v['v_tm']; "
                 "v['CBG_tuc']      = v['CBG_tu'] / v['v_tu']; "
                 "v['CBG_pc']       = v['CBG_p'] / v['v_p']; "
                 "v['CBG-Cort_tlc'] = v['CBG-Cort_tl'] / v['v_tl']; "
                 "v['CBG-Cort_tmc'] = v['CBG-Cort_tm'] / v['v_tm']; "
                 "v['CBG-Cort_tuc'] = v['CBG-Cort_tu'] / v['v_tu']; "
                 "v['CBG-Cort_pc']  = v['CBG-Cort_p'] / v['v_p'];"
                 "rhs = v"
                 )


# Lower Lymph CBG
TL_LymphCBGRate_head = "def tl_lymphcbgrate(self, t, x, v):\n\trhs = "
TL_LymphCBGRate_body = (
                        "((v['Klp'] * v['CBG_pc'])"
                        "+ ( (1 - v['Klp']) * v['CBG_tlc'] ))"
                        " * v['TL_LW_R']"
                       )


TL_LymphCCRate_head = "def tl_lymphccrate(self, t, x, v):\n\trhs = "
TL_LymphCCRate_body = (
                        "((v['Klp'] * v['CBG-Cort_pc'])"
                        "+ ( (1 - v['Klp']) * v['CBG-Cort_tlc'] ))"
                        " * v['TL_LW_R']"
                       )


# Middle Lymph CBG
TM_LymphCBGRate_head = "def tm_lymphcbgrate(self, t, x, v):\n\trhs = "
TM_LymphCBGRate_body = (
                        "((v['Klp'] * v['CBG_pc'])"
                        "+ ( (1 - v['Klp']) * v['CBG_tmc'] ))"
                        " * v['TM_LW_R']"
                       )

TM_LymphCCRate_head = "def tm_lymphccrate(self, t, x, v):\n\trhs = "
TM_LymphCCRate_body = (
                        "((v['Klp'] * v['CBG-Cort_pc'])"
                        "+ ( (1 - v['Klp']) * v['CBG-Cort_tmc'] ))"
                        " * v['TM_LW_R']"
                       )



# Upper Lymph CBG
TU_LymphCBGRate_head = "def tu_lymphcbgrate(self, t, x, v):\n\trhs = "
TU_LymphCBGRate_body = (
                        "((v['Klp'] * v['CBG_pc'])"
                        "+ ( (1 - v['Klp']) * v['CBG_tuc'] ))"
                        " * v['TU_LW_R']"
                       )

TU_LymphCCRate_head = "def tu_lymphccrate(self, t, x, v):\n\trhs = "
TU_LymphCCRate_body = (
                        "((v['Klp'] * v['CBG-Cort_pc'])"
                        "+ ( (1 - v['Klp']) * v['CBG-Cort_tuc'] ))"
                        " * v['TU_LW_R']"
                       )


# Lower Cap CBG
TL_CapCBGRate_head = "def tl_capcbgrate(self, t, x, v):\n\trhs = "
TL_CapCBGRate_body = (
                      "v['LAPerm'] * v['sigC'] "
                      " * (v['CBG_pc'] - v['CBG_tlc'])"
                     )


TL_CapCCRate_head = "def tl_capccrate(self, t, x, v):\n\trhs = "
TL_CapCCRate_body = (
                      "v['LAPerm'] * v['sigC'] "
                      " * (v['CBG-Cort_pc'] - v['CBG-Cort_tlc'])"
                     )

# Middle Cap CBG
TM_CapCBGRate_head = "def tm_capcbgrate(self, t, x, v):\n\trhs = "
TM_CapCBGRate_body = (
                      "v['MAPerm'] * v['sigC'] " 
                      "* (v['CBG_pc'] - v['CBG_tmc'])"
                     )


TM_CapCCRate_head = "def tm_capccrate(self, t, x, v):\n\trhs = "
TM_CapCCRate_body = (
                      "v['MAPerm'] * v['sigC'] " 
                      "* (v['CBG-Cort_pc'] - v['CBG-Cort_tmc'])"
                     )

# Upper Cap CBG
TU_CapCBGRate_head = "def tu_capcbgrate(self, t, x, v):\n\trhs = "
TU_CapCBGRate_body = (
                      "v['UAPerm'] * v['sigC']"
                      "* (v['CBG_pc'] - v['CBG_tuc'])"
                       )

TU_CapCCRate_head = "def tu_capccrate(self, t, x, v):\n\trhs = "
TU_CapCCRate_body = (
                      "v['UAPerm'] * v['sigC']"
                      "* (v['CBG-Cort_pc'] - v['CBG-Cort_tuc'])"
                       )





# Torso Lower Interstitial CBG
CBGtl_head = ("def cbg_tl(self, t, x, v, TL_CapCBGRate,"
              "TL_LymphCBGRate):\n\trhs = ")
CBGtl_body = (
              "TL_CapCBGRate - TL_LymphCBGRate"
              # On Binding
              "+ (- v['konC'] * (1/v['MWCort']) * v['CBG_tlc']"  
              "* v['Cort_ic'] " 
              # cbg off binding
              "+ v['koffC'] * v['CBG-Cort_tlc'] * (v['MWCBG']/v['MWCC'])"
              " * (1/v['mictonano'])) * v['v_tl']"
              )


CBGCorttl_head = ("def cbgcort_tl(self, t, x, v, TL_CapCCRate,"
              "TL_LymphCCRate):\n\trhs = ")
CBGCorttl_body = (
                  "TL_CapCCRate - TL_LymphCCRate"
                  # cbg-cort binding
                  "+( v['konC'] * (1/v['MWCBG']) * v['CBG_tlc']"  
                  "* v['Cort_ic'] * (v['MWCC']/v['MWCort'])  "
                  "* v['mictonano']"
                  # off binding
                  "- v['koffC'] * v['CBG-Cort_tlc']) * v['v_tl']"
                 )


# Torso Middle Interstitial Protein CBG
CBGtm_head = ("def cbg_tm(self, t, x, v, TM_CapCBGRate,"
              "TM_LymphCBGRate):\n\trhs = ")
CBGtm_body = (
              "TM_CapCBGRate - TM_LymphCBGRate"
              # On Binding
              " + (- v['konC'] * (1/v['MWCort']) * v['CBG_tmc']"  
              "* v['Cort_ic'] " 
              # cbgumin off binding
              "+ v['koffC'] * v['CBG-Cort_tmc'] * (v['MWCBG']/v['MWCC'])"
              " * (1/v['mictonano'])) * v['v_tm']"
              )


CBGCorttm_head = ("def cbgcort_tm(self, t, x, v, TM_CapCCRate,"
              "TM_LymphCCRate):\n\trhs = ")
CBGCorttm_body = (
                  "TM_CapCCRate - TM_LymphCCRate"
                  # cbg-cort binding
                  "+ (v['konC'] * (1/v['MWCBG']) * v['CBG_tmc']"  
                  "* v['Cort_ic'] * (v['MWCC']/v['MWCort']) "
                  "* v['mictonano']"
                  # off binding
                  "- v['koffC'] * v['CBG-Cort_tmc']) * v['v_tm']"
                  )


# Torso_Upper_InterstitialProtein.[CBG]
CBGtu_head = ("def cbg_tu(self, t, x, v, TU_CapCBGRate,"
              "TU_LymphCBGRate):\n\trhs = ")
CBGtu_body = (
              "TU_CapCBGRate - TU_LymphCBGRate"
              # On Binding
              "+ (- v['konC'] * (1/v['MWCort']) * v['CBG_tuc']"  
              "* v['Cort_ic'] " 
              # cbgumin off binding
              "+ v['koffC'] * v['CBG-Cort_tuc'] * (v['MWCBG']/v['MWCC'])"
              " * (1/v['mictonano'])) * v['v_tu']"
              )


CBGCorttu_head = ("def cbgcort_tu(self, t, x, v, TU_CapCCRate,"
              "TU_LymphCCRate):\n\trhs = ")
CBGCorttu_body = (
                  "TU_CapCCRate - TU_LymphCCRate"
                  # cbg-cort binding
                  "+ (v['konC'] * (1/v['MWCBG']) * v['CBG_tuc']"  
                  "* v['Cort_ic'] * (v['MWCC']/v['MWCort']) "
                  "* v['mictonano']"
                  # off binding
                  "- v['koffC'] * v['CBG-Cort_tuc']) * v['v_tu']"
                  )


# Plasma 
CBGp_head = ("def cbg_p(self, t, x, v, Total_LymphCBGRate, "
        "Total_CapCBGRate):\n\trhs = ")
CBGp_body = (
             "v['k5'] + Total_LymphCBGRate - Total_CapCBGRate "
             # half life rate
             "- (v['half'] * v['CBG_p'])"
             # On Binding
             "+ (- v['konC'] * (1/v['MWCort']) * v['CBG_pc']"  
             "* v['Cort_pc'] " 
             # cbgumin off binding
             "+ v['koffC'] * v['CBG-Cort_pc'] * (v['MWCBG']/v['MWCC'])"
             " * (1/v['mictonano'])) * v['v_p']"
             )


CBGCortp_head = ("def cbgcort_p(self, t, x, v, Total_LymphCCRate, "
                 "Total_CapCCRate):\n\trhs = ")
CBGCortp_body = (
                 " Total_LymphCCRate - Total_CapCCRate"
                 # On Binding
                 "+ (v['konC'] * (1/v['MWCBG']) * v['CBG_pc']"  
                  "* v['Cort_pc'] * (v['MWCC']/v['MWCort']) "
                  "* v['mictonano']"
                 # cbgumin off binding
                 "- v['koffC'] * v['CBG-Cort_pc']) * v['v_p']" 
                 )



def Unit_test(ToConc, CBGtl, CBGCorttl, CBGtm, CBGCorttm, CBGtu,
        CBGCorttu, CBGp, CBGCortp, rates):
    
    # Creat the vector of Cortisone state vars
    testvec = [val for val in DEM.stateVarsMass.values()][19:]
    
    # Calculate concentration and add to params
    ps = {**DEM.Params, **DEM.stateVars, **DEM.stateVarsMass}
    ps = ToConc.cbgtoconc(testvec, ps)
    
    # Calculate the rates 
    TL_LymphCBGRate    = rates[0].tl_lymphcbgrate(0, testvec, ps)
    TL_LymphCCRate     = rates[1].tl_lymphccrate(0, testvec, ps)
    TM_LymphCBGRate    = rates[2].tm_lymphcbgrate(0, testvec, ps)
    TM_LymphCCRate     = rates[3].tm_lymphccrate(0, testvec, ps)
    TU_LymphCBGRate    = rates[4].tu_lymphcbgrate(0, testvec, ps)
    TU_LymphCCRate     = rates[5].tu_lymphccrate(0, testvec, ps)
    TL_CapCBGRate      = rates[6].tl_capcbgrate(0, testvec, ps)
    TL_CapCCRate       = rates[7].tl_capccrate(0, testvec, ps)
    TM_CapCBGRate      = rates[8].tm_capcbgrate(0, testvec, ps)
    TM_CapCCRate       = rates[9].tm_capccrate(0, testvec, ps)
    TU_CapCBGRate      = rates[10].tu_capcbgrate(0, testvec, ps)
    TU_CapCCRate       = rates[11].tu_capccrate(0, testvec, ps)

    # CBG tl
    print("CBGtl:", CBGtl.cbg_tl(0, testvec, ps, 
        TL_CapCBGRate, TL_LymphCBGRate).to(DEM.u.ng/DEM.u.min))
    
    # CBG-Cort tl
    print("CBGCorttl:", CBGCorttl.cbgcort_tl(0, testvec, ps, 
        TL_CapCCRate, TL_LymphCCRate))
    
    # CBG tm
    print("CBGtm:", CBGtm.cbg_tm(0, testvec, ps, 
        TM_CapCBGRate, TM_LymphCBGRate))
    
    # CBG-Cort tm
    print("CBGCorttm:", CBGCorttm.cbgcort_tm(0, testvec, ps, 
        TM_CapCCRate, TM_LymphCCRate))
    
    # CBG tu
    print("CBGtu:", CBGtu.cbg_tu(0, testvec, ps, 
        TU_CapCBGRate, TU_LymphCBGRate))
    
    # CBG-Cort tu
    print("CBGCorttu:", CBGCorttu.cbgcort_tu(0, testvec, ps, 
        TU_CapCCRate, TU_LymphCCRate))
    
    # Calculate cbg totals
    Total_LymphCBGRate = (TL_LymphCBGRate + TM_LymphCBGRate 
                       + TU_LymphCBGRate)

    Total_CapCBGRate = (TL_CapCBGRate + TM_CapCBGRate 
                       + TU_CapCBGRate)
    
    # CBG 
    print("CBGp:", CBGp.cbg_p(0, testvec, ps, Total_LymphCBGRate,
        Total_CapCBGRate))
    
    # Calculate cbg cort totals
    Total_LymphCCRate = (TL_LymphCCRate + TM_LymphCCRate 
                       + TU_LymphCCRate)

    Total_CapCCRate = (TL_CapCCRate + TM_CapCCRate 
                       + TU_CapCCRate)

    # CBG-Cort tm
    print("CBGCortp:", CBGCortp.cbgcort_p(0, testvec, ps, 
        Total_LymphCCRate, Total_CapCCRate))

    return 0


def numerical_test(ToConc, CBGtl, CBGCorttl, CBGtm, CBGCorttm, CBGtu,
        CBGCorttu, CBGp, CBGCortp, rates):
    # Create the vector of Cortisone state vars
    testvec = [val for val in DEM.stateVarsMassu.values()][19:]
    # Calculate concentration and add to params
    ps = {**DEM.Paramsu, **DEM.stateVarsu, **DEM.stateVarsMassu}
    ps = ToConc.cbgtoconc(testvec, ps)

    # Calculate the rates 
    TL_LymphCBGRate    = rates[0].tl_lymphcbgrate(0, testvec, ps)
    TL_LymphCCRate     = rates[1].tl_lymphccrate(0, testvec, ps)
    TM_LymphCBGRate    = rates[2].tm_lymphcbgrate(0, testvec, ps)
    TM_LymphCCRate     = rates[3].tm_lymphccrate(0, testvec, ps)
    TU_LymphCBGRate    = rates[4].tu_lymphcbgrate(0, testvec, ps)
    TU_LymphCCRate     = rates[5].tu_lymphccrate(0, testvec, ps)
    TL_CapCBGRate      = rates[6].tl_capcbgrate(0, testvec, ps)
    TL_CapCCRate       = rates[7].tl_capccrate(0, testvec, ps)
    TM_CapCBGRate      = rates[8].tm_capcbgrate(0, testvec, ps)
    TM_CapCCRate       = rates[9].tm_capccrate(0, testvec, ps)
    TU_CapCBGRate      = rates[10].tu_capcbgrate(0, testvec, ps)
    TU_CapCCRate       = rates[11].tu_capccrate(0, testvec, ps)

    # CBG tl
    print("CBGtl:", CBGtl.cbg_tl(0, testvec, ps, 
        TL_CapCBGRate, TL_LymphCBGRate))
    
    # CBG-Cort tl
    print("CBGCorttl:", CBGCorttl.cbgcort_tl(0, testvec, ps, 
        TL_CapCCRate, TL_LymphCCRate))
    
    # CBG tm
    print("CBGtm:", CBGtm.cbg_tm(0, testvec, ps, 
        TM_CapCBGRate, TM_LymphCBGRate))
    
    # CBG-Cort tm
    print("CBGCorttm:", CBGCorttm.cbgcort_tm(0, testvec, ps, 
        TM_CapCCRate, TM_LymphCCRate))
    
    # CBG tu
    print("CBGtu:", CBGtu.cbg_tu(0, testvec, ps, 
        TU_CapCBGRate, TU_LymphCBGRate))
    
    # CBG-Cort tu
    print("CBGCorttu:", CBGCorttu.cbgcort_tu(0, testvec, ps, 
        TU_CapCCRate, TU_LymphCCRate))
    
    # Calculate cbg totals
    Total_LymphCBGRate = (TL_LymphCBGRate + TM_LymphCBGRate 
                       + TU_LymphCBGRate)

    Total_CapCBGRate = (TL_CapCBGRate + TM_CapCBGRate 
                       + TU_CapCBGRate)
    
    # CBG 
    print("CBGp:", CBGp.cbg_p(0, testvec, ps, Total_LymphCBGRate,
        Total_CapCBGRate))
    
    # Calculate cbg cort totals
    Total_LymphCCRate = (TL_LymphCCRate + TM_LymphCCRate 
                       + TU_LymphCCRate)

    Total_CapCCRate = (TL_CapCCRate + TM_CapCCRate 
                       + TU_CapCCRate)

    # CBG-Cort p
    print("CBGCortp:", CBGCortp.cbgcort_p(0, testvec, ps, 
        Total_LymphCCRate, Total_CapCCRate))    
    
    return 0


##########################################################################
# main                                                                   #
##########################################################################

def main():

    """ 
    # Test the solution
    Vars = {
            "v\['CBG_tl'\]": 'x[0]', "v\['CBG-Cort_tl'\]": 'x[1]',
            "v\['CBG_tm'\]": 'x[2]', "v\['CBG-Cort_tm'\]": 'x[3]',
            "v\['CBG_tu'\]": 'x[4]', "v\['CBG-Cort_tu'\]": 'x[5]',
            "v\['CBG_p'\]" : 'x[6]', "v\['CBG-Cort_p'\]" : 'x[7]',
           }


    CBGtoconc = DEM.DeEquation(CBGtoconc_body)
    CBGtoconc.make_fun('cbgtoconc', CBGtoconc_head, 
        DEM.ret_str, Vars)

    # Torso Lower Lymph 
    TL_LymphCBGRate = DEM.DeEquation(TL_LymphCBGRate_body)
    TL_LymphCBGRate.make_fun('tl_lymphcbgrate', 
        TL_LymphCBGRate_head, DEM.ret_str, Vars)

    TL_LymphCCRate = DEM.DeEquation(TL_LymphCCRate_body)
    TL_LymphCCRate.make_fun('tl_lymphccrate', 
        TL_LymphCCRate_head, DEM.ret_str, Vars)

    # Torso Middle Lymph
    TM_LymphCBGRate = DEM.DeEquation(TM_LymphCBGRate_body)
    TM_LymphCBGRate.make_fun('tm_lymphcbgrate', 
        TM_LymphCBGRate_head, DEM.ret_str, Vars)

    TM_LymphCCRate = DEM.DeEquation(TM_LymphCCRate_body)
    TM_LymphCCRate.make_fun('tm_lymphccrate', 
        TM_LymphCCRate_head, DEM.ret_str, Vars)

    # Torso Lower Lymph
    TU_LymphCBGRate = DEM.DeEquation(TU_LymphCBGRate_body)
    TU_LymphCBGRate.make_fun('tu_lymphcbgrate', 
        TU_LymphCBGRate_head, DEM.ret_str, Vars)

    TU_LymphCCRate = DEM.DeEquation(TU_LymphCCRate_body)
    TU_LymphCCRate.make_fun('tu_lymphccrate', 
        TU_LymphCCRate_head, DEM.ret_str, Vars)


    # Torso Lower Cap
    TL_CapCBGRate = DEM.DeEquation(TL_CapCBGRate_body)
    TL_CapCBGRate.make_fun('tl_capcbgrate', 
        TL_CapCBGRate_head, DEM.ret_str, Vars)

    TL_CapCCRate = DEM.DeEquation(TL_CapCCRate_body)
    TL_CapCCRate.make_fun('tl_capccrate', 
        TL_CapCCRate_head, DEM.ret_str, Vars)

    # Torso Middle Cap
    TM_CapCBGRate = DEM.DeEquation(TM_CapCBGRate_body)
    TM_CapCBGRate.make_fun('tm_capcbgrate', 
        TM_CapCBGRate_head, DEM.ret_str, Vars)

    TM_CapCCRate = DEM.DeEquation(TM_CapCCRate_body)
    TM_CapCCRate.make_fun('tm_capccrate', 
        TM_CapCCRate_head, DEM.ret_str, Vars)

    # Torso Upper Cap
    TU_CapCBGRate = DEM.DeEquation(TU_CapCBGRate_body)
    TU_CapCBGRate.make_fun('tu_capcbgrate', 
        TU_CapCBGRate_head, DEM.ret_str, Vars)

    TU_CapCCRate = DEM.DeEquation(TU_CapCCRate_body)
    TU_CapCCRate.make_fun('tu_capccrate', 
        TU_CapCCRate_head, DEM.ret_str, Vars)


    # Torso Lower Balance
    CBGtl = DEM.DeEquation(CBGtl_body)
    CBGtl.make_fun('cbg_tl', CBGtl_head, DEM.ret_str, Vars)

    CBGCorttl = DEM.DeEquation(CBGCorttl_body)
    CBGCorttl.make_fun('cbgcort_tl', 
        CBGCorttl_head, DEM.ret_str, Vars)

    # Torso Middle
    CBGtm = DEM.DeEquation(CBGtm_body)
    CBGtm.make_fun('cbg_tm', 
        CBGtm_head, DEM.ret_str, Vars)

    CBGCorttm = DEM.DeEquation(CBGCorttm_body)
    CBGCorttm.make_fun('cbgcort_tm', 
        CBGCorttm_head, DEM.ret_str, Vars)

    # Torso Upper
    CBGtu = DEM.DeEquation(CBGtu_body)
    CBGtu.make_fun('cbg_tu', 
        CBGtu_head, DEM.ret_str, Vars)

    CBGCorttu = DEM.DeEquation(CBGCorttu_body)
    CBGCorttu.make_fun('cbgcort_tu', 
        CBGCorttu_head, DEM.ret_str, Vars)

    # Plasma
    CBGp = DEM.DeEquation(CBGp_body)
    CBGp.make_fun('cbg_p', 
        CBGp_head, DEM.ret_str, Vars)

    CBGCortp = DEM.DeEquation(CBGCortp_body)
    CBGCortp.make_fun('cbgcort_p', 
        CBGCortp_head, DEM.ret_str, Vars)
    
    
    rates = [TL_LymphCBGRate, TL_LymphCCRate, TM_LymphCBGRate, 
            TM_LymphCCRate, TU_LymphCBGRate, TU_LymphCCRate, 
            TL_CapCBGRate, TL_CapCCRate, TM_CapCBGRate, TM_CapCCRate, 
            TU_CapCBGRate, TU_CapCCRate]
    
    
    # Test    
    Unit_test(CBGtoconc, CBGtl, CBGCorttl, CBGtm, CBGCorttm, CBGtu,
        CBGCorttu, CBGp, CBGCortp, rates)
    numerical_test(CBGtoconc, CBGtl, CBGCorttl, CBGtm, CBGCorttm, CBGtu,
        CBGCorttu, CBGp, CBGCortp, rates)
    """
    

    Vars = {"v\['k5'\]"    : 'x[0]', "v\['sigC'\]"       : 'x[1]',
            "v\['CBG_tl'\]": 'x[2]', "v\['CBG-Cort_tl'\]": 'x[3]',
            "v\['CBG_tm'\]": 'x[4]', "v\['CBG-Cort_tm'\]": 'x[5]',
            "v\['CBG_tu'\]": 'x[6]', "v\['CBG-Cort_tu'\]": 'x[7]',
            "v\['CBG-Cort_p'\]": 'x[8]'}

    ######################################################################
    # CBG eqs                                                            #
    ###################################################################### 
    
    CBGtoconc = DEM.DeEquation(CBGtoconc_body)
    CBGtoconc.make_fun('cbgtoconc', CBGtoconc_head, 
        DEM.ret_str, Vars)

    # Torso Lower Lymph 
    TL_LymphCBGRate = DEM.DeEquation(TL_LymphCBGRate_body)
    TL_LymphCBGRate.make_fun('tl_lymphcbgrate', 
        TL_LymphCBGRate_head, DEM.ret_str, Vars)

    TL_LymphCCRate = DEM.DeEquation(TL_LymphCCRate_body)
    TL_LymphCCRate.make_fun('tl_lymphccrate', 
        TL_LymphCCRate_head, DEM.ret_str, Vars)

    # Torso Middle Lymph
    TM_LymphCBGRate = DEM.DeEquation(TM_LymphCBGRate_body)
    TM_LymphCBGRate.make_fun('tm_lymphcbgrate', 
        TM_LymphCBGRate_head, DEM.ret_str, Vars)

    TM_LymphCCRate = DEM.DeEquation(TM_LymphCCRate_body)
    TM_LymphCCRate.make_fun('tm_lymphccrate', 
        TM_LymphCCRate_head, DEM.ret_str, Vars)

    # Torso Lower Lymph
    TU_LymphCBGRate = DEM.DeEquation(TU_LymphCBGRate_body)
    TU_LymphCBGRate.make_fun('tu_lymphcbgrate', 
        TU_LymphCBGRate_head, DEM.ret_str, Vars)

    TU_LymphCCRate = DEM.DeEquation(TU_LymphCCRate_body)
    TU_LymphCCRate.make_fun('tu_lymphccrate', 
        TU_LymphCCRate_head, DEM.ret_str, Vars)


    # Torso Lower Cap
    TL_CapCBGRate = DEM.DeEquation(TL_CapCBGRate_body)
    TL_CapCBGRate.make_fun('tl_capcbgrate', 
        TL_CapCBGRate_head, DEM.ret_str, Vars)

    TL_CapCCRate = DEM.DeEquation(TL_CapCCRate_body)
    TL_CapCCRate.make_fun('tl_capccrate', 
        TL_CapCCRate_head, DEM.ret_str, Vars)

    # Torso Middle Cap
    TM_CapCBGRate = DEM.DeEquation(TM_CapCBGRate_body)
    TM_CapCBGRate.make_fun('tm_capcbgrate', 
        TM_CapCBGRate_head, DEM.ret_str, Vars)

    TM_CapCCRate = DEM.DeEquation(TM_CapCCRate_body)
    TM_CapCCRate.make_fun('tm_capccrate', 
        TM_CapCCRate_head, DEM.ret_str, Vars)

    # Torso Upper Cap
    TU_CapCBGRate = DEM.DeEquation(TU_CapCBGRate_body)
    TU_CapCBGRate.make_fun('tu_capcbgrate', 
        TU_CapCBGRate_head, DEM.ret_str, Vars)

    TU_CapCCRate = DEM.DeEquation(TU_CapCCRate_body)
    TU_CapCCRate.make_fun('tu_capccrate', 
        TU_CapCCRate_head, DEM.ret_str, Vars)


    # Torso Lower Balance
    CBGtl = DEM.DeEquation(CBGtl_body)
    CBGtl.make_fun('cbg_tl', CBGtl_head, DEM.ret_str, Vars)

    CBGCorttl = DEM.DeEquation(CBGCorttl_body)
    CBGCorttl.make_fun('cbgcort_tl', 
        CBGCorttl_head, DEM.ret_str, Vars)

    # Torso Middle
    CBGtm = DEM.DeEquation(CBGtm_body)
    CBGtm.make_fun('cbg_tm', 
        CBGtm_head, DEM.ret_str, Vars)

    CBGCorttm = DEM.DeEquation(CBGCorttm_body)
    CBGCorttm.make_fun('cbgcort_tm', 
        CBGCorttm_head, DEM.ret_str, Vars)

    # Torso Upper
    CBGtu = DEM.DeEquation(CBGtu_body)
    CBGtu.make_fun('cbg_tu', 
        CBGtu_head, DEM.ret_str, Vars)

    CBGCorttu = DEM.DeEquation(CBGCorttu_body)
    CBGCorttu.make_fun('cbgcort_tu', 
        CBGCorttu_head, DEM.ret_str, Vars)

    # Plasma
    CBGp = DEM.DeEquation(CBGp_body)
    CBGp.make_fun('cbg_p', 
        CBGp_head, DEM.ret_str, Vars)

    CBGCortp = DEM.DeEquation(CBGCortp_body)
    CBGCortp.make_fun('cbgcort_p', 
        CBGCortp_head, DEM.ret_str, Vars)
    
    
    CBGrates = [TL_LymphCBGRate, TL_LymphCCRate, TM_LymphCBGRate, 
            TM_LymphCCRate, TU_LymphCBGRate, TU_LymphCCRate, 
            TL_CapCBGRate, TL_CapCCRate, TM_CapCBGRate, TM_CapCCRate, 
            TU_CapCBGRate, TU_CapCCRate]


    ######################################################################
    # Make sys                                                           #
    ######################################################################

    def sys(t, x, v):
    
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

        # Additional Constraint of interstitial CBG 75% of plasma
        InterCBG = (
                    # 75% of plasma    
                    0.75 * (v['CBG_pc'])
                    # is equal to the total 
                    - ((x[2]/v['v_tl'] + (x[4]/v['v_tm'])  
                        + (x[6]/v['v_tu'])))
                   )

        lhs = [cbgp, cbgtl, InterCBG, cbgcorttl, cbgtm, 
                cbgcorttm, cbgtu, cbgcorttu, cbgcortp]
        
        return lhs

    # Parameters
    full_dict = DEM.Paramsu
    full_dict['CBG_p'] = 37.675 * full_dict['v_p']
    full_dict['CBG_pc'] = 37.675
    full_dict['Cort_ic'] = 3.06
    full_dict['Cort_pc'] = 3.06
    # Create the solver function
    y = lambda x : sys(0, x, full_dict)
    # create initial Guess 
    guess = DEM.np.zeros(9)
    guess[0] = full_dict['k5']
    guess[1] = 0.65
    guess[2] = 3000
    guess[3] = 1000
    guess[4] = 1000
    guess[5] = 1000
    guess[6] = 1000
    guess[7] = 1000
    guess[8] = 1000

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


