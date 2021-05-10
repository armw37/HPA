import DEModel as DEM

##########################################################################
# Albumin                                                                #
##########################################################################

# albumin binding
#"- v['konA'] * (1/v['MWCort']) * v['Alb_pc']"  
#"* v['Cort_p'] " 
# albumin off binding
#"+ v['koffA'] * v['Alb-Cort_p'] * (v['MWAlb']/v['MWAC'])"
#" * (1/v['tonano'])"


# alb-cort binding
#"+ v['konA'] * (1/v['MWAlb']) * v['Alb_pc']"  
#"* v['Cort_p'] * (1/v['MWCort']) * (v['MWAC']) "
#"* v['tonano']"
# off binding
#"- v['koffA'] * v['Alb-Cort_p']"
  


# Convert to concentration
Albtoconc_head = "def albtoconc(self, x, v):\n\t"
Albtoconc_body = (
                 "v['Alb_tlc']      = v['Alb_tl'] / v['v_tl']; "
                 "v['Alb_tmc']      = v['Alb_tm'] / v['v_tm']; "
                 "v['Alb_tuc']      = v['Alb_tu'] / v['v_tu']; "
                 "v['Alb_pc']       = v['Alb_p'] / v['v_p']; "
                 "v['Alb-Cort_tlc'] = v['Alb-Cort_tl'] / v['v_tl']; "
                 "v['Alb-Cort_tmc'] = v['Alb-Cort_tm'] / v['v_tm']; "
                 "v['Alb-Cort_tuc'] = v['Alb-Cort_tu'] / v['v_tu']; "
                 "v['Alb-Cort_pc']  = v['Alb-Cort_p'] / v['v_p'];"
                 "rhs = v"
                 )


# Lower Lymph Alb
TL_LymphAlbRate_head = "def tl_lymphalbrate(self, t, x, v):\n\trhs = "
TL_LymphAlbRate_body = (
                        "((v['Klp'] * v['Alb_pc'])"
                        "+ ( (1 - v['Klp']) * v['Alb_tlc'] ))"
                        " * v['TL_LW_R']"
                       )


TL_LymphACRate_head = "def tl_lymphacrate(self, t, x, v):\n\trhs = "
TL_LymphACRate_body = (
                        "((v['Klp'] * v['Alb-Cort_pc'])"
                        "+ ( (1 - v['Klp']) * v['Alb-Cort_tlc'] ))"
                        " * v['TL_LW_R']"
                       )


# Middle Lymph Alb
TM_LymphAlbRate_head = "def tm_lymphalbrate(self, t, x, v):\n\trhs = "
TM_LymphAlbRate_body = (
                        "((v['Klp'] * v['Alb_pc'])"
                        "+ ( (1 - v['Klp']) * v['Alb_tmc'] ))"
                        " * v['TM_LW_R']"
                       )

TM_LymphACRate_head = "def tm_lymphacrate(self, t, x, v):\n\trhs = "
TM_LymphACRate_body = (
                        "((v['Klp'] * v['Alb-Cort_pc'])"
                        "+ ( (1 - v['Klp']) * v['Alb-Cort_tmc'] ))"
                        " * v['TM_LW_R']"
                       )



# Upper Lymph Alb
TU_LymphAlbRate_head = "def tu_lymphalbrate(self, t, x, v):\n\trhs = "
TU_LymphAlbRate_body = (
                        "((v['Klp'] * v['Alb_pc'])"
                        "+ ( (1 - v['Klp']) * v['Alb_tuc'] ))"
                        " * v['TU_LW_R']"
                       )

TU_LymphACRate_head = "def tu_lymphacrate(self, t, x, v):\n\trhs = "
TU_LymphACRate_body = (
                        "((v['Klp'] * v['Alb-Cort_pc'])"
                        "+ ( (1 - v['Klp']) * v['Alb-Cort_tuc'] ))"
                        " * v['TU_LW_R']"
                       )


# Lower Cap Alb
TL_CapAlbRate_head = "def tl_capalbrate(self, t, x, v):\n\trhs = "
TL_CapAlbRate_body = (
                      "v['LAPerm'] * v['sigA'] "
                      " * (v['Alb_pc'] - v['Alb_tlc'])"
                     )


TL_CapACRate_head = "def tl_capacrate(self, t, x, v):\n\trhs = "
TL_CapACRate_body = (
                      "v['LAPerm'] * v['sigA'] "
                      " * (v['Alb-Cort_pc'] - v['Alb-Cort_tlc'])"
                     )

# Middle Cap Alb
TM_CapAlbRate_head = "def tm_capalbrate(self, t, x, v):\n\trhs = "
TM_CapAlbRate_body = (
                      "v['MAPerm'] * v['sigA'] " 
                      "* (v['Alb_pc'] - v['Alb_tmc'])"
                     )


TM_CapACRate_head = "def tm_capacrate(self, t, x, v):\n\trhs = "
TM_CapACRate_body = (
                      "v['MAPerm'] * v['sigA'] " 
                      "* (v['Alb-Cort_pc'] - v['Alb-Cort_tmc'])"
                     )

# Upper Cap Alb
TU_CapAlbRate_head = "def tu_capalbrate(self, t, x, v):\n\trhs = "
TU_CapAlbRate_body = (
                      "v['UAPerm'] * v['sigA']"
                      "* (v['Alb_pc'] - v['Alb_tuc'])"
                       )

TU_CapACRate_head = "def tu_capacrate(self, t, x, v):\n\trhs = "
TU_CapACRate_body = (
                      "v['UAPerm'] * v['sigA']"
                      "* (v['Alb-Cort_pc'] - v['Alb-Cort_tuc'])"
                       )


# Albumin partial colliod osmotic
Albpcop_head = "def albpcop(self, x, v):\n\trhs = "
Albpcop_body = (
                 "((v['Alb_pc']/v['MWAlb'])"
                 " + (v['Alb-Cort_pc']/v['MWAC'])*1/v['tonano'])"
                 " * 1000 * 62.36367 * 2.1 * v['Temp_K'] * v['sigA']"
                )



# colloid osmotic pressure
COP_head = "def cop(self, x, v, Albpcop):\n\trhs = "
COP_body = (
            "Albpcop + v['IgpCOP'] + v['HetapCOP']"
            " + v['AGPpCOP']"
            )



# COPeffect 
COPE_head = "def cope(self, t, x, v, copeCurve, COP):\n\trhs = "
COPE_body = (
             "0.015 * (copeCurve(COP) - v['COPEffect'])"
            )


# set up the cop effect curve
copefCurve = DEM.cs([20, 28, 40],[1.15, 1, 0],[0, -0.07, 0])

def copeCurve(cop):
    if cop <= 20.0:
        return 1.15
    elif cop >= 40.0:
        return 0
    else:
        return copefCurve.__call__(cop)



# Torso Lower Interstitial Alb
Albtl_head = ("def alb_tl(self, t, x, v, TL_CapAlbRate,"
              "TL_LymphAlbRate):\n\trhs = ")
Albtl_body = (
              "TL_CapAlbRate - TL_LymphAlbRate"
              # On Binding
              "+ (- v['konA'] * (1/v['MWCort']) * v['Alb_tlc']"  
              "* v['Cort_ic'] " 
              # albumin off binding
              "+ v['koffA'] * v['Alb-Cort_tlc'] * (v['MWAlb']/v['MWAC'])"
              " * (1/v['tonano'])) * v['v_tl']"
              )


AlbCorttl_head = ("def albcort_tl(self, t, x, v, TL_CapACRate,"
              "TL_LymphACRate):\n\trhs = ")
AlbCorttl_body = (
                  "TL_CapACRate - TL_LymphACRate"
                  # alb-cort binding
                  "+( v['konA'] * (1/v['MWAlb']) * v['Alb_tlc']"  
                  "* v['Cort_ic'] * (v['MWAC']/v['MWCort'])  "
                  "* v['tonano']"
                  # off binding
                  "- v['koffA'] * v['Alb-Cort_tlc']) * v['v_tl']"
                 )


# Torso Middle Interstitial Protein Alb
Albtm_head = ("def alb_tm(self, t, x, v, TM_CapAlbRate,"
              "TM_LymphAlbRate):\n\trhs = ")
Albtm_body = (
              "TM_CapAlbRate - TM_LymphAlbRate"
              # On Binding
              " + (- v['konA'] * (1/v['MWCort']) * v['Alb_tmc']"  
              "* v['Cort_ic'] " 
              # albumin off binding
              "+ v['koffA'] * v['Alb-Cort_tmc'] * (v['MWAlb']/v['MWAC'])"
              " * (1/v['tonano'])) * v['v_tm']"
              )


AlbCorttm_head = ("def albcort_tm(self, t, x, v, TM_CapACRate,"
              "TM_LymphACRate):\n\trhs = ")
AlbCorttm_body = (
                  "TM_CapACRate - TM_LymphACRate"
                  # alb-cort binding
                  "+ (v['konA'] * (1/v['MWAlb']) * v['Alb_tmc']"  
                  "* v['Cort_ic'] * (v['MWAC']/v['MWCort']) "
                  "* v['tonano']"
                  # off binding
                  "- v['koffA'] * v['Alb-Cort_tmc']) * v['v_tm']"
                  )


# Torso_Upper_InterstitialProtein.[Alb]
Albtu_head = ("def alb_tu(self, t, x, v, TU_CapAlbRate,"
              "TU_LymphAlbRate):\n\trhs = ")
Albtu_body = (
              "TU_CapAlbRate - TU_LymphAlbRate"
              # On Binding
              "+ (- v['konA'] * (1/v['MWCort']) * v['Alb_tuc']"  
              "* v['Cort_ic'] " 
              # albumin off binding
              "+ v['koffA'] * v['Alb-Cort_tuc'] * (v['MWAlb']/v['MWAC'])"
              " * (1/v['tonano'])) * v['v_tu']"
              )


AlbCorttu_head = ("def albcort_tu(self, t, x, v, TU_CapACRate,"
              "TU_LymphACRate):\n\trhs = ")
AlbCorttu_body = (
                  "TU_CapACRate - TU_LymphACRate"
                  # alb-cort binding
                  "+ (v['konA'] * (1/v['MWAlb']) * v['Alb_tuc']"  
                  "* v['Cort_ic'] * (v['MWAC']/v['MWCort']) "
                  "* v['tonano']"
                  # off binding
                  "- v['koffA'] * v['Alb-Cort_tuc']) * v['v_tu']"
                  )


# Plasma 
Albp_head = ("def alb_p(self, t, x, v, Total_LymphAlbRate, "
        "Total_CapAlbRate):\n\trhs = ")
Albp_body = (
             "v['k4'] * v['COPEffect']" 
             " + Total_LymphAlbRate - v['Deg_a'] "
             "- Total_CapAlbRate "
             "- (v['Alb_pc']*(v['ModPerm']+v['AlbPerm'])*v['GFR']/1000)"
             # On Binding
             "+ (- v['konA'] * (1/v['MWCort']) * v['Alb_pc']"  
             "* v['Cort_pc'] " 
             # albumin off binding
             "+ v['koffA'] * v['Alb-Cort_pc'] * (v['MWAlb']/v['MWAC'])"
             " * (1/v['tonano'])) * v['v_p']"
             )


AlbCortp_head = ("def albcort_p(self, t, x, v, Total_LymphACRate, "
                 "Total_CapACRate):\n\trhs = ")
AlbCortp_body = (
                 " Total_LymphACRate - Total_CapACRate"
                 # Urine
                 "- (v['Alb-Cort_pc']*(v['ModPerm']+v['AlbPerm'])"
                 "*v['GFR']/1000)"
                 # On Binding
                 "+ (v['konA'] * (1/v['MWAlb']) * v['Alb_pc']"  
                  "* v['Cort_pc'] * (v['MWAC']/v['MWCort']) "
                  "* v['tonano']"
                 # albumin off binding
                 "- v['koffA'] * v['Alb-Cort_pc']) * v['v_p']"
                 
                )



def Unit_test(ToConc, Albtl, AlbCorttl, Albtm, AlbCorttm, Albtu,
        AlbCorttu, COPE, Albp, AlbCortp, rates):
    
    # Creat the vector of Cortisone state vars
    testvec = [val for val in DEM.stateVarsMass.values()][10:19]
    
    # Calculate concentration and add to params
    ps = {**DEM.Params, **DEM.stateVars, **DEM.stateVarsMass}
    ps = ToConc.albtoconc(testvec, ps)
    
    # Calculate the rates 
    TL_LymphAlbRate    = rates[0].tl_lymphalbrate(0, testvec, ps)
    TL_LymphACRate     = rates[1].tl_lymphacrate(0, testvec, ps)
    TM_LymphAlbRate    = rates[2].tm_lymphalbrate(0, testvec, ps)
    TM_LymphACRate     = rates[3].tm_lymphacrate(0, testvec, ps)
    TU_LymphAlbRate    = rates[4].tu_lymphalbrate(0, testvec, ps)
    TU_LymphACRate     = rates[5].tu_lymphacrate(0, testvec, ps)
    TL_CapAlbRate      = rates[6].tl_capalbrate(0, testvec, ps)
    TL_CapACRate       = rates[7].tl_capacrate(0, testvec, ps)
    TM_CapAlbRate      = rates[8].tm_capalbrate(0, testvec, ps)
    TM_CapACRate       = rates[9].tm_capacrate(0, testvec, ps)
    TU_CapAlbRate      = rates[10].tu_capalbrate(0, testvec, ps)
    TU_CapACRate       = rates[11].tu_capacrate(0, testvec, ps)
    Albpcop            = rates[12].albpcop(testvec, ps)
    COP                = rates[13].cop(testvec, ps, Albpcop)

    # Alb tl
    print("Albtl:", Albtl.alb_tl(0, testvec, ps, 
        TL_CapAlbRate, TL_LymphAlbRate).to(DEM.u.ng/DEM.u.min))
    
    # Alb-Cort tl
    print("AlbCorttl:", AlbCorttl.albcort_tl(0, testvec, ps, 
        TL_CapACRate, TL_LymphACRate))
    
    # Alb tm
    print("Albtm:", Albtm.alb_tm(0, testvec, ps, 
        TM_CapAlbRate, TM_LymphAlbRate))
    
    # Alb-Cort tm
    print("AlbCorttm:", AlbCorttm.albcort_tm(0, testvec, ps, 
        TM_CapACRate, TM_LymphACRate))
    
    # Alb tu
    print("Albtu:", Albtu.alb_tu(0, testvec, ps, 
        TU_CapAlbRate, TU_LymphAlbRate))
    
    # Alb-Cort tu
    print("AlbCorttu:", AlbCorttu.albcort_tu(0, testvec, ps, 
        TU_CapACRate, TU_LymphACRate))
    
    # COPE
    print('COP Effect:', COPE.cope(0, testvec, ps, 
        copeCurve, COP.magnitude))
    
    # Calculate alb totals
    Total_LymphAlbRate = (TL_LymphAlbRate + TM_LymphAlbRate 
                       + TU_LymphAlbRate)

    Total_CapAlbRate = (TL_CapAlbRate + TM_CapAlbRate 
                       + TU_CapAlbRate)
    
    # Alb 
    print("Albp:", Albp.alb_p(0, testvec, ps, Total_LymphAlbRate,
        Total_CapAlbRate))
    
    # Calculate alb cort totals
    Total_LymphACRate = (TL_LymphACRate + TM_LymphACRate 
                       + TU_LymphACRate)

    Total_CapACRate = (TL_CapACRate + TM_CapACRate 
                       + TU_CapACRate)

    # Alb-Cort tm
    print("AlbCortp:", AlbCortp.albcort_p(0, testvec, ps, 
        Total_LymphACRate, Total_CapACRate))

    return 0


def numerical_test(ToConc, Albtl, AlbCorttl, Albtm, AlbCorttm, Albtu,
        AlbCorttu, COPE, Albp, AlbCortp, rates):
    # Create the vector of Cortisone state vars
    testvec = [val for val in DEM.stateVarsMassu.values()][10:19]
    # Calculate concentration and add to params
    ps = {**DEM.Paramsu, **DEM.stateVarsu, **DEM.stateVarsMassu}
    ps = ToConc.albtoconc(testvec, ps)

    # Calculate the rates 
    TL_LymphAlbRate    = rates[0].tl_lymphalbrate(0, testvec, ps)
    TL_LymphACRate     = rates[1].tl_lymphacrate(0, testvec, ps)
    TM_LymphAlbRate    = rates[2].tm_lymphalbrate(0, testvec, ps)
    TM_LymphACRate     = rates[3].tm_lymphacrate(0, testvec, ps)
    TU_LymphAlbRate    = rates[4].tu_lymphalbrate(0, testvec, ps)
    TU_LymphACRate     = rates[5].tu_lymphacrate(0, testvec, ps)
    TL_CapAlbRate      = rates[6].tl_capalbrate(0, testvec, ps)
    TL_CapACRate       = rates[7].tl_capacrate(0, testvec, ps)
    TM_CapAlbRate      = rates[8].tm_capalbrate(0, testvec, ps)
    TM_CapACRate       = rates[9].tm_capacrate(0, testvec, ps)
    TU_CapAlbRate      = rates[10].tu_capalbrate(0, testvec, ps)
    TU_CapACRate       = rates[11].tu_capacrate(0, testvec, ps)
    Albpcop            = rates[12].albpcop(testvec, ps)
    COP                = rates[13].cop(testvec, ps, Albpcop)

    # Alb tl
    print("Albtl:", Albtl.alb_tl(0, testvec, ps, 
        TL_CapAlbRate, TL_LymphAlbRate))
    
    # Alb-Cort tl
    print("AlbCorttl:", AlbCorttl.albcort_tl(0, testvec, ps, 
        TL_CapACRate, TL_LymphACRate))
    
    # Alb tm
    print("Albtm:", Albtm.alb_tm(0, testvec, ps, 
        TM_CapAlbRate, TM_LymphAlbRate))
    
    # Alb-Cort tm
    print("AlbCorttm:", AlbCorttm.albcort_tm(0, testvec, ps, 
        TM_CapACRate, TM_LymphACRate))
    
    # Alb tu
    print("Albtu:", Albtu.alb_tu(0, testvec, ps, 
        TU_CapAlbRate, TU_LymphAlbRate))
    
    # Alb-Cort tu
    print("AlbCorttu:", AlbCorttu.albcort_tu(0, testvec, ps, 
        TU_CapACRate, TU_LymphACRate))
    
    # COPE
    print('COP Effect:', COPE.cope(0, testvec, ps, 
        copeCurve, COP))
    
    # Calculate alb totals
    Total_LymphAlbRate = (TL_LymphAlbRate + TM_LymphAlbRate 
                       + TU_LymphAlbRate)

    Total_CapAlbRate = (TL_CapAlbRate + TM_CapAlbRate 
                       + TU_CapAlbRate)
    
    # Alb 
    print("Albp:", Albp.alb_p(0, testvec, ps, Total_LymphAlbRate,
        Total_CapAlbRate))
    
    # Calculate alb cort totals
    Total_LymphACRate = (TL_LymphACRate + TM_LymphACRate 
                       + TU_LymphACRate)

    Total_CapACRate = (TL_CapACRate + TM_CapACRate 
                       + TU_CapACRate)

    # Alb-Cort p
    print("AlbCortp:", AlbCortp.albcort_p(0, testvec, ps, 
        Total_LymphACRate, Total_CapACRate))    
    
    return 0


##########################################################################
# main                                                                   #
##########################################################################

def main():

    """ 
    # Test the solution
    Vars = {
            "v\['Alb_tl'\]": 'x[0]', "v\['Alb-Cort_tl'\]": 'x[1]',
            "v\['Alb_tm'\]": 'x[2]', "v\['Alb-Cort_tm'\]": 'x[3]',
            "v\['Alb_tu'\]": 'x[4]', "v\['Alb-Cort_tu'\]": 'x[5]',
            "v\['COPEffect'\]" : 'x[6]', "v\['Alb_p'\]" : 'x[7]', 
            "v\['Alb-Cort_p'\]" : 'x[8]',
           }


    Albtoconc = DEM.DeEquation(Albtoconc_body)
    Albtoconc.make_fun('albtoconc', Albtoconc_head, 
        DEM.ret_str, Vars)

    # Torso Lower Lymph 
    TL_LymphAlbRate = DEM.DeEquation(TL_LymphAlbRate_body)
    TL_LymphAlbRate.make_fun('tl_lymphalbrate', 
        TL_LymphAlbRate_head, DEM.ret_str, Vars)

    TL_LymphACRate = DEM.DeEquation(TL_LymphACRate_body)
    TL_LymphACRate.make_fun('tl_lymphacrate', 
        TL_LymphACRate_head, DEM.ret_str, Vars)

    # Torso Middle Lymph
    TM_LymphAlbRate = DEM.DeEquation(TM_LymphAlbRate_body)
    TM_LymphAlbRate.make_fun('tm_lymphalbrate', 
        TM_LymphAlbRate_head, DEM.ret_str, Vars)

    TM_LymphACRate = DEM.DeEquation(TM_LymphACRate_body)
    TM_LymphACRate.make_fun('tm_lymphacrate', 
        TM_LymphACRate_head, DEM.ret_str, Vars)

    # Torso Lower Lymph
    TU_LymphAlbRate = DEM.DeEquation(TU_LymphAlbRate_body)
    TU_LymphAlbRate.make_fun('tu_lymphalbrate', 
        TU_LymphAlbRate_head, DEM.ret_str, Vars)

    TU_LymphACRate = DEM.DeEquation(TU_LymphACRate_body)
    TU_LymphACRate.make_fun('tu_lymphacrate', 
        TU_LymphACRate_head, DEM.ret_str, Vars)


    # Torso Lower Cap
    TL_CapAlbRate = DEM.DeEquation(TL_CapAlbRate_body)
    TL_CapAlbRate.make_fun('tl_capalbrate', 
        TL_CapAlbRate_head, DEM.ret_str, Vars)

    TL_CapACRate = DEM.DeEquation(TL_CapACRate_body)
    TL_CapACRate.make_fun('tl_capacrate', 
        TL_CapACRate_head, DEM.ret_str, Vars)

    # Torso Middle Cap
    TM_CapAlbRate = DEM.DeEquation(TM_CapAlbRate_body)
    TM_CapAlbRate.make_fun('tm_capalbrate', 
        TM_CapAlbRate_head, DEM.ret_str, Vars)

    TM_CapACRate = DEM.DeEquation(TM_CapACRate_body)
    TM_CapACRate.make_fun('tm_capacrate', 
        TM_CapACRate_head, DEM.ret_str, Vars)

    # Torso Upper Cap
    TU_CapAlbRate = DEM.DeEquation(TU_CapAlbRate_body)
    TU_CapAlbRate.make_fun('tu_capalbrate', 
        TU_CapAlbRate_head, DEM.ret_str, Vars)

    TU_CapACRate = DEM.DeEquation(TU_CapACRate_body)
    TU_CapACRate.make_fun('tu_capacrate', 
        TU_CapACRate_head, DEM.ret_str, Vars)

    # Partial pressure
    Albpcop = DEM.DeEquation(Albpcop_body)
    Albpcop.make_fun('albpcop', Albpcop_head, 
        DEM.ret_str, Vars)

    # Colloid Oncotic Pressure
    COP = DEM.DeEquation(COP_body)
    COP.make_fun('cop', COP_head, 
        DEM.ret_str, Vars)

    # COP effect
    COPE = DEM.DeEquation(COPE_body)
    COPE.make_fun('cope', COPE_head, 
        DEM.ret_str, Vars)

    # Torso Lower Balance
    Albtl = DEM.DeEquation(Albtl_body)
    Albtl.make_fun('alb_tl', Albtl_head, DEM.ret_str, Vars)

    AlbCorttl = DEM.DeEquation(AlbCorttl_body)
    AlbCorttl.make_fun('albcort_tl', 
        AlbCorttl_head, DEM.ret_str, Vars)

    # Torso Middle
    Albtm = DEM.DeEquation(Albtm_body)
    Albtm.make_fun('alb_tm', 
        Albtm_head, DEM.ret_str, Vars)

    AlbCorttm = DEM.DeEquation(AlbCorttm_body)
    AlbCorttm.make_fun('albcort_tm', 
        AlbCorttm_head, DEM.ret_str, Vars)

    # Torso Upper
    Albtu = DEM.DeEquation(Albtu_body)
    Albtu.make_fun('alb_tu', 
        Albtu_head, DEM.ret_str, Vars)

    AlbCorttu = DEM.DeEquation(AlbCorttu_body)
    AlbCorttu.make_fun('albcort_tu', 
        AlbCorttu_head, DEM.ret_str, Vars)

    # Plasma
    Albp = DEM.DeEquation(Albp_body)
    Albp.make_fun('alb_p', 
        Albp_head, DEM.ret_str, Vars)

    AlbCortp = DEM.DeEquation(AlbCortp_body)
    AlbCortp.make_fun('albcort_p', 
        AlbCortp_head, DEM.ret_str, Vars)
    
    
    rates = [TL_LymphAlbRate, TL_LymphACRate, TM_LymphAlbRate, 
            TM_LymphACRate, TU_LymphAlbRate, TU_LymphACRate, 
            TL_CapAlbRate, TL_CapACRate, TM_CapAlbRate, TM_CapACRate, 
            TU_CapAlbRate, TU_CapACRate, Albpcop, COP]
    
    
    # Test    
    Unit_test(Albtoconc, Albtl, AlbCorttl, Albtm, AlbCorttm, Albtu,
        AlbCorttu, COPE, Albp, AlbCortp, rates)
    numerical_test(Albtoconc, Albtl, AlbCorttl, Albtm, AlbCorttm, Albtu,
        AlbCorttu, COPE, Albp, AlbCortp, rates)
    """
    # Vars
    Vars = {"v\['Alb_tl'\]": 'x[0]' , "v\['Alb-Cort_tl'\]": 'x[1]' ,
            "v\['Alb_tm'\]": 'x[2]' , "v\['Alb-Cort_tm'\]": 'x[3]' ,
            "v\['Alb_tu'\]": 'x[4]' , "v\['Alb-Cort_tu'\]": 'x[5]' , 
            "v\['Alb_p'\]" : 'x[6]' , "v\['Alb-Cort_p'\]" : 'x[7]' ,
            "v\['konA'\]"  : 'x[8]', "v\['koffA'\]"      : 'x[9]', 
           }


    ######################################################################
    # albumin eqs                                                        #
    ######################################################################  

    Albtoconc = DEM.DeEquation(Albtoconc_body)
    Albtoconc.make_fun('albtoconc', Albtoconc_head, 
        DEM.ret_str, Vars)

    # Torso Lower Lymph 
    TL_LymphAlbRate = DEM.DeEquation(TL_LymphAlbRate_body)
    TL_LymphAlbRate.make_fun('tl_lymphalbrate', 
        TL_LymphAlbRate_head, DEM.ret_str, Vars)

    TL_LymphACRate = DEM.DeEquation(TL_LymphACRate_body)
    TL_LymphACRate.make_fun('tl_lymphacrate', 
        TL_LymphACRate_head, DEM.ret_str, Vars)

    # Torso Middle Lymph
    TM_LymphAlbRate = DEM.DeEquation(TM_LymphAlbRate_body)
    TM_LymphAlbRate.make_fun('tm_lymphalbrate', 
        TM_LymphAlbRate_head, DEM.ret_str, Vars)

    TM_LymphACRate = DEM.DeEquation(TM_LymphACRate_body)
    TM_LymphACRate.make_fun('tm_lymphacrate', 
        TM_LymphACRate_head, DEM.ret_str, Vars)

    # Torso Lower Lymph
    TU_LymphAlbRate = DEM.DeEquation(TU_LymphAlbRate_body)
    TU_LymphAlbRate.make_fun('tu_lymphalbrate', 
        TU_LymphAlbRate_head, DEM.ret_str, Vars)

    TU_LymphACRate = DEM.DeEquation(TU_LymphACRate_body)
    TU_LymphACRate.make_fun('tu_lymphacrate', 
        TU_LymphACRate_head, DEM.ret_str, Vars)


    # Torso Lower Cap
    TL_CapAlbRate = DEM.DeEquation(TL_CapAlbRate_body)
    TL_CapAlbRate.make_fun('tl_capalbrate', 
        TL_CapAlbRate_head, DEM.ret_str, Vars)

    TL_CapACRate = DEM.DeEquation(TL_CapACRate_body)
    TL_CapACRate.make_fun('tl_capacrate', 
        TL_CapACRate_head, DEM.ret_str, Vars)

    # Torso Middle Cap
    TM_CapAlbRate = DEM.DeEquation(TM_CapAlbRate_body)
    TM_CapAlbRate.make_fun('tm_capalbrate', 
        TM_CapAlbRate_head, DEM.ret_str, Vars)

    TM_CapACRate = DEM.DeEquation(TM_CapACRate_body)
    TM_CapACRate.make_fun('tm_capacrate', 
        TM_CapACRate_head, DEM.ret_str, Vars)

    # Torso Upper Cap
    TU_CapAlbRate = DEM.DeEquation(TU_CapAlbRate_body)
    TU_CapAlbRate.make_fun('tu_capalbrate', 
        TU_CapAlbRate_head, DEM.ret_str, Vars)

    TU_CapACRate = DEM.DeEquation(TU_CapACRate_body)
    TU_CapACRate.make_fun('tu_capacrate', 
        TU_CapACRate_head, DEM.ret_str, Vars)

    # Partial pressure
    Albpcop = DEM.DeEquation(Albpcop_body)
    Albpcop.make_fun('albpcop', Albpcop_head, 
        DEM.ret_str, Vars)

    # Colloid Oncotic Pressure
    COP = DEM.DeEquation(COP_body)
    COP.make_fun('cop', COP_head, 
        DEM.ret_str, Vars)

    # COP effect
    COPE = DEM.DeEquation(COPE_body)
    COPE.make_fun('cope', COPE_head, 
        DEM.ret_str, Vars)

    # Torso Lower Balance
    Albtl = DEM.DeEquation(Albtl_body)
    Albtl.make_fun('alb_tl', Albtl_head, DEM.ret_str, Vars)

    AlbCorttl = DEM.DeEquation(AlbCorttl_body)
    AlbCorttl.make_fun('albcort_tl', 
        AlbCorttl_head, DEM.ret_str, Vars)

    # Torso Middle
    Albtm = DEM.DeEquation(Albtm_body)
    Albtm.make_fun('alb_tm', 
        Albtm_head, DEM.ret_str, Vars)

    AlbCorttm = DEM.DeEquation(AlbCorttm_body)
    AlbCorttm.make_fun('albcort_tm', 
        AlbCorttm_head, DEM.ret_str, Vars)

    # Torso Upper
    Albtu = DEM.DeEquation(Albtu_body)
    Albtu.make_fun('alb_tu', 
        Albtu_head, DEM.ret_str, Vars)

    AlbCorttu = DEM.DeEquation(AlbCorttu_body)
    AlbCorttu.make_fun('albcort_tu', 
        AlbCorttu_head, DEM.ret_str, Vars)

    # Plasma
    Albp = DEM.DeEquation(Albp_body)
    Albp.make_fun('alb_p', 
        Albp_head, DEM.ret_str, Vars)

    AlbCortp = DEM.DeEquation(AlbCortp_body)
    AlbCortp.make_fun('albcort_p', 
        AlbCortp_head, DEM.ret_str, Vars)
    
    
    Albrates = [TL_LymphAlbRate, TL_LymphACRate, TM_LymphAlbRate, 
            TM_LymphACRate, TU_LymphAlbRate, TU_LymphACRate, 
            TL_CapAlbRate, TL_CapACRate, TM_CapAlbRate, TM_CapACRate, 
            TU_CapAlbRate, TU_CapACRate, Albpcop, COP]



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
        cope      = COPE.cope(0, x, v, copeCurve, ACOP)
        T_LAR     = TL_LAR + TM_LAR + TU_LAR
        T_CAR     = TL_CAR + TM_CAR + TU_CAR
        albp      = Albp.alb_p(0, x, v, T_LAR,T_CAR)
        T_LACR    = TL_LACR + TM_LACR + TU_LACR
        T_CACR    = TL_CACR + TM_CACR + TU_CACR
        albcortp  = AlbCortp.albcort_p(0, x, v, T_LACR, T_CACR)    


        # Addition constraints
        KDA       = x[8] / x[9] - 3e-3

        lhs = [albtl, albcorttl, albtm, albcorttm, albtu, albcorttu, 
               cope, albp, albcortp, KDA]
        
        return lhs

    # Parameters
    full_dict = {**DEM.Paramsu, **DEM.stateVarsu, **DEM.stateVarsMassu}
    # Create the solver function
    y = lambda x : sys(0, x, full_dict)
    # create initial Guess 
    guess = DEM.np.zeros(10)
    guess[0]  = full_dict['Alb_tl']
    guess[1]  = full_dict['Alb-Cort_tl']
    guess[2]  = full_dict['Alb_tm']
    guess[3]  = full_dict['Alb-Cort_tm']
    guess[4]  = full_dict['Alb_tu']
    guess[5]  = full_dict['Alb-Cort_tu']
    guess[6]  = full_dict['Alb_p']
    guess[7]  = full_dict['Alb-Cort_p']
    guess[8] = full_dict['konA']
    guess[9] = full_dict['koffA']
    print(y(guess))

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




