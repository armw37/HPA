import DEModel as DEM

"""
CRH ACTH and Cortisol
"""

##########################################################################
# CRF and ACTH functions                                                 #
##########################################################################

CRF_head = "def crf(self, t, x, v):\n\trhs = "
CRF_body = ("v['k1']  * "  
           # Positive hill feedback by cort_p
           "(1 + v['zeta'] * ((v['Cort_pc'] ** v['alpha'])"
           " / (v['Cort_pc'] ** v['alpha'] + v['c'] ** v['alpha']))" 
           # Negative hill feedback by cort_p
           " - (v['psi'] * ((v['Cort_pc'] ** v['gamma'])"
           "  / (v['Cort_pc'] ** v['gamma'] + v['c3'] ** v['gamma']))))"
           # Generic waste
           "- v['w1'] * v['CRH']"
           )



ACTH_head = "def acth(self, t, x, v):\n\trhs = "
ACTH_body = ("v['k2'] * "
            # Negative hill feedback by cort_p
            "(1 - v['rho'] * (v['Cort_pc'] ** v['alpha'])"
            " /(v['Cort_pc'] ** v['alpha'] + v['c'] ** v['alpha']))"
            # Saturable CRH Stimulation
            #"* ((v['CRH'] ** v['omega'])"
            #"/ (v['CRH'] ** v['omega'] + v['cc'] ** v['omega']))"
            "* v['CRHc'] * v['v_p']"
            # Generic waste
            "- v['w2'] * v['ACTH']"
            )



##########################################################################
# Plasma Cortisol functions                                              #
##########################################################################



# Numerical is good
Cortp_head = "def cort_p(self, t, x, v):\n\trhs = "
Cortp_body = (# stimulation
              "v['k3'] * v['ACTHc'] * v['v_p']"   
              # lymph
              "+ v['f_lymp_w'] * v['Cort_ic']"
              # capillary
              "- v['f_cap_w'] * v['Cort_pc']"  
              # diffusion to interstitium
              "- v['Dcort_pi'] * (v['Cort_pc'] - v['Cort_ic'])"
              # diffusion to the kidney
              "- v['Dcort_pk'] * (v['Cort_pc'] - v['Cort_kc'])"            
              # Urine
              "-  v['ucp'] * v['Cort_p']"
              # Albumin Binding
              #"+ (-v['konA'] *(1/v['MWAlb']) * v['Alb_pc'] * v['Cort_pc']"
              #" * v['tonano']"
              # albumin off binding
              #"+ v['koffA'] * v['Alb-Cort_pc'] * (v['MWCort']/v['MWAC'])"
              # CBG binding
              #"- v['konC'] * (1/v['MWCBG']) * v['CBG_pc'] * v['Cort_pc']"
              #" * v['mictonano']"
              # CBG off binding
              #"+ v['koffC'] * v['CBG-Cort_pc'] * (v['MWCort']/v['MWCC']))"
              #"* v['v_p']"
             )

##########################################################################
# Interstitial Cortisol functions                                        #
##########################################################################

# numerical done but derivative is very low 
Corti_head = "def cort_i(self, t, x, v):\n\trhs = "
Corti_body = (# lymph
              "- v['f_lymp_w'] * v['Cort_ic'] "
              # capillary
              "+ v['f_cap_w'] * v['Cort_pc'] "
              # diffusion
              "+ v['Dcort_pi'] * (v['Cort_pc'] - v['Cort_ic']) "
              # diffusion
              "- v['Dcort_ik'] * (v['Cort_ic'] - v['Cort_kc'])"
              # diffusion
              "- v['Dcort_io'] * (v['Cort_ic'] - v['Cort_oc']) "
              
              ##################### Albumin #####################

              # ON TL
              #"+ (- v['konA'] * (1/v['MWAlb']) * v['Alb_tlc']"  
              #"* v['Cort_ic'] * v['tonano']" 
              # OFF TL
              #"+ v['koffA'] * v['Alb-Cort_tlc'] * (v['MWCort']/v['MWAC']))"
              #" * v['v_tl']"
              
              # ON TM
              #" + (- v['konA'] * (1/v['MWAlb']) * v['Alb_tmc']"  
              #"* v['Cort_ic'] * v['tonano']" 
              # OFF TM
              #"+ v['koffA'] * v['Alb-Cort_tmc'] * (v['MWCort']/v['MWAC']))"
              #" * v['v_tm']"

              # ON TU
              #"+ (- v['konA'] * (1/v['MWAlb']) * v['Alb_tuc']"  
              #"* v['Cort_ic']  * v['tonano']" 
              # OFF TU
              #"+ v['koffA'] * v['Alb-Cort_tuc'] * (v['MWCort']/v['MWAC']))"
              #" * v['v_tu']"


              ##################### CBG #####################

              # ON TL
              #"+ (- v['konC'] * (1/v['MWCBG']) * v['CBG_tlc']"  
              #"* v['Cort_ic'] * v['mictonano']" 
              # OFF TL
              #"+ v['koffC'] * v['CBG-Cort_tlc'] * (v['MWCort']/v['MWCC']))"
              #" * v['v_tl']"

              # ON TM
              #" + (- v['konC'] * (1/v['MWCBG']) * v['CBG_tmc']"  
              #"* v['Cort_ic'] * v['mictonano'] " 
              # ON TM
              #"+ v['koffC'] * v['CBG-Cort_tmc'] * (v['MWCort']/v['MWCC']))"
              #" * v['v_tm']"

              # ON TU
              #"+ (- v['konC'] * (1/v['MWCBG']) * v['CBG_tuc']"  
              #"* v['Cort_ic'] * v['mictonano'] " 
              # OFF TU
              #"+ v['koffC'] * v['CBG-Cort_tuc'] * (v['MWCort']/v['MWCC']))"
              #" * v['v_tu']"
             )


##########################################################################
# Cortisol in kidney and other cells                                     #
##########################################################################


MR_head = "def mr_b(self, t, x, v):\n\trhs = "
MR_body = ("v['Cort_kc'] * (1/v['MWCort']) "
           "/ (v['Kd_mr'] + v['Cort_kc'] * (1/v['MWCort']) "
           "+ v['Ald'] * (1/v['MWAld']))"
          )


MRA_head = "def mra_b(self, t, x, v):\n\trhs = "
MRA_body = ("v['Ald'] * (1/v['MWAld']) "
           "/ (v['Kd_mr'] + v['Cort_kc'] * (1/v['MWCort']) "
           "+ v['Ald'] * (1/v['MWAld']))"
          )


# Numerical good
Cortk_head = "def cort_k(self, t, x, v):\n\trhs = "
Cortk_body = (# Diffusion from interstitium
              "v['Dcort_ik'] * (v['Cort_ic'] - v['Cort_kc'])"
              # Diffusion from plasma
              " + v['Dcort_pk'] * (v['Cort_pc'] - v['Cort_kc'])"
              # 11B-HSD2
              "- ( v['Vm2'] * v['Cort_kc']" 
              "/(v['Cort_kc'] * (1/v['MWCort']) + v['km2'])) * v['v_k']"
             )



# Numerical Good
Corto_head = "def cort_o(self, t, x, v):\n\trhs = "
Corto_body = ( # Diffusion
             "v['Dcort_io'] * (v['Cort_ic'] - v['Cort_oc']) " 
             # 11B-HSD1
             "+ v['v_o'] * (v['Vm1'] * v['MWCort'] "
             "* (v['Cortisone_oc'] * (1/v['MWsone']))"
             " / (v['Cortisone_oc'] * (1/v['MWsone']) + v['km1'])) "
             # General Metabolism
             #"- v['mcort'] * v['Cort_o']"
             ) 


ToConc_head = "def cortToConc(self, x, v):\n\t"
ToConc_body = (# interstitial to concentration
               "v['Cort_ic'] = v['Cort_i'] / v['v_i'];"
               # plasma to concentration
               "v['Cort_pc'] = v['Cort_p'] / v['v_p'];"
               # other to concentration
               "v['Cort_oc'] = v['Cort_o'] / v['v_o'];"
               # kidney to concentration
               "v['Cort_kc'] = v['Cort_k'] / v['v_k'];"
               # CRH to concentration
               "v['CRHc'] = v['CRH'] / v['v_p'];"
               # ACTH to concentration
               "v['ACTHc'] = v['ACTH'] / v['v_p'];"
               # return v
               "rhs = v"
               )




def Unit_test(ToConc, CRH, ACTH, Cort_i, Cort_p, Cort_o, Cort_k, MR_B):
    
    # Creat the vector of Cortisone state vars
    testvec = [val for val in DEM.stateVarsMass.values()][:6]
    # Calculate concentration and add to params
    ps = {**DEM.Params, **DEM.stateVars, **DEM.stateVarsMass}
    ps = ToConc.cortToConc(testvec, ps)
    ps['mrb'] = MR_B.mr_b(0, testvec, ps).to('').magnitude
    # CRH
    print("CRH:", CRH.crf(0, testvec, 
        ps))
    # ACTH
    print("ACTH:", ACTH.acth(0, testvec, 
        ps).to(DEM.u.pg / DEM.u.min))
    # Cort_p
    print("Cort_p:", Cort_p.cort_p(0, testvec, 
        ps).to(DEM.u.ng / DEM.u.min))
    # Cort_i
    print("Cort_i:", Cort_i.cort_i(0, testvec, ps).to(DEM.u.ng / DEM.u.min))
    # Cort_k
    print("Cort_k:", Cort_k.cort_k(0, testvec, ps))
    # Cort_o
    print("Cort_o:", Cort_o.cort_o(0, testvec, ps))
    
    return 0


def numerical_test(ToConc, CRH, ACTH, Cort_i, Cort_p, Cort_o, Cort_k, MR_B):
    # Create the vector of Cortisone state vars
    testvec = [val for val in DEM.stateVarsMassu.values()][:6]
    # Calculate concentration and add to params
    ps = {**DEM.Paramsu, **DEM.stateVarsu, **DEM.stateVarsMassu}
    ps = ToConc.cortToConc(testvec, ps)
    ps['mrb'] = MR_B.mr_b(0, testvec, ps)
    # CRH
    print("CRH:", CRH.crf(0, testvec, 
        ps))#.to(DEM.u.ng / DEM.u.min))
    # ACTH
    print("ACTH:", ACTH.acth(0, testvec, ps))
    # Cort_p
    print("Cort_p:", Cort_p.cort_p(0, testvec, ps))
    # Cort_i
    print("Cort_i:", Cort_i.cort_i(0, testvec, ps))
    # Cort_k
    print("Cort_k:", Cort_k.cort_k(0, testvec, ps))
    # Cort_o
    print("Cort_o:", Cort_o.cort_o(0, testvec, ps))
    
    return 0




##########################################################################
# main                                                                   #
##########################################################################

def main():

    
    # Test the solution
    Vars = {
            "v\['CRH'\]": 'x[0]', "v\['ACTH'\]": 'x[1]',
            "v\['Cort_p'\]": 'x[2]', "v\['Cort_i'\]": 'x[3]',
            "v\['Cort_k'\]": 'x[4]', "v\['Cort_o'\]": 'x[5]',
           }

    CRH = DEM.DeEquation(CRF_body)
    CRH.make_fun('crf', CRF_head, DEM.ret_str, Vars)

    ACTH = DEM.DeEquation(ACTH_body)
    ACTH.make_fun('acth', ACTH_head, DEM.ret_str, Vars)

    Cort_p = DEM.DeEquation(Cortp_body)
    Cort_p.make_fun('cort_p', Cortp_head, DEM.ret_str, Vars)

    Cort_i = DEM.DeEquation(Corti_body)
    Cort_i.make_fun('cort_i', Corti_head, DEM.ret_str,
        Vars)

    MR_B = DEM.DeEquation(MR_body)
    MR_B.make_fun('mr_b', MR_head, DEM.ret_str,
        Vars)

    Cort_k = DEM.DeEquation(Cortk_body)
    Cort_k.make_fun('cort_k', Cortk_head, DEM.ret_str, Vars)

    Cort_o = DEM.DeEquation(Corto_body)
    Cort_o.make_fun('cort_o', Corto_head, DEM.ret_str, Vars) 
    
    ToConc = DEM.DeEquation(ToConc_body)
    ToConc.make_fun('cortToConc', ToConc_head,
        DEM.ret_str, Vars)

    
    # DE System
    def sys(t, x, v):
        # Calculate concentration    
        v = ToConc.cortToConc(x, v)
        v['mrb'] = MR_B.mr_b(t, x, v)
        ddt = [
               CRH.crf(t, x, v),
               ACTH.acth(t, x, v),
               Cort_p.cort_p(t, x, v),
               Cort_i.cort_i(t, x, v),
               Cort_k.cort_k(t, x, v),
               Cort_o.cort_o(t, x, v)
              ]
        return ddt

    # ignore protein binding for now

    # Solution 
    vs = {**DEM.Paramsu, **DEM.stateVarsu, **DEM.stateVarsMassu}   

    # Play with parms
    vs['k3'] = 0.0005
    #vs['Dcort_pi'] = 10.8
    #vs['Dcort_pk'] = 18
    #vs['Dcort_ik'] = 10
    #vs['Dcort_io'] = 50
    vs['Vm2']       = 10
    vs['Vm1']       = 0.1
    vs['Cortisone_oc'] = 0.1        

    i0 = [DEM.stateVarsMassu['CRH'],                
          DEM.stateVarsMassu['ACTH'],
          DEM.stateVarsMassu['Cort_p'],
          DEM.stateVarsMassu['Cort_i'],
          DEM.stateVarsMassu['Cort_k'],
          DEM.stateVarsMassu['Cort_o']]


    sol = DEM.solve_ivp(sys, [0, 6000], i0, args=(vs,),
             method='LSODA')
    z = sol.y.T
    
    # To Conc
    z[:,0] = z[:,0] / vs['v_p']
    z[:,1] = z[:,1] / vs['v_p']
    z[:,2] = z[:,2] / vs['v_p']
    z[:,3] = z[:,3] / vs['v_i']
    z[:,4] = z[:,4] / vs['v_k']
    z[:,5] = z[:,5] / vs['v_o']


    fig, axes = DEM.plt.subplots(2,3, sharex=True) 
    axes = axes.ravel()
    axes[0].plot(sol.t, z[:,0]) 
    axes[0].set_title('CRH'); #axes[0].set_ylim([0, 50])
    axes[1].plot(sol.t, z[:,1]) 
    axes[1].set_title('ACTH'); #axes[1].set_ylim([0, 50])
    axes[2].plot(sol.t, z[:,2]) 
    axes[2].set_title('Plasma Cortisol'); #axes[2].set_ylim([0, 10])
    axes[3].plot(sol.t, z[:,3]) 
    axes[3].set_title('Interstitial Cortisol'); #axes[3].set_ylim([0, 50])
    axes[4].plot(sol.t, z[:,4]) 
    axes[4].set_title('Kidney Cortisol'); #axes[0].set_ylim([0, 50])
    axes[5].plot(sol.t, z[:,5]) 
    axes[5].set_title('Other Cortisol'); #axes[1].set_ylim([0, 50])

    DEM.plt.tight_layout() 
    DEM.plt.show(block=True) 
    """

    # Running Parameter Estimates, going 1 at a time
    eVars = {
            "v\['k1'\]"     : 'x[0]',
            "v\['k2'\]"     : 'x[1]',
            "v\['k3'\]"     : 'x[2]',
            "v\['Cort_o'\]" : 'x[3]',
            "v\['Vm2'\]"    : 'x[4]',
            }
    
    CRH = DEM.DeEquation(CRF_body)
    CRH.make_fun('crf', CRF_head, DEM.ret_str, eVars)

    ACTH = DEM.DeEquation(ACTH_body)
    ACTH.make_fun('acth', ACTH_head, DEM.ret_str, eVars)

    Cort_p = DEM.DeEquation(Cortp_body)
    Cort_p.make_fun('cort_p', Cortp_head, DEM.ret_str, eVars)

    Cort_i = DEM.DeEquation(Corti_body)
    Cort_i.make_fun('cort_i', Corti_head, DEM.ret_str,
       eVars)

    MR_B = DEM.DeEquation(MR_body)
    MR_B.make_fun('mr_b', MR_head, DEM.ret_str,
        eVars)

    Cort_k = DEM.DeEquation(Cortk_body)
    Cort_k.make_fun('cort_k', Cortk_head, DEM.ret_str, eVars)

    Cort_o = DEM.DeEquation(Corto_body)
    Cort_o.make_fun('cort_o', Corto_head, DEM.ret_str, eVars)

    ToConc = DEM.DeEquation(ToConc_body)
    ToConc.make_fun('cortToConc', ToConc_head,
        DEM.ret_str, eVars)
       
    # DE System
    def sys(t, x, v):
        # Calculate concentration    
        v = ToConc.cortToConc(x, v)
        v['mrb'] = MR_B.mr_b(t, x, v)
        ddt = [
               CRH.crf(t, x, v),
               ACTH.acth(t, x, v),
               Cort_p.cort_p(t, x, v),
               Cort_i.cort_i(t, x, v),
               Cort_k.cort_k(t, x, v),
              ]
        return ddt

    # Solve For the unknown parameter
    full_dict = {**DEM.Paramsu, **DEM.stateVarsu, **DEM.stateVarsMassu}
    # Change the cort_kc to equal the ald concentration
    full_dict['Cort_kc'] = full_dict['Ald']
    full_dict['Cort_k']  = full_dict['Cort_kc'] * full_dict['v_k']
    y = lambda x : sys(0, x, full_dict) 
    guess = [2000, 0.1, 0.001, 80000, 50]
    steady = DEM.fsolve(y, guess) 
    print('The Steady State solution is:')
    print(steady)
    

    # Update params_dict
    full_dict['k1']     = steady[0]
    full_dict['k2']     = steady[1]  
    full_dict['k3']     = steady[2]
    full_dict['Cort_o'] = steady[3]
    full_dict['Vm2']    = steady[4]

    # Test to make sure it is right
    Vars = {
            "v\['CRH'\]"    : 'x[0]',
            "v\['ACTH'\]"   : 'x[1]',
            "v\['Cort_p'\]" : 'x[2]',
            "v\['Cort_i'\]" : 'x[3]',
            "v\['Cort_k'\]" : 'x[4]'
           }

    CRH = DEM.DeEquation(CRF_body)
    CRH.make_fun('crf', CRF_head, DEM.ret_str, Vars)

    ACTH = DEM.DeEquation(ACTH_body)
    ACTH.make_fun('acth', ACTH_head, DEM.ret_str, Vars)

    Cort_p = DEM.DeEquation(Cortp_body)
    Cort_p.make_fun('cort_p', Cortp_head, DEM.ret_str, Vars)

    Cort_i = DEM.DeEquation(Corti_body)
    Cort_i.make_fun('cort_i', Corti_head, DEM.ret_str,
        Vars)

    MR_B = DEM.DeEquation(MR_body)
    MR_B.make_fun('mr_b', MR_head, DEM.ret_str,
        Vars)

    Cort_k = DEM.DeEquation(Cortk_body)
    Cort_k.make_fun('cort_k', Cortk_head, DEM.ret_str, Vars)


    ToConc = DEM.DeEquation(ToConc_body)
    ToConc.make_fun('cortToConc', ToConc_head,
        DEM.ret_str, Vars)

    # DE System
    def sys(t, x, v):
        # Calculate concentration    
        v = ToConc.cortToConc(x, v)
        v['mrb'] = MR_B.mr_b(t, x, v)
        ddt = [
               CRH.crf(t, x, v),
               ACTH.acth(t, x, v),
               Cort_p.cort_p(t, x, v),
               Cort_i.cort_i(t, x, v),
               Cort_k.cort_k(t, x, v)
              ]
        return ddt




    states = [val for val in DEM.stateVarsMassu.values()][:6]
    states[5] = steady[3]
    states[4] = full_dict['Cort_k']
    print('Test the DE')
    print(sys(0, states, full_dict))
    """


    return 0


if __name__ == '__main__':
    main()


