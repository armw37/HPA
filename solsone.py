import DEModel as DEM
import Cortisol as sol
import Cortisone as sone

##########################################################################
# Solve for unknowns of Cortisol and Cortisone                           #
##########################################################################

"""
Solution Params:
# CRH
'CRH'       
'zeta'      
'alpha'     
'c'        
'vi'       
'gamma'     
'c3'        
'w1'        
# ACTH
'ACTH'      
'rho'       
'w2'        
# Cortisol
'Cort_p'    
'f_lymp_w'  
'v_p'       
'v_i'       
'f_cap_w'   
'Dcort_pi'  
'Dcort_pk'  
'ucp'       
'konA'      
'koffA'     
'konC'      
'koffC'     
# Cort_i
'Cort_i'   
'Dcort_ik'  
'Dcort_io'  
# Cort_k
'v_k'        
'v_o'         
'km2'       
'Kd_mr'     
'Ald'       
# cort_o
'Vm1'       
'km1'       
'usp'       
'Cortisone_i' 
'Cortisone_p' 

Solution Variables:
'k1'        
'k2'       
'k3'       
'Cort_k'   
'Vm2'      
'Cort_o'    
'Vm1'      
'mcort'    
'Cortisone_k'  
'Cortisone_o'  
'msone'

EQS     = 10
Unkowns = 11


Extra Constraints:


"""

def main(sol, sone):
    
    # Vars
    Vars = {"v\['k1'\]"         : 'x[0]', "v\['k2'\]"         : 'x[1]',
            "v\['k3'\]"         : 'x[2]', "v\['Cort_k'\]"     : 'x[3]',
            "v\['Vm2'\]"        : 'x[4]', "v\['Cort_o'\]"     : 'x[5]', 
            "v\['Vm1'\]"        : 'x[6]', "v\['Cortisone_k'\]": 'x[7]', 
            "v\['Cortisone_o'\]": 'x[8]', "v\['msone'\]"      : 'x[9]'
           }

    Vars =  {
             "v\['CRH'\]": 'x[0]', "v\['ACTH'\]": 'x[1]',
             "v\['Cort_p'\]": 'x[2]', "v\['Cort_i'\]": 'x[3]',
             "v\['Cort_k'\]": 'x[4]', "v\['Cort_o'\]": 'x[5]',
             "v\['Cortisone_k'\]": 'x[6]'," v\['Cortisone_i'\]": 'x[7]',
             "v\['Cortisone_p'\]": 'x[8]', "v\['Cortisone_o'\]": 'x[9]'
            }

    ######################################################################
    # Cortisol eqs                                                       #
    ######################################################################  

    DEM
    CRH = DEM.DeEquation(sol.CRF_body)
    CRH.make_fun('crf', sol.CRF_head, DEM.ret_str, Vars)

    ACTH = DEM.DeEquation(sol.ACTH_body)
    ACTH.make_fun('acth', sol.ACTH_head, DEM.ret_str, Vars)

    Cort_p = DEM.DeEquation(sol.Cortp_body)
    Cort_p.make_fun('cort_p', sol.Cortp_head, DEM.ret_str, Vars)

    Cort_i = DEM.DeEquation(sol.Corti_body)
    Cort_i.make_fun('cort_i', sol.Corti_head, DEM.ret_str,
        Vars)

    MR_B = DEM.DeEquation(sol.MR_body)
    MR_B.make_fun('mr_b', sol.MR_head, DEM.ret_str,
        Vars)

    Cort_k = DEM.DeEquation(sol.Cortk_body)
    Cort_k.make_fun('cort_k', sol.Cortk_head, DEM.ret_str, Vars)

    Cort_o = DEM.DeEquation(sol.Corto_body)
    Cort_o.make_fun('cort_o', sol.Corto_head, DEM.ret_str, Vars) 
    
    SolToConc = DEM.DeEquation(sol.ToConc_body)
    SolToConc.make_fun('cortToConc', sol.ToConc_head,
        DEM.ret_str, Vars)

    ######################################################################
    # Cortisone eqs                                                      #
    ###################################################################### 

    Cortisone_k = DEM.DeEquation(sone.Cortisonek_body)
    Cortisone_k.make_fun('cortisone_k', sone.Cortisonek_head,
                     DEM.ret_str, Vars)

    Cortisone_i = DEM.DeEquation(sone.Cortisonei_body)
    Cortisone_i.make_fun('cortisone_i', sone.Cortisonei_head,
                     DEM.ret_str, Vars)

    Cortisone_p = DEM.DeEquation(sone.Cortisonep_body)
    Cortisone_p.make_fun('cortisone_p', sone.Cortisonep_head,
        DEM.ret_str, Vars)


    Cortisone_o = DEM.DeEquation(sone.Cortisoneo_body)
    Cortisone_o.make_fun('cortisone_o', sone.Cortisoneo_head,
        DEM.ret_str, Vars)

    SoneToConc = DEM.DeEquation(sone.ToConc_body)
    SoneToConc.make_fun('soneToConc', sone.ToConc_head,
        DEM.ret_str, Vars)


    ######################################################################
    # Make sys                                                           #
    ######################################################################

    def sys(t, x, v):
    
        # Calculate concentration for Cortisol components
        v = SolToConc.cortToConc(x, v)
   
        #print("Params after 1: ")
        #print('\n'.join([key + ': ' + str(val) for key, val in v.items()]))

        # Cortisol Base eqautions
        crh   = CRH.crf(0, x, v)
        acth  = ACTH.acth(0, x, v)
        cortp = Cort_p.cort_p(0, x, v)
        corti = Cort_i.cort_i(0, x, v)
        cortk = Cort_k.cort_k(0, x, v)
        corto = Cort_o.cort_o(0, x, v)

        # Calculate concentration for CBG components
        v = SoneToConc.soneToConc(x, v)

        #print("Params after 2: ")
        #print('\n'.join([key + ': ' + str(val) for key, val in v.items()]))

        # Cortisone Base equations
        sonek = Cortisone_k.cortisone_k(0, x, v)
        sonei = Cortisone_i.cortisone_i(0, x, v)
        sonep = Cortisone_p.cortisone_p(0, x, v)
        soneo = Cortisone_o.cortisone_o(0, x, v)

        # additonal Constraint: %30 of production is from 11B-HSD2
        #prod = (
        #        ((x[4] * ((1 - v['mrb']) * x[3]/v['v_k']) 
        #        /((1 - v['mrb']) * (x[3]/v['v_k']) * (1/v['MWCort'])
        #        + v['km2'])) * v['v_k']) 
        #        - 0.3 * (x[2] * v['ACTH'])
        #       )
        #eq11Bhsd = (# 11B-HSD2
        #            (x[4] * ((1 - v['mrb']) * x[3]/v['v_k']) 
        #            /((1 - v['mrb']) * (x[3]/v['v_k']) * (1/v['MWCort'])
        #            + v['km2'])) * v['v_k']
        #            # 11B-HSD1
        #            - v['v_o'] * (x[6] * v['MWCort']
        #            * (x[9]/v['v_o'] * (1/v['MWsone']))
        #            / (x[9]/v['v_o'] * (1/v['MWsone']) + v['km1']))
        #           )

                
        lhs = [crh, acth, cortp, corti, cortk, corto, 
                sonek, sonei, sonep, soneo]
        
        return lhs

    # Parameters
    full_dict = {**DEM.Paramsu, **DEM.stateVarsu, **DEM.stateVarsMassu}
    
    # Update the values for the Albumin and CBG
    full_dict['Alb_tlc']     = 5.42588639e+01 / full_dict['v_tl']
    full_dict['Alb-Cort_tlc']= 1.38170786e+06 / full_dict['v_tl']
    full_dict['Alb_tmc']     = 5.44027971e+01 / full_dict['v_tm']
    full_dict['Alb-Cort_tmc']= 1.38537314e+06 / full_dict['v_tm']
    full_dict['Alb_tuc']     = 4.38962198e+01 / full_dict['v_tu']
    full_dict['Alb-Cort_tuc']= 1.11782201e+06 / full_dict['v_tu']
    full_dict['Alb_pc']      = 1.50676010e+02 / full_dict['v_p']
    full_dict['Alb-Cort_pc'] = 3.83698096e+06 / full_dict['v_p']
    full_dict['CBG_tlc']     = 3.00639349e+04 / full_dict['v_tl']
    full_dict['CBG-Cort_tlc']= 1.01650433e+06 / full_dict['v_tl']
    full_dict['CBG_tmc']     = 2.71639552e+04 / full_dict['v_tm']
    full_dict['CBG-Cort_tmc']= 9.18451897e+05 / full_dict['v_tm']
    full_dict['CBG_tuc']     = 2.53549308e+04 / full_dict['v_tu']
    full_dict['CBG-Cort_tuc']= 8.57286213e+05 / full_dict['v_tu']
    full_dict['CBG_pc']      = 1.13832238e+05 / full_dict['v_p']
    full_dict['CBG-Cort_pc'] = 1.18000000e+03 / full_dict['v_p']
 
    
    # Vary Parameters

    # Play with parms
    full_dict['k3'] = 0.01
    #full_dict['Dcort_pi'] = 10.8
    #full_dict['Dcort_pk'] = 18
    #full_dict['Dcort_ik'] = 10
    #full_dict['Dcort_io'] = 50
    full_dict['Vm2']       = 10
    full_dict['Vm1']       = 1
    #full_dict['Cortisone_oc'] = 0.1 
    
    
    
    i0 = [DEM.stateVarsMassu['CRH'],                
          DEM.stateVarsMassu['ACTH'],
          DEM.stateVarsMassu['Cort_p'],
          DEM.stateVarsMassu['Cort_i'],
          DEM.stateVarsMassu['Cort_k'],
          DEM.stateVarsMassu['Cort_o'],
          DEM.stateVarsMassu['Cortisone_k'],                
          DEM.stateVarsMassu['Cortisone_i'],
          DEM.stateVarsMassu['Cortisone_p'],
          DEM.stateVarsMassu['Cortisone_o'],
          ]   
    
    #print("Params before: ")
    #print('\n'.join([key + ': ' + str(val) for key, val in full_dict.items()]))
    
    #print("Vector: ", i0)
    #print(sys(0, i0, full_dict))

    # Testvector
    #test = [24987.31050334082, 68495.15256950058,
    #    9979.367714267932, 37102.62031531384,
    #    0.3364453742701157, 2734.609949664298,
    #    19.497703036478587, 12166.747643672683,
    #    3248.7923396605634, 25285.184265632135]
   
    #print("Params before: ")
    #print('\n'.join([key + ': ' + str(val) for key, val in full_dict.items()]))
    
    #print("Vector: ", test)
    #print(sys(0, i0, full_dict))
    #sys(0, test, full_dict)


    # test ODE sol
    sol = DEM.solve_ivp(sys, [0, 20], i0, args=(full_dict,),
             method='LSODA')
    z = sol.y.T
    print(z)
    
    # To Conc
    #z[:,0] = z[:,0] / full_dict['v_p']
    #z[:,1] = z[:,1] / full_dict['v_p']
    #z[:,2] = z[:,2] / full_dict['v_p']
    #z[:,3] = z[:,3] / full_dict['v_i']
    #z[:,4] = z[:,4] / full_dict['v_k']
    #z[:,5] = z[:,5] / full_dict['v_o']
    #z[:,6] = z[:,6] / full_dict['v_k']
    #z[:,7] = z[:,7] / full_dict['v_i']
    #z[:,8] = z[:,8] / full_dict['v_p']
    #z[:,9] = z[:,9] / full_dict['v_o']

    
    fig, axes = DEM.plt.subplots(2,5, sharex=True) 
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
    axes[6].plot(sol.t, z[:,6]) 
    axes[6].set_title('Kidney Cortisone')#; axes[0].set_ylim([0, 1])
    axes[7].plot(sol.t, z[:,7])
    axes[7].set_title('Interstitial Cortisone')#; axes[1].set_ylim([0, 1])
    axes[8].plot(sol.t, z[:,8]) 
    axes[8].set_title('Plasma Cortisone')#; axes[2].set_ylim([0, 1])
    axes[9].plot(sol.t, z[:,9]) 
    axes[9].set_title('Other Cortisone')#; axes[2].set_ylim([0, 1])


    DEM.plt.tight_layout() 
    DEM.plt.show(block=True)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    """    
    
    # Create the solver function
    y = lambda x : sys(0, x, full_dict)
    # create initial Guess 
    guess = DEM.np.zeros(10)
    guess[0]  = full_dict['k1']
    guess[1]  = full_dict['k2']
    guess[2]  = full_dict['k3']
    guess[3]  = full_dict['Cort_k']
    guess[4]  = full_dict['Vm2']
    guess[5]  = full_dict['Cort_o']
    guess[6]  = full_dict['Vm1']
    guess[7]  = full_dict['Cortisone_k']
    guess[8]  = full_dict['Cortisone_o']   
    guess[9] = full_dict['msone']



    print('Testing the return of the system')
    row_format ="{:<11}" * 10
    print(row_format.format(* ['crh', 'acth', 'cortp', 'cortk', 
        'sonek', 'corti', 'corto', 'sonep', 'sonei', 'soneo']))

    nrow_format ="{:<11.3f}" * 10
    print(nrow_format.format(* y(guess)))


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

    """

    return 0

if __name__ == '__main__':
    main(sol, sone)
