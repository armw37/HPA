import DEModel as DEM

"""
Cortisone: Units work and the unit and numerical calculation match
"""


##########################################################################
# Cortisone  in kidney and other cells                                   #
##########################################################################


Cortisonek_head = "def cortisone_k(self, t, x, v):\n\trhs = "
Cortisonek_body = ( # Diffusion from interstitium
                   "v['Dcort_ik'] *(v['Cortisone_ic']"
                   " - v['Cortisone_kc']) "
                   # Diffusion to plasma
                   "- v['Dcort_pk'] *(v['Cortisone_kc']" 
                   " - v['Cortisone_pc'])"
                   # 11B-HSD2
                   "+ v['Vm2'] * v['Cort_kc']" 
                   "/(v['Cort_kc'] " 
                   "* (1/v['MWCort']) + v['km2']) * v['v_k']"
                   " * (v['MWsone']/v['MWCort'])"
                  )
    



# Numerical Good
Cortisonei_head = "def cortisone_i(self, t, x, v):\n\trhs = "
Cortisonei_body = ( # lymph
                   "- v['f_lymp_w'] * v['Cortisone_ic'] "
                   # capillary
                   "+ v['f_cap_w'] * v['Cortisone_pc']"
                   # Diffusion to plasma
                   "+ v['Dcort_pi'] * (v['Cortisone_pc'] - "
                   "v['Cortisone_ic'])"
                   # Diffusion to other
                   "- v['Dcort_io'] * (v['Cortisone_ic'] - "
                   "v['Cortisone_oc'])"
                   # Diffusion to kidney
                   "- v['Dcort_ik'] * (v['Cortisone_ic'] - "
                   "v['Cortisone_kc'])"
                   )


# Numerical Good
Cortisonep_head = "def cortisone_p(self, t, x, v):\n\trhs = "
Cortisonep_body = (
                   # lymph 
                   " + v['f_lymp_w'] * v['Cortisone_ic'] "
                   # cap
                   "- v['f_cap_w'] * v['Cortisone_pc'] "
                   # Diffusion to interstitium
                   "- v['Dcort_pi'] * (v['Cortisone_pc'] - "
                   "v['Cortisone_ic'])"
                    # Diffusion from kidney
                   "+ v['Dcort_pk'] *(v['Cortisone_kc']" 
                   " - v['Cortisone_pc'])"
                   # Urine
                   "- v['usp'] * v['Cortisone_p']"
                  )


Cortisoneo_head = "def cortisone_o(self, t, x, v):\n\trhs = "
Cortisoneo_body = (# Diffusion
                   "v['Dcort_io'] * (v['Cortisone_ic'] - "
                   "v['Cortisone_oc'])" 
                   # 11B-HSD1
                   "- v['v_o'] * (v['Vm1'] * (v['Cortisone_oc']"
                   " / (v['Cortisone_oc'] * (1/v['MWsone']) + v['km1'])))"
                   # General Metabolism
                   #"- v['msone'] * v['Cortisone_o']"     
                   )


ToConc_head = "def soneToConc(self, x, v):\n\t"
ToConc_body = (# interstitial to concentration
               "v['Cortisone_ic'] = v['Cortisone_i'] / v['v_i'];"
               # plasma to concentration
               "v['Cortisone_pc'] = v['Cortisone_p'] / v['v_p'];"
               # other to concentration
               "v['Cortisone_oc'] = v['Cortisone_o'] / v['v_o'];"
               # other to concentration
               "v['Cortisone_kc'] = v['Cortisone_k'] / v['v_k'];"              
               # return v
               "rhs = v"
               )



def Unit_test(ToConc, Cortisone_k, Cortisone_i, Cortisone_p, Cortisone_o):
    
    # Creat the vector of Cortisone state vars
    testvec = [val for val in DEM.stateVarsMass.values()]
    testvec = testvec[6:10]
    # Calculate concentration and add to params
    ps = {**DEM.Params, **DEM.stateVars, **DEM.stateVarsMass}
    ps = ToConc.soneToConc(testvec, ps)
    # Set a constant MR binding for testing 
    ps['mrb'] = 0.235169
    # Cortisone_k
    print("Cortisone_k:", Cortisone_k.cortisone_k(0, testvec, 
        ps).to(DEM.u.ng / DEM.u.min))
    # Cortisone_i
    print("Cortisone_i:", Cortisone_i.cortisone_i(0, testvec, ps))
    # Cortisone_p
    print("Cortisone_p:", Cortisone_p.cortisone_p(0, testvec, 
        ps).to(DEM.u.ng / DEM.u.min))
    # Cortisone_o
    print("Cortisone_o:", Cortisone_o.cortisone_o(0, testvec, ps))

    return 0


def numerical_test(ToConc, Cortisone_k ,Cortisone_i, Cortisone_p, Cortisone_o):
    # Create the vector of Cortisone state vars
    testvec = [val for val in DEM.stateVarsMassu.values()]
    testvec = testvec[6:10]
    # Calculate concentration and add to params
    ps = {**DEM.Paramsu, **DEM.stateVarsu, **DEM.stateVarsMassu}
    ps = ToConc.soneToConc(testvec, ps)
    # Set a constant MR binding for testing 
    ps['mrb'] = 0.235169
    # Cortisone_k
    print("Cortisone_k:", Cortisone_k.cortisone_k(0, testvec, ps))
    # Cortisone_i
    print("Cortisone_i:", Cortisone_i.cortisone_i(0, testvec, ps))
    # Cortisone_p
    print("Cortisone_p:", Cortisone_p.cortisone_p(0, testvec, ps))
    # Cortisone_o
    print("Cortisone_o:", Cortisone_o.cortisone_o(0, testvec, ps))
    
    return 0



##########################################################################
# main                                                                   #
##########################################################################

def main():

    
    # Test the solution
    Vars = {
            "v\['Cortisone_k'\]": 'x[0]'," v\['Cortisone_i'\]": 'x[1]',
            "v\['Cortisone_p'\]": 'x[2]', "v\['Cortisone_o'\]": 'x[3]'
           }

    Cortisone_k = DEM.DeEquation(Cortisonek_body)
    Cortisone_k.make_fun('cortisone_k', Cortisonek_head,
                     DEM.ret_str, Vars)

    Cortisone_i = DEM.DeEquation(Cortisonei_body)
    Cortisone_i.make_fun('cortisone_i', Cortisonei_head,
                     DEM.ret_str, Vars)

    Cortisone_p = DEM.DeEquation(Cortisonep_body)
    Cortisone_p.make_fun('cortisone_p', Cortisonep_head,
        DEM.ret_str, Vars)


    Cortisone_o = DEM.DeEquation(Cortisoneo_body)
    Cortisone_o.make_fun('cortisone_o', Cortisoneo_head,
        DEM.ret_str, Vars)

    ToConc = DEM.DeEquation(ToConc_body)
    ToConc.make_fun('soneToConc', ToConc_head,
        DEM.ret_str, Vars)


    # DE System
    def sys(t, x, v):
        # Calculate concentration    
        v = ToConc.soneToConc(x, v)
        ddt = [
               Cortisone_k.cortisone_k(t, x, v),              
               Cortisone_i.cortisone_i(t, x, v),
               Cortisone_p.cortisone_p(t, x, v),    
               Cortisone_o.cortisone_o(t, x, v)
              ]
        return ddt

    
    # Solution 
    vs = {**DEM.Paramsu, **DEM.stateVarsu, **DEM.stateVarsMassu}   

    # Play with parms
    #vs['Dcort_pi'] = 10.8
    #vs['Dcort_pk'] = 18
    #vs['Dcort_ik'] = 10
    #vs['Dcort_io'] = 50
    vs['Vm2']       = 20
    vs['Vm1']       = 0.1
    vs['Cort_kc']   = 0.1        

    i0 = [DEM.stateVarsMassu['Cortisone_k'],                
          DEM.stateVarsMassu['Cortisone_i'],
          DEM.stateVarsMassu['Cortisone_p'],
          DEM.stateVarsMassu['Cortisone_o'],
          ]


    sol = DEM.solve_ivp(sys, [0, 6000], i0, args=(vs,),
             method='LSODA')
    z = sol.y.T
    
    # To Conc
    z[:,0] = z[:,0] / vs['v_k']
    z[:,1] = z[:,1] / vs['v_i']
    z[:,2] = z[:,2] / vs['v_p']
    z[:,3] = z[:,3] / vs['v_o']

        
    fig, axes = DEM.plt.subplots(2,2, sharex=True) 
    axes = axes.ravel()
    axes[0].plot(sol.t, z[:,0]) 
    axes[0].set_title('Kidney Cortisone')#; axes[0].set_ylim([0, 1])
    axes[1].plot(sol.t, z[:,1])
    axes[1].set_title('Interstitial Cortisone')#; axes[1].set_ylim([0, 1])
    axes[2].plot(sol.t, z[:,2]) 
    axes[2].set_title('Plasma Cortisone')#; axes[2].set_ylim([0, 1])
    axes[3].plot(sol.t, z[:,3]) 
    axes[3].set_title('Other Cortisone')#; axes[2].set_ylim([0, 1])

    DEM.plt.tight_layout() 
    DEM.plt.show(block=True) 
    
    """
    # Running Parameter Estimates, going 1 at a time
    eVars = {
            "v\['Cortisone_k'\]"   : 'x[0]',
            "v\['Cortisone_o'\]"   : 'x[1]',
            #"v\['Vm1'\]"           : 'x[2]',
            #"v\['msone'\]"         : 'x[3]',
            }
    
    Cortisone_k = DEM.DeEquation(Cortisonek_body)
    Cortisone_k.make_fun('cortisone_k', Cortisonek_head,
                     DEM.ret_str, eVars)

    Cortisone_i = DEM.DeEquation(Cortisonei_body)
    Cortisone_i.make_fun('cortisone_i', Cortisonei_head,
                     DEM.ret_str, eVars)

    Cortisone_p = DEM.DeEquation(Cortisonep_body)
    Cortisone_p.make_fun('cortisone_p', Cortisonep_head,
        DEM.ret_str, eVars)


    Cortisone_o = DEM.DeEquation(Cortisoneo_body)
    Cortisone_o.make_fun('cortisone_o', Cortisoneo_head,
        DEM.ret_str, eVars)

    ToConc = DEM.DeEquation(ToConc_body)
    ToConc.make_fun('soneToConc', ToConc_head,
        DEM.ret_str, eVars)
 
    
    # DE System
    def sys(t, x, v):
        # Calculate concentration    
        v = ToConc.soneToConc(x, v)
        ddt = [
               Cortisone_k.cortisone_k(t, x, v),              
               Cortisone_i.cortisone_i(t, x, v),
               #Cortisone_p.cortisone_p(t, x, v),    
               #Cortisone_o.cortisone_o(t, x, v)
              ]
        return ddt


    # Solution 
    full_dict = {**DEM.Paramsu, **DEM.stateVarsu, **DEM.stateVarsMassu}
    # Assume a constant mr_binding
    full_dict['mrb'] = 0.235169
    full_dict['Cort_kc'] = full_dict['Ald']

    y = lambda x : sys(0, x, full_dict) 
    guess = [3,30000]
    steady = DEM.fsolve(y, guess) 
    print('The Steady State solution is:')
    print(steady)


    # Update params_dict
    #full_dict['k1']     = steady[0]
    #full_dict['k2']     = steady[1]  
    #full_dict['k3']     = steady[2]
    #full_dict['Cort_o'] = steady[3]

    # Test to make sure it is right
    Vars = {
            "v\['Cortisone_k'\]"    : 'x[0]',
            "v\['Cortisone_i'\]"    : 'x[1]',
            #"v\['Cort_p'\]" : 'x[2]',
            #"v\['Cort_i'\]" : 'x[3]',
           }

    Cortisone_k = DEM.DeEquation(Cortisonek_body)
    Cortisone_k.make_fun('cortisone_k', Cortisonek_head,
                     DEM.ret_str, Vars)

    Cortisone_i = DEM.DeEquation(Cortisonei_body)
    Cortisone_i.make_fun('cortisone_i', Cortisonei_head,
                     DEM.ret_str, Vars)

    ToConc = DEM.DeEquation(ToConc_body)
    ToConc.make_fun('soneToConc', ToConc_head,
        DEM.ret_str, Vars)


    states = [val for val in DEM.stateVarsMassu.values()][6:10]
    states[0] = steady[0];
    full_dict['Cortisone_o'] = steady[1]
    print('Test the DE')
    print(sys(0, states, full_dict))


   
    """ 
    return 0

if __name__ == '__main__':
    main()


