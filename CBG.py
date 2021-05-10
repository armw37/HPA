import DEModel as DEM

"""
Working well for no binding interpreation of the sigma c is weird

"""

##########################################################################
# plasma CBG and CBG bound Cortisol functions                            #
##########################################################################

# Model as mass

# Plasma CBG
CBGp_head = "def cbg_p(self, t, x, v):\n\trhs = "
CBGp_body = (# Stimulation
             "v['k5'] "
             # lymph
             "+ (v['f_lymp_w'] * (v['Klp'] * v['CBG_p']"
             " + (1 - v['Klp']) * v['CBG_i']))"  
             # cap
             "- (v['f_cap_w'] * v['sigC'] * (v['CBG_p'] - v['CBG_i']))"  
             # half life rate
             "- (v['half'] * v['CBG_p'])"  
            )



# Interstitial CBG gotta remember the on off reactions are concentration
CBGi_head = "def cbg_i(self, t, x, v):\n\trhs = "
CBGi_body = (# lymph
             "+ (v['f_lymp_w'] * (v['Klp'] * v['CBG_p']"
             " + (1 - v['Klp']) * v['CBG_i']))"  
             # cap
             "- (v['f_cap_w'] * v['sigC'] * (v['CBG_p'] - v['CBG_i']))"             
            )



CBGtoconc_head = "def cbgtoconc(self, x, v):\n\trhs = "
CBGtoconc_body = (
                 "[v['CBG_p'] / v['v_p'], "
                 " v['CBG_i'] / v['v_i']] "
                 )


##########################################################################
# main                                                                   #
##########################################################################

def main():

    # solve for the steady state guesses
    parmvars = {"v\['k5'\]": 'x[0]', "v\['sigC'\]": 'x[1]'}
    CBG_p = DEM.DeEquation(CBGp_body)
    CBG_p.make_fun('cbg_p', CBGp_head, DEM.ret_str, parmvars)

    CBG_i = DEM.DeEquation(CBGi_body)
    CBG_i.make_fun('cbg_i', CBGi_head, DEM.ret_str, parmvars)


    # DE System
    def CBG_sys(t, x, v):
        ddt = [CBG_p.cbg_p(t, x, v), CBG_i.cbg_i(t, x, v)]
        return ddt
   
    full_dict = {**DEM.Paramsu, **DEM.stateVarsu}
    y = lambda x : CBG_sys(0, x, full_dict) 
    guess = [3, 0.98]
    steady = DEM.fsolve(y, guess) 
    print(steady) 

    # update the parameters
    full_dict['k5'] = steady[0]; full_dict['sigC'] = steady[1]

    # Test the solution
    Vars = {"v\['CBG_p'\]": 'x[0]', "v\['CBG_i'\]": 'x[1]'}
    CBG_p.make_fun('cbg_p', CBGp_head, DEM.ret_str, Vars)

    CBG_i.make_fun('cbg_i', CBGi_head, DEM.ret_str, Vars)

    CBGtoconc = DEM.DeEquation(CBGtoconc_body)
    CBGtoconc.make_fun('cbgtoconc', CBGtoconc_head,
        DEM.ret_str, Vars)

    # DE System
    def CBG_sys(t, x, v):
        x = CBGtoconc.cbgtoconc(x, v)
        ddt = [CBG_p.cbg_p(t, x, v), CBG_i.cbg_i(t, x, v)]
        return ddt


    # Solution 
    i0 = [35 * DEM.Paramsu['v_p'], 35 * 0.6 * DEM.Paramsu['v_i']]
    sol = DEM.solve_ivp(CBG_sys, [0, 6000], i0, args=(full_dict,),
             method='LSODA')
    z = sol.y.T
   
    y = lambda x: CBGtoconc.cbgtoconc(x, full_dict)
    c = DEM.np.apply_along_axis(y, 1, z)
    
    fig, axes = DEM.plt.subplots(2,2, sharex=True) 
 
    axes[0,0].plot(sol.t, z[:,0]) 
    axes[0,0].set_title('Plasma Mass')
    axes[0,0].set_ylim([10000, 500000])
    axes[0,1].plot(sol.t, z[:,1]) 
    axes[0,1].set_title('Interstial Mass')
    axes[0,1].set_ylim([10000, 500000])   
    axes[1,0].plot(sol.t, c[:,0]) 
    axes[1,0].set_title('Plasma Concentration') 
    axes[1,0].set_ylim([0, 50])
    axes[1,1].plot(sol.t, c[:,1]) 
    axes[1,1].set_title('Interstial Concentration')
    axes[1,1].set_ylim([0, 50])

    DEM.plt.tight_layout() 
    DEM.plt.show(block=True) 
 
    return 0


if __name__ == '__main__':
    main() 





