""" 
    Parameter fitting options
    # Cortisol
    'Dcort_pi' : 10.8      * u.ml / u.min,
    # Cort_i
    'Dcort_ik' : 15        * u.ml / u.min, # need to add cell surface
    'Dcort_io' : 15        * u.ml / u.min, # need to add cell surface 
    # Cort_k
    'Vm2'      : 48        * u.nmol / (u.l * u.min), # no idea here
    # cort_o
    'Vm1'      : 48        * u.nmol / (u.l * u.min), # match / play around a bit
    cortisone_o
    'mcort'    : 0.5       * 1 / u.min,              # guesses
    'msone'    : 0.005     * 1 / u.min,              # guesses
    State Options
    'Cortisone_k' : 0.5 * u.ng / u.ml,              
    """

    # Just going to treat acth and crh as constant
    # 8 eqs and which are my 8 unknowns
    # I can add a constraint on cortk to be =ald or maybe 1/10, 1/100
    # assume at eq the cort_o and cortisone_o are at equilibrium with 
    # interstitial space

    """
    # solve for the steady state guesses
    parmvars = {"v\['Dcort_pi'\]" : 'x[0]', "v\['Dcort_ik'\]" : 'x[1]',
            "v\['Dcort_io'\]" : 'x[2]', "v\['Vm2'\]" : 'x[3]',
            "v\['Vm1'\]" : 'x[4]', "v\['mcort'\]" : 'x[5]', 
            "v\['msone'\]" : 'x[6]', "v\['Cortisone_k'\]" : 'x[7]'}
    
    Cort_p = DEM.DeEquation(Cortp_body)
    Cort_p.make_fun('cort_p', Cortp_head, DEM.ret_str, parmvars)


    Cort_i = DEM.DeEquation(Corti_body)
    Cort_i.make_fun('cort_i', Corti_head, DEM.ret_str,
        parmvars)


    MR_B = DEM.DeEquation(MR_body)
    MR_B.make_fun('mr_b', MR_head, DEM.ret_str,
        parmvars)


    Cort_k = DEM.DeEquation(Cortk_body)
    Cort_k.make_fun('cort_k', Cortk_head, DEM.ret_str, parmvars)


    Cort_o = DEM.DeEquation(Corto_body)
    Cort_o.make_fun('cort_o', Corto_head, DEM.ret_str, parmvars)

    # DE System
    def Cort_sys(t, x, v):
        v['mrb'] = MR_B.mr_b(t, x, v)
        ddt = [
               Cort_p.cort_p(t, x, v), Cort_i.cort_i(t, x, v),
               Cort_k.cort_k(t, x, v), Cort_o.cort_o(t, x, v),
               Cortisone_k.cortisone_k(t, x, v),
               Cortisone_i.cortisone_i(t, x, v),
               Cortisone_p.cortisone_p(t, x, v),    
               Cortisone_o.cortisone_o(t, x, v)
              ]
        return ddt
   
    full_dict = {**DEM.Paramsu, **DEM.stateVarsu}
    # Modify cortk
    full_dict['Cort_k'] = full_dict['Ald']
    y = lambda x : Cort_sys(0, x, full_dict) 
    guess = [10, 8, 20, 50, 50, 0.005, 0.005, 1]
    steady = DEM.fsolve(y, guess) 
    print(steady) 

    # update the parameters
    full_dict['Dcort_pi']    = steady[0]
    full_dict['Dcort_ik']    = steady[1]
    full_dict['Dcort_io']    = steady[2]
    full_dict['Vm2']         = steady[3]
    full_dict['Vm1']         = steady[4]
    full_dict['mcort']       = steady[5]
    full_dict['msone']       = steady[6]
    full_dict['Cortisone_k'] = steady[7]

