# testing what I think is the inconsistent part of the model

import DEModel as DEM

vs = {**DEM.Paramsu, **DEM.stateVarsu, **DEM.stateVarsMassu}



def sys(x, v):
    # calculate MRB
    mrb = (
           (x[0]/v['v_p'] * (1/v['MWCort'])) 
           / (v['Kd_mr'] + (x[0]/v['v_p'] * (1/v['MWCort']))
           + v['Ald'] * (1/v['MWAld']))
          )

    # Calc Cortk
    Cortk = (# Diffusion from interstitium
              v['Dcort_ik'] * (v['Cort_ic'] - (1-mrb)*(x[0]/v['v_p']))
             # Diffusion from plasma
             + v['Dcort_pk'] * (v['Cort_pc'] - (1-mrb)*(x[0]/v['v_p']))
             # 11B-HSD2
             - ( x[1] * ((1 - mrb) * (x[0]/v['v_p']))
             /((1 - mrb) * (x[0]/v['v_p']) * (1/v['MWCort'])
             + v['km2'])) * v['v_k']
            )



    Cortisonek = ( # Diffusion from interstitium
                  v['Dcort_ik'] * (v['Cortisone_ic'] - (x[2]/v['v_p'])) 
                  # Diffusion to plasma
                  + v['Dcort_pk'] * (v['Cortisone_pc'] - (x[2]/v['v_p'])) 
                  # 11B-HSD2
                  + (x[1] * (v['MWsone']/v['MWCort']) 
                  * ((1 - mrb) * (x[0]/v['v_p'])) 
                  /((1 - mrb) * (x[0]/v['v_p'])
                  * (1/v['MWCort']) + v['km2'])) * v['v_k']
                  )
    

    Cortisonep = (
                  # lymph 
                  + v['f_lymp_w'] * v['Cortisone_ic']
                  # cap
                  - v['f_cap_w'] * v['Cortisone_pc']
                  # Diffusion to interstitium
                  - v['Dcort_pi'] * (v['Cortisone_pc'] 
                  - v['Cortisone_ic'])
                  # Diffusion from kidney
                  - v['Dcort_pk'] *(v['Cortisone_pc'] - (x[2]/v['v_p']))
                  # Urine
                  - v['usp'] * v['Cortisone_p']
                 )


    return [Cortk, Cortisonek, Cortisonep]



i0 = [val for key, val in full_dict.items() if key in stateVariables]
sol = DEM.solve_ivp(sone_sys, [0, 6000], i0, args=(full_dict,),
    method='LSODA')
z = sol.y.T





# Create the solver function
y = lambda x : sys(x, vs)
# create initial Guess 
guess = DEM.np.zeros(3)
guess[0]  = vs['Cort_k']
guess[1]  = vs['Vm2'] 
guess[2]  = vs['Cortisone_k']


print(y(

# Solve
steady, info, ier, mesg  = DEM.fsolve(y, guess, full_output=True)
if ier == 1:
    print('The Steady State solution is:')
    print(steady)
    print('Checking the Sys: ')
    print(sys(steady, vs))
else:
    print('Solution was not found')
    print('Error:', mesg)

