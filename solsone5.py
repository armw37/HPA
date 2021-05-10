import DEModel as DEM

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
# cort_o
'Vm1'       
'km1'       
'usp'       
'Cortisone_i' 
'Cortisone_p' 

Solution Variables:
# CRH
    'k1'        
# ACTH
    'k2'       
# Cortp
    'k3'
# Corti
    DCortpi
# Cortk
    Vm2
#Corto
    Dcortio
# Cortisonek
    Dcortpk     
# Cortisonei
    Dcortik
# Cortisonep   
    kcp
# Cortisoneo
    'Vm1'      
     
  

EQS     = 10
Unkowns = 10


"""

def main():
    

    ######################################################################
    # Make sys                                                           #
    ######################################################################
    #x[0] : Dik
    #x[1] : Vm1
    #x[2] : Vm2
    #x[3] : Cortisoneo
    #x[4] : Dcort_io

    def sys(t, x, v): 
        CRF = (x[0]  * 
               # Positive hill feedback by cort_p
               (1 + v['zeta'] * (((v['Cortp']/v['v_p']) ** v['alpha'])
               / ((v['Cortp']/v['v_p']) ** v['alpha'] + v['c'] ** v['alpha'])) 
               # Negative hill feedback by cort_p
               - (v['psi'] * (((v['Cortp']/v['v_p']) ** v['gamma'])
               / ((v['Cortp']/v['v_p']) ** v['gamma'] + v['c3'] ** v['gamma']))))
               # Generic waste
               - v['w1'] * v['CRH']
              )



        ACTH = (x[1] * 
                # Negative hill feedback by cort_p
                (1 - v['rho'] * ((v['Cortp']/v['v_p']) ** v['alpha'])
                /((v['Cortp']/v['v_p']) ** v['alpha'] + v['c'] ** v['alpha']))
               * (v['CRH']/v['v_p']) * v['v_p']
               # Generic waste
               - v['w2'] * v['ACTH']
               )



        Cortp = (# stimulation
                 x[2] * (v['ACTH']/v['v_p']) * v['v_p']   
                 # diffusion to interstitium
                 - v['Dcort_pi'] * ((v['Cortp']/v['v_p']) - (v['Corti']/v['v_i']))
                 # diffusion to the kidney
                 - v['Dcort_pk'] * ((v['Cortp']/v['v_p']) - (v['Cortk']/v['v_k']))            
                 # Source from Vm1
                 +  (x[3] * v['MWCort'] * v['v_p'] 
                 * (v['Cortisonep']/(v['v_p'] * v['MWsone']))
                 / ((v['Cortisonep']/(v['v_p'] * v['MWsone'])) + v['km1']))
                 # Urine
                 - v['ucp'] * v['Cortp']
                )

        Cortisonep = (# Constant source 
                      x[4]
                      # Diffusion from kidney
                      - v['Dcort_pk'] *((v['Cortisonep']/v['v_p']) 
                      - (v['Cortisonek']/v['v_k']))
                      # Vm1 sink
                      +  (x[3] * v['MWCort'] * v['v_p'] 
                      * (v['Cortisonep']/(v['v_p'] * v['MWsone']))
                      / ((v['Cortisonep']/(v['v_p'] * v['MWsone'])) + v['km1']))
                      # Urine
                      - v['usp'] * v['Cortisonep']
                     )       
        
        
        # changed v['Cortisonep'] to high cort o   
        Corti = (
                 # diffusion
                 + v['Dcort_pi'] * ((v['Cortp']/v['v_p']) - (v['Corti']/v['v_i'])) 
                 # Add binding 
                 )


        Cortk = (
                 # Diffusion from plasma
                 + v['Dcort_pk'] * ((v['Cortp']/v['v_p']) - (v['Cortk']/v['v_k']))
                 # 11B-HSD2
                 - (x[5] * v['v_k'] * v['MWCort']
                 * ((v['Cortk']/(v['v_k'] * v['MWCort'])) 
                 /((v['Cortk']/(v['v_k'] * v['MWCort'])) + v['km2'])))
                )

        
        """
        Cortisonek = (
                      # Diffusion to plasma
                      + v['Dcort_pk'] *((v['Cortisonep']/v['v_p']) - (x[4]/v['v_k']))
                      # 11B-HSD2
                      + v['Vm2'] * (v['Cortk']/v['v_k']) 
                      /((v['Cortk']/v['v_k'])  
                      * (1/v['MWCort']) + v['km2']) * v['v_k']
                      * (v['MWsone']/v['MWCort'])
                     )
        """

        # Extra Constraint
        extra = (- (x[5] * v['v_k']
                 * ((v['Cortk']/(v['v_k'] * v['MWCort'])) 
                 /((v['Cortk']/(v['v_k'] * v['MWCort'])) + v['km2'])))
                 + (x[3] * v['MWCort'] * v['v_p'] 
                 * (v['Cortisonep']/(v['v_p'] * v['MWsone']))
                 / ((v['Cortisonep']/(v['v_p'] * v['MWsone'])) + v['km1']))
                )


        lhs = [CRF, ACTH, Cortp, extra, Cortisonep, Cortk]
        
        return lhs

    # Parameters
    full_dict = {**DEM.Paramsu}
    
    # Ensure State variable concentrations are right
    full_dict['CRH']         = 7.66 * full_dict['v_p'] 
    full_dict['ACTH']        = 21 * full_dict['v_p']
    full_dict['Cortp']       = 3.06 * full_dict['v_p']
    full_dict['Corti']       = 3.06 * full_dict['v_i']
    full_dict['Cortk']       = full_dict['Ald'] * full_dict['v_k'] * 0.1
    full_dict['Cortisonep']  = 3.06 * full_dict['v_p']
    full_dict['Cortisonek']  = (
                                (((full_dict['MWCort']/full_dict['MWsone']) 
                                * 3.06 + 3.06 - full_dict['Cortk']/full_dict['v_k'])
                                / (full_dict['MWCort']/full_dict['MWsone']))
                                * full_dict['v_k']
                               )


    # Create the solver function
    y = lambda x : sys(0, x, full_dict)
    # create initial Guess 
    guess = DEM.np.zeros(6)
    guess[0]  = 2000
    guess[1]  = 0.1
    guess[2]  = 0.01
    guess[3]  = 1
    guess[4]  = 20
    guess[5]  = 20

    y(guess)
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



    def sys(t, x, v): 
        CRF = (x[0]  * 
               # Positive hill feedback by cort_p
               (1 + v['zeta'] * (((v['Cortp']/v['v_p']) ** v['alpha'])
               / ((v['Cortp']/v['v_p']) ** v['alpha'] + v['c'] ** v['alpha'])) 
               # Negative hill feedback by cort_p
               - (v['psi'] * (((v['Cortp']/v['v_p']) ** v['gamma'])
               / ((v['Cortp']/v['v_p']) ** v['gamma'] + v['c3'] ** v['gamma']))))
               # Generic waste
               - v['w1'] * v['CRH']
              )



        ACTH = (x[1] * 
                # Negative hill feedback by cort_p
                (1 - v['rho'] * ((v['Cortp']/v['v_p']) ** v['alpha'])
                /((v['Cortp']/v['v_p']) ** v['alpha'] + v['c'] ** v['alpha']))
               * (v['CRH']/v['v_p']) * v['v_p']
               # Generic waste
               - v['w2'] * v['ACTH']
               )



        Cortp = (# stimulation
                 x[2] * (v['ACTH']/v['v_p']) * v['v_p']   
                 # diffusion to interstitium
                 - v['Dcort_pi'] * ((v['Cortp']/v['v_p']) - (v['Corti']/v['v_i']))
                 # diffusion to the kidney
                 - v['Dcort_pk'] * ((v['Cortp']/v['v_p']) - (v['Cortk']/v['v_k']))            
                 # Source from Vm1
                 +  (x[3] * v['MWCort'] * v['v_p'] 
                 * (v['Cortisonep']/(v['v_p'] * v['MWsone']))
                 / ((v['Cortisonep']/(v['v_p'] * v['MWsone'])) + v['km1']))
                 # Urine
                 - v['ucp'] * v['Cortp']
                )

        Cortisonep = (# Constant source 
                      x[4]
                      # Diffusion from kidney
                      - v['Dcort_pk'] *((v['Cortisonep']/v['v_p']) 
                      - (v['Cortisonek']/v['v_k']))
                      # Vm1 sink
                      +  (x[3] * v['MWCort'] * v['v_p'] 
                      * (v['Cortisonep']/(v['v_p'] * v['MWsone']))
                      / ((v['Cortisonep']/(v['v_p'] * v['MWsone'])) + v['km1']))
                      # Urine
                      - v['usp'] * v['Cortisonep']
                     )       
        
        
        # changed v['Cortisonep'] to high cort o   
        Corti = (
                 # diffusion
                 + v['Dcort_pi'] * ((v['Cortp']/v['v_p']) - (v['Corti']/v['v_i'])) 
                 # Add binding 
                 )


        Cortk = (
                 # Diffusion from plasma
                 + v['Dcort_pk'] * ((v['Cortp']/v['v_p']) - (v['Cortk']/v['v_k']))
                 # 11B-HSD2
                 - (x[5] * v['v_k'] * v['MWCort']
                 * ((v['Cortk']/(v['v_k'] * v['MWCort'])) 
                 /((v['Cortk']/(v['v_k'] * v['MWCort'])) + v['km2'])))
                )

        
        
        Cortisonek = (
                      # Diffusion to plasma
                      + v['Dcort_pk'] *((v['Cortisonep']/v['v_p']) 
                          - (v['Cortisonek']/v['v_k']))
                      # 11B-HSD2
                      + (x[5] * v['v_k'] * v['MWsone']
                      * ((v['Cortk']/(v['v_k'] * v['MWCort'])) 
                      /((v['Cortk']/(v['v_k'] * v['MWCort'])) + v['km2'])))
                      )
        

        # Extra Constraint

        lhs = [CRF, ACTH, Cortp, Corti, Cortk, Cortisonep, Cortisonek]

        return lhs
    
    print(sys(0, steady, full_dict))
    



    return 0

if __name__ == '__main__':
    main()
