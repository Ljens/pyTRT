""" 
Script containing versions of the Infinite Cylindrical Source Model based on: 

    R.A. Beier, Soil Thermal Conductivity Tests, Oklahoma State University, 2008. 
    http://www.igshpa.okstate.edu/research/papers/tc_testing_copyright.pdf.

"""

def ICSM_sim(ks,Rb,Cpvs,rb,q,T0,t):
    """    
    Based on original model formulation from Ingersol et al (1954).    
    Integral and bessel functions are solved numerically.
    
    Parameters:
    ks:     ground thermal conductivity (W/mK)
    Rb:     GHE thermal resistance (mK/W)
    Cpvs:   Soil Volumetric specific heat capacity (J/m3K)
    rb:     GHE radius
    q:      Heating power rate (W/m)
    T0:     Undisturbed Ground Temperature
    t:      time (hours)
    
    """
    import numpy as np
    import pandas as pd
    from scipy.integrate import quad
    
    # Import Bessel functions required in ICSM formulation
    from scipy.special import j0, j1, y0, y1
    
    #Set Constants
    r0 = rb
    mu = rb/r0

    def integrand(x, Fo):    
        return(((np.exp(-x**2 * Fo)-1)/(j1(x)**2 + y1(x)**2)) * ((j0(mu*x)*y1(x) - j1(x)*y0(mu*x))/x**2))

    #return(T0 + q/(np.pi**2 * ks) * quad(integrand, 0, np.inf, args=(Fo))[0] + q*Rb)

    if isinstance(t, float)==True or isinstance(t, int)==True:
        Fo = (ks*t)/(Cpvs*rb**2)*3600
        return(T0 + q/(np.pi**2 * ks) * quad(integrand, 0, np.inf, args=(Fo))[0] + q*Rb)      
            
    else:
        result=[]
        for i in t:
            Fo = (ks*i)/(Cpvs*rb**2)*3600
            result.append(T0 + q/(np.pi**2 * ks) * quad(integrand, 0, np.inf, args=(Fo))[0] + q*Rb)
        return(result)


def ICSM_ksRb(ks, Rb, Cpvs, t, temp, Power, H, T0, rb, t_start, t_end, method, verbose):
    """ Fit the ILS to measured data by curve fitting two parameters:
            ks - Soil Thermal Conductivity (W/m.K)
            Rb - Borehole Resistance (m.K/W)
            These are defined above as the initial guess.
        
        ILSM formula is from Hemmingway et al 2012. Could be updated to match a more reputable paper.
            
        Also need to define:
            t - time (Hours)
            temp - TRT temp response curve (degC)
            Power - TRT measured heating power curve (W)
            Cpvs - the volumetric heat capacity of the soil (J/m3.K)
            H - vertical length of pipe (m)
            T0 - initial ground temp (degC)
            rb - borehole radius (m)
            start_hour - the lower bound of the regression region in hours
            end_hour - the upper bound of the regression region in hours
            method - the optimzation algorithm to use, passed to scikit.minimize function
    """

    from .ICSM import ICSM_sim
    from scipy.optimize import minimize
    import numpy as np
    import pandas as pd
    
    # Function to find the nearest index to a value in an array    
    def find_nearest_idx(array,value):
        idx = (np.abs(array-value)).argmin()
        return(idx)    
        
    # Convert the start and end hour values to index values
    t1 = find_nearest_idx(t,t_start)
    t2 = find_nearest_idx(t,t_end)
    q = np.mean(Power[t1:t2])/H

    # Sum of Least Squares objective function to be minimised by the optimization algorithm      
    def RMSE(param):
        pred = ICSM_sim(param[0],param[1]/10,Cpvs,rb,q,T0,t[t1:t2])
        E = np.sqrt(np.mean((pred - temp[t1:t2])**2))
        if verbose == True: print('RMSE = %f' % (E))
        return(E)
    #RMSE([2.7,0.2*10])
    
    # Optimization using minimize function
    param=np.array([ks,Rb*10]) #array of rough guess of input parameters
    bounds=((0,4), (0,10))
    loc_min=minimize(RMSE, param, method=method, bounds=bounds)
    ks_result=loc_min.x[0]
    Rb_result=loc_min.x[1]/10
    print("\n" + "{}".format(loc_min.message))
    print("Final RMSE = %0.5f\n" % (loc_min.fun))

    return(dict(ks=ks_result,Rb=Rb_result))

def ICSM_ksRbCpvs(ks, Rb, alpha, t, temp, Power, H, T0, rb, t_start, t_end, method, verbose):
    """ Fit the ICSM to measured data by curve fitting three parameters:
            ks - Soil Thermal Conductivity (W/m.K)
            Rb - Borehole Resistance (m.K/W)
            Cpvs
            These are defined at function call as the initial guess.

            
        Also need to define:
            t - time (seconds)
            temp - TRT temp response curve (degC)
            Power - TRT measured heating power curve (W)
            Cpvs - the volumetric heat capacity of the soil (J/m3.K)
            H - vertical length of pipe (m)
            T0 - initial ground temp (degC)
            rb - borehole radius (m)
            start_hour - the lower bound of the regression region in hours
            end_hour - the upper bound of the regression region in hours
            method - the optimzation algorithm to use, passed to scikit.minimize function
    """
    from .ICSM import ICSM_sim  
    from scipy.optimize import minimize
    import numpy as np
    import pandas as pd
    
    # Function to find the nearest index to a value in an array    
    def find_nearest_idx(array,value):
        idx = (np.abs(array-value)).argmin()
        return(idx)    
        
    # Convert the start and end hour values to index values
    t1 = find_nearest_idx(t,t_start)
    t2 = find_nearest_idx(t,t_end)
    q = np.mean(Power[t1:t2])/H
    gama=0.5772156649 #eulers number

    # Sum of Least Squares objective function to be minimised by the optimization algorithm              
    def RMSE(param):
        pred = ICSM_sim(param[0],param[1]/10,param[2],rb,q,T0,t[t1:t2])
        E = np.sqrt(np.sum((pred - temp[t1:t2])**2)/len(pred))
        if verbose == True: print('RMSE = %f' % (E))
        return(E)
    #RMSE([ks,Rb*10, alpha*1e6])
    
    # Optimization using minimize function
    param=np.array([ks,Rb*10, Cpvs])
    bounds=((0,4), (0,10), (0, None))
    loc_min=minimize(RMSE, param, method=method, bounds=bounds)
    ks_result=loc_min.x[0]
    Rb_result=loc_min.x[1]/10
    alpha_result=loc_min.x[2]
    print("\n" + "{}".format(loc_min.message))
    print("Final RMSE = %0.5f\n" % (loc_min.fun))

    return(dict(ks=ks_result, Rb=Rb_result, alpha=alpha_result))