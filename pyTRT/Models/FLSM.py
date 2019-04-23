

def FLSM_sim(ks,Rb,alpha,H,rb,z,T0,q,t):  
    
    """ Base FLSM function.

        ks = soil thermal conductivity (W/m.k)
        a = soil thermal diffusivity (m2/s)
        H = depth of the borehole (m)
        rb = borehole radius (m)
        z = depth to evaluate the temperature at (m)
        T0 = undisturbed ground temp (degC)
        q = power applied to the GHE (W/m)
        Rb = GHE thermal resistance (m.K/W)
        t = time (hours)
    """
    
    import numpy as np
    import pandas as pd
    from scipy.integrate import quad
    from scipy.special import erfc

    
    # Calculate alpha from Cpvs and ks internally
    a=alpha
    
    # Define the FLS integral function over the length of the borehole H    
    def H_int(h,rb,z,a,t):
        return((erfc(np.sqrt(rb**2 + (z - h)**2)/(2*np.sqrt(a*t)))/np.sqrt(rb**2 + (z - h)**2)) - 
                (erfc(np.sqrt(rb**2 + (z + h)**2)/(2*np.sqrt(a*t)))/np.sqrt(rb**2 + (z + h)**2)))
    
    def get_Temp(z,ks,a,rb,t,H):
        return(T0 + (q/(4*np.pi*ks))*quad(H_int, 0, H, args=(rb,z,a,t), epsabs=1.49e-8, epsrel=1.49e-8, limit=1000, points=[z])[0])
        
    # If t is a single variable, calculate the integral directly output a float
    # If t is a pd.Series, iterate over each value in s and output an array
    if isinstance(t, float)==True or isinstance(t, int)==True:
        t=t*60*60
        Temp=float(T0 + (q/(4*np.pi*ks))*quad(H_int, 0, H, args=(rb,z,a,t))[0] + (q*Rb))
        
    else:
        t = pd.Series(t)*3600 # convert t to seconds
        Temp = np.array(t.apply(lambda x: get_Temp(z,ks,a,rb,x,H))) + (q*Rb)

    return(Temp)

#Test 
#FLSM_sim(2.7,0.2,1.6e-7,50,0.075,50/2,18,2000/50,np.arange(1,100,1))

def FLSM_sim_intmean(ks,Rb,alpha,H,rb,T0,q,t):  
    
    """ Base FLSM function.

        ks = soil thermal conductivity (W/m.k)
        a = soil thermal diffusivity (m2/s)
        H = depth of the borehole (m)
        rb = borehole radius (m)
        z = depth to evaluate the temperature at (m)
        T0 = undisturbed ground temp (degC)
        q = power applied to the GHE (W/m)
        Rb = GHE thermal resistance (m.K/W)
        t = time (hours)
    """
    
    import numpy as np
    import pandas as pd
    from scipy.integrate import quad
    from scipy.special import erfc

    
    # Calculate alpha from Cpvs and ks internally
    a=alpha
    
    # Define the FLS integral function over the length of the borehole H    
    def H_int(h,rb,z,a,t):
        return((erfc(np.sqrt(rb**2 + (z - h)**2)/(2*np.sqrt(a*t)))/np.sqrt(rb**2 + (z - h)**2)) - 
                (erfc(np.sqrt(rb**2 + (z + h)**2)/(2*np.sqrt(a*t)))/np.sqrt(rb**2 + (z + h)**2)))
    
    # Now get_temp adds another layer of inegration of the temps along the depth
    def depth_int(z,ks,a,rb,t,H):
        return(T0 + (q/(4*np.pi*ks))*quad(H_int, 0, H, args=(rb,z,a,t), epsabs=1.49e-8, epsrel=1.49e-8, limit=1000, points=[z])[0])
        
    def get_Temp(ks,a,rb,t,H):
        return(quad(depth_int, 0, 30, args=(ks,a,rb,t,H))[0]/H)
        
    # If t is a single variable, calculate the integral directly output a float
    # If t is a pd.Series, iterate over each value in s and output an array
    if isinstance(t, float)==True or isinstance(t, int)==True:
        t=t*60*60
        Temp=float(get_Temp(ks,a,rb,t,H) + (q*Rb))
        
    else:
        t = pd.Series(t)*3600 # convert t to seconds
        Temp = np.array(t.apply(lambda x: get_Temp(ks,a,rb,x,H))) + (q*Rb)

    return(Temp)

#Test 
#FLSM_sim_intmean(2.7,0.2,1.6e-7,50,0.075,50/2,18,2000/50,np.arange(1,100,1))

def FLSM_sim_Lamarche(ks,Rb,alpha,H,rb,T0,q,t):  
    
    """ Improved ILSM function from Lamarche and Beauchamp 2007

        ks = soil thermal conductivity (W/m.k)
        a = soil thermal diffusivity (m2/s)
        H = depth of the borehole (m)
        rb = borehole radius (m)
        T0 = undisturbed ground temp (degC)
        q = power applied to the GHE (W/m)
        Rb = GHE thermal resistance (m.K/W)
        t = time (hours)
    """
    
    import numpy as np
    import pandas as pd
    from scipy.integrate import quad
    from scipy.special import erfc

    
    # Function to evaluate full improved FLSM at single time t
    def get_Temp(ks,Rb,alpha,H,rb,T0,q,t):
        # Set up constants
        B = rb/H
        Fo = (alpha*t)/rb**2
        t_star = 9*B**2*Fo
        #g = (3/2)*np.sqrt(t_star)
        g = H/np.sqrt(4*alpha*t)
        
        # Define Da and Db components
        Da = ((np.sqrt(B**2 + 1)*erfc(g*np.sqrt(B**2 + 1))) - (B*erfc(g*B)) - 
              ((np.exp(-g**2 * (B**2 + 1)) - np.exp(-g**2 * B**2))/(g*np.sqrt(np.pi))))
        #Da = quad(lambda z: erfc(q*z), B, np.sqrt(B**2 + 1))[0]
        
        Db = (np.sqrt(B**2 + 1)*erfc(g*np.sqrt(B**2 + 1)) - 
              0.5*(B*erfc(g*B) + np.sqrt(B**2 + 4)*erfc(g*np.sqrt(B**2 + 4))) - 
              (np.exp(-g**2 * (B**2 + 1)) - 0.5*(np.exp(-g**2 * B**2) + np.exp(-g**2 * (B**2 + 4))))/(g*np.sqrt(np.pi)))
        
        #Db = 0.5*(quad(lambda z: erfc(g*z), B, np.sqrt(B**2 + 1))[0] + quad(lambda z: erfc(g*z), np.sqrt(B**2 + 1), np.sqrt(B**2 + 4))[0])
        # Define A, B integrals, they are the same with different bounds
        def AB_int(z, B, g):
            return(erfc(g*z)/np.sqrt(z**2 - B**2))
        
        # Put it together in a g func
        Gfunc = (quad(AB_int, B, np.sqrt(B**2 + 1), args=(B, g))[0] - Da) - (quad(AB_int, np.sqrt(B**2 + 1), np.sqrt(B**2 + 4), args=(B, g))[0] + Db)
        
        return((q/(2*np.pi*ks))*Gfunc)
    
        
    # If t is a single variable, calculate the integral directly output a float
    # If t is a pd.Series, iterate over each value in s and output an array
    if isinstance(t, float)==True or isinstance(t, int)==True:
        t=t*60*60
        Temp=float( T0 + get_Temp(ks,Rb,alpha,H,rb,T0,q,t)  + q*Rb)
        
    else:
        t = pd.Series(t)*3600 # convert t to seconds
        Temp = T0 + np.array(t.apply(lambda x: get_Temp(ks,Rb,alpha,H,rb,T0,q,x))) + q*Rb

    return(Temp)

#Test 
#FLSM_sim_Lamarche(2.7,0.2,1.6e-7,50,0.075,18,2000/50,np.arange(1,100,10))

def FLSM_ksRb(t, temp, Q, tstart, tend, method, ks, Rb, alpha, H, rb, z, T0, verbose, fluid_calc):
    """ Function to fit the FLSM model to set of TRT data using the scipy-optimize minimize function.
        Inputs:
            t = vector of time steps recorded in test (hours)
            temp = vector of measured mean fluid temp data (degC)
            Q = Power applid to the GHE (W)
            tstart = time to start fitting region (hours)
            tend = time to finish fitting region (hours)
            method = optimization method to pass to minimize function (string)
                     (Defuat is L-BFGS-B, bounds are already set)
            ks = soil thermal conductivity (W/m.K)
            Rb = GHE thermal resistance (m.K/W)
            a = alpha, thermal diffusivity (m2/s)
            H = GHE depth (m)
            rb = GHE radius (m)
            T0 = undisturbed ground temperature (degC)
            z = depth at which temperature should be evaluated
            plot = boolean, whether to plot the result against the measured data or not
                   (Default is False)
    """
    from .FLSM import FLSM_sim, FLSM_sim_intmean
    from scipy.optimize import minimize 
    import numpy as np
    import pandas as pd
    
    def find_nearest_idx(array,value):
        idx = (np.abs(array-value)).argmin()
        return(idx)     
    # Set start and end time indeces
    t1=find_nearest_idx(t,tstart)
    t2=find_nearest_idx(t,tend)
    
    # Set fitting regions of data
    t_fit = t[t1:t2]
    temp_fit = temp[t1:t2]
    q = np.mean(Q[t1:t2])/H
    
    # Objective function using RMSE        
    if fluid_calc == 'mid':
        def RMSE(x0):
            pred = FLSM_sim(x0[0],x0[1]/10,alpha,H,rb,z,T0,q,t_fit)
            E = np.sqrt(sum((pred - temp_fit)**2)/len(pred))
            if verbose == True: print('RMSE = %f' % (E))
            return(E)
    elif fluid_calc == 'intmean':
        def RMSE(x0):
            pred = FLSM_sim_intmean(x0[0],x0[1]/10,alpha,H,rb,T0,q,t_fit)
            E = np.sqrt(sum((pred - temp_fit)**2)/len(pred))
            if verbose == True: print('RMSE = %f' % (E))
            return(E)
    elif fluid_calc == 'Lamarche':
        def RMSE(x0):
            pred = FLSM_sim_Lamarche(x0[0],x0[1]/10,alpha,H,rb,T0,q,t_fit)
            E = np.sqrt(sum((pred - temp_fit)**2)/len(pred))
            if verbose == True: print('RMSE = %f' % (E))
            return(E)
    
    # Now minimise this objective function
    x0=[ks,Rb*10] #initialization parameter values
    bounds = ((0,4), (0,10))
    opt_result=minimize(RMSE, x0, method=method, bounds=bounds)
    message = 'Message:' + str(opt_result.message)
    success = 'Optimization Successful? ' + str(opt_result.success)
    ks_result=opt_result.x[0]
    Rb_result=opt_result.x[1]/10
    print(success)
    print(message)
    print("Final RMSE = %0.5f\n" % (opt_result.fun))

    return({'ks':ks_result, 'Rb':Rb_result})
        
def FLSM_ksRbAlpha(t, temp, Q, tstart, tend, method, ks, Rb, alpha, H, rb, z, T0, verbose):
    """ Function to fit the FLSM model to set of TRT data using the scipy-optimize minimize function.
        Inputs:
            t = vector of time steps recorded in test (hours)
            temp = vector of measured mean fluid temp data (degC)
            Q = Power applid to the GHE (W)
            tstart = time to start fitting region (hours)
            tend = time to finish fitting region (hours)
            method = optimization method to pass to minimize function (string)
                     (Defuat is L-BFGS-B, bounds are already set)
            ks = soil thermal conductivity (W/m.K)
            Rb = GHE thermal resistance (m.K/W)
            a = alpha, thermal diffusivity (m2/s)
            H = GHE depth (m)
            rb = GHE radius (m)
            T0 = undisturbed ground temperature (degC)
            z = depth at which temperature should be evaluated

    """
    from .FLSM import FLSM_sim
    from scipy.optimize import minimize 
    import numpy as np
    import pandas as pd
    
    def find_nearest_idx(array,value):
        idx = (np.abs(array-value)).argmin()
        return(idx)     
    # Set start and end time indeces
    t1=find_nearest_idx(t,tstart)
    t2=find_nearest_idx(t,tend)
    
    # Set fitting regions of data
    t_fit = t[t1:t2]
    temp_fit = temp[t1:t2]
    q = np.mean(Q[t1:t2])/H
    
    # Objective function using RMSE        
        
    def RMSE(x0):
        pred = FLSM_sim(x0[0],x0[1]/10,x0[2]/1e6,H,rb,z,T0,q,t_fit)
        E = np.sqrt(sum((pred - temp_fit)**2)/len(pred))
        if verbose == True: print('RMSE = %f' % (E))
        return(E)
    
    # Now minimise this objective function
    x0=[ks,Rb*10,alpha*1e6] #initialization parameter values
    bounds = ((0,4), (0,10), (0,10))
    opt_result=minimize(RMSE, x0, method=method, bounds=bounds)
    message = 'Message:' + str(opt_result.message)
    success = 'Optimization Successful? ' + str(opt_result.success)
    ks_result=opt_result.x[0]
    Rb_result=opt_result.x[1]/10
    alpha_result=opt_result.x[2]/1e6
    print(success)
    print(message)
    print("Final RMSE = %0.5f\n" % (opt_result.fun))

    return({'ks':ks_result, 'Rb':Rb_result, 'alpha': alpha_result})