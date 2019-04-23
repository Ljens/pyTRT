
import numpy as np
import pandas as pd

def ILSM_sim(ks,Rb,alpha,q,rb,t,T0):
    
    """ Base ILSM function.
        Takes all required inputs including a single time value in seconds are returns a corresponding value of temperature.
    
        ks = soil thermal conducitivity (W/mK)
        Rb = borehole thermal resistance (m.K/W)
        rb = borehole radius (m)
        alpha = soil thermal diffusivity (m/s^2)
        q = heating/cooling power (W/m)
        t = time/s at which temperature occurs (hrs). Can be int, float, list, or an 1D iterable (np.array, pd.Series)
        T0 = undisturbed ground temperature (degC)
        
    """
    import numpy as np
    import pandas as pd
    
    gama=0.5772156649 # gamma constant is Euler's number
    
    if isinstance(t, float)==True or isinstance(t, int)==True:
        s=t*3600 # convert t to seconds
        return(q/(4*np.pi*ks)*np.log(s) + q*(Rb + (1/(4*np.pi*ks))*(np.log((4*alpha)/(rb**2))-gama)) + T0)        

    elif isinstance(t, list)==True:
        s = [x * 3600 for x in t] # convert t to seconds
        result=[]
        for i in s:
            result.append(q/(4*np.pi*ks)*np.log(i) + q*(Rb + (1/(4*np.pi*ks))*(np.log((4*alpha)/(rb**2))-gama)) + T0)
        return(result)
            
    else:
        s=t*3600 # convert t to seconds
        result=[]
        for i in s:
            result.append(q/(4*np.pi*ks)*np.log(i) + q*(Rb + (1/(4*np.pi*ks))*(np.log((4*alpha)/(rb**2))-gama)) + T0)
        return(result)
        
def ILSM_EI_sim(ks,Rb,alpha,q,rb,t,T0):
    
    """ Base ILSM function using the EI() function instead of ln()
        Takes all required inputs including a single time value in seconds are returns a corresponding value of temperature.
    
        ks = soil thermal conducitivity (W/mK)
        Rb = borehole thermal resistance (m.K/W)
        rb = borehole radius (m)
        alpha = soil thermal diffusivity (m/s^2)
        q = heating/cooling power (W/m)
        t = time/s at which temperature occurs (hrs). Can be int, float, list, or an 1D iterable (np.array, pd.Series)
        T0 = undisturbed ground temperature (degC)
        
    """
    import numpy as np
    import pandas as pd
    from scipy.special import expi 
    from scipy.integrate import quad
    from mpmath import ei

    if isinstance(t, float)==True or isinstance(t, int)==True:
        s=t*3600 # convert t to seconds
        return(T0 + (q/(4*np.pi*ks))*-expi(-rb**2/(4*alpha*s)) + q*Rb)
        #return(T0 + (q/(4*np.pi*ks))*-ei(-rb**2/(4*alpha*s)) + q*Rb)  
        #return(T0 + (q/(4*np.pi*ks))*-quad(lambda t: np.exp(t)/t, -np.inf, -rb**2/(4*alpha*s))[0] + q*Rb, epsabs=1.49e-8, epsrel=1.49e-8, limit=1000, points=[t])        

    elif isinstance(t, list)==True:
        s = [x * 3600 for x in t] # convert t to seconds
        result=[]
        for i in s:
            result.append(T0 + (q/(4*np.pi*ks))*-expi(-rb**2/(4*alpha*i)) + q*Rb)
            #result.append(T0 + (q/(4*np.pi*ks))*-ei(-rb**2/(4*alpha*i)) + q*Rb)
            #result.append(T0 + (q/(4*np.pi*ks))*-quad(lambda t: np.exp(t)/t, -np.inf, -rb**2/(4*alpha*i))[0] + q*Rb, epsabs=1.49e-8, epsrel=1.49e-8, limit=1000, points=[t])
            
        return(result)
            
    else:
        s=t*3600 # convert t to seconds
        result=[]
        for i in s:
            result.append(T0 + (q/(4*np.pi*ks))*-expi(-rb**2/(4*alpha*i)) + q*Rb)
            #result.append(T0 + (q/(4*np.pi*ks))*-ei(-rb**2/(4*alpha*i)) + q*Rb)
            #result.append(T0 + (q/(4*np.pi*ks))*-quad(lambda t: np.exp(t)/t, -np.inf, -rb**2/(4*alpha*i))[0] + q*Rb,  epsabs=1.49e-8, epsrel=1.49e-8, limit=1000, points=[t])
        return(result)
    
    
# Test Vlaues
#ILSM_sim(3,0.2,1.4e-7,70,0.075,pd.Series([1,2,3]),18)


def ILSM_linreg(t,Temp,Power,t_start,t_end,H,T0,rb,Cpvs, mode='heating'):
    """ 
    Take a set of measured data and fit a linear model to it
    using the ILSM to calculate thermal conductivity and borehole resistance.
    
    Parameters:
    t:          time (hours)
    Temp:       Mean Fluid Temperature (°C)
    t_start:    starting time of regression window (hours)
    t_end:      end time of regression window (hours)
    mode:       either 'heating' or 'cooling'
    
    """
    from sklearn import linear_model
    import numpy as np
    import pandas as pd
    
    def find_nearest_idx(array,value):
        idx = (np.abs(array-value)).idxmin()
        return(idx)   
    
    if mode == 'heating':
        mode_sign = 1
    elif mode == 'cooling':
        mode_sign = -1
    else:
        print('ERROR: select mode as heating or cooling')
        
    # Convert t in hours to the log space
    t_log=np.log(t)
    
    # Convert times to indices
    t1=find_nearest_idx(t_log,np.log(t_start))
    t2=find_nearest_idx(t_log,np.log(t_end))
    
    # Use either an array for power or s single value
    if isinstance(Power, float) == True:
        Q=Power
    else:
        Q=np.mean(Power[t1:t2])
    

    # Calculate ks and Rb
    lm=linear_model.LinearRegression(fit_intercept=True) #create a linear model object
    lm.fit(np.array(t_log).reshape(-1, 1)[t1:t2],np.array(Temp).reshape(-1, 1)[t1:t2]) #fit it to a subset of the data
    m=mode_sign * float(lm.coef_) #the slope of the line
    c=float(lm.predict(np.array(0).reshape(1, -1))) #the y-intercept of the line
    ks=(Q/(4*np.pi*m*H))
    alpha=ks/(Cpvs)*3600 # Thermal diffusivity in m2/hr
    Rb=(1/(4*np.pi*ks)*(abs(c - T0)/m - np.log((4*alpha)/(1.78*rb**2))))
    return (dict(ks=ks,Rb=Rb, c=c, m=m))
    
def ILSM_linreg_uncert(t,Temp,Q,t_start,t_end,H,T0,rb,Cpvs, mode='heating'):
    """ 
    Take a set of measured data and fit a linear model to it
    using the ILSM to calculate thermal conductivity and borehole resistance.
    
    Parameters:
    t:          time (hours)
    Temp:       Mean Fluid Temperature (°C)
    t_start:    starting time of regression window (hours)
    t_end:      end time of regression window (hours)
    mode:       either 'heating' or 'cooling'
    
    """
    from sklearn import linear_model
    import numpy as np
    import pandas as pd
    
    def find_nearest_idx(array,value):
        idx = (np.abs(array-value)).argmin()
        return(idx)   
    
    if mode == 'heating':
        mode_sign = 1
    elif mode == 'cooling':
        mode_sign = -1
    else:
        print('ERROR: select mode as heating or cooling')
    
    
    # Convert t in hours to the log space
    t_log=np.log(t)
    
    # Convert times to indices
    t1=find_nearest_idx(t_log,np.log(t_start))
    t2=find_nearest_idx(t_log,np.log(t_end))
    
    # Calculate ks and Rb
    lm=linear_model.LinearRegression(fit_intercept=True) #create a linear model object
    lm.fit(np.array(t_log)[t1:t2,np.newaxis],np.array(Temp)[t1:t2]) #fit it to a subset of the data
    m=mode_sign * float(lm.coef_) #the slope of the line
    c=float(lm.predict(np.array(0).reshape(1, -1))) #the y-intercept of the line
    ks=(Q/(4*np.pi*m*H))
    alpha=ks/(Cpvs)*3600 # Thermal diffusivity in m2/hr
    #Rb=(1/(4*np.pi*ks)*((c - T0)/m - np.log((4*alpha)/(1.78*rb**2)))) 
    Rb=(1/(4*np.pi*ks)*(abs(c - T0)/m - np.log((4*alpha)/(1.78*rb**2)))) 
    return (dict(ks=ks,Rb=Rb, c=c, m=m))


def ILSM_linreg_multi(t,Temp,Power,t_start1,t_start2,t_end,res,H,T0,rb,Cpvs):
    """ 
        Take a set of measured data and fit a linear model to it
        using the ILSM to calculate thermal conductivity and borehole resistance.
        t = time in hours
        Temp = TRT temp response curve
        t_start1 = starting time of first regression window (hours)
        t_start2 = starting time of last regression window (hours)
        res = how many indices to increase the regression window after each run by
    """
    from sklearn import linear_model
    import numpy as np
    import pandas as pd
    
    def find_nearest_idx(array,value):
        idx = (np.abs(array-value)).argmin()
        return(idx)       

    # Convert t in hours to the log space
    t_log=np.log(t)
    
    #t_start1=0.25
    #t_start2=1.9
    # define the upper and lower time bounds in terms of indices
    hours = np.arange(t_start1,t_start2,res) # in hours
    idxs = [] 
    for i in hours:
        idxs.append(find_nearest_idx(t_log,np.log(i)))

    #t_end=len(t_log)
    t_end = find_nearest_idx(t, t_end)
    
    ks_out=[]
    Rb_out=[]
    m_out=[]
    c_out=[]
    for i in idxs:
        lm=linear_model.LinearRegression(fit_intercept=True) #create a linear model object
        lm.fit(np.array(t_log)[i:t_end,np.newaxis],np.array(Temp)[i:t_end]) #fit it to a subset of the data
        m=float(lm.coef_) #the slope of the line
        c=float(lm.predict(np.array(0).reshape(1, -1))) #the y-intercept of the line
        Q=np.mean(Power[i:t_end])
        ks=(Q/(4*np.pi*m*H))
        alpha=ks/(Cpvs)*3600 # Thermal diffusivity in m2/hr
        Rb=(1/(4*np.pi*ks)*(abs(c - T0)/m - np.log((4*alpha)/(1.78*rb**2)))) 
        ks_out.append(ks)
        Rb_out.append(Rb)
        m_out.append(m)
        c_out.append(c)
        
    df=pd.DataFrame()
    df['Hours']=t[idxs]
    df['ks']=ks_out
    df['Rb']=Rb_out
    df['m']=m_out
    df['c']=c_out
    df.reset_index(inplace=True,drop=True)
    return(df)
        
# Test
# ILS_linreg_multi(t,Temp,Power,t_start1,t_start2,res,H,T0,rb,Cpvs)  
#ILS_linreg_multi(heat_needle.data.Hours,heat_needle.data.Sensor1,heat_needle.data.Power,t_start1,t_start2,res,H,T0,rb,Cpvs)
        

def ILSM_linreg_multi_reverse(t,Temp,Power,tstart,res,H,T0,rb,Cpvs):
    """ 
        Take a set of measured data and fit a linear model to it
        using the ILSM to calculate thermal conductivity and borehole resistance.
        t = time in hours
        Temp = TRT temp response curve
        t_start1 = starting time of first regression window (hours)
        t_start2 = starting time of last regression window (hours)
        res = how many indices to increase the regression window after each run by
    """
    from sklearn import linear_model
    import numpy as np
    import pandas as pd
    
    def find_nearest_idx(array,value):
        idx = (np.abs(array-value)).argmin()
        return(idx)       

    # Convert t in hours to the log space
    t_log=np.log(t)

    #t_end=len(t_log)
    tend = max(t)
    
    #t_start1=0.25
    #t_start2=1.9
    # define the upper and lower time bounds in terms of indices
    hours = np.arange(tstart,tend,res) # in hours
    idxs = []
    for i in hours:
        idxs.append(find_nearest_idx(t_log,np.log(i)))

    
    ks_out=[]
    Rb_out=[]
    for i in idxs[1:]:
        lm=linear_model.LinearRegression(fit_intercept=True) #create a linear model object
        lm.fit(np.array(t_log)[idxs[0]:i,np.newaxis],np.array(Temp)[idxs[0]:i]) #fit it to a subset of the data
        m=float(lm.coef_) #the slope of the line
        c=float(lm.predict(np.array(0).reshape(1, -1))) #the y-intercept of the line
        Q=np.mean(Power[idxs[0]:i])
        ks=(Q/(4*np.pi*m*H))
        alpha=ks/(Cpvs)*3600 # Thermal diffusivity in m2/hr
        Rb=(1/(4*np.pi*ks)*(abs(c - T0)/m - np.log((4*alpha)/(1.78*rb**2)))) 
        ks_out.append(ks)
        Rb_out.append(Rb)
        #print(idxs[0])        
        #print(i)
        
    df=pd.DataFrame()
    df['Hours']=t[idxs[1:]]
    df['ks']=ks_out
    df['Rb']=Rb_out
    df.reset_index(inplace=True,drop=True)
    return(df)
        
# Test
#ILS_linreg_multi_Reverse(t,Temp,Power,10,5,H,T0,rb,Cpvs)    

def ILSM_ksRb(ks, Rb, alpha, t, temp, Power, H, T0, rb, t_start, t_end, method, verbose):
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

    from .ILSM import ILSM_sim
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
        pred = ILSM_sim(param[0],param[1]/10,alpha,q,rb,t[t1:t2],T0)
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
    
#ILS_2par(1, 0.5, alpha, t, temp, Power, H, T0, rb, t_start, t_end, 'L-BFGS-B')
    
def ILSM_EI_ksRb(ks, Rb, alpha, t, temp, Power, H, T0, rb, t_start, t_end, method, verbose):
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

    from .ILSM import ILSM_EI_sim
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
        pred = ILSM_EI_sim(param[0],param[1]/10,alpha,q,rb,t[t1:t2],T0)
        E = np.sqrt(np.mean((pred - temp[t1:t2])**2))
        if verbose == True: print('RMSE = %f' % (E))
        return(E)
    #RMSE([2.7,0.2*10])
    
    # Optimization using minimize function
    param=np.array([ks,Rb*10]) #array of rough guess of input parameters
    bounds=((0,4), (0,10))
    loc_min=minimize(RMSE, param, method=method, bounds=bounds, options={'gtol': 1e-06, 'ftol':1e-10, 'maxiter':20000, 'maxfun': 20000})
    ks_result=loc_min.x[0]
    Rb_result=loc_min.x[1]/10
    print("\n" + "{}".format(loc_min.message))
    print("Final RMSE = %0.5f\n" % (loc_min.fun))

    return(dict(ks=ks_result,Rb=Rb_result))
        
def ILSM_ksRbAlpha(ks, Rb, alpha, t, temp, Power, H, T0, rb, t_start, t_end, method, verbose):
    """ Fit the ILS to measured data by curve fitting two parameters:
            ks - Soil Thermal Conductivity (W/m.K)
            Rb - Borehole Resistance (m.K/W)
            These are defined above as the initial guess.
        
        ILSM formula is from Hemmingway et al 2012. Could be updated to match a more reputable paper.
            
        Also need to define:
            Seconds - time (seconds)
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
    from .ILSM import ILSM_sim  
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
        pred = ILSM_sim(param[0],param[1]/10,param[2]/1e6,q,rb,t[t1:t2],T0)
        E = np.sqrt(np.sum((pred - temp[t1:t2])**2)/len(pred))
        if verbose == True: print('RMSE = %f' % (E))
        return(E)
    #RMSE([ks,Rb*10, alpha*1e6])
    
    # Optimization using minimize function
    param=np.array([ks,Rb*10, alpha*1e6])
    bounds=((0,4), (0,10), (0, 10))
    loc_min=minimize(RMSE, param, method=method, bounds=bounds)
    ks_result=loc_min.x[0]
    Rb_result=loc_min.x[1]/10
    alpha_result=loc_min.x[2]/1e6
    print("\n" + "{}".format(loc_min.message))
    print("Final RMSE = %0.5f\n" % (loc_min.fun))

    return(dict(ks=ks_result, Rb=Rb_result, alpha=alpha_result))
    
#ILS_3par(ks, Rb, alpha, t, temp, Power, H, T0, rb, t_start, t_end, 'Powell')