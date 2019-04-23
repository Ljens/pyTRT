import numpy as np
import pandas as pd


def Gfunc_sim(ks,Cpvs,Rc,q,rb,t,T0,n,SDR,F,AR,pipe_placement='PNE'):
    """
    ks = soil thermal conductivity
    Cpvs = soil vollumetric specific heat capacity
    Rb = pile thermal resistance
    q = applied power (W/m)
    rb = pile radius
    t = time (s)
    T0 = Undisturbed ground temperature
    n = number of pipes in pile
    SDR = pipe standard dimensional ratio
    F = flow rate (L/s)
    AR = aspect ration, must be on of: 15, 25, 33 or 50
    pipe_placement = must be: 'PNE' or 'CP'
    """
    
    ##############################################################################
    #n = 6 #number of pipes
    #SDR = 13.6 
    #F = 0.25 #flow rate (L/s)

    do = 0.025 #pipe outer diameter
    p_thickness = do/SDR #pipe wall thickness (m)
    di = do-2*p_thickness #inner pipe diameter
    ri = di/2
    ro = do/2

    u_ave = F/(np.pi*ri**2) #average fluid velocity (m/s)

    rhof = 1000 #fluid density (kg/m3)
    mu = 0.001 #dynamic viscosity (Pa.s)
    cpf = 4180 #fluid specific heat capacity (J/kg.K)
    kf = 0.582 #fluid thermal conductivity (W/m.K)
    kp = 0.4 #pipe thermal conductivity (W/m.K)
    e = 0.0015 #pipe absolute roughness (mm)

    Re = (rhof*u_ave*di)/mu #Reynolds number
    Pr = (cpf*mu)/kf # Prandtl Number
    hi = (0.023*Re**0.8*Pr**0.35*kf)/(2*ri) #heat transfer coefficient

    # The Convective fluid resistance is..
    Rp_conv=1/(2*n*np.pi*ri*hi)
    #print(Rp_conv)

    # The conductive pipe wall resistance
    Rp_cond=np.log(ro/ri)/(2*n*np.pi*kp)
    #print(Rp_cond)

    ##############################################################################    
    # Import Gfunc coefs from file
    #Pile
    Gp_coefs = pd.read_csv(fileOut + 'Gfunc_pile_coefs.csv')
    #Gp_coefs = ',Coefficient,CP_Lower,CP_Upper,PNE_Lower,PNE_Upper\n0,a,-0.000101,3.5499999999999996e-05,-1.4400000000000001e-05,-2.99e-05\n1,b,-0.000234,6.02e-05,1.28e-05,-8.04e-06\n2,c,0.003037,-0.000603,0.000953,0.0008609999999999999\n3,d,0.0018030000000000001,0.0013,0.000131,-0.0011300000000000001\n4,e,-0.04339,-0.00744,-0.0245,-0.0109\n5,f,0.1029,0.0256,0.0757,0.0479\n6,g,0.9095,0.9690000000000001,0.9209999999999999,0.9390000000000001\n'
    #Gp_coefs = pd.read_csv(Gp_coefs)
    Gp_coefs.set_index('Coefficient', inplace=True)
    
    #Ground
    Gg_coefs = pd.read_csv(fileOut + 'Gfunc_ground_coefs2.csv')
    #Gg_coefs = ',Coefficient,AR15_Upper,AR25_Upper,AR33_Upper,AR50_Upper,AR15_Lower,AR25_Lower,AR33_Lower,AR50_Lower\n0,a,-4.837e-07,-3.796e-07,-2.192e-07,-5.1420000000000006e-08,2.68e-07,-6.107999999999999e-07,-8.984e-07,-8.741e-08\n1,b,6.5970000000000005e-06,6.441e-06,4.311e-06,8.756e-07,-1.306e-05,1.83e-05,3.137e-05,8.242999999999999e-06\n2,c,6.592e-05,4.129e-05,2.939e-05,3.233e-05,0.0001827,-0.0001942,-0.0003894,-0.0001835\n3,d,-0.0008843,-0.0008687,-0.0007328,-0.0005292,-9.15e-05,0.001366,0.002361,0.001894\n4,e,-0.004678,-0.003276,-0.002647,-0.00279,-0.01434,-0.01275,-0.01257,-0.01375\n5,f,0.03975,0.04415,0.0443,0.04284,0.05634,0.04932,0.04341,0.04905\n6,g,0.3018,0.3071,0.3076,0.3144,0.3722,0.3863,0.3928,0.3997\n7,h,0.5715,0.5819,0.5861,0.59,0.3989,0.4173,0.4245,0.4267\n'
    #Gg_coefs = pd.read_csv(Gg_coefs)
    Gg_coefs.set_index('Coefficient', inplace=True)
    
    ##############################################################################
    # calcualte Alpha
    alpha_g = ks/(Cpvs)

    Gp_Upper = (Gp_coefs.loc['a'][pipe_placement + '_Upper']*(np.log((alpha_g*t)/rb**2))**6 +
                    Gp_coefs.loc['b'][pipe_placement + '_Upper']*(np.log((alpha_g*t)/rb**2))**5 +
                    Gp_coefs.loc['c'][pipe_placement + '_Upper']*(np.log((alpha_g*t)/rb**2))**4 +
                    Gp_coefs.loc['d'][pipe_placement + '_Upper']*(np.log((alpha_g*t)/rb**2))**3 +
                    Gp_coefs.loc['e'][pipe_placement + '_Upper']*(np.log((alpha_g*t)/rb**2))**2 +
                    Gp_coefs.loc['f'][pipe_placement + '_Upper']*(np.log((alpha_g*t)/rb**2)) +
                    Gp_coefs.loc['g'][pipe_placement + '_Upper'])

    Gp_Lower = (Gp_coefs.loc['a'][pipe_placement + '_Lower']*(np.log((alpha_g*t)/rb**2))**6 +
                    Gp_coefs.loc['b'][pipe_placement + '_Lower']*(np.log((alpha_g*t)/rb**2))**5 +
                    Gp_coefs.loc['c'][pipe_placement + '_Lower']*(np.log((alpha_g*t)/rb**2))**4 +
                    Gp_coefs.loc['d'][pipe_placement + '_Lower']*(np.log((alpha_g*t)/rb**2))**3 +
                    Gp_coefs.loc['e'][pipe_placement + '_Lower']*(np.log((alpha_g*t)/rb**2))**2 +
                    Gp_coefs.loc['f'][pipe_placement + '_Lower']*(np.log((alpha_g*t)/rb**2)) +
                    Gp_coefs.loc['g'][pipe_placement + '_Lower'])

    Gg_Upper = (Gg_coefs.loc['a']['AR%d_Upper' % AR]*(np.log((alpha_g*t)/rb**2))**7 +
                    Gg_coefs.loc['b']['AR%d_Upper' % AR]*(np.log((alpha_g*t)/rb**2))**6 +
                    Gg_coefs.loc['c']['AR%d_Upper' % AR]*(np.log((alpha_g*t)/rb**2))**5 +
                    Gg_coefs.loc['d']['AR%d_Upper' % AR]*(np.log((alpha_g*t)/rb**2))**4 +
                    Gg_coefs.loc['e']['AR%d_Upper' % AR]*(np.log((alpha_g*t)/rb**2))**3 +
                    Gg_coefs.loc['f']['AR%d_Upper' % AR]*(np.log((alpha_g*t)/rb**2))**2 +
                    Gg_coefs.loc['g']['AR%d_Upper' % AR]*(np.log((alpha_g*t)/rb**2)) +
                    Gg_coefs.loc['h']['AR%d_Upper' % AR])

    Gg_Lower = (Gg_coefs.loc['a']['AR%d_Lower' % AR]*(np.log((alpha_g*t)/rb**2))**7 +
                    Gg_coefs.loc['b']['AR%d_Lower' % AR]*(np.log((alpha_g*t)/rb**2))**6 +
                    Gg_coefs.loc['c']['AR%d_Lower' % AR]*(np.log((alpha_g*t)/rb**2))**5 +
                    Gg_coefs.loc['d']['AR%d_Lower' % AR]*(np.log((alpha_g*t)/rb**2))**4 +
                    Gg_coefs.loc['e']['AR%d_Lower' % AR]*(np.log((alpha_g*t)/rb**2))**3 +
                    Gg_coefs.loc['f']['AR%d_Lower' % AR]*(np.log((alpha_g*t)/rb**2))**2 +
                    Gg_coefs.loc['g']['AR%d_Lower' % AR]*(np.log((alpha_g*t)/rb**2)) +
                    Gg_coefs.loc['h']['AR%d_Lower' % AR])

    """ Old version with .ix, deprecated in pandas V20.0 >
    Gp_Upper = (Gp_coefs.ix['a'][pipe_placement + '_Upper']*(np.log((alpha_g*t)/rb**2))**6 +
                    Gp_coefs.ix['b'][pipe_placement + '_Upper']*(np.log((alpha_g*t)/rb**2))**5 +
                    Gp_coefs.ix['c'][pipe_placement + '_Upper']*(np.log((alpha_g*t)/rb**2))**4 +
                    Gp_coefs.ix['d'][pipe_placement + '_Upper']*(np.log((alpha_g*t)/rb**2))**3 +
                    Gp_coefs.ix['e'][pipe_placement + '_Upper']*(np.log((alpha_g*t)/rb**2))**2 +
                    Gp_coefs.ix['f'][pipe_placement + '_Upper']*(np.log((alpha_g*t)/rb**2)) +
                    Gp_coefs.ix['g'][pipe_placement + '_Upper'])

    Gp_Lower = (Gp_coefs.ix['a'][pipe_placement + '_Lower']*(np.log((alpha_g*t)/rb**2))**6 +
                    Gp_coefs.ix['b'][pipe_placement + '_Lower']*(np.log((alpha_g*t)/rb**2))**5 +
                    Gp_coefs.ix['c'][pipe_placement + '_Lower']*(np.log((alpha_g*t)/rb**2))**4 +
                    Gp_coefs.ix['d'][pipe_placement + '_Lower']*(np.log((alpha_g*t)/rb**2))**3 +
                    Gp_coefs.ix['e'][pipe_placement + '_Lower']*(np.log((alpha_g*t)/rb**2))**2 +
                    Gp_coefs.ix['f'][pipe_placement + '_Lower']*(np.log((alpha_g*t)/rb**2)) +
                    Gp_coefs.ix['g'][pipe_placement + '_Lower'])

    Gg_Upper = (Gg_coefs.ix['a']['AR%d_Upper' % AR]*(np.log((alpha_g*t)/rb**2))**7 +
                    Gg_coefs.ix['b']['AR%d_Upper' % AR]*(np.log((alpha_g*t)/rb**2))**6 +
                    Gg_coefs.ix['c']['AR%d_Upper' % AR]*(np.log((alpha_g*t)/rb**2))**5 +
                    Gg_coefs.ix['d']['AR%d_Upper' % AR]*(np.log((alpha_g*t)/rb**2))**4 +
                    Gg_coefs.ix['e']['AR%d_Upper' % AR]*(np.log((alpha_g*t)/rb**2))**3 +
                    Gg_coefs.ix['f']['AR%d_Upper' % AR]*(np.log((alpha_g*t)/rb**2))**2 +
                    Gg_coefs.ix['g']['AR%d_Upper' % AR]*(np.log((alpha_g*t)/rb**2)) +
                    Gg_coefs.ix['h']['AR%d_Upper' % AR])

    Gg_Lower = (Gg_coefs.ix['a']['AR%d_Lower' % AR]*(np.log((alpha_g*t)/rb**2))**7 +
                    Gg_coefs.ix['b']['AR%d_Lower' % AR]*(np.log((alpha_g*t)/rb**2))**6 +
                    Gg_coefs.ix['c']['AR%d_Lower' % AR]*(np.log((alpha_g*t)/rb**2))**5 +
                    Gg_coefs.ix['d']['AR%d_Lower' % AR]*(np.log((alpha_g*t)/rb**2))**4 +
                    Gg_coefs.ix['e']['AR%d_Lower' % AR]*(np.log((alpha_g*t)/rb**2))**3 +
                    Gg_coefs.ix['f']['AR%d_Lower' % AR]*(np.log((alpha_g*t)/rb**2))**2 +
                    Gg_coefs.ix['g']['AR%d_Lower' % AR]*(np.log((alpha_g*t)/rb**2)) +
                    Gg_coefs.ix['h']['AR%d_Lower' % AR])
    """

    UU = q*(Rp_cond+Rp_conv) + q*Rc*Gp_Upper + q/(2*np.pi*ks)*Gg_Upper+T0
    UU = pd.DataFrame({'Temp': UU})
    UU['Seconds'] = t
    UU['Hours'] = t/3600

    LL = q*(Rp_cond+Rp_conv) + q*Rc*Gp_Lower + q/(2*np.pi*ks)*Gg_Lower+T0
    LL = pd.DataFrame({'Temp': LL})
    LL['Seconds'] = t
    LL['Hours'] = t/3600
    
    UL = q*(Rp_cond+Rp_conv) + q*Rc*Gp_Upper + q/(2*np.pi*ks)*Gg_Lower+T0
    UL = pd.DataFrame({'Temp': UL})
    UL['Seconds'] = t
    UL['Hours'] = t/3600

    LU = q*(Rp_cond+Rp_conv) + q*Rc*Gp_Lower + q/(2*np.pi*ks)*Gg_Upper+T0
    LU = pd.DataFrame({'Temp': LU})
    LU['Seconds'] = t
    LU['Hours'] = t/3600
    
    Gfunc_result = {'UU': UU, 'LL': LL, 'UL':UL, 'LU':LU}
    
    return(Gfunc_result)
    
#Gfunc_sim(2.7,2400*750,0.01,60,0.3,np.arange(0,100000,1000),18.35,6,13.6,0.25,50) 
    
    
def Gfunc_sim_VR(ks,Cpvs,Rc,q,rb,t,T0,n,SDR,F,AR,pipe_placement='PNE'):
    """
    ks = soil thermal conductivity
    Cpvs = soil vollumetric specific heat capacity
    Rb = pile thermal resistance
    q = applied power (W/m)
    rb = pile radius
    t = time (s)
    T0 = Undisturbed ground temperature
    n = number of pipes in pile
    SDR = pipe standard dimensional ratio
    F = flow rate (L/s)
    AR = aspect ration, must be on of: 15, 25, 33 or 50
    pipe_placement = must be: 'PNE' or 'CP'
    """
    
    ##############################################################################
    #n = 6 #number of pipes
    #SDR = 13.6 
    #F = 0.25 #flow rate (L/s)

    do = 0.025 #pipe outer diameter
    p_thickness = do/SDR #pipe wall thickness (m)
    di = do-2*p_thickness #inner pipe diameter
    ri = di/2
    ro = do/2

    u_ave = F/(np.pi*ri**2) #average fluid velocity (m/s)

    rhof = 1000 #fluid density (kg/m3)
    mu = 0.001 #dynamic viscosity (Pa.s)
    cpf = 4180 #fluid specific heat capacity (J/kg.K)
    kf = 0.582 #fluid thermal conductivity (W/m.K)
    kp = 0.4 #pipe thermal conductivity (W/m.K)
    e = 0.0015 #pipe absolute roughness (mm)

    Re = (rhof*u_ave*di)/mu #Reynolds number
    Pr = (cpf*mu)/kf # Prandtl Number
    hi = (0.023*Re**0.8*Pr**0.35*kf)/(2*ri) #heat transfer coefficient

    # The Convective fluid resistance is..
    Rp_conv=1/(2*n*np.pi*ri*hi)
    #print(Rp_conv)

    # The conductive pipe wall resistance
    Rp_cond=np.log(ro/ri)/(2*n*np.pi*kp)
    #print(Rp_cond)
    
    ##############################################################################    
    # Import Gfunc coefs from file
    #Pile
    Gp_coefs = pd.read_csv(fileOut + 'Gfunc_pile_coefs.csv')
    #Gp_coefs = ',Coefficient,CP_Lower,CP_Upper,PNE_Lower,PNE_Upper\n0,a,-0.000101,3.5499999999999996e-05,-1.4400000000000001e-05,-2.99e-05\n1,b,-0.000234,6.02e-05,1.28e-05,-8.04e-06\n2,c,0.003037,-0.000603,0.000953,0.0008609999999999999\n3,d,0.0018030000000000001,0.0013,0.000131,-0.0011300000000000001\n4,e,-0.04339,-0.00744,-0.0245,-0.0109\n5,f,0.1029,0.0256,0.0757,0.0479\n6,g,0.9095,0.9690000000000001,0.9209999999999999,0.9390000000000001\n'
    #Gp_coefs = pd.read_csv(Gp_coefs)
    Gp_coefs.set_index('Coefficient', inplace=True)
    
    #Ground
    Gg_coefs = pd.read_csv(fileOut + 'Gfunc_ground_coefs2.csv')
    #Gg_coefs = ',Coefficient,AR15_Upper,AR25_Upper,AR33_Upper,AR50_Upper,AR15_Lower,AR25_Lower,AR33_Lower,AR50_Lower\n0,a,-4.837e-07,-3.796e-07,-2.192e-07,-5.1420000000000006e-08,2.68e-07,-6.107999999999999e-07,-8.984e-07,-8.741e-08\n1,b,6.5970000000000005e-06,6.441e-06,4.311e-06,8.756e-07,-1.306e-05,1.83e-05,3.137e-05,8.242999999999999e-06\n2,c,6.592e-05,4.129e-05,2.939e-05,3.233e-05,0.0001827,-0.0001942,-0.0003894,-0.0001835\n3,d,-0.0008843,-0.0008687,-0.0007328,-0.0005292,-9.15e-05,0.001366,0.002361,0.001894\n4,e,-0.004678,-0.003276,-0.002647,-0.00279,-0.01434,-0.01275,-0.01257,-0.01375\n5,f,0.03975,0.04415,0.0443,0.04284,0.05634,0.04932,0.04341,0.04905\n6,g,0.3018,0.3071,0.3076,0.3144,0.3722,0.3863,0.3928,0.3997\n7,h,0.5715,0.5819,0.5861,0.59,0.3989,0.4173,0.4245,0.4267\n'
    #Gg_coefs = pd.read_csv(Gg_coefs)
    Gg_coefs.set_index('Coefficient', inplace=True)
    
    
    ##############################################################################
    # calcualte Alpha
    alpha_g = ks/(Cpvs)


    def Gp(t,alpha_g,rb, bound):
        """
        Function for calculating the pile Gfunc from Loveridge and Powerie 2013.

        t = time in seconds, must be a pd.Series (so that the pd.Series.apply function works)
        alpha_g = thermal diffusivity of the ground - ks/(rho_s*Cp_s)) (m/s)
        rb = pile radius
        bound = either 'Upper' or 'Lower', must be a string
        """

        Fo = (alpha_g*t)/rb**2

        def check_Fo(Fo):
            if Fo < 0.01:
                return(0)
            elif (Fo >= 0.01) & (Fo <= 10.00):
                return(Gp_coefs.ix['a'][pipe_placement + '_' + bound]*(np.log(Fo))**6 +
                            Gp_coefs.ix['b'][pipe_placement + '_' + bound]*(np.log(Fo))**5 +
                            Gp_coefs.ix['c'][pipe_placement + '_' + bound]*(np.log(Fo))**4 +
                            Gp_coefs.ix['d'][pipe_placement + '_' + bound]*(np.log(Fo))**3 +
                            Gp_coefs.ix['e'][pipe_placement + '_' + bound]*(np.log(Fo))**2 +
                            Gp_coefs.ix['f'][pipe_placement + '_' + bound]*(np.log(Fo)) +
                            Gp_coefs.ix['g'][pipe_placement + '_' + bound])
            elif Fo > 10.00:
               return(1)

        return(Fo.apply(lambda x: check_Fo(x)))

    def Gg_Upper(t,alpha_g,rb, AR):
        """
        Function for calculating the upper bound ground Gfunc from Loveridge and Powerie 2013.

        t = time in seconds, must be a pd.Series (so that the pd.Series.apply function works)
        alpha_g = thermal diffusivity of the ground - ks/(rho_s*Cp_s)) (m/s)
        rb = pile radius
        AR = pile aspect ratio, either 15, 25, 33, or 50
        """

        Fo = (alpha_g*t)/rb**2

        def check_Fo(Fo):
            if Fo < 0.1:
                return(0)
            elif Fo >= 0.1:
                return(Gg_coefs.ix['a']['AR%d_Upper' % AR]*(np.log(Fo))**7 +
                        Gg_coefs.ix['b']['AR%d_Upper' % AR]*(np.log(Fo))**6 +
                        Gg_coefs.ix['c']['AR%d_Upper' % AR]*(np.log(Fo))**5 +
                        Gg_coefs.ix['d']['AR%d_Upper' % AR]*(np.log(Fo))**4 +
                        Gg_coefs.ix['e']['AR%d_Upper' % AR]*(np.log(Fo))**3 +
                        Gg_coefs.ix['f']['AR%d_Upper' % AR]*(np.log(Fo))**2 +
                        Gg_coefs.ix['g']['AR%d_Upper' % AR]*(np.log(Fo)) +
                        Gg_coefs.ix['h']['AR%d_Upper' % AR])

        return(Fo.apply(lambda x: check_Fo(x)))
    
    def Gg_Lower(t,alpha_g,rb, AR):
        """
        Function for calculating the lower bound ground Gfunc from Loveridge and Powerie 2013.

        t = time in seconds, must be a pd.Series (so that the pd.Series.apply function works)
        alpha_g = thermal diffusivity of the ground - ks/(rho_s*Cp_s)) (m/s)
        rb = pile radius
        AR = pile aspect ratio, either 15, 25, 33, or 50
        """

        Fo = (alpha_g*t)/rb**2

        def check_Fo(Fo):
            if Fo < 0.25:
                return(0)
            elif Fo >= 0.25:
                return(Gg_coefs.ix['a']['AR%d_Lower' % AR]*(np.log(Fo))**7 +
                        Gg_coefs.ix['b']['AR%d_Lower' % AR]*(np.log(Fo))**6 +
                        Gg_coefs.ix['c']['AR%d_Lower' % AR]*(np.log(Fo))**5 +
                        Gg_coefs.ix['d']['AR%d_Lower' % AR]*(np.log(Fo))**4 +
                        Gg_coefs.ix['e']['AR%d_Lower' % AR]*(np.log(Fo))**3 +
                        Gg_coefs.ix['f']['AR%d_Lower' % AR]*(np.log(Fo))**2 +
                        Gg_coefs.ix['g']['AR%d_Lower' % AR]*(np.log(Fo)) +
                        Gg_coefs.ix['h']['AR%d_Lower' % AR])

        return(Fo.apply(lambda x: check_Fo(x)))



    UU = q*(Rp_cond+Rp_conv) + q*Rc*Gp(t,alpha_g,rb, 'Upper') + q/(2*np.pi*ks)*Gg_Upper(t,alpha_g,rb, AR)+T0
    UU = pd.DataFrame({'Temp': UU})
    UU['Seconds'] = t
    UU['Hours'] = t/3600

    LL = q*(Rp_cond+Rp_conv) + q*Rc*Gp(t,alpha_g,rb, 'Lower') + q/(2*np.pi*ks)*Gg_Lower(t,alpha_g,rb, AR)+T0
    LL = pd.DataFrame({'Temp': LL})
    LL['Seconds'] = t
    LL['Hours'] = t/3600
    
    UL = q*(Rp_cond+Rp_conv) + q*Rc*Gp(t,alpha_g,rb, 'Upper') + q/(2*np.pi*ks)*Gg_Lower(t,alpha_g,rb, AR)+T0
    UL = pd.DataFrame({'Temp': UL})
    UL['Seconds'] = t
    UL['Hours'] = t/3600

    LU = q*(Rp_cond+Rp_conv) + q*Rc*Gp(t,alpha_g,rb, 'Lower') + q/(2*np.pi*ks)*Gg_Upper(t,alpha_g,rb, AR)+T0
    LU = pd.DataFrame({'Temp': LU})
    LU['Seconds'] = t
    LU['Hours'] = t/3600
    
    Gfunc_result = {'UU': UU, 'LL': LL, 'UL':UL, 'LU':LU}
    
    return(Gfunc_result)


#Gfunc_sim_VR(2.7,2400*750,0.01,60,0.3,pd.Series(np.arange(0,100000,1000)),18.35,6,13.6,0.25,50) 


def Gfunc_ksRb(ks,Cpvs,Rc,Power,rb,t,T0,n,SDR,F,AR,pipe_placement,H,temp,t_start,t_end,method, verbose):
    
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
    
    Gfunc_bounds = ['UU', 'LL', 'UL', 'LU']
    
    Gfunc_results = []
    for i in Gfunc_bounds:
        def RMSE(param, Gfunc_bound):
            pred = Gfunc_sim(param[0],Cpvs,param[1],q,rb,t[t1:t2]*3600,T0,n,SDR,F,AR,pipe_placement)[Gfunc_bound].Temp
            E = np.sqrt(np.mean((pred - temp[t1:t2])**2))
            if verbose == True: print('RMSE = %f' % (E))
            return(E)
        #RMSE([2.7,0.2*10])

        # Optimization using minimize function
        param=np.array([ks,Rc]) #array of rough guess of input parameters
        bounds=((0,4), (0,10))
        loc_min=minimize(RMSE, param, args = (i), method=method, bounds=bounds)
        ks_result=loc_min.x[0]
        Rc_result=loc_min.x[1]
        print("\n" + "{}".format(loc_min.message))
        print("Final RMSE = %0.5f\n" % (loc_min.fun))
        Gfunc_results.append({i: [ks_result, Rc_result]})
    

    return(Gfunc_results)




def Gfunc_ksRb_VR(ks,Cpvs,Rc,Power,rb,t,T0,n,SDR,F,AR,pipe_placement,H,temp,t_start,t_end,method, verbose):
    
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
    
    Gfunc_bounds = ['UU', 'LL', 'UL', 'LU']
    
    Gfunc_results = []
    for i in Gfunc_bounds:
        def RMSE(param, Gfunc_bound):
            pred = Gfunc_sim_VR(param[0],Cpvs,param[1],q,rb,t[t1:t2]*3600,T0,n,SDR,F,AR,pipe_placement)[Gfunc_bound].Temp
            E = np.sqrt(np.mean((pred - temp[t1:t2])**2))
            if verbose == True: print('RMSE = %f' % (E))
            return(E)
        #RMSE([2.7,0.2*10])

        # Optimization using minimize function
        param=np.array([ks,Rc]) #array of rough guess of input parameters
        bounds=((0,4), (0,10))
        loc_min=minimize(RMSE, param, args = (i), method=method, bounds=bounds)
        ks_result=loc_min.x[0]
        Rc_result=loc_min.x[1]
        print("\n" + "{}".format(loc_min.message))
        print("Final RMSE = %0.5f\n" % (loc_min.fun))
        Gfunc_results.append({i: [ks_result, Rc_result]})
    

    return(Gfunc_results)






