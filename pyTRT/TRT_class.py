""" Base TRT class """

import pandas as pd
import numpy as np

class TRT:
    """
    The TRT class allows the creation, import, and analysis of TRT data for a range of different models.

    Class initialisation requires setting of several key parameters:
        name:   TRT class object name, eg. TRT1
        H:      GHE length (m)
        Cpvs:   Volumetric Specific heat capacity of the ground aroung GHE (J/m^3K)
        rb:     GHE radius (m)
        T0:     Undisturbed Ground Temperature (degC)
        ks:     Ground Thermal Conductivity (W/mK)
        Rb:     GHE Thermal Resistance (mK/W)
    
    These are the minimum required to simulate a TRT with the Infinite Line Source Model, and are used as the
    default first guess for parameter estimation methods.
    
    The following parameters are calcualted on class init:
        alpha:  Soil thermal diffusivity (m/s^2)
    """

    def __init__(self,name,params):
        """ Initial TRT class with most basic parameter settings, enough to generate a plot of the ILSM.
            ks, Cpvs and Rb parameters will be used as first guess in curve fitting.

            params is a dict containing key value pairs
        """
        self.name=name
        self.params=params

        # Get individual params from params dict
        self.H=self.params['H']
        self.Cps=self.params['Cps']
        self.rhos=self.params['rhos']
        self.rb=self.params['rb'] 
        self.T0=self.params['T0'] 
        self.ks=self.params['ks'] 
        self.Rb=self.params['Rb']

        # Calculate additional params
        self.Cpvs=self.rhos*self.Cps
        self.alpha=self.ks/self.Cpvs

        
    ################################################################################
    # Import TRT Data      
    ################################################################################
    def import_msrd(self,data_file):
        """ Import function for measured data. Essentially just for Stuarts data. Column names need to match TRT19
            Power (W) = Flow (L/s) * (Tin-Tout) * Cpv.w (J/LK)
        """
        self.data=pd.read_csv(data_file)
        self.data['Hours']=abs(self.data.Seconds/3600)
        self.data['logHours']=np.log(self.data.Hours)
        
        # Need to convert flow from L/min to L/s
        self.data.Flow=self.data.Flow/60
        
        # Calculate the power applied to the borehole
        self.data['Power']=abs(self.data.Flow*(self.data.Temp_in - self.data.Temp_out)*4200)
        self.data['DeltaT'] = self.data.Temp_in - self.data.Temp_out
        
    def import_Zprof(self,data_file):        
        # function to read in Zprof data and attach to TRT object
        self.zdata=pd.read_csv(data_loc + 'season_mod/constq_data_from_matlab/' + data_file)
        self.zdata.rename(columns={'depth_Tin':'depth'}, inplace=True)
          
    def import_TRTrig(self,data_file,skipfooter=1, calib = True):
        """ Method to import TRT rig data.
            data_file = file name of dataset.
            skiprows = number of rows to skip (int)
            skipfooter = number of rows at bottom to skip (int)
        """
        self.data=pd.read_csv(data_file, skiprows=2, skipfooter=skipfooter, engine='python', header=None)
        
        # Make sure there are only 23 columns
        if len(self.data.columns)==24:
            self.data.drop(23, 1, inplace=True)
        elif len(self.data.columns)==23:
            self.data=self.data
            
        # change the column names to something useful
        #### Old export file ####
        """
        self.data.columns = ['Date', 'Time', 'Temp_in','Temp_out', 'Control_Tank_Temp',
                             'Return_Tank_Temp', 'Supply_Voltage',
                             'Supply_Current', 'Supply_Watts', 'Supply_VoltAmps',
                             'Compressor_Speed', 'Pump_Speed',
                             'Coil_Off_Air_Temp', 'Coil_On_Air_Temp',
                             'Heat_Mode_Flow_Rate', 'Cool_Mode_Flow_Rate',
                             'Heat_Mode_DeltaT_Setpoint', 'Cool_Mode_DeltaT_Setpoint',
                             'Flowrate', 'DeltaT', 'Power', 'Cpv_w']
        """
        # change the column names to something useful 
        # (adds the acummulative power column on the end)
        #"""
        self.data.columns = ['Date', 'Time', 'Temp_in','Temp_out', 'Control_Tank_Temp',
                             'Return_Tank_Temp', 'Supply_Voltage',
                             'Supply_Current', 'Supply_Watts', 'Supply_VoltAmps',
                             'Compressor_Speed', 'Pump_Speed',
                             'Coil_Off_Air_Temp', 'Coil_On_Air_Temp',
                             'Heat_Mode_Flow_Rate', 'Cool_Mode_Flow_Rate',
                             'Heat_Mode_DeltaT_Setpoint', 'Cool_Mode_DeltaT_Setpoint',
                             'Flow', 'DeltaT', 'Power', 'Cpv_w', 'Accum_Flow']
        #"""
        # Convert string times into pd.Timestamp objects
        self.data['Timestamp']=pd.to_datetime(self.data.Date+' '+self.data.Time, dayfirst=True)
                #NOTE:dayfirst has been set to true since some sets of dates seem to confuse the month and day in the string. 
        
        # Get seconds since start of test
        seconds = [0]
        days = [0]
        for i in range(1,len(self.data.Timestamp)):
            timedelta=self.data.Timestamp[i]-self.data.Timestamp[0]
            seconds.append(timedelta.seconds)
            days.append(timedelta.days)
        self.data['Days']=days
        self.data['Seconds']=seconds + self.data['Days']*86400
        self.data['Hours']=self.data.Seconds/3600 
        
        # save the test duration in hours as a variable
        self.test_end=max(self.data.Hours)
        # add logHours]
        self.data['logHours']=np.log(self.data.Hours)
        
        print("Imported File: " + data_file)        
        
        if calib == True:
            # Copy the the uncalibrated features and make new columns named *_raw
            self.data['Temp_in_raw'] = self.data.Temp_in
            self.data['Temp_out_raw'] = self.data.Temp_out
            self.data['DeltaT_raw'] = self.data.DeltaT
            self.data['Power_raw'] = self.data.Power
            
            self.data['Temp_in'] = (self.data['Temp_in_raw'] + 0.0342)/1.0156
            self.data['Temp_out'] = (self.data['Temp_out_raw'] - 0.2428)/1.0087
            self.data['DeltaT'] = self.data.Temp_in - self.data.Temp_out
            self.data['Power'] = self.data.Flow*self.data.DeltaT*4200
            self.data['Temp_ave']=(self.data.Temp_in +self.data.Temp_out)/2
            
        elif calib==False:
            self.data['Temp_ave']=(self.data.Temp_in +self.data.Temp_out)/2
            print('No calibration applied to this data.\n')

        
    def import_TRTrig_Yale(self,data_file):
        """ Import function for measured data. Essentially just for Stuarts data. Column names need to match TRT19
            Power (W) = Flow (L/s) * (Tin-Tout) * Cpv.w (J/LK)
        """
        self.data=pd.read_csv(data_file)
        self.data['Hours']=abs(self.data.Seconds/3600)
        self.data['logHours']=np.log(self.data.Hours)
        self.data['DeltaT'] = ((self.data.Temp_In1 - self.data.Temp_Out1)+(self.data.Temp_In2 - self.data.Temp_Out2))/2
        
        # Need to convert flow from L/min to L/s
        #self.data.Flow=self.data.Flow/60
        
        # Calculate the power applied to the borehole
        #self.data['Power']=abs(self.data.Flow*(self.data.Temp_in - self.data.Temp_out)*4200)
        
    def import_geocube(self, data_file, skipfooter=0, calib = False):
        """ Method to import Geocube data.
            data_file = file name of dataset. Automatically includes data_loc variable at start of string.
            skipfooter = number of rows at bottom to skip (int)
        """
        self.data=pd.read_csv(data_file, skiprows=1, skipfooter=skipfooter, engine='python', header=0)
            
        # change the column names to something useful
        self.data.columns = ['#', 'Time', 'Voltage', 'Current', 'Counts', 
                             'Temp_in1', 'Temp_in2', 'Temp_out1', 'Temp_out2']
                             
        # Convert temps in degF to degC
        self.data.Temp_in1 = (self.data.Temp_in1 - 32) * (5/9)
        self.data.Temp_in2 = (self.data.Temp_in2 - 32) * (5/9)
        self.data.Temp_out1 = (self.data.Temp_out1 - 32) * (5/9)
        self.data.Temp_out2 = (self.data.Temp_out2 - 32) * (5/9)
                             
        # Convert string times into pd.Timestamp objects
        self.data['Timestamp']=pd.to_datetime(self.data.Time, dayfirst=False)
        
        # Turn timestamp into seconds, hours, days
        seconds = [0]
        days = [0]
        for i in range(1,len(self.data.Timestamp)):
            timedelta=self.data.Timestamp[i]-self.data.Timestamp[0]
            seconds.append(timedelta.seconds)
            days.append(timedelta.days)
        self.data['Days']=days
        self.data['Seconds']=seconds + self.data['Days']*86400
        self.data['Hours']=self.data.Seconds/3600 
        
        # Do some more stuff
        self.test_end=max(self.data.Hours)
        self.data['logHours']=np.log(self.data.Hours)
        self.data['Temp_ave']=(self.data.Temp_in1 + self.data.Temp_in2 + self.data.Temp_out1 + self.data.Temp_out1)/4
        self.data['DeltaT'] = ((self.data.Temp_in1 - self.data.Temp_out1)+(self.data.Temp_in2 - self.data.Temp_out2))/2
                
        
        # Convert Counts to flow rate
        ## First need to get the timestep in seconds
        timesteps = [self.data.Seconds[0]]
        for ix, i in enumerate(self.data.Seconds[:len(self.data.Seconds)-1]):
            timesteps.append(self.data.Seconds[ix+1] - self.data.Seconds[ix])
            
        self.data['timesteps'] = timesteps
            
        ### check to make sure all the timestep entries are the same
        def all_same(items):
            return all(x == items[0] for x in items)
        
        if all_same(timesteps) == True:
            timestep=np.mean(timesteps)
        elif all_same(timesteps) == False:
            print("WARNING: flow rate is not consistent, check and make sure it does not impact on power calculation")
            timestep=np.mean(timesteps)
        
        ## calculate the flow rate in L/s
        self.flow_coef = 3.785 #L/count
        self.data['Flow2'] = (self.data.Counts.cumsum() * self.flow_coef)/self.data.Seconds # L/count for geocube flow meter calculated based on GLD
        self.data['Flow'] = (self.data.Counts * self.flow_coef)/self.data.timesteps # L/count for geocube flow meter calculated based on GLD
        ## Calculate the power applied to the borehole
        self.data['Power2']=abs(self.data.Flow2*self.data.DeltaT*4200)
        self.data['Power']=abs(self.data.Flow*self.data.DeltaT*4200)
        ## Calculate the true power
        self.data['TruePower'] = self.data.Voltage * self.data.Current
         
         
    def import_comsol(self, data_file):
        """ Method to import comsol data from the matlab code.
            data_file = file name of dataset. Automatically includes data_loc variable at start of string.
            
            Includes TRT temp response and the Rg already calculated.
        """    
        import numpy as np
        self.data=pd.read_csv(''.join((self.data_loc,data_file)), skiprows=4)
        
        # Rename Columns        
        self.data.columns = ['Seconds', 'Temp_out', 'Temp_in', 'Power']
        
        # Add additional variables
        self.data['Hours']    = self.data.Seconds/3600
        self.data['Temp_ave'] = (self.data.Temp_out + self.data.Temp_in)/2
        self.data['logHours'] = np.log(self.data.Hours)
        
        # save the test duration in hours as a variable
        self.test_end = max(self.data.Hours)
          
    def import_comsol_new(self, data_file):
        """ Method to import comsol data from the matlab code.
            data_file = file name of dataset. Automatically includes data_loc variable at start of string.
            
            Includes TRT temp response and the Rg already calculated.
        """    
        import numpy as np
        
        self.data=pd.read_csv(''.join((self.data_loc,'season_mod/constq_data_from_matlab/',data_file)))
        #self.data=pd.read_csv(''.join((self.data_loc,data_file)))
        
        # Add additional variables
        self.data['Seconds']  = self.data.time
        self.data['Hours']    = self.data.Seconds/3600
        self.data['Temp_ave'] = (self.data.Temp_out + self.data.Temp_in)/2
        self.data.Power       = -self.data.Power
        self.data['logHours'] = np.log(self.data.Hours)
        self.data['Temp_bhw_ave'] = self.data.Temp_bhw/(2*np.pi*self.rb*self.H) - 273.15 #get the average BHW wall temp in degC 
        
        # save the test duration in hours as a variable
        self.test_end = max(self.data.Hours)
            
    ################################################################################       
    def import_comsol_old(self,data_file, Rb='no'):
        """ Method for importing the old COMOSL TRT data.
            data_file = file name of dataset. Automatically includes data_loc variable at start of string.
            
            For each TRT data file, there is a corresponding 3 Rb data files stored in the Rb subdirectory.
            Rb - If Rb = yes, the function will look in the Rb subdirectory, import the Rb data, and caluclate Rb, 
                         then add each to the TRT.data attribute.
                 If Rb = no (deafult), the import function will just import the TRT temp data.
        """
        self.data=pd.read_csv(''.join((self.data_loc,'season_mod/constq_data/',data_file)))
        
        # Rename columns
        self.data.columns = ['Seconds', 'Temp_out', 'Temp_in', 'Temp_ave', 'Power']
        
        # Add additional Variables
        self.data['Hours'] = self.data.Seconds/3600
        self.data.Power = -self.data.Power
        
        # save the test duration in hours as a variable
        self.test_end = max(self.data.Hours)
        
        # Now make a df of just the temp_ave quick access
        df = pd.DataFrame()
        df['Hours'] = self.data.Hours
        df['Temp'] = self.data.Temp_ave
        self.msrd_data['data'] = df
        
        ### Thermal Resistance ###
        # Define Rb file list
        Rb_file='Rb' + data_file[3:-4] + '*.csv' 
        
        # Import Rb data
        if Rb == 'yes':
            # Get list of files matching TRT data with glob
            Rb_files = glob.glob(self.data_loc + 'season_mod/constq_data/Rb/' + Rb_file)
            
            # Import the 3 Rb variables and add to the TRT.data attribute
            self.data['Temp_bhw_ext'] = pd.read_csv(Rb_files[0], skiprows = 4).ix[:,1]
            self.data['Temp_pipe_ext'] = pd.read_csv(Rb_files[1], skiprows = 4).ix[:,1]
            self.data['Power_pipe_wall'] = pd.read_csv(Rb_files[2], skiprows = 4).ix[:,1]
            
            # Calculate Rg and add to TRT.data
            # Set constants
            self.Hb=30.066 #length of borehole (m)
            self.spacing = 0.04 #spacing between pipes (m)
            self.H_total = 2*self.H + self.spacing
            self.Abhw = 2*np.pi*self.rb*self.Hb #Surface are of sides of the borehole wall (m2)
            
            self.spacing = 0.04 #spacing between pipes (m)
            self.bottom_gap = 0.36602950000000006
            self.Hb = self.H + self.bottom_gap #length of borehole (m)
            self.H_total = 2*self.H + self.spacing
            self.Abhw = 2*np.pi*0.125/2*self.Hb #Surface are of sides of the borehole wall (m2)
            
            # Calculate the grout resistance.
            # NOTE: I assume that the Rg from COMSOL = Rb from analytical solutions because there isnt any pipe materials in the 3D models.
            #       This may be wrong however, need to investigate further.
            self.data['Rg']=((self.data.Temp_pipe_ext/self.H_total) - (self.data.Temp_bhw_ext/self.Abhw))/(-self.data.Power_pipe_wall/self.H_total)
            
        elif Rb == 'no':
            'do nothing'
        
        else:
            print('ERROR: Rb must be <yes> or <no>')

    ################################################################################
    # Simulate TRT data
    ################################################################################
        
    def ILSM_sim(self,q,t):
        """ Basic ILSM simulation function: ILSM(q,times)
            Takes times as any iterable in hours, outputs temperature based on the parameters set at TRT class initialization
            q:      power (W/m)
            t = time/s at which temperature occurs (hrs). Can be int, float, list, or an 1D iterable (np.array, pd.Series)
        """
        # Import the ILSM function
        from pyTRT.Models.ILSM import ILSM_sim
        
        temps=ILSM_sim(self.ks,self.Rb,self.alpha,q,self.rb,t,self.T0)
        
        if isinstance(t, float)==True or isinstance(t, int)==True:
            return(temps)
            
        else:
            df=pd.DataFrame()
            df['Hours']=t
            df['Temp']=temps
            self.ILSM_sim_result=df
            return(df)
            
    def ILSM_EI_sim(self,q,t):
        """ Basic ILSM simulation function: ILSM(q,times)
            Takes times as any iterable in hours, outputs temperature based on the parameters set at TRT class initialization
            q:      power (W/m)
            t = time/s at which temperature occurs (hrs). Can be int, float, list, or an 1D iterable (np.array, pd.Series)
        """
        # Import the ILSM function
        from pyTRT.Models.ILSM import ILSM_EI_sim
        
        temps=ILSM_EI_sim(self.ks,self.Rb,self.alpha,q,self.rb,t,self.T0)
        
        if isinstance(t, float)==True or isinstance(t, int)==True:
            return(temps)
            
        else:
            df=pd.DataFrame()
            df['Hours']=t
            df['Temp']=temps
            self.ILSM_EI_sim_result=df
            return(df)
        
    def FLSM_sim(self,z,q,t):
        """ Basic FLSM simulation function: FLSM_sim(ks,a,H,rb,z,T0,q,Rb,t)
        
            Takes time inputs and outputs temperature at some depth along the GHE z.
            
            z = depth at which to evaluate the fluid temperature (m) NOTE: H/2 ~ mean fluid temperature
            q = power (W/m)
            times = int, float, series, or np.array (hours)
        """
        from pyTRT.Models.FLSM import FLSM_sim
        
        temps=FLSM_sim(self.ks,self.Rb,self.alpha,self.H,self.rb,z,self.T0,q,t)
        
        if isinstance(t, float)==True or isinstance(t, int)==True:
            return(temps)
            
        else:
            df=pd.DataFrame()
            df['Hours']=t
            df['Temp']=temps
            self.FLSM_sim_result=df
            return(df)
            
    def FLSM_sim_intmean(self,q,t):
        """ Basic FLSM simulation function: FLSM_sim(ks,a,H,rb,z,T0,q,Rb,t)
        
            Takes time inputs and outputs temperature at some depth along the GHE z.
            
            z = depth at which to evaluate the fluid temperature (m) NOTE: H/2 ~ mean fluid temperature
            q = power (W/m)
            times = int, float, series, or np.array (hours)
        """
        from pyTRT.Models.FLSM import FLSM_sim_intmean
        
        temps=FLSM_sim_intmean(self.ks,self.Rb,self.alpha,self.H,self.rb,self.T0,q,t)
        
        if isinstance(t, float)==True or isinstance(t, int)==True:
            return(temps)
            
        else:
            df=pd.DataFrame()
            df['Hours']=t
            df['Temp']=temps
            self.FLSM_sim_result=df
            return(df)
            
    def FLSM_sim_Lamarche(self,q,t):
        """ Basic FLSM simulation function: FLSM_sim(ks,a,H,rb,z,T0,q,Rb,t)
        
            Takes time inputs and outputs temperature at some depth along the GHE z.
            
            q = power (W/m)
            times = int, float, series, or np.array (hours)
        """
        from pyTRT.Models.FLSM import FLSM_sim_Lamarche
        
        temps=FLSM_sim_Lamarche(self.ks,self.Rb,self.alpha,self.H,self.rb,self.T0,q,t)
        
        if isinstance(t, float)==True or isinstance(t, int)==True:
            return(temps)
            
        else:
            df=pd.DataFrame()
            df['Hours']=t
            df['Temp']=temps
            self.FLSM_sim_result=df
            return(df)
            
            
    def FLSM_sim_ave(self,z,q,t):
        """ Basic FLSM simulation function: FLSM_sim(ks,a,H,rb,z,T0,q,Rb,t)
        
            Takes time inputs and outputs temperature at some depth along the GHE z.
            
            z = depth at which to evaluate the fluid temperature (m) NOTE: H/2 ~ mean fluid temperature
            q = power (W/m)
            times = int, float, series, or np.array (hours)
        """
        from pyTRT.Models.FLSM import FLSM_sim
        
        temps=FLSM_sim(self.ks,self.Rb,self.alpha,self.H,self.rb,z,self.T0,q,t)
        
        if isinstance(t, float)==True or isinstance(t, int)==True:
            return(temps)
            
        else:
            df=pd.DataFrame()
            df['Hours']=t
            df['Temp']=temps
            self.FLSM_sim_result=df
            return(df)
            
    def ICSM_sim(self,q,t):
        """ Basic ICSM simulation function: ICSM_sim(ks,Rb,Cpvs,rb,q,T0,t)
            Takes time inputs and outputs mean fluid temperature.
            
            Parameters:
            q = power (W/m)
            t = int, float, or iterable in hours
        """
        from pyTRT.Models.ICSM import ICSM_sim
        
        temps=ICSM_sim(self.ks,self.Rb,self.Cpvs,self.rb,q,self.T0,t)
        
        if isinstance(t, float)==True or isinstance(t, int)==True:
            return(temps)
            
        else:
            df=pd.DataFrame()
            df['Hours']=t
            df['Temp']=temps
            self.ICSM_sim_result=df
            return(df)
        
    ################################################################################
    # ILSM Linereg Functions
    ################################################################################
    
    def ILSM_linreg(self,t_start,t_end,mode='heating',plot=False):
        """ 
        Fit a linear model with linear regression to the measured TRT data and
        use the ILSM to calculate thermal conductivity and borehole resistance.
        
        The linear model is fitted between t_start and t_end, which define the fitting region, or data window.
        Generally, t_end should be taken as the last data point in the dataset.
        
        Parameters:

        t_start:    starting time of regression window (hours)
        t_end:      end time of regression window (hours)
        mode:       either 'heating' (default) or 'cooling'
        plot:       True or False (default), should the function call plot the linear model against the data?
        
        """
        
        from pyTRT.Models.ILSM import ILSM_linreg
        import matplotlib.pyplot as plt
        
        self.results=ILSM_linreg(self.data.Hours,self.data.Temp_ave,
                                self.data.Power,t_start,t_end,
                                self.H,self.T0,self.rb,self.Cpvs, mode)
        
        if mode == 'heating':
            mode_sign = 1
        elif mode == 'cooling':
            mode_sign = -1
        else:
            print('ERROR: select mode as heating or cooling')
        
        # Creat results output dictionary    
        self.ILSM_linreg_result = {}                            
        self.ILSM_linreg_result['ks']=self.results['ks']
        self.ILSM_linreg_result['Rb']=self.results['Rb']
        self.ILSM_linreg_result['m']=self.results['m']
        self.ILSM_linreg_result['c']=self.results['c']
        
        # OPTIONAL: print the ILSM fit over the meausred data
        if plot==True:
            f, ax1 = plt.subplots(1,1, figsize=(8,5))
            ax1.plot(self.data.logHours, self.data.Temp_ave, color='y', label=self.name)
            ax1.set_xlabel('log(Hours)')
            ax1.set_ylabel('Temperature (°C)')
            ax1.legend(loc=4)
            ax1.set_title(self.name)
            ax1.grid(b=False)
            ax1.text(0.1,0.9,
                             'ks = %0.2f\nRb = %0.2f' % (self.results['ks'], self.results['Rb']), 
                             fontsize=12, weight= 'normal', horizontalalignment='center', color = 'k',
                             verticalalignment='center',transform = ax1.transAxes, 
                             bbox={'facecolor':'white','alpha':0, 'pad':1})
            
            # plot the linear fit between two chosen x-values
            x=np.log(np.arange(t_start,t_end))
            y=mode_sign*self.results['m']*x+ mode_sign*self.results['c']
            ax1.plot(x, y, color='k', linewidth=1)
            plt.shot()
            
            return(self.ILSM_linreg_result)
        
        elif plot==False:
            return(self.ILSM_linreg_result)
            

            
    def ILSM_linreg_multi(self,t_start1,t_start2,t_end,res):
        """ 
        Take a set of measured data and fit several linear models to it,
        using the ILSM to calculate thermal conductivity and borehole resistance.
        t = time in hours
        Temp = TRT temp response curve
        t_start1 = starting time of first regression window (hours)
        t_start2 = starting time of last regression window (hours)
        res = how many indices to increase the regression window after each run by
        """
        
        from pyTRT.Models.ILSM import ILSM_linreg_multi
        
        self.ILSM_linreg_multi_results=ILSM_linreg_multi(self.data.Hours,self.data.Temp_ave,
                                                            self.data.Power,t_start1,t_start2,t_end,res,
                                                            self.H,self.T0,self.rb,self.Cpvs)
        return(self.ILSM_linreg_multi_results)
        
    def ILSM_linreg_multi_reverse(self,tstart,res):
        
        from pyTRT.Models.ILSM import ILSM_linreg_multi_reverse
        
        self.ILSM_linreg_multi_reverse_results=ILSM_linreg_multi_reverse(self.data.Hours,self.data.Temp_ave,
                                                                         self.data.Power,tstart, res,
                                                                         self.H,self.T0,self.rb,self.Cpvs)
        return(self.ILSM_linreg_multi_reverse_results)
        
    ################################################################################
    # ILSM Fitting Functions
    ################################################################################
    def ILSM_ksRb(self, tstart, tend, method, verbose=False, plot=False):
        """ Fit the ILS to measured data by curve fitting two parameters:
                ks - Soil Thermal Conductivity (W/m.K)
                Rb - Borehole Resistance (m.K/W)
                These are defined at TRT class init as the initial guess.
            
            ILSM formula is from Hemmingway et al 2012. Could be updated to match a more reputable paper.
                
            Also need to define:
                tstart:     start time (hours) of fitting region
                tend:       end time (hours) of fitting region
                method:     (string) the optimzation algorithm to use, 
                            see scikit.minimize function, usually 'L-BFGS-B'
                verbose:    Print the RMSE of each iteration of the optimization algorithm
                plot:       print the fitted function against the measured data?
        """
            
    
        from pyTRT.Models.ILSM import ILSM_ksRb, ILSM_sim
        import matplotlib.pyplot as plt
        
        def find_nearest_idx(array,value):
            idx = (np.abs(array-value)).argmin()
            return(idx)    
        
        t = self.data.Hours
        temp = self.data.Temp_ave
        q = self.data.Power
        
        #if multiStart == False:
        self.ILSM_ksRb_results = ILSM_ksRb(self.ks, self.Rb, self.alpha, 
                                               t, temp, q, self.H, self.T0, self.rb, 
                                               tstart, tend, method, verbose)
        
        if plot==True:
            f, ax1 = plt.subplots(1,1, figsize=(8,5))
            ax1.plot(self.data.Hours, self.data.Temp_ave, color='y', label=self.name)
            ax1.set_xlabel('Hours')
            ax1.set_ylabel('Temperature (°C)')
            ax1.grid(b=False)
            
            # plot the fitted model
            t1 = find_nearest_idx(t,tstart)
            t2 = find_nearest_idx(t,tend)
            pred=ILSM_sim(self.ILSM_ksRb_results['ks'],self.ILSM_ksRb_results['Rb'], self.alpha, 
                              np.mean(q[t1:t2])/self.H, self.rb,t[t1:t2], self.T0)  
            ax1.plot(t[t1:t2], pred, color='k', linewidth=1, label='ILSM_2par Fit')
            plt.legend(loc='best')
            
            return(self.ILSM_ksRb_results)
        
        elif plot==False:
            return(self.ILSM_ksRb_results)
            
            
    def ILSM_KsRb_multi(self, tstarts, tend, method, verbose=False, plot=False):
        from pyTRT.Models.ILSM import ILSM_ksRb, ILSM_sim
        
        def find_nearest_idx(array,value):
            idx = (np.abs(array-value)).argmin()
            return(idx)    
        
        t = self.data.Hours
        temp = self.data.Temp_ave
        q = self.data.Power
        
        ks=[]
        Rb=[]
        for i in tstarts:
            result = ILSM_ksRb(self.ks, self.Rb, self.alpha, 
                               t, temp, q, self.H, self.T0, self.rb, 
                               i, tend, method, verbose)
            ks.append(result['ks'])
            Rb.append(result['Rb'])
        
        result = pd.DataFrame({'Hours':tstarts})
        result['ks'] = ks
        result['Rb'] = Rb
       
        self.ILSM_ksRb_multi_results = result
        return(self.ILSM_ksRb_multi_results)

    def ILSM_EI_ksRb(self, tstart, tend, method, verbose=False, plot=False):
        """ Fit the ILS to measured data by curve fitting two parameters:
                ks - Soil Thermal Conductivity (W/m.K)
                Rb - Borehole Resistance (m.K/W)
                These are defined at TRT class init as the initial guess.
            
            ILSM formula is from Hemmingway et al 2012. Could be updated to match a more reputable paper.
                
            Also need to define:
                tstart:     start time (hours) of fitting region
                tend:       end time (hours) of fitting region
                method:     (string) the optimzation algorithm to use, 
                            see scikit.minimize function, usually 'L-BFGS-B'
                verbose:    Print the RMSE of each iteration of the optimization algorithm
                plot:       print the fitted function against the measured data?
        """
            
    
        from pyTRT.Models.ILSM import ILSM_EI_ksRb, ILSM_EI_sim
        import matplotlib.pyplot as plt
        
        def find_nearest_idx(array,value):
            idx = (np.abs(array-value)).argmin()
            return(idx)    
        
        t = self.data.Hours
        temp = self.data.Temp_ave
        q = self.data.Power
        
        #if multiStart == False:
        self.ILSM_ksRb_results = ILSM_EI_ksRb(self.ks, self.Rb, self.alpha, 
                                               t, temp, q, self.H, self.T0, self.rb, 
                                               tstart, tend, method, verbose)
        
        if plot==True:
            f, ax1 = plt.subplots(1,1, figsize=(8,5))
            ax1.plot(self.data.Hours, self.data.Temp_ave, color='y', label=self.name)
            ax1.set_xlabel('Hours')
            ax1.set_ylabel('Temperature (°C)')
            ax1.grid(b=False)
            
            # plot the fitted model
            t1 = find_nearest_idx(t,tstart)
            t2 = find_nearest_idx(t,tend)
            pred=ILSM_EI_sim(self.ILSM_ksRb_results['ks'],self.ILSM_ksRb_results['Rb'], self.alpha, 
                              np.mean(q[t1:t2])/self.H, self.rb,t[t1:t2], self.T0)  
            ax1.plot(t[t1:t2], pred, color='k', linewidth=1, label='ILSM_EI_2par Fit')
            plt.legend(loc='best')
            
            return(self.ILSM_ksRb_results)
        
        elif plot==False:
            return(self.ILSM_ksRb_results) 
            
    def ILSM_EI_KsRb_multi(self, tstarts, tend, method, verbose=False, plot=False):
        from pyTRT.Models.ILSM import ILSM_EI_ksRb, ILSM_EI_sim
        
        def find_nearest_idx(array,value):
            idx = (np.abs(array-value)).argmin()
            return(idx)    
        
        t = self.data.Hours
        temp = self.data.Temp_ave
        q = self.data.Power
        
        ks=[]
        Rb=[]
        for i in tstarts:
            result = ILSM_EI_ksRb(self.ks, self.Rb, self.alpha, 
                               t, temp, q, self.H, self.T0, self.rb, 
                               i, tend, method, verbose)
            ks.append(result['ks'])
            Rb.append(result['Rb'])
        
        result = pd.DataFrame({'Hours':tstarts})
        result['ks'] = ks
        result['Rb'] = Rb
       
        self.ILSM_EI_ksRb_multi_results = result
        return(self.ILSM_EI_ksRb_multi_results)
            
    def ILSM_ksRbAlpha(self, tstart, tend, method, verbose=False, plot=False):
        
        from pyTRT.Models.ILSM import ILSM_ksRbAlpha, ILSM_sim
        import matplotlib.pyplot as plt
        
        def find_nearest_idx(array,value):
            idx = (np.abs(array-value)).argmin()
            return(idx)    
        
        t = self.data.Hours
        temp = self.data.Temp_ave
        q = self.data.Power
        
        self.ILSM_ksRbAlpha_results = ILSM_ksRbAlpha(self.ks, self.Rb, self.alpha, 
                                               t, temp, q, self.H, self.T0, self.rb, 
                                               tstart, tend, method, verbose)
        
        if plot==True:
            f, ax1 = plt.subplots(1,1, figsize=(8,5))
            ax1.plot(self.data.Hours, self.data.Temp_ave, color='y', label=self.name)
            ax1.set_xlabel('Hours')
            ax1.set_ylabel('Temperature (°C)')
            ax1.grid(b=False)
            
            # plot the fitted model
            t1 = find_nearest_idx(t,tstart)
            t2 = find_nearest_idx(t,tend)
            pred=ILSM_sim(self.ILSM_ksRbAlpha_results['ks'],self.ILSM_ksRbAlpha_results['Rb'], self.ILSM_ksRbAlpha_results['alpha'], 
                              np.mean(q[t1:t2])/self.H, self.rb,t[t1:t2], self.T0)  
            ax1.plot(t[t1:t2], pred, color='k', linewidth=1, label='ILSM_3par Fit')
            plt.legend(loc='best')
            
            return(self.ILSM_ksRbAlpha_results)
        
        elif plot==False:
            return(self.ILSM_ksRbAlpha_results)
    
   
        

    """
    def ILSM_RbAlpha(self, tstart, tend, method, ks='ILSM_LR', verbose=False, plot=False):
        import fit.ILS_opt as ILSM
        import ILS.ILS_base as ILS_base
        
        def find_nearest_idx(array,value):
            idx = (np.abs(array-value)).argmin()
            return(idx)    
        
        t = self.data.Hours
        temp = self.data.Temp_ave
        q = self.data.Power
        
        if ks =='ILSM_LR':
            ks = self.ILS_fit_result['ks']
        else:
            ks=ks
        
        self.ILSM_RbAlpha_results = ILSM.ILS_RbAlpha(ks, self.Rb, self.alpha, 
                                               t, temp, q, self.H, self.T0, self.rb, 
                                               tstart, tend, method, verbose)
        
        if plot==True:
            f, ax1 = plt.subplots(1,1, figsize=(8,5))
            ax1.plot(self.data.Hours, self.data.Temp_ave, color='y', label=self.name)
            ax1.set_xlabel('Hours')
            ax1.set_ylabel('Temperature (°C)')
            ax1.grid(b=False)
            
            # plot the fitted model
            t1 = find_nearest_idx(t,tstart)
            t2 = find_nearest_idx(t,tend)
            pred=ILS_base.ILS(self.ILS_fit_result['ks'],self.ILSM_RbAlpha_results['Rb'], self.ILSM_RbAlpha_results['alpha'], 
                              np.mean(q[t1:t2])/H, self.rb,t[t1:t2], self.T0)  
            ax1.plot(t[t1:t2], pred, color='k', linewidth=1, label='ILSM_3par Fit')
            plt.legend(loc='best')
            
            return(self.ILSM_RbAlpha_results)
        
        elif plot==False:
            return(self.ILSM_RbAlpha_results)    
        """
    
    ################################################################################
    # FLSM Fitting Functions
    ################################################################################
    def FLSM_ksRb(self, tstart, tend, method, verbose=False, plot=False, fluid_calc='mid'):
        
        from pyTRT.Models.FLSM import FLSM_ksRb, FLSM_sim, FLSM_sim_intmean
        import matplotlib.pyplot as plt
        
        def find_nearest_idx(array,value):
            idx = (np.abs(array-value)).argmin()
            return(idx)    
        
        t = self.data.Hours
        temp = self.data.Temp_ave
        Q = self.data.Power
        
        self.FLSM_ksRb_results = FLSM_ksRb(t, temp, Q, tstart, tend, method, 
                                           self.ks, self.Rb, self.alpha, 
                                           self.H, self.rb, self.H/2, self.T0, 
                                           verbose, fluid_calc = fluid_calc)
        
        if plot==True:
            f, ax1 = plt.subplots(1,1, figsize=(8,5))
            ax1.plot(self.data.Hours, self.data.Temp_ave, color='y', label=self.name)
            ax1.set_xlabel('Hours')
            ax1.set_ylabel('Temperature (°C)')
            ax1.grid(b=False)
            
            # plot the fitted model
            t1 = find_nearest_idx(t,tstart)
            t2 = find_nearest_idx(t,tend)
            
            if fluid_calc == 'mid':
                pred = FLSM_sim(self.FLSM_ksRb_results['ks'],self.FLSM_ksRb_results['Rb'],self.alpha,
                               self.H,self.rb,self.H/2,self.T0,np.mean(Q[t1:t2])/self.H,t[t1:t2])
            elif fluid_calc == 'intmean':
                pred = FLSM_sim_intmean(self.FLSM_ksRb_results['ks'],self.FLSM_ksRb_results['Rb'],self.alpha,
                               self.H,self.rb,self.H/2,self.T0,np.mean(Q[t1:t2])/self.H,t[t1:t2])
            
            ax1.plot(t[t1:t2], pred, color='k', linewidth=1, label='FLSM_ksRb')
            plt.legend(loc='best')
            
            return(self.FLSM_ksRb_results)
        
        elif plot==False:
            return(self.FLSM_ksRb_results)
            
            
    def FLSM_ksRb_multi(self, tstart, tend, method, verbose=False, plot=False):
        
        from pyTRT.Models.FLSM import FLSM_ksRb, FLSM_sim
        import matplotlib.pyplot as plt
        
        def find_nearest_idx(array,value):
            idx = (np.abs(array-value)).argmin()
            return(idx)    
        
        t = self.data.Hours
        temp = self.data.Temp_ave
        Q = self.data.Power
        
        self.FLSM_ksRb_results = FLSM_ksRb(t, temp, Q, tstart, tend, method, 
                                           self.ks, self.Rb, self.alpha, 
                                           self.H, self.rb, self.H/2, self.T0, 
                                           verbose)
        
        if plot==True:
            f, ax1 = plt.subplots(1,1, figsize=(8,5))
            ax1.plot(self.data.Hours, self.data.Temp_ave, color='y', label=self.name)
            ax1.set_xlabel('Hours')
            ax1.set_ylabel('Temperature (°C)')
            ax1.grid(b=False)
            
            # plot the fitted model
            t1 = find_nearest_idx(t,tstart)
            t2 = find_nearest_idx(t,tend)
            pred = FLSM_sim(self.FLSM_ksRb_results['ks'],self.FLSM_ksRb_results['Rb'],self.alpha,
                           self.H,self.rb,self.H/2,self.T0,np.mean(Q[t1:t2])/self.H,t[t1:t2])
            
            ax1.plot(t[t1:t2], pred, color='k', linewidth=1, label='FLSM_ksRb')
            plt.legend(loc='best')
            
            return(self.FLSM_ksRb_results)
        
        elif plot==False:
            return(self.FLSM_ksRb_results) 
        
    
    def FLSM_ksRbAlpha(self, tstart, tend, method, verbose=False, plot=False):
        
        from pyTRT.Models.FLSM import FLSM_ksRbAlpha, FLSM_sim
        import matplotlib.pyplot as plt
        
        def find_nearest_idx(array,value):
            idx = (np.abs(array-value)).argmin()
            return(idx)    
            
        t = self.data.Hours
        temp = self.data.Temp_ave
        Q = self.data.Power
        
        self.FLSM_ksRbAlpha_results = FLSM_ksRbAlpha(t, temp, Q, tstart, tend, method, 
                                           self.ks, self.Rb, self.alpha, 
                                           self.H, self.rb, self.H/2, self.T0, 
                                           verbose)

        if plot==True:
            f, ax1 = plt.subplots(1,1, figsize=(8,5))
            ax1.plot(self.data.Hours, self.data.Temp_ave, color='y', label=self.name)
            ax1.set_xlabel('Hours')
            ax1.set_ylabel('Temperature (°C)')
            ax1.grid(b=False)
            
            # plot the fitted model
            t1 = find_nearest_idx(t,tstart)
            t2 = find_nearest_idx(t,tend)
            pred = FLSM_sim(self.FLSM_ksRbAlpha_results['ks'], self.FLSM_ksRbAlpha_results['Rb'], self.FLSM_ksRbAlpha_results['alpha'],
                           self.H, self.rb, self.H/2, self.T0, np.mean(Q[t1:t2])/self.H, t[t1:t2])
            
            ax1.plot(t[t1:t2], pred, color='k', linewidth=1, label='FLSM_ksRbAlpha')
            plt.legend(loc='best')
            
            return(self.FLSM_ksRbAlpha_results)
        
        elif plot==False:
            return(self.FLSM_ksRbAlpha_results) 
            
    ################################################################################
    # ICSM Fitting Functions
    ################################################################################
    def ICSM_ksRb(self, tstart, tend, method, verbose=False, plot=False):
        """ Fit the ILS to measured data by curve fitting two parameters:
                ks - Soil Thermal Conductivity (W/m.K)
                Rb - Borehole Resistance (m.K/W)
                These are defined at TRT class init as the initial guess.
            
            ILSM formula is from Hemmingway et al 2012. Could be updated to match a more reputable paper.
                
            Also need to define:
                tstart:     start time (hours) of fitting region
                tend:       end time (hours) of fitting region
                method:     (string) the optimzation algorithm to use, 
                            see scikit.minimize function, usually 'L-BFGS-B'
                verbose:    Print the RMSE of each iteration of the optimization algorithm
                plot:       print the fitted function against the measured data?
        """
            
    
        from pyTRT.Models.ICSM import ICSM_ksRb, ICSM_sim
        import matplotlib.pyplot as plt
        
        def find_nearest_idx(array,value):
            idx = (np.abs(array-value)).argmin()
            return(idx)    
        
        t = self.data.Hours
        temp = self.data.Temp_ave
        q = self.data.Power
        
        #if multiStart == False:
        self.ICSM_ksRb_results = ICSM_ksRb(self.ks, self.Rb, self.Cpvs, 
                                               t, temp, q, self.H, self.T0, self.rb, 
                                               tstart, tend, method, verbose)
        
        if plot==True:
            f, ax1 = plt.subplots(1,1, figsize=(8,5))
            ax1.plot(self.data.Hours, self.data.Temp_ave, color='y', label=self.name)
            ax1.set_xlabel('Hours')
            ax1.set_ylabel('Temperature (°C)')
            ax1.grid(b=False)
            
            # plot the fitted model
            t1 = find_nearest_idx(t,tstart)
            t2 = find_nearest_idx(t,tend)
            pred=ICSM_sim(self.ICSM_ksRb_results['ks'],self.ICSM_ksRb_results['Rb'], 
                          self.Cpvs,self.rb, np.mean(q[t1:t2])/self.H, self.T0, t[t1:t2])  
            ax1.plot(t[t1:t2], pred, color='k', linewidth=1, label='ICSM_ksRb')
            plt.legend(loc='best')
            
            return(self.ICSM_ksRb_results)
        
        elif plot==False:
            return(self.ICSM_ksRb_results)      
    
    """
    ################################################################################
    # ILSM Transient analysis functions
    ################################################################################
    def ILS_linreg_trans(self, tstarts, tskip, tend, res=0.5, plot=False):
        ""
        Transient analysis with the ILSM and linreg.
        Provide list of fixed regression window start times as tstarts (hours).
        Returns TRT results calcualted over data windows increasing from tstarts 
        to the end of the data by res (hours).
        
        To plot the results in the function call, set plot = 'yes'.
        If not, set to plot = 'no' (default).
        ""
        import matplotlib.pyplot as plt
        import fit.ILS_linreg as ILS_linreg
        import itertools
        self.ILS_linreg_results_trans=ILS_linreg.ILS_linreg_trans(self.data.Hours,self.data.Temp_ave,
                                                            self.data.Power,tstarts,tskip,tend,res,
                                                            self.H,self.T0,self.rb,self.Cpvs)
        if plot == True:
            f, (ax1, ax2)= plt.subplots(1,2, figsize=(6,5), dpi=300)
            palette = itertools.cycle(sns.color_palette('gnuplot'))
            for key, grp in self.ILS_linreg_results_trans.groupby('start_time', axis=0):
                grp.plot('stop_time','ks', ax=ax1, label='Start Time = %d ' % (key), 
                       color=next(palette))
            ax1.set_xlim([min(self.data.Hours), max(self.data.Hours)])    
            ax1.set_ylabel('Thermal Conductivity (W/m.K)')
            ax1.set_title('Thermal Conductivity')
            
            for key, grp in self.ILS_linreg_results_trans.groupby('start_time', axis=0):
                grp.plot('stop_time','Rb', ax=ax2, label='Start Time = %d ' % (key), 
                       color=next(palette))
            ax2.set_xlim([min(self.data.Hours), max(self.data.Hours)])    
            ax2.set_ylabel('GHE Thermal Resistance (m.K/W)')
            ax2.set_title('GHE Thermal Resistance')
            plt.tight_layout()    
            
        elif plot == False:
            return(self.ILS_linreg_results_trans)
        else:
            return('Set plot to True or False')
            
            
    def ILS_trans_hist(self, minRsqrd, ks_hist_range, Rb_hist_range, Rsqrd_hist_range, plot=True):
        import fit.ILS_linreg as ILS_linreg
        
        if self.ILS_linreg_results_trans.empty == True:
            return('ERROR: need to run TRT_obj.ILS_linreg_trans first')
        else:
            self.ILS_trans_hist_results = ILS_linreg.ILS_trans_hist(self.ILS_linreg_results_trans, minRsqrd, ks_hist_range, Rb_hist_range, Rsqrd_hist_range, self.name, plot=plot)
            return(self.ILS_trans_hist_results)
    """            

    
        
        

