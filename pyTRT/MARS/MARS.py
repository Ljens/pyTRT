
def MARS_ILSM(TRTobj, startIndex=0):
    """ 
    Function that calls the MARS R script and applies it to TRT data.
    The cuts are then used to fit the ILSM model in segments to the data, with the thermal properties derived.
    
    Imputs:
        TRTobj = TRT class object with imported data
        startIndex = choose which row in the data to start fitting. Rows before are left out.
    """
    import subprocess
    import os
    
    ### ONLY WORKS ON MAC ###
    # Save the TRT data to an in memory text string in .csv format
    #data_file_InMemory = TRTobj.data[['Hours','Temp_ave']][startIndex:].to_csv()
    
    # Save to actual CSV then import
    # This saves it to the current working directory
    TRTobj.data[['Hours','Temp_ave']][startIndex:].to_csv('MARS_script_data_file.csv')
    data_file = 'MARS_script_data_file.csv' # This defines a variable the location of the data_file on the disk

    ### Run the R scripts to get the cuts and fitted temp values
    command = 'Rscript'
    path2CUTS = 'getMARScuts.R'
    path2FITTED = 'getMARSfitted.R'

    cmd_CUTS = [command, path2CUTS, data_file]
    cmd_CUTS = [command, path2CUTS]
    cmd_FITTED = [command, path2FITTED, data_file]

    # check_output will run the command and store to result
    MARScuts = subprocess.getoutput(cmd_CUTS, universal_newlines=True, shell=True)
    MARScuts = subprocess.Popen(cmd_CUTS)
    MARScuts = subprocess.check_output(cmd_CUTS, universal_newlines=True, stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL)
    MARScuts = subprocess.check_output(cmd_CUTS, universal_newlines=True, stdin=subprocess.DEVNULL)
    MARScuts = subprocess.check_output(cmd_CUTS, universal_newlines=True, stderr=subprocess.STDOUT, shell=True)
    MARSfitted = subprocess.check_output(cmd_FITTED, universal_newlines=True, stderr=subprocess.STDOUT, shell=True, bufsize=-5000)
    
    substring = 'Loading required package:'
    substring = 'Loading required package: class\nLoaded mda 0.4-7\n\n'
    # Once the data has been read in by the R scripts, delete the data_file
    #os.remove('MARS_script_data_file.csv')

    # Turn getMARS... output string into array
    MARScuts = MARScuts.split(' ')
    MARSfitted = MARSfitted.split(' ')

    # convert the string entries to numeric
    MARScuts = pd.DataFrame({'logHours':list(map(float, MARScuts))})
    MARSfitted = pd.DataFrame({'MARStemp':list(map(float, MARSfitted))})[startIndex:]
    
    # Add the logHours column to the MARSfitted df
    MARSfitted['logHours'] = TRTobj.data.logHours[startIndex:]
    
    #### Update the MARScuts dataframe and and ILSM results ####
    # Find indicies of the cuts in logHours from the logHours array and add to dataframe
    MARScuts['startIndex'] = MARScuts.logHours.apply(lambda x: find_nearest_idx(TRTobj.data.logHours, x))

    # Convert each of the logHours to hours and add to the dataframe labelled startHour
    MARScuts['startHour'] = MARScuts.startIndex.apply(lambda x: TRTobj.data.Hours[x])

    # Find the end index for each segment
    MARScuts['endIndex'] = MARScuts.startIndex[1:].append(pd.Series(len(TRTobj.data.Hours))).reset_index(drop=True)
    MARScuts

    # Subset the data between the start and end indices and fit LM to each
    from sklearn import linear_model

    ks_out=[]
    Rb_out=[]
    gradient = []
    intercept = []
    meanQ = []
    for ix, i in enumerate(MARScuts.startIndex):
        lm=linear_model.LinearRegression(fit_intercept=True) #create a linear model object
        lm.fit(np.array(TRTobj.data.logHours)[i:MARScuts.endIndex[ix],np.newaxis],
               np.array(TRTobj.data.Temp_ave)[i:MARScuts.endIndex[ix]]) #fit it to a subset of the data
        m=float(lm.coef_) #the slope of the line
        c=float(lm.predict(0)) #the y-intercept of the line
        Q=np.mean(TRTobj.data.Power[i:MARScuts.endIndex[ix]])
        ks=(Q/(4*np.pi*m*TRTobj.H))
        alpha=ks/(TRTobj.Cpvs)*3600 # Thermal diffusivity in m2/hr
        Rb=(1/(4*np.pi*ks)*(abs(c - TRTobj.T0)/m - np.log((4*alpha)/(1.78*TRTobj.rb**2)))) 
        ks_out.append(ks)
        Rb_out.append(Rb)
        gradient.append(m)
        intercept.append(c)
        meanQ.append(Q)

    MARScuts['ks'] = ks_out
    MARScuts['Rb'] = Rb_out
    MARScuts['m'] = gradient
    MARScuts['c'] = intercept
    MARScuts['meanQ'] = meanQ

    return({'Cuts': MARScuts, 'Fitted': MARSfitted})
