# Required Packages
import os
#import numpy as np
#import pandas as pd
#import matplotlib as mpl
#import matplotlib.pyplot as plt
#import seaborn as sns; sns.set_style("whitegrid")
print(os.getcwd())

# Test importing pyTRT as a module
from pyTRT.TRT_class import TRT

# Read in params to make the TRT class
# They are in a .txt file, so need to parse it
params = {}
with open("pyTRT/tests/test_param_input.txt") as myfile:
    for line in myfile:
        name, var = line.partition("=")[::2]
        params[name.strip()] = float(var)

# Creat TRT object with input params
TRTobj=TRT('TRT test',params)

# Import data
TRTobj.import_geocube('pyTRT/tests/TRT1_Footscray.csv') # Add measured data

# Run the ILSM_linreg method
# This fits the ILSM with linear regression to the data and returns a dictionary with predictions for
#   - Thermal conducitivity - ks
#   - GHE Thermal Resistance - Rb
# along with the fitting paramters.
print(TRTobj.ILSM_linreg(5,50))


