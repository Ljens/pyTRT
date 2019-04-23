# Test the ILSM_sim method
import sys

# Requires two command line arguments to be passed
# e.g. >>> $ python3 test_sum.py 50 100
# returns the borehole fluid temperature at 100 seconds assuming input power of 50W/m
q = float(sys.argv[1]) # Input power in W/m, constant value
t = float(sys.argv[2]) # # Input times in seconds

# Import the pyTRT class
from pyTRT.TRT_class import TRT

#print(os.getcwd())

# Read in params to make the TRT class
# They are in a .txt file, so need to parse it
params = {}
with open("pyTRT/tests/test_param_input.txt") as myfile:
    for line in myfile:
        name, var = line.partition("=")[::2]
        params[name.strip()] = float(var)

# put a unittest here to check the key names are correc
# Creat TRT object with input params
TRTobj=TRT('TRT test',params)

# Return the simulation result
print(TRTobj.ILSM_sim(q,t))

