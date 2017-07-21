#!/usr/bin/env python

'''This function takes in a three column text file demarcating the 
pointing limits of a telescope. Format of the text file is
alt, az_lower, alt_upper. It returns two 1-D interpolation functions,
one for the lower limits and one for the upper limits.'''

import numpy as np
from scipy.interpolate import interp1d

def read_scope_limits(infile):
    data = np.loadtxt(infile, comments = '#').T
    az = data[0]
    alt_l = data[1]
    alt_u = data[2]
    interp_l = interp1d(az, alt_l)
    interp_u = interp1d(az, alt_u)
    scope_limits = [interp_l, interp_u]
    return scope_limits
    
    
## test case    
#il, iu = read_scope_limits('keck2_limits.txt')
#print(iu([0,100,200,300]))