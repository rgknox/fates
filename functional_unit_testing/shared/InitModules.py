import numpy as np
import matplotlib.pyplot as plt
import json
import sys
sys.path.insert(0,'py_src')
from PyF90Utils import c8, ci, cchar, c8_arr, ci_arr, ccharnb
from CtypesInit import f90_modules # Directly import f90_modules from CtypesInit
from ctypes import *
import code  # For development: code.interact(local=dict(globals(), **locals()))


f90 = f90_modules('bld/')

parameterfile = '../../parameter_files/fates_params_default.json'

with open(parameterfile, 'r') as file:
    params = json.load(file)
    
iret = f90.json_setinval(c8(9.9e32+10))
iret = f90.json_setloginit(ci(6))
iret = f90.wrapjson_read(*ccharnb(parameterfile.strip()))

iret = f90.getparams()


