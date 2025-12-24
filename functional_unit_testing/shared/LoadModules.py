import numpy as np
import matplotlib.pyplot as plt
import json
import sys
sys.path.insert(0,'py_src')
from PyF90Utils import c8, ci, cchar, c8_arr, ci_arr, ccharnb
import CtypesInit
exec(open("py_src/CtypesInit.py").read())

f90 = f90_modules('bld/')
