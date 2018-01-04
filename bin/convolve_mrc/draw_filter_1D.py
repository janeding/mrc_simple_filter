#!/usr/bin/env python
import sys
import matplotlib.pyplot as plt
import numpy as np
from math import *

class InputError(Exception):
    """ A generic exception object containing a string for error reporting.
        (Raising this exception implies that the caller has provided
         a faulty input file or argument.)

    """
    def __init__(self, err_msg):
        self.err_msg = err_msg
    def __str__(self):
        return self.err_msg
    def __repr__(self):
        return str(self)

A = 2.7
B = -1.7
a = 4.0
b = 5.0

if sys.argv[1] == '-gauss':
    if len(sys.argv) <= 4:
        raise InputError('Error: Expected numbers arguments following \"'+sys.argv[1]+'\"\n')
    A = float(sys.argv[1])
    a = float(sys.argv[2])
    m = 2.0
    n = 2.0
    B = 0.0
    b = 1.0

elif sys.argv[1] == '-gaussgen':
    if len(sys.argv) <= 4:
        raise InputError('Error: Expected numbers arguments following \"'+sys.argv[1]+'\"\n')
    A = float(sys.argv[1])
    a = float(sys.argv[2])
    m = 2.0
    n = 2.0
    B = 0.0
    b = 1.0
    # By default, we assume the user wants us to plot:
    #   h(r) = A*exp(-0.5*(r/a)^2) - B*exp(-0.5*(r/b)^2)
    # Since this script instead plots
    #   h(r) = A*exp(-(r/a)^m) - B*exp(-(r/b)^n)
    # so we absorb the factor of sqrt(0.5) into "a" and "b":
    a *= sqrt(2.0)
    if len(sys.argv) >= 4:
        # if the user manually specified (m and n)
        # then we assume the user wants us to plot this instead:
        #   h(r) = A*exp(-(r/a)^m) - B*exp(-(r/b)^n)
        m = float(sys.argv[3])
        n = float(sys.argv[4])

elif sys.argv[1] == '-dog':
    if len(sys.argv) <= 4:
        raise InputError('Error: Expected numbers arguments following \"'+sys.argv[1]+'\"\n')
    A = float(sys.argv[1])
    B = float(sys.argv[2])
    a = float(sys.argv[3])
    b = float(sys.argv[4])
    m = 2.0
    n = 2.0

elif sys.argv[1] == '-doggen':
    if len(sys.argv) <= 4:
        raise InputError('Error: Expected numbers arguments following \"'+sys.argv[1]+'\"\n')
    A = float(sys.argv[1])
    B = -float(sys.argv[2])
    a = float(sys.argv[3])
    b = float(sys.argv[4])
    m = 2.0
    n = 2.0
    # By default, we assume the user wants us to plot:
    #   h(r) = A*exp(-0.5*(r/a)^2) - B*exp(-0.5*(r/b)^2)
    # Since this script instead plots
    #   h(r) = A*exp(-(r/a)^m) - B*exp(-(r/b)^n)
    # so we absorb the factor of sqrt(0.5) into "a" and "b":
    a *= sqrt(2.0)
    b *= sqrt(2.0)
    if len(sys.argv) == 7:
        # if the user manually specified (m and n)
        # then we assume the user wants us to plot this instead:
        #   h(r) = A*exp(-(r/a)^m) - B*exp(-(r/b)^n)
        m = float(sys.argv[5])
        n = float(sys.argv[6])


width = max(a,b)
r_i = np.arange(-ceil(4*width), ceil(4*width), 1.0)
p_i = np.array([A*exp(-(abs(r)/a)**m)-B*exp(-(abs(r)/b)**n) for r in r_i])
#plt.plot(r_i, p_i)
plt.step(r_i+0.5, p_i)
plt.show()
