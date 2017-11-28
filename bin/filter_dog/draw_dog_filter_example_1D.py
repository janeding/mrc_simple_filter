#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
from math import *
A = 2.7
B = -1.7
s = 4.0
t = 5.0
a = 4.0
b = 4.0
r_i = np.arange(-20.0, 21.0, 1.0)
p_i = np.array([A*exp(-(abs(r)/s)**a)+B*exp(-(abs(r)/t)**b) for r in r_i])
#plt.plot(r_i, p_i)
plt.step(r_i+0.5, p_i)
plt.show()
