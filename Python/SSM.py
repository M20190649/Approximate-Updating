#!/usr/bin/env python

#import edward as ed
#from edward.models import Normal
import matplotlib.pyplot as plt
#import matplotlib.gridspec as gridspec
import numpy as np
#import os
import tensorflow as tf
import numpy.random as rn

T = 50
y = np.empty(T)
x = np.empty(T+1)

xmu = 0
xsd = 1
ymu = 2
xar = 0.5
ysd = 1

x[0] = rn.normal(xmu, xsd)
for t in range(T):
    x[t+1] = rn.normal(xmu + xar*x[t], xsd)
    y[t] = rn.normal(x[t+1], ysd)

