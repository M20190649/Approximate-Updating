#!/home/nltom2/.virtualenvs/Python/bin/python3


import edward as ed
from edward.models import Normal, Uniform, InverseGamma
import matplotlib.pyplot as plt
#import matplotlib.gridspec as gridspec
import numpy as np
#import os
import tensorflow as tf
import numpy.random as rn


#Setting up parameters for simulated data
T = 50
xsdTrue = 1
bTrue = 2
arTrue = 0.5
ysdTrue = 1

#Simulating data
yData = np.empty(T)
xData = np.empty(T+1)
xData[0] = rn.normal(loc=0, scale=xsdTrue/np.sqrt(1 - arTrue**2))
for t in range(T):
    xData[t+1] = rn.normal(loc=arTrue*xData[t], scale=xsdTrue)
    yData[t] = rn.normal(loc=xData[t+1] + bTrue, scale=ysdTrue)


#Create Model
b = Normal(mu=0.0, sigma=5.0)
ar = Uniform(-1.0, 1.0)
xsd = InverseGamma(alpha=1.0, beta=1.0)
ysd = InverseGamma(alpha=1.0, beta=1.0)

x = [0] * (T)
x0 = Normal(mu=0.0, sigma=xsd / tf.sqrt(1-ar**2))
for t in range(T):
    if(t == 0):
        x[t] = Normal(mu=ar*x0, sigma=xsd)
    else:
        x[t] = Normal(mu=ar*x[t-1], sigma=xsd)

y = Normal(mu=x+b, sigma=ysd)


#Set up approximation
qar = Normal(mu = tf.Variable(tf.random_normal([1])),
            sigma=tf.nn.softplus(tf.Variable(tf.random_normal([1]))))
qb = Normal(mu = tf.Variable(tf.random_normal([1])),
            sigma=tf.nn.softplus(tf.Variable(tf.random_normal([1]))))
qysd = InverseGamma(
    alpha=tf.nn.softplus(tf.Variable(tf.random_normal([1]))),
    beta=tf.nn.softplus(tf.Variable(tf.random_normal([1]))))
qxsd = InverseGamma(
    alpha=tf.nn.softplus(tf.Variable(tf.random_normal([1]))),
    beta=tf.nn.softplus(tf.Variable(tf.random_normal([1]))))
qx0 = Normal(mu = tf.Variable(tf.random_normal([1])),
            sigma=tf.nn.softplus(tf.Variable(tf.random_normal([1]))))
qx = Normal(mu = tf.Variable(tf.random_normal([T])),
            sigma=tf.nn.softplus(tf.Variable(tf.random_normal([T]))))

inference = ed.KLqp({ar: qar, b: qb, xsd: qxsd, ysd: qysd, x0: qx0, x: qx}, data={y: yData})
inference.run(n_iter=1000, n_samples=5)
