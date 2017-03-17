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
y = [0] * (T)
x0 = Normal(mu=0.0, sigma=xsd / tf.sqrt(1-ar**2))
for t in range(T):
    if(t == 0):
        x[t] = Normal(mu=ar*x0, sigma=xsd)
    else:
        x[t] = Normal(mu=ar*x[t-1], sigma=xsd)
    y[t] = Normal(mu=x[t]+b, sigma=ysd)

#Set up approximation
qAr = Normal(mu = tf.Variable(tf.random_normal([])),
            sigma=tf.nn.softplus(tf.Variable(tf.random_normal([]))))
qB = Normal(mu = tf.Variable(tf.random_normal([])),
            sigma=tf.nn.softplus(tf.Variable(tf.random_normal([]))))
qYsd = InverseGamma(
    alpha=tf.nn.softplus(tf.Variable(tf.random_normal([]))),
    beta=tf.nn.softplus(tf.Variable(tf.random_normal([]))))
qXsd = InverseGamma(
    alpha=tf.nn.softplus(tf.Variable(tf.random_normal([]))),
    beta=tf.nn.softplus(tf.Variable(tf.random_normal([]))))
qX0 = Normal(mu = tf.Variable(tf.random_normal([])),
            sigma=tf.nn.softplus(tf.Variable(tf.random_normal([]))))
qX = [0] * (T)
for t in range(T):
    qX[t] = Normal(mu = tf.Variable(tf.random_normal([])),
            sigma=tf.nn.softplus(tf.Variable(tf.random_normal([]))))

# Inference
parameters = {
    ar: qAr,
    b: qB,
    ysd: qYsd,
    xsd: qXsd,
    x0: qX0
}
states = {x: qX for x, qX in zip(x, qX)}
parameters.update(states)
data = {yt: yData for yt, yData in zip(y, yData)}
inference = ed.KLqp(parameters, data)
inference.run()
