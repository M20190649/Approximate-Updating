#!/home/nltom2/.virtualenvs/Python/bin/python3
import tensorflow as tf
import edward as ed
import numpy as np
from edward.models import Normal, InverseGamma, Uniform

# Model Parameters
T = 50
sigmaTrue = 0.5
gammaTrue = 0.2
phiTrue = 0.8

# Simulate Data
hTrue = [0] * T
yTrue = [0] * T
for t in range(T):
    if (t == 1):
        hTrue[t] = np.random.normal(loc=gammaTrue / (1 - phiTrue), scale=sigmaTrue / np.sqrt(1 - phiTrue ** 2))
    else:
        hTrue[t] = np.random.normal(loc=gammaTrue + phiTrue * hTrue[t - 1], scale=sigmaTrue)
    yTrue[t] = np.random.normal(loc=0, scale=np.exp(hTrue[t] / 2))

# Prior distributions
gamma = Normal(mu=0.0, sigma=5.0)
phi = Uniform(-1.0, 1.0)
sigma = InverseGamma(alpha=1.0, beta=1.0)
# Create h and y as lists
h = [0] * T
y = [0] * T
for t in range(T):
    if (t == 0):
        h[0] = Normal(mu=0.0, sigma=2.0)
    else:
        h[t] = Normal(mu=gamma + phi * h[t - 1], sigma=sigma)
    y[t] = Normal(mu=0.0, sigma=tf.exp(h[t] / 2))

# Create Approximating Distribution
qGamma = Normal(mu=tf.Variable(tf.random_normal([])),
                sigma=tf.nn.softplus(tf.Variable(tf.random_normal([]))))
qPhi = Normal(mu=tf.Variable(tf.random_normal([])),
              sigma=tf.nn.softplus(tf.Variable(tf.random_normal([]))))
qSigma = InverseGamma(
    alpha=tf.nn.softplus(tf.Variable(tf.random_normal([]))),
    beta=tf.nn.softplus(tf.Variable(tf.random_normal([]))))
qH = [0] * T
for t in range(T):
    qH[t] = Normal(mu=tf.Variable(tf.random_normal([])),
                   sigma=tf.nn.softplus(tf.Variable(tf.random_normal([]))))

# Inference
parameters = {
    gamma: qGamma,
    phi: qPhi,
    sigma: qSigma,
}
states = {h: qH for h, qH in zip(h, qH)}
parameters.update(states)
data = {yt: yTrue for yt, yTrue in zip(y, yTrue)}
inference = ed.KLqp(parameters, data)
inference.run()