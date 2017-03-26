import tensorflow as tf
import edward as ed
import numpy as np
from edward.models import Normal

a = tf.train.GradientDescentOptimizer
a.compute_gradients()

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
        hTrue[t] = np.random.normal(loc=gammaTrue / (1 - phiTrue),
                                   scale=sigmaTrue / np.sqrt(1 - phiTrue ** 2))
    else:
        hTrue[t] = np.random.normal(loc=gammaTrue + phiTrue * hTrue[t - 1],
                                    scale=sigmaTrue)
    yTrue[t] = np.random.normal(loc=0, scale=np.exp(hTrue[t] / 2))

# Prior distributions
gamma = Normal(mu=0.0, sigma=5.0)
phi = Normal(mu=0.0, sigma=0.5)
lsigma = Normal(mu=0.0, sigma=1.0)
sigma = tf.exp(lsigma)
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
qlSigma = Normal(mu=tf.Variable(tf.random_normal([])),
              sigma=tf.nn.softplus(tf.Variable(tf.random_normal([]))))
qH = [0] * T
for t in range(T):
    qH[t] = Normal(mu=tf.Variable(tf.random_normal([])),
                   sigma=tf.nn.softplus(tf.Variable(tf.random_normal([]))))

# Inference
parameters = {
    gamma: qGamma,
    phi: qPhi,
    lsigma: qlSigma,
}
states = {h: qH for h, qH in zip(h, qH)}
parameters.update(states)
data = {yt: yTrue for yt, yTrue in zip(y, yTrue)}
inference = ed.KLqp(parameters, data=data)
inference.run()
