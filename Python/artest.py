#!/home/nltom2/.virtualenvs/Python/bin/python3
from edward.models import Normal, InverseGamma
import tensorflow as tf
import edward as ed
import numpy as np

mu_true = 0.0
beta_true = 0.9
noise_obs = 0.1
T = 32
N = 50
x_true = np.random.randn(T,N)*noise_obs
for t in range(1, T):
    x_true[t] += beta_true*x_true[t-1]


mu = Normal(mu=0., sigma=10.0)
beta = Normal(mu=0., sigma=2.0)
noise_proc = tf.constant(0.1) #InverseGamma(alpha=1.0, beta=1.0)
noise_obs = tf.constant(0.1) #InverseGamma(alpha=1.0, beta=1.0)
x = [0] * T
x[0] = Normal(mu=tf.ones([N])*mu, sigma=tf.ones([N])*10.0)
for n in range(1, T):
    x[n] = Normal(mu=tf.ones([N])*mu + beta * x[n-1], sigma=tf.ones([N])*noise_proc)


qmu = Normal(mu=tf.Variable(tf.random_normal([])),
             sigma=tf.nn.softplus(tf.Variable(tf.random_normal([]))))
qbeta = Normal(mu=tf.Variable(tf.random_normal([])),
               sigma=tf.nn.softplus(tf.Variable(tf.random_normal([]))))
inference_vb = ed.KLqp({beta: qbeta, mu:qmu},
                       data={xt: xt_true for xt, xt_true in zip(x, x_true)})
inference_vb.run()
