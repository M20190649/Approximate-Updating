library(tensorflow)
sess = tf$Session()

file_writer = tf$summary$FileWriter('logs', sess$graph)
error <- rnorm(100)
mu <- 3
y_data <- mu + 0.1 * error
Mu <- tf$Variable(tf$zeros(shape(1L)))
Sd <- tf$Variable(tf$ones(shape(1L)))
X <- tf$distributions$Normal(loc=tf$zeros(shape(1L)), scale=tf$ones(shape(1L)))

eps <- tf$random_normal(shape(1L), mean=0, stddev=1)
epsVal <- eps$eval(session = sess)
y <- Mu + Sd * epsVal * tf$ones(100) + 0.1* error
logP <- - tf$divide(tf$reduce_sum(tf$pow(tf$subtract(y_data, y), 2)),  0.02)
logJ <- tf$log(Sd)
logQ <- - epsVal^2 / 2
elbo <- tf$subtract(tf$add(logP, logJ), logQ)

optimizer <- tf$train$GradientDescentOptimizer(0.1)
train <- optimizer$minimize(-elbo)  
sess <- tf$Session()
sess$run(tf$global_variables_initializer())

# Fit the line (Learns best fit is W: 0.1, b: 0.3)
for (step in 1:15) {
  sess$run(train)
  if (step %% 2 == 0)
    cat(step, "-", sess$run(Mu), sess$run(Sd), "\n")
}

theta <- tf$contrib$bayesflow$stochastic_tensor(tf$contrib$distributions)
