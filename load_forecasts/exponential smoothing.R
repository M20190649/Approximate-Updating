library(ggplot2)
sim_exp_smoothing = function(params, lags, initial, T, sigmasq){
  m = length(params)
  dim = sum(lags) - m + 1
  if(length(lags) != m){
    stop("Parameter and lag vectors must be the same length")
  }
  if(length(initial) != dim){
    stop("Length of initial state vector must equal sum of lag vector minus number of latent states plus one")
  }
  if(any(params >= 1 | params <= 0)){
    stop("All parameters must be between zero and one")
  }
  states = matrix(0, ncol = T+1, nrow = dim)
  states[,1] = initial
  y = vector(length = T+1)
  newstates = c(1, 1+cumsum(lags))[-(m+1)]
  for(t in 2:(T+1)){
    y[t] = rnorm(1, 0, sqrt(sigmasq)) + states[1, t-1] - sum(states[(1:sum(lags))[-cumsum(lags)],t-1])
    states[1, t] = (1-params[1])*states[1, t-1] + params[1]*(y[t] + sum(states[(1:sum(lags))[-cumsum(lags)],t-1])) 
    
    
    for(j in 2:m){
      states[newstates[j], t] = -sum(states[(cumsum(lags)[j-1]+1):(cumsum(lags)[j]-1),t-1]) +
                                       params[j]*(y[t] - states[1, t-1] + sum(states[(1:sum(lags))[-cumsum(lags)], t-1]))
    }
    for(j in (1:(sum(lags)-m+1))[-newstates]){
      states[j, t] = states[j-1, t-1]
    }
  }
  out = list(y = y[2:(T+1)], states = states[newstates,2:(T+1)])
  return(out)
}

alpha = c(0.3, 0.5)
lags = c(1, 7)
initial = c(0, 1, 2, 6, 1, 1, -3)
T = 250
sigmasq = 1
exp_data = sim_exp_smoothing(alpha, lags, initial, T, sigmasq)
ggplot() + geom_line(aes(1:T, exp_data$y), colour = "red") + 
  geom_line(aes(1:T, exp_data$states[1,]), colour = "black") + geom_line(aes(1:T, exp_data$states[2, ]), colour = "blue")

##MCMC Setup
m = length(lags)
dim = sum(lags) - m +1
y = exp_data$y
x = c(1, rep(-1, dim-1))
Trans = matrix(0, dim, dim)
Trans[1,1] = 1
newstates = c(1, 2+cumsum(lags-1)[-m])
for(i in 2:m){
  if(i < m){
    Trans[newstates[i], newstates[i]:(newstates[i+1]-1)] = -1
  } else {
    Trans[newstates[m], newstates[m]:dim] = -1
  }
}
for(i in (1:dim)[-newstates]){
  Trans[i, i-1] = 1
}


rep = 5000
test = MCMC_exp_smoothing(y, rep, x, Trans, lags, 0.05)


biunif = function(x) {x[1]*x[2]}

resolution = 11
support = seq(0, 1, length.out = resolution)
grid = expand.grid(support, support)
grid$dens = apply(grid, 1, cond_density_only, y = y, Trans = Trans, x = x, lags = lags, sigmasq = sigmasq, initial = initial)
grid$dens = grid$dens / sum(grid$dens)
grid$expdens = exp(grid$dens)

ggplot(grid) + geom_tile(aes(Var1, Var2, fill = expdens))
grid[which.max(grid$dens),]

densitymatrix = matrix(grid$dens, resolution)
marginalx = colSums(densitymatrix)
XCDF = cumsum(marginalx)


draws = replicate(10000, drawing())
draws = t(draws)

rownames(densitymatrix) = support
colnames(densitymatrix) = support




x_tilde = matrix(0, T, dim)
y_tilde = rep(0, T)
b_bar = matrix(0, dim, T)
alpha_full = rep(0, dim)
csum_lags = c(1, 2)
alpha_full[csum_lags] = alpha
D = Trans - alpha_full %*% t(x)
x_tilde[1,] = t(x)
b_bar[,1] = alpha_full * y[1]
y_tilde[1] = y[1]
for(t in 2:T){
  b_bar[,t] = D %*% b_bar[,t-1] + alpha_full * y[t]
  x_tilde[t,] = x_tilde[t-1,] %*% D 
  y_tilde[t] = y[t] - t(x) %*% b_bar[,t-1] 
}
xtxinv = solve(t(x_tilde) %*% x_tilde)
beta0hat = xtxinv %*% t(x_tilde) %*% y
ssquared = (t(y_tilde - x_tilde %*% beta0hat) %*% (y_tilde - x_tilde %*% beta0hat)) / (T - size)
log(sqrt(det(xtxinv)))   -(T - size)/2 * log(ssquared)





#Progress:
#Fix problem in eigenvalue of D matrix (coding error)
#Write simulation and MCMC code (in progress - alpha draws)
#Discover problem with XTX being singular
#Reparameterise seasonal effects to have equivalent to XTX that is not singluar
#Fix old code to work with new parameterisation
#Wrote inverse CDF sampler for two and three dimensional alpha

#To do:
#Handle extremely small densities - via conditional and MH is not working - gets stuck very easily
#Run and see what happens
#Hope to fit marginals to easy to use densities
#Figure out copula structure and fit
#Write Variational Bayes algorithm (and hope distribution allows reparameterisation)

