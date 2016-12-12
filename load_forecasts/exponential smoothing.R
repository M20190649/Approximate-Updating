library(ggplot2)
library(gridExtra)
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

resolution = 500
xsup = seq(0.2, 0.45, length.out = resolution)
ysup = seq(0.4, 0.65, length.out = resolution)
grid = expand.grid(xsup, ysup)
grid$dens = apply(grid, 1, density_only, y = y, Trans = Trans, x = x, lags = lags)
grid$dens[grid$dens == "NaN"] = 0
grid$dens = grid$dens / sum(grid$dens*(xsup[2]-xsup[1])*(ysup[2]-ysup[1]))

densityplot = ggplot(grid) + geom_tile(aes(Var1, Var2, fill = dens))
densityplot
grid[which.max(grid$dens),]

densitymatrix = matrix(grid$dens, resolution, byrow = TRUE)
rownames(densitymatrix) = ysup
colnames(densitymatrix) = xsup
marginalx = colSums(densitymatrix)
XCDF = cumsum(marginalx)

densarray = array(0, dim = c(500, 500, 1))
densarray[,,1] = densitymatrix
testdraws = replicate(10000, draw_alpha(XCDF, xsup, ysup, densarray))
draws = as.data.frame(matrix(testdraws, ncol = 2, byrow = TRUE))
drawsplot = ggplot(draws, aes(V1, V2)) + geom_density2d()
grid.arrange(densityplot, drawsplot, ncol = 2)

draws = data.frame()
store = MC_exp_smoothing(y, 100, x, Trans, lags, XCDF, xsup, ysup, densarray, 1:10)
draws = rbind(draws, store)
colnames(draws) = c("Alpha", "Delta", "SigmaSq", "In1", "In2", "In3", "In4", "In5", "In6", "In7")
apply(draws, 2, mean)

#Progress:
#Fix problem in eigenvalue of D matrix (coding error)
#Write simulation and MCMC code (in progress - alpha draws)
#Discover problem with XTX being singular
#Reparameterise seasonal effects to have equivalent to XTX that is not singluar
#Fix old code to work with new parameterisation
#Wrote inverse CDF sampler for two and three dimensional alpha
#Fixed bug leading to very small densities
#Calculate CDFs

#To do:
#Run and see what happens
#Hope to fit marginals to easy to use densities
#Figure out copula structure and fit
#Write Variational Bayes algorithm (and hope distribution allows reparameterisation)

