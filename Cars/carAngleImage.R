x <- seq(1, 3.5, 0.001)
y <- -(x - 3)^2 + 5
y[x > 3] <- 5
x[y < 4.5] <- 3 - sqrt(0.5)

xa1 <- seq(2.11, 2.15, 0.001)
ya1 <- sqrt(0.15^2 - (xa1 - 2)^2) + 4.25

xa2 <- seq(2.905, 3.1, 0.0001)
ya2 <- sqrt(0.1^2 - (xa2 - 3)^2) + 4.625

library(ggplot2)
ggplot() + geom_path(aes(x, y), size = 1.25) + 
  geom_abline(aes(slope = 1, intercept = 2.25), linetype = 2, size = 1.25) + 
  geom_abline(aes(slope = -0.25, intercept = 5.375), linetype = 3, size = 1.25) + 
  geom_point(aes(2.5, 4.75), colour = 'red', size = 7) + 
  geom_point(aes(2.3, 4.8), colour = 'blue', size = 7) + 
  geom_line(aes(x = c(2, 2.15), y = c(4.25, 4.25)), linetype = 2) + 
  geom_line(aes(xa1, ya1), linetype = 2) + 
  geom_text(aes(x = 2.1, y = 4.28, label = "lambda"), size = 8,  parse = TRUE) + 
  geom_line(aes(x = c(3, 3.1), y = c(4.625, 4.625)), linetype = 3) + 
  geom_line(aes(xa2, ya2), linetype = 3) + 
  geom_text(aes(x = 3, y = 4.66, label = 'psi'), size = 8, parse = TRUE) + 
  ylim(4.2, 5) + xlim(2, 3.3) + 
  theme_bw() + 
  theme(axis.text = element_blank(),
        axis.ticks = element_blank())


xa3 <- seq(1.105, 1.15, 0.001)
ya3 <- sqrt(0.15^2 - (xa3 - 1)^2) + 1

ggplot() + geom_path(aes(x = c(1, 2, 2, 1), y = c(1, 1, 2, 1))) + 
  geom_point(aes(x = 1, y = 1), colour = 'red', size = 5) + 
  geom_point(aes(x = 2, y = 2), colour = 'blue', size = 5) + 
  geom_text(aes(x = 0.85, y = 0.8, label = 'x[t-1]'), size = 10, parse = TRUE) + 
  geom_text(aes(x = 0.93, y = 0.8, label = ','), size = 10) +
  geom_text(aes(x = 1.02, y = 0.8, label = 'y[t-1]'), size = 10, parse = TRUE) +
  geom_text(aes(x = 2.1, y = 2, label = 'x[t]'), size = 10, parse = TRUE) + 
  geom_text(aes(x = 2.15, y = 2, label = ','), size = 10) +
  geom_text(aes(x = 2.23, y = 2, label = 'y[t]'), size = 10, parse = TRUE) +
  geom_line(aes(xa3, ya3), linetype = 2) + 
  geom_text(aes(x = 1.18, y = 1.09, label = "delta[t]"), size = 10, parse = TRUE) +
  geom_text(aes(x = 1.4, y = 1.5, label = "v[t]"), size = 10, parse = TRUE) +
  geom_text(aes(x = c(0.85, 1, 2.14, 2.26), y = c(0.84, 0.85, 2.05, 2.05), label = "*"), size = 10) + 
  ylim(0.75, 2.1) + xlim(0.6, 2.3) + 
  theme_bw() + 
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title=element_text(size=30)) + 
  labs(x = "x", y = "y")


