setwd("/Users/HPL/Documents/GitHub/hawkes-process/")
source("simulate_marked_hawkes_v2.R")

mu <- c(0.1, 0.2)
alpha <- matrix(c(3, 5, 5, 5), nrow = 2, byrow = TRUE)/10
beta <- matrix(c(7, 5, 5, 6), nrow = 2, byrow = TRUE)
#mu <- 0.1
#alpha <- 3/10
#beta <- 7/10
t_max <- 100

set.seed(1111)

aa <- simulate_marked_hawkes(t_max = t_max, mu = mu, alpha = alpha, beta = beta, show_info = FALSE)
lambda <- aa$lambda
times <- aa$event_times

plot(1: t_max, lambda[1, ], type = "l", xlab = "Time", ylab = "Intensity 1")
plot(1: t_max, lambda[2, ], type = "l", xlab = "Time", ylab = "Intensity 2")

################
##### Test #####
################
setwd("/Users/HPL/Documents/GitHub/hawkes-process/")
source("simulate_marked_hawkes_v2.R")

t_max <- 100
p <- 0.5
n_mc <- 200
d <- 1

if (d == 1) {
  alpha <- 0.3
  beta <- 7
  mu <- 0.1
} else if (d == 2) {
  alpha <- matrix(c(0.3, 0.1, 0.1, 0.5), nrow = 2, byrow = TRUE)
  beta <- matrix(c(7, 8, 8, 5), nrow = 2, byrow = TRUE)
  mu <- matrix(c(0.1, 0.3), nrow = 2)
}

counts1 <- numeric(n_mc)
counts2 <- numeric(n_mc)

set.seed(1111)

for (i in 1: n_mc) {
  sim1 <- simulate_marked_hawkes(t_max = t_max, mu = mu, alpha = alpha, beta = beta)
  sim2 <- simulate_marked_hawkes(t_max = t_max, mu = mu, alpha = alpha, beta = beta)
  counts1[i] <- length(sim1$event_times[[1]])
  counts2[i] <- length(sim2$event_times[[1]])
}

stat1 <- t_max^p * ( mean(counts1/t_max) - solve(diag(d) - alpha/beta) %*% mu )
stat2 <- t_max^p * ( mean(counts2/t_max) - solve(diag(d) - alpha/beta) %*% mu )
round(stat1, 3)
round(stat2, 3)
mu / (1 - alpha/beta)

# example


