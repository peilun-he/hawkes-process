setwd("/Users/HPL/Documents/GitHub/hawkes-process/")
source("simulate_marked_hawkes.R")

kernel1 <- function(t, times, a=3, b=7) {
  tau <- t - times
  tau <- tau[tau > 0]
  return(sum(a * exp(-b * tau)))
}

kernel2 <- function(t, times, a=5, b=6) {
  tau <- t - times
  tau <- tau[tau > 0]
  return(sum(a * exp(-b * tau)))
}

kernel12 <- function(t, times, a=5, b=5) {
  tau <- t - times
  tau <- tau[tau > 0]
  return(sum(a * exp(-b * tau)))
}

kernel21 <- function(t, times, a=5, b=5) {
  tau <- t - times
  tau <- tau[tau > 0]
  return(sum(a * exp(-b * tau)))
}

kernel <- matrix(list(kernel1, NA, NA, kernel1), nrow = 2, byrow = TRUE)
mu <- c(0.1, 0.2)
t_max <- 100

set.seed(1111)

aa <- simulate_marked_hawkes(t_max = t_max, mu = mu, kernel = kernel, show_info = TRUE)
lambda <- aa$lambda
times <- aa$event_times

plot(1: t_max, lambda[[1]], type = "l", xlab = "Time", ylab = "Intensity 1", ylim = c(0, max(lambda[[1]])))
plot(1: t_max, lambda[[2]], type = "l", xlab = "Time", ylab = "Intensity 2", ylim = c(0, max(lambda[[2]])))

################
##### Test #####
################
setwd("/Users/HPL/Documents/GitHub/hawkes-process/")
source("simulate_marked_hawkes.R")

alpha <- 3
beta <- 7
mu <- 0.1
t_max <- 100
p <- 0
n_mc <- 200

kernel1 <- function(t, times, a=alpha, b=beta) {
  tau <- t - times
  tau <- tau[tau > 0]
  return(sum(a * exp(-b * tau)))
}

counts1 <- numeric(n_mc)
counts2 <- numeric(n_mc)

set.seed(1111)

for (i in 1: n_mc) {
  sim1 <- simulate_marked_hawkes(t_max = t_max, mu = mu, kernel = kernel1)
  sim2 <- simulate_marked_hawkes(t_max = t_max, mu = mu, kernel = kernel1)
  counts1[i] <- length(sim1$event_times[[1]])
  counts2[i] <- length(sim2$event_times[[1]])
}

stat1 <- t_max^p * ( mean(counts1/t_max) - mu / (1 - alpha/beta) )
stat2 <- t_max^p * ( mean(counts2/t_max) - mu / (1 - alpha/beta) )
round(stat1, 3)
round(stat2, 3)
mu / (1 - alpha/beta)



