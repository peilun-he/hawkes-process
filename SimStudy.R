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

kernel <- matrix(list(kernel1, kernel12, kernel21, kernel2), nrow = 2, byrow = TRUE)
mu <- c(0.1, 0.2)
t_max <- 10

set.seed(1111)

aa <- simulate_marked_hawkes(t_max = t_max, mu = mu, kernel = kernel)
lambda <- aa$lambda
times <- aa$event_times

plot(1: t_max, lambda[[1]], type = "l", xlab = "Time", ylab = "Intensity 1", ylim = c(0, max(lambda[[1]])))
plot(1: t_max, lambda[[2]], type = "l", xlab = "Time", ylab = "Intensity 2", ylim = c(0, max(lambda[[2]])))


