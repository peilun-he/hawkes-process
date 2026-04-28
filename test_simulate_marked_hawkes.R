#############################################
##### Test 1: univariate Hawkes process #####
#############################################
kernel <- function(t, times, a=5, b=6) {
  tau <- t - times
  tau <- tau[tau > 0]
  return(sum(a * exp(-b * tau)))
}

mu <- 0.2

t_max <- 10

set.seed(1111)

results <- simulate_marked_hawkes(t_max = t_max, mu = mu, kernel = kernel)
results$lambda <- lapply(results$lambda, round, digits = 4)
results$event_times <- lapply(results$event_times, round, digits = 4)

exp_results <- list(lambda = list(c(0.2000, 0.2000, 0.2000, 0.2000, 2.8923, 
                                    0.2824, 0.2002, 0.2000, 0.2000, 0.2000)), 
                    event_times = list(c(4.8968, 5.3016))) # expected results 

tinytest::expect_equal(results, exp_results)

#############################################################
##### Test 2: multivariate independent Hawkes processes #####
#############################################################
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

kernel <- matrix(list(kernel1, NA, NA, kernel2), nrow = 2, byrow = TRUE)
mu <- c(0.1, 0.2)
t_max <- 10

set.seed(1111)

results <- simulate_marked_hawkes(t_max = t_max, mu = mu, kernel = kernel)
results$lambda <- lapply(results$lambda, round, digits = 4)
results$event_times <- lapply(results$event_times, round, digits = 4)

exp_results <- list(lambda = list(c(0.1000, 0.1000, 0.1000, 0.2327, 0.1001,
                                    0.1000, 0.1000, 0.1000, 0.1000, 0.1000),
                                  c(0.2000, 0.2000, 0.2000, 0.2000, 5.0248, 
                                    0.2120, 0.2667, 0.2002, 0.2000, 0.2000)), 
                    event_times = list(c(3.2645, 3.5344), 
                                       c(4.9941, 6.2805))
                    ) # expected results 

tinytest::expect_equal(results, exp_results)


#################################################
##### Test 3: multivariate Hawkes processes #####
#################################################
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

results <- simulate_marked_hawkes(t_max = t_max, mu = mu, kernel = kernel)
results$lambda <- lapply(results$lambda, round, digits = 4)
results$event_times <- lapply(results$event_times, round, digits = 4)

exp_results <- list(lambda = list(c(0.1000, 0.1000, 0.1000, 0.2327, 4.9536,
                                    0.1327, 0.2372, 0.1009, 0.1000, 0.1000),
                                  c(0.2000, 0.2000, 0.2000, 0.8139, 5.0289, 
                                    0.2120, 0.2667, 0.2002, 0.2000, 0.2000)), 
                    event_times = list(c(3.2645, 3.5344), 
                                       c(4.9941, 6.2805))
) # expected results 

tinytest::expect_equal(results, exp_results)




