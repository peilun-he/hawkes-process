simulate_marked_hawkes <- function(t_max, mu, alpha, beta, show_info = FALSE) {
  # Simulate multivariate marked Hawkes process: 
  #   lambda_{i}(t) = mu_{i} + \sum_{j=1}^{p} \sum_{t_k < t} alpha * exp(-beta * (t - t_k))
  #
  # Inputs: 
  #   t_max: scalar specifying the maximum time horizon.
  #   mu: numeric vector of baseline intensities.
  #   alpha: numeric matrix of excitation parameters. 
  #   beta: numeric matrix of decay rates. 
  #   show_info: logical; if TRUE, prints the probability of each event type.
  #
  # Outputs: 
  #   lambda: matrix of intensity trajectories for each dimension.
  #   event_times: list of event times for each dimension.
  
  p <- length(mu) # dimension of Hawkes process
  
  if (p == 1) {
    alpha <- matrix(alpha)
    beta <- matrix(beta)
  }
  
  # Check 
  if ( !is.numeric(mu) ) stop("mu must be a vector of numbers.")
  if ( !is.numeric(alpha) ) stop("alpha must be a matrix of numbers.")
  if ( !is.numeric(beta) ) stop("beta must be a matrix of numbers.")
  if ( nrow(alpha) != ncol(alpha) ) stop("alpha must be a square matrix.")
  if ( nrow(beta) != ncol(beta) ) stop("beta must be a square matrix.")
  if ( nrow(alpha) != p | nrow(beta) != p ) stop("The dimensions of mu, alpha, and beta are different.")
  if ( any(is.na(alpha)) ) alpha[ which(is.na(alpha)) ] <- 0
  if ( any(alpha < 0) ) stop("Excitation parameters must be non-negative. ")
  if ( any(beta <= 0) ) stop("Decay rates must be positive. ")
  if ( any(mu < 0) ) stop("Baseline intensities must be non-negative. ")
  t_max <- floor(t_max)
  
  t <- 0
  event_times <- vector("list", p) # when events occur
  lambda <- matrix(0, nrow = p, ncol = t_max) # intensities 
  
  intensity <- function(t) {
    lambda_temp <- mu
    for (i in 1: p) {
      for (j in 1: p) {
        if ( length(event_times[[j]]) > 0 ) {
          past_times <- event_times[[j]][ event_times[[j]] <= t ]
          lambda_temp[i] <- lambda_temp[i] + sum( alpha[i,j] * exp( -beta[i,j] * (t - past_times) ) )
        }
      }
    }
    return(lambda_temp)
  }
  
  rt <- numeric(0)
  
  while (t <= t_max) {
    lambda_old <- intensity(t)
    M <- sum(lambda_old) + max(colSums(alpha)) # upper bound 
    
    t_new <- t + rexp(1, rate = M)
    
    if (t_new > t_max) break
    
    lambda_new <- intensity(t_new)
    accept <- runif(1) <= sum(lambda_new) / M
    t <- t_new
    
    rt <- c(rt, sum(lambda_new) / M)
    
    if (accept) {
      type <- sample( 1: p, size = 1, prob = lambda_new / sum(lambda_new) )
      event_times[[type]] <- c(event_times[[type]], t)
      if (show_info) {
        print(paste("An event is accepted at time ", t, ".", sep = ""))
        print("The probability of each type is:")
        print(round(lambda_new / sum(lambda_new), 4))
        print(paste("This event is finally allocated to type ", type, ".", sep = ""))
        print("")
      }
    }
  }
  
  for (t in 1: t_max) {
    lambda[, t] <- intensity(t)
  }
  
  return(list(lambda = lambda, 
              event_times = event_times, 
              rt = rt))
}