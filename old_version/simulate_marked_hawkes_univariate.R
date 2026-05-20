simulate_marked_hawkes_univariate <- function(t_max, p, intensity) {
  # Simulate multivariate marked Hawkes process: 
  #   lambda_{i}(t) = mu_{i} + \sum_{j=1}^{p} b_{ij} \sum_{t_k < t} phi_{ij}(t - t_k)
  # Inputs: 
  
  # Outputs:
  
  if (p > 1 & !is.list(intensity)) stop("For the multivariate Hawkes process, intensity must be a list of functions.")
  if (p > 1 & any( !sapply(intensity, is.function) )) stop("For the multivariate Hawkes process, intensity must be a list of functions.")
  if (p == 1 & !is.function(intensity)) stop("Intensity must be a function.")
  
  t <- 0
  times <- c() # a vector indicating when events occur
  lambda <- c() # a vector of intensities 
  
  while (t <= t_max) {
    M <- intensity(t, times)
    t_new <- t + rexp(1, rate = M)
    
    lambda_new <- intensity(t_new, times)
    accept <- runif(1) <= lambda_new / M
    
    if (t_new > t_max) break

    t <- t_new
    
    if (accept) {
      type <- sample(1: p, size = 1, prob = lambda_new / lambda_new)
      times <- c(times, t)
    }
  }
  
  lambda <- sapply(1: t_max, intensity, times)
  
  return(list(lambda = lambda, 
              event_times = times))
}