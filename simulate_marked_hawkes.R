simulate_marked_hawkes <- function(t_max, mu, kernel, show_info = FALSE) {
  # Simulate multivariate marked Hawkes process: 
  #   lambda_{i}(t) = mu_{i} + \sum_{j=1}^{p} \sum_{t_k < t} phi_{ij}(t - t_k)
  # Inputs: 
  
  # Outputs:
  
  if (is.function(kernel)) {
    p <- 1 # dimension of Hawkes process
    kernel <- matrix(list(kernel))
  } else if (is.matrix(kernel)) {
    p <- nrow(kernel)
    if ( nrow(kernel) != ncol(kernel) ) stop("Kernel must be a square matrix.")
    if ( !is.list(kernel) ) stop("For the multivariate Hawkes process, kernel must be a matrix of functions.")  
    if ( any(!sapply(kernel, function(x) is.function(x) | is.na(x))) ) stop("For the multivariate Hawkes process, kernel must be a matrix of functions.")  
    
    if ( any(is.na(kernel)) ) kernel[ which(is.na(kernel)) ] <- list(function(t, times) 0)
  } else {
    stop("Kernel must be a function or a matrix of functions.")
  }
  
  if ( !is.numeric(mu) ) stop("mu must be a vector of numbers.")
  if ( length(mu) != p ) stop("The dimensions of mu and kernel are different.")
  
  t <- 0
  event_times <- vector("list", p) # when events occur
  lambda <- lapply(mu, rep, times = t_max) # intensities 
  
  kernel_self <- diag(kernel)
  
  while (t <= t_max) {
    M <- sum(mu) + sum( mapply(function(f, x) f(t, x), kernel_self, event_times) )
    t_new <- t + rexp(1, rate = M)
    
    if (t_new > t_max) break
    
    lambda_local <- mu + mapply(function(f, x) f(t_new, x), kernel_self, event_times)
    lambda_new <- sum(lambda_local)
    accept <- runif(1) <= lambda_new / M
    
    t <- t_new
    
    if (accept) {
      type <- sample(1: p, size = 1, prob = lambda_local / lambda_new)
      
      event_times[[type]] <- c(event_times[[type]], t)
      
      if (show_info) {
        print(paste("An event is accepted at time ", t_new, ".", sep = ""))
        print("The probability of each type is:")
        print(round(lambda_local / lambda_new), 4)
        print(paste("This event is finally allocated to type ", type, ".", sep = ""))
        print("")
      }
    }
  }
  
  for (t in 1: t_max) {
    for (j in 1: p) {
      lambda[[j]][t] <- lambda[[j]][t] + sum( mapply(function(f, x) f(t, x), kernel[j, ], event_times) )
    }
    #excitation_self <- mapply(function(f, x) f(t, x), kernel, event_times)
    #for (i in 1: length(excitation_self)) {
    #  lambda[[i]] <- c(lambda[[i]], excitation_self[i])
    #}
  }
  
  return(list(lambda = lambda, 
              event_times = event_times))
}