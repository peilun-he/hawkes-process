# ------------------------------------------------------------
# Marked Hawkes simulator (Ogata thinning) with:
#   mu(t) = exp(beta0 + beta^T c(t))  [c(t) piecewise constant]
#   kappa(m) = exp(alpha0 + alpha^T m)
#   g(tau) = sum_k w_k exp(-omega_k tau)
# ------------------------------------------------------------

# Remarks:
# Covariates: edit cov_times and cov_values to reflect a real geochem time series as a step function (or discretized time grid).

# Baseline sensitivity: beta0, beta (controls mean activity and its response to chemistry).

# Mark effect: alpha0, alpha and mark_sampler (controls triggering productivity).

# Kernel time scales: omega (large = short memory; small = long memory) and w (weights of each time scale).

simulate_marked_hawkes <- function(
    T_end,
    beta0, beta,
    cov_times,
    cov_values,
    alpha0, alpha,
    w, omega,
    mark_sampler = function(n) rnorm(n, mean = 0, sd = 1),
    marks_given = NULL,
    max_events = 1e5,
    seed = NULL,
    verbose = FALSE,
    show_progress = TRUE,
    progress_every = 200,
    eta_smooth = 0.2
) {
  
  if (!is.null(seed)) set.seed(seed)
  
  # ---- checks
  K <- length(w)
  stopifnot(length(omega) == K, all(w >= 0), all(omega > 0))
  
  p <- length(beta)
  stopifnot(is.numeric(cov_times), is.matrix(cov_values), ncol(cov_values) == p)
  stopifnot(length(cov_times) == nrow(cov_values) + 1)
  stopifnot(abs(cov_times[1] - 0) < 1e-12)
  stopifnot(abs(tail(cov_times, 1) - T_end) < 1e-12)
  
  # ---- covariates
  cov_at <- function(t) {
    if (t >= T_end) return(cov_values[nrow(cov_values), , drop = TRUE])
    r <- findInterval(t, cov_times, all.inside = TRUE)
    cov_values[r, , drop = TRUE]
  }
  
  mu_at <- function(t) exp(beta0 + sum(beta * cov_at(t)))
  
  kernel_sum_at <- function(t, times, kappas) {
    if (length(times) == 0) return(0)
    
    dt <- t - times
    dt <- dt[dt > 0]
    if (length(dt) == 0) return(0)
    
    E <- exp(-outer(dt, omega))
    as.numeric(sum(kappas[seq_along(dt)] * (E %*% w)))
  }
  
  kappa_of_m <- function(m) exp(alpha0 + alpha * m)
  
  # ---- storage
  times  <- numeric(0)
  marks  <- numeric(0)
  kappas <- numeric(0)
  
  t <- 0
  n_events <- 0
  n_steps  <- 0
  
  # ---- progress + ETA
  if (show_progress) {
    
    pb <- utils::txtProgressBar(min = 0, max = T_end, style = 3)
    
    start_time <- Sys.time()
    eta_ema <- NA_real_
    
    fmt_eta <- function(seconds) {
      if (!is.finite(seconds) || seconds < 0) return("ETA --:--:--")
      
      h <- floor(seconds / 3600)
      m <- floor((seconds - 3600*h) / 60)
      s <- floor(seconds - 3600*h - 60*m)
      
      sprintf("ETA %02d:%02d:%02d", h, m, s)
    }
    
    print_status <- function(force = FALSE) {
      
      frac <- max(min(t / T_end, 1), 1e-12)
      
      elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
      
      eta_raw <- elapsed * (1/frac - 1)
      
      if (is.na(eta_ema)) {
        eta_ema <<- eta_raw
      } else {
        eta_ema <<- (1 - eta_smooth) * eta_ema + eta_smooth * eta_raw
      }
      
      acc_rate <- if (n_steps > 0) n_events / n_steps else 0
      
      status <- sprintf(
        "  %s | t=%.2f/%.2f | N=%d | proposals=%d | acc=%.3f",
        fmt_eta(eta_ema), t, T_end, n_events, n_steps, acc_rate
      )
      
      cat("\r", status, "    ", sep = "")
      
      if (force) cat("\n")
      
      flush.console()
    }
    
    on.exit({
      try(utils::setTxtProgressBar(pb, T_end), silent = TRUE)
      try(close(pb), silent = TRUE)
      print_status(force = TRUE)
    }, add = TRUE)
    
  }
  
  # ---- simulation loop (Ogata thinning)
  while (TRUE) {
    
    mu_t   <- mu_at(t)
    self_t <- kernel_sum_at(t, times, kappas)
    
    M <- mu_t + self_t
    
    if (!is.finite(M) || M <= 0) break
    
    t_prop <- t + rexp(1, rate = M)
    
    if (t_prop > T_end) break
    
    lambda_prop <- mu_at(t_prop) + kernel_sum_at(t_prop, times, kappas)
    
    accept <- runif(1) <= lambda_prop / M
    
    t <- t_prop
    n_steps <- n_steps + 1
    
    if (accept) {
      
      n_events <- n_events + 1
      
      if (n_events > max_events) stop("Exceeded max_events")
      
      times <- c(times, t)
      
      m_i <- if (is.null(marks_given)) mark_sampler(1) else marks_given[n_events]
      
      marks <- c(marks, m_i)
      
      kappas <- c(kappas, kappa_of_m(m_i))
      
    }
    
    if (show_progress) {
      if (accept || (n_steps %% progress_every == 0)) {
        utils::setTxtProgressBar(pb, min(t, T_end))
        print_status()
      }
    }
    
  }
  
  # ---- intensity function
  lambda_fun <- function(tgrid) {
    sapply(tgrid, function(tt) mu_at(tt) + kernel_sum_at(tt, times, kappas))
  }
  
  structure(
    list(
      T_end = T_end,
      times = times,
      marks = marks,
      kappas = kappas,
      lambda = lambda_fun
    ),
    class = "marked_hawkes_sim"
  )
}

plot_marked_hawkes <- function(sim, n_grid = 2000) {
  stopifnot(inherits(sim, "marked_hawkes_sim"))
  tg <- seq(0, sim$T_end, length.out = n_grid)
  lam <- sim$lambda(tg)
  
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  
  par(mfrow = c(2,1), mar = c(4,4,2,1))
  
  plot(tg, lam, type = "l", xlab = "time", ylab = "lambda(t)",
       main = "Simulated marked Hawkes: intensity and events")
  if (length(sim$times) > 0) rug(sim$times)
  
  if (length(sim$times) > 0) {
    plot(sim$times, seq_along(sim$times), type = "s",
         xlab = "time", ylab = "N(t)", main = "Cumulative event count")
  } else {
    plot(0, 0, type = "n", xlab = "time", ylab = "N(t)",
         xlim = c(0, sim$T_end), ylim = c(0, 1),
         main = "Cumulative event count (no events)")
  }
}


T_end <- 50

cov_times  <- c(0, 15, 25, 50)
cov_values <- matrix(c(0.0, 1.0, -0.5), ncol = 1)

sim <- simulate_marked_hawkes(
  T_end = T_end,
  beta0 = -2.0, beta = 1.2,
  cov_times = cov_times, cov_values = cov_values,
  alpha0 = -0.9, alpha = 0.05,
  w = c(0.7, 0.3),
  omega = c(1.0, 0.08),
  seed = 42,
  show_progress = TRUE,
  progress_every = 300
)

cat("Simulated events:", length(sim$times), "\n")
plot_marked_hawkes(sim)


#####################################################################

#####################################################################
# Estimation via Algorithm 1 - robust, end-to-end solution
# Fixes ALL dimension/NULL issues for c_i (and marks), including n=0 or n=1 events.
# Fixes NA/NaN/Inf objective issues by:
#   - using safe_exp everywhere exp(.) appears
#   - guarding lambda/loglik/objective against non-finite values
#   - robust stopping criterion that never evaluates NA in if()
#   - step-size damping when objective becomes non-finite
#####################################################################

# Safe exp to avoid overflow/underflow
safe_exp <- function(x, max_x = 50, min_x = -50) {
  exp(pmin(pmax(x, min_x), max_x))
}

# ------------------------------------------------------------
# Utilities: soft-threshold, elastic-net proximal
# ------------------------------------------------------------
soft_thresh <- function(x, lam) sign(x) * pmax(abs(x) - lam, 0)

prox_elastic_net <- function(x, step, lam1, lam2) {
  soft_thresh(x, step * lam1) / (1 + step * lam2)
}

# ------------------------------------------------------------
# Robust covariate builder: c_i = c(t_i) from piecewise-constant covariate path
# ALWAYS returns an n x p_c matrix, even when n=0 or n=1 or p_c=1.
# ------------------------------------------------------------
build_c_i_from_piecewise <- function(t, cov_times, cov_values, T_end) {
  cov_values <- as.matrix(cov_values)
  p_c <- ncol(cov_values)
  n <- length(t)
  
  if (n == 0) return(matrix(numeric(0), nrow = 0, ncol = p_c))
  
  cov_at <- function(tt) {
    if (tt >= T_end) return(cov_values[nrow(cov_values), , drop = TRUE])
    r <- findInterval(tt, cov_times, all.inside = TRUE)
    cov_values[r, , drop = TRUE]
  }
  
  out <- t(vapply(t, cov_at, FUN.VALUE = numeric(p_c)))
  out <- as.matrix(out) # forces 2D even when n=1/p_c=1
  
  # enforce exact n x p_c shape
  if (nrow(out) != n || ncol(out) != p_c) {
    out <- matrix(as.vector(out), nrow = n, ncol = p_c, byrow = TRUE)
  }
  out
}

# ------------------------------------------------------------
# Forward pass: compute loglik and gradients using recursions
# Robust handling:
#   - If c_i is NULL or mismatched, rebuild from cov_times/cov_values
#   - If m is scalar vs matrix, handles both
#   - Uses safe_exp to prevent overflow
#   - Guards loglik/lambda against non-finite values
# ------------------------------------------------------------
forward_pass_marked_hawkes <- function(
    t, m, c_i,
    cov_times, cov_values,
    alpha0, alpha, beta0, beta,
    w, omega,
    T_end,
    eps = 1e-12
) {
  # sort events and align marks
  n <- length(t)
  
  cov_values <- as.matrix(cov_values)
  p_c <- ncol(cov_values)
  
  # handle n=0 quickly
  if (n == 0) {
    mu_int <- 0
    for (r in 1:nrow(cov_values)) {
      mu_r_lin <- beta0 + sum(beta * cov_values[r, ])
      mu_int <- mu_int + safe_exp(mu_r_lin) * (cov_times[r + 1] - cov_times[r])
    }
    return(list(
      loglik = -mu_int,
      lambda_i = numeric(0),
      kappa = numeric(0),
      mu_i = numeric(0),
      R = matrix(0, nrow = 0, ncol = length(w)),
      Ck = rep(0, length(w)),
      grad_w = rep(0, length(w)),
      grad_alpha_tilde = rep(0, 1 + length(alpha)),
      grad_beta_tilde = -c(mu_int, rep(0, p_c)),
      t_sorted = t,
      m_sorted = m,
      c_i_sorted = matrix(numeric(0), nrow = 0, ncol = p_c)
    ))
  }
  
  ord <- order(t)
  t <- t[ord]
  if (is.null(dim(m))) {
    m <- m[ord]
  } else {
    m <- m[ord, , drop = FALSE]
  }
  
  # robust c_i: rebuild if NULL or wrong rows or not a matrix
  if (is.null(c_i) || !is.matrix(c_i) || nrow(c_i) != n) {
    c_i <- build_c_i_from_piecewise(t, cov_times, cov_values, T_end)
  } else {
    c_i <- c_i[ord, , drop = FALSE]
  }
  
  # basic checks
  stopifnot(ncol(c_i) == p_c)
  stopifnot(abs(cov_times[1] - 0) < 1e-12, abs(tail(cov_times, 1) - T_end) < 1e-12)
  
  K <- length(w)
  p_m <- length(alpha)
  
  # kappa_i = exp(alpha0 + alpha^T m_i)
  if (is.null(dim(m))) {
    stopifnot(p_m == 1)
    lin_k <- alpha0 + alpha[1] * m
    z_mat <- cbind(1, m)
  } else {
    stopifnot(nrow(m) == n, ncol(m) == p_m)
    lin_k <- alpha0 + as.vector(m %*% alpha)
    z_mat <- cbind(1, m)
  }
  kappa <- safe_exp(lin_k)
  kappa[!is.finite(kappa)] <- 0
  
  # baseline mu_i = exp(beta0 + beta^T c_i)
  stopifnot(length(beta) == p_c)
  mu_lin <- beta0 + as.vector(c_i %*% beta)
  mu_i <- safe_exp(mu_lin)
  mu_i[!is.finite(mu_i)] <- 0
  
  # Recursions
  R <- matrix(0, nrow = n, ncol = K)
  p_alpha_tilde <- 1 + p_m
  A <- array(0, dim = c(n, K, p_alpha_tilde))
  
  if (n >= 2) {
    for (i in 2:n) {
      dt <- t[i] - t[i - 1]
      for (k in 1:K) {
        ed <- exp(-omega[k] * dt)
        R[i, k] <- ed * (R[i - 1, k] + kappa[i - 1])
        A[i, k, ] <- ed * (A[i - 1, k, ] + kappa[i - 1] * z_mat[i - 1, ])
      }
    }
  }
  
  lambda_i <- mu_i + as.vector(R %*% w)
  lambda_i[!is.finite(lambda_i)] <- eps
  lambda_i <- pmax(lambda_i, eps)
  
  # integral baseline mu over time (piecewise constant)
  mu_int <- 0
  for (r in 1:nrow(cov_values)) {
    mu_r_lin <- beta0 + sum(beta * cov_values[r, ])
    mu_int <- mu_int + safe_exp(mu_r_lin) * (cov_times[r + 1] - cov_times[r])
  }
  
  # Ck
  Ck <- numeric(K)
  for (k in 1:K) {
    Ck[k] <- sum(kappa * (1 - exp(-omega[k] * (T_end - t))) / omega[k])
  }
  Ck[!is.finite(Ck)] <- 0
  
  loglik <- sum(log(lambda_i)) - mu_int - sum(w * Ck)
  if (!is.finite(loglik)) loglik <- -Inf
  
  # gradients
  grad_w <- numeric(K)
  for (k in 1:K) grad_w[k] <- sum(R[, k] / lambda_i) - Ck[k]
  grad_w[!is.finite(grad_w)] <- 0
  
  # grad alpha~
  term1 <- rep(0, p_alpha_tilde)
  for (i in 1:n) {
    tmp <- rep(0, p_alpha_tilde)
    for (k in 1:K) tmp <- tmp + w[k] * A[i, k, ]
    term1 <- term1 + (1 / lambda_i[i]) * tmp
  }
  
  term2 <- rep(0, p_alpha_tilde)
  for (k in 1:K) {
    weights_jk <- (1 - exp(-omega[k] * (T_end - t))) / omega[k]
    zj_weighted_sum <- colSums((kappa * weights_jk) * z_mat)
    term2 <- term2 + w[k] * zj_weighted_sum
  }
  grad_alpha_tilde <- term1 - term2
  grad_alpha_tilde[!is.finite(grad_alpha_tilde)] <- 0
  
  # grad beta~
  c_tilde <- cbind(1, c_i)  # n x (1+p_c)
  grad_beta_tilde_events <- colSums((mu_i / lambda_i) * c_tilde)
  
  grad_beta_tilde_int <- rep(0, 1 + p_c)
  for (r in 1:nrow(cov_values)) {
    Delta <- cov_times[r + 1] - cov_times[r]
    mu_r_lin <- beta0 + sum(beta * cov_values[r, ])
    mu_r <- safe_exp(mu_r_lin)
    grad_beta_tilde_int <- grad_beta_tilde_int + mu_r * Delta * c(1, cov_values[r, ])
  }
  grad_beta_tilde <- grad_beta_tilde_events - grad_beta_tilde_int
  grad_beta_tilde[!is.finite(grad_beta_tilde)] <- 0
  
  list(
    loglik = loglik,
    lambda_i = lambda_i,
    kappa = kappa,
    mu_i = mu_i,
    R = R,
    Ck = Ck,
    grad_w = grad_w,
    grad_alpha_tilde = grad_alpha_tilde,
    grad_beta_tilde  = grad_beta_tilde,
    t_sorted = t,
    m_sorted = m,
    c_i_sorted = c_i
  )
}

# ------------------------------------------------------------
# Algorithm 1: Block-coordinate proximal estimation (robust)
# Adds:
#   - robust convergence test (finite checks)
#   - step damping when objective becomes non-finite
#   - safe_exp in kappa_est for branching ratio
# ------------------------------------------------------------
fit_algo1_prox <- function(
    t, m, c_i = NULL,
    cov_times, cov_values,
    omega,
    T_end = max(t),
    lam1 = 0.0, lam2 = 0.0,
    step_alpha = 1e-3,
    step_beta  = 1e-3,
    step_w     = 1e-2,
    max_iter = 300,
    tol = 1e-6,
    shrink_w_if_unstable = 0.9,
    init = NULL,
    verbose = TRUE,
    # recovery controls
    damp_on_nonfinite = TRUE,
    damp_factor = 0.5,
    min_step = 1e-10
) {
  cov_values <- as.matrix(cov_values)
  p_c <- ncol(cov_values)
  K <- length(omega)
  p_m <- if (is.null(dim(m))) 1 else ncol(m)
  
  if (is.null(init)) {
    alpha0 <- 0; alpha <- rep(0, p_m)
    beta0  <- 0; beta  <- rep(0, p_c)
    w      <- rep(0, K)
  } else {
    alpha0 <- init$alpha0; alpha <- init$alpha
    beta0  <- init$beta0;  beta  <- init$beta
    w      <- init$w
  }
  
  obj_prev <- NA_real_
  history <- data.frame(iter = integer(0), obj = numeric(0), loglik = numeric(0), n_branch = numeric(0))
  
  for (it in 1:max_iter) {
    fp <- forward_pass_marked_hawkes(
      t = t, m = m, c_i = c_i,
      cov_times = cov_times, cov_values = cov_values,
      alpha0 = alpha0, alpha = alpha,
      beta0 = beta0, beta = beta,
      w = w, omega = omega,
      T_end = T_end
    )
    
    # keep aligned, sorted data
    t <- fp$t_sorted
    m <- fp$m_sorted
    c_i <- fp$c_i_sorted
    
    loglik <- fp$loglik
    
    theta_pen <- c(alpha0, alpha, beta0, beta, w)
    obj <- -loglik + lam1 * sum(abs(theta_pen)) + (lam2/2) * sum(theta_pen^2)
    
    # If objective non-finite, damp and retry next iteration
    if (!is.finite(obj)) {
      if (damp_on_nonfinite) {
        step_alpha <- max(step_alpha * damp_factor, min_step)
        step_beta  <- max(step_beta  * damp_factor, min_step)
        step_w     <- max(step_w     * damp_factor, min_step)
        w <- 0.9 * w
        if (verbose) cat(sprintf("iter %d | objective non-finite -> damping steps (a=%.2e b=%.2e w=%.2e)\n",
                                 it, step_alpha, step_beta, step_w))
        next
      } else {
        stop(sprintf("Objective became non-finite at iter %d.", it))
      }
    }
    
    # w update (prox + nonnegativity)
    grad_w_obj <- -fp$grad_w + lam2 * w
    grad_w_obj[!is.finite(grad_w_obj)] <- 0
    w_tmp <- w - step_w * grad_w_obj
    w_new <- pmax(0, soft_thresh(w_tmp, step_w * lam1))
    
    # alpha~ update
    grad_alpha_obj <- -fp$grad_alpha_tilde + lam2 * c(alpha0, alpha)
    grad_alpha_obj[!is.finite(grad_alpha_obj)] <- 0
    alpha_tilde <- c(alpha0, alpha)
    alpha_tilde_new <- prox_elastic_net(alpha_tilde - step_alpha * grad_alpha_obj, step_alpha, lam1, lam2)
    alpha0_new <- alpha_tilde_new[1]
    alpha_new  <- alpha_tilde_new[-1]
    
    # beta~ update
    grad_beta_obj <- -fp$grad_beta_tilde + lam2 * c(beta0, beta)
    grad_beta_obj[!is.finite(grad_beta_obj)] <- 0
    beta_tilde <- c(beta0, beta)
    beta_tilde_new <- prox_elastic_net(beta_tilde - step_beta * grad_beta_obj, step_beta, lam1, lam2)
    beta0_new <- beta_tilde_new[1]
    beta_new  <- beta_tilde_new[-1]
    
    # branching ratio stability (empirical Ekappa)
    if (length(t) == 0) {
      n_branch <- 0
    } else {
      if (is.null(dim(m))) {
        kappa_est <- safe_exp(alpha0_new + alpha_new[1] * m)
      } else {
        kappa_est <- safe_exp(alpha0_new + as.vector(m %*% alpha_new))
      }
      Ekappa <- mean(kappa_est)
      n_branch <- Ekappa * sum(w_new / omega)
      if (!is.finite(n_branch)) n_branch <- Inf
      
      if (n_branch >= 1) {
        w_new <- w_new * shrink_w_if_unstable
        n_branch <- Ekappa * sum(w_new / omega)
        if (!is.finite(n_branch)) n_branch <- 0
      }
    }
    
    # commit
    w <- w_new
    alpha0 <- alpha0_new; alpha <- alpha_new
    beta0  <- beta0_new;  beta  <- beta_new
    
    history <- rbind(history, data.frame(iter = it, obj = obj, loglik = loglik, n_branch = n_branch))
    
    if (verbose && (it %% 25 == 0 || it == 1)) {
      cat(sprintf("iter %d | obj %.6f | loglik %.6f | n %.3f | steps(a,b,w)=(%.1e,%.1e,%.1e) | nnz(w)=%d\n",
                  it, obj, loglik, n_branch, step_alpha, step_beta, step_w, sum(w > 1e-10)))
    }
    
    # robust convergence test (never evaluates NA)
    if (is.finite(obj_prev) && is.finite(obj)) {
      rel_change <- abs(obj_prev - obj) / (abs(obj_prev) + 1e-12)
      if (is.finite(rel_change) && rel_change < tol) break
    }
    obj_prev <- obj
  }
  
  list(
    alpha0 = alpha0, alpha = alpha,
    beta0 = beta0, beta = beta,
    w = w, omega = omega,
    history = history
  )
}

#####################################################################
# Synthetic test case (robust): uses the simulator you already have
# NOTE: Requires simulate_marked_hawkes() and plot_marked_hawkes()
# from earlier in your session (the version that returns $times and $marks).
#####################################################################

# ------------------------------------------------------------
# TEST CASE
# ------------------------------------------------------------
T_end <- 60
cov_times  <- c(0, 5, 10, T_end)
cov_values <- matrix(c(0.0, 0.1, -0.5), ncol = 1)

beta0_true <- -2
beta_true  <-  1.2
alpha0_true <- -0.9
alpha_true  <-  0.05

omega_grid <- c(1.0, 0.2, 0.05)
w_true     <- c(0.40, 0.25, 0.10)

sim <- simulate_marked_hawkes(
  T_end = T_end,
  beta0 = beta0_true, beta = beta_true,
  cov_times = cov_times, cov_values = cov_values,
  alpha0 = alpha0_true, alpha = alpha_true,
  w = c(0.7, 0.3),
  omega = c(1.0, 0.08),
  seed = 42,
  show_progress = TRUE,
  progress_every = 300
)

cat("Simulated events:", length(sim$times), "\n")
plot_marked_hawkes(sim)

# If too few events, bump baseline/horizon (otherwise estimation can be unstable)
if (length(sim$times) < 5) {
  message("Too few events; re-simulating with longer horizon and higher baseline.")
  T_end <- 80
  cov_times  <- c(0, 20, 40, 80)
  cov_values <- matrix(c(0.0, 0.1, -0.5), ncol = 1)
  
  sim <- simulate_marked_hawkes(
    T_end = T_end,
    beta0 = beta0_true + 0.8, beta = beta_true,
    cov_times = cov_times, cov_values = cov_values,
    alpha0 = alpha0_true, alpha = alpha_true,
    w = c(0.7, 0.3),
    omega = c(1.0, 0.08),
    seed = 43,
    show_progress = TRUE,
    progress_every = 300
  )
  cat("Simulated events (retry):", length(sim$times), "\n")
}

t_obs <- sim$times
m_obs <- sim$marks

# Fit WITHOUT manually building c_i (the fitter will build it robustly)
fit <- fit_algo1_prox(
  t = t_obs,
  m = m_obs,
  c_i = NULL,                 # <-- key: let the estimator rebuild c_i robustly
  cov_times = cov_times,
  cov_values = cov_values,
  omega = omega_grid,
  T_end = T_end,
  lam1 = 1e-3, lam2 = 1e-3,
  step_alpha = 5e-4,
  step_beta  = 5e-4,
  step_w     = 5e-3,
  max_iter = 400,
  tol = 1e-6,
  verbose = TRUE
)

stopifnot(nrow(fit$history) >= 2)
stopifnot(tail(fit$history$obj, 1) <= head(fit$history$obj, 1) + 1e-6)

cat("\nTRUE vs EST:\n")
cat(sprintf("beta0:  true %.3f  | est %.3f\n", beta0_true, fit$beta0))
cat(sprintf("beta :  true %.3f  | est %.3f\n", beta_true,  fit$beta))
cat(sprintf("alpha0: true %.3f  | est %.3f\n", alpha0_true, fit$alpha0))
cat(sprintf("alpha : true %.3f  | est %.3f\n", alpha_true,  fit$alpha))
cat("w true:", paste(round(w_true, 3), collapse = ", "), "\n")
cat("w est :", paste(round(fit$w, 3), collapse = ", "), "\n")
cat(sprintf("branching ratio estimate (last iter): %.3f\n", tail(fit$history$n_branch, 1)))



