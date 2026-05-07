#==============================================================================
# AMAZON RAINFOREST TIPPING POINT - BAYESIAN CALIBRATION
#==============================================================================

rm(list = ls())
#setwd("")
set.seed(1234)

#==============================================================================
# LIBRARIES
#==============================================================================
libraries <- c(
  "dplyr", "tidyr", "purrr", "readxl", "writexl",
  "knitr", "kableExtra", "readr",
  "ggplot2", "deSolve", "parallel", "doParallel",
  "lhs", "patchwork", "truncnorm", "digest")
lapply(libraries, require, character.only = TRUE)

#define direction
mydirection= file.path(getwd(),"analysis_figures/")
mydirection_numericalmodel= file.path(getwd(),"numerical_model/")

#==============================================================================
# PARALLEL CLUSTER
#==============================================================================
n_cores <- parallel::detectCores() / 2
cl      <- makeCluster(n_cores)
registerDoParallel(cl)
cat("Cluster started with", n_cores, "cores\n")

#==============================================================================
# LOAD DATA
#==============================================================================
epsilons <- read.table(
  file.path(mydirection, "output/epsilons.csv"),
  sep = ";", header = FALSE, stringsAsFactors = FALSE)
epsilons    <- setNames(epsilons$V2, epsilons$V1)
epsilon_max <- epsilons[["epsilon_max"]]

# Time steps (5-year periods over 1980-2200 ~ 37 periods)
Times <- seq(1, 37, by = 1)

#==============================================================================
# EXPERT ELICITATION - Kriegler et al. (2009)
# Model each expert interval as uniform, fit Beta by method of moments
#==============================================================================
proba <- read.csv(file.path(mydirection, "data/kriegler_expertproba.csv"))

sample_expert_distribution <- function(proba_data, corridor_name, 
                                        n_samples = 10000) {
  experts  <- proba_data %>% filter(corridor == corridor_name)
  n_experts <- nrow(experts)
  if (n_experts == 0) stop(paste("No experts found for corridor:", corridor_name))
  
  samples <- numeric(n_samples)
  for (i in 1:n_samples) {
    idx        <- sample(1:n_experts, 1)
    samples[i] <- runif(1, experts$lower[idx], experts$upper[idx])
  }
  return(samples)
}

set.seed(123)
expert_samples_low    <- sample_expert_distribution(proba, "low")
expert_samples_medium <- sample_expert_distribution(proba, "medium")
expert_samples_high   <- sample_expert_distribution(proba, "high")

expert_samples_df <- data.frame(
  value    = c(expert_samples_low, expert_samples_medium, expert_samples_high),
  corridor = factor(
    rep(c("low", "medium", "high"), each = 10000),
    levels = c("low", "medium", "high")
  )
)

#==============================================================================
# FIT BETA DISTRIBUTIONS BY METHOD OF MOMENTS
#==============================================================================
fit_beta_from_moments <- function(mean_val, var_val) {
  if (mean_val <= 0 || mean_val >= 1) {
    warning(paste("Mean", mean_val, "out of (0,1), using Beta(1,1)"))
    return(list(alpha = 1, beta = 1))
  }
  if (var_val <= 0) {
    warning("Variance <= 0, using Beta(1,1)")
    return(list(alpha = 1, beta = 1))
  }
  common_term <- mean_val * (1 - mean_val) / var_val - 1
  if (common_term <= 0) {
    warning("Invalid Beta parameters (variance too large), using Beta(1,1)")
    return(list(alpha = 1, beta = 1))
  }
  alpha <- mean_val * common_term
  beta  <- (1 - mean_val) * common_term
  if (alpha <= 0 || beta <= 0 || !is.finite(alpha) || !is.finite(beta)) {
    warning(paste("Alpha=", alpha, "Beta=", beta, "invalid, using Beta(1,1)"))
    return(list(alpha = 1, beta = 1))
  }
  return(list(alpha = alpha, beta = beta))
}

expert_distributions <- data.frame(
  corridor = c("low", "medium", "high")
) %>%
  rowwise() %>%
  mutate(
    mean_empirical = mean(get(paste0("expert_samples_", corridor))),
    var_empirical  = var(get(paste0("expert_samples_", corridor))),
    beta_params    = list(fit_beta_from_moments(mean_empirical, var_empirical)),
    alpha_fitted   = beta_params$alpha,
    beta_fitted    = beta_params$beta
  ) %>%
  ungroup() %>%
  select(corridor, mean_empirical, var_empirical, alpha_fitted, beta_fitted)

cat("Fitted Beta distributions:\n")
print(expert_distributions)

#==============================================================================
# LOAD KAPPA AND TEMPERATURE TRAJECTORIES
#==============================================================================
parameters <- read.csv(file.path(mydirection_numericalmodel, "parameters/param_model.csv"), sep = ";")
colnames(parameters) <- c("variable", "value", "description")

kappa_raw <- read.csv(file.path(mydirection_numericalmodel, "parameters/param_scenariosEU.csv"),
                      header = FALSE)
kappa <- parameters$value[parameters$variable == "ratio_deg"] *
         parameters$value[parameters$variable == "deltaT"] *
         kappa_raw$V1 /
         parameters$value[parameters$variable == "total_area"] / 1000

# SSP temperature paths scaled to Kriegler corridors
temp_2200    <- read.csv(file.path(mydirection,"data/temp_ssp_ipcc2200.csv"), stringsAsFactors = FALSE)
temp         <- cbind(temp_2200[, 6], temp_2200[, 7], temp_2200[, 8])
temp_initial <- temp

# Scale temperature deviations to match Kriegler corridor endpoints
# and convert to regional Amazon warming (regional/global TCRE = 2/1.73)
temp[, 1] <- (temp[, 1] - temp[1, 1]) * 1.5 / (temp[nrow(temp), 1] - temp[1, 1]) * (2 / 1.73)
temp[, 2] <- (temp[, 2] - temp[1, 2]) * 3.0 / (temp[nrow(temp), 2] - temp[1, 2]) * (2 / 1.73)
temp[, 3] <- (temp[, 3] - temp[1, 3]) * 6.0 / (temp[nrow(temp), 3] - temp[1, 3]) * (2 / 1.73)

traj_to_corridor <- c("1" = "low", "2" = "medium", "3" = "high")

# Helper: convert (mu, phi) Beta parametrization to (alpha, beta)
mu_phi_to_alpha_beta <- function(mu_eps, phi_eps) {
  list(alpha = mu_eps * phi_eps,
       beta  = (1 - mu_eps) * phi_eps)
}

#==============================================================================
# PRECOMPUTE TEMPERATURE INTERPOLATION FUNCTIONS
# Each corridor has n_temp_points temperature trajectories uniformly spaced
# over the corridor bounds, drawn uniformly as in Kriegler (2009)
#==============================================================================
corridor_bounds <- list(
  low    = c(1.0, 2.0),
  medium = c(2.0, 4.0),
  high   = c(4.0, 8.0)
)

precompute_globals <- function(n_temp_points = 5) {
  temp_interp <<- list()
  for (corridor_name in c("low", "medium", "high")) {
    traj_idx <- which(c("low", "medium", "high") == corridor_name)
    T_min    <- corridor_bounds[[corridor_name]][1]
    T_max    <- corridor_bounds[[corridor_name]][2]
    T_finals <- seq(T_min, T_max, length.out = n_temp_points)
    
    interp_list <- lapply(T_finals, function(T_final) {
      temp_base <- temp_initial[, traj_idx]
      temp_dev  <- temp_base - temp_base[1]
      if (max(abs(temp_dev)) > 0) {
        temp_scaled <- temp_dev * T_final / temp_dev[nrow(temp_initial)] * (2 / 1.73)
      } else {
        temp_scaled <- rep(0, nrow(temp_initial))
      }
      approxfun(1:nrow(temp_initial), temp_scaled, rule = 2)
    })
    temp_interp[[corridor_name]] <<- interp_list
  }
  kappa_interp <<- approxfun(1:length(kappa), kappa, rule = 2)
  cat(sprintf(
    "Precomputed %d temperature trajectories per corridor ",
    n_temp_points
  ))
  cat(sprintf(
    "([low:%.1f-%.1fC] [med:%.1f-%.1fC] [high:%.1f-%.1fC])\n",
    corridor_bounds$low[1],    corridor_bounds$low[2],
    corridor_bounds$medium[1], corridor_bounds$medium[2],
    corridor_bounds$high[1],   corridor_bounds$high[2]
  ))
}

precompute_globals(n_temp_points = 10)

#==============================================================================
# FOREST ODE - compiled for speed
#==============================================================================
ab <- function(time, state, parms) {
  with(as.list(c(state, parms)), {
    x        <- min(max(x, 0), 1)
    time_idx <- max(1, min(floor(time), length(epsilon)))
    eps_t    <- epsilon[time_idx]
    traj_name <- c("low", "medium", "high")[traj]
    temp_val  <- temp_interp[[traj_name]][[temp_idx]](time)
    kappa_val <- kappa_interp(time)
    check <- g0 * (1 - drop * ((feedback * (1 - x)) / beta0)^eta) * x * (1 - x) -
             (eps_t * log(temp_val + 1)) * x - x * kappa_val
    return(list(c(check)))
  })
}
ab <- compiler::cmpfun(ab, options = list(optimize = 3))

#==============================================================================
# NUMERICAL INTEGRATION - one trajectory
#==============================================================================
process <- function(epsilon_chain, fun, growth, Upsilon, beta0, eta, drop,
                    traj, temp_idx = 1) {
  yini <- c(x = 1)
  pars <- list(
    epsilon  = epsilon_chain, g0 = growth, feedback = Upsilon,
    beta0    = beta0,         eta = eta,   drop     = drop,
    traj     = traj,          temp_idx = temp_idx
  )
  deSolve::lsoda(y = yini, times = Times, func = fun, parms = pars,
                 rtol = 1e-4, atol = 1e-6)
}
process <- compiler::cmpfun(process, options = list(optimize = 3))

#==============================================================================
# PARALLEL MONTE CARLO SIMULATION
# Averages over n_temp_points temperature trajectories within each corridor
# Uses stratified sampling to reduce Monte Carlo variance
#==============================================================================
simulate_trajectories_fast <- function(mu_eps, phi_eps, g0, calibration,
                                       Upsilon, beta0, eta, iteration,
                                       n_traj, traj, drop,
                                       use_stratified = TRUE, cl,
                                       n_temp_points  = 10) {
  clusterSetRNGStream(cl, iseed = sample.int(.Machine$integer.max, 1))
  ab_alpha <- mu_eps * phi_eps
  ab_beta  <- (1 - mu_eps) * phi_eps
  clusterExport(
    cl,
    varlist = c("ab_alpha", "ab_beta", "epsilon_max", "Times",
                "g0", "Upsilon", "beta0", "eta", "drop", "traj",
                "process", "ab", "n_temp_points"),
    envir = environment()
  )
  sim_one <- function(i) {
    set.seed(as.integer(Sys.time()) + i * 1000 + traj * 10000)
    if (use_stratified) {
      n_strata   <- 10
      strata_idx <- ((i - 1) %% n_strata) + 1
      u_base     <- runif(length(Times),
                          (strata_idx - 1) / n_strata,
                          strata_idx / n_strata)
      eps <- epsilon_max * qbeta(u_base, ab_alpha, ab_beta)
    } else {
      eps <- epsilon_max * rbeta(length(Times), ab_alpha, ab_beta)
    }
    # Average tipping probability over temperature corridor
    results_temp <- sapply(1:n_temp_points, function(t_idx) {
      sim <- tryCatch(
        process(eps, ab, g0, Upsilon, beta0, eta, drop, traj, temp_idx = t_idx),
        error = function(e) NULL
      )
      if (!is.null(sim)) {
        x_trimmed <- sim[-c(1, 2), "x"]
        return(max(as.numeric(x_trimmed < 0.5)))
      }
      return(NA)
    })
    valid_temp <- !is.na(results_temp)
    if (sum(valid_temp) == 0) return(NA)
    return(mean(results_temp[valid_temp]))
  }
  results <- parallel::parSapply(cl, 1:n_traj, sim_one)
  valid   <- !is.na(results)
  return(mean(results[valid]))
}

#==============================================================================
# LIKELIHOOD CACHE
#==============================================================================
likelihood_cache <- new.env()

clear_likelihood_cache <- function() {
  rm(list = ls(envir = likelihood_cache), envir = likelihood_cache)
  cat("Likelihood cache cleared\n")
}

#==============================================================================
# COMPUTE LIKELIHOOD
# shape parameter controls width of credible interval:
#   shape=0.50 -> [q25:q75], shape=0.30 -> [q35:q65]
#==============================================================================
compute_likelihood_progressive <- function(params, proba_mean, n_traj = 1000,
                                           drop = 1, verbose = FALSE,
                                           use_cache = TRUE, cl,
                                           shape = 0.5,
                                           penalty_strength = 0) {
  if (use_cache) {
    param_hash <- digest::digest(list(params, shape, penalty_strength),
                                 algo = "xxhash64")
    if (exists(param_hash, envir = likelihood_cache)) {
      if (verbose) cat("[CACHE HIT]\n")
      return(get(param_hash, envir = likelihood_cache))
    }
  }
  
  traj_results <- tryCatch({
    sapply(1:3, function(traj) {
      simulate_trajectories_fast(
        params$mu_eps, params$phi_eps, params$g0, 1.0,
        params$Upsilon, params$beta0, params$eta,
        1, n_traj, traj, drop, FALSE, cl
      )
    })
  }, error = function(e) rep(NA, 3))
  
  if (any(is.na(traj_results)) || any(!is.finite(traj_results))) {
    result <- list(log_lik = -Inf, mae = Inf, mae_raw = Inf,
                   diff1 = NA, diff2 = NA, diff3 = NA,
                   traj_low = NA, traj_medium = NA, traj_high = NA,
                   diffs_normalized = rep(NA, 3), n_in_bounds = 0, shape = shape)
    if (use_cache) assign(param_hash, result, envir = likelihood_cache)
    return(result)
  }
  
  lower_q <- (1 - shape) / 2
  upper_q <- 1 - lower_q
  
  diffs            <- numeric(3)
  diffs_normalized <- numeric(3)
  log_lik          <- 0
  n_in_bounds      <- 0
  
  if (verbose) cat(sprintf("\nCI [q%.0f:q%.0f]:\n", lower_q * 100, upper_q * 100))
  
  for (j in 1:3) {
    corridor_name <- traj_to_corridor[as.character(j)]
    row_idx       <- which(proba_mean$corridor == corridor_name)
    alpha_exp     <- proba_mean$alpha_fitted[row_idx]
    beta_exp      <- proba_mean$beta_fitted[row_idx]
    mean_exp      <- proba_mean$mean_empirical[row_idx]
    sd_exp        <- sqrt(alpha_exp * beta_exp /
                          ((alpha_exp + beta_exp)^2 * (alpha_exp + beta_exp + 1)))
    result_j      <- max(0.001, min(0.999, traj_results[j]))
    
    log_lik_corridor <- dbeta(result_j, alpha_exp, beta_exp, log = TRUE)
    log_lik          <- log_lik + log_lik_corridor
    
    diff                <- result_j - mean_exp
    diffs[j]            <- diff
    diffs_normalized[j] <- diff / max(sd_exp, 1e-8)
    
    ic_lower <- qbeta(lower_q, alpha_exp, beta_exp)
    ic_upper <- qbeta(upper_q, alpha_exp, beta_exp)
    in_ic    <- (result_j >= ic_lower) && (result_j <= ic_upper)
    if (in_ic) n_in_bounds <- n_in_bounds + 1
    
    if (verbose) {
      cat(sprintf(
        "  %s: sim=%.3f, CI=[%.3f,%.3f] %s, diff=%+.3f (%.2f sd), loglik=%.2f\n",
        corridor_name, result_j, ic_lower, ic_upper,
        ifelse(in_ic, "OK", "X"), diff, diffs_normalized[j], log_lik_corridor
      ))
    }
  }
  
  mae_raw <- mean(abs(diffs))
  mae     <- mean(abs(diffs_normalized))
  
  if (verbose) {
    cat(sprintf("loglik=%.3f | MAE_raw=%.4f | MAE_norm=%.4f | in_CI: %d/3\n",
                log_lik, mae_raw, mae, n_in_bounds))
  }
  
  result <- list(
    log_lik          = log_lik,
    mae              = mae,
    mae_raw          = mae_raw,
    diff1            = diffs[1],
    diff2            = diffs[2],
    diff3            = diffs[3],
    traj_low         = traj_results[1],
    traj_medium      = traj_results[2],
    traj_high        = traj_results[3],
    diffs_normalized = diffs_normalized,
    n_in_bounds      = n_in_bounds,
    shape            = shape
  )
  if (use_cache) assign(param_hash, result, envir = likelihood_cache)
  return(result)
}

#==============================================================================
# BISECTION: find tightest credible interval (smallest shape) with 3/3 in bounds
#==============================================================================
find_best_shape <- function(params, proba_mean, n_traj, drop, cl,
                            shape_min = 0.10, shape_max = 0.90, tol = 0.05) {
  ll_hi <- compute_likelihood_progressive(
    params, proba_mean, n_traj, drop,
    shape = shape_max, verbose = FALSE, use_cache = FALSE, cl = cl
  )
  if (ll_hi$n_in_bounds < 3) {
    return(list(best_shape = NA, n_in_bounds = ll_hi$n_in_bounds,
                mae_raw = ll_hi$mae_raw, ll = ll_hi))
  }
  lo <- shape_min
  hi <- shape_max
  while ((hi - lo) > tol) {
    mid    <- (lo + hi) / 2
    ll_mid <- compute_likelihood_progressive(
      params, proba_mean, n_traj, drop,
      shape = mid, verbose = FALSE, use_cache = FALSE, cl = cl
    )
    if (ll_mid$n_in_bounds == 3) hi <- mid else lo <- mid
  }
  final_ll <- compute_likelihood_progressive(
    params, proba_mean, n_traj, drop,
    shape = hi, verbose = FALSE, use_cache = FALSE, cl = cl
  )
  return(list(best_shape = hi, n_in_bounds = final_ll$n_in_bounds,
              mae_raw = final_ll$mae_raw, ll = final_ll))
}

#==============================================================================
# PRIOR PARAMETERS
#==============================================================================
prior_params <- list(
  mu_eps  = list(lower = 0.05, upper = 0.95),
  phi_eps = list(lower = 0.01, upper = 1.00),
  g0      = list(mu = 0.0175, sigma = 0.008, lower = 0.005, upper = 0.05),
  eta     = list(lower = 1.5, upper = 4.0, type = "uniform"),
  beta0   = list(mu = 6.0, sigma = 2/3, lower = 4.0, upper = 8.0),
  Upsilon = list(mu = 5.52, sigma = 0.51, lower = 3.98, upper = 7.06)
)

log_prior <- function(params) {
  lp <- 0
  # mu_eps: uniform
  if (params$mu_eps  < prior_params$mu_eps$lower  ||
      params$mu_eps  > prior_params$mu_eps$upper)  return(-Inf)
  lp <- lp + log(1 / (prior_params$mu_eps$upper - prior_params$mu_eps$lower))
  # phi_eps: uniform
  if (params$phi_eps < prior_params$phi_eps$lower ||
      params$phi_eps > prior_params$phi_eps$upper) return(-Inf)
  lp <- lp + log(1 / (prior_params$phi_eps$upper - prior_params$phi_eps$lower))
  # g0: truncated normal
  if (params$g0 < prior_params$g0$lower || params$g0 > prior_params$g0$upper) return(-Inf)
  z_g0     <- (params$g0 - prior_params$g0$mu) / prior_params$g0$sigma
  log_Z_g0 <- log(pnorm((prior_params$g0$upper - prior_params$g0$mu) / prior_params$g0$sigma) -
                  pnorm((prior_params$g0$lower - prior_params$g0$mu) / prior_params$g0$sigma))
  lp <- lp + dnorm(z_g0, log = TRUE) - log(prior_params$g0$sigma) - log_Z_g0
  # eta: uniform
  if (params$eta < prior_params$eta$lower || params$eta > prior_params$eta$upper) return(-Inf)
  lp <- lp + log(1 / (prior_params$eta$upper - prior_params$eta$lower))
  # beta0: truncated normal
  if (params$beta0 < prior_params$beta0$lower || params$beta0 > prior_params$beta0$upper) return(-Inf)
  z_b     <- (params$beta0 - prior_params$beta0$mu) / prior_params$beta0$sigma
  log_Z_b <- log(pnorm((prior_params$beta0$upper - prior_params$beta0$mu) / prior_params$beta0$sigma) -
                 pnorm((prior_params$beta0$lower - prior_params$beta0$mu) / prior_params$beta0$sigma))
  lp <- lp + dnorm(z_b, log = TRUE) - log(prior_params$beta0$sigma) - log_Z_b
  # Upsilon: truncated normal
  if (params$Upsilon < prior_params$Upsilon$lower || params$Upsilon > prior_params$Upsilon$upper) return(-Inf)
  z_u     <- (params$Upsilon - prior_params$Upsilon$mu) / prior_params$Upsilon$sigma
  log_Z_u <- log(pnorm((prior_params$Upsilon$upper - prior_params$Upsilon$mu) / prior_params$Upsilon$sigma) -
                 pnorm((prior_params$Upsilon$lower - prior_params$Upsilon$mu) / prior_params$Upsilon$sigma))
  lp <- lp + dnorm(z_u, log = TRUE) - log(prior_params$Upsilon$sigma) - log_Z_u
  return(lp)
}

sample_from_prior <- function(n = 1) {
  data.frame(
    mu_eps  = runif(n, prior_params$mu_eps$lower,  prior_params$mu_eps$upper),
    phi_eps = runif(n, prior_params$phi_eps$lower, prior_params$phi_eps$upper),
    g0      = rtruncnorm(n, a = prior_params$g0$lower, b = prior_params$g0$upper,
                         mean = prior_params$g0$mu, sd = prior_params$g0$sigma),
    eta     = runif(n, prior_params$eta$lower, prior_params$eta$upper),
    beta0   = rtruncnorm(n, a = prior_params$beta0$lower, b = prior_params$beta0$upper,
                         mean = prior_params$beta0$mu, sd = prior_params$beta0$sigma),
    Upsilon = rtruncnorm(n, a = prior_params$Upsilon$lower, b = prior_params$Upsilon$upper,
                         mean = prior_params$Upsilon$mu, sd = prior_params$Upsilon$sigma)
  )
}

#==============================================================================
# ADAPTIVE LHS EXPLORATION 
# Generates n_per_step candidates via improved LHS, evaluates each,
# applies bisection on 3/3 candidates, shrinks search bounds adaptively
#==============================================================================
adaptive_lhs_exploration <- function(proba_mean, n_traj, drop, cl,
                                     n_per_step    = 5000,
                                     top_k         = 10,
                                     n_steps       = 1,
                                     shrink_factor = 0.5,
                                     shape_min     = 0.10,
                                     shape_max     = 0.90,
                                     tol_bisect    = 0.05) {
  all_candidates <- NULL
  current_bounds <- NULL
  current_shape  <- shape_max

  for (step in 1:n_steps) {
    cat(sprintf("\n--- Step %d/%d (reference shape=%.2f) ---\n",
                step, n_steps, current_shape))

    if (is.null(current_bounds)) {
      bounds <- data.frame(
        lower = c(prior_params$mu_eps$lower, prior_params$phi_eps$lower,
                  prior_params$g0$lower,     prior_params$eta$lower,
                  prior_params$beta0$lower,  prior_params$Upsilon$lower),
        upper = c(prior_params$mu_eps$upper, prior_params$phi_eps$upper,
                  prior_params$g0$upper,     prior_params$eta$upper,
                  prior_params$beta0$upper,  prior_params$Upsilon$upper)
      )
    } else {
      bounds <- current_bounds
    }

    lhs_sample      <- improvedLHS(n_per_step, 6)
    candidates_step <- data.frame(
      mu_eps  = qunif(lhs_sample[, 1], bounds$lower[1], bounds$upper[1]),
      phi_eps = qunif(lhs_sample[, 2], bounds$lower[2], bounds$upper[2]),
      g0      = qtruncnorm(lhs_sample[, 3],
                           a    = prior_params$g0$lower,
                           b    = prior_params$g0$upper,
                           mean = prior_params$g0$mu,
                           sd   = prior_params$g0$sigma),
      eta     = qunif(lhs_sample[, 4], bounds$lower[4], bounds$upper[4]),
      beta0   = qunif(lhs_sample[, 5], bounds$lower[5], bounds$upper[5]),
      Upsilon = qunif(lhs_sample[, 6], bounds$lower[6], bounds$upper[6])
    )

    best_shapes        <- rep(1.0, n_per_step)
    n_in_bounds_scores <- numeric(n_per_step)
    mae_scores         <- numeric(n_per_step)
    candidates_step$traj_low    <- NA_real_
    candidates_step$traj_medium <- NA_real_
    candidates_step$traj_high   <- NA_real_
    candidates_step$diff_low    <- NA_real_
    candidates_step$diff_medium <- NA_real_
    candidates_step$diff_high   <- NA_real_

    log_file <- file.path(mydirection, "output/candidates_log.csv")


write.table(
  data.frame(step=integer(), i=integer(),
             mu_eps=numeric(), phi_eps=numeric(), g0=numeric(),
             eta=numeric(), beta0=numeric(), Upsilon=numeric(),
             diff_low=numeric(), diff_medium=numeric(), diff_high=numeric(),
             traj_low=numeric(), traj_medium=numeric(), traj_high=numeric(),
             mae_raw=numeric()),
  log_file, sep=",", row.names=FALSE, col.names=TRUE
)

    for (i in 1:n_per_step) {
      # Quick pre-filter: discard hopeless candidates cheaply
      pf_low <- tryCatch(
        simulate_trajectories_fast(
          candidates_step$mu_eps[i], candidates_step$phi_eps[i],
          candidates_step$g0[i], 1.0,
          candidates_step$Upsilon[i], candidates_step$beta0[i],
          candidates_step$eta[i], 1, 500, 1, drop, FALSE, cl
        ), error = function(e) NA)
      if (is.na(pf_low) || pf_low < 0.05) {
        n_in_bounds_scores[i] <- 0; mae_scores[i] <- 99; next
      }
      pf_high <- tryCatch(
        simulate_trajectories_fast(
          candidates_step$mu_eps[i], candidates_step$phi_eps[i],
          candidates_step$g0[i], 1.0,
          candidates_step$Upsilon[i], candidates_step$beta0[i],
          candidates_step$eta[i], 1, 500, 3, drop, FALSE, cl
        ), error = function(e) NA)
      if (is.na(pf_high) || pf_high > 0.9) {
        n_in_bounds_scores[i] <- 0; mae_scores[i] <- 99; next
      }

      # Full likelihood evaluation
      ll_quick <- compute_likelihood_progressive(
        candidates_step[i, ], proba_mean,
        n_traj = n_traj, drop = drop, shape = current_shape,
        verbose = FALSE, use_cache = FALSE, cl = cl
      )
      n_in_bounds_scores[i]          <- ll_quick$n_in_bounds
      mae_scores[i]                  <- ll_quick$mae_raw
      candidates_step$traj_low[i]    <- ll_quick$traj_low
      candidates_step$traj_medium[i] <- ll_quick$traj_medium
      candidates_step$traj_high[i]   <- ll_quick$traj_high
      candidates_step$diff_low[i]    <- ll_quick$diff1
      candidates_step$diff_medium[i] <- ll_quick$diff2
      candidates_step$diff_high[i]   <- ll_quick$diff3

      # Bisection only on 3/3 candidates
      if (ll_quick$n_in_bounds == 3) {
        log_row <- data.frame(
          step = step, i = i,
          mu_eps = candidates_step$mu_eps[i], phi_eps = candidates_step$phi_eps[i],
          g0 = candidates_step$g0[i], eta = candidates_step$eta[i],
          beta0 = candidates_step$beta0[i], Upsilon = candidates_step$Upsilon[i],
          diff_low = ll_quick$diff1, diff_medium = ll_quick$diff2,
          diff_high = ll_quick$diff3, traj_low = ll_quick$traj_low,
          traj_medium = ll_quick$traj_medium, traj_high = ll_quick$traj_high,
          mae_raw = ll_quick$mae_raw
        )
        write.table(log_row, log_file, append = TRUE, sep = ",",
                    row.names = FALSE, col.names = !file.exists(log_file))

        res            <- find_best_shape(
          candidates_step[i, ], proba_mean,
          n_traj = n_traj, drop = drop, cl = cl,
          shape_min = shape_min, shape_max = current_shape, tol = tol_bisect
        )
        best_shapes[i] <- ifelse(is.na(res$best_shape), current_shape, res$best_shape)
      }

      if (i %% 50 == 0) {
        valid <- best_shapes[1:i] < 1.0
        cat(sprintf("  %d/%d | 3/3: %d | best shape=%.3f | MAE=%.4f\n",
                    i, n_per_step, sum(valid),
                    ifelse(any(valid), min(best_shapes[1:i][valid]), NA),
                    min(mae_scores[1:i], na.rm = TRUE)))
        if (any(valid)) {
          best_i <- which(best_shapes[1:i] == min(best_shapes[1:i][valid]))[1]
          cat(sprintf(
            "    >> Best candidate (#%d): mu_eps=%.4f phi_eps=%.4f g0=%.5f eta=%.3f beta0=%.3f Upsilon=%.3f\n",
            best_i,
            candidates_step$mu_eps[best_i], candidates_step$phi_eps[best_i],
            candidates_step$g0[best_i],     candidates_step$eta[best_i],
            candidates_step$beta0[best_i],  candidates_step$Upsilon[best_i]
          ))
          cat(sprintf(
            "       diff: low=%+.4f med=%+.4f high=%+.4f | traj: low=%.3f(t:0.244) med=%.3f(t:0.437) high=%.3f(t:0.616)\n",
            candidates_step$diff_low[best_i], candidates_step$diff_medium[best_i],
            candidates_step$diff_high[best_i],
            candidates_step$traj_low[best_i], candidates_step$traj_medium[best_i],
            candidates_step$traj_high[best_i]
          ))
        }
      }
    } # end loop over candidates

    candidates_step$best_shape  <- best_shapes
    candidates_step$n_in_bounds <- n_in_bounds_scores
    candidates_step$mae         <- mae_scores
    candidates_step$score       <- -best_shapes
    all_candidates              <- rbind(all_candidates, candidates_step)

    cat(sprintf("  Step %d | 3/3: %d/%d | best shape=%.3f\n",
                step, sum(n_in_bounds_scores == 3), n_per_step,
                min(best_shapes[best_shapes < 1.0], na.rm = TRUE)))

    # Update reference shape
    best_shape_step <- min(best_shapes[best_shapes < 1.0], na.rm = TRUE)
    if (is.finite(best_shape_step)) {
      current_shape <- max(best_shape_step, shape_min)
      cat(sprintf("  New reference shape: %.3f\n", current_shape))
    }

    # Shrink bounds around top_k for next step
    if (step < n_steps) {
      top_idx     <- order(candidates_step$score, decreasing = TRUE)[1:min(top_k, n_per_step)]
      top_rows    <- candidates_step[top_idx, 1:6]
      lower_prior <- c(prior_params$mu_eps$lower, prior_params$phi_eps$lower,
                       prior_params$g0$lower,     prior_params$eta$lower,
                       prior_params$beta0$lower,  prior_params$Upsilon$lower)
      upper_prior <- c(prior_params$mu_eps$upper, prior_params$phi_eps$upper,
                       prior_params$g0$upper,     prior_params$eta$upper,
                       prior_params$beta0$upper,  prior_params$Upsilon$upper)
      current_range <- bounds$upper - bounds$lower
      new_lower     <- pmax(apply(top_rows, 2, min) - shrink_factor * current_range, lower_prior)
      new_upper     <- pmin(apply(top_rows, 2, max) + shrink_factor * current_range, upper_prior)
      current_bounds <- data.frame(lower = new_lower, upper = new_upper)
      cat("  New bounds:\n")
      print(data.frame(
        param = c("mu_eps", "phi_eps", "g0", "eta", "beta0", "Upsilon"),
        lower = round(new_lower, 4), upper = round(new_upper, 4)
      ))
    }
  } # end loop over steps

  all_candidates <- all_candidates[order(all_candidates$score, decreasing = TRUE), ]
  cat(sprintf("\nTotal explored: %d | 3/3: %d | best shape: %.3f\n",
              nrow(all_candidates),
              sum(all_candidates$n_in_bounds == 3),
              min(all_candidates$best_shape[all_candidates$best_shape < 1.0], na.rm = TRUE)))

  # Save good candidates
  good_step1 <- all_candidates[all_candidates$n_in_bounds == 3 &
                                 all_candidates$best_shape < 1.0, ]
  good_step1 <- good_step1[order(good_step1$best_shape, good_step1$mae), ]
  cat(sprintf("\n=== GOOD CANDIDATES (n_in_bounds==3): %d total ===\n", nrow(good_step1)))
  print(good_step1[, c("mu_eps", "phi_eps", "g0", "eta", "beta0", "Upsilon",
                        "best_shape", "mae", "traj_low", "traj_medium", "traj_high")])
  write.csv(good_step1,
            file.path(mydirection, "output/good_candidates_step1.csv"),
            row.names = TRUE)

  return(all_candidates)
}

#==============================================================================
# MAIN OPTIMIZATION FUNCTION 
#==============================================================================
optimize_MAP_smart <- function(proba_mean,
                               n_traj_initial = 3000,
                               n_traj_final   = 10000,
                               drop           = 1,
                               shape_min      = 0.10,
                               shape_max      = 0.90,
                               tol_bisect     = 0.05,
                               n_per_step     = 5000,
                               top_k          = 10,
                               n_steps        = 1,
                               shrink_factor  = 0.5,
                               candidates_file = NULL,
                               cl) {

  param_names <- c("mu_eps", "phi_eps", "g0", "eta", "beta0", "Upsilon")

  cat(rep("=", 70), "\n", sep = "")
  cat("BAYESIAN CALIBRATION - PHASE 1 LHS EXPLORATION\n")
  cat(rep("=", 70), "\n\n", sep = "")

  # Load existing candidates or run LHS
  if (!is.null(candidates_file) && file.exists(candidates_file)) {
    cat("Loading candidates from file:", candidates_file, "\n")
    all_candidates       <- read.csv(candidates_file, row.names = 1)
    all_candidates$score <- -all_candidates$best_shape
    cat(sprintf("  Loaded %d candidates | 3/3: %d | best shape: %.3f\n",
                nrow(all_candidates),
                sum(all_candidates$n_in_bounds == 3),
                min(all_candidates$best_shape[all_candidates$best_shape < 1.0],
                    na.rm = TRUE)))
  } else {
    cat("Running LHS exploration\n")
    all_candidates <- adaptive_lhs_exploration(
      proba_mean    = proba_mean,
      n_traj        = n_traj_initial,
      drop          = drop,
      cl            = cl,
      n_per_step    = n_per_step,
      top_k         = top_k,
      n_steps       = n_steps,
      shrink_factor = shrink_factor,
      shape_max     = shape_max,
      shape_min     = shape_min,
      tol_bisect    = tol_bisect
    )
  }

  all_candidates <- all_candidates[order(all_candidates$score, decreasing = TRUE), ]

  # Print top 5
  best_idx <- order(all_candidates$score, decreasing = TRUE)[1:min(5, nrow(all_candidates))]
  cat("\nTop 5 candidates:\n")
  print(all_candidates[best_idx,
                        c("mu_eps", "phi_eps", "g0", "eta", "beta0", "Upsilon",
                          "best_shape", "mae", "traj_low", "traj_medium", "traj_high")])

  # Final evaluation of best candidate
  best_row    <- all_candidates[1, param_names]
  best_params <- as.list(best_row)

  cat(rep("=", 70), "\n", sep = "")
  cat("FINAL EVALUATION\n")
  cat(rep("=", 70), "\n", sep = "")

  final_ll <- compute_likelihood_progressive(
    best_params, proba_mean,
    n_traj = n_traj_final, drop = drop,
    shape = 0.50, verbose = TRUE, use_cache = FALSE, cl = cl
  )

  map_estimate <- data.frame(
    mu_eps            = best_row$mu_eps,
    phi_eps           = best_row$phi_eps,
    g0                = best_row$g0,
    eta               = best_row$eta,
    beta0             = best_row$beta0,
    Upsilon           = best_row$Upsilon,
    neg_log_posterior = -final_ll$log_lik,
    alpha_implied     = best_row$mu_eps * best_row$phi_eps,
    beta_implied      = (1 - best_row$mu_eps) * best_row$phi_eps
  )
  print(map_estimate)

  return(list(
    estimate            = map_estimate,
    final_ll            = final_ll,
    final_mae           = final_ll$mae,
    final_mae_raw       = final_ll$mae_raw,
    final_n_in_bounds   = final_ll$n_in_bounds,
    best_shape          = final_ll$shape,
    exploration_results = list(candidates = all_candidates)
  ))
}

#==============================================================================
# EXPORT TO CLUSTER
#==============================================================================
clusterExport(cl, varlist = c(
  "epsilon_max", "Times", "temp_initial", "kappa",
  "process", "ab", "temp_interp", "kappa_interp",
  "corridor_bounds", "mu_phi_to_alpha_beta", "traj_to_corridor",
  "simulate_trajectories_fast"
), envir = environment())

#==============================================================================
# RUN CALIBRATION
#==============================================================================
cat("\n")
cat(rep("#", 70), "\n", sep = "")
cat("BAYESIAN CALIBRATION\n")
cat("n_cores =", n_cores, "\n")
cat(rep("#", 70), "\n", sep = "")

clear_likelihood_cache()

# Dans le nouveau code propre
map_result <- optimize_MAP_smart(
  proba_mean      = expert_distributions,
  n_traj_initial  = 3000,
  n_traj_final    = 10000,
  drop            = 1,
  shape_min       = 0.10,
  shape_max       = 0.50,  
  tol_bisect      = 0.05,
  n_per_step      = 1000,
  top_k           = 10,
  n_steps         = 1,
  shrink_factor   = 0.5,
  candidates_file = NULL,
  cl              = cl
)


cands_log <- read.csv(file.path(mydirection, "output/candidates_log.csv"),
                      header = FALSE)
names(cands_log) <- c("step", "i", "mu_eps", "phi_eps", "g0", "eta",
                       "beta0", "Upsilon", "diff_low", "diff_medium",
                       "diff_high", "traj_low", "traj_medium", "traj_high",
                       "mae_raw")
cands_log$mae         <- cands_log$mae_raw
cands_log$n_in_bounds <- 3
cands_log$best_shape  <- 0.30
cands_log$traj_low    <- cands_log$traj_low
cands_log$traj_medium <- cands_log$traj_medium
cands_log$traj_high   <- cands_log$traj_high


good_candidates = cands_log

#==============================================================================
# BUILD THREE PARAMETER SETS FOR WELFARE ANALYSIS
#
# (1) Baseline: best MAE among candidates with tightest CI (shape=0.30)
# (2) Robust-low: minimizes error on low corridor (most policy-relevant)
# (3) Robust-distant: maximally distant from baseline in normalized param space
#     among candidates with shape <= 0.50 and 3/3 in bounds
#==============================================================================
#good_candidates <- read.csv(
#  file.path(mydirection, "output/good_candidates_step1.csv"),
#  row.names = 1
#)

# (1) Baseline: tightest CI, then best MAE
cand_best_shape <- good_candidates[
  good_candidates$best_shape == min(good_candidates$best_shape), ]
baseline <- cand_best_shape[which.min(cand_best_shape$mae), ]

# (2) Robust-low: minimizes |diff_low| among all 3/3 candidates
good_candidates$abs_diff_low <- abs(good_candidates$diff_low)
robust_low <- good_candidates[which.min(good_candidates$abs_diff_low), ]

# (3) Robust-distant: furthest from baseline in normalized parameter space
#     restricted to shape <= 0.50 to stay within credible region
cand_valid   <- good_candidates[good_candidates$best_shape <= 0.50, ]
scaled       <- scale(cand_valid[, c("mu_eps", "phi_eps", "g0",
                                     "eta", "beta0", "Upsilon")])
baseline_v   <- scaled[rownames(cand_valid) == rownames(baseline), , drop = FALSE]
if (nrow(baseline_v) == 0) baseline_v <- matrix(colMeans(scaled), nrow = 1)
dists        <- apply(scaled, 1, function(row) sqrt(sum((row - baseline_v)^2)))
dists[rownames(cand_valid) == rownames(baseline)] <- -Inf
robust_distant <- cand_valid[which.max(dists), ]

# Assemble and save
cols_keep <- c("mu_eps", "phi_eps", "g0", "eta", "beta0", "Upsilon",
               "best_shape", "n_in_bounds", "mae",
               "diff_low", "diff_medium", "diff_high",
               "traj_low", "traj_medium", "traj_high")

sensitivity_sets <- rbind(
  data.frame(baseline[,       cols_keep], role = "baseline"),
  data.frame(robust_low[,     cols_keep], role = "robust_low"),
  data.frame(robust_distant[, cols_keep], role = "robust_distant")
)

cat("\n=== THREE PARAMETER SETS ===\n")
print(sensitivity_sets[, c("role", "mu_eps", "g0", "eta", "beta0", "Upsilon",
                             "best_shape", "mae",
                             "diff_low", "diff_medium", "diff_high")])

ab_values <- mu_phi_to_alpha_beta(sensitivity_sets$mu_eps,
                                 sensitivity_sets$phi_eps)

sensitivity_sets$beta_alpha <- ab_values$alpha
sensitivity_sets$beta_beta  <- ab_values$beta

#write.csv(sensitivity_sets,
#          file.path(mydirection, "output/sensitivity_sets_final.csv"),
#          row.names = TRUE)

sensitivity_sets=sensitivity_sets[,c("beta_alpha","beta_beta","g0","beta0","eta","Upsilon","role")]
colnames(sensitivity_sets)=c("beta_alpha","beta_beta","growth0","beta0","eta","Upsilon","calibration_id")

sensitivity_sets$calibration_id[sensitivity_sets$calibration_id=="baseline"] = 1
sensitivity_sets$calibration_id[sensitivity_sets$calibration_id=="robust_low"] = 2
sensitivity_sets$calibration_id[sensitivity_sets$calibration_id=="robust_distant"] = 3

write.csv(sensitivity_sets,
          file.path(mydirection_numericalmodel, "parameters/calibration_kriegler.csv"),
          row.names = TRUE)


#==============================================================================
# EXPERT VS BETA FIT PLOT
# Produced after calibration so best_shape is known and can be annotated
#==============================================================================
best_shape_found <- min(good_candidates$best_shape)

df_beta <- expert_distributions

df_beta$sd_fitted <- sqrt(
  df_beta$alpha_fitted * df_beta$beta_fitted /
  ((df_beta$alpha_fitted + df_beta$beta_fitted)^2 *
   (df_beta$alpha_fitted + df_beta$beta_fitted + 1))
)

# Credible interval bands at multiple levels
ranges      <- 100 * best_shape_found
#ranges      <- c(40, 30, 20, 10, 5)
y_levels    <- c(2.7)
range_rects <- lapply(seq_along(ranges), function(i) {
  r    <- ranges[i]
  qmin <- qbeta((1 - r/100)/2, df_beta$alpha_fitted, df_beta$beta_fitted)
  qmax <- qbeta(1 - (1 - r/100)/2, df_beta$alpha_fitted, df_beta$beta_fitted)
  data.frame(
    corridor = df_beta$corridor,
    xmin = qmin, xmax = qmax,
    ymin = y_levels[1], ymax = y_levels[1] + 0.1,
    label    = paste0("[q", 100 * (1 - r/100)/2,
                      ":q", 100 * (1 - (1 - r/100)/2), "]"),
    y_text   = y_levels[1] + 0.05
  )
}) %>% bind_rows()

x_grid       <- seq(0, 1, length.out = 500)
beta_densities <- lapply(1:3, function(i) {
  corridor_name <- c("low", "medium", "high")[i]
  row           <- expert_distributions[i, ]
  data.frame(
    x        = x_grid,
    density  = dbeta(x_grid, row$alpha_fitted, row$beta_fitted),
    corridor = corridor_name
  )
}) %>% bind_rows()
beta_densities$corridor <- factor(beta_densities$corridor,
                                   levels = c("low", "medium", "high"))

x_text <- max(range_rects$xmax) + 0.01
y_top  <- max(range_rects$ymax)

# Annotate best shape on the plot
best_shape_label <- sprintf("Best shape (calibration): %.2f -> [q%.0f:q%.0f]",
                             best_shape_found,
                             100 * (1 - best_shape_found) / 2,
                             100 * (1 - (1 - best_shape_found) / 2))

ggplot() +
  geom_density(data = expert_samples_df,
               aes(x = value, fill = corridor), alpha = 0.3) +
  geom_line(data = beta_densities,
            aes(x = x, y = density, color = corridor), linewidth = 1.5) +
  geom_segment(data = df_beta,
               aes(x = mean_empirical, xend = mean_empirical,
                   y = 0, yend = y_top, color = corridor), linewidth = 1.2) +
  geom_rect(data = range_rects,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = corridor),
            alpha = 0.35, color = NA) +
  geom_text(data = range_rects %>% distinct(label, y_text),
            aes(x = x_text, y = y_text, label = label),
            hjust = 0, size = 5) +
  annotate("text", x = 0.02, y = 2.95, label = best_shape_label,
           hjust = 0, size = 4.5, color = "black", fontface = "italic") +
  scale_fill_manual(values  = c("low" = "#4DAF4A", "medium" = "#FF7F00",
                                 "high" = "#E41A1C"), name = "Corridor") +
  scale_color_manual(values = c("low" = "#4DAF4A", "medium" = "#FF7F00",
                                 "high" = "#E41A1C"), name = "Corridor") +
  coord_cartesian(ylim = c(0, 3), xlim = c(0, 1.05)) +
  labs(title    = "Expert opinions: empirical mixture vs fitted Beta",
       subtitle = paste("Rectangles = credible intervals |",
                        "Vertical lines = mean |", best_shape_label),
       x = "Probability of partial dieback", y = "Density") +
  theme_minimal() +
  theme(legend.position = "top",
        axis.title  = element_text(size = 16),
        axis.text   = element_text(size = 14),
        plot.title  = element_text(size = 18, face = "bold"),
        plot.subtitle = element_text(size = 13))

ggsave(file.path(mydirection, "figure/expert_vs_beta_fit.png"),
       width = 10, height = 6)

#==============================================================================
# VALIDATION: simulated dieback frequencies vs Kriegler targets
# Averages over n_temp_points temperature draws within each corridor
# (consistent with calibration procedure)
#==============================================================================
targets_kriegler <- setNames(expert_distributions$mean_empirical,
                              expert_distributions$corridor)

n_sim            <- 10000
plots_validation <- list()
set.seed(123)

role_names  <- c("baseline", "robust_low", "robust_distant")
role_labels <- c("Baseline", "Robust-low", "Robust-distant")

for (traj_id in 1:3) {
  corridor_name <- c("low", "medium", "high")[traj_id]
  target        <- targets_kriegler[corridor_name]

  for (r in seq_along(role_names)) {
    data_val  <- sensitivity_sets[sensitivity_sets$role == role_names[r], ]
    ab_params <- mu_phi_to_alpha_beta(data_val$mu_eps, data_val$phi_eps)
    results_list <- vector("list", n_sim)

    for (i in 1:n_sim) {
      eps        <- epsilon_max * rbeta(length(Times),
                                        ab_params$alpha, ab_params$beta)
      temp_idx_i <- sample(1:10, 1)
      res        <- tryCatch(
        process(eps, ab,
                growth  = data_val$g0, Upsilon = data_val$Upsilon,
                beta0   = data_val$beta0, eta   = data_val$eta,
                drop    = 1, traj = traj_id, temp_idx = temp_idx_i),
        error = function(e) NULL)
      if (!is.null(res)) results_list[[i]] <- res[, c("time", "x")]
    }

    df <- as.data.frame(do.call(rbind, lapply(seq_along(results_list), function(i) {
      if (!is.null(results_list[[i]])) cbind(results_list[[i]], Sim = i)
    })))
    colnames(df) <- c("Time", "x", "Sim")

    dieback_sims <- df %>%
      filter(Time > 2) %>%
      group_by(Sim) %>%
      summarise(has_dieback = any(x < 0.5), .groups = "drop")

    pct_dieback <- round(100 * mean(dieback_sims$has_dieback), 1)

    df         <- left_join(df, dieback_sims, by = "Sim")
    summary_df <- df %>%
      group_by(Time) %>%
      summarise(mean_x = mean(x), q10 = quantile(x, 0.10),
                q90 = quantile(x, 0.90), .groups = "drop")

    cat(sprintf("Corridor %-6s | %-15s | dieback: %.1f%% | target: %.1f%%\n",
                corridor_name, role_names[r], pct_dieback, 100 * target))

    # Index: row = role (1-3), col = corridor (1-3)
    # Layout: role varie en ligne, corridor en colonne
    plot_idx <- (r - 1) * 3 + traj_id

    plots_validation[[plot_idx]] <- ggplot() +
      geom_line(
        data  = df %>% filter(!has_dieback) %>%
                filter(Sim %in% sample(unique(Sim), min(200, n_sim))),
        aes(x = Time, y = x, group = Sim),
        color = "grey70", alpha = 0.15, linewidth = 0.3) +
      geom_line(
        data  = df %>% filter(has_dieback),
        aes(x = Time, y = x, group = Sim),
        color = "red", alpha = 0.4, linewidth = 0.4) +
      geom_ribbon(
        data = summary_df,
        aes(x = Time, ymin = q10, ymax = q90),
        fill = "blue", alpha = 0.15) +
      geom_line(
        data  = summary_df,
        aes(x = Time, y = mean_x),
        color = "blue", linewidth = 1.2) +
      geom_hline(yintercept = 0.5, linetype = "dotted",
                 color = "darkred", linewidth = 1.5) +
      ylim(0, 1) +
      labs(
        y     = "x(t)",
        title = paste0(role_labels[r], " | Corridor ", corridor_name,
                       " | dieback: ", pct_dieback,
                       "% | target: ", round(100 * target, 1), "%")
      ) +
      theme_minimal() +
      theme(plot.title = element_text(size = 10))
  }
}

combined_validation <- (plots_validation[[1]] | plots_validation[[2]] | plots_validation[[3]]) /
                       (plots_validation[[4]] | plots_validation[[5]] | plots_validation[[6]]) /
                       (plots_validation[[7]] | plots_validation[[8]] | plots_validation[[9]])

ggsave(file.path(mydirection, "figure/validation_dieback.pdf"),
       plot = combined_validation, width = 20, height = 18)


# PLOT 2 : TEMPERATURE PATHS SSP  (correction : epsilon varie dans le temps)
# On utilise le jeu de paramètres "baseline" pour les simulations SSP
#==============================================================================
cat("=== Generating temperature_paths_ssp.pdf ===\n")

#raw temperature data (not centered on kriegler)
temp_ssp <- read.csv(file.path(mydirection, "data/temp_ssp_ipcc2200.csv"),
                     stringsAsFactors = FALSE)

#anomaly wrt first year
temp_ssp <- temp_ssp %>%
  mutate(across(everything(), ~ . - first(.)))

#scaling to Amazon temperature (cf. Leduc et al. 2016)
temp_ssp <- temp_ssp * (2 / 1.65)

ssp_names <- c("SSP1.2.6", "SSP4.6.0", "SSP5.8.5")

#parameters
data_baseline <- sensitivity_sets[sensitivity_sets$role == "baseline", ]
ab_params_ssp <- mu_phi_to_alpha_beta(data_baseline$mu_eps,
                                       data_baseline$phi_eps)

#local version ode with temperature path
ab_ssp <- function(time, state, parms) {
  with(as.list(c(state, parms)), {
    x         <- min(max(x, 0), 1)
    time_idx  <- max(1, min(floor(time), length(epsilon)))
    eps_t     <- epsilon[time_idx]
    temp_val  <- temp_vec[time_idx]          # vecteur de température SSP
    kappa_val <- kappa_interp(time)
    check <- g0 * (1 - drop * ((feedback * (1 - x)) / beta0)^eta) * x * (1 - x) -
             (eps_t * log(max(temp_val, 0) + 1)) * x - x * kappa_val
    return(list(c(check)))
  })
}
ab_ssp <- compiler::cmpfun(ab_ssp, options = list(optimize = 3))

process_ssp <- function(epsilon_chain, temp_vec, growth, Upsilon,
                         beta0, eta, drop) {
  yini <- c(x = 1)
  pars <- list(
    epsilon  = epsilon_chain, temp_vec = temp_vec,
    g0       = growth,        feedback = Upsilon,
    beta0    = beta0,         eta      = eta,
    drop     = drop
  )
  deSolve::lsoda(y = yini, times = Times, func = ab_ssp, parms = pars,
                 rtol = 1e-4, atol = 1e-6)
}
process_ssp <- compiler::cmpfun(process_ssp, options = list(optimize = 3))

n_sim_ssp <- 1000
all_df    <- data.frame()

for (traj_id in ssp_names) {
  col_idx  <- which(colnames(temp_ssp) == traj_id)
  temp_vec_ssp <- temp_ssp[[col_idx]][1:length(Times)]

  results_list <- vector("list", n_sim_ssp)

  for (i in 1:n_sim_ssp) {
    if (i %% 100 == 0)
      cat(sprintf("  traj = %s | sim = %d/%d\n", traj_id, i, n_sim_ssp))

    epsilon_chain <- epsilon_max * rbeta(length(Times),
                                         ab_params_ssp$alpha,
                                         ab_params_ssp$beta)

    res <- tryCatch(
      process_ssp(epsilon_chain, temp_vec_ssp,
                  growth  = data_baseline$g0,
                  Upsilon = data_baseline$Upsilon,
                  beta0   = data_baseline$beta0,
                  eta     = data_baseline$eta,
                  drop    = 1),
      error = function(e) NULL
    )
    if (!is.null(res)) results_list[[i]] <- res[, "x"]
  }

  valid_sims <- Filter(Negate(is.null), results_list)
  df_wide    <- as.data.frame(do.call(cbind, valid_sims))
  df_wide$Time <- Times

  df_long <- pivot_longer(df_wide, cols = -Time,
                           names_to = "Sim", values_to = "x")

  summary_df <- df_long %>%
    group_by(Time) %>%
    summarise(
      mean_x  = mean(x),
      minimum = min(x),
      maximum = max(x),
      q05     = quantile(x, 0.05),
      q10     = quantile(x, 0.10),
      q90     = quantile(x, 0.90),
      q95     = quantile(x, 0.95),
      .groups = "drop"
    ) %>%
    mutate(traj_id = factor(traj_id, levels = ssp_names))

  all_df <- bind_rows(all_df, summary_df)
}

#75GtC
all_df <- all_df %>%
  mutate(
    mean_x  = 75 - 75 * mean_x,
    minimum = 75 - 75 * minimum,
    maximum = 75 - 75 * maximum,
    q10     = 75 - 75 * q10,
    q90     = 75 - 75 * q90
  )

combined_plot <- ggplot(all_df,
                        aes(x       = 2015 + Time * 5,
                            color   = traj_id,
                            fill    = traj_id)) +

  geom_ribbon(aes(ymin  = minimum,
                  ymax  = maximum,
                  alpha = "Uncertainty (1 000 stochastic paths)"),
              color = NA) +

  geom_line(aes(y       = mean_x,
                linetype = "Average"),
            linewidth = 1.3) +

  geom_hline(yintercept = 75 * 0.5,
             linetype = "dashed", color = "black", linewidth = 0.8) +
  geom_hline(yintercept = 75,
             linetype = "solid",  color = "black", linewidth = 0.8) +

  scale_color_brewer(palette = "Set1", name = "Trajectory") +
  scale_fill_brewer(palette  = "Set1", name = "Trajectory") +
  scale_linetype_manual(name   = "",
                        values = c("Average" = "solid")) +
  scale_alpha_manual(name   = "",
                     values = c("Uncertainty (1 000 stochastic paths)" = 0.2)) +

  coord_cartesian(ylim = c(0, 75)) +

  labs(y = "Net cumulative carbon losses (GtC)", x = "Year") +

  theme_minimal(base_size = 20) +
  theme(
    axis.text        = element_text(size = 20),
    axis.title       = element_text(size = 22, face = "bold"),
    legend.title     = element_text(size = 20, face = "bold"),
    legend.text      = element_text(size = 18),
    legend.position  = "top",
    legend.key.size  = unit(1.2, "cm"),
    panel.grid.minor = element_blank()
  )

ggsave(file.path(mydirection, "figure/temperature_paths_ssp.pdf"),
       plot = combined_plot, width = 16, height = 10)
cat("  -> output/temperature_paths_ssp.pdf saved\n")


cat("\nCalibration complete. Outputs saved to output/\n")
cat("  - good_candidates_step1.csv\n")
cat("  - sensitivity_sets_final.csv\n")
cat("  - expert_vs_beta_fit.png\n")
cat("  - validation_dieback.pdf\n")

#write a table with parameter values

df <- read_csv(file.path(mydirection_numericalmodel, "parameters/calibration_kriegler.csv")) %>%
  dplyr::select(calibration_id, beta_alpha, beta_beta, growth0, beta0, eta, Upsilon) %>%
  arrange(as.numeric(calibration_id))

# Arrondi
df_rounded <- df %>%
  mutate(across(where(is.numeric), ~ signif(., 4)))

# Noms de colonnes en LaTe
col_names <- c(
  "ID",
  "$\\alpha$",
  "$\\beta$",
  "$g_0$",
  "$\\beta_0$",
  "$\\eta$",
  "$\\Upsilon$"
)
# Génération de la table LaTeX
latex_table <- df_rounded %>%
  kable(
    format   = "latex",
    booktabs = TRUE,
    col.names = col_names,
    align    = c("c", rep("r", 6)),
    caption  = "Parameter values",
    label    = "calibration_kriegler",
    escape   = FALSE
  ) %>%
  kable_styling(latex_options = c("hold_position"))

# Sauvegarde
output_path <- file.path(mydirection, "output/calibration_kriegler.tex")
writeLines(as.character(latex_table), output_path)


stopCluster(cl)

