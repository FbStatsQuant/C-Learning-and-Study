library(glmgen)
library(ggplot2)
library(tidyr)
library(dplyr)
library(zoo)

# Function to compute mode (most frequent value using density estimate)
compute_mode <- function(x) {
  d <- density(x)
  d$x[which.max(d$y)]
}

# =============================================================================
# DATA GENERATION (SAME AS ORIGINAL)
# =============================================================================

set.seed(100)
n_data <- 1000
x <- seq(1, n_data, 1)
#baseline <- sin(2 * pi * x / n_data)
#baseline <- sin(10 * pi * x / n_data) + cos(15 * pi * x / n_data)
baseline <- sqrt(x/n_data)
sd_noise <- 0.4
y <- baseline + rnorm(n_data, 0, sd_noise)

# Parameters (SAME AS ORIGINAL)
order <- 2
eps <- 0.01
delta <- 0.02
gamma <- 0.1

# Prior parameters
prior_lambda_mean <- 3 * n_data 
prior_lambda_sd <- 1*n_data    
a <- (prior_lambda_mean / prior_lambda_sd)^2
b <- prior_lambda_mean / (prior_lambda_sd^2)

# Gibbs sampler settings
gibbs_iter <- 5000
gibbs_burn <- 2000

cat("=== QUANTILE RMSE EVALUATION ACROSS τ = 0.1, 0.2, ..., 0.9 ===\n")
cat("Data size:", n_data, "\n")
cat("Gibbs iterations:", gibbs_iter, "| Burn-in:", gibbs_burn, "\n")
cat("Baseline function: sqrt(x/n)\n")
cat("Noise SD:", sd_noise, "\n\n")

# =============================================================================
# ALL FUNCTIONS FROM FAST_GIBBS_EXTENDED
# =============================================================================

D_matrix <- function(n, order) {
  D <- diag(n)
  for (i in 1:(order + 1)) {
    D <- diff(D)
  }
  D
}

D <- D_matrix(n_data, order = order)

rho_tau <- function(u, tau) ifelse(u >= 0, tau*u, (tau-1)*u)

true_check_loss <- function(u, tau) {
  loss <- rho_tau(u, tau)
  (loss)
}

smooth_check_loss <- function(u, tau, eps = eps) {
  tau * u + eps * log(1 + exp(-u / eps))
}

grad_smooth_check <- function(y, theta, tau, eps = eps) {
  u <- y - theta
  sigmoid <- 1 / (1 + exp(u / eps))
  -(tau - sigmoid)
}

neg_log_posterior_true <- function(theta, y, tau, lambda, D, eps) {
  resid <- y - theta
  f_theta <- sum(true_check_loss(resid, tau))
  g_theta <- lambda * sum(abs(D %*% theta))
  f_theta + g_theta
}

prox_g <- function(z, D, gamma, lambda, order = 2) {
  eff <- gamma * lambda
  fit <- glmgen::trendfilter(x = seq_along(z), y = z, k = order, lambda = eff)
  as.numeric(fit$beta)
}

grad_f_plus_g <- function(y, theta, tau, eps, gamma, D, lambda, order = 2) {
  grad_f_theta <- grad_smooth_check(y, theta, tau, eps)
  x_tilde <- theta - gamma * grad_f_theta
  prox_step <- prox_g(x_tilde, D, gamma, lambda, order = order)
  (theta - prox_step) / gamma
}

# P-MALA SINGLE SAMPLE (no burn-in, returns proposal and acceptance decision)
pmala_single_sample <- function(y, tau = 0.5, lambda = 50, gamma = 0.1,
                                D, delta = 0.1, eps = 0.01, 
                                theta_current = NULL, order = 2) {
  n <- length(y)
  
  if (is.null(theta_current)) {
    theta_current <- rep(0, n)
  }
  
  # Generate proposal
  grad <- grad_smooth_check(y, theta_current, tau, eps)
  z <- theta_current - delta * grad + sqrt(2 * delta) * rnorm(n)
  theta_proposal <- prox_g(z, D, gamma, lambda, order = order)
  
  # Compute acceptance probability
  log_posterior_current <- -neg_log_posterior_true(theta_current, y, tau, lambda, D, eps)
  log_posterior_proposal <- -neg_log_posterior_true(theta_proposal, y, tau, lambda, D, eps)
  log_posterior_ratio <- log_posterior_proposal - log_posterior_current
  
  grad_psi_current <- grad_f_plus_g(y, theta_current, tau, eps, gamma, D, lambda, order)
  mu_current <- theta_current - delta * grad_psi_current
  
  grad_psi_proposal <- grad_f_plus_g(y, theta_proposal, tau, eps, gamma, D, lambda, order)
  mu_proposal <- theta_proposal - delta * grad_psi_proposal
  
  log_q_proposal_given_current <- -sum((theta_proposal - mu_current)^2) / (4 * delta)
  log_q_current_given_proposal <- -sum((theta_current - mu_proposal)^2) / (4 * delta)
  log_proposal_ratio <- log_q_current_given_proposal - log_q_proposal_given_current
  
  log_alpha <- min(0, log_posterior_ratio + log_proposal_ratio)
  
  # Accept/reject
  accepted <- (log(runif(1)) < log_alpha)
  
  if (accepted) {
    theta_new <- theta_proposal
  } else {
    theta_new <- theta_current
  }
  
  list(theta = theta_new, accepted = accepted)
}

# GIBBS SAMPLER WITH GAMMA PRIOR
gibbs_gamma_prior <- function(y, tau = 0.5, a = 2, b = 0.002,
                              gamma = 0.1, D, delta = 0.1,
                              gibbs_iter = 30, gibbs_burn = 10,
                              eps = 0.01, theta_init = NULL, order = 2) {
  
  n <- length(y)
  
  # Storage for Gibbs chain
  gibbs_post <- gibbs_iter - gibbs_burn
  theta_gibbs_samples <- matrix(NA, nrow = gibbs_post, ncol = n)
  lambda_gibbs_samples <- numeric(gibbs_post)
  
  # Storage for full chain (including burn-in)
  theta_gibbs_chain <- matrix(NA, nrow = gibbs_iter, ncol = n)
  lambda_gibbs_chain <- numeric(gibbs_iter)
  
  # Acceptance tracking
  acceptance_count <- 0
  
  # Initialize
  if (is.null(theta_init)) {
    theta_current <- y
  } else {
    theta_current <- theta_init
  }
  lambda_current <- 1000  # Fixed initial lambda
  
  for (i in 1:gibbs_iter) {
    # STEP 1: Sample theta given lambda (Single P-MALA sample)
    pmala_result <- pmala_single_sample(
      y = y,
      tau = tau,
      lambda = lambda_current,
      gamma = gamma,
      D = D,
      delta = delta,
      eps = eps,
      theta_current = theta_current,
      order = order
    )
    
    theta_current <- pmala_result$theta
    
    # Track acceptance
    if (pmala_result$accepted) {
      acceptance_count <- acceptance_count + 1
    }
    
    # STEP 2: Sample lambda given theta (Gamma posterior)
    S <- sum(abs(D %*% theta_current))
    lambda_current <- rgamma(1, shape = a+1, rate = b + S)
    
    # Store samples
    theta_gibbs_chain[i, ] <- theta_current
    lambda_gibbs_chain[i] <- lambda_current
    
    if (i > gibbs_burn) {
      sample_idx <- i - gibbs_burn
      theta_gibbs_samples[sample_idx, ] <- theta_current
      lambda_gibbs_samples[sample_idx] <- lambda_current
    }
  }
  
  # Final estimates - Mean, Median, and Mode
  theta_mean_estimate <- colMeans(theta_gibbs_samples)
  theta_median_estimate <- apply(theta_gibbs_samples, 2, median)
  theta_mode_estimate <- apply(theta_gibbs_samples, 2, compute_mode)
  lambda_final_estimate <- mean(lambda_gibbs_samples)
  
  list(
    theta_estimate = theta_mean_estimate,  # Keep this name for backward compatibility
    theta_mean_estimate = theta_mean_estimate,
    theta_median_estimate = theta_median_estimate,
    theta_mode_estimate = theta_mode_estimate,
    lambda_estimate = lambda_final_estimate,
    theta_samples = theta_gibbs_samples,
    lambda_samples = lambda_gibbs_samples,
    theta_chain = theta_gibbs_chain,
    lambda_chain = lambda_gibbs_chain,
    acceptance_rate = acceptance_count / gibbs_iter
  )
}

# Initialize smart starting values
y_smooth <- rollmedian(y, k=21, fill="extend")
noise_sd_est <- mad(diff(y)) / sqrt(2)

# =============================================================================
# MAIN EVALUATION LOOP
# =============================================================================

# Define quantiles to evaluate
quantiles <- seq(0.1, 0.9, 0.1)
n_quantiles <- length(quantiles)

# Initialize results storage
results_table <- data.frame(
  quantile = quantiles,
  rmse_mean = numeric(n_quantiles),
  rmse_median = numeric(n_quantiles),
  rmse_mode = numeric(n_quantiles),
  lambda_estimate = numeric(n_quantiles),
  acceptance_rate = numeric(n_quantiles),
  best_estimator = character(n_quantiles),
  stringsAsFactors = FALSE
)

cat("Starting evaluation loop...\n\n")

# Loop over all quantiles
for (i in 1:n_quantiles) {
  tau <- quantiles[i]
  
  cat(sprintf("=== EVALUATING QUANTILE τ = %.1f (%d/%d) ===\n", tau, i, n_quantiles))
  
  # Compute true theoretical quantile for current tau
  Q_true <- baseline + qnorm(tau) * sd_noise
  
  # Smart initialization for this quantile
  theta_init_smart <- y_smooth + qnorm(tau) * noise_sd_est
  
  # Run Gibbs sampler for this quantile
  result <- gibbs_gamma_prior(
    y = y,
    tau = tau,
    a = a,
    b = b,
    gamma = gamma,
    D = D,
    delta = delta,
    gibbs_iter = gibbs_iter,
    gibbs_burn = gibbs_burn,
    eps = eps,
    theta_init = theta_init_smart,
    order = order
  )
  
  # Compute RMSE for all three estimates
  rmse_mean <- sqrt(mean((Q_true - result$theta_mean_estimate)^2))
  rmse_median <- sqrt(mean((Q_true - result$theta_median_estimate)^2))
  rmse_mode <- sqrt(mean((Q_true - result$theta_mode_estimate)^2))
  
  # Identify best estimator
  rmse_values <- c(rmse_mean, rmse_median, rmse_mode)
  estimator_names <- c("Mean", "Median", "Mode")
  best_idx <- which.min(rmse_values)
  best_estimator <- estimator_names[best_idx]
  
  # Store results
  results_table[i, "rmse_mean"] <- rmse_mean
  results_table[i, "rmse_median"] <- rmse_median
  results_table[i, "rmse_mode"] <- rmse_mode
  results_table[i, "lambda_estimate"] <- result$lambda_estimate
  results_table[i, "acceptance_rate"] <- result$acceptance_rate
  results_table[i, "best_estimator"] <- best_estimator
  
  # Progress report
  cat(sprintf("  RMSE - Mean: %.4f | Median: %.4f | Mode: %.4f\n", 
              rmse_mean, rmse_median, rmse_mode))
  cat(sprintf("  Best: %s (%.4f) | λ: %.1f | Accept: %.3f\n\n", 
              best_estimator, min(rmse_values), result$lambda_estimate, result$acceptance_rate))
}

# =============================================================================
# FINAL RESULTS TABLE
# =============================================================================

cat("=== COMPREHENSIVE RMSE EVALUATION COMPLETE ===\n\n")

cat("FAST GIBBS QUANTILE EVALUATION TABLE:\n")
cat("=====================================\n")
print(results_table, digits = 4, row.names = FALSE)

cat("\n=== SUMMARY STATISTICS ===\n")

# Overall performance summary
mean_performance <- colMeans(results_table[, c("rmse_mean", "rmse_median", "rmse_mode")])
cat("Average RMSE across all quantiles:\n")
cat(sprintf("  Mean estimator:   %.4f\n", mean_performance["rmse_mean"]))
cat(sprintf("  Median estimator: %.4f\n", mean_performance["rmse_median"]))
cat(sprintf("  Mode estimator:   %.4f\n", mean_performance["rmse_mode"]))

# Best overall estimator
overall_best <- names(mean_performance)[which.min(mean_performance)]
cat(sprintf("\nOverall best estimator: %s (avg RMSE = %.4f)\n", 
            gsub("rmse_", "", overall_best), min(mean_performance)))

# Count wins by estimator
estimator_wins <- table(results_table$best_estimator)
cat("\nWins by estimator:\n")
print(estimator_wins)

# Quantile-specific insights
extreme_quantiles <- results_table[results_table$quantile %in% c(0.1, 0.9), ]
middle_quantiles <- results_table[results_table$quantile %in% c(0.4, 0.5, 0.6), ]

cat("\nPerformance by quantile region:\n")
cat("Extreme quantiles (0.1, 0.9) - Avg RMSE:\n")
cat(sprintf("  Mean: %.4f | Median: %.4f | Mode: %.4f\n",
            mean(extreme_quantiles$rmse_mean),
            mean(extreme_quantiles$rmse_median),
            mean(extreme_quantiles$rmse_mode)))

cat("Middle quantiles (0.4, 0.5, 0.6) - Avg RMSE:\n")
cat(sprintf("  Mean: %.4f | Median: %.4f | Mode: %.4f\n",
            mean(middle_quantiles$rmse_mean),
            mean(middle_quantiles$rmse_median),
            mean(middle_quantiles$rmse_mode)))

# Lambda behavior
cat(sprintf("\nLambda estimates - Mean: %.1f | Range: %.1f to %.1f\n",
            mean(results_table$lambda_estimate),
            min(results_table$lambda_estimate),
            max(results_table$lambda_estimate)))

cat(sprintf("Acceptance rates - Mean: %.3f | Range: %.3f to %.3f\n",
            mean(results_table$acceptance_rate),
            min(results_table$acceptance_rate),
            max(results_table$acceptance_rate)))

cat("\n=== EVALUATION COMPLETE ===\n")
cat("Results stored in 'results_table' data frame.\n")
cat("Ready for further analysis or export.\n")

# Save results for potential later use
write.csv(results_table, "fast_gibbs_quantile_rmse_evaluation.csv", row.names = FALSE)
cat("Results also saved to: fast_gibbs_quantile_rmse_evaluation.csv\n")
