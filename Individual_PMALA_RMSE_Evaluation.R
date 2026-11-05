library(glmgen)
library(ggplot2)
library(tidyr)
library(dplyr)

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
#baseline <- sin(4 * pi * x / n_data) + cos(15 * pi * x / n_data)
baseline <- sqrt(x/n_data)
sd_noise <- 0.4
y <- baseline + rnorm(n_data, 0, sd_noise)

# Parameters (SAME AS ORIGINAL INDIVIDUAL P-MALA)
order <- 2
lambda <- 1000  # FIXED lambda (key difference from Fast Gibbs)
eps <- 0.01
gamma <- 0.1
delta <- 0.1

# P-MALA sampler settings
pmala_iter <- 4000
pmala_burn <- 1000

cat("=== INDIVIDUAL P-MALA QUANTILE RMSE EVALUATION ===\n")
cat("Data size:", n_data, "\n")
cat("P-MALA iterations:", pmala_iter, "| Burn-in:", pmala_burn, "\n")
cat("Fixed lambda:", lambda, "(no adaptation)\n")
cat("Baseline function: sqrt(x/n)\n")
cat("Noise SD:", sd_noise, "\n\n")

# =============================================================================
# ALL FUNCTIONS FROM INDIVIDUAL_PMALA_EXTENDED
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
  fit <- glmgen::trendfilter(x = seq_along(z), y = z, k = order, lambda = eff,
                             control = trendfilter.control.list(max_iter = 1000L))
  as.numeric(fit$beta)
}

grad_f_plus_g <- function(y, theta, tau, eps, gamma, D, lambda, order = 2) {
  grad_f_theta <- grad_smooth_check(y, theta, tau, eps)
  x_tilde <- theta - gamma * grad_f_theta
  prox_step <- prox_g(x_tilde, D, gamma, lambda, order = order)
  (theta - prox_step) / gamma
}

# P-MALA SAMPLER (EXTENDED VERSION WITH ALL THREE ESTIMATES)
pmala_sampler <- function(y, tau = 0.5, lambda = 50, gamma = 0.1,
                          D, delta = 0.1, n_iter = 5000,
                          eps = 0.01, burn_in = 1000, theta_init = NULL, 
                          order = 2, verbose = FALSE) {
  n <- length(y)
  theta_chain <- matrix(NA, nrow = n_iter, ncol = n)
  samples <- matrix(NA, nrow = n_iter - burn_in, ncol = n)

  if (is.null(theta_init)) {
    theta_current <- rep(0, n)
  } else {
    theta_current <- theta_init
  }

  acceptance_count <- 0

  for (i in 1:n_iter) {
    grad <- grad_smooth_check(y, theta_current, tau, eps)
    z <- theta_current - delta * grad + sqrt(2 * delta) * rnorm(n)
    theta_proposal <- prox_g(z, D, gamma, lambda, order = order)

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

    if (log(runif(1)) < log_alpha) {
      theta_current <- theta_proposal
      if (i > burn_in) acceptance_count <- acceptance_count + 1
    }

    theta_chain[i, ] <- theta_current

    if (i > burn_in) {
      samples[i - burn_in, ] <- theta_current
    }

    if (verbose && i %% 1000 == 0) {
      if (i <= burn_in) {
        cat("  P-MALA burn-in:", i, "\n")
      } else {
        cat("  P-MALA sampling:", i - burn_in,
            "Acceptance:", round(acceptance_count/(i - burn_in), 3), "\n")
      }
    }
  }

  theta_post <- theta_chain[(burn_in + 1):n_iter, , drop = FALSE]

  # Compute all three estimates - Mean, Median, and Mode
  theta_mean_estimate <- colMeans(theta_post)
  theta_median_estimate <- apply(theta_post, 2, median)
  theta_mode_estimate <- apply(theta_post, 2, compute_mode)

  list(samples = theta_post,
       estimate = theta_mean_estimate,  # Keep this name for backward compatibility
       theta_mean_estimate = theta_mean_estimate,
       theta_median_estimate = theta_median_estimate,
       theta_mode_estimate = theta_mode_estimate,
       acceptance_rate = acceptance_count / (n_iter - burn_in))
}

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
  acceptance_rate = numeric(n_quantiles),
  best_estimator = character(n_quantiles),
  stringsAsFactors = FALSE
)

cat("Starting Individual P-MALA evaluation loop...\n\n")

# Loop over all quantiles
for (i in 1:n_quantiles) {
  tau <- quantiles[i]
  
  cat(sprintf("=== EVALUATING QUANTILE τ = %.1f (%d/%d) ===\n", tau, i, n_quantiles))
  
  # Compute true theoretical quantile for current tau
  Q_true <- baseline + qnorm(tau) * sd_noise
  
  # Run P-MALA sampler for this quantile (fixed lambda)
  result <- pmala_sampler(
    y = y,
    tau = tau,
    lambda = lambda,  # FIXED - no adaptation
    gamma = gamma,
    D = D,
    delta = delta,
    n_iter = pmala_iter,
    eps = eps,
    burn_in = pmala_burn,
    theta_init = y,
    order = order,
    verbose = FALSE
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
  results_table[i, "acceptance_rate"] <- result$acceptance_rate
  results_table[i, "best_estimator"] <- best_estimator
  
  # Progress report
  cat(sprintf("  RMSE - Mean: %.4f | Median: %.4f | Mode: %.4f\n", 
              rmse_mean, rmse_median, rmse_mode))
  cat(sprintf("  Best: %s (%.4f) | λ: %d (fixed) | Accept: %.3f\n\n", 
              best_estimator, min(rmse_values), lambda, result$acceptance_rate))
}

# =============================================================================
# FINAL RESULTS TABLE
# =============================================================================

cat("=== INDIVIDUAL P-MALA QUANTILE EVALUATION COMPLETE ===\n\n")

cat("INDIVIDUAL P-MALA QUANTILE EVALUATION TABLE:\n")
cat("============================================\n")
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

# Acceptance rate behavior (since lambda is fixed)
cat(sprintf("\nAcceptance rates - Mean: %.3f | Range: %.3f to %.3f\n",
            mean(results_table$acceptance_rate),
            min(results_table$acceptance_rate),
            max(results_table$acceptance_rate)))

cat("\n=== KEY DIFFERENCES FROM FAST GIBBS ===\n")
cat("• Fixed lambda = ", lambda, " (no adaptation)\n")
cat("• Pure P-MALA sampling (no Gibbs component)\n") 
cat("• Consistent penalty across all quantiles\n")
cat("• Different computational profile\n")

cat("\n=== EVALUATION COMPLETE ===\n")
cat("Results stored in 'results_table' data frame.\n")
cat("Ready for comparison with Fast Gibbs results.\n")

# Save results for potential later use
write.csv(results_table, "individual_pmala_quantile_rmse_evaluation.csv", row.names = FALSE)
cat("Results also saved to: individual_pmala_quantile_rmse_evaluation.csv\n")

# =============================================================================
# COMPARISON PREPARATION
# =============================================================================

cat("\n=== COMPARISON ANALYSIS PREP ===\n")
cat("To compare with Fast Gibbs results:\n")
cat("1. Load both CSV files\n")
cat("2. Compare average RMSE by estimator type\n") 
cat("3. Analyze fixed vs adaptive lambda performance\n")
cat("4. Check acceptance rate stability\n")
cat("5. Identify quantile-specific patterns\n\n")

cat("Individual P-MALA evaluation with fixed λ =", lambda, "complete!\n")
