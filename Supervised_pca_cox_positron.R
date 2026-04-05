# Collective Estimation Simulation with High-dim Scalers Data in Survival Analysis
# Load necessary libraries
library(survival)
library(tidyverse)

#' Supervised Dimension Reduction Cox Model
#' 
#' @param Z Numeric matrix of scalar covariates (n x p)
#' @param time Numeric vector of survival times (n x 1)
#' @param status Numeric vector of event indicators (n x 1) (1 for event, 0 for censoring)
#' @param K Integer, number of latent dimensions to extract
#' @param sigma2 Numeric, variance of the data reconstruction (tuning parameter)
#' @param gamma Numeric, learning rate for Gradient Ascent on U
#' @param max_iter Integer, maximum number of alternating optimization iterations
#' @param tol Numeric, tolerance for convergence

supervised_pca_cox <- function(Z, time, status, K, 
                               sigma2 = 1.0, gamma = 0.001, 
                               max_iter = 100, tol = 1e-4) {
    
  n <- nrow(Z)
  p <- ncol(Z)
  
  # Center the data for PCA stability
  Z_centered <- scale(Z, center = TRUE, scale = FALSE)

  # ==========================================
  # STEP 1: INITIALIZATION
  # ==========================================

  # Initialize U using standard unsupervised PCA on the centered data
  svd_Z <- svd(Z_centered)
  U <- svd_Z$v[, 1:K, drop = FALSE]  # U is currently p x K and orthonormal

  # Initialize alpha by fitting Cox model on the initial PCA scores
  Z_proj <- Z_centered %*% U  # n x K
  df <- data.frame(time = time, status = status, Z_proj)
  cox_init <- coxph(Surv(time, status) ~ ., data = df)
  alpha <- matrix(coef(cox_init), ncol = 1)  # K x 1

  cat("Starting alternating optimization...\n")
    
  # ==========================================
  # STEP 2: ALTERNATING OPTIMIZATION LOOP
  # ==========================================

  for (iter in 1:max_iter){
    U_old <- U
    alpha_old <- alpha

    # --- Block A: Update Alpha (Newton-Raphson) --- 
    Z_proj <- Z_centered %*% U  # n x K
    df <- data.frame(time = time, status = status, Z_proj)

    # Fit Cox model holding U fixed
    cox_fit <- coxph(Surv(time, status) ~ ., data = df)
    alpha <- matrix(coef(cox_fit), ncol = 1)  # K x 1

    # Extract Martingale residuals for gradient calculation (Cox gradient)
    d <- matrix(residuals(cox_fit, type = "martingale"), ncol = 1)  # n x 1

    # --- Block B: Update U (Gradient Ascent) ---
    # Compute gradient of the Cox partial likelihood w.r.t. U
    # 1. Data Gradient: (1/sigma2) * Z^T *Z * U
    G_data <- (1/sigma2) * crossprod(Z_centered) %*% U  # p x K

    # 2. Cox Gradient: Z^T * d * alpha^T
    G_cox <- crossprod(Z_centered, d) %*% t(alpha)  # p x K

    # Total Gradient
    G_total <- G_data + G_cox  # p x K
    # Update U with a gradient ascent step
    U_temp <- U + gamma * G_total  # p x K

    # Orthogonalize using QR decomposition to enforce U^T U = I
    qr_res <- qr(U_temp)
    U <- qr.Q(qr_res)  # p x K, orthonormal columns
    R <- qr.R(qr_res)  # K x K, upper triangular

    # Standardize the signs of Q (QR decomposition can arbitrarily flip signs between iterations)
    sign_diag <- sign(diag(R))
    U <- scale(Q, center = FALSE, scale = sign_diag)

    # --- Block C: Check for Convergence ---
    diff_alpha <- max(abs(alpha - alpha_old))
    diff_U <- max(abs(U - U_old))

    if (iter %% 10 == 0 || iter == 1) {
      cat(sprintf("Iter %d: max delta_alpha = %.6f, max delta_U = %.6f\n", 
                  iter, diff_alpha, diff_U))
    }
    
    if (diff_alpha < tol && diff_U < tol) {
      cat(sprintf("=> Converged successfully at iteration %d\n", iter))
      break
    }
  }

  if (iter == max_iter) cat("=> Reached maximum iterations without full convergence.\n")
  
  return(list(
    u = U,
    alpha = alpha,
    iterations = iter,
    Z_proj = Z_centered %*% U,
    final_model = cox_fit
  ))
}


# ==========================================
# STEP 3: SIMULATION TEST BLOCK
# ==========================================
set.seed(42)

# Simulation parameters
n_sim <- 300 # number of patients
p_sim <- 50 # number of scalar covariates (High-dimensional)
K_sim <- 3 # number of latent dimensions (True latent biological phenotypes)

# 1. Generate True Orthogonal Loading Matrix U_true (p x K)
U_true <- qr.Q(qr(matrix(rnorm(p_sim * K_sim), nrow = p_sim, ncol = K_sim)))

# 2. Generate Latent Patient Scores (c_i) (n x K)
C_scores <- matrix(rnorm(n_sim * K_sim), nrow = n_sim, ncol = K_sim)

# 3. Generate Scalar Covariates Z (n x p) with noise (Z = C * U^T + noise)
Z_sim <- C_scores %*% t(U_true) + matrix(rnorm(n_sim * p_sim, sd = 0.5), nrow = n_sim, ncol = p_sim)

# 4. Generate Survival Times using a Cox model with the latent scores
# Let's make the 3rd shape highly lethal, the 1st protective, and 2nd neutral
alpha_true <- c(-1.5, 0, 2.0) 
linear_predictor <- C_scores %*% alpha_true

# Simulate exponential survival times
hazard <- exp(linear_predictor)
time_sim <- rexp(n_sim, rate = hazard)  # Survival times

# Add random censoring
censor_time <- rexp(n_sim, rate = 0.5)  # Censoring times
status_sim <- ifelse(time_sim <= censor_time, 1, 0)  # Event indicator
time_sim <- pmin(time_sim, censor_time)  # Observed times

cat("Running supervised PCA Cox model on simulated data...\n")
results <- supervised_pca_cox(Z = Z_sim, time = time_sim, status = status_sim, K = K_sim, 
                              sigma2 = 1, gamma = 0.001, max_iter = 100) 


cat("\nFinal Estimated Alpha Weights:\n")
print(results$alpha)
cat("\nTrue Alpha Weights were: -1.5, 0.0, 2.0\n")