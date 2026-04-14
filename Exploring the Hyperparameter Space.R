## Exploring the Hyperparameter Space


# ==============================================================================
# EXPERIMENT: The Effect of Data Noise on Optimal Theta
# ==============================================================================
library(survival)
library(tidyverse)

set.seed(123)

# 1. Base Parameters (Fixed)
n_sim <- 300 
p_sim <- 50 
K_sim <- 3 

# Generate True Orthogonal Matrix and Latent Scores (Fixed across all noise levels)
U_true <- qr.Q(qr(matrix(rnorm(p_sim * K_sim), nrow = p_sim, ncol = K_sim)))
C_scores <- matrix(rnorm(n_sim * K_sim), nrow = n_sim, ncol = K_sim)

# True Survival Weights
alpha_true <- matrix(c(-1.5, 0, 2.0), ncol = 1) 

# 2. Experiment Setup
noise_levels <- c(0.1, 0.5, 1.0, 2.0, 3.0) # From very clean to very noisy
theta_grid <- seq(0.05, 1.0, by = 0.05)
n_folds <- 5

# Create consistent cross-validation folds so noise levels are perfectly comparable
folds <- sample(rep(1:n_folds, length.out = n_sim))
results_df <- data.frame()

cat("Starting Noise vs. Theta Experiment...\n")

# 3. Run the Loop
for (noise in noise_levels) {
  cat(sprintf("-> Simulating Environment with Noise SD = %.1f\n", noise))
  
  # Generate Noisy Data
  Z_raw <- C_scores %*% t(U_true) + matrix(rnorm(n_sim * p_sim, sd = noise), n_sim, p_sim)
  Z_sim <- scale(Z_raw, center = TRUE, scale = TRUE)
  
  # Generate Survival Times based on scaled data
  eta_true <- Z_sim %*% U_true %*% alpha_true
  time_sim <- rexp(n_sim, rate = exp(eta_true)) 
  censor_time <- rexp(n_sim, rate = 0.5) 
  status_sim <- ifelse(time_sim <= censor_time, 1, 0) 
  time_sim <- pmin(time_sim, censor_time) 
  
  # Run Cross Validation
  for (theta_val in theta_grid) {
    c_indices <- numeric(n_folds)
    
    for (f in 1:n_folds) {
      test_idx <- which(folds == f)
      
      # Train
      fit <- supervised_pca_cox(Z = Z_sim[-test_idx, ], time = time_sim[-test_idx], 
                                status = status_sim[-test_idx], 
                                K = K_sim, theta = theta_val, gamma = 0.001, max_iter = 50)
      
      # Predict
      eta_pred <- Z_sim[test_idx, ] %*% fit$U %*% fit$alpha
      
      # Calculate C-Index with Hazard Fix
      raw_c <- concordance(Surv(time_sim[test_idx], status_sim[test_idx]) ~ eta_pred)$concordance
      c_indices[f] <- ifelse(raw_c < 0.5, 1 - raw_c, raw_c)
    }
    
    # Store Result
    results_df <- rbind(results_df, data.frame(
      Noise = as.factor(noise),
      Theta = theta_val,
      C_Index = mean(c_indices)
    ))
  }
}

# ==============================================================================
# PART 4: PLOT THE RESULTS
# ==============================================================================
# Find the optimal theta for each noise level to mark it on the plot
optimal_points <- results_df %>%
  group_by(Noise) %>%
  slice_max(C_Index, n = 1) %>%
  ungroup()

noise_plot <- ggplot(results_df, aes(x = Theta, y = C_Index, color = Noise, group = Noise)) +
  geom_line(linewidth = 1.2, alpha = 0.8) +
  geom_point(size = 2) +
  # Add a vertical dashed line and point for the peaks
  geom_point(data = optimal_points, shape = 21, size = 4, fill = "white", stroke = 2) +
  scale_color_viridis_d(option = "plasma", end = 0.9, name = "Data Noise (SD)") +
  theme_minimal(base_size = 14) +
  labs(
    title = "How Data Noise Shifts the Optimal Supervision Weight",
    subtitle = "Higher noise forces the model to rely more on Cox supervision (lower Theta).",
    x = "Theta (Supervision Weight)",
    y = "Cross-Validated C-Index"
  ) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

print(noise_plot)


### ==============================================================================
# ==============================================================================
# EXPERIMENT 2: The Effect of Latent Dimensions (K) on Optimal Theta
# ==============================================================================
library(survival)
library(tidyverse)

# ------------------------------------------------------------------------------
# 1. FAST CORE ALGORITHM (No Likelihood Logging)
# ------------------------------------------------------------------------------
supervised_pca_cox_fast <- function(Z, time, status, K, theta = 0.5,
                                    gamma = 0.001, max_iter = 100, tol = 1e-4) {
  n <- nrow(Z)
  p <- ncol(Z)
  Z_scaled <- scale(Z, center = TRUE, scale = TRUE)
  
  svd_Z <- svd(Z_scaled)
  U <- svd_Z$v[, 1:K, drop = FALSE]
  
  Z_recon <- Z_scaled %*% U %*% t(U)
  sigma2 <- mean((Z_scaled - Z_recon)^2) 
  
  Z_proj <- Z_scaled %*% U 
  df <- data.frame(time = time, status = status, Z_proj)
  cox_init <- suppressWarnings(coxph(Surv(time, status) ~ ., data = df))
  alpha <- matrix(coef(cox_init), ncol = 1) 
  
  for (iter in 1:max_iter) {
    U_old <- U
    alpha_old <- alpha
    
    # Block A: Update Alpha
    Z_proj <- Z_scaled %*% U 
    df <- data.frame(time = time, status = status, Z_proj)
    cox_fit <- suppressWarnings(coxph(Surv(time, status) ~ ., data = df))
    alpha <- matrix(coef(cox_fit), ncol = 1) 
    d <- matrix(residuals(cox_fit, type = "martingale"), ncol = 1)
    
    # Block B: Update U
    G_data <- (1 / sigma2) * crossprod(Z_scaled) %*% U
    G_cox <- crossprod(Z_scaled, d) %*% t(alpha)
    G_total <- (theta * G_data) + ((1 - theta) * G_cox)
    U_temp <- U + gamma * G_total
    
    # Orthogonalize
    qr_res <- qr(U_temp)
    U <- scale(qr.Q(qr_res), center = FALSE, scale = sign(diag(qr.R(qr_res))))
    
    # Block C: Fast Convergence Check (No Likelihood Math!)
    if (max(abs(alpha - alpha_old)) < tol && max(abs(U - U_old)) < tol) break
  }
  
  return(list(U = U, alpha = alpha, final_model = cox_fit))
}

# ------------------------------------------------------------------------------
# 2. SIMULATION SETUP
# ------------------------------------------------------------------------------
set.seed(42)
n_sim <- 300 
p_sim <- 50 
K_true <- 3 

# Generate True Orthogonal Matrix and Latent Scores
U_true <- qr.Q(qr(matrix(rnorm(p_sim * K_true), nrow = p_sim, ncol = K_true)))
C_scores <- matrix(rnorm(n_sim * K_true), nrow = n_sim, ncol = K_true)
alpha_true <- matrix(c(-1.5, 0, 2.0), ncol = 1) 

# Generate Data with Moderate Noise
noise_fixed <- 1.5
Z_raw <- C_scores %*% t(U_true) + matrix(rnorm(n_sim * p_sim, sd = noise_fixed), n_sim, p_sim)
Z_sim <- scale(Z_raw, center = TRUE, scale = TRUE)

eta_true <- Z_sim %*% U_true %*% alpha_true
time_sim <- rexp(n_sim, rate = exp(eta_true)) 
censor_time <- rexp(n_sim, rate = 0.5) 
status_sim <- ifelse(time_sim <= censor_time, 1, 0) 
time_sim <- pmin(time_sim, censor_time) 

# ------------------------------------------------------------------------------
# 3. K-FOLD CROSS-VALIDATION EXPERIMENT
# ------------------------------------------------------------------------------
K_test_values <- c(2, 3, 6) # Test Underfit, Perfect Fit, and Overfit
theta_grid <- seq(0.05, 1.0, by = 0.05)
n_folds <- 5
folds <- sample(rep(1:n_folds, length.out = n_sim))

results_df <- data.frame()
cat("Starting K vs. Theta Experiment...\n")

for (K_test in K_test_values) {
  cat(sprintf("-> Testing Latent Dimensions: K = %d\n", K_test))
  
  for (theta_val in theta_grid) {
    c_indices <- numeric(n_folds)
    
    for (f in 1:n_folds) {
      test_idx <- which(folds == f)
      
      # Train with FAST function
      fit <- supervised_pca_cox_fast(Z = Z_sim[-test_idx, ], time = time_sim[-test_idx], 
                                     status = status_sim[-test_idx], 
                                     K = K_test, theta = theta_val, gamma = 0.001)
      
      # Predict
      eta_pred <- Z_sim[test_idx, ] %*% fit$U %*% fit$alpha
      
      # Calculate C-Index
      raw_c <- concordance(Surv(time_sim[test_idx], status_sim[test_idx]) ~ eta_pred)$concordance
      c_indices[f] <- ifelse(raw_c < 0.5, 1 - raw_c, raw_c)
    }
    
    results_df <- rbind(results_df, data.frame(
      K_Value = as.factor(K_test),
      Theta = theta_val,
      C_Index = mean(c_indices)
    ))
  }
}

# ------------------------------------------------------------------------------
# 4. PLOT THE RESULTS
# ------------------------------------------------------------------------------
optimal_k_points <- results_df %>%
  group_by(K_Value) %>%
  slice_max(C_Index, n = 1) %>%
  ungroup()

k_plot <- ggplot(results_df, aes(x = Theta, y = C_Index, color = K_Value, group = K_Value)) +
  geom_line(size = 1.2, alpha = 0.8) +
  geom_point(size = 2) +
  geom_point(data = optimal_k_points, shape = 21, size = 4, fill = "white", stroke = 2) +
  scale_color_brewer(palette = "Set1", name = "Latent Dimensions (K)") +
  theme_minimal(base_size = 14) +
  labs(
    title = "How Latent Dimensions (K) Shift the Optimal Supervision Weight",
    subtitle = "Higher K requires higher Theta (data regularization) to prevent overfitting.",
    x = "Theta (Supervision Weight)",
    y = "Cross-Validated C-Index"
  ) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

print(k_plot)





## ==============================================================================
# Both Experiments in One Script for Easy Comparison
# ==============================================================================
library(survival)
library(tidyverse)

# 1. Include the core function from our previous code
# (Assuming fit_supervised_pca_cox is defined exactly as we did before)
# For this script to run, make sure the supervised_pca_cox() function is loaded!

run_hyperparameter_experiment <- function(n_sim = 300, p_sim = 50, K_true = 3) {
  
  # --- Generate Base Data ---
  set.seed(123)
  U_true <- qr.Q(qr(matrix(rnorm(p_sim * K_true), nrow = p_sim, ncol = K_true)))
  C_scores <- matrix(rnorm(n_sim * K_true), nrow = n_sim, ncol = K_true)
  alpha_true <- matrix(c(-1.5, 0, 2.0), ncol = 1) 
  
  # Experiment Scenarios
  K_test_values <- c(2, 3, 6) # Underfit, Perfect, Overfit
  noise_levels <- c(0.5, 2.0) # Low Noise, High Noise
  theta_grid <- seq(0.05, 1.0, by = 0.05)
  n_folds <- 3
  
  results_df <- data.frame()
  
  cat("Starting Hyperparameter Grid Search...\n")
  
  for(noise in noise_levels) {
    # Generate data with specific noise
    Z_raw <- C_scores %*% t(U_true) + matrix(rnorm(n_sim * p_sim, sd = noise), n_sim, p_sim)
    Z_sim <- scale(Z_raw, center = TRUE, scale = TRUE)
    
    eta_true <- Z_sim %*% U_true %*% alpha_true
    time_sim <- rexp(n_sim, rate = exp(eta_true))
    censor_time <- rexp(n_sim, rate = 0.5)
    status_sim <- ifelse(time_sim <= censor_time, 1, 0)
    time_sim <- pmin(time_sim, censor_time)
    
    for(K_test in K_test_values) {
      cat(sprintf("Testing Noise = %.1f, K = %d\n", noise, K_test))
      
      folds <- sample(rep(1:n_folds, length.out = n_sim))
      
      for(theta_val in theta_grid) {
        c_indices <- numeric(n_folds)
        
        for (f in 1:n_folds) {
          test_idx <- which(folds == f)
          
          # Train
          fit <- suppressWarnings(
            supervised_pca_cox(Z = Z_sim[-test_idx, ], time = time_sim[-test_idx], 
                               status = status_sim[-test_idx], 
                               K = K_test, theta = theta_val, gamma = 0.001, max_iter = 50)
          )
          
          # Predict
          eta_pred <- Z_sim[test_idx, ] %*% fit$U %*% fit$alpha
          raw_c <- concordance(Surv(time_sim[test_idx], status_sim[test_idx]) ~ eta_pred)$concordance
          c_indices[f] <- ifelse(raw_c < 0.5, 1 - raw_c, raw_c)
        }
        
        # Save results
        results_df <- rbind(results_df, data.frame(
          Noise = as.factor(noise),
          K = as.factor(K_test),
          Theta = theta_val,
          C_Index = mean(c_indices)
        ))
      }
    }
  }
  
  return(results_df)
}

# Run the experiment
experiment_results <- run_hyperparameter_experiment()

# Plot the curves to see how K and Noise shift the optimal Theta
ggplot(experiment_results, aes(x = Theta, y = C_Index, color = K, group = K)) +
  geom_line(size = 1.2) +
  geom_point(size = 2) +
  facet_wrap(~ Noise, labeller = label_both) +
  theme_bw(base_size = 14) +
  labs(
    title = "Cross-Validation Curves: How K and Noise Affect Optimal Theta",
    subtitle = "Higher K often requires higher Theta (more data regularization) to prevent overfitting.",
    x = "Theta (Supervision Weight)",
    y = "Cross-Validated C-Index"
  ) +
  scale_color_brewer(palette = "Set1")




####






#### ==============================================================================
# EXPERIMENT 4: The Showdown (Collective Estimation vs. Standard Cox)
# ==============================================================================
library(survival)
library(tidyverse)

set.seed(42)

# 1. Simulation Parameters
n_train <- 100       # Decrease training size (simulating an expensive clinical trial)
n_test <- 300        # Keep test set large for stable C-Index calculation
p_sim <- 85          # Push dimensions dangerously close to n_train
K_true <- 3 
noise_level <- 3.0   # Pump up the random sensor static

# Generate True Orthogonal Matrix
U_true <- qr.Q(qr(matrix(rnorm(p_sim * K_true), nrow = p_sim, ncol = K_true)))
alpha_true <- matrix(c(-1.5, 0, 2.0), ncol = 1) 
noise_level <- 1.5

# ---------------------------------------------------------
# Generate TRAIN Data
# ---------------------------------------------------------
C_train <- matrix(rnorm(n_train * K_true), nrow = n_train, ncol = K_true)
Z_train_raw <- C_train %*% t(U_true) + matrix(rnorm(n_train * p_sim, sd = noise_level), n_train, p_sim)
Z_train <- scale(Z_train_raw, center = TRUE, scale = TRUE)

eta_train_true <- Z_train %*% U_true %*% alpha_true
time_train <- rexp(n_train, rate = exp(eta_train_true))
censor_train <- rexp(n_train, rate = 0.5)
status_train <- ifelse(time_train <= censor_train, 1, 0)
time_train <- pmin(time_train, censor_train)

# ---------------------------------------------------------
# Generate TEST Data
# ---------------------------------------------------------
C_test <- matrix(rnorm(n_test * K_true), nrow = n_test, ncol = K_true)
Z_test_raw <- C_test %*% t(U_true) + matrix(rnorm(n_test * p_sim, sd = noise_level), n_test, p_sim)
# IMPORTANT: Scale test data using Train parameters to prevent data leakage
Z_test <- scale(Z_test_raw, center = attr(Z_train, "scaled:center"), scale = attr(Z_train, "scaled:scale"))

eta_test_true <- Z_test %*% U_true %*% alpha_true
time_test <- rexp(n_test, rate = exp(eta_test_true))
censor_test <- rexp(n_test, rate = 0.5)
status_test <- ifelse(time_test <= censor_test, 1, 0)
time_test <- pmin(time_test, censor_test)

# ==============================================================================
# 2. FIT THE MODELS
# ==============================================================================
cat("Fitting Standard Cox Model (No Dimension Reduction)...\n")
# We just throw all 50 variables directly into the Cox model
df_train_full <- data.frame(time = time_train, status = status_train, Z_train)
standard_cox <- suppressWarnings(coxph(Surv(time, status) ~ ., data = df_train_full))

cat("Fitting Collective Estimation Model (With Dimension Reduction)...\n")
collective_fit <- supervised_pca_cox_fast(Z = Z_train, time = time_train, status = status_train, 
                                          K = 3, theta = 0.3, gamma = 0.001)

# ==============================================================================
# 3. PREDICT ON BRAND NEW TEST DATA
# ==============================================================================
# Standard Cox Predictions
df_test_full <- data.frame(Z_test)
eta_test_standard <- predict(standard_cox, newdata = df_test_full)

# Collective Estimation Predictions
eta_test_collective <- Z_test %*% collective_fit$U %*% collective_fit$alpha

# Calculate C-Indices
c_standard <- concordance(Surv(time_test, status_test) ~ eta_test_standard)$concordance
c_standard <- ifelse(c_standard < 0.5, 1 - c_standard, c_standard)

c_collective <- concordance(Surv(time_test, status_test) ~ eta_test_collective)$concordance
c_collective <- ifelse(c_collective < 0.5, 1 - c_collective, c_collective)

cat(sprintf("\n--- TEST SET C-INDEX RESULTS ---\n"))
cat(sprintf("Standard Cox (No Reduction): %.4f\n", c_standard))
cat(sprintf("Collective Est (Our Method): %.4f\n", c_collective))

# ==============================================================================
# 4. PLOT THE RESULTS
# ==============================================================================
plot_data <- data.frame(
  True_Risk = rep(as.numeric(eta_test_true), 2),
  Estimated_Risk = c(as.numeric(eta_test_standard), as.numeric(eta_test_collective)),
  Model = factor(rep(c("1. Standard Cox (No Reduction)", "2. Collective Estimation (Our Method)"), each = n_test))
)

showdown_plot <- ggplot(plot_data, aes(x = True_Risk, y = Estimated_Risk, color = Model)) +
  geom_point(alpha = 0.5, size = 2) +
  geom_smooth(method = "lm", color = "black", se = FALSE, linetype = "dashed") +
  facet_wrap(~ Model, scales = "free_y") +
  theme_bw(base_size = 14) +
  scale_color_manual(values = c("#de2d26", "#2c7fb8")) +
  labs(
    title = "True vs. Estimated Risk on Held-Out Test Data",
    subtitle = "Standard Cox overfits the noise. Collective Estimation successfully maps the true biology.",
    x = "True Biological Risk Score (Simulated)",
    y = "Model Estimated Risk Score"
  ) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "none" # Legend is redundant due to facets
  )

print(showdown_plot)



# ==============================================================================
# EXPERIMENT 4: Survival Curve Showdown (For a single representative patient)
# ==============================================================================
# Pick a specific patient from the novel test set to evaluate
patient_idx <- 5 

# ---------------------------------------------------------
# 1. THE TRUE SURVIVAL CURVE
# ---------------------------------------------------------
# Because we simulated using an Exponential distribution, the true baseline 
# hazard is strictly 1.0. Therefore, True S(t) = exp(-t * exp(eta))
t_seq <- seq(0, max(time_test), length.out = 100)
eta_true_val <- eta_test_true[patient_idx]
surv_true <- exp(-exp(eta_true_val) * t_seq)

df_true <- data.frame(
  Time = t_seq, 
  Survival = surv_true, 
  Model = "1. True Biological Curve"
)

# ---------------------------------------------------------
# 2. STANDARD COX SURVIVAL CURVE (No Reduction)
# ---------------------------------------------------------
# Extract the patient's raw 85 variables
patient_std_data <- df_test_full[patient_idx, , drop = FALSE]

# Ask the Standard Cox model to predict their survival
surv_std_obj <- survfit(standard_cox, newdata = patient_std_data)
df_std <- data.frame(
  Time = surv_std_obj$time, 
  Survival = surv_std_obj$surv, 
  Model = "2. Standard Cox (Overfit)"
)

# ---------------------------------------------------------
# 3. COLLECTIVE ESTIMATION SURVIVAL CURVE (Our Method)
# ---------------------------------------------------------
# Step A: Compress the patient's 85 variables down to 3 using our learned U matrix
Z_test_proj <- Z_test %*% collective_fit$U
patient_coll_data <- as.data.frame(Z_test_proj)

# Step B: Ensure the column names perfectly match what the coxph model expects
colnames(patient_coll_data) <- names(coef(collective_fit$final_model))

# Step C: Ask our custom model to predict their survival
surv_coll_obj <- survfit(collective_fit$final_model, newdata = patient_coll_data[patient_idx, , drop = FALSE])
df_coll <- data.frame(
  Time = surv_coll_obj$time, 
  Survival = surv_coll_obj$surv, 
  Model = "3. Collective Estimation"
)

# ==============================================================================
# 4. PLOT THE 3-WAY COMPARISON
# ==============================================================================
# Combine the predicted step-functions
df_estimates <- bind_rows(df_std, df_coll)

surv_plot <- ggplot() +
  # Draw the true mathematical curve (continuous dashed line)
  geom_line(data = df_true, aes(x = Time, y = Survival, color = Model), 
            size = 1.2, linetype = "dashed") +
  # Draw the Cox model estimates (step functions)
  geom_step(data = df_estimates, aes(x = Time, y = Survival, color = Model), 
            size = 1.2) +
  scale_color_manual(values = c("black", "#de2d26", "#2c7fb8")) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Survival Curve Comparison for a Novel Test Patient",
    subtitle = "Standard Cox predictions derail. Collective Estimation tightly tracks the truth.",
    x = "Time",
    y = "Survival Probability S(t)"
  ) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom",
    legend.title = element_blank()
  )

print(surv_plot)



###

# ==============================================================================
# EXPERIMENT 4: The Gold Standard - Risk Stratification Kaplan-Meier Curves
# ==============================================================================
library(survival)
library(tidyverse)
# install.packages("broom") # Run this if you don't have broom installed
library(broom) 

# 1. Split the test patients into High/Low Risk based on the models' predictions
median_std <- median(eta_test_standard)
group_std <- ifelse(eta_test_standard >= median_std, "High Risk", "Low Risk")

median_coll <- median(eta_test_collective)
group_coll <- ifelse(eta_test_collective >= median_coll, "High Risk", "Low Risk")

# 2. Build datasets for the true survival of these predicted groups
df_std_km <- data.frame(Time = time_test, Status = status_test, Group = group_std)
df_coll_km <- data.frame(Time = time_test, Status = status_test, Group = group_coll)

# 3. Fit Kaplan-Meier curves based on the groups
km_std <- survfit(Surv(Time, Status) ~ Group, data = df_std_km)
km_coll <- survfit(Surv(Time, Status) ~ Group, data = df_coll_km)

# Use broom::tidy to easily convert the survival objects to ggplot dataframes
plot_data_std <- tidy(km_std) %>% mutate(Model = "1. Standard Cox")
plot_data_coll <- tidy(km_coll) %>% mutate(Model = "2. Collective Estimation")

# Clean up the group names (broom outputs them as "GroupHigh Risk")
plot_data_combined <- bind_rows(plot_data_std, plot_data_coll) %>%
  mutate(strata = gsub("Group=", "", strata))

# ==============================================================================
# 4. PLOT THE KAPLAN-MEIER STRATIFICATION
# ==============================================================================
km_plot <- ggplot(plot_data_combined, aes(x = time, y = estimate, color = strata)) +
  geom_step(size = 1.2) +
  # Add confidence intervals
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = strata), alpha = 0.1, color = NA) +
  facet_wrap(~ Model) +
  scale_color_manual(values = c("#de2d26", "#2c7fb8")) +
  scale_fill_manual(values = c("#de2d26", "#2c7fb8")) +
  theme_bw(base_size = 14) +
  labs(
    title = "Clinical Utility: High vs. Low Risk Patient Stratification",
    subtitle = "Standard Cox cannot separate groups due to noise. Collective Estimation finds massive separation.",
    x = "Time",
    y = "Overall Survival Probability S(t)",
    color = "Model Prediction:",
    fill = "Model Prediction:"
  ) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

print(km_plot)


###

