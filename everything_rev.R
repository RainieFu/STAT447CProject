library(tmle)
library(SuperLearner)
library(rstan)
library(tibble)
library(posterior)
library(bayesplot)
library(ggplot2)
library(gridExtra)
library(bartMachine)
library(dplyr)

plots_dir <- "bayesian_tmle_plots"
if (!dir.exists(plots_dir)) {
  dir.create(plots_dir)
}

# SuperLearner setup
SL.glm.DCDR <- function(...) SL.glm(...)
SL.gam4.DCDR <- function(...) SL.gam(..., deg.gam=4)
SL.gam6.DCDR <- function(...) SL.gam(..., deg.gam=6)
SL.nnet.DCDR <- function(...) SL.nnet(..., size=4)
SL.mean.DCDR <- function(...) SL.mean(...)
SL.library <- c("SL.glm.DCDR", "SL.gam4.DCDR", "SL.gam6.DCDR", "SL.nnet.DCDR", "SL.mean.DCDR")

# Set seed for reproducibility
set.seed(447)

# Define true ATE from simulation study
TRUE_ATE <- -0.129

# Directly adopted from my previous research work
full_cvtmle_single <- function(data, covars, K = 5) {
  tryCatch({
    cat("Starting full_cvtmle_single\n")
    folds <- sample(rep(1:K, length.out = nrow(data)))

    # Initialize storage for fold estimates and influence curves
    fold_estimates <- numeric(K)
    fold_ey1 <- numeric(K)
    fold_ey0 <- numeric(K)
    all_ic <- numeric(0)

    cat("Number of observations:", nrow(data), "\n")
    cat("Number of covariates:", length(covars), "\n")

    for (k in 1:K) {
      cat("Processing fold", k, "\n")

      # Split data into training and validation sets
      train_idx <- folds != k
      train_data <- data[train_idx, ]
      test_data  <- data[!train_idx, ]

      cat("Train size:", nrow(train_data), "Test size:", nrow(test_data), "\n")

      ## --- Fit the initial estimators on the training data ---

      # Q-model: fit outcome regression using treatment and covariates
      cat("Fitting Q model\n")
      Q_fit <- SuperLearner(
        Y = train_data$Y,
        X = train_data[, c("statin", covars), drop = FALSE],
        family = binomial(),
        SL.library = SL.library
      )
      cat("Q model fit complete\n")

      # g-model: fit propensity score using covariates (excluding treatment)
      cat("Fitting g model\n")
      g_fit <- SuperLearner(
        Y = train_data$statin,
        X = train_data[, covars, drop = FALSE],
        family = binomial(),
        SL.library = SL.library
      )
      cat("g model fit complete\n")

      ## --- Obtain predictions on the validation set ---

      # For g: predict probability of treatment (statin == 1)
      cat("Making g predictions\n")
      g_pred <- predict(g_fit, newdata = test_data[, covars, drop = FALSE])$pred
      g_pred <- pmax(pmin(g_pred, 0.99), 0.01)  # Bound predictions

      # For Q: create two versions of test data (set statin to 0 and 1)
      test_data_A0 <- test_data
      test_data_A1 <- test_data
      test_data_A0$statin <- 0
      test_data_A1$statin <- 1

      cat("Making Q predictions for A=0\n")
      Q0_pred <- predict(Q_fit, newdata = test_data_A0[, c("statin", covars), drop = FALSE])$pred
      cat("Making Q predictions for A=1\n")
      Q1_pred <- predict(Q_fit, newdata = test_data_A1[, c("statin", covars), drop = FALSE])$pred

      # Bound Q predictions if necessary (here assuming probabilities)
      Q0_pred <- pmax(pmin(Q0_pred, 0.99), 0.01)
      Q1_pred <- pmax(pmin(Q1_pred, 0.99), 0.01)
      cat("Q predictions A=0 range:", range(Q0_pred), "\n")
      cat("Q predictions A=1 range:", range(Q1_pred), "\n")

      # Create Q-matrix with two columns: first for A=0, second for A=1.
      Q_matrix <- cbind(Q0_pred, Q1_pred)

      ## --- Run the targeting (TMLE) step on the validation data ---
      cat("Running TMLE on validation data\n")
      tmle_fit <- tmle(
        Y = test_data$Y,
        A = test_data$statin,
        W = test_data[, covars, drop = FALSE],
        Q = Q_matrix,
        g1W = g_pred,
        family = "binomial"
      )
      cat("TMLE complete\n")

      # Extract the TMLE estimates
      fold_estimates[k] <- tmle_fit$estimates$ATE$psi
      fold_ey1[k] <- tmle_fit$estimates$EY1$psi
      fold_ey0[k] <- tmle_fit$estimates$EY0$psi

      cat("Fold estimates - ATE:", fold_estimates[k],
          "EY1:", fold_ey1[k],
          "EY0:", fold_ey0[k], "\n")

      ## --- Compute the efficient influence curve (EIC) ---
      # The literature gives the EIC for ATE as:
      #   D*(O) = (A/g(W) - (1-A)/(1-g(W)))*(Y - Q*(A,W)) + (Q*(1,W) - Q*(0,W)) - ATE
      #
      # Here Q* is the targeted (updated) outcome regression.
      #
      # Retrieve the targeted predictions:
      Qstar <- tmle_fit$Qstar  # Matrix with two columns (A=0 and A=1)
      Q0star <- Qstar[, 1]
      Q1star <- Qstar[, 2]

      # For each observation in the validation set, select the targeted prediction
      Qstar_obs <- ifelse(test_data$statin == 1, Q1star, Q0star)

      # Compute the clever covariate h(W):
      h_a <- test_data$statin / g_pred - (1 - test_data$statin) / (1 - g_pred)

      # Residual term:
      resid <- test_data$Y - Qstar_obs

      # Compute the influence curve for each observation:
      fold_ic <- h_a * resid + (Q1star - Q0star) - tmle_fit$estimates$ATE$psi

      cat("Influence Curve range - ATE:", range(fold_ic), "\n")

      # Append the fold's influence curve
      all_ic <- c(all_ic, fold_ic)
    }

    cat("Computing final estimates\n")
    final_est <- mean(fold_estimates)
    ey1_est <- mean(fold_ey1)
    ey0_est <- mean(fold_ey0)

    # Variance estimate: empirical variance of the EIC divided by the sample size
    var_est <- var(all_ic) / length(all_ic)

    cat("Final estimates - ATE:", final_est,
        "EY1:", ey1_est,
        "EY0:", ey0_est, "\n")

    results <- tibble(
      r1 = ey1_est,
      r0 = ey0_est,
      rd = final_est,
      v1 = var_est,
      v0 = var_est,
      vd = var_est
    )

    cat("Single run complete\n")
    return(results)

  }, error = function(e) {
    cat("Error occurred:", e$message, "\n")
    cat("Error call:", deparse(e$call), "\n")
    warning(paste("Error in full CV-TMLE:", e$message))
    return(tibble(r1 = NA, r0 = NA, rd = NA, v1 = NA, v0 = NA, vd = NA))
  })
}

# Function to run CV-TMLE
run_cvtmle <- function(data, covars, treatment_var, outcome_var, K = 5) {
  # Ensure data has appropriate column names
  data_renamed <- data
  names(data_renamed)[names(data_renamed) == treatment_var] <- "statin"
  names(data_renamed)[names(data_renamed) == outcome_var] <- "Y"

  # Run CV-TMLE
  cv_results <- full_cvtmle_single(
    data = data_renamed,
    covars = covars,
    K = K
  )

  # Extract results and calculate CI
  ate_est <- cv_results$rd
  ate_se <- sqrt(cv_results$vd)
  ate_CI <- c(ate_est - 1.96*ate_se, ate_est + 1.96*ate_se)
  ate_bias <- ate_est - TRUE_ATE
  ate_mse <- ate_bias^2 + ate_se^2
  ate_coverage <- ate_CI[1] <= TRUE_ATE && TRUE_ATE <= ate_CI[2]

  return(list(
    model_type = "cv-tmle",
    ate_mean = ate_est,
    ate_median = ate_est,  # Same as mean for frequentist methods
    ate_mode = ate_est,    # Same as mean for frequentist methods
    ate_se = ate_se,
    ate_bias = ate_bias,
    ate_mse = ate_mse,
    ate_CI = ate_CI,
    ate_coverage = ate_coverage,
    ey1 = cv_results$r1,
    ey0 = cv_results$r0
  ))
}

# Function to get TMLE inputs using SuperLearner for both Q and g
get_tmle_inputs <- function(data, covars, treatment_var, outcome_var) {
  W <- data[, covars, drop = FALSE]
  tmle_fit <- tmle(Y = data[[outcome_var]],
                   A = data[[treatment_var]],
                   W = W,
                   Q.SL.library = SL.library,
                   g.SL.library = SL.library,
                   family = "binomial",
                   cvQinit = FALSE)
  QAW_pred <- tmle_fit$Qstar[cbind(1:nrow(data), data[[treatment_var]] + 1)]
  Q1_pred <- tmle_fit$Qstar[, 2]
  Q0_pred <- tmle_fit$Qstar[, 1]
  gn <- tmle_fit$g$g1W
  H <- (data[[treatment_var]] / gn) - ((1 - data[[treatment_var]]) / (1 - gn))

  return(list(
    QAW = QAW_pred,
    Q1 = Q1_pred,
    Q0 = Q0_pred,
    H = H,
    Y = data[[outcome_var]],
    A = data[[treatment_var]],
    tmle_fit = tmle_fit
  ))
}

# Function to get inputs using BART for Q model and SuperLearner for g model
get_tmle_inputs_bart <- function(data, covars, treatment_var, outcome_var) {
  # Set up data
  W <- data[, covars, drop = FALSE]
  Y <- data[[outcome_var]]
  A <- data[[treatment_var]]

  # Initialize bartMachine
  set_bart_machine_num_cores(2)

  # Fit BART models for Q(1,W) and Q(0,W)
  cat("\nFitting BART for treated group (A=1)...\n")
  bart_1 <- bartMachine(X = W[A == 1, , drop = FALSE],
                        y = Y[A == 1],
                        num_trees = 50,
                        serialize = TRUE,
                        verbose = FALSE)

  cat("Fitting BART for control group (A=0)...\n")
  bart_0 <- bartMachine(X = W[A == 0, , drop = FALSE],
                        y = Y[A == 0],
                        num_trees = 50,
                        serialize = TRUE,
                        verbose = FALSE)

  # Get standard predictions for initial Q values
  Q1 <- predict(bart_1, new_data = W)
  Q0 <- predict(bart_0, new_data = W)
  QAW <- ifelse(A == 1, Q1, Q0)

  # Use SuperLearner for propensity score model g(W)
  sl_g <- SuperLearner(Y = A,
                       X = W,
                       SL.library = SL.library,
                       family = binomial())
  g1W <- predict(sl_g)$pred

  # Bound g to avoid extreme propensity scores
  g1W <- pmin(pmax(g1W, 0.025), 0.975)

  # Construct clever covariate H
  H <- (A / g1W) - ((1 - A) / (1 - g1W))
  H <- as.vector(H)

  # Determine number of posterior samples available
  # Check one posterior sample to determine how many are available
  test_posterior <- bart_machine_get_posterior(bart_1, new_data = W[1:5,])
  n_posterior <- ncol(test_posterior$y_hat_posterior_samples)

  # Use a reasonable number that won't cause memory issues
  n_posterior <- min(n_posterior, 100)

  # Return list with models but without storing all predictions
  return(list(
    QAW = QAW,
    Q1 = Q1,
    Q0 = Q0,
    H = H,
    Y = Y,
    A = A,
    W = W,  # Include W for use in posterior sampling
    g_fit = sl_g,
    bart_1 = bart_1,
    bart_0 = bart_0,
    n_posterior = n_posterior
  ))
}

# Stan model for Bayesian TMLE targeting step
# Integrate stan into R script so I don't need to create separate .stan files for each prior candidates
write_stan_model <- function(prior_sd, file_path = "bayes_tmle.stan") {
  stan_code <- sprintf('
data {
  int<lower=1> N;
  vector<lower=0,upper=1>[N] QAW;
  vector[N] H;
  int<lower=0,upper=1> Y[N];
}
parameters {
  real epsilon;
}
model {
  // Prior
  epsilon ~ normal(0, %.2f);

  // Likelihood
  for (i in 1:N) {
    real p = inv_logit(logit(QAW[i]) + epsilon * H[i]);
    Y[i] ~ bernoulli(p);
  }
}
', prior_sd)
  writeLines(stan_code, con = file_path)
}

# Function to perform the exact invariance test on the Bayesian TMLE
run_exact_invariance_test <- function(tmle_inputs, prior_sd, model_type = "standard", n_samples = 1000, verbose = TRUE) {
  # Extract key inputs
  QAW_bounded <- pmin(pmax(tmle_inputs$QAW, 0.001), 0.999)
  H <- tmle_inputs$H
  Y <- tmle_inputs$Y

  # Create Stan model file for this test
  write_stan_model(prior_sd)

  # Forward simulator for Bayesian TMLE
  forward <- function() {
    # Sample epsilon from prior
    epsilon <- rnorm(1, 0, prior_sd)

    # Calculate probabilities with this epsilon
    logit_QAW <- qlogis(QAW_bounded)
    probs <- plogis(logit_QAW + epsilon * H)

    # Sample Y from the model using these probabilities
    Y_sim <- rbinom(length(probs), 1, probs)

    # Return both the parameter and data
    return(list(
      epsilon = epsilon,
      Y = Y_sim
    ))
  }

  # Joint log-likelihood function
  joint_loglik <- function(epsilon, Y_data) {
    # Prior
    prior_loglik <- dnorm(epsilon, 0, prior_sd, log = TRUE)

    # Likelihood
    logit_QAW <- qlogis(QAW_bounded)
    probs <- plogis(logit_QAW + epsilon * H)
    data_loglik <- sum(dbinom(Y_data, 1, probs, log = TRUE))

    # Joint
    return(prior_loglik + data_loglik)
  }

  # Simple Metropolis-Hastings for Bayesian TMLE
  stationary_mcmc <- function(n_iterations) {
    # Initialize from the forward simulation
    initialization <- forward()
    Y_data <- initialization$Y

    if (n_iterations == 0) {
      return(initialization$epsilon)
    } else {
      current_epsilon <- initialization$epsilon

      for (i in 1:n_iterations) {
        # Propose new epsilon
        proposed_epsilon <- current_epsilon + rnorm(1)

        # Calculate acceptance ratio
        ratio <- exp(joint_loglik(proposed_epsilon, Y_data) -
                       joint_loglik(current_epsilon, Y_data))

        # Accept or reject
        if (runif(1) < ratio) {
          current_epsilon <- proposed_epsilon
        }
      }

      return(current_epsilon)
    }
  }

  # Arrays to store samples
  if (verbose) cat("Running exact invariance test with", n_samples, "samples...\n")

  # Generate forward-only samples
  if (verbose) cat("Generating forward-only samples...\n")
  forward_samples <- replicate(n_samples, stationary_mcmc(0))

  # Generate with-MCMC samples
  if (verbose) cat("Generating with-MCMC samples...\n")
  with_mcmc <- replicate(n_samples, stationary_mcmc(200))

  # Run Kolmogorov-Smirnov test
  ks_result <- ks.test(forward_samples, with_mcmc)

  if (verbose) {
    cat("\nExact Invariance Test Results:\n")
    cat("  KS statistic:", round(ks_result$statistic, 4), "\n")
    cat("  p-value:", format.pval(ks_result$p.value), "\n")

    # Generate diagnostic plots
    test_plot <- ggplot(data.frame(
      Sample = c(rep("Forward", n_samples), rep("MCMC", n_samples)),
      Epsilon = c(forward_samples, with_mcmc)
    ), aes(x = Epsilon, fill = Sample)) +
      geom_density(alpha = 0.5) +
      labs(title = paste("Exact Invariance Test - Prior: N(0,", prior_sd, ")"),
           subtitle = paste("KS p-value =", format.pval(ks_result$p.value)),
           x = "Epsilon", y = "Density") +
      theme_minimal()
    print(test_plot)
  }

  return(list(
    forward_samples = forward_samples,
    mcmc_samples = with_mcmc,
    ks_test = ks_result,
    plot = test_plot
  ))
}

# Function to run Bayesian TMLE with SuperLearner initial estimator
run_bayes_tmle <- function(prior_sd, tmle_inputs, model_type = "standard") {
  # Create the Stan model
  write_stan_model(prior_sd)
  bayes_tmle_model <- stan_model("bayes_tmle.stan")

  # Bound QAW due to error
  QAW_bounded <- pmin(pmax(tmle_inputs$QAW, 0.001), 0.999)

  # Ensure all inputs are proper vectors
  tmle_inputs$H <- as.vector(tmle_inputs$H)
  tmle_inputs$QAW <- as.vector(tmle_inputs$QAW)
  tmle_inputs$Y <- as.vector(tmle_inputs$Y)
  stan_data <- list(
    N = length(tmle_inputs$Y),
    QAW = QAW_bounded,
    H = tmle_inputs$H,
    Y = tmle_inputs$Y
  )

  # Prior predictive check
  eps_prior <- rnorm(1000, 0, prior_sd)
  qstar_sim <- plogis(qlogis(stan_data$QAW) + sample(eps_prior, size=length(tmle_inputs$Y), replace=TRUE) * stan_data$H)
  y_rep_prior <- rbinom(length(qstar_sim), 1, qstar_sim)

  df_priorcheck <- data.frame(
    y_rep = y_rep_prior,
    y_obs = stan_data$Y
  )

  # Summarize prior predictions and observed outcomes for better visualization
  prior_prop_table <- data.frame(
    Source = c("Prior predictive", "Observed data"),
    Mean = c(mean(y_rep_prior), mean(stan_data$Y)),
    SD = c(sd(y_rep_prior), sd(stan_data$Y))
  )

  prior_plot <- ggplot(prior_prop_table, aes(x = Source, y = Mean, fill = Source)) +
    geom_bar(stat = "identity") +
    geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.2) +
    labs(title = paste("Prior Predictive Check — N(0,", prior_sd, ")"),
         y = "Proportion of Y=1", x = "") +
    theme_minimal()
  print(prior_plot)
  ggsave(file.path(plots_dir, paste0("prior_predictive_check_", model_type, "_sd", prior_sd, ".png")),
         prior_plot, width = 8, height = 6)

  # Sample from posterior
  fit <- sampling(
    object = bayes_tmle_model,
    data = stan_data,
    chains = 4,
    iter = 2000,
    warmup = 1000,
    seed = 447,
    refresh = 500
  )

  # MCMC diagnostics
  cat("\nMCMC Diagnostics (Prior SD =", prior_sd, "):\n")
  cat("R-hat for epsilon:", summary(fit)$summary["epsilon", "Rhat"], "\n")
  cat("n_eff for epsilon:", summary(fit)$summary["epsilon", "n_eff"], "\n")

  # Check divergences
  sampler_params <- get_sampler_params(fit, inc_warmup=FALSE)
  divergences <- sum(sapply(sampler_params, function(x) sum(x[, "divergent__"])))
  cat("Divergences:", divergences, "\n")

  # Plot trace of epsilon parameter
  trace_plot <- mcmc_trace(as_draws_df(fit), pars = "epsilon") +
    labs(title = paste("Trace Plot - Prior: N(0,", prior_sd, ")"))
  print(trace_plot)
  ggsave(file.path(plots_dir, paste0("trace_plot_", model_type, "_sd", prior_sd, ".png")),
         trace_plot, width = 10, height = 6)

  # Posterior predictive check
  eps_samples <- extract(fit)$epsilon

  # Creating a matrix where each row is an observation and each column is a posterior sample
  H_mat <- matrix(rep(stan_data$H, length(eps_samples)), ncol=length(eps_samples))
  QAW_mat <- matrix(rep(qlogis(stan_data$QAW), length(eps_samples)), ncol=length(eps_samples))
  prob_mat <- plogis(QAW_mat + H_mat * matrix(eps_samples, nrow=length(stan_data$Y),
                                              ncol=length(eps_samples), byrow=TRUE))

  # Simulate from posterior predictive
  Y_rep <- matrix(rbinom(prod(dim(prob_mat)), 1, prob_mat), nrow=nrow(prob_mat))

  # Calculate mean prediction for each observation
  Y_rep_mean <- rowMeans(Y_rep)

  post_check_df <- data.frame(
    Source = c(rep("Posterior predictive", length(Y_rep_mean)), rep("Observed data", length(stan_data$Y))),
    Value = c(Y_rep_mean, stan_data$Y)
  )

  post_plot <- ggplot(post_check_df, aes(x = Value, fill = Source)) +
    geom_density(alpha = 0.5) +
    labs(title = paste("Posterior Predictive Check - Prior: N(0,", prior_sd, ")"),
         x = "Y (or probability of Y=1)", y = "Density") +
    theme_minimal()
  print(post_plot)
  ggsave(file.path(plots_dir, paste0("posterior_check_", model_type, "_sd", prior_sd, ".png")),
         post_plot, width = 8, height = 6)

  # ATE calculation with BART posterior samples - memory efficient version!
  if (model_type == "bart" && !is.null(tmle_inputs$bart_1) && !is.null(tmle_inputs$bart_0)) {
    # Extract W from tmle_inputs
    W <- tmle_inputs$W

    n_obs <- length(tmle_inputs$Y)
    n_eps_samples <- length(eps_samples)

    # Define batch size (because of memory constraints)
    batch_size <- 50
    n_posterior <- tmle_inputs$n_posterior
    n_batches <- ceiling(n_posterior / batch_size)

    # Initialize storage for combined ATE samples
    all_ate_samples <- numeric(0)

    # Process BART samples in batches
    for (batch in 1:n_batches) {
      start_idx <- (batch - 1) * batch_size + 1
      end_idx <- min(batch * batch_size, n_posterior)
      actual_batch_size <- end_idx - start_idx + 1

      cat(sprintf("Processing BART posterior samples %d to %d...\n", start_idx, end_idx))

      # Get batch of BART posterior samples
      Q1_post <- bart_machine_get_posterior(tmle_inputs$bart_1, new_data = W)
      Q0_post <- bart_machine_get_posterior(tmle_inputs$bart_0, new_data = W)

      # Extract the posterior samples matrix (rows are observations, columns are posterior samples)
      Q1_batch <- Q1_post$y_hat_posterior_samples[, 1:actual_batch_size, drop=FALSE]
      Q0_batch <- Q0_post$y_hat_posterior_samples[, 1:actual_batch_size, drop=FALSE]

      # For each BART sample in this batch
      for (b in 1:actual_batch_size) {
        # Extract the b-th BART posterior sample
        Q1_b <- pmin(pmax(Q1_batch[, b], 0.001), 0.999)
        Q0_b <- pmin(pmax(Q0_batch[, b], 0.001), 0.999)

        # Transform to logit scale
        Q1_logit_b <- qlogis(Q1_b)
        Q0_logit_b <- qlogis(Q0_b)

        # Apply targeting for each epsilon sample and calculate ATE
        batch_ate_samples <- numeric(n_eps_samples)

        for (e in 1:n_eps_samples) {
          Q1_star_b_e <- plogis(Q1_logit_b + eps_samples[e])
          Q0_star_b_e <- plogis(Q0_logit_b + eps_samples[e])
          batch_ate_samples[e] <- mean(Q1_star_b_e) - mean(Q0_star_b_e)
        }

        # Append to all ATE samples
        all_ate_samples <- c(all_ate_samples, batch_ate_samples)
      }

      # Clean up to free memory
      rm(Q1_batch, Q0_batch, Q1_post, Q0_post)
      gc()
    }

    # Use the combined samples
    ate_samples <- all_ate_samples

  } else {
    # Original code for non-BART methods
    # Prepare Q1 and Q0 probabilities
    Q1_bounded <- pmin(pmax(tmle_inputs$Q1, 0.001), 0.999)
    Q0_bounded <- pmin(pmax(tmle_inputs$Q0, 0.001), 0.999)

    # Transform to logit scale
    Q1_logit <- qlogis(Q1_bounded)
    Q0_logit <- qlogis(Q0_bounded)

    # Create matrices for efficient calculation
    Q1_star_logit <- outer(Q1_logit, eps_samples, "+")
    Q0_star_logit <- outer(Q0_logit, eps_samples, "+")

    # Transform back to probability scale
    Q1_star <- plogis(Q1_star_logit)
    Q0_star <- plogis(Q0_star_logit)

    # Calculate ATE samples
    ate_samples <- colMeans(Q1_star) - colMeans(Q0_star)
  }

  # Calculate evaluation metrics
  ate_mean <- mean(ate_samples)
  ate_median <- median(ate_samples)
  ate_mode <- density(ate_samples)$x[which.max(density(ate_samples)$y)]
  ate_se <- sd(ate_samples)
  ate_bias <- ate_mean - TRUE_ATE
  ate_mse <- mean((ate_samples - TRUE_ATE)^2)
  ate_CI <- quantile(ate_samples, c(0.025, 0.975))
  ate_coverage <- ate_CI[1] <= TRUE_ATE && TRUE_ATE <= ate_CI[2]

  # Calculate TMLE metrics for comparison - different handling for standard vs BART
  if (model_type == "standard" && !is.null(tmle_inputs$tmle_fit)) {
    tmle_ate <- tmle_inputs$tmle_fit$estimates$ATE$psi
    tmle_se <- sqrt(tmle_inputs$tmle_fit$estimates$ATE$var.psi)
    tmle_CI <- c(tmle_ate - 1.96*tmle_se, tmle_ate + 1.96*tmle_se)
  } else {
    # For BART, calculate naive ATE estimate
    tmle_ate <- mean(tmle_inputs$Q1 - tmle_inputs$Q0)
    tmle_se <- NA
    tmle_CI <- c(NA, NA)
  }

  tmle_bias <- tmle_ate - TRUE_ATE
  tmle_mse <- tmle_bias^2 + (if(is.na(tmle_se)) 0 else tmle_se^2)
  tmle_coverage <- if(any(is.na(tmle_CI))) NA else (tmle_CI[1] <= TRUE_ATE && TRUE_ATE <= tmle_CI[2])

  # Plot ATE posterior with true value and TMLE estimate
  ate_plot <- ggplot(data.frame(ate = ate_samples), aes(x = ate)) +
    geom_density(fill = "orange", alpha = 0.5) +
    geom_vline(xintercept = TRUE_ATE, color = "red", linetype = "dashed", size = 1) +
    geom_vline(xintercept = tmle_ate, color = "blue", linetype = "dashed", size = 1) +
    geom_vline(xintercept = ate_mean, color = "green3", linetype = "solid", size = 1) +
    labs(title = paste("ATE Posterior - Prior: N(0,", prior_sd, ")"),
         subtitle = paste("True ATE (red) =", round(TRUE_ATE, 4),
                          "| TMLE est. (blue) =", round(tmle_ate, 4),
                          "| Bayes est. (green) =", round(ate_mean, 4)),
         x = "Average Treatment Effect", y = "Density") +
    theme_minimal()
  print(ate_plot)
  ggsave(file.path(plots_dir, paste0("ate_posterior_", model_type, "_sd", prior_sd, ".png")),
         ate_plot, width = 8, height = 6)

  return(list(
    prior_sd = prior_sd,
    model_type = model_type,
    fit = fit,
    ate_mean = ate_mean,
    ate_median = ate_median,
    ate_mode = ate_mode,
    ate_se = ate_se,
    ate_bias = ate_bias,
    ate_mse = ate_mse,
    ate_CI = ate_CI,
    ate_coverage = ate_coverage,
    tmle_ate = tmle_ate,
    tmle_se = tmle_se,
    tmle_bias = tmle_bias,
    tmle_mse = tmle_mse,
    tmle_CI = tmle_CI,
    tmle_coverage = tmle_coverage,
    eps_samples = eps_samples,
    prior_plot = prior_plot,
    post_plot = post_plot,
    ate_plot = ate_plot,
    trace_plot = trace_plot
  ))
}

# === EXECUTE SIMULATION STUDY ===
# Load simulation data
data <- read.csv("/Users/rainiefu/Downloads/STAT447CProject/Data/statin_sim_data_full_rep.csv")
covars <- c("age", "ldl_log", "diabetes", "risk_score")

# === TEST MCMC IMPLEMENTATION USING EXACT INVARIANCE TEST ===
cat("\n\n=========== TESTING MCMC IMPLEMENTATION WITH EXACT INVARIANCE TEST ===========\n")

# Test SuperLearner TMLE implementation with each prior
invariance_results_sl <- lapply(prior_sds, function(sd) {
  cat("\nTesting SuperLearner Bayesian TMLE with prior SD =", sd, "\n")
  run_exact_invariance_test(tmle_inputs_sl, sd, n_samples = 500)
})

# Test BART TMLE implementation with each prior
invariance_results_bart <- lapply(prior_sds, function(sd) {
  cat("\nTesting BART Bayesian TMLE with prior SD =", sd, "\n")
  run_exact_invariance_test(tmle_inputs_bart, sd, n_samples = 500)
})

# Save invariance test plots
for (i in 1:length(prior_sds)) {
  ggsave(file.path(plots_dir, paste0("invariance_test_sl_sd", prior_sds[i], ".png")),
         invariance_results_sl[[i]]$plot, width = 8, height = 6)
  ggsave(file.path(plots_dir, paste0("invariance_test_bart_sd", prior_sds[i], ".png")),
         invariance_results_bart[[i]]$plot, width = 8, height = 6)
}

# === APPROACH 1: SuperLearner for Q initial, Bayesian for targeting ===
cat("\n\n=========== APPROACH 1: SuperLearner + Bayesian targeting ===========\n")
tmle_inputs_sl <- get_tmle_inputs(
  data = data,
  covars = covars,
  treatment_var = "statin",
  outcome_var = "Y"
)

prior_sds <- c(1, 2, 5)
sl_results <- lapply(prior_sds, function(sd) {
  run_bayes_tmle(sd, tmle_inputs_sl, model_type = "standard")
})

# === APPROACH 2: BART for Q initial, Bayesian for targeting ===
cat("\n\n=========== APPROACH 2: BART + Bayesian targeting ===========\n")
tmle_inputs_bart <- get_tmle_inputs_bart(
  data = data,
  covars = covars,
  treatment_var = "statin",
  outcome_var = "Y"
)

bart_results <- lapply(prior_sds, function(sd) {
  run_bayes_tmle(sd, tmle_inputs_bart, model_type = "bart")
})

# === APPROACH 3: CV-TMLE ===
cat("\n\n=========== APPROACH 3: Cross-Validated TMLE ===========\n")
cvtmle_result <- run_cvtmle(
  data = data,
  covars = covars,
  treatment_var = "statin",
  outcome_var = "Y",
  K = 5  # 5-fold cross-validation
)
# === COMPILE RESULTS FOR COMPARISON ===
# Create a comprehensive evaluation metrics table
create_eval_table <- function(results_list, method_prefix) {
  do.call(rbind, lapply(results_list, function(res) {
    data.frame(
      Method = paste0(method_prefix, " (Prior SD=", res$prior_sd, ")"),
      Mean_Est = round(res$ate_mean, 4),
      Median_Est = round(res$ate_median, 4),
      Mode_Est = round(res$ate_mode, 4),
      Std_Error = round(res$ate_se, 4),
      Bias = round(res$ate_bias, 4),
      MSE = round(res$ate_mse, 5),
      CI_Lower = round(res$ate_CI[1], 4),
      CI_Upper = round(res$ate_CI[2], 4),
      CI_Width = round(res$ate_CI[2] - res$ate_CI[1], 4),
      Coverage = res$ate_coverage
    )
  }))
}

sl_eval <- create_eval_table(sl_results, "SuperLearner + Bayes")
bart_eval <- create_eval_table(bart_results, "BART + Bayes")

# Add TMLE results to the table
tmle_row <- data.frame(
  Method = "Vanilla-TMLE (frequentist)",
  Mean_Est = round(sl_results[[1]]$tmle_ate, 4),
  Median_Est = round(sl_results[[1]]$tmle_ate, 4),
  Mode_Est = round(sl_results[[1]]$tmle_ate, 4),
  Std_Error = round(sl_results[[1]]$tmle_se, 4),
  Bias = round(sl_results[[1]]$tmle_bias, 4),
  MSE = round(sl_results[[1]]$tmle_mse, 5),
  CI_Lower = round(sl_results[[1]]$tmle_CI[1], 4),
  CI_Upper = round(sl_results[[1]]$tmle_CI[2], 4),
  CI_Width = round(sl_results[[1]]$tmle_CI[2] - sl_results[[1]]$tmle_CI[1], 4),
  Coverage = sl_results[[1]]$tmle_coverage
)

# Add CV-TMLE to evaluation table
cvtmle_row <- data.frame(
  Method = "CV-TMLE (frequentist)",
  Mean_Est = round(cvtmle_result$ate_mean, 4),
  Median_Est = round(cvtmle_result$ate_median, 4),
  Mode_Est = round(cvtmle_result$ate_mode, 4),
  Std_Error = round(cvtmle_result$ate_se, 4),
  Bias = round(cvtmle_result$ate_bias, 4),
  MSE = round(cvtmle_result$ate_mse, 5),
  CI_Lower = round(cvtmle_result$ate_CI[1], 4),
  CI_Upper = round(cvtmle_result$ate_CI[2], 4),
  CI_Width = round(cvtmle_result$ate_CI[2] - cvtmle_result$ate_CI[1], 4),
  Coverage = cvtmle_result$ate_coverage
)

evaluation_table <- rbind(sl_eval, bart_eval, tmle_row, cvtmle_row)

# Print the evaluation metrics table
cat("\n===== EVALUATION METRICS =====\n")
cat("TRUE ATE =", TRUE_ATE, "\n\n")
print(evaluation_table)

# Create a visualization comparing all methods
comparison_df <- data.frame(
  Method = c(
    paste0("SuperLearner + Bayes (SD=", prior_sds, ")"),
    paste0("BART + Bayes (SD=", prior_sds, ")"),
    "Vanilla-TMLE (frequentist)",
    "CV-TMLE (frequentist)"
  ),
  Estimate = c(
    sapply(sl_results, function(r) r$ate_mean),
    sapply(bart_results, function(r) r$ate_mean),
    sl_results[[1]]$tmle_ate,
    cvtmle_result$ate_mean
  ),
  Lower = c(
    sapply(sl_results, function(r) r$ate_CI[1]),
    sapply(bart_results, function(r) r$ate_CI[1]),
    sl_results[[1]]$tmle_CI[1],
    cvtmle_result$ate_CI[1]
  ),
  Upper = c(
    sapply(sl_results, function(r) r$ate_CI[2]),
    sapply(bart_results, function(r) r$ate_CI[2]),
    sl_results[[1]]$tmle_CI[2],
    cvtmle_result$ate_CI[2]
  )
)

compare_plot <- ggplot(comparison_df, aes(x = reorder(Method, Estimate), y = Estimate, color = Method)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2) +
  geom_hline(yintercept = TRUE_ATE, linetype = "dashed", color = "red") +
  coord_flip() +
  labs(title = "Comparison of ATE Estimates",
       subtitle = paste("True ATE =", TRUE_ATE, "(red dashed line)"),
       y = "Average Treatment Effect", x = "") +
  theme_minimal() +
  theme(legend.position = "none")
print(compare_plot)

# Prior sensitivity analysis
sensitivity_df <- data.frame(
  prior_sd = rep(prior_sds, 2),
  method = rep(c("SuperLearner", "BART"), each = length(prior_sds)),
  ate_mean = c(sapply(sl_results, function(r) r$ate_mean),
               sapply(bart_results, function(r) r$ate_mean)),
  ate_se = c(sapply(sl_results, function(r) r$ate_se),
             sapply(bart_results, function(r) r$ate_se)),
  mse = c(sapply(sl_results, function(r) r$ate_mse),
          sapply(bart_results, function(r) r$ate_mse)),
  coverage = c(sapply(sl_results, function(r) as.numeric(r$ate_coverage)),
               sapply(bart_results, function(r) as.numeric(r$ate_coverage)))
)

# Prior sensitivity visualization
sensitivity_plot <- ggplot(sensitivity_df, aes(x = prior_sd, y = ate_mean, color = method, group = method)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  geom_ribbon(aes(ymin = ate_mean - ate_se, ymax = ate_mean + ate_se, fill = method), alpha = 0.2) +
  geom_hline(yintercept = TRUE_ATE, linetype = "dashed", color = "red") +
  labs(title = "Prior Sensitivity Analysis",
       subtitle = paste("True ATE =", TRUE_ATE, "(red dashed line)"),
       x = "Prior SD", y = "ATE Estimate") +
  theme_minimal()
print(sensitivity_plot)

# MSE comparison
mse_plot <- ggplot(sensitivity_df, aes(x = prior_sd, y = mse, color = method, group = method)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  labs(title = "MSE by Prior Standard Deviation",
       x = "Prior SD", y = "Mean Squared Error") +
  theme_minimal()
print(mse_plot)

# Create distributions of epsilon (fluctuation parameter) for comparing priors
eps_df <- do.call(rbind, lapply(1:length(prior_sds), function(i) {
  sl_eps <- data.frame(
    prior_sd = prior_sds[i],
    method = "SuperLearner",
    epsilon = sl_results[[i]]$eps_samples
  )
  bart_eps <- data.frame(
    prior_sd = prior_sds[i],
    method = "BART",
    epsilon = bart_results[[i]]$eps_samples
  )
  rbind(sl_eps, bart_eps)
}))

eps_plot <- ggplot(eps_df, aes(x = epsilon, fill = as.factor(prior_sd))) +
  geom_density(alpha = 0.5) +
  facet_grid(method ~ .) +
  labs(title = "Posterior Distributions of Fluctuation Parameter (ε)",
       x = "Epsilon Value", y = "Density",
       fill = "Prior SD") +
  theme_minimal()
print(eps_plot)

# Final summary: combine into a single plot
grid.arrange(
  compare_plot,
  sensitivity_plot,
  mse_plot,
  eps_plot,
  ncol = 2,
  top = "Bayesian TMLE Study: Key Results"
)
