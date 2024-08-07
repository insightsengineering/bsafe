# This script is meant to create the thresholds for the B-Safe App, that will
# be used for testing, the intention is to use only RBesT and none of the other
# features of bsafe


# Setup the code ----------------------------------------------------------


# Set your work directory to your file location, all other path go out from there
setwd("~/bsafe/tests")
source("testthat/thresholds/setup_proportions.R")
source("testthat/thresholds/setup_rates.R")


# Load necessary library
library(purrr)

# Define the path where the CSV files are stored
path <- "testthat/scenarios"

# Get a list of all CSV files in the directory
files_scenarios <- list.files(path, pattern = "*.csv")

# Read all CSV files into a list
scenarios_list <- files_scenarios %>%
  purrr::map(~ read.csv(file.path(path, .)))

# Number of Scenarios
n_scenarios <- length(scenarios_list)

# Numbers for rounding
deci <- 4

# For testing purpsoe and to set to the same settings as in B-Safe
testing <- FALSE
if (testing == FALSE) {
  # for exact method look stan/jags up, but if no rounding etc n_mcmc = (n_iter-n_warm)*n_chain/n_thin
  n_iter <- 6000
  n_burn <- 2000
  n_chains <- 4
  n_thin <- 4

  # Number of Runs for the thresholds
  n_runs <- 1000
} else {
  n_iter <- 1000
  n_burn <- 100
  n_chains <- 2
  n_thin <- 1

  # Number of Runs for the thresholds
  n_runs <- 10
}

# boundaries for credible intervals
crilb <- 0.025
criub <- 0.975


# Thresholds - Setup Arrays  ------------------------------------------------
# +1 for the index, because the first index is the any scenario

th_prop_final <- testing_list_props[[1 + 1]]$treshholds
th_prop_NA <- th_prop_final[, 1:6]
th_prop_NA[, ] <- NA
array_prop_th <- array(unlist(replicate(n_runs, th_prop_NA, simplify = FALSE)),
  dim = c(nrow(th_prop_NA), ncol(th_prop_NA), n_runs)
)

dimnames(array_prop_th) <- list(
  rownames(th_prop_NA),
  c("mean", "sd", "med", "cri025", "cri975", "ess"),
  paste0("Run ", 1:n_runs)
)


th_rate_final <- testing_list_rates[[1 + 1]]$treshholds
th_rate_NA <- th_rate_final[, 1:6]
th_rate_NA[, ] <- NA
array_rate_th <- array(unlist(replicate(n_runs, th_prop_NA, simplify = FALSE)),
  dim = c(nrow(th_rate_NA), ncol(th_rate_NA), n_runs)
)

dimnames(array_rate_th) <- list(
  rownames(th_rate_NA),
  c("mean", "sd", "med", "cri025", "cri975", "ess"),
  paste0("Run ", 1:n_runs)
)


# each scenario the values will now be calculated:
# For loop scen ----------------------------------------------------------------
# still issue with scen 4
for (scen in c(1:11)[-4]) {
  print(paste0("Scenario: ", scen, " begun at ", Sys.time()))

  if (scen == 1) {
    starttime <- Sys.time()
  }

  # read in the data
  raw_data <- scenarios_list[[scen]]

  # prepare the data
  # filter for historic information and correct group
  input_data <-
    raw_data %>% dplyr::filter(HIST == 1 &
      ARM == as.character(testing_list_props[[scen + 1]]$parameters["group"]))

  # define seed for seeds
  seed <- as.numeric(testing_list_rates[[scen + 1]]$parameters["seed"])

  # Tau, as proposed by Neuenschwander
  hist_borrow <- as.character(testing_list_props[[scen + 1]]$parameters["heterog"])
  if (hist_borrow == "Small") {
    tau <- 0.125
  } else if (hist_borrow == "Moderate") {
    tau <- 0.25
  } else if (hist_borrow == "Substantial") {
    tau <- 0.5
  } else if (hist_borrow == "Large") {
    tau <- 1
  } else if (hist_borrow == "Very Large") {
    tau <- 2
  }

  # in between heterogeneity distirbution
  tau_dist <- as.character(testing_list_props[[scen + 1]]$parameters["tau"])

  # define seeds
  set.seed(seed)
  seeds <- sample.int(9999999, size = n_runs, replace = FALSE)

  for (run in 1:n_runs) {
    if (run %% 100 == 0) {
      print(paste("Iteration:", run, "Time:", Sys.time(), "Scenario:", scen))
    }

    # MAP prior proportions ---------------------------------------------------

    set.seed(seeds[run])
    map_prop <- RBesT::gMAP(
      cbind(N_WITH_AE, N - N_WITH_AE) ~ 1 | STUDYID,
      data = input_data,
      tau.dist = tau_dist, # HalfNormal
      tau.prior = tau,
      beta.prior = 2,
      iter = getOption("RBesT.MC.iter", n_iter),
      warmup = getOption("RBesT.MC.warmup", n_burn),
      thin = getOption("RBesT.MC.thin", n_thin),
      chains = getOption("RBesT.MC.chains", n_chains),
      family = "binomial"
    )

    sample_prop <- rstan::extract(map_prop$fit)$theta_resp_pred

    array_prop_th["MAP Prior", "mean", run] <- mean(sample_prop)
    array_prop_th["MAP Prior", "sd", run] <- sd(sample_prop)
    array_prop_th["MAP Prior", "med", run] <- median(sample_prop)
    array_prop_th["MAP Prior", "cri025", run] <- quantile(sample_prop, probs = crilb)
    array_prop_th["MAP Prior", "cri975", run] <- quantile(sample_prop, probs = criub)

    fit_prop <- RBesT::automixfit(map_prop, Nc = seq(3, 3))

    # FIXME: Scen 4 fix me, infinite size, with testing everything is fine, with proper we get bad results

    array_prop_th["MAP Prior", "ess", run] <- RBesT::ess(fit_prop, method = "elir")


    # Robust MAP proportions --------------------------------------------------

    rob_map_prop <- RBesT::robustify(fit_prop, weight = as.numeric(testing_list_props[[scen + 1]]$parameters["rob_weight"]))


    ana_obj <- summary(rob_map_prop,
      quantile.type = 7, digits = deci, na.rm = TRUE,
      names = TRUE, type = c("quartiles"),
      probs = c(crilb, 0.5, criub)
    )

    array_prop_th["Robustified", "mean", run] <- ana_obj["mean"]
    array_prop_th["Robustified", "sd", run] <- ana_obj["sd"]
    array_prop_th["Robustified", "med", run] <- ana_obj["50.0%"]
    array_prop_th["Robustified", "cri025", run] <- ana_obj[3]
    array_prop_th["Robustified", "cri975", run] <- ana_obj[5]
    array_prop_th["Robustified", "ess", run] <- RBesT::ess(rob_map_prop)


    # Likelihood proportions --------------------------------------------------

    new_r <- testing_list_props[[scen + 1]]$parameters["nta_event"]
    new_n <- testing_list_props[[scen + 1]]$parameters["nta_npat"]


    # Only for the display of the likelihood, there occur errors for a,b <1 or a=b
    if (new_r == 0) {
      alpha <- new_r + 1
      beta <- new_n - (new_r + 1)
    } else if (new_r == new_n) {
      alpha <- new_r - 1
      beta <- new_n - (new_r - 1)
    } else {
      alpha <- new_r
      beta <- new_n - new_r
    }


    lik_prop <- RBesT::mixbeta(lik = c(1, as.numeric(alpha), as.numeric(beta)))

    ana_obj <- summary(lik_prop,
      quantile.type = 7, digits = deci, na.rm = TRUE,
      names = TRUE, type = c("quartiles"),
      probs = c(crilb, 0.5, criub)
    )

    array_prop_th["Likelihood", "mean", run] <- ana_obj["mean"]
    array_prop_th["Likelihood", "sd", run] <- ana_obj["sd"]
    array_prop_th["Likelihood", "med", run] <- ana_obj["50.0%"]
    array_prop_th["Likelihood", "cri025", run] <- ana_obj[3]
    array_prop_th["Likelihood", "cri975", run] <- ana_obj[5]


    # Posterior proportions ---------------------------------------------------

    post_prop <- RBesT::postmix(rob_map_prop,
      n = as.numeric(new_n),
      r = as.numeric(new_r)
    )


    ana_obj <- summary(post_prop,
      quantile.type = 7, digits = deci, na.rm = TRUE,
      names = TRUE, type = c("quartiles"),
      probs = c(crilb, 0.5, criub)
    )

    array_prop_th["Posterior", "mean", run] <- ana_obj["mean"]
    array_prop_th["Posterior", "sd", run] <- ana_obj["sd"]
    array_prop_th["Posterior", "med", run] <- ana_obj["50.0%"]
    array_prop_th["Posterior", "cri025", run] <- ana_obj[3]
    array_prop_th["Posterior", "cri975", run] <- ana_obj[5]

    # Potential for parralelization
    # MAP prior rates ---------------------------------------------------------

    hist_borrow <- as.character(testing_list_rates[[scen + 1]]$parameters["heterog"])

    if (hist_borrow == "Small") {
      tau <- 0.0625
    } else if (hist_borrow == "Moderate") {
      tau <- 0.125
    } else if (hist_borrow == "Substantial") {
      tau <- 0.25
    } else if (hist_borrow == "Large") {
      tau <- 0.5
    } else if (hist_borrow == "Very Large") {
      tau <- 1
    }

    set.seed(seeds[run])
    map_rate <- RBesT::gMAP(N_WITH_AE ~ 1 + offset(log(TOT_EXP)) | STUDYID,
      data = input_data,
      tau.dist = tau_dist,
      tau.prior = tau,
      beta.prior = 1,
      iter = getOption("RBesT.MC.iter", n_iter),
      warmup = getOption("RBesT.MC.warmup", n_burn),
      thin = getOption("RBesT.MC.thin", n_thin),
      chains = getOption("RBesT.MC.chains", n_chains),
      family = "poisson"
    )

    sample_rate_log <- rstan::extract(map_rate$fit)$theta_pred
    sample_rate_exp <- rstan::extract(map_rate$fit)$theta_resp_pred

    array_rate_th["log MAP Prior", "mean", run] <- mean(sample_rate_log)
    array_rate_th["log MAP Prior", "sd", run] <- sd(sample_rate_log)
    array_rate_th["log MAP Prior", "med", run] <- median(sample_rate_log)
    array_rate_th["log MAP Prior", "cri025", run] <- quantile(sample_rate_log, probs = crilb)
    array_rate_th["log MAP Prior", "cri975", run] <- quantile(sample_rate_log, probs = criub)

    array_rate_th["exp MAP Prior", "mean", run] <- mean(sample_rate_exp)
    array_rate_th["exp MAP Prior", "sd", run] <- sd(sample_rate_exp)
    array_rate_th["exp MAP Prior", "med", run] <- median(sample_rate_exp)
    array_rate_th["exp MAP Prior", "cri025", run] <- quantile(sample_rate_exp, probs = crilb)
    array_rate_th["exp MAP Prior", "cri975", run] <- quantile(sample_rate_exp, probs = criub)


    theta_fit <- as.data.frame(map_rate$fit)
    theta_pred <- theta_fit$theta_pred
    fit_rate <- RBesT::automixfit(theta_pred, Nc = seq(3, 3))

    array_rate_th["log MAP Prior", "ess", run] <- RBesT::ess(fit_rate, method = "elir", sigma = 1)


    # Robust MAP rates -------------------------------------------------------------

    rob_map_rate <- RBesT::robustify(fit_rate,
      weight = as.numeric(testing_list_rates[[scen + 1]]$parameters["rob_weight"]),
      mean = log(as.numeric(testing_list_rates[[scen + 1]]$parameters["rob_mean"])), sigma = 1
    )

    set.seed(seeds[run])
    rob_rate_exp <- exp(RBesT::rmix(rob_map_rate, 100000))
    rob_rate_log <- summary(rob_map_rate,
      quantile.type = 7, digits = deci, na.rm = TRUE,
      names = TRUE, type = c("quartiles"),
      probs = c(crilb, 0.5, criub)
    )

    array_rate_th["log Robustified", "mean", run] <- rob_rate_log["mean"]
    array_rate_th["log Robustified", "sd", run] <- rob_rate_log["sd"]
    array_rate_th["log Robustified", "med", run] <- rob_rate_log["50.0%"]
    array_rate_th["log Robustified", "cri025", run] <- rob_rate_log[3]
    array_rate_th["log Robustified", "cri975", run] <- rob_rate_log[5]
    array_rate_th["log Robustified", "ess", run] <- RBesT::ess(rob_map_rate, method = "elir", sigma = 1)

    array_rate_th["exp Robustified", "mean", run] <- mean(rob_rate_exp)
    array_rate_th["exp Robustified", "sd", run] <- sd(rob_rate_exp)
    array_rate_th["exp Robustified", "med", run] <- median(rob_rate_exp)
    array_rate_th["exp Robustified", "cri025", run] <- quantile(rob_rate_exp, probs = crilb)
    array_rate_th["exp Robustified", "cri975", run] <- quantile(rob_rate_exp, probs = criub)



    # Likelihood rates --------------------------------------------------------

    nta_events <- testing_list_rates[[scen + 1]]$parameters["nta_event"]
    nta_time <- testing_list_rates[[scen + 1]]$parameters["nta_time"]

    new_hzrate <- as.numeric(nta_events / nta_time)
    new_se <- as.numeric(sqrt(1 / nta_events))


    lik_rate <- RBesT::mixnorm(lik = c(1, log(new_hzrate), sigma = new_se))

    lik_ana <- summary(lik_rate,
      quantile.type = 7, digits = deci, na.rm = TRUE,
      names = TRUE, type = c("quartiles"),
      probs = c(crilb, 0.5, criub)
    )


    array_rate_th["log Likelihood", "mean", run] <- lik_ana["mean"]
    array_rate_th["log Likelihood", "sd", run] <- lik_ana["sd"]
    array_rate_th["log Likelihood", "med", run] <- lik_ana["50.0%"]
    array_rate_th["log Likelihood", "cri025", run] <- lik_ana[3]
    array_rate_th["log Likelihood", "cri975", run] <- lik_ana[5]

    # Posterior rates

    post_rate <- RBesT::postmix(rob_map_rate, m = log(new_hzrate), se = new_se)


    post_rate_log <- summary(post_rate,
      quantile.type = 7, digits = deci, na.rm = TRUE,
      names = TRUE, type = c("quartiles"),
      probs = c(crilb, 0.5, criub)
    )

    array_rate_th["log Posterior", "mean", run] <- post_rate_log["mean"]
    array_rate_th["log Posterior", "sd", run] <- post_rate_log["sd"]
    array_rate_th["log Posterior", "med", run] <- post_rate_log["50.0%"]
    array_rate_th["log Posterior", "cri025", run] <- post_rate_log[3]
    array_rate_th["log Posterior", "cri975", run] <- post_rate_log[5]

    set.seed(seeds[run])
    post_rate_exp <- exp(RBesT::rmix(post_rate, 100000))

    array_rate_th["exp Posterior", "mean", run] <- mean(post_rate_exp)
    array_rate_th["exp Posterior", "sd", run] <- sd(post_rate_exp)
    array_rate_th["exp Posterior", "med", run] <- median(post_rate_exp)
    array_rate_th["exp Posterior", "cri025", run] <- quantile(post_rate_exp, probs = crilb)
    array_rate_th["exp Posterior", "cri975", run] <- quantile(post_rate_exp, probs = criub)

    if (scen == n_scenarios && run == n_runs) {
      endtime <- Sys.time()
    }
  }

  # Building the treshholds
  th_rate_final[, c("mean_lb", "sd_lb", "median_lb", "cri_lb_lb", "cri_ub_lb", "ess_lb")] <- apply(
    array_rate_th, c(1, 2), function(x) quantile(x, probs = c(0.01), na.rm = TRUE)
  )
  th_rate_final[, c("mean_min", "sd_min", "median_min", "cri_lb_min", "cri_ub_min", "ess_min")] <- apply(
    array_rate_th, c(1, 2), function(x) min(x, na.rm = TRUE)
  )

  th_rate_final[, c("mean_ub", "sd_ub", "median_ub", "cri_lb_ub", "cri_ub_ub", "ess_ub")] <- apply(
    array_rate_th, c(1, 2), function(x) quantile(x, probs = c(0.99), na.rm = TRUE)
  )
  th_rate_final[, c("mean_max", "sd_max", "median_max", "cri_lb_max", "cri_ub_max", "ess_max")] <- apply(
    array_rate_th, c(1, 2), function(x) max(x, na.rm = TRUE)
  )

  th_prop_final[, c("mean_lb", "sd_lb", "median_lb", "cri_lb_lb", "cri_ub_lb", "ess_lb")] <- apply(
    array_prop_th, c(1, 2), function(x) quantile(x, probs = c(0.01), na.rm = TRUE)
  )
  th_prop_final[, c("mean_min", "sd_min", "median_min", "cri_lb_min", "cri_ub_min", "ess_min")] <- apply(
    array_prop_th, c(1, 2), function(x) min(x, na.rm = TRUE)
  )

  th_prop_final[, c("mean_ub", "sd_ub", "median_ub", "cri_lb_ub", "cri_ub_ub", "ess_ub")] <- apply(
    array_prop_th, c(1, 2), function(x) quantile(x, probs = c(0.99), na.rm = TRUE)
  )
  th_prop_final[, c("mean_max", "sd_max", "median_max", "cri_lb_max", "cri_ub_max", "ess_max")] <- apply(
    array_prop_th, c(1, 2), function(x) max(x, na.rm = TRUE)
  )

  th_rate_final[th_rate_final == Inf | th_rate_final == -Inf] <- NA
  th_prop_final[th_prop_final == Inf | th_prop_final == -Inf] <- NA

  # Be careful not to override
  if (scen < 10) {
    scen_chr <- paste0(0, scen)
  } else {
    scen_chr <- scen
  }

  saveRDS(array_rate_th, paste0("th_scen/th_rate_scen", scen_chr, ".rds"))
  saveRDS(array_prop_th, paste0("th_scen/th_prop_scen", scen_chr, ".rds"))

  testing_list_props[[scen + 1]]$treshholds <- th_prop_final
  testing_list_rates[[scen + 1]]$treshholds <- th_rate_final
}

# Be careful not to override
# setwd("~/bsafe/tests/testthat/thresholds")
# saveRDS(testing_list_rates, "testing_list_rates.rds")
# saveRDS(testing_list_props, "testing_list_props.rds")

# [1] "Scenario: 1 begun at 2024-07-29 15:10:59.723102"
# [1] "Iteration: 400 Time: 2024-07-30 07:00:58.302488 Scenario: 10"
# [1] "Iteration: 1000 Time: 2024-07-30 07:57:19.108441 Scenario: 10"
# [1] "Scenario: 11 begun at 2024-07-30 07:57:24.341301"
# [1] "Iteration: 500 Time: 2024-07-30 08:44:22.789593 Scenario: 11"
# [1] "Iteration: 600 Time: 2024-07-30 08:53:49.109465 Scenario: 11"
# [1] "Iteration: 1000 Time: 2024-07-30 09:31:39.830038 Scenario: 11"

# Issues with scenario 4 and 12
# Error in if (is.finite(x[xmax])) return(log1p(sum(exp(x[-xmax] - x[xmax]))) +  :
# argument is of length zero
# Final MCMC sample equivalent to less than 1000 independent draws.
# Please consider increasing the MCMC simulation size.
