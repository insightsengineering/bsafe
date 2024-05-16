#' @title Conjugate Posterior Analysis
#'
#' @description Calculates the posterior distribution for given prior robust_map_prior,
#' where the prior is a mixture of conjugate distributions. The posterior is then also
#' a mixture of conjugate distributions.
#' posterior MCMC samples from EXNEX analysis for exposure-adjusted analysis
#'
#' @param select_analysis Incidence proportion and Exposure-adjusted AE rate
#' @param robust_map_prior mixture distribution output from robust_map()
#' @param new_v1 first new variable #n for porportion and #events for rates
#' @param new_v2 second new variable #r for proprotion and observed exposure time for rates
#' @param input_data dataframe from data_table_prep() including the new trial
#' @param adj_tau numeric value of tau from tau_adjust()
#' @param p_exch probabilitiy of exchangeable compononent (use 1-robust_weight)
#' @param nex_mean Mean for Nex part (use robust_mean)
#' @param nex_sd  Standard Deviation for Nex part
#' @param seed   for reproduceability
#' @param testing for faster MCMC
#' @param explore for exploring a prior data conflict set to TRUE, to combine the data
#'
#' @return posterior mixture distribution for incidence proportion analysis and
#' @export
posterior_dist <- function(select_analysis,
                           robust_map_prior,
                           new_v1,
                           new_v2,
                           input_data,
                           explore = FALSE,
                           p_exch = NULL,
                           nex_mean = NULL,
                           nex_sd = NULL,
                           adj_tau,
                           seed = NULL,
                           testing = FALSE) {
  if (select_analysis == "Incidence proportion") {
    RBesT::postmix(robust_map_prior, n = new_v1, r = new_v2)
  } else if (select_analysis == "Exposure-adjusted AE rate") {
    RBesT::postmix(robust_map_prior, m = round(new_v1 / new_v2), se = 1 / new_v1)
  }
}
#' @title Robustify Mixture Priors
#'
#' @description Add a non-informative component to a mixture prior
#'
#' @param select_analysis Incidence proportion or Exposure-adjusted AE rate
#' @param param_approx the best fitting mixture model from parametric_approx()
#' @param robust_weight weight given to the non-informative component (0 < weight < 1)
#' @param robust_mean mean of the non-informative component
#' @param input_data read in data
#' @param adj_tau in between heterogeneity
#' @param seed seed for reprduceability
#' @param testing faster mcmc for testing
#'
#' @return new mixture distribution with an extra non-informative component
#' @export
robust_map <-
  function(select_analysis,
           param_approx,
           input_data,
           robust_weight,
           robust_mean,
           adj_tau,
           seed = NULL,
           testing = FALSE) {
    if (select_analysis == "Incidence proportion") {
      RBesT::robustify(param_approx, weight = robust_weight, mean = robust_mean)
    } else if (select_analysis == "Exposure-adjusted AE rate") {
      RBesT::robustify(param_approx, weight = robust_weight, mean = log(robust_mean), sigma = 1)
    }
  }


#' @title Automatic Fitting of Mixtures of Conjugate Distributions to a Sample
#'
#' @description Fitting a series of mixtures of conjugate distributions to a sample
#' using Expectation-Maximization. First, a Nc[1] component mixture is fitted, then a
#' Nc[2] component mixture, and so on. The mixture providing the best AIC value is selected.
#'
#' @param select_analysis Incidence proportion or Exposure-adjusted AE rate
#' @param map_prior an S3 object (list) of type gMAP for the incidence proportion analysis
#' or a vector of the posterior MCMC samples for the exposure-adjusted analysis
#'
#' @return the best fitting mixture model, i.e., the model with the lowest AIC
#' @export
parametric_approx <- function(select_analysis, map_prior) {
  if (select_analysis == "Incidence proportion") {
    RBesT::automixfit(map_prior, Nc = seq(3, 3))
  } else if (select_analysis == "Exposure-adjusted AE rate") {
    # Rbest gives the model alreday out for exp(theta_pred), but we need
    # theta_pred for further analysis
    theta_fit <- as.data.frame(map_prior$fit)
    theta_pred <- theta_fit$theta_pred
    RBesT::automixfit(theta_pred, Nc = seq(3, 3))
  }
}
#' @title Meta-Analytic-Predictive Analysis for Generalized Linear Models
#'
#' @description MAP prior calculation based on the selected analysis.
#' RBesT is used for the incidence proportion analysis. Alexander Stemke's code
#' is used for the exposure-adjusted analysis.
#' MAP prior should be computed using only historical data
#'
#' @param input_data dataframe from data_table_prep()
#' @param select_analysis Incidence proportion or Exposure-adjusted AE rate
#' @param tau_dist assumed distribution of tau
#' @param adj_tau numeric value from tau_adjust
#' @param seed a seed as input for reproducability
#' @param testing for testing purposes faster MCMC
#'
#' @return an S3 object (list) of type gMAP for the incidence proportion analysis
#' or a vector of the posterior MCMC samples for the exposure-adjusted analysis
#' @export
map_prior_func <-
  function(input_data,
           select_analysis,
           tau_dist,
           adj_tau,
           seed = NULL,
           testing = FALSE) {
    input_data <- input_data %>% dplyr::filter(HIST == 1)

    if (testing == FALSE) {
      # for exact method look stan/jags up, but if no rounding etc n_mcmc = (n_iter-n_warm)*n_chain/n_thin
      n_iter <- 6000
      n_burn <- 2000
      n_chains <- 4
      n_thin <- 4
    } else {
      n_iter <- 400
      n_burn <- 100
      n_chains <- 2
      n_thin <- 2
    }

    if (is.null(seed)) {
      seed <- as.numeric(Sys.time())
    }

    if (select_analysis == "Incidence proportion") {
      set.seed(seed)
      RBesT::gMAP(
        cbind(N_WITH_AE, N - N_WITH_AE) ~ 1 | STUDYID,
        data = input_data,
        tau.dist = tau_dist,
        tau.prior = adj_tau,
        beta.prior = 2,
        iter = getOption("RBesT.MC.iter", n_iter),
        warmup = getOption("RBesT.MC.warmup", n_burn),
        thin = getOption("RBesT.MC.thin", n_thin),
        chains = getOption("RBesT.MC.chains", n_chains),
        family = "binomial"
      )
    } else if (select_analysis == "Exposure-adjusted AE rate") {
      set.seed(seed)
      RBesT::gMAP(N_WITH_AE ~ 1 + offset(log(TOT_EXP)) | STUDYID,
        data = input_data,
        tau.dist = tau_dist,
        tau.prior = adj_tau, ## assuming moderate heterogeniety
        beta.prior = 1,
        iter = getOption("RBesT.MC.iter", n_iter),
        warmup = getOption("RBesT.MC.warmup", n_burn),
        thin = getOption("RBesT.MC.thin", n_thin),
        chains = getOption("RBesT.MC.chains", n_chains),
        family = "poisson"
      )
    }
  }
#' @title Tau Adjustment for Amount of Historical Borrowing
#'
#' @description Change the value of tau to control the amount of historical borrowing
#'
#' @param select_analysis Incidence proportion or Exposure-adjusted AE rate
#' @param hist_borrow the amount of historical borrowing (ranges from small to very large)
#'
#' @return the numeric value of tau based on the selected amount of historical borrowing
#' @export
tau_adjust <- function(select_analysis, hist_borrow) {
  if (select_analysis == "Incidence proportion") {
    # beta is chosen to be 2 therefor tau/beta should be as recommended by Neuenschwander & Co.
    # see RBesT documentation

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
  } else if (select_analysis == "Exposure-adjusted AE rate") {
    # sigma is chosen to be 1 therefor tau/sigma should be ...

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
  }
}
