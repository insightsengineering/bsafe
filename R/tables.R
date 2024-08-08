#' @title Display input data
#'
#' @description Display input data in the historical data tab
#'
#' @param data Uploaded Data
#' @param select_analysis Incidence proportion or Exposure-adjusted AE rate
#' @param saf_topic Selected safety topic to analyze/the adverse event of interest
#' @export
input_data_display <- function(data, select_analysis, saf_topic) {
  tab <- data
  rownames(tab) <- NULL
  if (select_analysis == "Incidence proportion") {
    colnames(tab) <- c(
      "STUDY ID", "Number of Patients in Arm",
      paste0("Number of Patients with ", saf_topic, " in Arm"), "Historical"
    )
  } else if (select_analysis == "Exposure-adjusted AE rate") {
    colnames(tab) <- c(
      "STUDY ID", "Number of Patients in Arm",
      paste0("Number of Patients with ", saf_topic, " in Arm"),
      "Total Exposure Time", "Historical"
    )
  }
  tab
}


#' @title MAP Prior table
#'
#' @description Function to display descriptive characteristics of the MAP prior,
#' mainly in the MAP Prior tab
#'
#' @param map_object MAP Prior object
#' @param select_analysis Incidence proportion or Exposure-adjusted AE rate
#' @param param_approx map prior; best fitting mixture model from parametric_approx()
#' @param ess_method ESS Method, currently only ELIR available
#' @param numerical TRUE or FALSE, for return values in Dataframe or to be displayed
#'
#' @export
model_summary_display <- function(map_object, select_analysis,
                                  param_approx, ess_method, numerical = FALSE) {
  if (is.null(map_object)) {
    return(NULL)
  } else if (select_analysis == "Incidence proportion") {
    sample_prop <- rstan::extract(map_object$fit)$theta_resp_pred
    stats_mat_prop <- data.frame(mcmc_desc(sample_prop),
      row.names = c("MAP Prior")
    )

    ESS <- RBesT::ess(param_approx, method = ess_method)

    if (numerical == TRUE) {
      stats_mat_prop$ESS <- ESS
      return(stats_mat_prop)
    } else {
      disp_mat_prop <- text_prop(stats_mat_prop = stats_mat_prop)
      disp_mat_prop$ESS <- round(ESS, 1)
      disp_mat_prop %>%
        dplyr::rename(Mean = mean, SD = sd, Median = median, "95% CrI" = cri)
    }
  } else if (select_analysis == "Exposure-adjusted AE rate") {
    sample_rate_log <- rstan::extract(map_object$fit)$theta_pred
    sample_rate_exp <- rstan::extract(map_object$fit)$theta_resp_pred

    stats_mat_rate <- data.frame(
      rbind(
        mcmc_desc(sample_rate_log),
        mcmc_desc(sample_rate_exp)
      ),
      row.names = c(
        "MAP Prior: log(hazard)",
        "MAP Prior: hazard"
      )
    )

    ESS <- RBesT::ess(param_approx, method = ess_method, sigma = 1)

    if (numerical == TRUE) {
      stats_mat_rate$ESS <- c(ESS, NA)
      return(stats_mat_rate)
    } else {
      disp_mat_rate <- round(stats_mat_rate, 4)
      disp_mat_rate <- cri_char(disp_mat_rate)

      disp_mat_rate$ESS <- c(round(ESS, 1), "Not applicable.")

      disp_mat_rate %>%
        dplyr::rename(Mean = mean, SD = sd, Median = median, "95% CrI" = cri)
    }
  }
}


#' @title Display Summary Stats
#'
#' @description Display summary stats of robust MAP prior and MAP prior in the Robust MAP tab
#'
#' @param map_object MAP Prior object
#' @param param_approx map prior; best fitting mixture model from parametric_approx()
#' @param ess_method Method to calculate ESS: Currently only 1 available ELIR
#' @param robust_map_object map prior; mixture distribution with a non-informative component from robust_map()
#' @param rob_ess_method Method to calculate ESS: Currently only 1 available ELIR
#' @param select_analysis Chosen Analysis proportion vs. rates
#' @param seed a seed to reproduce the results
#' @param numerical TRUE or FALSE, for return values in Dataframe or to be displayed
#'
#' @export
summary_stats_robust_map_prior_display <- function(
    map_object,
    select_analysis,
    param_approx,
    ess_method,
    robust_map_object,
    rob_ess_method,
    numerical = FALSE,
    seed) {
  if (select_analysis == "Incidence proportion") {
    # Summary statistics for MAP prior

    sample_prop <- rstan::extract(map_object$fit)$theta_resp_pred
    stats_mat_prop <- data.frame(
      rbind(
        mcmc_desc(sample_prop),
        rmix_desc(robust_map_object)
      ),
      row.names = c("MAP Prior", "Robustified MAP Prior")
    )

    ESS_MAP <- RBesT::ess(param_approx, method = ess_method)
    ESS_ROB <- RBesT::ess(robust_map_object, method = rob_ess_method)

    if (numerical == TRUE) {
      stats_mat_prop$ESS <- c(ESS_MAP, ESS_ROB)
      return(stats_mat_prop)
    } else {
      disp_mat_prop <- text_prop(stats_mat_prop = stats_mat_prop)
      disp_mat_prop$ESS <- round(c(ESS_MAP, ESS_ROB), 1)
      disp_mat_prop %>%
        dplyr::rename(Mean = mean, SD = sd, Median = median, "95% CrI" = cri)
    }
  } else if (select_analysis == "Exposure-adjusted AE rate") {
    # Rbest gives the model alreday out for exp(theta_pred), but we need
    # theta_pred for further analysis

    sample_rate_log <- rstan::extract(map_object$fit)$theta_pred
    sample_rate_exp <- rstan::extract(map_object$fit)$theta_resp_pred

    stats_mat_rate <- data.frame(
      rbind(
        mcmc_desc(sample_rate_log),
        rmix_desc(robust_map_object),
        mcmc_desc(sample_rate_exp),
        mcmc_desc(robust_map_object,
          trans = TRUE, seed = seed
        )
      ),
      row.names = c(
        "MAP Prior: log(hazard)",
        "robust MAP Prior: log(hazard)",
        "MAP Prior: hazard",
        "robust MAP Prior: hazard"
      )
    )

    ESS_MAP <- RBesT::ess(param_approx, method = ess_method, sigma = 1)
    ESS_ROB <- RBesT::ess(robust_map_object, method = ess_method, sigma = 1)

    if (numerical == TRUE) {
      stats_mat_rate$ESS <- c(ESS_MAP, ESS_ROB, NA, NA)
      return(stats_mat_rate)
    } else {
      disp_mat_rate <- round(stats_mat_rate, 4)
      disp_mat_rate <- cri_char(disp_mat_rate)

      disp_mat_rate$ESS <- c(
        round(c(ESS_MAP, ESS_ROB), 1), "Not applicable.", "Not applicable."
      )

      disp_mat_rate %>%
        dplyr::rename(Mean = mean, SD = sd, Median = median, "95% CrI" = cri)
    }
  }
}


#' @title Summary statistics for robust prior, likelihood, and posterior
#'
#' @description Summary statistics for robust prior, likelihood, and posterior in the New trial analysis tab
#'
#' @param select_analysis Incidence proportion or Exposure-adjusted AE rate
#' @param robust_map_object map prior; mixture distribution with a non-informative component from robust_map()
#' @param ess_method ESS Method, currently only ELIR available
#' @param current_trial_data information whether it is a historical (0) or current (1) trial
#' @param select_analysis Incidence proportion or Exposure-adjusted AE rate
#' @param post_dist posterior mixture distribution or MCMC samples from posterior_dist()
#' @param seed a seed
#' @param numerical TRUE or FALSE, for return values in Dataframe or to be displayed
#'
#' @export
summary_stat_all_display <- function(
    select_analysis,
    robust_map_object,
    ess_method,
    current_trial_data,
    post_dist,
    numerical = FALSE,
    seed) {
  if (is.na(seed)) {
    seed <- as.numeric(Sys.time())
  }

  # assign overall variables
  if (select_analysis == "Incidence proportion") {
    new_n <- current_trial_data[["new_v1"]]
    new_r <- current_trial_data[["new_v2"]]
  } else if (select_analysis == "Exposure-adjusted AE rate") {
    new_hzrate <- current_trial_data[["new_v1"]] / current_trial_data[["new_v2"]]
    new_unit_sd <- sqrt(1 / current_trial_data[["new_v1"]])
  }

  if (select_analysis == "Incidence proportion") {
    # RBesT cannot handle a ratio of 0 or 1, therefor:
    beta_param <- check_minmax(r = new_r, n = new_n)
    alpha <- beta_param["alpha"]
    beta <- beta_param["beta"]

    lik <- RBesT::mixbeta(lik = c(1, alpha, beta))

    stats_mat_prop <- data.frame(
      rbind(
        rmix_desc(robust_map_object),
        rmix_desc(lik),
        rmix_desc(post_dist)
      ),
      row.names = c("Robust MAP Prior", "Likelihood", "Posterior")
    )


    ESS_ROB <- RBesT::ess(robust_map_object)
    ESS_LIK <- NA # Likelihood does not have an effective sample size
    ESS_POST <- NA

    if (numerical == TRUE) {
      stats_mat_prop$ESS <- c(ESS_ROB, ESS_LIK, ESS_POST)
      return(stats_mat_prop)
    } else {
      disp_mat_prop <- text_prop(stats_mat_prop = stats_mat_prop)
      disp_mat_prop$ESS <- c(round(ESS_ROB, 1), "Not applicable.", "Not applicable.")
      disp_mat_prop %>%
        dplyr::rename(Mean = mean, SD = sd, Median = median, "95% CrI" = cri)
    }
  } else if (select_analysis == "Exposure-adjusted AE rate") {
    lik <- RBesT::mixnorm(lik = c(1, log(new_hzrate), new_unit_sd))

    post_dist <-
      posterior_dist(
        select_analysis = "Exposure-adjusted AE rate",
        robust_map_prior = robust_map_object,
        new_v1 = current_trial_data[["new_v1"]],
        new_v2 = current_trial_data[["new_v2"]],
        seed = seed
      )

    stats_mat_rate <- data.frame(
      rbind(
        rmix_desc(lik),
        rmix_desc(robust_map_object),
        rmix_desc(post_dist),
        mcmc_desc(robust_map_object,
          trans = TRUE, seed = seed
        ),
        mcmc_desc(post_dist,
          trans = TRUE, seed = seed
        )
      ),
      row.names = c(
        "Likelihood",
        "robust MAP Prior: log(hazard)", "Posterior: log(hazard)",
        "robust MAP Prior: hazard", "Posterior: hazard"
      )
    )

    ESS_ROB <- RBesT::ess(robust_map_object, sigma = 1)

    if (numerical == TRUE) {
      stats_mat_rate$ESS <- c(
        NA, ESS_ROB,
        NA, NA, NA
      )
      return(stats_mat_rate)
    } else {
      disp_mat_rate <- round(stats_mat_rate, 4)
      disp_mat_rate <- cri_char(disp_mat_rate)

      disp_mat_rate$ESS <- c(
        "Not applicable.", round((ESS_ROB), 1),
        "Not applicable.", "Not applicable.", "Not applicable."
      )

      disp_mat_rate %>%
        dplyr::rename(Mean = mean, SD = sd, Median = median, "95% CrI" = cri)
    }
  }
}

#' @title preset table for inference
#'
#' @description Table of preset statistical inference statement
#'
#' @param mix mixture distribution object
#' @param saf_topic Selected safety topic to analyze/the adverse event of interest
#' @param select_analysis Incidence proportion or Exposure-adjusted AE rate
#' @export
preset_stat_table <- function(
    mix,
    saf_topic,
    select_analysis) {
  certainty90 <- round(100 * RBesT::qmix(mix, 0.10, lower.tail = TRUE))
  certainty95 <- round(100 * RBesT::qmix(mix, 0.05, lower.tail = TRUE))
  certainty99 <- round(100 * RBesT::qmix(mix, 0.01, lower.tail = TRUE))

  postfix <- ""
  denominator <- 1

  if (select_analysis == "Incidence proportion") {
    postfix <- "%."
    denominator <- 1
    midfix <- " proportion of patients with "
  } else if (select_analysis == "Exposure-adjusted AE rate") {
    postfix <- "."
    denominator <- 100
    midfix <- " incidence rate of patients with "
  }

  inf_mat <- as.data.frame(
    matrix(
      c(
        paste0(
          "We are at least 90% certain that the", midfix, saf_topic,
          " is greater than ", certainty90 / denominator, postfix
        ),
        paste0(
          "We are at least 95% certain that the", midfix, saf_topic,
          " is greater than ", certainty95 / denominator, postfix
        ),
        paste0(
          "We are at least 99% certain that the", midfix, saf_topic,
          " is greater than ", certainty99 / denominator, postfix
        )
      ),
      nrow = 3, ncol = 1
    )
  )

  colnames(inf_mat) <- ""

  inf_mat
}


# TODO: Adding further helpers
# FIXME: Check correct display
