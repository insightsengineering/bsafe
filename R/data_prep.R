#' Title
#'
#' @param input_data
#' @export
trt_data_wrangler <- function(input_data) {
  choices_trt <- input_data %>%
    dplyr::distinct(input_data[, "ARM"])
  colnames(choices_trt) <- "Arms"
  return(choices_trt)
}

# ae events input filter
#' Title
#'
#' @param input_data
#' @param selected_trt
#' @export
ae_events_wrangler <- function(input_data, selected_trt) {
  safety_topics <- as.character(unlist(input_data[, "SAF_TOPIC"]))
  choices_ae <- safety_topics[as.character(unlist(input_data[, "ARM"])) == selected_trt]
  return(choices_ae)
}
# MAP Prior, Robust MAP Prior, Likelihood, or Posterior Distribution Samples in the Decision Making Tab
#' Title
#'
#' @param select_dist
#' @param param_approx
#' @param select_analysis
#' @param new_trial_analysis
#'
#' @export
sampling_all_plot <- function(select_analysis, select_dist, param_approx, new_trial_analysis) {
  if (select_analysis == "Incidence proportion") {
    x <- seq(0.0001, 1, length = 500)
  } else if (select_analysis == "Exposure-adjusted AE rate") {
    a <- RBesT::qmix(param_approx, 0.02)
    b <- RBesT::qmix(param_approx, 0.98)
    x <- seq(a, b, length = 500)
    rm(a, b)
  }

  if (select_dist == "MAP Prior") {
    # LS
    # MAP prior
    mixture_mat <- data.frame(param_approx[, 1:ncol(param_approx)])
    str_vec <- vector(length = ncol(mixture_mat))

    if (select_analysis == "Incidence proportion") {
      for (j in 1:ncol(mixture_mat)) {
        if (j != ncol(mixture_mat)) {
          str_vec[j] <- paste0(mixture_mat[1, j], "*", "dbeta(x, shape1 = ", mixture_mat[2, j], ", shape2 = ", mixture_mat[3, j], ")", " + ")
        } else {
          str_vec[ncol(mixture_mat)] <- paste0(mixture_mat[1, ncol(mixture_mat)], "*", "dbeta(x, shape1 = ", mixture_mat[2, ncol(mixture_mat)], ", shape2 = ", mixture_mat[3, ncol(mixture_mat)], ")")
        }
      }

      # Prior density function

      Prior <- eval(parse(text = paste(str_vec, collapse = "")))

      # Create dataframe for ggplot
      df <- data.frame(
        Probability = x,
        Density = rep("MAP Prior", each = 500),
        Value = Prior
      )
    } else if (select_analysis == "Exposure-adjusted AE rate") {
      for (j in 1:ncol(mixture_mat)) {
        if (j != ncol(mixture_mat)) {
          str_vec[j] <-
            paste0(
              mixture_mat[1, j],
              "*",
              "dnorm(x, mean = ",
              mixture_mat[2, j],
              ", sd = ",
              mixture_mat[3, j],
              ")",
              " + "
            )
        } else {
          str_vec[ncol(mixture_mat)] <-
            paste0(
              mixture_mat[1, ncol(mixture_mat)],
              "*",
              "dnorm(x, mean = ",
              mixture_mat[2, ncol(mixture_mat)],
              ", sd = ",
              mixture_mat[3, ncol(mixture_mat)],
              ")"
            )
        }
      }

      # Prior density function

      Prior <- eval(parse(text = paste(str_vec, collapse = "")))

      # Create dataframe for ggplot
      df <- data.frame(
        Probability = x,
        Density = rep("MAP Prior", each = 500),
        Value = Prior
      )
    }
  } else if (select_dist == "Robust MAP Prior") {
    df <- new_trial_analysis %>% dplyr::filter(Density == "Robust MAP Prior")
  } else if (select_dist == "Posterior") {
    df <- new_trial_analysis %>% dplyr::filter(Density == "Posterior")
  } else if (select_dist == "Likelihood") {
    df <- new_trial_analysis %>% dplyr::filter(Density == "Likelihood")
  }
}


# Mixture distributions for MAP Prior, Robust MAP Prior, Likelihood, or Posterior Distribution
#' Title
#'
#' @param current_trial_data
#' @param select_dist
#' @param param_approx
#' @param robust_map_object
#' @param select_analysis
#' @param post_dist
#'
#' @export
mix_distribution_all <- function(select_analysis, current_trial_data, select_dist, param_approx, robust_map_object, post_dist) {
  # assign overall variables
  if (select_analysis == "Incidence proportion") {
    new_n <- current_trial_data[["new_v1"]]
    new_r <- current_trial_data[["new_v2"]]
  } else if (select_analysis == "Exposure-adjusted AE rate") {
    new_n_with_ae <- current_trial_data[["new_v1"]]
    new_tot_exp <- current_trial_data[["new_v2"]]
  }

  if (select_dist == "MAP Prior") {
    param_approx
  } else if (select_dist == "Robust MAP Prior") {
    if (select_analysis == "Incidence proportion") {
      robust_map_object
    } else if (select_analysis == "Exposure-adjusted AE rate") {
      robust_map_object
    }
  } else if (select_dist == "Posterior") {
    if (select_analysis == "Incidence proportion") {
      post_dist
    } else if (select_analysis == "Exposure-adjusted AE rate") {
      post_dist
    }
  } else if (select_dist == "Likelihood") {
    if (select_analysis == "Incidence proportion") {
      RBesT::mixbeta(inf = c(1, new_r + 1, new_n - new_r + 1))
    } else if (select_analysis == "Exposure-adjusted AE rate") {
      theta <- new_n_with_ae / new_tot_exp
      like_sample <- rnorm(1000, mean = log(theta), sd = 1 / new_n_with_ae)
      RBesT::automixfit(like_sample, Nc = seq(3, 3))

      # RBesT::mixnorm(inf = c(1, log(new_n_with_ae/new_tot_exp), sqrt(1/new_n_with_ae)))
    }
  }
}
#' @title Dataframe preparation for plotting robust map prior, likelihood, and posterior distributions
#'
#' @description Creates a dataframe for plotting robust map prior, likelihood, and posterior distributions
#'
#' @param select_analysis Incidence proportion or Exposure-adjusted AE rate
#' @param robust_map_prior robust map prior; mixture distribution with a non-informative component from robust_map()
#' @param new_v1 #new_n sample size in the new trial for the selected safety topic and background treatment
#' @param new_v2 #new_r number of patients with at least one occurrence of the AE of interest
#' @param post_dist posterior mixture distribution or MCMC samples from posterior_dist()
#'
#' @return dataframe for plotting robust map prior, likelihood, and posterior distributions
#' @export
new_trial_compare <-
  function(select_analysis,
           robust_map_prior,
           new_v1,
           new_v2,
           post_dist) {
    # assign overall variables
    if (select_analysis == "Incidence proportion") {
      x <- seq(0.0001, 1, length = 500)
      new_n <- new_v1
      new_r <- new_v2
    } else if (select_analysis == "Exposure-adjusted AE rate") {
      a <- RBesT::qmix(robust_map_prior, 0.02)
      b <- RBesT::qmix(robust_map_prior, 0.98)
      x <- seq(a, b, length = 500)
      rm(a, b)
      new_n_with_ae <- new_v1
      new_tot_exp <- new_v2
    }

    # Robust MAP prior
    rob_mixture_mat <-
      data.frame(robust_map_prior[, 1:ncol(robust_map_prior)])
    rob_str_vec <- vector(length = ncol(rob_mixture_mat))


    # Robust MAP Prior density function
    if (select_analysis == "Incidence proportion") {
      for (j in 1:ncol(rob_mixture_mat)) {
        if (j != ncol(rob_mixture_mat)) {
          rob_str_vec[j] <-
            paste0(
              rob_mixture_mat[1, j],
              "*",
              "dbeta(x, shape1 = ",
              rob_mixture_mat[2, j],
              ", shape2 = ",
              rob_mixture_mat[3, j],
              ")",
              " + "
            )
        } else {
          rob_str_vec[ncol(rob_mixture_mat)] <-
            paste0(
              rob_mixture_mat[1, ncol(rob_mixture_mat)],
              "*",
              "dbeta(x, shape1 = ",
              rob_mixture_mat[2, ncol(rob_mixture_mat)],
              ", shape2 = ",
              rob_mixture_mat[3, ncol(rob_mixture_mat)],
              ")"
            )
        }
      }
    } else if (select_analysis == "Exposure-adjusted AE rate") {
      for (j in 1:ncol(rob_mixture_mat)) {
        if (j != ncol(rob_mixture_mat)) {
          rob_str_vec[j] <-
            paste0(
              rob_mixture_mat[1, j],
              "*",
              "dnorm(x, mean = ",
              rob_mixture_mat[2, j],
              ", sd = ",
              rob_mixture_mat[3, j],
              ")",
              " + "
            )
        } else {
          rob_str_vec[ncol(rob_mixture_mat)] <-
            paste0(
              rob_mixture_mat[1, ncol(rob_mixture_mat)],
              "*",
              "dnorm(x, mean = ",
              rob_mixture_mat[2, ncol(rob_mixture_mat)],
              ", sd = ",
              rob_mixture_mat[3, ncol(rob_mixture_mat)],
              ")"
            )
        }
      }
    }


    rob_Prior <-
      eval(parse(text = paste(rob_str_vec, collapse = "")))

    # Likelihood
    if (select_analysis == "Incidence proportion") {
      Likelihood <-
        dbeta(
          x = x,
          shape1 = new_r + 1,
          shape2 = new_n - new_r + 1
        )
    } else if (select_analysis == "Exposure-adjusted AE rate") {
      # Distribution of log lambda
      # http://www.math.wm.edu/~leemis/chart/UDR/PDFs/GammaLoggamma.pdf

      theta <- new_n_with_ae / new_tot_exp
      Likelihood <- dnorm(x, mean = log(theta), sd = 1 / new_n_with_ae)
    }
    # Update prior distribution - compute the conditional distribution given the data and the prior
    # Compute the posterior distribution

    post_mat <- data.frame(post_dist[, 1:ncol(post_dist)])
    post_str_vec <- vector(length = ncol(post_mat))


    # Posterior density function
    if (select_analysis == "Incidence proportion") {
      for (j in 1:ncol(post_mat)) {
        if (j != ncol(post_mat)) {
          post_str_vec[j] <-
            paste0(
              post_mat[1, j],
              "*",
              "dbeta(x, shape1 = ",
              post_mat[2, j],
              ", shape2 = ",
              post_mat[3, j],
              ")",
              " + "
            )
        } else {
          post_str_vec[ncol(post_mat)] <-
            paste0(
              post_mat[1, ncol(post_mat)],
              "*",
              "dbeta(x, shape1 = ",
              post_mat[2, ncol(post_mat)],
              ", shape2 = ",
              post_mat[3, ncol(post_mat)],
              ")"
            )
        }
      }
    } else if (select_analysis == "Exposure-adjusted AE rate") {
      for (j in 1:ncol(post_mat)) {
        if (j != ncol(post_mat)) {
          post_str_vec[j] <-
            paste0(
              post_mat[1, j],
              "*",
              "dnorm(x, mean = ",
              post_mat[2, j],
              ", sd = ",
              post_mat[3, j],
              ")",
              " + "
            )
        } else {
          post_str_vec[ncol(post_mat)] <-
            paste0(
              post_mat[1, ncol(post_mat)],
              "*",
              "dnorm(x, mean = ",
              post_mat[2, ncol(post_mat)],
              ", sd = ",
              post_mat[3, ncol(post_mat)],
              ")"
            )
        }
      }
    }

    Posterior <-
      eval(parse(text = paste(post_str_vec, collapse = "")))

    # Create dataframe for ggplot
    df <- data.frame(
      Probability = rep(x, 3),
      Density = rep(c(
        "Robust MAP Prior", "Likelihood", "Posterior"
      ), each = 500),
      Value = c(rob_Prior, Likelihood, Posterior)
    )
    df$Density <-
      factor(df$Density,
        levels = c("Robust MAP Prior", "Likelihood", "Posterior")
      )
    return(df)
  }
#' @title robust_compare
#'
#' @description Creates the dataframe for plotting the comparison between the map prior and robust map prior
#'
#' @param select_analysis Incidence proportion or Exposure-adjusted AE rate
#' @param robust_map_prior robust map prior; mixture distribution with a non-informative component from robust_map()
#' @param param_approx map prior; best fitting mixture model from parametric_approx()
#'
#' @return dataframe of comparison of map prior and robust map prior for ggplot
#' @export
robust_compare <- function(select_analysis,
                           robust_map_prior,
                           param_approx) {
  rob_mixture_mat <-
    data.frame(robust_map_prior[, 1:ncol(robust_map_prior)])
  rob_str_vec <- vector(length = ncol(rob_mixture_mat))

  # Robust MAP Prior density function
  if (select_analysis == "Incidence proportion") {
    x <- seq(0.0001, 1, length = 500)
    for (j in 1:ncol(rob_mixture_mat)) {
      if (j != ncol(rob_mixture_mat)) {
        rob_str_vec[j] <-
          paste0(
            rob_mixture_mat[1, j],
            "*",
            "dbeta(x, shape1 = ",
            rob_mixture_mat[2, j],
            ", shape2 = ",
            rob_mixture_mat[3, j],
            ")",
            " + "
          )
      } else {
        rob_str_vec[ncol(rob_mixture_mat)] <-
          paste0(
            rob_mixture_mat[1, ncol(rob_mixture_mat)],
            "*",
            "dbeta(x, shape1 = ",
            rob_mixture_mat[2, ncol(rob_mixture_mat)],
            ", shape2 = ",
            rob_mixture_mat[3, ncol(rob_mixture_mat)],
            ")"
          )
      }
    }
  } else if (select_analysis == "Exposure-adjusted AE rate") {
    a <- RBesT::qmix(param_approx, 0.02)
    b <- RBesT::qmix(param_approx, 0.98)
    x <- seq(a, b, length = 500)
    for (j in 1:ncol(rob_mixture_mat)) {
      if (j != ncol(rob_mixture_mat)) {
        rob_str_vec[j] <-
          paste0(
            rob_mixture_mat[1, j],
            "*",
            "dnorm(x, mean = ",
            rob_mixture_mat[2, j],
            ", sd = ",
            rob_mixture_mat[3, j],
            ")",
            " + "
          )
      } else {
        rob_str_vec[ncol(rob_mixture_mat)] <-
          paste0(
            rob_mixture_mat[1, ncol(rob_mixture_mat)],
            "*",
            "dnorm(x, mean = ",
            rob_mixture_mat[2, ncol(rob_mixture_mat)],
            ", sd = ",
            rob_mixture_mat[3, ncol(rob_mixture_mat)],
            ")"
          )
      }
    }
  }

  rob_Prior <-
    eval(parse(text = paste(rob_str_vec, collapse = "")))

  # MAP prior
  mixture_mat <- data.frame(param_approx[, 1:ncol(param_approx)])
  str_vec <- vector(length = ncol(mixture_mat))

  # MAP Prior density function
  if (select_analysis == "Incidence proportion") {
    for (j in 1:ncol(mixture_mat)) {
      if (j != ncol(mixture_mat)) {
        str_vec[j] <-
          paste0(
            mixture_mat[1, j],
            "*",
            "dbeta(x, shape1 = ",
            mixture_mat[2, j],
            ", shape2 = ",
            mixture_mat[3, j],
            ")",
            " + "
          )
      } else {
        str_vec[ncol(mixture_mat)] <-
          paste0(
            mixture_mat[1, ncol(mixture_mat)],
            "*",
            "dbeta(x, shape1 = ",
            mixture_mat[2, ncol(mixture_mat)],
            ", shape2 = ",
            mixture_mat[3, ncol(mixture_mat)],
            ")"
          )
      }
    }
  } else if (select_analysis == "Exposure-adjusted AE rate") {
    for (j in 1:ncol(mixture_mat)) {
      if (j != ncol(mixture_mat)) {
        str_vec[j] <-
          paste0(
            mixture_mat[1, j],
            "*",
            "dnorm(x, mean = ",
            mixture_mat[2, j],
            ", sd = ",
            mixture_mat[3, j],
            ")",
            " + "
          )
      } else {
        str_vec[ncol(mixture_mat)] <-
          paste0(
            mixture_mat[1, ncol(mixture_mat)],
            "*",
            "dnorm(x, mean = ",
            mixture_mat[2, ncol(mixture_mat)],
            ", sd = ",
            mixture_mat[3, ncol(mixture_mat)],
            ")"
          )
      }
    }
  }

  Prior <- eval(parse(text = paste(str_vec, collapse = "")))

  # Create dataframe for ggplot
  df <- data.frame(
    Probability = rep(x, 2),
    Density = rep(c("MAP Prior", "Robust MAP Prior"), each = 500),
    Value = c(Prior, rob_Prior)
  )
  df$Density <-
    factor(df$Density, levels = c("MAP Prior", "Robust MAP Prior"))

  return(df)
}
#' @title Summary Table Preparation
#'
#' @description Returns a dataframe for the selected safety topic with the study ID,
#' the number of patients in the selected arm, the number of patients in the arm with at
#' least one occurrence of the adverse event, and the total exposure in patient years
#' if the exposure-adjusted analysis is chosen.
#'
#' @param input_data the raw summary level data
#' @param select_analysis Incidence proportion or Exposure-adjusted AE rate
#' @param saf_topic the adverse event of interest
#' @param current_trial information whether it is a historical (0) or current (1) trial
#' @param select_btrt selected treatment arms
#'
#' @return a dataframe for the selected safety topic with the study ID,
#' the number of patients in the selected arm, the number of patients in the arm with at
#' least one occurrence of the adverse event, and the total exposure in patient years
#' if the exposure-adjusted analysis is chosen
#' @export
data_table_prep <-
  function(input_data,
           select_analysis,
           saf_topic,
           select_btrt,
           bool_pooled = FALSE,
           current_trial = FALSE) {
    # Rename columns
    dat <- input_data
    dat <- dat[order(dat$STUDYID), ]

    # Filter for the selected background treatment and safety topic
    # Return total exposure as well for exposure-adjusted analysis
    if (select_analysis == "Incidence proportion") {
      dat <- dat %>%
        dplyr::filter(ARM %in% select_btrt &
          SAF_TOPIC == saf_topic) %>%
        dplyr::select(STUDYID, N, N_WITH_AE, HIST) %>%
        na.omit()
    } else if (select_analysis == "Exposure-adjusted AE rate") {
      dat <- dat %>%
        dplyr::filter(ARM %in% select_btrt &
          SAF_TOPIC == saf_topic) %>%
        dplyr::select(STUDYID, N, N_WITH_AE, TOT_EXP, HIST) %>%
        na.omit()
    }

    # test this

    if (bool_pooled == TRUE) {
      if (select_analysis == "Incidence proportion") {
        dat <- dat %>%
          dplyr::group_by(STUDYID, HIST) %>%
          dplyr::summarise(N = sum(N), N_WITH_AE = sum(N_WITH_AE), .groups = "drop")
        dat <- dat %>% dplyr::select(STUDYID, N, N_WITH_AE, HIST)
      }

      if (select_analysis == "Exposure-adjusted AE rate") {
        dat <- dat %>%
          dplyr::group_by(STUDYID, HIST) %>%
          dplyr::summarise(N = sum(N), N_WITH_AE = sum(N_WITH_AE), TOT_EXP = sum(TOT_EXP), .groups = "drop")
        dat <- dat %>% dplyr::select(STUDYID, N, N_WITH_AE, TOT_EXP, HIST)
      }
    }

    return(dat)
  }
