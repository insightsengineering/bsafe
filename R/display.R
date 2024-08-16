#' @title Parametric Approximation density display
#'
#' @description Display parametric approximation mixture density function in the MAP Prior Tab
#'
#' @param param_approx parametric approximation
#' @param select_analysis Incidence proportion or Exposure-adjusted AE rate
#'
#' @export
map_prior_function_display <- function(param_approx, select_analysis) {
  mixture_mat <- round(data.frame(param_approx[, seq_len(ncol(param_approx))]), 2)

  str_vec <- vector(length = ncol(mixture_mat))

  if (select_analysis == "Incidence proportion") {
    for (j in seq_len(ncol(mixture_mat))) {
      if (j == 1) {
        str_vec[j] <- paste0(
          "$$", mixture_mat[1, j], " \\cdot ", "Beta(", mixture_mat[2, j],
          ",", mixture_mat[3, j], ")", " + "
        )
      } else if (j != ncol(mixture_mat)) {
        str_vec[j] <- paste0(
          mixture_mat[1, j], " \\cdot ", "Beta(", mixture_mat[2, j],
          ",", mixture_mat[3, j], ")", " + "
        )
      } else {
        str_vec[ncol(mixture_mat)] <- paste0(
          mixture_mat[1, ncol(mixture_mat)], " \\cdot ", "Beta(",
          mixture_mat[2, ncol(mixture_mat)], ",",
          mixture_mat[3, ncol(mixture_mat)], ")", "$$"
        )
      }
    }
  } else if (select_analysis == "Exposure-adjusted AE rate") {
    for (j in seq_len(ncol(mixture_mat))) {
      if (j == 1) {
        str_vec[j] <- paste0(
          "$$", mixture_mat[1, j], " \\cdot ", "Normal(",
          mixture_mat[2, j], ",", mixture_mat[3, j], ")", " + "
        )
      } else if (j != ncol(mixture_mat)) {
        str_vec[j] <- paste0(
          mixture_mat[1, j], " \\cdot ", "Normal(", mixture_mat[2, j],
          ",", mixture_mat[3, j], ")", " + "
        )
      } else {
        str_vec[ncol(mixture_mat)] <- paste0(
          mixture_mat[1, ncol(mixture_mat)], " \\cdot ", "Normal(",
          mixture_mat[2, ncol(mixture_mat)], ",",
          mixture_mat[3, ncol(mixture_mat)], ")", "$$"
        )
      }
    }
  }
  shiny::withMathJax(
    shiny::h4(paste(str_vec, collapse = ""))
  )
}

#' @title Robust MAP density display
#'
#' @description Display robust MAP prior mixture density function in the robust MAP Prior Tab
#'
#' @param robust_map_object map prior; mixture distribution with a non-informative component from robust_map()
#' @param select_analysis Incidence proportion or Exposure-adjusted AE rate
#' @export
robust_map_prior_mix_dens_display <- function(
    robust_map_object,
    select_analysis) {
  rob_mixture_mat <- round(data.frame(robust_map_object[, seq_len(ncol(robust_map_object))]), 2)
  rob_str_vec <- vector(length = ncol(rob_mixture_mat))

  if (select_analysis == "Incidence proportion") {
    dist_label <- "Beta("
  } else {
    dist_label <- "Normal("
  }

  for (j in seq_len(ncol(rob_mixture_mat))) {
    rob_str_vec[j] <- paste0(
      rob_mixture_mat[1, j], " \\cdot ", dist_label, rob_mixture_mat[2, j],
      ",", rob_mixture_mat[3, j], ")", " + "
    )
  }

  rob_str_vec <- paste(rob_str_vec, collapse = "")
  rob_str_vec <- gsub(" \\+ $", "", rob_str_vec) # replace last + sign with empty str # nolint
  rob_str_vec <- paste0("$$", rob_str_vec, "$$")

  shiny::withMathJax(
    shiny::h4(rob_str_vec)
  )
}

#' @title Area under the curve
#'
#' @description Interpret area under the curve
#'
#' @param ae_prop Vector of quantiles
#' @param mix mixture distribution object
#' @param saf_topic Selected safety topic to analyze/the adverse event of interest
#' @param select_analysis Incidence proportion or Exposure-adjusted AE rate
#'
#' @export
area_under_the_curve <- function(ae_prop, mix, saf_topic, select_analysis) {
  # Interpret area under the curve

  if (select_analysis == "Incidence proportion") {
    certainty <- round(
      100 * (RBesT::pmix(mix, ae_prop[2], lower.tail = TRUE) - RBesT::pmix(mix, ae_prop[1], lower.tail = TRUE)),
      digits = 2
    ) # nolint

    # Bound certainty from 1% to 99%
    if (certainty >= 99.99) {
      certainty <- 99.99
    } else if (certainty <= 0.01) {
      certainty <- 0.01
    } else {
      certainty <- certainty
    }
    # Statistical inference based on certainty and quantile of incidence proportion

    if ((ae_prop[1] > 0) & (ae_prop[2] == 1)) {
      out <- paste0(
        "We are at least ", certainty,
        "% certain that the proportion of patients with ", saf_topic,
        " is greater than ", round(100 * ae_prop[1], 2), "%."
      )
    } else if (
      ((ae_prop[1] > 0) & (ae_prop[2] < 1)) | ((ae_prop[1] == 0) & (ae_prop[2] == 1))) { # nolint
      out <- paste0(
        "We are at least ", certainty,
        "% certain that the proportion of patients with ", saf_topic,
        " is between ", round(100 * ae_prop[1], 2), "% and ",
        round(100 * ae_prop[2], 2), "%."
      )
    } else if ((ae_prop[1] == 0) & (ae_prop[2] < 1)) {
      out <- paste0(
        "We are at least ", certainty,
        "% certain that the proportion of patients with ", saf_topic,
        " is less than ", round(100 * ae_prop[2], 2), "%."
      )
    } else if (ae_prop[1] == ae_prop[2]) {
      out <- paste0("Please ensure that the slider represents an interval.")
    }
    # hotfix
  } else if (select_analysis == "Exposure-adjusted AE rate") {
    ae_prop <- ae_prop * 100

    certainty <- round(
      100 * (RBesT::pmix(mix, ae_prop[2], lower.tail = TRUE) - RBesT::pmix(mix, ae_prop[1], lower.tail = TRUE)),
      digits = 2
    ) # nolint

    if (ae_prop[1] == ae_prop[2]) {
      out <- paste0("Please ensure that the slider represents an interval.")
    } else {
      out <- paste0(
        "We are at least ", certainty,
        "% certain that the exposure adjusted density rate of patients with ", saf_topic,
        " is between ", round(exp(ae_prop[1]), 2), " exp(", round(ae_prop[1], 2), ") and ",
        round(exp(ae_prop[2]), 2), " exp(", round(ae_prop[2], 2), ")."
      )
    }
  }

  return(out)
}
