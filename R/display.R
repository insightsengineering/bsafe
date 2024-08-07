#' @title Parametric Approximation density display
#'
#' @description Display parametric approximation mixture density function in the MAP Prior Tab
#'
#' @param map prior; best fitting mixture model from parametric_approx()
#' @param select_analysis Incidence proportion or Exposure-adjusted AE rate
#' @export
map_prior_function_display <- function(param_approx, select_analysis) {
  mixture_mat <- round(data.frame(param_approx[, 1:ncol(param_approx)]), 2)

  str_vec <- vector(length = ncol(mixture_mat))

  if (select_analysis == "Incidence proportion") {
    for (j in 1:ncol(mixture_mat)) {
      if (j == 1) {
        str_vec[j] <- paste0("$$", mixture_mat[1, j], " \\cdot ", "Beta(", mixture_mat[2, j], ",", mixture_mat[3, j], ")", " + ")
      } else if (j != ncol(mixture_mat)) {
        str_vec[j] <- paste0(mixture_mat[1, j], " \\cdot ", "Beta(", mixture_mat[2, j], ",", mixture_mat[3, j], ")", " + ")
      } else {
        str_vec[ncol(mixture_mat)] <- paste0(mixture_mat[1, ncol(mixture_mat)], " \\cdot ", "Beta(", mixture_mat[2, ncol(mixture_mat)], ",", mixture_mat[3, ncol(mixture_mat)], ")", "$$")
      }
    }
  } else if (select_analysis == "Exposure-adjusted AE rate") {
    for (j in 1:ncol(mixture_mat)) {
      if (j == 1) {
        str_vec[j] <- paste0("$$", mixture_mat[1, j], " \\cdot ", "Normal(", mixture_mat[2, j], ",", mixture_mat[3, j], ")", " + ")
      } else if (j != ncol(mixture_mat)) {
        str_vec[j] <- paste0(mixture_mat[1, j], " \\cdot ", "Normal(", mixture_mat[2, j], ",", mixture_mat[3, j], ")", " + ")
      } else {
        str_vec[ncol(mixture_mat)] <- paste0(mixture_mat[1, ncol(mixture_mat)], " \\cdot ", "Normal(", mixture_mat[2, ncol(mixture_mat)], ",", mixture_mat[3, ncol(mixture_mat)], ")", "$$")
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
robust_map_prior_mix_dens_display <- function(robust_map_object, select_analysis) {
  rob_mixture_mat <- round(data.frame(robust_map_object[, 1:ncol(robust_map_object)]), 2)
  rob_str_vec <- vector(length = ncol(rob_mixture_mat))

  if (select_analysis == "Incidence proportion") {
    dist_label <- "Beta("
  } else {
    dist_label <- "Normal("
  }

  for (j in 1:ncol(rob_mixture_mat)) {
    rob_str_vec[j] <- paste0(rob_mixture_mat[1, j], " \\cdot ", dist_label, rob_mixture_mat[2, j], ",", rob_mixture_mat[3, j], ")", " + ")
  }

  rob_str_vec <- paste(rob_str_vec, collapse = "")
  rob_str_vec <- gsub(" \\+ $", "", rob_str_vec) # replace last + sign with empty str
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
#' @export
area_under_the_curve <- function(ae_prop, mix, saf_topic) {
  # Interpret area under the curve
  certainty <- round(100 * (RBesT::pmix(mix, ae_prop[2], lower.tail = TRUE) - RBesT::pmix(mix, ae_prop[1], lower.tail = TRUE)))

  # Bound certainty from 1% to 99%
  if (certainty >= 99) {
    certainty <- 99
  } else if (certainty <= 1) {
    certainty <- 1
  } else {
    certainty <- certainty
  }
  # Statistical inference based on certainty and quantile of incidence proportion
  if ((ae_prop[1] > 0) & (ae_prop[2] == 1)) {
    paste0("We are at least ", certainty, "% certain that the proportion of patients with ", saf_topic, " is greater than ", round(100 * ae_prop[1]), "%.")
  } else if (((ae_prop[1] > 0) & (ae_prop[2] < 1)) | ((ae_prop[1] == 0) & (ae_prop[2] == 1))) {
    paste0("We are at least ", certainty, "% certain that the proportion of patients with ", saf_topic, " is between ", round(100 * ae_prop[1]), "% and ", round(100 * ae_prop[2]), "%.")
  } else if ((ae_prop[1] == 0) & (ae_prop[2] < 1)) {
    paste0("We are at least ", certainty, "% certain that the proportion of patients with ", saf_topic, " is less than ", round(100 * ae_prop[2]), "%.")
  } else if (ae_prop[1] == ae_prop[2]) {
    paste0("Please ensure that the slider represents an interval.")
  }
}
