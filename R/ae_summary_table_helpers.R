#' @title Data availability
#'
#' @description Gives as output 4 T/F values whether input data is available, ix stands for index
#'
#' @param hist_trt Historical treatment data
#' @param hist_ctr Historical control data
#' @param trt_current current treatment data
#' @param ctr_current current control data
#'
#' @return A vector to inform where information is available
#' @export
data_available <- function(hist_trt, hist_ctr, trt_current, ctr_current) {
  data_check <- c(
    hist_trt = nrow(hist_trt) > 0,
    hist_ctr = nrow(hist_ctr) > 0,
    trt_current = nrow(trt_current) > 0,
    ctr_current = nrow(ctr_current) > 0
  )

  return(data_check)
}


# Function to add a new row
#' @title Update warning text data frame
#'
#' @description Adds another row, if data is missing in a comparison
#'
#' @param df Data frame with existing warning texts
#' @param data_check Output from data_available()
#' @param topic Safety topics
#' @param group Comparison number
#' @param analysis Type of analysis
#'
#' @return Another row for a dataframe returning warning texts about missing information and its consequence.
#' @export
add_row <- function(df, data_check, topic, group, analysis) {
  dc <- sum(data_check)

  if (dc == 0) {
    issue <- "Warning"
    message <- paste0(c(" No data available for ", topic, " in any Arm ", "."))
  } else if (dc < 4) {
    misidx <- which(data_check == FALSE)
    issue <- "Warning"
    for (idx in misidx) {
      message <- "! "
      if (idx == 1) {
        message <- paste0(c(
          message, " No historic data available for ",
          topic, " in ", "Arm A", ".", "An uniformative prior was used."
        ))
      }

      if (idx == 2) {
        message <- paste0(c(
          message, " No historic data available for ",
          topic, " in ", "Arm B", ".", "An uniformative prior was used."
        ))
      }

      if (idx == 3) {
        message <- paste0(c(
          message, " No current data available for ",
          topic, " in ", "Arm A", ".", "The prior information was used as posterior."
        ))
      }

      if (idx == 3) {
        message <- paste0(c(
          message, " No current data available for ",
          topic, " in ", "Arm B", ".", "The prior information was used as posterior."
        ))
      }
    }
  }

  if (sum(data_check) < 4) {
    new_row <- data.frame(
      Issue = issue,
      Analysis = analysis,
      Comparison = group,
      Topic = topic,
      Message = message,
      stringsAsFactors = FALSE
    )
    df <- rbind(df, new_row)
  }

  return(df)
}



#' @title Naive incidence
#'
#' @description Gives the naive estimation outputs for the whole data set
#'
#' @param data Historic and Current data
#' @param array_inci The array setup for incidence proportions
#' @param arm "Arm A" (treatment) or "Arm B" (control)
#' @param topic Safety Topic
#' @param group Comparison number
#'
#' @return Returns an array with naive incidence rates
#' @export
inci_naiv <- function(data, array_inci, arm, topic, group) {
  array_inci[topic, "r", arm, group] <- sum(data$N_WITH_AE)
  array_inci[topic, "n", arm, group] <- sum(data$N)
  array_inci[topic, "%", arm, group] <- array_inci[topic, "r", arm, group] / array_inci[topic, "n", arm, group]

  return(array_inci)
}


#' @title  Array rmix
#'
#' @description Assigning values from a rmix object from RBesT to the array_inci
#'
#' @param rmix_obj A S3 class matrix from R RBesT with a mixture distribution
#' @param array The array setup
#' @param arm "Arm A" (treatment) or "Arm B" (control)
#' @param topic Safety Topic
#' @param group Comparison number
#' @param lb Lower bound for credible interval
#' @param ub Upper bound for credible interval
#'
#' @export
array_rmix <- function(rmix_obj, array, arm, topic, group, lb = lb, ub = ub) {
  ana <- rmix_desc(rmix_obj = rmix_obj, crilb = lb, criub = ub)

  array[topic, "cri_l", arm, group] <- ana$crilb
  array[topic, "cri_u", arm, group] <- ana$criub
  array[topic, "post_median", arm, group] <- ana$median
  array[topic, "post_mean", arm, group] <- ana$mean

  return(array)
}


#' @title  Array Comparison
#'
#' @description Assigning values from a sample to array_inci_comp
#'
#' @param array_comp The array setup for comparisons
#' @param comp "Risk Diff" or "Risk Ratio"
#' @param topic Variable to be analyzed
#' @param group The current comparison group
#' @param pr_sample posterior sample for difference or ratio
#' @param crilb Cri lower bound
#' @param criub Cri upper bound
#' @param array_ana Setup for analyzed array
#' @param ana Type of analysis incidence proportions or rates
#'
#' @export
array_comp <- function(ana, array_comp = NA, comp, pr_sample,
                       topic = topic, group, crilb = lb, criub = ub, array_ana = NA) {
  if (ana == "inci") {
    if (comp == "Risk Diff") {
      array_comp[topic, "naive", "Risk Diff", group] <- as.numeric(
        array_ana[topic, "%", "Arm A", group]
      ) - as.numeric(array_ana[topic, "%", "Arm B", group])
    } else if (comp == "Risk Ratio") {
      array_comp[topic, "naive", "Risk Ratio", group] <- as.numeric(
        array_ana[topic, "%", "Arm A", group]
      ) / as.numeric(array_ana[topic, "%", "Arm B", group])
    }
  } else if (ana == "rate") {
    if (comp == "Risk Diff") {
      array_comp[topic, "naive", "Risk Diff", group] <- as.numeric(
        array_ana[topic, "naive", "Arm A", group]
      ) - as.numeric(array_ana[topic, "naive", "Arm B", group])
    } else if (comp == "Risk Ratio") {
      array_comp[topic, "naive", "Risk Ratio", group] <- as.numeric(
        array_ana[topic, "naive", "Arm A", group]
      ) / as.numeric(array_ana[topic, "naive", "Arm B", group])
    }
  }

  array_comp[topic, "post_mean", comp, group] <- mean(pr_sample)
  array_comp[topic, "post_median", comp, group] <- median(pr_sample)
  array_comp[topic, "cri_l", comp, group] <- quantile(pr_sample, crilb)
  array_comp[topic, "cri_u", comp, group] <- quantile(pr_sample, criub)
  array_comp[topic, "A>B%", comp, group] <- sum(pr_sample > 0) / length(pr_sample)

  return(array_comp)
}
