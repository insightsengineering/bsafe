#' @title Describe MCMC sample
#'
#' @description  Helper for descriptive statistics of sample
#'
#' @param mcmc_obj A mcmc object, or in general a sample vector
#' @param crilb Lower Bound for the credible interval
#' @param criub Upper Bound for the credible interval
#' @param seed  Seed for reproducibility
#' @param trans If a Matrix with a mixture is given, it will be sampled from it.
#' Needed for interpretation of log(hazard)
#'
#' @return A data frame including descriptive statistics
#' @export
mcmc_desc <- function(mcmc_obj = NA, crilb = 0.025, criub = 0.975,
                      seed = NA, trans = FALSE) {
  if (is.na(seed)) {
    seed <- as.numeric(Sys.time())
  }



  if (trans == TRUE) {
    set.seed(seed)
    mcmc_obj <- exp(RBesT::rmix(mcmc_obj, 100000))
  }


  res <- data.frame(
    mean   = mean(mcmc_obj),
    sd     = sd(mcmc_obj),
    median = median(mcmc_obj),
    crilb  = quantile(mcmc_obj, probs = crilb),
    criub  = quantile(mcmc_obj, probs = criub)
  )

  return(res)
}


#' @title Describe Mixture Distributions
#'
#' @description Helper for descriptive statistics of a matrix with mixture
#' distribution, as object from RBesT
#'
#' @param rmix_obj A dataframe containing a mixture distribution
#' @param crilb Lower Bound for the credible interval
#' @param criub Upper Bound for the credible interval
#' @param deci Number of decimals for rounding
#'
#' @return A data frame including descriptive statistics
#' @export
#'
rmix_desc <- function(rmix_obj = NA, crilb = 0.025, criub = 0.975, deci = 4) {
  ana_obj <- summary(rmix_obj,
    quantile.type = 7, digits = deci, na.rm = TRUE,
    names = TRUE, type = c("quartiles"),
    probs = c(crilb, 0.5, criub)
  )

  res <- data.frame(
    mean   = ana_obj["mean"],
    sd     = ana_obj["sd"],
    median = ana_obj["50.0%"],
    crilb  = ana_obj[3],
    criub  = ana_obj[5]
  )

  return(res)
}



#' @title Helper for incidence proportion display
#'
#' @description  Data Frame to display proportions in %
#'
#' @param stats_mat_prop Matrix with statistical information for the proportional case
#'
#' @return Text for display
#' @export
text_prop <- function(stats_mat_prop = NA) {
  text_mat_prop <- round(stats_mat_prop * 100, digits = 4)
  text_mat_prop[] <- lapply(
    text_mat_prop,
    function(x) paste0(as.character(x), "%")
  )
  text_mat_prop <- cri_char(text_mat_prop)

  return(text_mat_prop)
}


#' @title Helper for credible interval display
#'
#' @description  Put Credible Interval in one Character
#'
#' @param df A dataframe that is expected to include the columns crilb and criub
#'
#' @return a data frame including the credible intervals
#' @export
cri_char <- function(df = NA) {
  df$cri <- mapply(function(x, y) paste0("(", x, ", ", y, ")"), df$crilb, df$criub)

  df$crilb <- NULL
  df$criub <- NULL

  return(df)
}

#' @title Beta Dist Update check
#'
#' @description Many unnecessary error occur with RBesT when in a Beta Dist
#' the parameters are <1
#'
#' @param r Response = alpha
#' @param n Total = alpha + beta
#'
#' @return A vector containing alpha and beta for the beta distribution
#' @export
#'
check_minmax <- function(r, n) {
  # RBesT cannot handle a ratio of 0 or 1, therefor:
  if (r == 0) {
    alpha <- r + 1
    beta <- n - (r + 1)
  } else if (r == n) {
    alpha <- r - 1
    beta <- n - (r - 1)
  } else {
    alpha <- r
    beta <- n - r
  }

  return(c(alpha = alpha, beta = beta))
}
