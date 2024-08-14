#' Simulate Test Data Set
#'
#' @param SimStudy_nPat
#' @param SimStudy_hz
#' @param SimStudy_dropout
#' @param SimStudy_accr
#' @param SimStudy_accr_method
#' @param SimStudy_surv_method
#' @param SimStudy_intensity
#' @param SimStudy_accr_timepoint
#' @param SimStudy_time_cutoff
#' @param SimStudy_NObsEvt
#' @param SimStudy_censor_type
#' @param nStudy Number
#' @param tau
#' @param prior_data_conflict
#' @param SAF_TOPIC Selected safety topic to analyze/the adverse event of interest
#' @param pdc_hz
#' @param diff_trt_length
#' @param seed
#'
#' @return
#' @export
#'
#' @examples
SimTestData <- function(
    SimStudy_nPat = c(g1 = 50, g2 = 100),
    SimStudy_hz = c(g1 = 0.1, g2 = 0.2),
    SimStudy_dropout = c(rate = 0.05, time = 18),
    SimStudy_accr = 6,
    SimStudy_accr_method = "Uniform",
    SimStudy_surv_method = "Exponential",
    SimStudy_intensity = c(2, 4, 6),
    SimStudy_accr_timepoint = c(0, 2, 4, 6),
    SimStudy_time_cutoff = 18,
    SimStudy_NObsEvt = 100,
    SimStudy_censor_type = 1,
    nStudy = 5,
    tau = 0,
    prior_data_conflict = FALSE,
    diff_trt_length = FALSE,
    pdc_hz = c(g1 = 0.05, g2 = 0.5),
    SAF_TOPIC = "Example",
    seed = 123) {
  res <- array(
    data = NA, dim = c(nStudy, 5, 2),
    dimnames = list(
      STUDYID = c(1:nStudy),
      c("HIST", "ARM", "N", "N_WITH_AE", "TOT_EXP"),
      c("g1", "g2")
    )
  )

  res[1:(nStudy - 1), "HIST", ] <- 1
  res[nStudy, "HIST", ] <- 0

  res[, "ARM", "g1"] <- 1
  res[, "ARM", "g2"] <- 2

  res[, "N", "g1"] <- SimStudy_nPat["g1"]
  res[, "N", "g2"] <- SimStudy_nPat["g2"]

  # intiualize the list to save the data
  res_SimStudy <- list()

  #
  if (!is.na(seed)) {
    set.seed(seed)
  }

  # For prior Data conflict, simulate n-1 similar and 1 different trial
  if (prior_data_conflict == TRUE) {
    nStudy <- nStudy - 1
  }

  # Simulate sutdies
  SimStudy_time_cutoff_set <- (c(0.5, 0.5, 1, 1, 1.5, 1.5) + 1) * 12
  for (i in 1:nStudy) {
    if (tau > 0) {
      SimStudy_hz <- exp(log(SimStudy_hz) + rnorm(2, mean = 0, sd = tau))
    }
    if (diff_trt_length == FALSE) {
      res_SimStudy[[i]] <- SimStudy(
        nPat = SimStudy_nPat,
        hz = SimStudy_hz,
        dropout = SimStudy_dropout,
        accr = SimStudy_accr,
        accr_method = SimStudy_accr_method,
        surv_method = SimStudy_surv_method,
        intensity = SimStudy_intensity,
        accr_timepoint = SimStudy_accr_timepoint,
        time_cutoff = SimStudy_time_cutoff_set,
        NObsEvt = SimStudy_NObsEvt,
        censor_type = SimStudy_censor_type
      )
    } else {
      # treatment_length = runif(1, 0.5, 1.5)
      # SimStudy_time_cutoff = (treatment_length + 1) * 12
      # SimStudy_time_cutoff_set = c(SimStudy_time_cutoff_set, treatment_length)
      res_SimStudy[[i]] <- SimStudy(
        nPat = SimStudy_nPat,
        hz = SimStudy_hz,
        dropout = SimStudy_dropout,
        accr = SimStudy_accr,
        accr_method = SimStudy_accr_method,
        surv_method = SimStudy_surv_method,
        intensity = SimStudy_intensity,
        accr_timepoint = SimStudy_accr_timepoint,
        time_cutoff = SimStudy_time_cutoff_set[i],
        NObsEvt = SimStudy_NObsEvt,
        censor_type = SimStudy_censor_type
      )
    }
  }

  # Simulate the different trial
  if (prior_data_conflict == TRUE) {
    nStudy <- nStudy + 1

    res_SimStudy[[nStudy]] <- SimStudy(
      nPat = SimStudy_nPat,
      hz = pdc_hz,
      dropout = SimStudy_dropout,
      accr = SimStudy_accr,
      accr_method = SimStudy_accr_method,
      surv_method = SimStudy_surv_method,
      intensity = SimStudy_intensity,
      accr_timepoint = SimStudy_accr_timepoint,
      time_cutoff = SimStudy_time_cutoff,
      NObsEvt = SimStudy_NObsEvt,
      censor_type = SimStudy_censor_type
    )
  }

  for (s in 1:nStudy) {
    for (g in 1:2) {
      res[s, "TOT_EXP", g] <-
        sum(res_SimStudy[[s]][res_SimStudy[[s]]$gID == g, ]$ObsTime)

      res[s, "N_WITH_AE", g] <-
        sum(res_SimStudy[[s]][res_SimStudy[[s]]$gID == g, ]$EventIndicator)
    }
  }

  res_df <- as.data.frame(rbind(res[, , 1], res[, , 2]))
  row.names(res_df) <- c(paste0(c(1:nStudy), "_g1"), paste0(c(1:nStudy), "_g2"))
  res_df$STUDYID <- c(paste0("Study#", 1:nStudy), paste0("Study#", 1:nStudy))
  res_df[res_df$ARM == 1, "ARM"] <- "g1"
  res_df[res_df$ARM == 2, "ARM"] <- "g2"
  res_df$SAF_TOPIC <- SAF_TOPIC
  res_df <- res_df[, c(
    "STUDYID", "HIST", "ARM", "N",
    "SAF_TOPIC", "N_WITH_AE", "TOT_EXP"
  )]


  if (diff_trt_length == TRUE) {
    res_df$LENGTH <- NA
    for (i in 1:nStudy) {
      res_df[res_df$STUDYID == paste0("Study#", i), ]$LENGTH <- round(SimStudy_time_cutoff_set[i] / 12 * 365)
    }
  }

  res_df$TREAT <- SAF_TOPIC

  return(res_df)
}
