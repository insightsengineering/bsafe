# Function to simulate 1 study
# nPat = Number of patients in each group
# g1 = group 1 (treatment); g2 = group 2 (control)
# dropout = 0_05: 5% dropout after time units of measure
# accr = accrual time, is to be in regards to the hazard
# NObsEvt = type 2 censoring, censor after NObsEvt number of events, probability of observing the event if < 1_
# accr_timepoint should include 0 and total accrual time_
# Pre-specify the censor type ahead of time


#' Simulate a study
#'
#' @param nPat
#' @param hz
#' @param dropout
#' @param accr
#' @param NObsEvt
#' @param accr_method
#' @param surv_method
#' @param intensity
#' @param accr_timepoint
#' @param censor_type
#' @param time_cutoff
#'
#' @return
#' @export
#'
#' @examples
SimStudy <- function(nPat = c(g1 = 100, g2 = 100),
                     hz = c(g1 = 0.1, g2 = 0.2),
                     dropout = c(rate = 0.05, time = 12),
                     accr = 6,
                     NObsEvt = 0.5,
                     accr_method = "Uniform",
                     surv_method = "Exponential",
                     intensity = c(2, 6, 10),
                     accr_timepoint = c(0, 2, 4, 6),
                     censor_type = 1,
                     time_cutoff = 18) {
  N <- sum(nPat)
  # Observed events either proportional ( < 100) or as absolute numbers
  if (NObsEvt < 1) {
    NObsEvt <- sum(nPat) * NObsEvt
  }

  # res: variable which stores the result output
  # gID: 1 treatment 2 control
  # ID: Subject
  # Entry: Entry time according to accrual
  # EventTime: Simulated Eventtime + Entry time
  # ObsTime: Time observed (min(EventTime, CensorTime)-Entry)
  # StudyTime: Timepoint in Study
  # Eventindicator: 1 event observed, 0 censored
  res <- matrix(
    data = NA, nrow = N, ncol = 8,
    dimnames = list(
      ID = 1:N,
      c(
        "gID", "ID", "Entry", "EventTime",
        "ObsTime", "CensorTime",
        "StudyTime", "EventIndicator"
      )
    )
  )
  # ID and gID just from 1 to number of patients in each group
  res[, "ID"] <- 1:N
  res[, "gID"] <- c(rep(1, nPat["g1"]), rep(2, nPat["g2"]))

  # Different methods for generating Enrollment Time
  if (accr_method == "Uniform") {
    res[, "Entry"] <- runif(N, 0, accr)
  }

  # Poisson accrual times
  if (accr_method == "Poisson") {
    rtlist <- lapply(intensity, function(x) rexp(N, x))
    recruit_time <- c()
    for (i in 1:length(intensity)) {
      recruit_time_new <- c(accr_timepoint[i] + cumsum(rtlist[[i]][(accr_timepoint[i] + cumsum(rtlist[[i]])) < accr_timepoint[i + 1]]))
      recruit_time <- c(recruit_time, recruit_time_new)
    }
    if (length(recruit_time) < N) {
      enrollment <- c(recruit_time, runif((N - length(recruit_time)), min(accr_timepoint), max(accr_timepoint)))
    } else {
      enrollment <- recruit_time[1:N]
    }
    res[, "Entry"] <- enrollment
  }

  # Piecewise Uniform accrual times
  if (accr_method == "Piecewise Uniform") {
    recruit_time <- c()
    for (i in 1:length(intensity)) {
      n_part <- intensity[i] * diff(accr_timepoint)[i]
      recruit_time_new <- runif(n_part, accr_timepoint[i], accr_timepoint[i + 1])
      recruit_time <- c(recruit_time, recruit_time_new)
    }
    if (length(recruit_time) < N) {
      enrollment <- c(recruit_time, runif((N - length(recruit_time)), min(accr_timepoint), max(accr_timepoint)))
    } else {
      enrollment <- recruit_time[1:N]
    }
    res[, "Entry"] <- enrollment
  }

  # Method for generating Survival Time
  if (surv_method == "Exponential") {
    for (i in 1:length(nPat)) {
      SurvTimesG <- rexp(nPat[i], hz[i])
      if (i == 1) {
        SurvTimes <- SurvTimesG
      } else {
        SurvTimes <- c(SurvTimes, SurvTimesG)
      }
    }
  }


  # Event Times
  res[, "EventTime"] <- res[, "Entry"] + SurvTimes

  # Get rate parameter for exponential distributed censoring times
  CensorRate <- if (dropout["rate"] > 0) {
    -log(1 - dropout["rate"]) / dropout["time"]
  } else {
    0
  }

  # Censoring times for all individuals, infinity if no censoring is applied
  CensorTime <- if (dropout["rate"] > 0) {
    rexp(N, CensorRate)
  } else {
    rep(Inf, N)
  }

  res[, "CensorTime"] <- CensorTime + res[, "Entry"]

  # Censor type 1, administrative censoring after cutoff time
  if (censor_type == 1) {
    evt_ind <- which(res[, "EventTime"] < res[, "CensorTime"] & res[, "EventTime"] < time_cutoff)
    non_evt_ind <- which(!(res[, "EventTime"] < res[, "CensorTime"] & res[, "EventTime"] < time_cutoff))
    res[evt_ind, "EventIndicator"] <- 1
    res[non_evt_ind, "EventIndicator"] <- 0
    res[evt_ind, "ObsTime"] <- res[evt_ind, "EventTime"] - res[evt_ind, "Entry"]
    res[non_evt_ind, "ObsTime"] <- ifelse(res[non_evt_ind, "CensorTime"] < time_cutoff,
      res[non_evt_ind, "CensorTime"] - res[non_evt_ind, "Entry"],
      time_cutoff - res[non_evt_ind, "Entry"]
    )
    res[, "StudyTime"] <- res[, "ObsTime"] + res[, "Entry"]
  }


  # Type 2 censoring, censoring after number of observed events
  if (censor_type == 2) {
    # Introduce censoring indices
    evt_ind <- which(res[, "EventTime"] < res[, "CensorTime"])
    non_evt_ind <- which(res[, "EventTime"] >= res[, "CensorTime"])
    res[evt_ind, "EventIndicator"] <- 1
    res[non_evt_ind, "EventIndicator"] <- 0
    res[evt_ind, "ObsTime"] <- res[evt_ind, "EventTime"] - res[evt_ind, "Entry"]
    res[non_evt_ind, "ObsTime"] <- res[non_evt_ind, "CensorTime"] - res[non_evt_ind, "Entry"]
    res[, "StudyTime"] <- res[, "ObsTime"] + res[, "Entry"]

    type2_censortime <- sort(res[, "StudyTime"], decreasing = FALSE)[NObsEvt]
    type2_censorind <- which(res[, "StudyTime"] > type2_censortime)
    res[type2_censorind, "StudyTime"] <- type2_censortime
    res[type2_censorind, "EventIndicator"] <- 0

    new_censored_row_idx <- which(res[, "StudyTime"] == type2_censortime)

    res[new_censored_row_idx, "ObsTime"] <- type2_censortime - res[new_censored_row_idx, "Entry"]
  }

  res <- as.data.frame(res)
  return(res)
}
