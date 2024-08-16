# Overview ----------------------------------------------------------------

testing_list_rates <- list()


# set up testing parameters for model summary display Exposure-adjusted AE rate
tp_rate_Any <- data.frame(
  csv        = "Scen01.csv",
  group      = "g1",
  analysis   = "Exposure-adjusted AE rate",
  saf_topic  = "Scen1",
  seed       = 1699874539,
  pool       = TRUE,
  tau        = "HalfNormal",
  heterog    = "Small",
  ESS        = "elir",
  rob_weight = 0.05,
  rob_mean   = 0.1,
  nta_event  = 100,
  nta_time   = 1000
)

thresholds_rate_Any <- data.frame(
  mean_min = c(rep(-Inf, 4), rep(0, 4)),
  mean_lb = c(rep(-Inf, 4), rep(0, 4)),
  mean_ub = c(rep(Inf, 4), rep(Inf, 4)),
  mean_max = c(rep(Inf, 4), rep(Inf, 4)),
  sd_min = c(rep(0, 4), rep(0, 4)),
  sd_lb = c(rep(0, 4), rep(0, 4)),
  sd_ub = c(rep(Inf, 4), rep(Inf, 4)),
  sd_max = c(rep(Inf, 4), rep(Inf, 4)),
  median_min = c(rep(-Inf, 4), rep(0, 4)),
  median_lb = c(rep(-Inf, 4), rep(0, 4)),
  median_ub = c(rep(Inf, 4), rep(Inf, 4)),
  median_max = c(rep(Inf, 4), rep(Inf, 4)),
  cri_lb_min = c(rep(-Inf, 4), rep(0, 4)),
  cri_lb_lb = c(rep(-Inf, 4), rep(0, 4)),
  cri_lb_ub = c(rep(Inf, 4), rep(Inf, 4)),
  cri_lb_max = c(rep(Inf, 4), rep(Inf, 4)),
  cri_ub_min = c(rep(-Inf, 4), rep(0, 4)),
  cri_ub_lb = c(rep(-Inf, 4), rep(0, 4)),
  cri_ub_ub = c(rep(Inf, 4), rep(Inf, 4)),
  cri_ub_max = c(rep(Inf, 4), rep(Inf, 4)),
  ess_min = c(rep(0, 2), rep(NA, 6)),
  ess_lb = c(rep(0, 2), rep(NA, 6)),
  ess_ub = c(rep(Inf, 2), rep(NA, 6)),
  ess_max = c(rep(Inf, 2), rep(NA, 6)),
  row.names = c("log MAP Prior", "log Robustified", "log Likelihood", "log Posterior", "exp MAP Prior", "exp Robustified", "exp Likelihood", "exp Posterior")
)

thresholds_rate_NA <- thresholds_rate_Any
thresholds_rate_NA[] <- NA

testing_list_rates[[1]] <- list(
  parameters = tp_rate_Any,
  treshholds = thresholds_rate_Any
)


# initializing all treshhold data frames
for (i in 1:14) {
  if (i < 10) {
    i_chr <- paste0(0, i)
  } else {
    i_chr <- i
  }
  assign(paste0("thresholds_rate_Scen", i_chr), thresholds_rate_NA)
}

# Deleting the initial ones
rm(tp_rate_Any, thresholds_rate_Any, thresholds_rate_NA)


# Scen1 -------------------------------------------------------------------

tp_rate_Scen01 <- data.frame(
  csv        = "Scen01.csv",
  group      = "g1",
  analysis   = "Exposure-adjusted AE rate",
  saf_topic  = "Scen01",
  seed       = 1699874539,
  pool       = TRUE,
  tau        = "HalfNormal",
  heterog    = "Small",
  ESS        = "elir",
  rob_weight = 0.05,
  rob_mean   = 0.1,
  nta_event  = 100,
  nta_time   = 1000,
  row.names  = c("Best Case Scenario")
)

testing_list_rates[[2]] <- list(
  parameters = tp_rate_Scen01,
  treshholds = thresholds_rate_Scen01
)

# Scen 2 ------------------------------------------------------------------

tp_rate_Scen02 <- data.frame(
  csv        = "Scen02.csv",
  group      = "g2",
  analysis   = "Exposure-adjusted AE rate",
  saf_topic  = "Scen02",
  seed       = 1701611344,
  pool       = TRUE,
  tau        = "HalfNormal",
  heterog    = "Small",
  ESS        = "elir",
  rob_weight = 0.8,
  rob_mean   = 0.3854,
  nta_event  = 200,
  nta_time   = 518,
  row.names  = c("Strong Prior Data Conflict")
)

testing_list_rates[[3]] <- list(
  parameters = tp_rate_Scen02,
  treshholds = thresholds_rate_Scen02
)

# Scen 3 ------------------------------------------------------------------

tp_rate_Scen03 <- data.frame(
  csv        = "Scen03.csv",
  group      = "g1",
  analysis   = "Exposure-adjusted AE rate",
  saf_topic  = "Scen03",
  seed       = 1701621384,
  pool       = TRUE,
  tau        = "HalfNormal",
  heterog    = "Substantial",
  ESS        = "elir",
  rob_weight = 0.25,
  rob_mean   = 0.0944,
  nta_event  = 31,
  nta_time   = 328.47,
  row.names  = c("Realistic scenario")
)

testing_list_rates[[4]] <- list(
  parameters = tp_rate_Scen03,
  treshholds = thresholds_rate_Scen03
)


# Scen 4 ------------------------------------------------------------------

tp_rate_Scen04 <- data.frame(
  csv        = "Scen04.csv",
  group      = "g1",
  analysis   = "Exposure-adjusted AE rate",
  saf_topic  = "Scen04",
  seed       = 1701626683,
  pool       = TRUE,
  tau        = "HalfNormal",
  heterog    = "Very Large",
  ESS        = "elir",
  rob_weight = 0.25,
  rob_mean   = 0.0539,
  nta_event  = 31,
  nta_time   = 328.47,
  row.names  = c("worst Case Scenario")
)


testing_list_rates[[5]] <- list(
  parameters = tp_rate_Scen04,
  treshholds = thresholds_rate_Scen04
)


# Scen 5 ------------------------------------------------------------------

tp_rate_Scen05 <- data.frame(
  csv        = "Scen05.csv",
  group      = "g1",
  analysis   = "Exposure-adjusted AE rate",
  saf_topic  = "Scen05",
  seed       = 1701628373,
  pool       = TRUE,
  tau        = "HalfNormal",
  heterog    = "Large",
  ESS        = "elir",
  rob_weight = 0.2,
  rob_mean   = 0.0865,
  nta_event  = 25,
  nta_time   = 289,
  row.names  = c("Heterogeneity Data (Medium)")
)

testing_list_rates[[6]] <- list(
  parameters = tp_rate_Scen05,
  treshholds = thresholds_rate_Scen05
)


# Scen 6 ------------------------------------------------------------------

tp_rate_Scen06 <- data.frame(
  csv        = "Scen06.csv",
  group      = "g1",
  analysis   = "Exposure-adjusted AE rate",
  saf_topic  = "Scen06",
  seed       = 1701628373,
  pool       = TRUE,
  tau        = "HalfNormal",
  heterog    = "Moderate",
  ESS        = "elir",
  rob_weight = 0.14,
  rob_mean   = 0.1204,
  nta_event  = 31,
  nta_time   = 257,
  row.names  = c("High Dropout")
)

testing_list_rates[[7]] <- list(
  parameters = tp_rate_Scen06,
  treshholds = thresholds_rate_Scen06
)


# Scen 7 ------------------------------------------------------------------

tp_rate_Scen07 <- data.frame(
  csv        = "Scen07.csv",
  group      = "g1",
  analysis   = "Exposure-adjusted AE rate",
  saf_topic  = "Scen07",
  seed       = 1701416989,
  pool       = TRUE,
  tau        = "HalfNormal",
  heterog    = "Very Large",
  ESS        = "elir",
  rob_weight = 0.2,
  rob_mean   = 0.19,
  nta_event  = 35,
  nta_time   = 200,
  row.names  = c("High Heterogenity")
)


testing_list_rates[[8]] <- list(
  parameters = tp_rate_Scen07,
  treshholds = thresholds_rate_Scen07
)

# Scen 8 ------------------------------------------------------------------

tp_rate_Scen08 <- data.frame(
  csv        = "Scen08.csv",
  group      = "g1",
  analysis   = "Exposure-adjusted AE rate",
  saf_topic  = "Scen08",
  seed       = 1701652217,
  pool       = TRUE,
  tau        = "HalfNormal",
  heterog    = "Large",
  ESS        = "elir",
  rob_weight = 0.2,
  rob_mean   = 0.0741,
  nta_event  = 35,
  nta_time   = 338,
  row.names  = c("Bad Scenario")
)

testing_list_rates[[9]] <- list(
  parameters = tp_rate_Scen08,
  treshholds = thresholds_rate_Scen08
)

# Scen 9 ------------------------------------------------------------------

tp_rate_Scen09 <- data.frame(
  csv        = "Scen09.csv",
  group      = "g1",
  analysis   = "Exposure-adjusted AE rate",
  saf_topic  = "Scen09",
  seed       = 1701655293,
  pool       = TRUE,
  tau        = "HalfNormal",
  heterog    = "Small",
  ESS        = "elir",
  rob_weight = 0.05,
  rob_mean   = 0.0926,
  nta_event  = 92,
  nta_time   = 1000,
  row.names  = c("Good Scenario")
)

testing_list_rates[[10]] <- list(
  parameters = tp_rate_Scen09,
  treshholds = thresholds_rate_Scen09
)

# Scen 10 ------------------------------------------------------------------

tp_rate_Scen10 <- data.frame(
  csv        = "Scen10.csv",
  group      = "g1",
  analysis   = "Exposure-adjusted AE rate",
  saf_topic  = "Scen10",
  seed       = 1701673095,
  pool       = TRUE,
  tau        = "HalfNormal",
  heterog    = "Small",
  ESS        = "elir",
  rob_weight = 0.6,
  rob_mean   = 0.2472,
  nta_event  = 150,
  nta_time   = 200,
  row.names  = c("Favored Control")
)

testing_list_rates[[11]] <- list(
  parameters = tp_rate_Scen10,
  treshholds = thresholds_rate_Scen10
)

# Scen 11 ------------------------------------------------------------------

tp_rate_Scen11 <- data.frame(
  csv        = "Scen11.csv",
  group      = "g1",
  analysis   = "Exposure-adjusted AE rate",
  saf_topic  = "Scen11",
  seed       = 1701876972,
  pool       = TRUE,
  tau        = "HalfNormal",
  heterog    = "Small",
  ESS        = "elir",
  rob_weight = 0.05,
  rob_mean   = 0.0952,
  nta_event  = 95,
  nta_time   = 1000,
  row.names  = c("Continued study duration with Realistic Setting")
)

testing_list_rates[[12]] <- list(
  parameters = tp_rate_Scen11,
  treshholds = thresholds_rate_Scen11
)


# Scen 12 ------------------------------------------------------------------

tp_rate_Scen12 <- data.frame(
  csv        = "Scen12.csv",
  group      = "g1",
  analysis   = "Exposure-adjusted AE rate",
  saf_topic  = "Scen12",
  seed       = 1701878308,
  pool       = TRUE,
  tau        = "HalfNormal",
  heterog    = "Large",
  ESS        = "elir",
  rob_weight = 0.5,
  rob_mean   = 0.2,
  nta_event  = 200,
  nta_time   = 1000,
  row.names  = c("Continued study duration with Worst Setting")
)

testing_list_rates[[13]] <- list(
  parameters = tp_rate_Scen12,
  treshholds = thresholds_rate_Scen12
)

# Scen 13 ------------------------------------------------------------------

tp_rate_Scen13 <- data.frame(
  csv        = "Scen13.csv",
  group      = "g1",
  analysis   = "Exposure-adjusted AE rate",
  saf_topic  = "Scen13",
  seed       = 1718356066,
  pool       = TRUE,
  tau        = "HalfNormal",
  heterog    = "Large",
  ESS        = "elir",
  rob_weight = 0.1,
  rob_mean   = 0.2,
  nta_event  = 186,
  nta_time   = 1681,
  row.names  = c("Different Study length")
)

testing_list_rates[[14]] <- list(
  parameters = tp_rate_Scen13,
  treshholds = thresholds_rate_Scen13
)


# Deleting obsolete data frames -------------------------------------------

for (i in 1:13) {
  if (i < 10) {
    i_chr <- paste0(0, i)
  } else {
    i_chr <- i
  }
  rm(list = paste0("tp_rate_Scen", i_chr))
  rm(list = paste0("thresholds_rate_Scen", i_chr))
}
rm(i)
rm(i_chr)
