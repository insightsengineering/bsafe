# Overview ----------------------------------------------------------------

testing_list_prop <- list()

# any Scenario ------------------------------------------------------------

# set up testing parameters for model summary display proportions
tp_prop_Any <- data.frame(
  csv        = "any.csv",
  group      = "any",
  analysis   = "Incidence proportion",
  saf_topic  = "any",
  seed       = 1699874539,
  pool       = TRUE,
  tau        = "HalfNormal",
  heterog    = "any",
  ESS        = "elir",
  rob_weight = 0.2, # needs change
  nta_event  = 50, # initially NA
  nta_npat   = 100 # initially NA
)

thresholds_prop_Any <- data.frame(
  mean_min = rep(0, 4),
  mean_lb = rep(0, 4),
  mean_ub = rep(1, 4),
  mean_max = rep(1, 4),
  sd_min = rep(0, 4),
  sd_lb = rep(0, 4),
  sd_ub = rep(Inf, 4),
  sd_max = rep(Inf, 4),
  median_min = rep(0, 4),
  median_lb = rep(0, 4),
  median_ub = rep(1, 4),
  median_max = rep(1, 4),
  cri_lb_min = rep(0, 4),
  cri_lb_lb = rep(0, 4),
  cri_lb_ub = rep(1, 4),
  cri_lb_max = rep(1, 4),
  cri_ub_min = rep(0, 4),
  cri_ub_lb = rep(0, 4),
  cri_ub_ub = rep(1, 4),
  cri_ub_max = rep(1, 4),
  ess_min = c(0, 4, NA, NA),
  ess_lb = c(0, 4, NA, NA),
  ess_ub = c(Inf, Inf, NA, NA),
  ess_max = c(Inf, Inf, NA, NA),
  row.names = c("MAP Prior", "Robustified", "Likelihood", "Posterior")
)


thresholds_prop_NA <- thresholds_prop_Any
thresholds_prop_NA[] <- NA

testing_list_prop[[1]] <- list(
  parameters = tp_prop_Any,
  treshholds = thresholds_prop_Any
)

# initializing all treshhold data frames
for (i in 1:13) {
  if (i < 10) {
    i_chr <- paste0(0, i)
  } else {
    i_chr <- i
  }
  assign(paste0("thresholds_prop_Scen", i_chr), thresholds_prop_NA)
}

# Deleting the initial ones
rm(tp_prop_Any, thresholds_prop_Any, thresholds_prop_NA)



# Scen 1 ------------------------------------------------------------------

tp_prop_Scen01 <- data.frame(
  csv        = "Scen01.csv",
  group      = "g1",
  analysis   = "Incidence proportion",
  saf_topic  = "Scen01",
  seed       = 1699874539,
  pool       = TRUE,
  tau        = "HalfNormal",
  heterog    = "Small",
  ESS        = "elir",
  rob_weight = 0.05,
  nta_event  = 194,
  nta_npat   = 200,
  row.names  = c("Best case scenario")
)

testing_list_prop[[2]] <- list(
  parameters = tp_prop_Scen01,
  treshholds = thresholds_prop_Scen01
)

# Scen 2 ------------------------------------------------------------------

tp_prop_Scen02 <- data.frame(
  csv        = "Scen02.csv",
  group      = "g1",
  analysis   = "Incidence proportion",
  saf_topic  = "Scen02",
  seed       = 1701611344,
  pool       = TRUE,
  tau        = "HalfNormal",
  heterog    = "Moderate",
  ESS        = "elir",
  rob_weight = 0.8,
  nta_event  = 199,
  nta_npat   = 200,
  row.names  = c("Strong Prior Data Conflict")
)


testing_list_prop[[3]] <- list(
  parameters = tp_prop_Scen02,
  treshholds = thresholds_prop_Scen02
)

# Scen 3 ------------------------------------------------------------------

tp_prop_Scen03 <- data.frame(
  csv        = "Scen03.csv",
  group      = "g1",
  analysis   = "Incidence proportion",
  saf_topic  = "Scen03",
  seed       = 1701621384,
  pool       = TRUE,
  tau        = "HalfNormal",
  heterog    = "Substantial",
  ESS        = "elir",
  rob_weight = 0.25,
  nta_event  = 31,
  nta_npat   = 200,
  row.names  = c("Realisitic Scenarios")
)


testing_list_prop[[4]] <- list(
  parameters = tp_prop_Scen03,
  treshholds = thresholds_prop_Scen03
)


# Scen 4 ------------------------------------------------------------------

tp_prop_Scen04 <- data.frame(
  csv        = "Scen04.csv",
  group      = "g1",
  analysis   = "Incidence proportion",
  saf_topic  = "Scen04",
  seed       = 1701626683,
  pool       = TRUE,
  tau        = "HalfNormal",
  heterog    = "Very Large",
  ESS        = "elir",
  rob_weight = 0.99,
  nta_event  = 27,
  nta_npat   = 50,
  row.names  = c("Worst Case Scenario")
)

testing_list_prop[[5]] <- list(
  parameters = tp_prop_Scen04,
  treshholds = thresholds_prop_Scen04
)


# Scen 5 ------------------------------------------------------------------

tp_prop_Scen05 <- data.frame(
  csv        = "Scen05.csv",
  group      = "g1",
  analysis   = "Incidence proportion",
  saf_topic  = "Scen05",
  seed       = 1701628373,
  pool       = TRUE,
  tau        = "HalfNormal",
  heterog    = "Large",
  ESS        = "elir",
  rob_weight = 0.4,
  nta_event  = 25,
  nta_npat   = 200,
  row.names  = c("Heterogenous Data ")
)

testing_list_prop[[6]] <- list(
  parameters = tp_prop_Scen05,
  treshholds = thresholds_prop_Scen05
)

# Scen 6 ------------------------------------------------------------------

tp_prop_Scen06 <- data.frame(
  csv        = "Scen06.csv",
  group      = "g1",
  analysis   = "Incidence proportion",
  saf_topic  = "Scen06",
  seed       = 1701628373,
  pool       = TRUE,
  tau        = "HalfNormal",
  heterog    = "Moderate",
  ESS        = "elir",
  rob_weight = 0.14,
  nta_event  = 31,
  nta_npat   = 200,
  row.names  = c("High Dropout")
)

testing_list_prop[[7]] <- list(
  parameters = tp_prop_Scen06,
  treshholds = thresholds_prop_Scen06
)

# Scen 7 ------------------------------------------------------------------

tp_prop_Scen07 <- data.frame(
  csv        = "Scen07.csv",
  group      = "g1",
  analysis   = "Incidence proportion",
  saf_topic  = "Scen07",
  seed       = 1701416989,
  pool       = TRUE,
  tau        = "HalfNormal",
  heterog    = "Very Large",
  ESS        = "elir",
  rob_weight = 0.2,
  nta_event  = 35,
  nta_npat   = 200,
  row.names  = c("High Heterogeneity")
)

testing_list_prop[[8]] <- list(
  parameters = tp_prop_Scen07,
  treshholds = thresholds_prop_Scen07
)

# Scen 8 ------------------------------------------------------------------

tp_prop_Scen08 <- data.frame(
  csv        = "Scen08.csv",
  group      = "g1",
  analysis   = "Incidence proportion",
  saf_topic  = "Scen08",
  seed       = 1701652217,
  pool       = TRUE,
  tau        = "HalfNormal",
  heterog    = "Large",
  ESS        = "elir",
  rob_weight = 0.2,
  nta_event  = 25,
  nta_npat   = 200,
  row.names  = c("Bad Scenario")
)

testing_list_prop[[9]] <- list(
  parameters = tp_prop_Scen08,
  treshholds = thresholds_prop_Scen08
)

# Scen 9 ------------------------------------------------------------------

tp_prop_Scen09 <- data.frame(
  csv        = "Scen09.csv",
  group      = "g1",
  analysis   = "Incidence proportion",
  saf_topic  = "Scen09",
  seed       = 1701655293,
  pool       = TRUE,
  tau        = "HalfNormal",
  heterog    = "Small",
  ESS        = "elir",
  rob_weight = 0.05,
  nta_event  = 175,
  nta_npat   = 200,
  row.names  = c("Good Scenario")
)

testing_list_prop[[10]] <- list(
  parameters = tp_prop_Scen09,
  treshholds = thresholds_prop_Scen09
)

# Scen 10 ------------------------------------------------------------------

tp_prop_Scen10 <- data.frame(
  csv        = "Scen10.csv",
  group      = "g1",
  analysis   = "Incidence proportion",
  saf_topic  = "Scen10",
  seed       = 1701673095,
  pool       = TRUE,
  tau        = "HalfNormal",
  heterog    = "Small",
  ESS        = "elir",
  rob_weight = 0.6,
  nta_event  = 175,
  nta_npat   = 200,
  row.names  = c("Favored Control")
)

testing_list_prop[[11]] <- list(
  parameters = tp_prop_Scen10,
  treshholds = thresholds_prop_Scen10
)

# Scen 11 ------------------------------------------------------------------

tp_prop_Scen11 <- data.frame(
  csv        = "Scen11.csv",
  group      = "g1",
  analysis   = "Incidence proportion",
  saf_topic  = "Scen11",
  seed       = 1701876972,
  pool       = TRUE,
  tau        = "HalfNormal",
  heterog    = "Small",
  ESS        = "elir",
  rob_weight = 0.05,
  nta_event  = 170,
  nta_npat   = 200,
  row.names  = c("Continued Study Duration with Realistic Setting")
)

testing_list_prop[[12]] <- list(
  parameters = tp_prop_Scen11,
  treshholds = thresholds_prop_Scen11
)

# Scen 12 ------------------------------------------------------------------

tp_prop_Scen12 <- data.frame(
  csv        = "Scen12.csv",
  group      = "g1",
  analysis   = "Incidence proportion",
  saf_topic  = "Scen12",
  seed       = 1701878308,
  pool       = TRUE,
  tau        = "HalfNormal",
  heterog    = "Large", # maybe should be Small
  ESS        = "elir",
  rob_weight = 0.5,
  nta_event  = 30,
  nta_npat   = 200,
  row.names  = c("Continued Study Duration with Worst Setting")
)

testing_list_prop[[13]] <- list(
  parameters = tp_prop_Scen12,
  treshholds = thresholds_prop_Scen12
)

# Scen 13 ------------------------------------------------------------------

tp_prop_Scen13 <-  data.frame(
  csv        = "Scen13.csv",
  group      = "g1",
  analysis   = "Incidence proportion",
  saf_topic  = "Scen13",
  seed       = 1718356066,
  pool       = TRUE,
  tau        = "HalfNormal",
  heterog    = "Large",
  ESS        = "elir",
  rob_weight = 0.1,
  nta_event  = 186,
  nta_npat   = 200,
  row.names  = c("Different treatment length")
)

testing_list_prop[[14]] <- list(parameters = tp_prop_Scen13,
                                      treshholds = thresholds_prop_Scen13)

testing_list_props <- testing_list_prop
# Deleting obsolote data frames -------------------------------------------

for (i in 1:13) {
  if (i < 10) {
    i_chr <- paste0(0, i)
  } else {
    i_chr <- i
  }
  rm(list = paste0("tp_prop_Scen", i_chr))
  rm(list = paste0("thresholds_prop_Scen", i_chr))
}
rm(i)
rm(i_chr)
rm(testing_list_prop)
