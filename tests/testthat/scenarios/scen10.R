# No censoring
# No Noise
# All Events observed
# Homogeneous Historical Data
# Heavy Prior Data conflict
# Hazard ratio in fovor of control group


# Scen10
scen10 <- SimTestData(
  SimStudy_nPat = c(g1 = 200, g2 = 200),
  SimStudy_hz = c(g1 = 0.2, g2 = 0.1),
  SimStudy_dropout = c(rate = 0.05, time = 12),
  SimStudy_accr = 6,
  SimStudy_accr_method = "Uniform",
  SimStudy_surv_method = "Exponential",
  SimStudy_intensity = NA,
  SimStudy_accr_timepoint = NA,
  SimStudy_time_cutoff = 18,
  SimStudy_NObsEvt = 93,
  SimStudy_censor_type = 2,
  nStudy = 6,
  tau = 0.02,
  prior_data_conflict = TRUE,
  pdc_hz = 1.2,
  SAF_TOPIC = "Scen10",
  seed = 1701673095
)


write.csv(scen10, file = "Scen10.csv")
